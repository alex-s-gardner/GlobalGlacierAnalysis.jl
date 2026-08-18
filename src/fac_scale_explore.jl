# Explore methods to scale FAC (Firn Air Content) as a function of change in melt rate.
# Develops a surrogate model: FAC_scaled = f(mscale, reference_state, location, climate)
# Trained on natural variation across GEMB elevation_delta and pscale perturbations.

# =============================================================================
# MODELS TESTED (see git history for full implementation)
# =============================================================================
# 01_Linear: Multiple linear regression
#    - Test R² = 0.601, RMSE = 4.02 m (baseline, underfits)
# 02_Poly2_Ridge: Polynomial degree 2 + Ridge regularization
#    - Test R² = 0.718, RMSE = 3.38 m (better but still underfits)
# 03_Loess_univar: Local regression on melt_rate only
#    - Test R² = 0.675, RMSE = 3.63 m (univariate, misses covariates)
# 04_RandomForest: Random Forest with hyperparameter tuning
#    - Test R² = 0.967, RMSE = 1.18 m (BEST - selected)
# 05_TimeResolved: Random Forest on time-resolved data
#    - Test R² = 0.94 (high accuracy but more complex, not needed)
# 06_PiecewiseLinear: Piecewise linear by melt quantiles
#    - Test R² = 0.848, RMSE = 2.48 m (interpretable but RF better)
# 07_RatioModel: Model dFAC/dMelt ratio directly
#    - Test R² = 0.24 on ratio prediction (too indirect)
# 08_GradientBoosting: EvoTrees gradient boosting
#    - Skipped due to API issues
#
# SELECTED MODEL: Random Forest (R² = 0.965 on ΔFAC)
# =============================================================================

begin
    import GlobalGlacierAnalysis as GGA
    using Dates
    using FileIO
    using ProgressMeter
    using DimensionalData
    using CairoMakie
    using Statistics
    using LinearAlgebra
    using Random
    using StatsBase
    using JLD2
    using DecisionTree

    # run parameters
    gemb_run_id = 5;
    minimum_land_coverage_fraction = 0.70;

    # Subset for faster iteration
    use_subset = true  # Set to false for full dataset
    n_locations_subset = 1000  # Number of unique locations to use

    vars2extract = ["fac", "acc", "refreeze", "melt", "rain", "ec"]
    dims2extract = ["latitude", "longitude", "date", "height"]

    gembinfo = GGA.gemb_info(; gemb_run_id);

    # define date and height binning ranges
    date_range, date_center = GGA.project_date_bins()

    # expand daterange to 1970
    date_end_new = Date(1970, 1, 1)
    Δd = 30
    date_range = reverse(last(date_range):-Day(Δd):date_end_new)
    date_center = date_range[1:end-1] .+ Day(Δd / 2)
    ddate = Dim{:date}(date_center)
    length(date_range) == (length(date_center) + 1) || error("date_range and date_center length mismatch")

    # Cache file for extracted training features
    cache_file = joinpath(dirname(gembinfo.filename_gemb_combined), "fac_surrogate_training_data_delta.jld2")
end;

# =============================================================================
# DATA LOADING: Use cache if available, otherwise load from raw .mat files
# =============================================================================
if isfile(cache_file)
    println("Loading cached training data from: $cache_file")
    training_data = load(cache_file)
    X_fac_all = training_data["X_fac"]
    X_refreeze_all = training_data["X_refreeze"]
    y_fac_all = training_data["y_fac"]
    y_refreeze_all = training_data["y_refreeze"]
    feature_names_fac = training_data["feature_names_fac"]
    feature_names_refreeze = training_data["feature_names_refreeze"]
    lat_all = training_data["lat"]
    lon_all = training_data["lon"]
    elev_delta_all = training_data["elev_delta"]
    pscale_all = training_data["pscale"]
    println("  Loaded $(size(X_fac_all, 1)) samples with 2 targets and separate feature sets")
else
    println("No cache found. Loading raw GEMB data (this takes ~10-40 min depending on I/O)...")
    gemb_files = vcat(GGA.allfiles.(gembinfo.gemb_folder; subfolders=false, fn_endswith=".mat", fn_contains=gembinfo.file_uniqueid)...)

    # read_gemb_files uses Threads.@threads internally for parallel file loading
    gembX = GGA.read_gemb_files(gemb_files, gembinfo; vars2extract=vcat(dims2extract, vars2extract), date_range, date_center, path2land_sea_mask=GGA.pathlocal.era5_land_sea_mask, minimum_land_coverage_fraction)

    # Remove refrozen rain from rain; remove non-refreezing rain from accumulation
    RefrozenRain = min.(gembX["rain"], gembX["refreeze"])
    gembX["rain"] .-= RefrozenRain
    gembX["acc"] .-= gembX["rain"]

    # =========================================================================
    # FEATURE EXTRACTION
    # =========================================================================
    println("Extracting training features...")

    # Compute rates from cumulative variables (per 30-day bin)
    melt_rate = diff(gembX["melt"], dims=2)
    acc_rate = diff(gembX["acc"], dims=2)
    refreeze_rate = diff(gembX["refreeze"], dims=2)
    rain_rate = diff(gembX["rain"], dims=2)
    ec_rate = diff(gembX["ec"], dims=2)
    fac_rate = diff(gembX["fac"], dims=2)

    # --- Time-averaged features (one sample per grid point) ---
    # NaN-safe mean along dim 2
    function nanmean_rows(X::Matrix)
        n = size(X, 1)
        result = zeros(n)
        Threads.@threads for i in 1:n
            row = @view X[i, :]
            valid = .!isnan.(row)
            if any(valid)
                result[i] = mean(row[valid])
            else
                result[i] = NaN
            end
        end
        return result
    end

    mean_melt_rate = nanmean_rows(melt_rate)
    mean_acc_rate = nanmean_rows(acc_rate)
    mean_refreeze_rate = nanmean_rows(refreeze_rate)
    mean_rain_rate = nanmean_rows(rain_rate)
    mean_ec_rate = nanmean_rows(ec_rate)
    mean_fac_rate = nanmean_rows(fac_rate)
    mean_fac_absolute = nanmean_rows(gembX["fac"][:, 2:end])  # Mean absolute FAC

    # --- BUILD REFERENCE LOOKUP ---
    println("\nBuilding reference lookup (mscale=1 equivalent: elevation_delta=0, pscale=1)...")

    # Identify reference samples
    ref_mask = (gembX["elevation_delta"] .== 0.0) .& (gembX["precipitation_scale"] .== 1.0)
    n_total = length(mean_fac_rate)

    # Create location key (lat, lon) for matching
    loc_keys = [(gembX["latitude"][i], gembX["longitude"][i]) for i in 1:n_total]

    # Build reference lookup: location -> reference state
    ref_lookup = Dict{Tuple{Float64,Float64}, @NamedTuple{
        fac::Float64, melt_rate::Float64, fac_rate::Float64,
        acc_rate::Float64, refreeze_rate::Float64, rain_rate::Float64, ec_rate::Float64,
        idx::Int
    }}()

    for i in 1:n_total
        if ref_mask[i] && !isnan(mean_fac_absolute[i]) && !isnan(mean_melt_rate[i])
            ref_lookup[loc_keys[i]] = (
                fac = mean_fac_absolute[i],
                melt_rate = mean_melt_rate[i],
                fac_rate = mean_fac_rate[i],
                acc_rate = mean_acc_rate[i],
                refreeze_rate = mean_refreeze_rate[i],
                rain_rate = mean_rain_rate[i],
                ec_rate = mean_ec_rate[i],
                idx = i
            )
        end
    end

    println("  Found $(length(ref_lookup)) locations with valid reference cases")

    # --- CREATE TRAINING SAMPLES ---
    println("\nCreating training samples with reference-only features...")

    # Features: use RATIO features (better for ΔFAC rate prediction)
    # Don't use current refreeze/ec because they're unknown at production time
    melt_ratio_vec = Float64[]           # melt_rate / ref_melt_rate (= mscale in production)
    delta_melt_rate_vec = Float64[]      # ADDED: Δmelt_rate captures temperature change!
    delta_fac_rate_vec = Float64[]       # TARGET 1: ΔFAC rate = current - reference
    delta_refreeze_rate_vec = Float64[]  # TARGET 2: Δrefreeze rate = current - reference
    delta_acc_rate_vec = Float64[]       # Keep deltas for refreeze model
    delta_rain_rate_vec = Float64[]
    acc_ratio_vec = Float64[]            # acc_rate / ref_melt_rate (for ΔFAC model)
    rain_ratio_vec = Float64[]           # rain_rate / ref_melt_rate (for ΔFAC model)
    ref_fac_vec = Float64[]
    ref_melt_rate_vec = Float64[]
    ref_fac_rate_vec = Float64[]
    ref_refreeze_rate_vec = Float64[]    # Use reference rate
    ref_ec_rate_vec = Float64[]          # Use reference rate
    lat_vec = Float64[]
    lon_vec = Float64[]
    elev_vec = Float64[]
    temp_anom_vec = Float64[]
    pscale_vec = Float64[]

    for i in 1:n_total
        loc = loc_keys[i]
        if !haskey(ref_lookup, loc)
            continue
        end

        ref_state = ref_lookup[loc]

        # Skip if this IS the reference (no perturbation to learn from)
        if i == ref_state.idx
            continue
        end

        # Skip if reference melt is too small
        if abs(ref_state.melt_rate) < 1e-6
            continue
        end

        # Skip if current sample has NaN
        if isnan(mean_fac_absolute[i]) || isnan(mean_melt_rate[i]) || isnan(mean_acc_rate[i])
            continue
        end

        # Compute RATIOS for ΔFAC model (better predictors than deltas)
        ref_melt_safe = max(abs(ref_state.melt_rate), 1e-10)
        push!(melt_ratio_vec, mean_melt_rate[i] / ref_state.melt_rate)
        push!(acc_ratio_vec, mean_acc_rate[i] / ref_melt_safe)
        push!(rain_ratio_vec, mean_rain_rate[i] / ref_melt_safe)

        # ADDED: Δmelt_rate captures temperature/energy change (drives compaction)
        push!(delta_melt_rate_vec, mean_melt_rate[i] - ref_state.melt_rate)

        # Compute TARGETS
        push!(delta_fac_rate_vec, mean_fac_rate[i] - ref_state.fac_rate)  # TARGET 1: ΔFAC rate
        push!(delta_refreeze_rate_vec, mean_refreeze_rate[i] - ref_state.refreeze_rate)  # TARGET 2: Δrefreeze rate

        # Compute DELTAS for refreeze model (uses different features)
        push!(delta_acc_rate_vec, mean_acc_rate[i] - ref_state.acc_rate)
        push!(delta_rain_rate_vec, mean_rain_rate[i] - ref_state.rain_rate)

        # Reference state features (all known at mscale=1)
        push!(ref_fac_vec, ref_state.fac)
        push!(ref_melt_rate_vec, ref_state.melt_rate)
        push!(ref_fac_rate_vec, ref_state.fac_rate)
        push!(ref_refreeze_rate_vec, ref_state.refreeze_rate)  # Use reference, not current
        push!(ref_ec_rate_vec, ref_state.ec_rate)              # Use reference, not current

        # Location and climate
        push!(lat_vec, abs(gembX["latitude"][i]))
        push!(lon_vec, gembX["longitude"][i])
        push!(elev_vec, gembX["height"][i])
        push!(temp_anom_vec, gembX["elevation_delta"][i] * (-6.5 / 1000))
        push!(pscale_vec, gembX["precipitation_scale"][i])
    end

    println("  Created $(length(delta_fac_rate_vec)) training samples")

    # Assemble feature matrices (different features for each model)
    feature_names_fac = [
        "Δmelt_rate",         # ADDED: captures temperature/energy change (drives compaction)
        "melt_ratio",         # melt_rate / ref_melt_rate (= mscale)
        "acc_ratio",
        "rain_ratio",
        "ref_fac",
        "ref_melt_rate",
        "ref_fac_rate",
        "ref_refreeze_rate",
        "ref_ec_rate",
        "|latitude|",
        "elevation",
        "temp_anomaly",
        "pscale"
    ]

    feature_names_refreeze = [
        "Δacc_rate",          # PRIMARY for refreeze: delta features work well
        "Δrain_rate",
        "ref_fac",
        "ref_melt_rate",
        "ref_fac_rate",
        "ref_refreeze_rate",
        "ref_ec_rate",
        "|latitude|",
        "elevation",
        "temp_anomaly",
        "pscale"
    ]

    X_fac_all = hcat(
        delta_melt_rate_vec,  # ADDED: temperature/energy change
        melt_ratio_vec,
        acc_ratio_vec,
        rain_ratio_vec,
        ref_fac_vec,
        ref_melt_rate_vec,
        ref_fac_rate_vec,
        ref_refreeze_rate_vec,
        ref_ec_rate_vec,
        lat_vec,
        elev_vec,
        temp_anom_vec,
        pscale_vec
    )

    X_refreeze_all = hcat(
        delta_acc_rate_vec,
        delta_rain_rate_vec,
        ref_fac_vec,
        ref_melt_rate_vec,
        ref_fac_rate_vec,
        ref_refreeze_rate_vec,
        ref_ec_rate_vec,
        lat_vec,
        elev_vec,
        temp_anom_vec,
        pscale_vec
    )
    y_fac_all = delta_fac_rate_vec         # TARGET 1: ΔFAC rate
    y_refreeze_all = delta_refreeze_rate_vec  # TARGET 2: Δrefreeze rate

    lat_all = lat_vec
    lon_all = lon_vec
    elev_delta_all = temp_anom_vec ./ (-6.5 / 1000)
    pscale_all = pscale_vec

    # Filter NaN
    valid_fac = vec(.!(any(isnan.(X_fac_all), dims=2)) .& .!isnan.(y_fac_all))
    valid_refreeze = vec(.!(any(isnan.(X_refreeze_all), dims=2)) .& .!isnan.(y_refreeze_all))

    if sum(.!valid_fac) > 0
        println("  Filtered $(sum(.!valid_fac)) NaN samples from FAC dataset")
        X_fac_all = X_fac_all[valid_fac, :]
        y_fac_all = y_fac_all[valid_fac]
    end

    if sum(.!valid_refreeze) > 0
        println("  Filtered $(sum(.!valid_refreeze)) NaN samples from refreeze dataset")
        X_refreeze_all = X_refreeze_all[valid_refreeze, :]
        y_refreeze_all = y_refreeze_all[valid_refreeze]
    end

    # Use FAC valid indices for location metadata (both should be same)
    lat_all = lat_all[valid_fac]
    lon_all = lon_all[valid_fac]
    elev_delta_all = elev_delta_all[valid_fac]
    pscale_all = pscale_all[valid_fac]

    # Save cache
    println("Saving training cache to: $cache_file")
    jldsave(cache_file;
        X_fac=X_fac_all, X_refreeze=X_refreeze_all,
        y_fac=y_fac_all, y_refreeze=y_refreeze_all,
        feature_names_fac=feature_names_fac, feature_names_refreeze=feature_names_refreeze,
        lat=lat_all, lon=lon_all, elev_delta=elev_delta_all, pscale=pscale_all
    )
    println("  Saved $(size(X_fac_all, 1)) samples")
end;

# =============================================================================
# SUBSET DATA (optional - for faster iteration)
# =============================================================================
if use_subset
    loc_keys = [(lat_all[i], lon_all[i]) for i in 1:length(lat_all)]
    unique_locs = unique(loc_keys)
    println("\nTotal unique locations: $(length(unique_locs))")

    Random.seed!(42)
    n_locs = min(n_locations_subset, length(unique_locs))
    subset_locs = Set(unique_locs[randperm(length(unique_locs))[1:n_locs]])

    subset_mask = [loc_keys[i] in subset_locs for i in 1:length(loc_keys)]
    X_fac_all = X_fac_all[subset_mask, :]
    X_refreeze_all = X_refreeze_all[subset_mask, :]
    y_fac_all = y_fac_all[subset_mask]
    y_refreeze_all = y_refreeze_all[subset_mask]
    lat_all = lat_all[subset_mask]
    lon_all = lon_all[subset_mask]
    elev_delta_all = elev_delta_all[subset_mask]
    pscale_all = pscale_all[subset_mask]

    println("Subset: $n_locs locations → $(size(X_fac_all, 1)) samples")
end;

# =============================================================================
# TRAIN / TEST / VALIDATION SPLIT (stratified by location)
# =============================================================================
begin
    Random.seed!(42)

    loc_keys = [(lat_all[i], lon_all[i]) for i in 1:length(lat_all)]
    unique_locs = unique(loc_keys)
    n_locs = length(unique_locs)

    loc_perm = randperm(n_locs)
    n_train_locs = round(Int, 0.7 * n_locs)
    n_val_locs = round(Int, 0.15 * n_locs)

    train_locs = Set(unique_locs[loc_perm[1:n_train_locs]])
    val_locs = Set(unique_locs[loc_perm[n_train_locs+1:n_train_locs+n_val_locs]])
    test_locs = Set(unique_locs[loc_perm[n_train_locs+n_val_locs+1:end]])

    train_idx = findall(i -> loc_keys[i] in train_locs, 1:length(loc_keys))
    val_idx = findall(i -> loc_keys[i] in val_locs, 1:length(loc_keys))
    test_idx = findall(i -> loc_keys[i] in test_locs, 1:length(loc_keys))

    X_fac_train = X_fac_all[train_idx, :]
    X_fac_val = X_fac_all[val_idx, :]
    X_fac_test = X_fac_all[test_idx, :]

    X_refreeze_train = X_refreeze_all[train_idx, :]
    X_refreeze_val = X_refreeze_all[val_idx, :]
    X_refreeze_test = X_refreeze_all[test_idx, :]

    y_fac_train = y_fac_all[train_idx]
    y_fac_val = y_fac_all[val_idx]
    y_fac_test = y_fac_all[test_idx]

    y_refreeze_train = y_refreeze_all[train_idx]
    y_refreeze_val = y_refreeze_all[val_idx]
    y_refreeze_test = y_refreeze_all[test_idx]

    n_train = length(train_idx)
    n_val = length(val_idx)
    n_test = length(test_idx)

    println("\n" * "=" ^ 70)
    println("DATASET SUMMARY (location-stratified split)")
    println("=" ^ 70)
    println("Locations: $n_train_locs train, $n_val_locs val, $(length(test_locs)) test")
    println("Total samples: $(n_train + n_val + n_test)")
    println("Train: $n_train ($(round(100*n_train/(n_train+n_val+n_test), digits=1))%) | Val: $n_val ($(round(100*n_val/(n_train+n_val+n_test), digits=1))%) | Test: $n_test ($(round(100*n_test/(n_train+n_val+n_test), digits=1))%)")

    println("\nFAC model features (training set):")
    for (i, name) in enumerate(feature_names_fac)
        μ = round(mean(X_fac_train[:, i]), sigdigits=3)
        σ = round(std(X_fac_train[:, i]), sigdigits=3)
        println("  $(rpad(name, 20)): μ=$μ, σ=$σ")
    end
    println("  $(rpad("ΔFAC rate (target)", 20)): μ=$(round(mean(y_fac_train), sigdigits=3)), σ=$(round(std(y_fac_train), sigdigits=3))")

    println("\nRefreeze model features (training set):")
    for (i, name) in enumerate(feature_names_refreeze)
        μ = round(mean(X_refreeze_train[:, i]), sigdigits=3)
        σ = round(std(X_refreeze_train[:, i]), sigdigits=3)
        println("  $(rpad(name, 20)): μ=$μ, σ=$σ")
    end
    println("  $(rpad("Δrefreeze rate (target)", 20)): μ=$(round(mean(y_refreeze_train), sigdigits=3)), σ=$(round(std(y_refreeze_train), sigdigits=3))")
end;

# =============================================================================
# EVALUATION UTILITIES
# =============================================================================
function metrics(y_true, y_pred)
    res = y_true .- y_pred
    ss_res = sum(res .^ 2)
    ss_tot = sum((y_true .- mean(y_true)) .^ 2)
    r2 = 1.0 - ss_res / ss_tot
    rmse = sqrt(mean(res .^ 2))
    mae = mean(abs.(res))
    return (; r2=round(r2, digits=5), rmse=round(rmse, sigdigits=4), mae=round(mae, sigdigits=4))
end;

# =============================================================================
# RANDOM FOREST MODELS (hyperparameter search for both targets)
# =============================================================================
begin
    println("\n" * "=" ^ 70)
    println("RANDOM FOREST MODEL 1: ΔFAC rate (hyperparameter grid search)")
    println("=" ^ 70)

    best_rf_fac = nothing
    best_rf_fac_val_r2 = -Inf
    best_rf_fac_params = (n_trees=0, max_depth=0, n_sub=-1)

    configs = [(nt, md, ns) for nt in [100, 300, 500], md in [6, 10, 14, 18, 22], ns in [-1, 4]]

    println("  Testing $(length(configs)) configurations for ΔFAC rate...")
    for (nt, md, ns) in configs
        rf = build_forest(y_fac_train, X_fac_train, ns, nt, 0.7, md)
        val_r2 = metrics(y_fac_val, apply_forest(rf, X_fac_val)).r2
        if val_r2 > best_rf_fac_val_r2
            global best_rf_fac_val_r2 = val_r2
            global best_rf_fac_params = (n_trees=nt, max_depth=md, n_sub=ns)
            global best_rf_fac = rf
        end
    end

    y_fac_pred_test = apply_forest(best_rf_fac, X_fac_test)
    m_fac_train = metrics(y_fac_train, apply_forest(best_rf_fac, X_fac_train))
    m_fac_val = metrics(y_fac_val, apply_forest(best_rf_fac, X_fac_val))
    m_fac_test = metrics(y_fac_test, y_fac_pred_test)

    println("\n  Best params: n_trees=$(best_rf_fac_params.n_trees), max_depth=$(best_rf_fac_params.max_depth), n_sub=$(best_rf_fac_params.n_sub)")
    println("  Train R²=$(m_fac_train.r2) | Val R²=$(m_fac_val.r2) | Test R²=$(m_fac_test.r2)")
    println("  Test RMSE=$(m_fac_test.rmse), MAE=$(m_fac_test.mae)")

    if m_fac_test.r2 >= 0.9
        println("  ✓ SUCCESS: R² on ΔFAC rate ≥ 0.9 target")
    else
        println("  ✗ WARNING: R² on ΔFAC rate < 0.9 target")
    end

    gap_fac = m_fac_train.r2 - m_fac_test.r2
    if gap_fac > 0.1
        println("  ⚠ Potential overfitting: train-test R² gap = $(round(gap_fac, digits=3))")
    else
        println("  ✓ Minimal overfitting: gap = $(round(gap_fac, digits=3))")
    end

    # =========================================================================
    # MODEL 2: Δrefreeze rate
    # =========================================================================
    println("\n" * "=" ^ 70)
    println("RANDOM FOREST MODEL 2: Δrefreeze rate (hyperparameter grid search)")
    println("=" ^ 70)

    best_rf_refreeze = nothing
    best_rf_refreeze_val_r2 = -Inf
    best_rf_refreeze_params = (n_trees=0, max_depth=0, n_sub=-1)

    println("  Testing $(length(configs)) configurations for Δrefreeze rate...")
    for (nt, md, ns) in configs
        rf = build_forest(y_refreeze_train, X_refreeze_train, ns, nt, 0.7, md)
        val_r2 = metrics(y_refreeze_val, apply_forest(rf, X_refreeze_val)).r2
        if val_r2 > best_rf_refreeze_val_r2
            global best_rf_refreeze_val_r2 = val_r2
            global best_rf_refreeze_params = (n_trees=nt, max_depth=md, n_sub=ns)
            global best_rf_refreeze = rf
        end
    end

    y_refreeze_pred_test = apply_forest(best_rf_refreeze, X_refreeze_test)
    m_refreeze_train = metrics(y_refreeze_train, apply_forest(best_rf_refreeze, X_refreeze_train))
    m_refreeze_val = metrics(y_refreeze_val, apply_forest(best_rf_refreeze, X_refreeze_val))
    m_refreeze_test = metrics(y_refreeze_test, y_refreeze_pred_test)

    println("\n  Best params: n_trees=$(best_rf_refreeze_params.n_trees), max_depth=$(best_rf_refreeze_params.max_depth), n_sub=$(best_rf_refreeze_params.n_sub)")
    println("  Train R²=$(m_refreeze_train.r2) | Val R²=$(m_refreeze_val.r2) | Test R²=$(m_refreeze_test.r2)")
    println("  Test RMSE=$(m_refreeze_test.rmse), MAE=$(m_refreeze_test.mae)")

    if m_refreeze_test.r2 >= 0.9
        println("  ✓ SUCCESS: R² on Δrefreeze rate ≥ 0.9 target")
    else
        println("  ⚠ Note: R² on Δrefreeze rate < 0.9")
    end

    gap_refreeze = m_refreeze_train.r2 - m_refreeze_test.r2
    if gap_refreeze > 0.1
        println("  ⚠ Potential overfitting: train-test R² gap = $(round(gap_refreeze, digits=3))")
    else
        println("  ✓ Minimal overfitting: gap = $(round(gap_refreeze, digits=3))")
    end
end;

# =============================================================================
# FEATURE IMPORTANCE (permutation-based)
# =============================================================================
begin
    println("\n" * "=" ^ 70)
    println("FEATURE IMPORTANCE (permutation-based)")
    println("=" ^ 70)

    function compute_importance(model, X_test0, y_test0; n_repeats=5)
        base_r2 = metrics(y_test0, apply_forest(model, X_test0)).r2
        imp = zeros(size(X_test0, 2))
        for j in 1:size(X_test0, 2)
            scores = Float64[]
            for _ in 1:n_repeats
                Xp = copy(X_test0)
                Xp[:, j] = shuffle(Xp[:, j])
                push!(scores, metrics(y_test0, apply_forest(model, Xp)).r2)
            end
            imp[j] = base_r2 - mean(scores)
        end
        return imp
    end

    # Feature importance for ΔFAC rate model
    println("\n  MODEL 1: ΔFAC rate")
    rf_fac_importance = compute_importance(best_rf_fac, X_fac_test, y_fac_test)
    sorted_idx_fac = sortperm(rf_fac_importance, rev=true)

    for i in sorted_idx_fac
        bar_len = max(1, round(Int, rf_fac_importance[i] / maximum(rf_fac_importance) * 40))
        bar = repeat("█", bar_len)
        println("  $(rpad(feature_names_fac[i], 20)) $(rpad(round(rf_fac_importance[i], digits=5), 8)) $bar")
    end

    # Feature importance for Δrefreeze rate model
    println("\n  MODEL 2: Δrefreeze rate")
    rf_refreeze_importance = compute_importance(best_rf_refreeze, X_refreeze_test, y_refreeze_test)
    sorted_idx_refreeze = sortperm(rf_refreeze_importance, rev=true)

    for i in sorted_idx_refreeze
        bar_len = max(1, round(Int, rf_refreeze_importance[i] / maximum(rf_refreeze_importance) * 40))
        bar = repeat("█", bar_len)
        println("  $(rpad(feature_names_refreeze[i], 20)) $(rpad(round(rf_refreeze_importance[i], digits=5), 8)) $bar")
    end
end;

# =============================================================================
# SAVE MODELS FOR PRODUCTION USE
# =============================================================================
begin
    model_file = joinpath(dirname(cache_file), "fac_surrogate_random_forest_delta.jld2")
    println("\n" * "=" ^ 70)
    println("SAVING PRODUCTION MODELS")
    println("=" ^ 70)

    jldsave(model_file;
        # Model 1: ΔFAC rate
        model_fac=best_rf_fac,
        fac_n_trees=best_rf_fac_params.n_trees,
        fac_max_depth=best_rf_fac_params.max_depth,
        fac_n_subfeatures=best_rf_fac_params.n_sub,
        fac_train_r2=m_fac_train.r2,
        fac_val_r2=m_fac_val.r2,
        fac_test_r2=m_fac_test.r2,
        fac_test_rmse=m_fac_test.rmse,
        fac_test_mae=m_fac_test.mae,

        # Model 2: Δrefreeze rate
        model_refreeze=best_rf_refreeze,
        refreeze_n_trees=best_rf_refreeze_params.n_trees,
        refreeze_max_depth=best_rf_refreeze_params.max_depth,
        refreeze_n_subfeatures=best_rf_refreeze_params.n_sub,
        refreeze_train_r2=m_refreeze_train.r2,
        refreeze_val_r2=m_refreeze_val.r2,
        refreeze_test_r2=m_refreeze_test.r2,
        refreeze_test_rmse=m_refreeze_test.rmse,
        refreeze_test_mae=m_refreeze_test.mae,

        # Shared metadata
        feature_names_fac=feature_names_fac,
        feature_names_refreeze=feature_names_refreeze,
        description="Two Random Forest models: (1) ΔFAC rate (ratio features) and (2) Δrefreeze rate (delta features). All features known at mscale=1."
    )

    println("  Saved to: $model_file")
    println("\n  MODEL 1 - ΔFAC rate:")
    println("    Test R² = $(round(m_fac_test.r2, digits=5)) ← PRIMARY METRIC")
    println("    Test RMSE = $(round(m_fac_test.rmse, digits=5))")
    println("\n  MODEL 2 - Δrefreeze rate:")
    println("    Test R² = $(round(m_refreeze_test.r2, digits=5))")
    println("    Test RMSE = $(round(m_refreeze_test.rmse, digits=5))")
    println("\n  Feature sets:")
    println("    MODEL 1 (ΔFAC): melt_ratio, acc_ratio, rain_ratio + reference rates")
    println("    MODEL 2 (Δrefreeze): Δacc_rate, Δrain_rate + reference rates")
    println("    Both: All features known at mscale=1 (production time)")
end;

println("\n" * "=" ^ 70)
println("TRAINING COMPLETE")
println("=" ^ 70)
println("Next steps:")
println("  1. Run src/fac_surrogate_model_export.jl to export final model")
println("  2. Run src/fac_surrogate_validate_change.jl to validate R² on ΔFAC")
println("  3. Run src/fac_surrogate_validation_plots.jl to generate diagnostic plots")
