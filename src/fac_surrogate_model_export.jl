# Export the trained Random Forest FAC surrogate model for integration into utilities_gemb.jl
# This script loads the cached training data and retrains the best model, then exports it

begin
    import GlobalGlacierAnalysis as GGA
    using JLD2
    using DecisionTree
    using Statistics
    using Random

    # Load cached training data
    cache_file = "/mnt/bylot-r3/data/gemb/raw/fac_surrogate_training_data_ratio.jld2"
    println("Loading training data from: $cache_file")
    training_data = load(cache_file)

    X_all = training_data["X"]
    y_all = training_data["y"]
    feature_names = training_data["feature_names"]
    lat_all = training_data["lat"]
    lon_all = training_data["lon"]

    println("  Loaded $(size(X_all, 1)) samples with $(size(X_all, 2)) features")
    println("  Feature names: ", feature_names)
end

# Subset to 1000 locations for faster model export (same as training)
use_subset = true
n_locations_subset = 1000

if use_subset
    loc_keys = [(lat_all[i], lon_all[i]) for i in 1:length(lat_all)]
    unique_locs = unique(loc_keys)
    Random.seed!(42)
    selected_locs = Set(shuffle(unique_locs)[1:min(n_locations_subset, length(unique_locs))])
    subset_mask = [loc_keys[i] in selected_locs for i in 1:length(lat_all)]

    X_all = X_all[subset_mask, :]
    y_all = y_all[subset_mask]
    lat_all = lat_all[subset_mask]
    lon_all = lon_all[subset_mask]

    println("\nSubset to $n_locations_subset locations: $(size(X_all, 1)) samples")
end

# Location-stratified train/val/test split (70/15/15)
loc_keys = [(lat_all[i], lon_all[i]) for i in 1:length(lat_all)]
unique_locs = unique(loc_keys)
n_locs = length(unique_locs)

Random.seed!(42)
shuffled_locs = shuffle(unique_locs)
n_train_locs = round(Int, 0.70 * n_locs)
n_val_locs = round(Int, 0.15 * n_locs)

train_locs = Set(shuffled_locs[1:n_train_locs])
val_locs = Set(shuffled_locs[n_train_locs+1:n_train_locs+n_val_locs])
test_locs = Set(shuffled_locs[n_train_locs+n_val_locs+1:end])

train_mask = [loc_keys[i] in train_locs for i in 1:length(lat_all)]
val_mask = [loc_keys[i] in val_locs for i in 1:length(lat_all)]
test_mask = [loc_keys[i] in test_locs for i in 1:length(lat_all)]

X_train = X_all[train_mask, :]
y_train = y_all[train_mask]
X_val = X_all[val_mask, :]
y_val = y_all[val_mask]
X_test = X_all[test_mask, :]
y_test = y_all[test_mask]

println("\nDataset splits:")
println("  Train: $(size(X_train, 1)) samples ($(length(train_locs)) locations)")
println("  Val: $(size(X_val, 1)) samples ($(length(val_locs)) locations)")
println("  Test: $(size(X_test, 1)) samples ($(length(test_locs)) locations)")

# Train best Random Forest model
println("\n" * "="^70)
println("TRAINING FINAL RANDOM FOREST MODEL")
println("="^70)

# Best hyperparameters from grid search
n_trees = 500
max_depth = 18
n_subfeatures = 4

println("  Hyperparameters: n_trees=$n_trees, max_depth=$max_depth, n_subfeatures=$n_subfeatures")
println("  Training...")

Random.seed!(42)
# build_forest(labels, features, n_subfeatures, n_trees, partial_sampling, max_depth)
model = build_forest(y_train, X_train, n_subfeatures, n_trees, 1.0, max_depth)

# Evaluate
y_train_pred = apply_forest(model, X_train)
y_val_pred = apply_forest(model, X_val)
y_test_pred = apply_forest(model, X_test)

function r2_score(y_true, y_pred)
    ss_res = sum((y_true .- y_pred).^2)
    ss_tot = sum((y_true .- mean(y_true)).^2)
    return 1 - ss_res / ss_tot
end

train_r2 = r2_score(y_train, y_train_pred)
val_r2 = r2_score(y_val, y_val_pred)
test_r2 = r2_score(y_test, y_test_pred)

println("  Train R² = $(round(train_r2, digits=5))")
println("  Val R² = $(round(val_r2, digits=5))")
println("  Test R² = $(round(test_r2, digits=5))")

# Export model
model_file = "/mnt/bylot-r3/data/gemb/raw/fac_surrogate_random_forest.jld2"
println("\n" * "="^70)
println("EXPORTING MODEL")
println("="^70)

jldsave(model_file;
    model=model,
    feature_names=feature_names,
    n_trees=n_trees,
    max_depth=max_depth,
    n_subfeatures=n_subfeatures,
    train_r2=train_r2,
    val_r2=val_r2,
    test_r2=test_r2,
    description="Random Forest surrogate for FAC scaling. Predicts absolute FAC given melt_ratio and reference state.",
    usage_notes="""
    To use this model in utilities_gemb.jl:

    1. Load model at module level:
       fac_surrogate = load("$model_file")
       fac_model = fac_surrogate["model"]
       fac_features = fac_surrogate["feature_names"]

    2. In process_gemb_geotiles, at line ~1902 where FAC is assigned:
       # Build feature matrix
       melt_ratio = mscale  # since melt is scaled by mscale
       acc_ratio = acc_rate_scaled / (abs(melt_ref) + 1e-10)
       refreeze_ratio = refreeze_rate_scaled / (abs(melt_ref) + 1e-10)
       rain_ratio = rain_rate_scaled / (abs(melt_ref) + 1e-10)
       ec_ratio = ec_rate_scaled / (abs(melt_ref) + 1e-10)

       # ref_fac, ref_melt_rate, ref_fac_rate are from mscale=1 case
       # latitude, elevation, temp_anomaly (from elevation_delta), pscale

       features = [melt_ratio, acc_ratio, refreeze_ratio, rain_ratio, ec_ratio,
                   ref_fac, ref_melt_rate, ref_fac_rate,
                   abs(latitude), elevation, temp_anomaly, pscale]

       # Predict FAC
       fac_scaled = apply_forest(fac_model, reshape(features, 1, :))
       v0[:, :] = max.(fac_scaled, 0.0)  # Enforce FAC >= 0
    """
)

println("  Saved model to: $model_file")
println("  Model type: Random Forest with $n_trees trees")
println("  Test performance: R² = $(round(test_r2, digits=5))")
println("\nModel export complete!")
