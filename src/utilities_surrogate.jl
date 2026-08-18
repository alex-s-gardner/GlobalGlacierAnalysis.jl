# Surrogate models for GEMB variables
# Provides trained models to predict variable responses under perturbations

using JLD2
using DecisionTree

"""
    FAC_SURROGATE

Singleton container for FAC surrogate model, loaded once at module initialization.

The model predicts absolute FAC given melt scaling and reference climatology:
- Input: 12 features including melt_ratio, ref_fac, acc_ratio, etc.
- Output: Absolute FAC (meters of air)
- Performance: Test R² = 0.96498 on ΔFAC (change from reference)

Features (must be in this exact order):
1. melt_ratio: melt_scaled / melt_reference (dimensionless)
2. acc_ratio: accumulation_rate / ref_melt_rate (dimensionless)
3. refreeze_ratio: refreeze_rate / ref_melt_rate (dimensionless)
4. rain_ratio: rain_rate / ref_melt_rate (dimensionless)
5. ec_ratio: evap_cond_rate / ref_melt_rate (dimensionless)
6. ref_fac: FAC at reference conditions (mscale=1) [m]
7. ref_melt_rate: Melt rate at reference [m/day]
8. ref_fac_rate: FAC rate at reference [m/day]
9. |latitude|: Absolute latitude [degrees]
10. elevation: Surface elevation [m]
11. temp_anomaly: Temperature perturbation [K]
12. pscale: Precipitation scaling factor [dimensionless]

Returns:
- Nothing if model file not found
- NamedTuple with model and metadata if loaded successfully
"""
const FAC_SURROGATE = let
    model_filename = "fac_surrogate_random_forest.jld2"

    # Search for model in multiple possible locations
    search_paths = [
        joinpath(dirname(@__FILE__), "..", "data", "gemb", "raw", model_filename),
        joinpath("/mnt/bylot-r3/data/gemb/raw", model_filename),
        joinpath(dirname(@__FILE__), "..", "..", "data", "gemb", "raw", model_filename)
    ]

    model_path = findfirst(isfile, search_paths)

    if isnothing(model_path)
        @warn """
        FAC surrogate model not found. Searched:
        $(join(search_paths, "\n  "))
        FAC will not be scaled with melt in :mscale method.
        """
        nothing
    else
        try
            model_data = load(search_paths[model_path])
            @info "Loaded FAC surrogate model from $(search_paths[model_path])"
            @info "  Model: $(model_data["n_trees"]) trees, max_depth=$(model_data["max_depth"])"
            @info "  Performance: Test R²=$(round(model_data["test_r2"], digits=5))"

            (
                model = model_data["model"],
                feature_names = model_data["feature_names"],
                n_trees = model_data["n_trees"],
                max_depth = model_data["max_depth"],
                n_subfeatures = model_data["n_subfeatures"],
                test_r2 = model_data["test_r2"],
                description = model_data["description"]
            )
        catch e
            @warn "Failed to load FAC surrogate model" exception=(e, catch_backtrace())
            nothing
        end
    end
end


"""
    predict_fac_scaled(mscale, ref_state, rates, location, climate_params)

Predict FAC at scaled melt conditions using the trained surrogate model.

# Arguments
- `mscale::Float64`: Melt scaling factor (e.g., 2.0 doubles melt)
- `ref_state::NamedTuple`: Reference state at mscale=1
  - `fac::Float64`: Reference FAC [m]
  - `melt_rate::Float64`: Reference melt rate [m/day]
  - `fac_rate::Float64`: Reference FAC rate [m/day]
- `rates::NamedTuple`: Scaled flux rates at current perturbation
  - `acc_rate::Float64`: Accumulation rate [m/day]
  - `refreeze_rate::Float64`: Refreezing rate [m/day]
  - `rain_rate::Float64`: Rain rate [m/day]
  - `ec_rate::Float64`: Evaporation/condensation rate [m/day]
- `location::NamedTuple`: Location parameters
  - `latitude::Float64`: Latitude [degrees, -90 to 90]
  - `elevation::Float64`: Surface elevation [m]
- `climate_params::NamedTuple`: Climate perturbation parameters
  - `temp_anomaly::Float64`: Temperature anomaly [K]
  - `pscale::Float64`: Precipitation scaling factor

# Returns
- `Float64`: Predicted FAC at scaled conditions [m], always ≥ 0

# Example
```julia
fac_predicted = predict_fac_scaled(
    2.0,  # mscale: double melt
    (fac=1.5, melt_rate=0.2, fac_rate=-0.0001),  # reference state
    (acc_rate=0.08, refreeze_rate=0.015, rain_rate=0.04, ec_rate=0.001),  # rates
    (latitude=65.0, elevation=1500.0),  # location
    (temp_anomaly=0.0, pscale=1.0)  # climate
)
```
"""
function predict_fac_scaled(mscale::Real,
                            ref_state::NamedTuple,
                            rates::NamedTuple,
                            location::NamedTuple,
                            climate_params::NamedTuple)

    # Check if model is loaded
    if isnothing(FAC_SURROGATE)
        @warn "FAC surrogate model not loaded - returning reference FAC"
        return max(ref_state.fac, 0.0)
    end

    # Safety: avoid division by very small reference melt
    ref_melt_safe = max(abs(ref_state.melt_rate), 1e-10)

    # Build feature vector in exact training order
    features = Float64[
        mscale,                                    # 1. melt_ratio
        rates.acc_rate / ref_melt_safe,           # 2. acc_ratio
        rates.refreeze_rate / ref_melt_safe,      # 3. refreeze_ratio
        rates.rain_rate / ref_melt_safe,          # 4. rain_ratio
        rates.ec_rate / ref_melt_safe,            # 5. ec_ratio
        ref_state.fac,                             # 6. ref_fac
        ref_state.melt_rate,                       # 7. ref_melt_rate
        ref_state.fac_rate,                        # 8. ref_fac_rate
        abs(location.latitude),                    # 9. |latitude|
        location.elevation,                        # 10. elevation
        climate_params.temp_anomaly,               # 11. temp_anomaly
        climate_params.pscale                      # 12. pscale
    ]

    # Predict (reshape to 2D: 1 sample × 12 features)
    fac_pred = apply_forest(FAC_SURROGATE.model, reshape(features, 1, :))[1]

    # Enforce physical constraint: FAC ≥ 0
    return max(fac_pred, 0.0)
end


"""
    predict_fac_scaled_batch(mscale, ref_states, rates_batch, locations, climate_params)

Vectorized batch prediction for multiple (date, height) points simultaneously.

More efficient than calling `predict_fac_scaled` in a loop.

# Arguments
All arguments are vectors/matrices where element `i` corresponds to point `i`:
- `mscale::Float64`: Scalar melt scaling factor (same for all points)
- `ref_states`: Vector of NamedTuples with (fac, melt_rate, fac_rate)
- `rates_batch`: Vector of NamedTuples with (acc_rate, refreeze_rate, rain_rate, ec_rate)
- `locations`: Vector of NamedTuples with (latitude, elevation)
- `climate_params::NamedTuple`: Scalar params (temp_anomaly, pscale) same for all

# Returns
- `Vector{Float64}`: Predicted FAC for each point [m], all ≥ 0
"""
function predict_fac_scaled_batch(mscale::Real,
                                  ref_states::Vector,
                                  rates_batch::Vector,
                                  locations::Vector,
                                  climate_params::NamedTuple)

    n_points = length(ref_states)
    @assert length(rates_batch) == n_points
    @assert length(locations) == n_points

    # Check if model is loaded
    if isnothing(FAC_SURROGATE)
        @warn "FAC surrogate model not loaded - returning reference FAC values"
        return max.([ref_states[i].fac for i in 1:n_points], 0.0)
    end

    # Build feature matrix: n_points × 12 features
    X = zeros(n_points, 12)

    for i in 1:n_points
        ref_melt_safe = max(abs(ref_states[i].melt_rate), 1e-10)

        X[i, 1] = mscale
        X[i, 2] = rates_batch[i].acc_rate / ref_melt_safe
        X[i, 3] = rates_batch[i].refreeze_rate / ref_melt_safe
        X[i, 4] = rates_batch[i].rain_rate / ref_melt_safe
        X[i, 5] = rates_batch[i].ec_rate / ref_melt_safe
        X[i, 6] = ref_states[i].fac
        X[i, 7] = ref_states[i].melt_rate
        X[i, 8] = ref_states[i].fac_rate
        X[i, 9] = abs(locations[i].latitude)
        X[i, 10] = locations[i].elevation
        X[i, 11] = climate_params.temp_anomaly
        X[i, 12] = climate_params.pscale
    end

    # Predict all points at once (vectorized)
    fac_predictions = apply_forest(FAC_SURROGATE.model, X)

    # Enforce physical constraint
    return max.(fac_predictions, 0.0)
end


"""
    is_fac_surrogate_available()

Check if FAC surrogate model is loaded and ready to use.

# Returns
- `Bool`: true if model is loaded, false otherwise
"""
is_fac_surrogate_available() = !isnothing(FAC_SURROGATE)


"""
    get_fac_surrogate_info()

Get information about the loaded FAC surrogate model.

# Returns
- `NamedTuple` with model metadata or `nothing` if not loaded
"""
function get_fac_surrogate_info()
    if isnothing(FAC_SURROGATE)
        return nothing
    end

    return (
        n_trees = FAC_SURROGATE.n_trees,
        max_depth = FAC_SURROGATE.max_depth,
        n_features = length(FAC_SURROGATE.feature_names),
        feature_names = FAC_SURROGATE.feature_names,
        test_r2 = FAC_SURROGATE.test_r2,
        description = FAC_SURROGATE.description
    )
end

# Export public API
export predict_fac_scaled, predict_fac_scaled_batch
export is_fac_surrogate_available, get_fac_surrogate_info
