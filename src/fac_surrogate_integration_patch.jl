# FAC Surrogate Model Integration Patch
# This file contains the code changes needed in utilities_gemb.jl to integrate the FAC surrogate model

# ============================================================================
# STEP 1: Add at module level (near top of utilities_gemb.jl, after imports)
# ============================================================================

"""
Load FAC surrogate model for melt scaling predictions.
Model predicts absolute FAC given melt_ratio and reference climatology.
"""
const FAC_SURROGATE_MODEL = let
    model_path = joinpath(dirname(@__FILE__), "..", "..", "data", "gemb", "raw", "fac_surrogate_random_forest.jld2")
    if isfile(model_path)
        try
            model_data = JLD2.load(model_path)
            @info "Loaded FAC surrogate model: R²=$(round(model_data["test_r2"], digits=4))"
            model_data
        catch e
            @warn "Failed to load FAC surrogate model from $model_path" exception=e
            nothing
        end
    else
        @warn "FAC surrogate model not found at $model_path - FAC will not be scaled with melt"
        nothing
    end
end

# ============================================================================
# STEP 2: Replace lines 1894-1925 in utilities_gemb.jl
# ============================================================================

# Original code structure:
# for k in Symbol.(vars)
#     for pscale in dpscale
#         v = gemb1[k][date=daterange, pscale=At(pscale)]
#         v0 = zeros(ddate, dheight)
#         v0 = v0[date=daterange]
#         for mscale in dmscale
#             if elevation_classes_method == :mscale
#                 if k == :melt
#                     v0[:, :] = v[:, index_height] .* mscale
#                 else
#                     v0[:, :] = v[:, index_height]  # ← PROBLEM HERE for FAC
#                 end
#             end
#             # ... store results ...
#         end
#     end
# end

# NEW CODE (replace the above section):

for k in Symbol.(vars)
    for pscale in dpscale
        v = gemb1[k][date=daterange, pscale=At(pscale)]
        v0 = zeros(ddate, dheight)
        v0 = v0[date=daterange]

        # For FAC surrogate model: store reference values (mscale=1) first pass
        v_ref_dict = if elevation_classes_method == :mscale && k == :fac && !isnothing(FAC_SURROGATE_MODEL)
            Dict{Symbol, Array}()
        else
            nothing
        end

        # First pass: extract reference values if needed
        if !isnothing(v_ref_dict)
            # Get all variables at this pscale for reference case
            for k_ref in Symbol.(vars)
                v_ref_temp = gemb1[k_ref][date=daterange, pscale=At(pscale)]
                v_ref_dict[k_ref] = v_ref_temp[:, index_height].data
            end
        end

        for mscale in dmscale
            if elevation_classes_method == :Δelevation
                Δheight = round((mscale * (1000 / -6.5)) / dh1, digits=0) * dh1
                height_effective = dheight1.val .- Δheight
                index_height_effective = (height_effective .>= dheight.val[1]) .& (height_effective .<= dheight.val[end])
                v0[:, :] = v[:, index_height_effective]

            elseif elevation_classes_method == :mscale

                if k == :melt
                    # Scale melt by mscale factor
                    v0[:, :] = v[:, index_height] .* mscale

                elseif k == :fac && !isnothing(FAC_SURROGATE_MODEL) && mscale != 1.0
                    # Use surrogate model to predict FAC at scaled melt

                    # Extract features
                    fac_model = FAC_SURROGATE_MODEL["model"]

                    # Reference values (mscale=1)
                    ref_fac = v_ref_dict[:fac]  # [n_dates, n_heights]
                    ref_melt = v_ref_dict[:melt]
                    ref_acc = v_ref_dict[:acc]
                    ref_refreeze = v_ref_dict[:refreeze]
                    ref_rain = v_ref_dict[:rain]
                    ref_ec = v_ref_dict[:ec]

                    # Compute reference FAC rate (dFAC/dt)
                    ref_fac_rate = zeros(size(ref_fac))
                    ref_fac_rate[2:end, :] = diff(ref_fac, dims=1) ./ 30.0  # per day
                    ref_fac_rate[1, :] = ref_fac_rate[2, :]  # extend to first timestep

                    # Compute reference melt rate
                    ref_melt_rate = zeros(size(ref_melt))
                    ref_melt_rate[2:end, :] = diff(ref_melt, dims=1) ./ 30.0
                    ref_melt_rate[1, :] = ref_melt_rate[2, :]

                    # Scaled rates (current perturbation)
                    acc_rate_scaled = zeros(size(ref_acc))
                    acc_rate_scaled[2:end, :] = diff(v_ref_dict[:acc], dims=1) ./ 30.0
                    acc_rate_scaled[1, :] = acc_rate_scaled[2, :]

                    refreeze_rate_scaled = zeros(size(ref_refreeze))
                    refreeze_rate_scaled[2:end, :] = diff(v_ref_dict[:refreeze], dims=1) ./ 30.0
                    refreeze_rate_scaled[1, :] = refreeze_rate_scaled[2, :]

                    rain_rate_scaled = zeros(size(ref_rain))
                    rain_rate_scaled[2:end, :] = diff(v_ref_dict[:rain], dims=1) ./ 30.0
                    rain_rate_scaled[1, :] = rain_rate_scaled[2, :]

                    ec_rate_scaled = zeros(size(ref_ec))
                    ec_rate_scaled[2:end, :] = diff(v_ref_dict[:ec], dims=1) ./ 30.0
                    ec_rate_scaled[1, :] = ec_rate_scaled[2, :]

                    # Build feature matrix for all (date, height) points
                    n_dates, n_heights = size(ref_fac)
                    n_points = n_dates * n_heights

                    # Feature order: ["melt_ratio", "acc_ratio", "refreeze_ratio", "rain_ratio", "ec_ratio",
                    #                 "ref_fac", "ref_melt_rate", "ref_fac_rate",
                    #                 "|latitude|", "elevation", "temp_anomaly", "pscale"]

                    X_features = zeros(n_points, 12)

                    # Get geotile-specific values (need to access from parent scope)
                    # These should be available in the function scope
                    lat_val = abs(geotile_row.lat_center)  # or appropriate latitude value
                    elev_vals = dheight.val  # elevation vector
                    temp_anomaly = 0.0  # Since we're using mscale method, no elevation delta
                    pscale_val = pscale.val[1]  # current pscale value

                    idx = 1
                    for iheight in 1:n_heights
                        for idate in 1:n_dates
                            # Safety: avoid division by very small ref_melt
                            ref_melt_safe = max(abs(ref_melt_rate[idate, iheight]), 1e-10)

                            X_features[idx, 1] = mscale  # melt_ratio
                            X_features[idx, 2] = acc_rate_scaled[idate, iheight] / ref_melt_safe  # acc_ratio
                            X_features[idx, 3] = refreeze_rate_scaled[idate, iheight] / ref_melt_safe  # refreeze_ratio
                            X_features[idx, 4] = rain_rate_scaled[idate, iheight] / ref_melt_safe  # rain_ratio
                            X_features[idx, 5] = ec_rate_scaled[idate, iheight] / ref_melt_safe  # ec_ratio
                            X_features[idx, 6] = ref_fac[idate, iheight]  # ref_fac
                            X_features[idx, 7] = ref_melt_rate[idate, iheight]  # ref_melt_rate
                            X_features[idx, 8] = ref_fac_rate[idate, iheight]  # ref_fac_rate
                            X_features[idx, 9] = lat_val  # |latitude|
                            X_features[idx, 10] = elev_vals[iheight]  # elevation
                            X_features[idx, 11] = temp_anomaly  # temp_anomaly (0 for mscale method)
                            X_features[idx, 12] = pscale_val  # pscale

                            idx += 1
                        end
                    end

                    # Predict FAC for all points at once (vectorized)
                    fac_predictions = DecisionTree.apply_forest(fac_model, X_features)

                    # Reshape back to (n_dates, n_heights)
                    fac_predicted_2d = reshape(fac_predictions, n_dates, n_heights)

                    # Enforce physical constraint: FAC >= 0
                    fac_predicted_2d = max.(fac_predicted_2d, 0.0)

                    # Assign to output
                    v0[:, :] = fac_predicted_2d

                else
                    # Default: use reference values (mscale=1) for all non-melt, non-FAC variables
                    v0[:, :] = v[:, index_height]
                end
            else
                error("unknown elevation_classes_method: $(elevation_classes_method)")
            end

            dv = @d v0 .* geotile_hyps_area_km2 ./ 1000
            dv[isnan.(dv)] .= 0
            dv = sum(dv; dims=:height)
            gemb_dv0[k][geotile=At(geotile_row.id), date=daterange, mscale=At(mscale), pscale=At(pscale)] = dropdims(dv, dims=:height)
        end
    end
end

# ============================================================================
# NOTES ON INTEGRATION
# ============================================================================

# IMPORTANT CONSIDERATIONS:
#
# 1. **Variable Access**: The code assumes access to:
#    - geotile_row.lat_center (or appropriate latitude field)
#    - dheight.val (elevation bins)
#    - pscale (current precipitation scale)
#    Need to verify these variable names match actual code
#
# 2. **Performance**: Prediction for n_dates × n_heights points per call
#    - For 100 dates × 50 heights = 5000 predictions per (pscale, mscale) combo
#    - Random Forest inference is fast but may add ~seconds per geotile
#    - Consider profiling and optimization if bottleneck
#
# 3. **Cumulative vs Rate Variables**:
#    - GEMB stores cumulative values for acc, melt, refreeze, rain, ec
#    - Need diff() to get rates
#    - FAC is a state variable, not cumulative
#    - Check if this is handled correctly in read_gemb_files
#
# 4. **Feature Order Critical**:
#    - MUST match training order exactly
#    - See FAC_SURROGATE_MODEL["feature_names"] to verify
#
# 5. **Error Handling**:
#    - If model fails to load, code falls back to reference values
#    - Consider adding try-catch around prediction for robustness
#
# 6. **Testing Strategy**:
#    - Test on single geotile first
#    - Compare FAC with/without surrogate (set FAC_SURROGATE_MODEL = nothing)
#    - Validate monotonicity: FAC should generally decrease with mscale
#    - Check for NaN/Inf in predictions
#    - Verify FAC >= 0 constraint enforced

# ============================================================================
# VALIDATION QUERY TO RUN AFTER INTEGRATION
# ============================================================================

# julia> using JLD2
# julia> using Statistics
# julia>
# julia> # Load output
# julia> gemb_result = load("path/to/gemb_output_with_surrogate.jld2")
# julia> fac_data = gemb_result["fac"]
# julia>
# julia> # Check for issues
# julia> println("Any NaN: ", any(isnan.(fac_data)))
# julia> println("Any negative: ", any(fac_data .< 0))
# julia> println("FAC range: [$(minimum(fac_data)), $(maximum(fac_data))]")
# julia>
# julia> # Check monotonicity vs mscale
# julia> fac_by_mscale = [mean(fac_data[mscale=At(m)]) for m in 0.5:0.5:2.0]
# julia> println("FAC vs mscale: ", fac_by_mscale)
# julia> println("Is decreasing? ", issorted(fac_by_mscale, rev=true))

