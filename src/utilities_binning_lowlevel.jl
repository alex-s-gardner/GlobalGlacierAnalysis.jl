# define models
begin
    # -------------------------------
    # Model definitions for binned data
    # -------------------------------

    # General polynomial + seasonal model to fit to all data binned by hypsometry.
    # Input:
    #   x: matrix with columns [time, elevation]
    #   p: parameter vector of length 9
    # Output:
    #   vector of modeled values
    # Parameters:
    #   p[1]: intercept
    #   p[2], p[3]: linear/quadratic dependence on time (x[:,1])
    #   p[4], p[5]: linear/quadratic dependence on elevation (x[:,2])
    #   p[6]: phase shift for annual sine terms
    #   p[7], p[8], p[9]: coefficients for annual sine (possibly elevation-modulated)
    model1(x, p) =
        p[1] .+                                       # Intercept
        p[2] .* x[:, 1] .+                            # Linear term in time
        p[3] .* x[:, 1] .^ 2 .+                       # Quadratic term in time
        p[4] .* x[:, 2] .+                            # Linear term in elevation
        p[5] .* x[:, 2] .^ 2 .+                       # Quadratic term in elevation
        p[7] .* sin.(2 .* pi .* (x[:, 1] .+ p[6])) .+ # Sine (annual, phase shifted)
        p[8] .* sin.(2 .* pi .* (x[:, 1] .+ p[6])) .* x[:, 2] .+           # Sine * elevation
        p[9] .* sin.(2 .* pi .* (x[:, 1] .+ p[6])) .* x[:, 2] .^ 2         # Sine * elevation^2

    # Initial parameters, lower and upper bounds for model1
    const p1 = zeros(9)
    const lb1 = [-10.0, -3.0, -2.0, -0.05, -0.0001, -1.0, -7.0, -0.05, -0.001]
    const ub1 = [+10.0, +3.0, +2.0, +0.05, +0.0001, +1.0, +7.0, +0.05, +0.001]

    # Seasonal-only model (annual sine with phase modulation and seasonal amplitude as quadratic in elevation)
    # p[1]: phase shift
    # p[2], p[3], p[4]: coefficients for amplitude (constant, linear, quadratic in elevation)
    model1_seasonal(x, p) =
        p[2] .* sin.(2 .* pi .* (x[:, 1] .+ p[1])) .+            # Sine term (base amplitude)
        p[3] .* sin.(2 .* pi .* (x[:, 1] .+ p[1])) .* x[:, 2] .+ # Sine * elevation
        p[4] .* sin.(2 .* pi .* (x[:, 1] .+ p[1])) .* x[:, 2] .^ 2 # Sine * elevation^2

    # Simple linear trend model: fit in elevation only
    # p[1]: intercept, p[2]: slope
    model1_trend(h, p) = p[1] .+ p[2] .* h
    const p_trend = zeros(2)

    # NOTE: Including quadratic for seasonal did not improve std(anom) in practical use.
    # (rejected alternate parameters commented below)
    #(p[6] .* cos.(2 .* pi .* x[:, 1]) .+  p[7].* sin.(2 .* pi .* x[:, 1])) .* (1 .+ p[8] .* x[:, 2] .+ p[9] .* x[:, 2] .^ 2)

    # Model for all geotiles for a region for a given year: quadratic in elevation
    # p[1]: intercept, p[2]: linear, p[3]: quadratic
    model2(h, p) = p[1] .+ p[2] .* h .+ p[3] .* h .^ 2
    const p2 = zeros(3)
    const lb2 = [-30.0, -0.1, -0.01]
    const ub2 = [+30.0, +0.1, 0.01]

    # Polynomial + seasonal model in time t
    # p = [offset, slope, acceleration, amplitude, phase]
    model3(t, p) = p[1] .+ p[2] .* t .+ p[3] .* t .^ 2 .+ p[4] .* sin.(2 .* pi .* (t .+ p[5]))
    const p3 = zeros(5)

    # Simple seasonal (sine) annual cycle model in time t
    # p = [offset, amplitude, phase]
    model4(t, p) = p[1] .+ p[2] .* sin.(2 .* pi .* (t .+ p[3]))
    const p4 = zeros(3)

    # Offset + trend + seasonal (sinusoidal) terms model for time series offsets
    # p[1]: offset, p[2]: trend, p[3]: seasonal amplitude, p[4]: phase
    offset_trend_seasonal(t, p) =
        p[1] .+                                     # Offset
        p[2] .* t .+                                # Linear trend
        p[3] .* sin.(2 .* pi .* (t .+ p[4]))        # Sinusoidal annual cycle

    # Offset + trend + full seasonality (cosine and sine terms); more flexible phase representation
    # p[1]: offset, p[2]: trend, p[3]: cos coefficient, p[4]: sin coefficient
    offset_trend_seasonal2(t, p) =
        p[1] .+                              # Offset
        p[2] .* t .+                         # Linear trend
        p[3] .* cos.(2π .* t) .+             # Cosine (annual cyclic)
        p[4] .* sin.(2π .* t)                # Sine (annual cyclic)

    # Offset + linear trend + acceleration + full seasonal cycle (cosine and sine terms)
    # p[1]: offset, p[2]: linear, p[3]: quadratic, p[4]: cos, p[5]: sin
    offset_trend_acceleration_seasonal2(t, p) =
        p[1] .+                            # Offset
        p[2] .* t .+                       # Linear trend
        p[3] .* t .^ 2 .+                  # Acceleration (quadratic term)
        p[4] .* cos.(2π .* t) .+           # Cosine (annual)
        p[5] .* sin.(2π .* t)              # Sine (annual)

    # 10th order polynomial 
    polynomial10(x, p) =
        p[1] .+ p[2] .* x .+ p[3] .* x .^ 2 .+ p[4] .* x .^ 3 .+ p[5] .* x .^ 4 .+ p[6] .* x .^ 5 .+ p[7] .* x .^ 6 .+ p[8] .* x .^ 7 .+ p[9] .* x .^ 8 .+ p[10] .* x .^ 9 .+ p[11] .* x .^ 10
    const p10 = zeros(11)
    
    # Initial parameters for offset/trend/seasonal fitting
    const p_offset_trend_seasonal = zeros(4)
end

offset_trend(t, p) = p[1] .+ p[2] .* t;
offset_trend_p = zeros(2);

"""
    replace_with_model!(dh, nobs, geotiles2replace::AbstractArray; mission2replace="hugonnet", missions2align2, missions2update)

Replace elevation change data for specified geotiles with model-fitted values.

# Arguments
- `dh`: Dictionary of DimensionalArrays containing elevation change data by mission
- `nobs`: Dictionary of DimensionalArrays containing observation counts by mission
- `geotiles2replace`: Array of geotile IDs to replace with model-fitted values
- `mission2replace`: Mission whose data will be replaced (default: "hugonnet")
-  `missions2align2`: List of missions to align to the reference mission

# Returns
- Modified `dh` dictionary with replaced values for specified geotiles

# Description
Fits elevation change models to reference mission data and uses these models to replace
values in the target mission for specified geotiles. Uses a combination of reference
missions to ensure data coverage, with the primary reference taking precedence.
"""
function replace_with_model!(dh, nobs, geotiles2replace; missions2replace="hugonnet", missions2align2=["icesat2", "icesat"])

    dgeotile = dims(dh[first(keys(dh))], :geotile)
    dheight = dims(dh[first(keys(dh))], :height)
    ddate = dims(dh[first(keys(dh))], :date)
    dmissions = Dim{:mission}(missions2align2)

    if isnothing(geotiles2replace) || isempty(geotiles2replace) || isempty(missions2replace)
        return dh, nobs
    end

    geotiles2replace = intersect(geotiles2replace, dgeotile)

    Threads.@threads for geotile in geotiles2replace
        for mission2replace in missions2replace
            valid0 = falses(dmissions, ddate, dheight)

            # in order for the data to be replaced with model, all missions2align2 must have data
            for mission in missions2align2
                valid1 = .!isnan.(dh[mission][geotile=At(geotile)])
                if !any(valid1)
                    @warn "No data for $(mission_proper_name(mission)) for geotile: $(geotile), observations not replaced with model"
                    return dh
                end
                valid0[At(mission), :, :] = valid1
            end

            # first missmissions2align2 is the perfered mission and will overwrite any overlaping data
            _, _, vheight = validrange(valid0)
            vdates, _, = validrange(.!isnan.(dh[mission2replace][At(geotile), :, vheight]))

            dh0 = fill(NaN, ddate, dheight)
            nobs0 = fill(NaN, ddate, dheight)
            for mission in missions2align2 # first mission is the perfered mission
                not_valid = isnan.(dh0)
                dh0[not_valid] = dh[mission][At(geotile), :, :][not_valid]
                nobs0[not_valid] = nobs[mission][At(geotile), :, :][not_valid]
            end
            valid0 = .!isnan.(dh0)

            # replace data with model fit
            t0 = decimalyear.(ddate)
            t0 = repeat(t0, 1, length(dheight))

            h0 = val(dheight)'
            h0 = repeat(h0, length(ddate), 1)

            t0 .-= mean(decimalyear.(parent(ddate)[vdates]))
            h0 .-= mean(val(dheight[vheight]))

            # some tiles have nobs == 0 where there is not icesat or iceast-2 data and the 
            # geotile has been filled with with neighbor values.. for this reason we need to add
            # one to all nobs
            nobs0 .+= 1

            # fit global model 
            if length(vheight) == 1
                fit1 = curve_fit(model3, t0[valid0], dh0[valid0], nobs0[valid0], p3)
                dh[mission2replace][At(geotile), vdates, vheight] = model3(vec(t0[vdates, vheight]), fit1.param)
            else
                fit1 = curve_fit(model1, hcat(t0[valid0], h0[valid0]), dh0[valid0], nobs0[valid0], p1; lower=lb1, upper=ub1)
                dh[mission2replace][At(geotile), vdates, vheight] = model1(hcat(vec(t0[vdates, vheight]), vec(h0[vdates, vheight])), fit1.param)
            end

            nobs[mission2replace][At(geotile), vdates, vheight] .= 9999
        end
    end
    return dh, nobs
end


"""
    hyps_model_fill!(dh1, nobs1, params; 
                     bincount_min=5, 
                     model1_nmad_max=5, 
                     smooth_n=9, 
                     smooth_h2t_length_scale=800,
                     variogram_range_ratio=false, 
                     show_times=false) -> Union{Nothing, Dict}

Fill gaps in elevation change data using a spatiotemporal model.

# Arguments
- `dh1`: Dictionary of elevation change DimArrays by mission
- `nobs1`: Dictionary of observation count DimArrays by mission
- `params`: Dictionary of parameter DataFrames by mission
- `bincount_min`: Minimum number of observations required for a valid bin (default: 5)
- `model1_nmad_max`: Maximum MAD normalization threshold for outlier filtering (default: 5)
- `smooth_n`: Number of nearest neighbors for smoothing (default: 9)
- `smooth_h2t_length_scale`: Scaling factor for height relative to time (default: 800)
- `variogram_range_ratio`: Whether to calculate and return variogram range ratios (default: false)
- `show_times`: Whether to display timing information (default: false)

# Returns
- If `variogram_range_ratio=true`, returns a dictionary of range ratios by mission
- Otherwise, modifies `dh1` and `nobs1` in place and returns nothing

# Description
This function fills gaps in elevation change data by:
1. Filtering data based on observation count and outlier detection
2. Fitting a spatiotemporal model to the filtered data
3. Smoothing residuals using k-nearest neighbors
4. Filling gaps with model predictions plus smoothed residuals
"""
function hyps_model_fill!(dh1, nobs1, params; bincount_min=5, model1_nmad_max=5, smooth_n=9, smooth_h2t_length_scale=800, show_times=false, missions2update=nothing)

    if smooth_h2t_length_scale < 1
        error("smooth_h2t_length_scale is < 1, should be in the range 2000 to 1, sypically 800")
    end

    t = decimalyear.(dims(dh1[first(keys(dh1))], :date))
    t = repeat(t, 1, length(dims(dh1[first(keys(dh1))], :height)))

    h = val(dims(dh1[first(keys(dh1))], :height))'
    h = repeat(h, length(dims(dh1[first(keys(dh1))], :date)), 1)

    dgeotiels = dims(dh1[first(keys(dh1))], :geotile)
    foo = DimArray(fill(NaN, length(dgeotiels)), dgeotiels)

    if isnothing(missions2update)
        missions2update = keys(dh1)
    end

    for mission in missions2update

        valid_all_mission = .!isnan.(dh1[mission])
        if !any(valid_all_mission)
            continue
        end

        # find valid range for entire mission so that all filled geotiles cover the same date range
        rrange, = validrange(vec(any(valid_all_mission, dims=(1, 3))))

        Threads.@threads for geotile in dims(dh1[mission], :geotile)

            show_times ? t1 = time() : nothing

            k = findfirst(params[mission].geotile .== geotile)
            df = @view params[mission][k, :]

            dh0 = dh1[mission][geotile=At(geotile)]
            nobs0 = nobs1[mission][geotile=At(geotile)]
            df.nobs_raw = sum(nobs0)
            df.nbins_raw = sum(nobs0 .> 0)

            ###################################### FILTER 1 ################################
            valid1 = .!isnan.(dh0) .& (nobs0 .> bincount_min[mission]) .& (abs.(dh0) .< 200)
            ################################################################################

            dh0[collect(.!valid1)] .= NaN
            dh1[mission][geotile=At(geotile)] = dh0
            nobs0[.!valid1] .= 0
            nobs1[mission][geotile=At(geotile)] = nobs0

            # if there are not enough points to fit a model then set all to NaNs
            va = sum(valid1)
            if va <= (length(p1) + 2)
                dh1[mission][geotile=At(geotile)] .= NaN
                nobs1[mission][geotile=At(geotile)] .= 0
                continue
            end

            # determine valid range of data
            (_, crange) = validrange(valid1)
            valid0 = valid1[rrange, crange]

            # if there are not enough points to fit a model then set all to NaNs
            va = sum(valid0)
            if va <= (length(p1) + 2)
                dh1[mission][At(geotile), :, :] .= NaN
                nobs1[mission][At(geotile), :, :] .= 0
                continue
            end

            dh0 = dh0[rrange, crange]
            nobs0 = nobs0[rrange, crange]
            t0 = t[rrange, crange]
            h0 = h[rrange, crange]
    
            show_times ? t2 = time() : nothing
            show_times && printstyled("initial data selection $mission: $(round(t2 - t1; digits=2))s\n", color=:yellow)

            # center predictors and observations
            t0_mean = round(mean(t0)) # remove an exact integer to keep phase from shifting
            df.t0 = t0_mean
            t0 = t0 .- t0_mean

            h0_mean = mean(h0[:])
            df.h0 = h0_mean
            h0 = h0 .- h0_mean

            dh0_median = median(dh0[valid0]) + 0.0000001 # add a small offset to prevent numerical instability
            if isnan(dh0_median)
                error("dh0_median is NaN")
            end

            df.dh0 = dh0_median
            dh0 = dh0 .- dh0_median

            df.bin_std = std(dh0[valid0])

            show_times ? t3 = time() : nothing
            show_times && printstyled("statistics $mission: $(round(t3 - t2; digits=2))s\n", color=:yellow)
            # fit global model 

            fit1 = []
            try
                fit1 = curve_fit(model1, hcat(t0[valid0], h0[valid0]), dh0[valid0], nobs0[valid0], p1; lower=lb1, upper=ub1)
            catch
                fit1 = curve_fit(model1, hcat(t0[valid0], h0[valid0]), dh0[valid0], nobs0[valid0], p1)
            end

            show_times ? t4 = time() : nothing
            show_times && printstyled("model fit $mission: $(round(t4 - t3; digits=2))s\n", color=:yellow)

            dh0_mdl = model1(hcat(t0[valid0], h0[valid0]), fit1.param)
            dh0_anom = dh0[valid0] .- dh0_mdl

            ###################################### FILTER 2 ####################################
            # filter model1_nmad_max sigma outliers
            valid0[valid0] = nmad(dh0_anom) .<= model1_nmad_max
            vb = sum(valid0)
            df.nbins_filt1 = vb

            # if there are not enough points to fit a model then set all to NaNs
            if vb <= (length(p1) + 2)
                dh1[mission][geotile = At(geotile)] .= NaN
                nobs1[mission][geotile = At(geotile)] .= 0
                continue

            end

            if vb < va
                dh0[.!valid0] .= NaN
                nobs0[.!valid0] .= 0
                fit1 = curve_fit(model1, hcat(t0[valid0], h0[valid0]), dh0[valid0], nobs0[valid0], p1; lower=lb1, upper=ub1)

                dh0_mdl = model1(hcat(t0[valid0], h0[valid0]), fit1.param)
                dh0_anom = dh0[valid0] .- dh0_mdl
            end

            df.bin_anom_std = std(dh0_anom)
            nobs1[mission][At(geotile), rrange, crange] = nobs0

            # add final parameters to DataFrame
            df.nobs_final = sum(nobs0)
            df.param_m1 = fit1.param

            show_times ? t5 = time() : nothing
            show_times && printstyled("filter 2 $mission: $(round(t5 - t4; digits=2))s\n", color=:yellow)


            #### THIS CODE BLOCK TAKES THE MOST TIME ####
            # take the median of the x closest neighbors
            if sum(valid0) < smooth_n[mission]
                anom_smooth = zeros(length(dh0))
            else

                # scale height distance relative to time (i.e. length-scale) 
                pts = hcat(t0[valid0], h0[valid0] / smooth_h2t_length_scale)'

                kdtree = KDTree(pts)
                idxs, _ = knn(kdtree, pts, smooth_n[mission])

                anom0 = map(ind -> median(dh0_anom[ind]), idxs)

                # extrema(anom0)
                # interpolate anomalies using weighted distance (Shepard(2))
                itp = ScatteredInterpolation.interpolate(Shepard(2), pts, anom0)
                pts = hcat(t0[:], h0[:] / smooth_h2t_length_scale)'
                anom_smooth = vec(ScatteredInterpolation.evaluate(itp, pts))
            end


            # fill out valid range (no extraploation beyond (rrange,crange) of geotile with the model only
            dh1[mission][At(geotile), rrange, crange] = model1(hcat(t0[:], h0[:]), fit1.param) .+ dh0_median .+ anom_smooth
            # println("granule interp: $(mission) - $(geotile)")
            show_times ? t6 = time() : nothing
            show_times && printstyled("granule interp $mission: $(round(t6 - t5; digits=2))s\n", color=:yellow)
        end
    end

    return dh1, nobs1, params
end

"""
    hyps_remove_land_surface_trend!(dh1; missions2update=nothing, remove_land_surface_trend=nothing)

Remove erroneous land surface trends from elevation change data.

# Arguments
- `dh1`: Dictionary of elevation change DimArrays indexed by mission
- `missions2update`: Missions to process (default: all missions)
- `remove_land_surface_trend`: Array of trend values to remove per mission (default: nothing)

# Description
For each mission, removes a linear trend from the elevation change data if specified.
The trend is centered around the mission's temporal midpoint to minimize impact on
the overall elevation change signal. This correction is typically applied to remove
artificial trends that may arise from land surface processing artifacts.

# Returns
- Modified `dh1` with land surface trends removed
"""
function hyps_remove_land_surface_trend!(dh1; missions2update=nothing, remove_land_surface_trend=nothing)

    dgeotile = dims(dh1[first(keys(dh1))], :geotile)
    dheight = dims(dh1[first(keys(dh1))], :height)
    ddate = dims(dh1[first(keys(dh1))], :date)
    decyear = decimalyear.(ddate)
    
    for mission in missions2update

        if !isnothing(remove_land_surface_trend) && (remove_land_surface_trend[At(mission)] != 0)
    
            # center date around mission center date
            valid1 = .!isnan.(dh1[mission])
            if !any(valid1)
                continue
            end

            _, date_range, _ = validrange(valid1)
            mid_date = mean(decyear[date_range])
            delta = (decyear .- mid_date) .* remove_land_surface_trend[At(mission)]

            for geotile in dgeotile
                for height in dheight
                    dh1[mission][geotile=At(geotile), height=At(height)] .-= delta
                end
            end
        end
    end
end


"""
    hyps_amplitude_normalize!(dh1, params, params_reference)

Normalize seasonal amplitude of elevation change data between missions.

# Arguments
- `dh1`: Dictionary of elevation change DimArrays indexed by mission
- `params`: Dictionary of parameter DataFrames indexed by mission  
- `params_reference`: Dictionary of parameter DataFrames for reference mission

# Description
For each mission, calculates the difference in seasonal components between
the mission's model and the reference mission's model. Adds this difference to normalize
the seasonal amplitude of the target mission's data to match the reference mission.
"""
function hyps_amplitude_normalize!(dh1, params, params_reference)

    t = decimalyear.(dims(dh1, :date))
    t = repeat(t, 1, length(dims(dh1, :height)))

    h = val(dims(dh1, :height))'
    h = repeat(h, length(dims(dh1, :date)), 1)

    Threads.@threads for geotile in dims(dh1, :geotile)
        k = findfirst(params.geotile .== geotile)
        df0 = params[k, :]
        dfr = params_reference[k, :]

        if any(isnan.(df0.param_m1)) || any(isnan.(dfr.param_m1))
            continue
        end

        dh0 = dh1[At(geotile), :, :]
        valid = .!isnan.(dh0)

        (rrange, crange) = validrange(valid)

        dh0 = dh0[rrange, crange]
        valid = valid[rrange, crange]

        # only last 6 fit parameters are related to seasonal cycle
        p0 = df0.param_m1[6:end]
        p_ref = dfr.param_m1[6:end]

        # parameters are tied to a relative height and time for eartlier fitting normalization
        t0 = t[rrange, crange] .- df0.t0
        h0 = h[rrange, crange] .- df0.h0
        model0 = DimArray(reshape(model1_seasonal(hcat(t0[:], h0[:]), p0), size(dh0)), dims(dh0))

        t0 = t[rrange, crange] .- dfr.t0
        h0 = h[rrange, crange] .- dfr.h0
        model_ref = DimArray(reshape(model1_seasonal(hcat(t0[:], h0[:]), p_ref), size(dh0)), dims(dh0))

        delta = model_ref .- model0

        if !any(isnan.(delta)) # there seem to be rare cases where model1_seasonal returns nans.
            dh1[At(geotile), rrange, crange] = dh0 .+ delta
        end
    end
end

"""
    hyps_align_dh!(dh, nobs, params, area_km2; missions2align2=["icesat2", "icesat"], missions2update=nothing)

Align elevation change data between different altimetry missions by calculating and applying weighted offsets.

# Arguments
- `dh`: Dict of elevation change DimArrays by mission
- `nobs`: Dict of observation count DimArrays by mission  
- `params`: Dict of parameter DataFrames by mission
- `area_km2`: DimArray of area values per geotile
- `missions2align2`: Reference missions to align others to (default: ["icesat2", "icesat"])

# Returns
- Modified `dh` with aligned elevation changes
- Modified `params` with offset parameters
"""
function hyps_align_dh!(dh, nobs, params, area_km2; missions2align2=["icesat2", "icesat"], missions2update=nothing)
    # Get dimensions from first mission
    ddate = dims(dh[first(keys(dh))], :date)
    dgeotile = dims(dh[first(keys(dh))], :geotile)

    # Initialize area-averaged elevation changes and observation counts
    area_avg_dh = Dict()
    nobs_sum = Dict()

    if isnothing(missions2update)
        missions2update = keys(dh)
    end

    for mission in keys(dh)
        area_avg_dh[mission] = fill(NaN, dgeotile, ddate)
        nobs_sum[mission] = zeros(Int, dgeotile, ddate)
    end

    # Calculate area-averaged elevation changes per geotile
    Threads.@threads for geotile in dgeotile
        for mission in keys(dh)
            area_avg_dh[mission][geotile=At(geotile)] = dh_area_average(dh[mission][geotile=At(geotile)], area_km2[geotile=At(geotile)])
            nobs_sum[mission][geotile=At(geotile)] = sum(nobs[mission][At(geotile), :, :], dims=:height)
        end
    end

    # Align non-reference missions to reference missions
    missions2align = setdiff(missions2update, missions2align2)

    Threads.@threads for geotile in dgeotile
        for mission in missions2align
            index_table = findfirst(params[mission][!, :geotile] .== geotile)
            offset0 = 0.0
            weight0 = 0

            # Calculate weighted offset from each reference mission
            for mission_ref in missions2align2
                offest1 = area_avg_dh[mission][geotile=At(geotile)] .- area_avg_dh[mission_ref][geotile=At(geotile)]
                index = .!isnan.(offest1)

                if any(index)
                    # Store offset statistics
                    params[mission][index_table, "offset_$mission_ref"] = median(offest1[index])
                    params[mission][index_table, "offset_nmad_$mission_ref"] = mean(nmad(offest1[index]))
                    params[mission][index_table, "offset_nobs_$mission_ref"] = sum(nobs[mission][geotile=At(geotile)])

                    # I don't think nmad if very relevant here as the data being adjusted can be pretty noisy
                    # weight1 = 1 / ((params[mission][index_table, "offset_nmad_$mission_ref"] / (sqrt(params[mission][index_table, "offset_nobs_$mission_ref"]))) .^ 2)
                    weight1 = sqrt(params[mission_ref][index_table, "offset_nobs_$mission_ref"] + params[mission][index_table, "offset_nobs_$mission_ref"])
                    offset0 += params[mission][index_table, "offset_$mission_ref"] * weight1
                    weight0 += weight1
                end
            end

            # Apply weighted average offset
            if weight0 > 0
                params[mission][index_table, "offset"] = offset0 / weight0

                # NOTE: offset is subtracted from height anomalies
                dh[mission][geotile=At(geotile)] .-= params[mission][index_table, "offset"]
            end
        end
    end

    return dh, params
end


"""
    hyps_fill_empty!(dh1, params, geotile_extent, area_km2; missions2update, mask=:glacier)

Fill empty geotiles with data interpolated from nearby geotiles.

# Arguments
- `dh1`: Dictionary of elevation change DimArrays by mission
- `params`: Dictionary of parameter DataFrames by mission  
- `geotile_extent`: Array of geotile extents for spatial calculations
- `area_km2`: DimArray containing area information for each geotile
- `mask`: Symbol specifying which area mask to use (default: `:glacier`)

# Description
For each mission and geotile:
1. Identifies geotiles with no data but containing the specified mask type (e.g., glacier)
2. Finds the nearest 5 geotiles with valid data at overlapping elevation ranges
3. Interpolates missing data using the median values from neighboring geotiles
4. Fills any remaining gaps along elevation profiles using linear interpolation
5. Sets geotiles with no mask area to zero

NOTE: if valid data extends beyond elevation range of surface_mask then extents of valid output data can differ.. this is not a problem

# Returns
- Modified `dh1` dictionary with filled geotiles
"""
function hyps_fill_empty!(dh1, params, geotile_extent, area_km2; missions2update=nothing)

    # ensure that dims match
    dgeotile = dims(dh1[first(keys(dh1))], :geotile)

    geotile_extent = geotile_extent[geotile=At(collect(dgeotile))]
    area_km2 = area_km2[geotile=At(collect(dgeotile))]

    geotile_rectangles = extent2rectangle.(geotile_extent)
    lonlat = GO.centroid.(geotile_rectangles)
    mission_specs = project_products(; project_id=:v01)

    if isnothing(missions2update)
        missions2update = keys(dh1)
    end

    for mission in missions2update
        # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><
        #mission = "hugonnet"
        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

        # if there are not enough points to fit a model then set all to NaNs
        valid1 = .!isnan.(dh1[mission])
        if !any(valid1)
            continue
        end

        # find valid range for entire mission so that all filled geotiles cover the same date range
        rrange, = validrange(vec(any(valid1, dims=(1, 3))))
        dh0_median = getindex.(params[mission][:, :param_m1], 1) .+ params[mission][:, :dh0]

        # copy dh as it gets morphed inside of parallel loop
        dh2 = deepcopy(dh1[mission])

        # rectangle of valid mission bounds
        mission_extent = (X=mission_specs[Symbol(mission)].longitude_limits, Y=mission_specs[Symbol(mission)].latitude_limits)
        mission_rectangle = extent2rectangle(mission_extent)

        within_mission_rectangles = GO.intersects.(geotile_rectangles, Ref(mission_rectangle))

        Threads.@threads for geotile in (dims(dh1[mission], :geotile)[within_mission_rectangles])
           
            # fill with median of nearest neighbors if no data
            has_ice = area_km2[geotile=At(geotile)] .> 0

            if !any(has_ice)
                #printstyled("$geotile $mask area is zero: setting dh .= 0 \n"; color=:blue) 
                dh0 = @view dh1[mission][At(geotile), rrange, :]
                dh0 .= 0
                continue
            else
                crange, = validrange(has_ice)
            end

            dh0 = @view dh1[mission][At(geotile), rrange, crange]

            if all(isnan.(dh0))

                # find distance between goetiles
                dist2geotile = haversine.(Ref(lonlat[geotile=At(geotile)]), lonlat, 6371000)

                # find closest X valid geotiles
                nnearest = 5

                #sorted index
                idx = sortperm(dist2geotile)

                # has data in desired elevation range
                has_data = vec(any(.!isnan.(dh2[:, rrange, crange]), dims=(2, 3)))
                has_data = has_data[idx]

                # check that there are enough geotiles with overlapping elevation ranges
                nnearest0 = sum(has_data)
                if nnearest0 < nnearest
                    if nnearest0 == 0
                        printstyled("$mission $geotile: no geotiles with overlapping elevation ranges to fill geotile, setting dh .= 0 \n"; color=:red, bold=true)
                        dh0 .= 0
                        continue
                    else
                        printstyled("$mission $geotile: less than $nnearest geotiles with overlapping elevation ranges to fill geotile, using closest $nnearest0 geotiles instead\n"; color=:light_red, bold=false)
                        nnearest = nnearest0
                    end
                end

                idx = idx[has_data][1:nnearest]
                f1 = copy(dh2[idx, rrange, crange])

                # remove mean offset to normalize between regions
                dh0_median0 = dh0_median[idx]
                for g in eachindex(dh0_median0)
                    f1[g, :, :] .-= dh0_median0[g]
                end

                for i in 1:length(dims(dh0, :date))
                    for j in 1:length(dims(dh0, :height))
                        f2 = vec(f1[:, i, j])
                        valid00 = .!isnan.(f2)
                        if any(valid00)
                            dh0[i, j] = median(f2[valid00])
                        end
                    end
                end

                # if there are any gaps along an elevation profile then they need to be filled
                (rrange0, crange0) = validrange(.!isnan.(dh0))
                if any(isnan.(dh0[rrange0, crange0]))
                    x = 1:size(dh0, 2)
                    for i in rrange0
                        y = vec(dh0[i, :])
                        valid = .!isnan.(y)
                        if any(.!valid)
                            itp = DataInterpolations.LinearInterpolation(y[valid], x[valid]; extrapolation=ExtrapolationType.Linear)
                            y[.!valid] = itp(x[.!valid])
                            dh0[i, :] = y
                        end
                    end
                end

                # add back median model offset 
                dh0 .+= median(dh0_median0)
            end
        end
    end
    return dh1
end


"""
    hyps_fill_updown!(dh1, area_km2; missions2update=nothing)

Fill elevation gaps at the lowest and highest elevations by extending the nearest valid values.

# Arguments
- `dh1`: Dictionary of elevation change DimArrays by mission
- `area_km2`: DimArray of glacier area

# Returns
- Modified `dh1` with filled elevation gaps

# Description
For each geotile, fills gaps at the lowest and highest elevations by extending the nearest valid 
elevation values outward. This ensures complete elevation coverage for each time step.
"""
function hyps_fill_updown!(dh1, area_km2; missions2update=nothing)

    if isnothing(missions2update)
        missions2update = keys(dh1)
    end

    for mission in missions2update
        # if there are not enough points to fit a model then set all to NaNs
        valid1 = .!isnan.(dh1[mission])
        if !any(valid1)
            continue
        end

        rrange, = validrange(vec(any(valid1, dims=(:geotile, :height))))
        dgeotile = dims(dh1[mission], :geotile)

        Threads.@threads for geotile in dgeotile

            # fill with median of nearest neighbors if no data
            valid = area_km2[geotile=At(geotile)] .> 0

            if !any(valid)
                continue
            else
                crange, = validrange(valid)
            end

            dh0 = @view dh1[mission][At(geotile), rrange, crange]

            # fill lowest and highest elevations with first and last values
            valid = .!isnan.(dh0)
            if .!any(valid) || all(valid)
                continue
            end

            for i in 1:length(dims(dh0, 1))
                if any(valid[i, :])
                    f = findfirst(vec(valid[i, :]))
                    l = findlast(vec(valid[i, :]))
                    dh0[i, 1:f] .= dh0[i, f]
                    dh0[i, l:end] .= dh0[i, l]
                end
            end
        end
    end
    return dh1
end

"""
    binned_filepath(;binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct)

Generate the filepath for binned elevation change data.

# Arguments
- `binned_folder`: Directory where binned files are stored
- `surface_mask`: Type of surface mask applied (e.g., "ice", "land")
- `dem_id`: Identifier for the DEM used
- `binning_method`: Method used for binning data
- `project_id`: Project identifier
- `curvature_correct`: Boolean indicating if curvature correction was applied

# Returns
- Full filepath to the binned data file with appropriate naming convention
"""
function binned_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct)
    if curvature_correct
        runid = "$(surface_mask)_dh_$(dem_id)_cc_$(binning_method)_$(project_id)"
    else
        runid = "$(surface_mask)_dh_$(dem_id)_$(binning_method)_$(project_id)"
    end

    binned_file = joinpath(binned_folder, "$(runid).jld2")

    return binned_file
end

"""
    binned_filled_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, fill_param)

Generate filepath for filled binned elevation change data and corresponding figure suffix.

# Arguments
- `binned_folder`: Path to the binned folder
- `surface_mask`: Type of surface mask applied (e.g., "ice", "land")
- `dem_id`: Identifier for the DEM used
- `binning_method`: Method used for binning data
- `project_id`: Project identifier
- `curvature_correct`: Boolean indicating if curvature correction was applied
- `amplitude_correct`: Boolean indicating if amplitude correction was applied
- `fill_param`: Parameter set identifier used for filling

# Returns
- `binned_filled_file`: Full filepath to the filled binned data file
- `figure_suffix`: Suffix string for related figure filenames
"""
function binned_filled_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, fill_param)

    if curvature_correct
        runid = "$(surface_mask)_dh_$(dem_id)_cc_$(binning_method)_$(project_id)"
    else
        runid = "$(surface_mask)_dh_$(dem_id)_$(binning_method)_$(project_id)"
    end

    if amplitude_correct
        binned_filled_file = joinpath(binned_folder, "$(runid)_filled_ac_p$(fill_param)_aligned.jld2")
    else
        binned_filled_file = joinpath(binned_folder, "$(runid)_filled_p$(fill_param)_aligned.jld2")
    end

    figure_suffix = splitpath(binned_filled_file)
    figure_suffix = figure_suffix[end]
    figure_suffix = replace(figure_suffix, ".jld2" => "")
    figure_suffix = replace(figure_suffix, "dh" => "dm")

    return binned_filled_file
end

function binned2filled_filepath(;binned_file, amplitude_correct, fill_param)

    binned_file0 = replace(binned_file, ".jld2" => "")
    
    if amplitude_correct
        binned_filled_file = "$(binned_file0)_filled_ac_p$(fill_param)_aligned.jld2"
    else
        binned_filled_file = "$(binned_file0)_filled_p$(fill_param)_aligned.jld2"
    end

    figure_suffix = splitpath(binned_filled_file)
    figure_suffix = figure_suffix[end]
    figure_suffix = replace(figure_suffix, ".jld2" => "")
    figure_suffix = replace(figure_suffix, "dh" => "dm")

    return binned_filled_file, figure_suffix
end


"""
    binned_synthesized_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, fill_param)

Generate filepath for synthesized binned elevation change data.

# Arguments
- `binned_folder`: Path to the binned folder
- `surface_mask`: Type of surface mask applied (e.g., "ice", "land")
- `dem_id`: Identifier for the DEM used
- `binning_method`: Method used for binning data
- `project_id`: Project identifier
- `curvature_correct`: Boolean indicating if curvature correction was applied
- `amplitude_correct`: Boolean indicating if amplitude correction was applied
- `fill_param`: Parameter set identifier used for filling

# Returns
- `binned_synthesized_file`: Full filepath to the synthesized binned data file
"""
function binned_synthesized_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, fill_param)
    binned_filled_file = binned_filled_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, fill_param)
    binned_synthesized_file = replace(binned_filled_file, ".jld2" => "_synthesized.jld2")
    return binned_synthesized_file
end


"""
    binned_filled_fileparts(binned_filled_file) -> NamedTuple

Extract metadata components from a binned filled file path.

# Arguments
- `binned_filled_file`: Path to a binned filled file

# Returns
A NamedTuple containing the following components:
- `surface_mask`: Type of surface mask applied (e.g., "ice", "land")
- `dem`: Digital elevation model identifier
- `curvature_correct`: Boolean indicating if curvature correction was applied
- `binning_method`: Method used for binning data
- `project_id`: Project identifier
- `amplitude_correct`: Boolean indicating if amplitude correction was applied
- `fill_param`: Parameter set identifier used for filling

# Description
Parses a binned filled filename to extract metadata components based on the standardized
naming convention used in the project.
"""
function binned_filled_fileparts(binned_filled_file)

    binned_folder = splitpath(replace(binned_filled_file, pathlocal[:data_dir] => ""))[1]

    geotile_width = parse(Int, binned_filled_file[findfirst("deg/", binned_filled_file)[1]-1])


    binned_filled_file = splitpath(binned_filled_file)[end]

    # plit file name at "_"
    param = split(binned_filled_file, "_")

    p = 1
    surface_mask = param[p]

    # skip "dm" in file name
    p += 1

    if !(param[p] == "dh")
        surface_mask = "$(surface_mask)_$(param[p])"
        p += 1
    end

    p += 1

    dem = param[p]

    if dem == "cop30"
        dem = "cop30_v2"
        p += 1
    end

    p += 1
    if param[p] == "cc"
        curvature_correct = true
        p += 1
    else
        curvature_correct = false
    end

    binning_method = param[p]

    p += 1
    project_id = param[p]

    # skip "filled"
    p += 2
    if param[p] == "ac"
        amplitude_correct = true
        p += 1
    else
        amplitude_correct = false
    end

    fill_param = parse(Int16, param[p][2])

    out = (; binned_folder, surface_mask, dem, curvature_correct, binning_method, project_id, amplitude_correct, fill_param, geotile_width)
    return out
end

"""
    geotile_bin2d(
        df; 
        var2bin="dh",
        dims_edges=("decyear" => 1990:(30/365):2026, "height_reference" => 0.:100.:10000.),
        binfunction::T = binningfun_define(binning_method)
    ) where {T <: Function} -> Tuple{Union{Nothing, DimArray}, Union{Nothing, DimArray}}

Bin data in a DataFrame by date and elevation into a 2D grid.

# Arguments
- `df`: DataFrame containing data to bin
- `var2bin`: Column name of the variable to bin (default: "dh")
- `dims_edges`: Tuple of dimension names and bin edges (default: decyear and height_reference)
- `binfunction`: Function to apply within each bin (default: defined by binning_method)

# Returns
- `var0`: DimArray of binned values or nothing if no valid data
- `nobs0`: DimArray of observation counts or nothing if no valid data

# Description
Creates a 2D grid of binned data based on the specified dimensions (typically time and elevation).
Filters data to the range of available values, applies the binning function to each cell,
and returns both the binned values and observation counts as DimArrays.
"""
function geotile_bin2d(
    df;
    var2bin="dh",
    dims_edges=("decyear" => 1990:(30/365):2026, "height_reference" => 0.:100.:10000.),
    binfunction::T=binningfun_define(binning_method)
) where {T<:Function}

    # bin data by date and elevation
    minmax_date = extrema(df[:, dims_edges[1][1]])
    minmax_height = extrema(df[:, dims_edges[2][1]])

    Δd = dims_edges[1][2][2] - dims_edges[1][2][1]
    Δh = dims_edges[2][2][2] - dims_edges[2][2][1]

    date_ind = (dims_edges[1][2] .>= minmax_date[1] - Δd) .&
               (dims_edges[1][2] .<= (minmax_date[2] + Δd))

    date_ind_center = findall(date_ind)[1:end-1]

    height_ind = (dims_edges[2][2] .>= minmax_height[1] - Δh) .&
                 (dims_edges[2][2] .<= (minmax_height[2] + Δh))

    height_ind_center = findall(height_ind)[1:end-1]

    nobs0 = nothing
    var0 = nothing

    # check bounds: binstats will throw an error if no data is passed to median()
    if !any(date_ind)
        return var0, nobs0
    end

    if !vector_overlap(df[!, dims_edges[1][1]], dims_edges[1][2][date_ind]) ||
       !vector_overlap(df[!, dims_edges[2][1]], dims_edges[2][2][height_ind])

        return var0, nobs0
    end

    df = binstats(df, [getindex.(dims_edges, 1)...], [getindex.(dims_edges, 2)...],
        var2bin; col_function=[binfunction], missing_bins=true)

    gdf = DataFrames.groupby(df, dims_edges[1][1])

    dd1 = Dim{Symbol(dims_edges[1][1])}(sort((dims_edges[1][2][1:end-1] .+ dims_edges[1][2][2:end]) ./ 2))
    dd2 = Dim{Symbol(dims_edges[2][1])}(sort((dims_edges[2][2][1:end-1] .+ dims_edges[2][2][2:end]) ./ 2))

    nobs0 = fill(0, (dd1, dd2))
    var0 = fill(NaN, (dd1, dd2)) # this should really match the type of the input data ... but hey... this is easier right now

    # get index into sorted array
    p = sortperm(gdf[1][:, dims_edges[2][1]])

    for (i, df) in enumerate(gdf)

        isval = .!ismissing.(df[p, "dh_function"])
        var2 = @view var0[i, :]
        nobs2 = @view nobs0[i, :]
        if any(isval)
            var2[isval] = df[p, "dh_function"][isval]
            nobs2[isval] = df.nrow[p][isval]
        end
    end

    return var0, nobs0
end

"""
    dh_area_average(dh, area0) -> dh_area_avg

Calculate area-weighted average of elevation change (dh) across all dimensions except date.

# Arguments
- `dh`: DimArray of elevation changes with dimensions including :date
- `area0`: Array of areas corresponding to the spatial dimensions of dh

# Returns
- `dh_area_avg`: DimArray of area-weighted average elevation changes with dimension :date

# Description
For each date, calculates the area-weighted mean of elevation change values where valid (non-NaN) 
data exists. The weighting is based on the corresponding areas in area0.
"""
function dh_area_average(dh, area0)

    ddate = dims(dh, :date);
    dh_area_avg = fill(NaN, dims(dh, :date))

    if !all(isnan.(dh))
        for date in ddate
            dh0 = dh[date = At(date)]
            valid0 = .!isnan.(dh0)
            if any(valid0)
                dh_area_avg[At(date)] = sum(dh0[valid0] .* (area0[valid0])) / sum(area0[valid0])
            end
        end
    end

    return dh_area_avg
end


"""
    binned_filled_filepaths(params)

Generate file paths for binned and filled data files based on parameter configurations.

# Arguments
- `params`: Vector of parameter dictionaries containing configuration settings for each data file

# Returns
- `path2runs`: Vector of valid file paths that exist in the filesystem

# Description
Iterates through parameter configurations to generate file paths for binned and filled data.
Only includes paths for files that actually exist in the filesystem. Each parameter set
is used to construct a file path using the `binned_filled_filepath` function.
"""
function binned_filled_filepaths(params; include_existing_files_only=false)
    path2runs = fill("", length(params))

    for i in eachindex(params)
        path2runs[i] = binned_filled_filepath(; params[i]...)
    end

    if include_existing_files_only
        file_exists = isfile.(path2runs)
        return path2runs[file_exists], params[file_exists]
    else
        return path2runs, params
    end
end

function binned_filled_filepaths(; project_id, surface_masks, dem_ids, curvature_corrects, amplitude_corrects, binning_methods, fill_params, binned_folders, include_existing_files_only=false)
    
    if !isa(project_id, Vector)
        project_id = [project_id]
    end

    param_nt = (; project_id=project_id, surface_mask=surface_masks, dem_id=dem_ids, curvature_correct=curvature_corrects, amplitude_correct=amplitude_corrects, binning_method=binning_methods, fill_param=fill_params, binned_folder=binned_folders)

    params = ntpermutations(param_nt)

    # only include files that exist
    path2binned_filled_files, params = binned_filled_filepaths(params; include_existing_files_only)

    return path2binned_filled_files, params
end