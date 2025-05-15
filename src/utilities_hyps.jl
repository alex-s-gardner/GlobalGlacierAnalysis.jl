
# define models
begin
    # Define model that will be fit to all data binned by hypsometry
    #model::Function = model(t, h; t_intercept=2010) = hcat(ones(size(t)), (t .- t_intercept), h, h.^2, cos.(2 * pi * t), sin.(2 * pi * t))
    model1::Function = model1(x, p) =
        p[1] .+
        p[2] .* x[:, 1] .+
        p[3] .* x[:, 1] .^ 2 .+
        p[4] .* x[:, 2] .+
        p[5] .* x[:, 2] .^ 2 .+
        sin.(2 .* pi .* (x[:, 1] .+ p[6])) .* (p[7] .+ p[8] .* x[:, 2] .* p[9] .* x[:, 2] .^ 2)

    const p1 = zeros(9);
    const lb1 = [-10.0, -3.0, -2.0, -0.05, -0.0001, -1.0, -7.0, -0.05, -0.001];
    const ub1 = [+10.0, +3.0, +2.0, +0.05, +0.0001, +1.0, +7.0, +0.05, +0.001];


    # seasonal only [amplitude is a quadratic function of elevation]
    model1_seasonal::Function = model1_seasonal(x, p) = 
        sin.(2 .* pi .* (x[:, 1] .+ p[6])) .* 
        (p[7] .+ p[8] .* x[:, 2] .* p[9] .* x[:, 2] .^ 2)

    # linear trend
    model_trend::Function = model1_trend(h, p) = p[1] .+ p[2] .* h;
    const p_trend = zeros(2);

    # including quadratic for seasonal does not improve std(anom)
    #(p[6] .* cos.(2 .* pi .* x[:, 1]) .+  p[7].* sin.(2 .* pi .* x[:, 1])) .* (1 .+ p[8] .* x[:, 2] .+ p[9] .* x[:, 2] .^ 2)
    #p1 = zeros(9);
    #lb1 = [-10., -3., -2., -.05, -0.0001, -10., -10., -0.01, -0.0001];
    #ub1 = [+10., +3., +2., +.05, +0.0001, +10., +10., +0.01, +0.0001];

    # model fit across all geotiles for a region for a given year
    model2::Function = model2(h, p) = p[1] .+ p[2] .* h .+ p[3] .* h .^ 2;
    const p2 = zeros(3);
    const lb2 = [-30.0, -0.1, -0.01];
    const ub2 = [+30.0, +0.1, 0.01];

    model3::Function = model3(t, p) = p[1] .+ p[2] .* t .+ p[3] .* t .^ 2 .+ p[4] .* sin.(2 .* pi .* (t .+ p[5]))
    const p3 = zeros(5)

    model4::Function = model4(t, p) = p[1] .+ p[2] .* sin.(2 .* pi .* (t .+ p[3]))
    const p4 = zeros(3)

    offset_trend_seasonal::Function =
        offset_trend_seasonal(t, p) =
        p[1] .+ 
        p[2] .* t .+ 
        p[3] .* sin.(2 .* pi .* (t .+ p[4]))


    offset_trend_seasonal2::Function =
        offset_trend_seasonal2(t, p) =
            p[1] .+
            p[2] .* t .+
            p[3] .* sin.(2π .* t) .+
            p[4] .* cos.(2π .* t)
            

    offset_trend_acceleration_seasonal2::Function =
        offset_trend_acceleration_seasonal2(t, p) =
            p[1] .+
            p[2] .* t .+
            p[3] .* t.^2 .+
            p[4] .* sin.(2π .* t) .+
            p[5] .* cos.(2π .* t)

    const p_offset_trend_seasonal = zeros(4)
end

"""
    geotiles_mutually_exclusive_rgi!(geotiles) -> (geotiles, reg)

Make RGI (Randolph Glacier Inventory) regions mutually exclusive in the geotiles DataFrame.

# Arguments
- `geotiles`: DataFrame containing geotile information with columns for RGI regions (prefixed with "rgi")

# Returns
- `geotiles`: Modified DataFrame with mutually exclusive RGI regions
- `reg`: Vector of column names corresponding to RGI regions

# Description
For each geotile, identifies the RGI region with the largest overlap value and sets all other
region values to zero, ensuring each geotile is assigned to exactly one RGI region.
"""
function geotiles_mutually_exclusive_rgi!(geotiles) 
    reg = names(geotiles);
    reg = reg[startswith.(reg, "rgi")];
    rgi_id = replace.(reg, Ref("rgi" => ""));
    rgi_id = parse.(Int64, rgi_id);

    # make regions mutually exclusive by assigning geotiles to region of largest overlap
    for geotile in eachrow(geotiles[!, reg])
        maxind = findfirst(maximum(geotile) .== collect(geotile))
        for r in eachindex(reg)
            if r == maxind
                continue
            else
                geotile[r] = 0
            end
        end
    end
    return geotiles, reg
end


"""
    hyps_model_fill!(dh1, nobs1, params; 
                     bincount_min=5, 
                     model1_madnorm_max=5, 
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
- `model1_madnorm_max`: Maximum MAD normalization threshold for outlier filtering (default: 5)
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
function hyps_model_fill!(dh1, nobs1, params; bincount_min=5, model1_madnorm_max=5, smooth_n=9, smooth_h2t_length_scale=800, variogram_range_ratio = false, show_times=false)

    if smooth_h2t_length_scale < 1
        error("smooth_h2t_length_scale is < 1, should be in the range 2000 to 1, sypically 800")
    end
   
    t = Altim.decimalyear.(dims(dh1[first(keys(dh1))], :date))
    t = repeat(t, 1, length(dims(dh1[first(keys(dh1))], :height)))

    h = val(dims(dh1[first(keys(dh1))], :height))'
    h = repeat(h, length(dims(dh1[first(keys(dh1))], :date)), 1)
    
    dgeotiels = dims(dh1[first(keys(dh1))], :geotile)
    foo = DimArray(fill(NaN, length(dgeotiels)), dgeotiels)

    if variogram_range_ratio
        range_ratio = Dict();
        for mission in keys(dh1)
            push!(range_ratio, String(mission) => copy(foo))
        end
    end

    for mission in keys(dh1)
        # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><>
        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        #mission = "hugonnet"

        # find valid range for entire mission so that all filled geotiles cover the same date range
        rrange, = Altim.validrange(vec(any(.!isnan.(dh1[mission]), dims=(1, 3))))

        Threads.@threads for geotile in dims(dh1[mission], :geotile)
            show_times ? t1 = time() : nothing
            
            #for geotile in dims(dh1[mission], :geotile)
            # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><>
            #geotile =       "lat[+80+82]lon[+058+060]"
            # geotile = first(dims(dh1[mission], :geotile))

            # geotile = geotiles[findfirst((geotiles.rgi1 .> 0.) .& (geotiles.glacier_frac .> 0.3)),:]

            # geotile = "lat[-72-70]lon[-014-012]"
            # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

            # println("$mission: $geotile")

            k = findfirst(params[mission].geotile .== geotile)
            df = @view params[mission][k, :]

            dh0 = dh1[mission][At(geotile), :, :]
            nobs0 = nobs1[mission][At(geotile), :, :]
            df.nobs_raw = sum(nobs0)
            df.nbins_raw = sum(nobs0 .> 0)

            ###################################### FILTER 1 ################################
            valid1 = .!isnan.(dh0) .& (nobs0 .> bincount_min[mission]) .& (abs.(dh0) .< 200)
            ################################################################################

            dh0[.!valid1] .= NaN
            dh1[mission][At(geotile), :, :] = dh0
            nobs0[.!valid1] .= 0
            nobs1[mission][At(geotile), :, :] = nobs0

            # if there are not enough points to fit a model then set all to NaNs
            va = sum(valid1)
            if va <= (length(p1) + 2)
                dh1[mission][At(geotile), :, :] .= NaN
                nobs1[mission][At(geotile), :, :] .= 0
                continue
            end

            # determine valid range of data
            (_, crange) = Altim.validrange(valid1)
            valid0 = valid1[rrange, crange]

            # if there are not enough points to fit a model then set all to NaNs
            va = sum(valid0)
            if va <= (length(p1) + 2)
                dh1[mission][At(geotile), :, :] .= NaN
                nobs1[mission][At(geotile), :, :] .= 0
                continue
            end

            dh0 = dh0[rrange, crange];
            nobs0 = nobs0[rrange, crange];
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
            # filter model1_madnorm_max sigma outliers
            valid0[valid0] = Altim.madnorm(dh0_anom) .<= model1_madnorm_max
            vb = sum(valid0)
            df.nbins_filt1 = vb

            # if there are not enough points to fit a model then set all to NaNs
            if vb <= (length(p1) + 2)
                dh1[mission][At(geotile), :, :] .= NaN
                nobs1[mission][At(geotile), :, :] .= 0
                continue

            end

            if vb < va
                dh0[.!valid0] .= NaN
                nobs0[.!valid0] .= 0
                fit1 = curve_fit(model1, hcat(t0[valid0], h0[valid0]), dh0[valid0], nobs0[valid0], p1; lower=lb1, upper=ub1)

                dh0_mdl = model1(hcat(t0[valid0], h0[valid0]), fit1.param)
                dh0_anom = dh0[valid0] .- dh0_mdl
            end

            ####################################################################################
            if variogram_range_ratio
                coord = [(a,b) for (a,b) in zip(t0[valid0], h0[valid0])]
                table = (; Z=vec(dh0_anom))
                geotable = georef(table, coord)

                try
                    vg1 = GeoStats.DirectionalVariogram((1.0, 0.0), geotable, :Z)
                    f1 = GeoStats.fit(SphericalVariogram, vg1)
                    vg2 = GeoStats.DirectionalVariogram((0.0, 1.0), geotable, :Z)
                    f2 = GeoStats.fit(SphericalVariogram, vg2)
                    range_ratio[mission][At(geotile)] = round(Int, (f2.ball.radii ./ f1.ball.radii)[1])
                    println("elevation / height range: $(range_ratio[mission][At(geotile)])")

                catch
                end

                continue
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

    if variogram_range_ratio
        return range_ratio
    end
end


"""
    hyps_amplitude_normalize!(dh1, params; mission_reference = "icesat2")

Normalize seasonal amplitude of elevation change data across different missions.

# Arguments
- `dh1`: Dictionary of elevation change DimArrays by mission
- `params`: Dictionary of parameter DataFrames by mission
- `mission_reference`: Reference mission to normalize against (default: "icesat2")

# Description
This function adjusts the seasonal amplitude of elevation change data from different
missions to match a reference mission (default: ICESat-2). For each geotile:
1. Extracts model parameters for both the target and reference missions
2. Calculates the difference in seasonal components between the models
3. Adds this difference to the target mission's data to normalize its seasonal amplitude
4. Skips geotiles where model parameters are not available for either mission

The function modifies `dh1` in place, adjusting only the seasonal component while
preserving other aspects of the elevation change signal.
"""
function hyps_amplitude_normalize!(dh1, params; mission_reference = "icesat2")

    t = Altim.decimalyear.(dims(dh1[first(keys(dh1))], :date))
    t = repeat(t, 1, length(dims(dh1[first(keys(dh1))], :height)))

    h = val(dims(dh1[first(keys(dh1))], :height))'
    h = repeat(h, length(dims(dh1[first(keys(dh1))], :date)), 1)

    for mission in keys(dh1)
        #mission = "hugonnet"

        if mission == mission_reference
            continue
        end

        for geotile in dims(dh1[mission], :geotile)
            #geotile = "lat[-72-70]lon[-014-012]"
            #geotile = first(dims(dh1[mission], :geotile))
            k = findfirst(params[mission].geotile .== geotile)
            df0 = params[mission][k, :]
            dfr = params[mission_reference][k, :]

            if any(isnan.(df0.param_m1)) || any(isnan.(dfr.param_m1))
                continue
            end

            dh0 = dh1[mission][At(geotile), :, :]
            valid = .!isnan.(dh0)
        
            (rrange, crange) = Altim.validrange(valid)

            dh0 = dh0[rrange, crange]
            valid = valid[rrange, crange]
            t0 = t[rrange, crange] .- df0.t0
            h0x = h[rrange, crange] .- df0.h0

            p0 = df0.param_m1
            pr = dfr.param_m1
            pr[1:5] = p0[1:5] # coeffients 1 to 5 are unrelated to seasonal cycle [pr[1:5] are not used but resetting here for clarity]

            delta = model1_seasonal(hcat(t0[:], h0x[:]), pr) .- model1_seasonal(hcat(t0[:], h0x[:]), p0)

            if !any(isnan.(delta)) # there seem to be rare cases where model1_seasonal returns nans.
                dh1[mission][At(geotile), rrange, crange] = dh0[:] .+ delta
            end
        end
    end
end

"""
    hyps_fill_empty!(dh1, params, geotiles; mask = :glacier)

Fill empty geotiles with data interpolated from nearby geotiles.

# Arguments
- `dh1`: Dictionary of elevation change DimArrays by mission
- `params`: Dictionary of parameter DataFrames by mission
- `geotiles`: DataFrame containing geotile information
- `mask`: Symbol specifying which area mask to use (default: `:glacier`)

# Description
For each mission and geotile:
1. Identifies geotiles with no data but containing the specified mask type (e.g., glacier)
2. Finds the nearest 5 geotiles with valid data at overlapping elevation ranges
3. Interpolates missing data using the median values from neighboring geotiles
4. Fills any remaining gaps along elevation profiles using linear interpolation
5. Sets geotiles with no mask area to zero

# Returns
- Modified `dh1` dictionary with filled geotiles
"""
function hyps_fill_empty!(dh1, params, geotiles; mask = :glacier)

    lon = mean.(getindex.(geotiles.extent, :X))
    lat = mean.(getindex.(geotiles.extent, :Y))
    ll = [(a, b) for (a, b) in zip(lon, lat)]

    for mission in keys(dh1)
        # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><
        #mission = "hugonnet"
        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

         # find valid range for entire mission so that all filled geotiles cover the same date range
        rrange, = Altim.validrange(vec(any(.!isnan.(dh1[mission]), dims=(1, 3))))
        dh0_median = getindex.(params[mission][:, :param_m1], 1) .+ params[mission][:, :dh0]
        
        # copy dh as it gets morphed inside of parallel loop
        dh2 = copy(dh1[mission])

        for geotile in (dims(dh1[mission], :geotile))
            #for geotile in dims(dh1[mission], :geotile)
            # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><>
            #geotile = first(dims(dh1[mission], :geotile)[empty_geotiles])
            #geotile = first(dims(dh1[mission], :geotile))

            # geotile = geotiles[findfirst((geotiles.rgi1 .> 0.) .& (geotiles.glacier_frac .> 0.3)),:]
            #geotile = "lat[-80-78]lon[+166+168]"
            # <><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

            # fill with median of nearest neighbors if no data
            k = findfirst(geotiles.id .== geotile)
            has_ice = geotiles[k, "$(mask)_area_km2"] .> 0

            if !any(has_ice)
                #printstyled("$geotile $mask area is zero: setting dh .= 0 \n"; color=:blue) 
                dh0 = @view dh1[mission][At(geotile), rrange, :]
                dh0 .= 0
                continue
            else
                crange, = Altim.validrange(has_ice)
            end

            dh0 = @view dh1[mission][At(geotile), rrange, crange]

            if all(isnan.(dh0))

                # find distance between goetiles
                dist2geotile = haversine.(Ref(ll[k]), ll, 6371000)

                # find closest X valid geotiles
                nnearest = 5;
                
                #sorted index
                idx = sortperm(dist2geotile)

                # has data in desired elevation range
                has_data = vec(any(.!isnan.(dh2[:, rrange, crange]), dims=(2, 3)));
                has_data = has_data[idx]

                # check that there are enough geotiles with overlapping elevation ranges
                nnearest0 = sum(has_data)
                if nnearest0 < nnearest
                    if nnearest0 == 0
                        printstyled("$mission $geotile: no geotiles with overlapping elevation ranges to fill geotile, setting dh .= 0 \n"; color=:red, bold=true)
                        dh0 .= 0;
                        continue
                    else
                        printstyled("$mission $geotile: less than $nnearest geotiles with overlapping elevation ranges to fill geotile, using closest $nnearest0 geotiles instead\n"; color=:light_red, bold=false)
                        nnearest = nnearest0
                    end
                end

                idx = idx[has_data][1:nnearest]
                f1 = dh2[idx, rrange, crange]

                # remove mean offset to normalize between regions
                dh0_median0 = dh0_median[idx]
                for g in eachindex(dh0_median0)
                    f1[g, :, :] .-= dh0_median0[g]
                end

                for i in 1:length(dims(dh0, :date))
                    for j in 1:length(dims(dh0, :height))
                        f2 = vec(f1[:,i, j])
                        valid00 = .!isnan.(f2)
                        if any(valid00)
                            dh0[i, j] = median(f2[valid00])
                        end
                    end
                end

                # if there are any gaps along an elevation profile then they need to be filled
                (rrange0, crange0) = Altim.validrange(.!isnan.(dh0))
                if any(isnan.(dh0[rrange0, crange0]))
                    x = 1:size(dh0,2)
                    for i in rrange0
                        y = vec(dh0[i,:]);
                        valid = .!isnan.(y)
                        if any(.!valid)
                            itp = DataInterpolations.LinearInterpolation(y[valid], x[valid]; extrapolation=ExtrapolationType.Linear)
                            y[.!valid] = itp(x[.!valid])
                            dh0[i,:] = y;
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
    hyps_fill_updown!(dh1, geotiles; mask = :glacier) -> dh1

Fill missing elevation data at the lowest and highest elevations within each geotile.

# Arguments
- `dh1`: Dictionary of elevation change DimArrays by mission
- `geotiles`: DataFrame containing geotile information
- `mask`: Symbol specifying which area mask to use (default: `:glacier`)

# Returns
- `dh1`: The modified input dictionary with filled elevation data

# Description
For each geotile with valid data, this function extrapolates elevation change values
to the lowest and highest elevation bands by using the first and last valid values
in each time series. This prevents gaps at the extremes of the elevation profile.
The function operates in-place, modifying the input arrays directly.
"""
function hyps_fill_updown!(dh1, geotiles; mask = :glacier)

    for mission in keys(dh1)
        # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><
        # mission = "hugonnet"
        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

        rrange, = Altim.validrange(vec(any(.!isnan.(dh1[mission]), dims=(:geotile, :height))))
        dgeotile = dims(dh1[mission], :geotile);

        Threads.@threads for geotile in dgeotile

            # fill with median of nearest neighbors if not data
            k = findfirst(geotiles.id .== geotile)
            if isnothing(mask)
                valid = geotiles[k, "area_km2"] .> 0
            else
                valid = geotiles[k, "$(mask)_area_km2"] .> 0
            end

            if !any(valid)
                continue
            else
                crange, = Altim.validrange(valid)
            end
            
            dh0 = @view dh1[mission][At(geotile), rrange, crange]

            # fill lowest hand highest elevations with first and last values
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
    hyps_geotile_aggrigate(dv, geotiles, reg; fun=sum::Function) -> DimArray

Aggregate elevation change data from geotiles to regional (RGI) level.

# Arguments
- `dv`: DimArray with dimensions (mission, geotile, date)
- `geotiles`: DataFrame containing geotile information with RGI region columns
- `reg`: Vector of RGI region column names
- `fun`: Aggregation function to apply across geotiles (default: sum)

# Returns
- DimArray with dimensions (mission, rgi, date) containing aggregated values

# Description
For each mission and RGI region, identifies relevant geotiles and applies the 
specified aggregation function to combine their data into regional values.
"""
function hyps_geotile_aggrigate(dv, geotiles, reg; fun=sum::Function)
    dmission, dgeotile, ddate = dims(dv)
    drgi = Dim{:rgi}(reg)

    dv_reg = DimArray(fill(NaN, length(dmission), length(reg), length(ddate)), (dmission, drgi, ddate))
    
    for mission in dmission
        for rgi in drgi
            rgi_ind = geotiles[:, rgi] .> 0
            geotile_ids = geotiles[rgi_ind, :].id
            dv_reg[At(mission), At(rgi), :] = vec(fun(dv[At(mission), At(geotile_ids), :], dims=1))
        end
    end

    return dv_reg
end

"""
    read_zemp2019(; datadir=setpaths().zemp_2019) -> NamedTuple

Read and process glacier mass balance data from Zemp et al. 2019.

# Arguments
- `datadir`: Directory containing Zemp 2019 data files (default: from setpaths())

# Returns
- NamedTuple containing:
  - `dm_gt`: DimArray of cumulative mass change in Gt by RGI region and date
  - `err_gt`: DimArray of uncertainty values in Gt by RGI region and date
  - `all`: Dictionary of raw data by RGI region

# Description
Reads CSV files containing regional glacier mass balance data from Zemp et al. 2019,
extracts RGI region identifiers from filenames, and organizes the data into DimArrays
with dimensions for RGI regions and dates (1950-2016).
"""
function read_zemp2019(;datadir= setpaths().zemp_2019)
    # Read Zemp 2019 data
    fn_startswith = "Zemp_etal_results_region_"
    fn_endswith = ".csv"
    files = allfiles(datadir; fn_startswith, fn_endswith)

    all = Dict()
    for file in files
        # pull RGI number form file name
        foo = file[end-10:(end-7)]
        s = findfirst('_', foo)
        e = findlast('_', foo)
        rgi = "rgi$(foo[(s+1):(e-1)])"
        data = DataFrame((CSV.File(file; header=28, skipto=29, stripwhitespace=true)))
        push!(all, rgi => data)
    end

    date = DateTime.(1950:2016)
    rgi = "rgi" .* string.(1:19)

    drgi = Dim{:rgi}(rgi)
    ddate = Dim{:date}(date)

    dm_gt = DimArray(fill(NaN, (length(drgi), length(ddate))), (drgi, ddate))
    err_gt = DimArray(fill(NaN, (length(drgi), length(ddate))), (drgi, ddate))

    for r in rgi
        d = DateTime.(all[r].Year)
        dm_gt[At(r),At(d)] = cumsum(all[r].INT_Gt)
        err_gt[At(r),At(d)] = all[r].sig_Total_Gt
    end

    return (; dm_gt, err_gt, all)
end


"""
    read_marzeion2020(; datadir=setpaths().marzeion_2020) -> NamedTuple

Read and process glacier mass change data from Marzeion et al. 2020.

# Arguments
- `datadir`: Directory containing Marzeion 2020 data files (default: from setpaths())

# Returns
- NamedTuple containing:
  - `dm_gt`: DimArray of mass change in Gt with dimensions for RGI regions, dates, 
             climate models, glacier models, and scenarios

# Description
Reads NetCDF data containing global glacier mass change projections from Marzeion et al. 2020,
and organizes the data into a multidimensional DimArray with appropriate dimensions.
"""
function read_marzeion2020(;datadir= setpaths().marzeion_2020)

    foo = Dataset(datadir)
    rgi = "rgi" .* string.(Int.(foo["Region"]))
    date = DateTime.(foo["Time"])
    climate_model= getindex.(collect(foo["Climate_Model"].attrib), 2)
    glacier_model= getindex.(collect(foo["Glacier_Model"].attrib), 2)
    scenario = getindex.(collect(foo["Scenario"].attrib), 2)

    drgi = Dim{:rgi}(rgi)
    ddate = Dim{:date}(date)
    dclimate_model = Dim{:climate_model}(climate_model)
    dglacier_model = Dim{:glacier_model}(glacier_model)
    dscenario = Dim{:scenario}(scenario)

    foo = collect(foo["Mass"].var);
    dm_gt = DimArray(foo, (drgi, ddate, dclimate_model, dglacier_model, dscenario))

    return (;dm_gt)
end

"""
    read_marzeion2012(; datadir=setpaths().marzeion_2012) -> NamedTuple

Read and process glacier mass change data from Marzeion et al. 2012.

# Arguments
- `datadir`: Directory containing Marzeion 2012 data files (default: from setpaths())

# Returns
- NamedTuple containing:
  - `dm_gt`: DimArray of mass change in Gt with dimensions for RGI regions, dates, 
             climate models, and scenarios
  - `err_gt`: DimArray of mass change errors in Gt with the same dimensions

# Description
Reads MATLAB data containing global glacier mass change projections from Marzeion et al. 2012,
and organizes the data into multidimensional DimArrays with appropriate dimensions.
"""
function read_marzeion2012(; datadir=setpaths().marzeion_2012)
    foo = matread(datadir)

    scenario = collect(keys(foo))

    # for s in scenario
    s = first(scenario)
    date = vec(DateTime.(foo[s]["yr"]))
    climate_model = vec(foo[s]["model"])
    rgi = vec("rgi" .* string.(Int.(foo[s]["rgiid"])))
    gt2sle = foo[s]["gt2sle"] .* -1

    drgi = Dim{:rgi}(rgi)
    ddate = Dim{:date}(date)
    dclimate_model = Dim{:climate_model}(climate_model)
    dscenario = Dim{:scenario}(scenario)


    dm_gt = fill(NaN, (length(rgi), length(date), length(climate_model), length(scenario)))
    err_gt = fill(NaN, (length(rgi), length(date), length(climate_model), length(scenario)))

    dm_gt = DimArray(dm_gt, (drgi, ddate, dclimate_model, dscenario))
    err_gt = DimArray(err_gt, (drgi, ddate, dclimate_model, dscenario))

    for (i, s) in enumerate(scenario)
        climate_model0 = vec(foo[s]["model"])
        foo1 = foo[s]["mb"] ./ gt2sle
        foo1 = permutedims(foo1, (2, 1, 3))
        dm_gt[:, :, At(climate_model0), i] = foo1

        foo1 = foo[s]["mb_err"] ./ gt2sle
        foo1 = permutedims(foo1, (2, 1, 3))
        err_gt[:, :, At(climate_model0), i] = foo1
    end

    foo = (;dm_gt, err_gt)

    return foo
end



"""
    read_hock2019(; datadir=setpaths().hock_2019) -> NamedTuple

Load and organize glacier mass change projections from Hock et al. 2019.

# Arguments
- `datadir`: Path to the Hock 2019 dataset file (default: from setpaths())

# Returns
- NamedTuple containing `dm_gt`: DimArray of glacier volume change with dimensions for 
  RGI regions, dates, climate models, glacier models, and scenarios

# Description
Reads NetCDF data containing global glacier mass change projections from Hock et al. 2019,
and organizes the data into a multidimensional DimArray with appropriate dimensions.
The function handles the complex combination of scenarios, glacier models, and climate models.
"""
function read_hock2019(; datadir=setpaths().hock_2019)
    #datadir= setpaths().hock_2019
    foo = Dataset(datadir)
    rgi = "rgi" .* string.(Int.(foo["region"]))
    date = DateTime.(foo["time"])

    # from metadata
    scenario = ["A1B", "RCP26", "RCP45", "RCP60", "RCP85"]
    glacier_model = ["BM", "GieOer", "GloGEM", "GloGEMebal", "HYOGA2", "RHB", "VASlangen"]
    climate_model = ["ACCESS1-0", "BCC-CSM1-1", "BNU-ESM", "CCSM4", "CGCM3-1(T63)", "CNRM-CM3", "CNRM-CM5", "CSIRO-Mk3-0", "CSIRO-Mk3-6-0", "CanESM2", "ECHAM5-MPI-OM", "GFDL-CM2-0", "GFDL-CM2-1", "GFDL-CM3", "GFDL-ESM2G", "GFDL-ESM2M", "GISS-E2-R", "HadGEM2-AO", "HadGEM2-CC", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-LR", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC-ESM-CHEM", "MIROC5", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-CGCM3", "NorESM1-M", "NorESM1-ME", "PCM", "UKMO-HadCM3"]

    drgi = Dim{:rgi}(rgi)
    ddate = Dim{:date}(date)
    dclimate_model = Dim{:climate_model}(climate_model)
    dglacier_model = Dim{:glacier_model}(glacier_model)
    dscenario = Dim{:scenario}(scenario)

    coord = (drgi, ddate, dclimate_model, dglacier_model, dscenario)
    dm_gt = DimArray(fill(NaN, length.(coord)), coord)

    for i in 1:length(dclimate_model)
        for j in 1:length(dglacier_model)
            for k in 1:length(dscenario)
                idx = (foo["scenario"] .== k) .& (foo["glaciermodel"] .== j) .& (foo["forcingmodel"] .== i)
                f = foo["volume"].var[:, :, findall(idx)]
                if isempty(f)
                    continue
                else
                    dm_gt[:, :, i, j, k][:] = f
                end
            end
        end
    end

    return (; dm_gt)
end

"""
    plot_height_time(dh1; geotile, fig_suffix, fig_folder, figure_suffix, mask=:glacier, mission, showplots=false)

Create a visualization of elevation change data across height and time dimensions.

# Arguments
- `dh1`: DimArray containing elevation change data
- `geotile`: Geotile information containing ID and area data
- `fig_suffix`: Suffix for the figure filename
- `fig_folder`: Directory where figures will be saved
- `figure_suffix`: Additional suffix for the figure filename
- `mask`: Type of area mask to use (default: `:glacier`)
- `mission`: Mission identifier for the data source
- `showplots`: Whether to display plots in addition to saving them (default: false)

# Description
Generates a two-panel figure with:
1. A heatmap of elevation change data across time (x-axis) and elevation (y-axis)
2. A bar chart showing area distribution by elevation
Both panels share the same x-axis scale for direct comparison.
"""
function plot_height_time(dh1; geotile, fig_suffix, fig_folder, figure_suffix, mask=:glacier, mission, showplots=false)
    height = dims(dh1, :height)
    valid = geotile["$(mask)_area_km2"] .> 0
    crange, = Altim.validrange(geotile["$(mask)_area_km2"] .> 0)

    # find valid range for entire mission so that all filled geotiles cover the same date range
    rrange, = Altim.validrange(vec(any(.!isnan.(dh1), dims=(1, 3))))

    # i need to plot twice, once to get colorbar and once to align x-axis.. makie does not yet support catigorical axis.. but coming soon
    color = Plots.cgrad(:balance, rev=true);
    clim = (-10, 10)

    fname = joinpath(fig_folder, "height_time_colormap.png")
    f = Plots.heatmap(rand(10, 10); clim, color, rightmargin=5Plots.mm)

    savefig(f, fname)

    for legend = [false]
        f = Plots.plot(layout=grid(2, 1, heights=[0.8, 0.2]), link=:x, framestyle=:box, legend=legend)

        Plots.heatmap!(dh1[At(geotile.id), rrange, crange]; subplot=1, legend=false, framestyle=:box, clim, color)
        Plots.xlabel!("")
        Plots.ylabel!("")

        Plots.bar!(val(height[crange]), geotile["$(mask)_area_km2"][crange], subplot=2, legend=false, framestyle=:box, link=:x)
        Plots.ylabel!("area [km²]", subplot=2)
        Plots.xlabel!("height [m]")

        if legend
            fname = joinpath(fig_folder, "$(figure_suffix)_$(mission)_$(geotile.id)_binned_$(fig_suffix)_colorbar.png")
        else
            fname = joinpath(fig_folder, "$(figure_suffix)_$(mission)_$(geotile.id)_binned_$(fig_suffix).png")
        end
        savefig(f, fname)
        showplots && display(f)
    end
end

"""
    plot_height_time(dh1::Dict; geotile, fig_suffix, fig_folder, figure_suffix, mask=:glacier, showplots=false)

Generate height-time plots for all missions in the provided dictionary.

# Arguments
- `dh1`: Dictionary of elevation change DimArrays by mission
- `geotile`: Geotile information containing ID and area data
- `fig_suffix`: Suffix for the figure filename
- `fig_folder`: Directory where figures will be saved
- `figure_suffix`: Additional suffix for the figure filename
- `mask`: Type of area mask to use (default: `:glacier`)
- `showplots`: Whether to display plots in addition to saving them (default: false)

# Description
Iterates through all missions in the dictionary and calls the single-mission version
of plot_height_time for each one, generating separate figures for each mission.
"""
function plot_height_time(dh1::Dict; geotile, fig_suffix, fig_folder, figure_suffix, mask=:glacier, showplots=false)
    for mission in keys(dh1)
        plot_height_time(dh1[mission]; geotile, fig_suffix, fig_folder, figure_suffix, mask, mission, showplots)
    end
end

"""
    read_ipccar6(; datadir=setpaths()[:ipcc_ar6], start2007=false)

Read IPCC AR6 glacier mass change data from CSV files.

# Arguments
- `datadir`: Directory containing IPCC AR6 data files (default: from setpaths())
- `start2007`: Whether to read data starting from 2007 (default: false)

# Returns
A NamedTuple containing:
- `dm_gt`: DimArray of glacier mass change values in Gt by RGI region, date, and scenario
- `err_gt`: DimArray of corresponding error values

# Description
Reads glacier mass change projections from IPCC AR6 Figure 9.21 data files.
The function processes both the main data files and corresponding error files,
organizing them into multidimensional arrays indexed by RGI region, date, and scenario.

# rsync -r devon:/mnt/bylot-r3/data/binned/2deg/figures /Users/gardnera/Research/20_01_GlobalGlacierChange/version\ 2/
"""
function read_ipccar6(;datadir=setpaths()[:ipcc_ar6], start2007=false)
    # Data from IPCC AR6 Figure 9.21
    fn_endswith = ".csv"
    files = allfiles(datadir; fn_endswith)

    f2007 = contains.(files, Ref("2007"))
    if start2007
        files = files[f2007]
    else
        files = files[.!f2007]
    end

    err_files = contains.(files, Ref("error"));
    files_err = files[err_files]
    files = files[.!err_files]
    all = Dict()

    if start2007
        scenarios = [f[end-13:end-9] for f in files]
    else
        scenarios = [f[end-8:end-4] for f in files]
    end

    dscenario = Dim{:scenario}(scenarios)

    data = DataFrame((CSV.File(files[1]; header=1, skipto=2, stripwhitespace=true)))
    ridx = contains.(names(data), Ref("RGI"))

    rgi = "rgi" .* [n[5:end] for n in names(data)[ridx]]
    drgi = Dim{:rgi}(rgi)

    date = DateTime.(data[:, 1])
    ddate = Dim{:date}(date)

    dm_gt = DimArray(fill(NaN, length(drgi), length(ddate), length(dscenario)), (drgi, ddate, dscenario))
    err_gt = copy(dm_gt)

    for k in eachindex(files)
        # pull RGI number form file name
        data = DataFrame((CSV.File(files[k]; header=1, skipto=2, stripwhitespace=true)))
        for j = 1:19
            dm_gt[j,:,k] = data[:,j+1]
        end

        data = DataFrame((CSV.File(files_err[k]; header=1, skipto=2, stripwhitespace=true)))
        for j = 1:19
            err_gt[j, :, k] = data[:, j+1]
        end
    end

    foo = (; dm_gt, err_gt)

    return foo
end

"""
    geotile_binarea!(geotile, ras, feature, bin_edges; invert=false, excludefeature=nothing, var_name)

Calculate area distribution by elevation bins for a specific feature within a geotile.

# Arguments
- `geotile`: Geotile object to be modified with binned area data
- `ras`: Raster containing elevation data
- `feature`: Polygon feature to calculate area for
- `bin_edges`: Elevation bin edges for area calculation
- `invert`: If true, calculate area outside the feature instead (default: false)
- `excludefeature`: Optional feature to exclude from area calculation (default: nothing)
- `var_name`: Name of the variable in geotile to store the binned area results

# Returns
- Modified geotile with binned area data stored in the specified variable

# Description
This function calculates the area distribution by elevation bins for a specific feature
within a geotile. It crops the raster to the geotile extent, creates a mask for the feature,
calculates the area per cell accounting for latitude-dependent cell size, and bins the
area by elevation. The results are stored in the geotile object under the specified variable name.
"""
function geotile_binarea!(geotile, ras, feature, bin_edges; invert=false, excludefeature=nothing, var_name)
    t1 = time()
    ras0 = read(Rasters.crop(ras, to=geotile.extent));

    # using rasterize is 15 times faster than reading its_live masks !!!!
    mask0 = Rasters.rasterize!(count, zeros(ras0.dims), feature; threaded=false, shape=:polygon, progress=false, verbose=false) .> 0;

    if invert
        mask0 = .!mask0
    end

    if .!isnothing(excludefeature)
        excludemask = Rasters.rasterize!(count, zeros(ras0.dims), excludefeature; threaded=false, shape=:polygon, progress=false, verbose=false) .> 0
        mask0 = mask0 .& .!excludemask
    end

    # calculate area per cell
    lon = lookup(ras0, X)
    lat = lookup(ras0, Y)
    d = Altim.meters2lonlat_distance.(Ref(1), lat)
    a = abs.((1 ./ getindex.(d, 2) * (lat[2] .- lat[1])) .* (1 / d[1][1] * (lon[2] - lon[1])))
    area_m2 = repeat(a', outer = [length(lon), 1])


    feature_area_km2 = mask0 .* area_m2 / (1000 * 1000)
    
    foo = geotile[var_name]

    df = DataFrame(ras = ras0[mask0], feature_area_km2 = feature_area_km2[mask0])

    dfbs = BinStatistics.binstats(df, :ras, bin_edges, :feature_area_km2, col_function = sum, missing_bins = true)

    valid = .!(ismissing.(dfbs.feature_area_km2_sum))
    foo[valid] .= dfbs.feature_area_km2_sum[valid]

    t = round(time() - t1)
    printstyled("    -> $(geotile.id) hypsometry calculated: $(t)s\n"; color=:light_black)
    
    return geotile
end

"""
    geotiles_mask_hyps(surface_mask, geotile_width) -> DataFrame

Load and prepare geotiles with hypsometry data for a specific surface mask.

# Arguments
- `surface_mask`: String identifier for the surface mask type (e.g., "glacier", "land")
- `geotile_width`: Width of geotiles in degrees

# Returns
- DataFrame containing geotile data with hypsometry information and geometry

# Description
Loads geotile hypsometry data from an Arrow file, processes the extent information,
and ensures geometry data is properly included by joining with projected geotiles
if necessary.
"""
function geotiles_mask_hyps(surface_mask, geotile_width)

    binned_folder = analysis_paths(; geotile_width).binned
    gt_file = joinpath(binned_folder, "geotile_$(surface_mask)_hyps.arrow");
    geotiles = DataFrame(Arrow.Table(gt_file));
    geotiles.extent = Extent.(getindex.(geotiles.extent, 1));
    geotiles = copy(geotiles);


    # As of now I can not figure out how to restore geometry using GI.Polygon
    df2 = Altim.project_geotiles(; geotile_width);
    if !any(Base.contains.(names(geotiles), "geometry"))
        geotiles = innerjoin(geotiles, df2[!, [:id, :geometry]], on = [:id])
    else
        geotiles = innerjoin(geotiles[!, Not(:geometry)], df2[!, [:id, :geometry]], on=[:id])
    end

    # println(gt_file)
    return geotiles
end
"""
    plot_dh(dh_reg, w; title="", xlims=(DateTime(2000), DateTime(2024)), ylims=nothing)

Plot elevation change time series with fitted model parameters.

# Arguments
- `dh_reg`: DimArray containing elevation change data with mission and date dimensions
- `w`: Weights for model fitting
- `title`: Plot title (default: "")
- `xlims`: X-axis limits as tuple of DateTime objects (default: (DateTime(2000), DateTime(2024)))
- `ylims`: Y-axis limits (default: nothing)

# Returns
- `p`: Plots object with elevation change time series
- `fit_param`: DimArray of fitted model parameters for each mission

# Description
Fits a model to elevation change data for each mission and creates a plot showing
time series with labels indicating mean, trend, acceleration, and amplitude values.
"""
function plot_dh(dh_reg, w; title="", xlims=(DateTime(2000), DateTime(2024)), ylims=nothing)
    dmission = dims(dh_reg, :mission)
    dmetric = Dim{:metric}(["mean", "trend", "acceleration", "amplitude", "date_intercept"])
    fit_param = fill(NaN, (dmission, dmetric))
    ddate = dims(dh_reg, :date)
    decdate = Altim.decimalyear.(ddate)
    for mission in dmission
        dh0 = dh_reg[At(mission), :]
        valid = .!isnan.(dh0)
        if any(valid)
            fit_param[At(mission), At("date_intercept")] = round(Int, mean(decdate[valid]))
            dh_fit = curve_fit(model3, decdate[valid] .- fit_param[At(mission), At("date_intercept")], dh0[valid], collect(w[At(mission), :][valid]), p3;)
            fit_param[At(mission), At("mean")] = dh_fit.param[1]
            fit_param[At(mission), At("trend")] = dh_fit.param[2]
            fit_param[At(mission), At("acceleration")] = dh_fit.param[3]
            fit_param[At(mission), At("amplitude")] = dh_fit.param[4]
        end
    end

    label = collect(dmission)
    prec = 2;
    label = replace.(label, "hugonnet" => "aster")
    
    dh_mean = round.(fit_param[:, At("mean")], digits=prec)
    trend = round.(fit_param[:, At("trend")], digits=prec)
    acc = round.(fit_param[:, At("acceleration")], digits=prec)
    amp = abs.(round.(fit_param[:, At("amplitude")], digits=prec))

    for i in eachindex(label)
        label[i] = "$(label[i]): $(dh_mean[i]) m $(trend[i]) m yr⁻¹ $(acc[i]) m yr⁻² $(amp[i]) m"
    end

    label = reverse(permutedims(label))
    p = Plots.plot(
        dh_reg';
        linewidth=2,
        ylabel="height anomaly [m]",
        label,
        legend=:topright,
        margin=5Plots.mm,
        title,
        legendtitle="mission: mean trend acc amp",
        xlims,
        ylims
    )

    Plots.plot!(size=(1000, 700))
    return p, fit_param
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
function binned_filepath(;binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct)
    if curvature_correct
        runid = "$(surface_mask)_dh_$(dem_id)_cc_$(binning_method)_$(project_id)"
    else
        runid = "$(surface_mask)_dh_$(dem_id)_$(binning_method)_$(project_id)"
    end

    binned_file = joinpath(binned_folder, "$(runid).jld2");

    return binned_file
end

"""
    binned_filled_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)

Generate filepath for filled binned elevation change data and corresponding figure suffix.

# Arguments
- `binned_folder`: Path to the binned folder
- `surface_mask`: Type of surface mask applied (e.g., "ice", "land")
- `dem_id`: Identifier for the DEM used
- `binning_method`: Method used for binning data
- `project_id`: Project identifier
- `curvature_correct`: Boolean indicating if curvature correction was applied
- `amplitude_correct`: Boolean indicating if amplitude correction was applied
- `paramater_set`: Parameter set identifier used for filling

# Returns
- `binned_filled_file`: Full filepath to the filled binned data file
- `figure_suffix`: Suffix string for related figure filenames
"""
function binned_filled_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)

    if curvature_correct
        runid = "$(surface_mask)_dh_$(dem_id)_cc_$(binning_method)_$(project_id)"
    else
        runid = "$(surface_mask)_dh_$(dem_id)_$(binning_method)_$(project_id)"
    end
    
    if amplitude_correct
        binned_filled_file = joinpath(binned_folder, "$(runid)_filled_ac_p$(paramater_set).jld2")
    else
        binned_filled_file = joinpath(binned_folder, "$(runid)_filled_p$(paramater_set).jld2")
    end

    figure_suffix = splitpath(binned_filled_file)
    figure_suffix = figure_suffix[end];
    figure_suffix = replace(figure_suffix, ".jld2"=>"")
    figure_suffix = replace(figure_suffix, "dh" => "dm")

    return binned_filled_file, figure_suffix 
end

"""
    binned_aligned_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)

Generate filepath for aligned binned elevation change data.

# Arguments
- `binned_folder`: Path to the binned folder
- `surface_mask`: Type of surface mask applied (e.g., "ice", "land")
- `dem_id`: Identifier for the DEM used
- `binning_method`: Method used for binning data
- `project_id`: Project identifier
- `curvature_correct`: Boolean indicating if curvature correction was applied
- `amplitude_correct`: Boolean indicating if amplitude correction was applied
- `paramater_set`: Parameter set identifier used for filling

# Returns
- `binned_aligned_file`: Full filepath to the aligned binned data file
"""
function binned_aligned_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)
    binned_filled_file, _ = Altim.binned_filled_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)
    binned_aligned_file = replace(binned_filled_file, ".jld2" => "_aligned.jld2")
    return binned_aligned_file
end

"""
    binned_synthesized_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)

Generate filepath for synthesized binned elevation change data.

# Arguments
- `binned_folder`: Path to the binned folder
- `surface_mask`: Type of surface mask applied (e.g., "ice", "land")
- `dem_id`: Identifier for the DEM used
- `binning_method`: Method used for binning data
- `project_id`: Project identifier
- `curvature_correct`: Boolean indicating if curvature correction was applied
- `amplitude_correct`: Boolean indicating if amplitude correction was applied
- `paramater_set`: Parameter set identifier used for filling

# Returns
- `binned_synthesized_file`: Full filepath to the synthesized binned data file
"""
function binned_synthesized_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)
    binned_filled_file, _ = Altim.binned_filled_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)
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

    out = (; surface_mask, dem, curvature_correct, binning_method, project_id, amplitude_correct, fill_param)
    return out
end


"""
    geotile_bin2d(
        df; 
        var2bin="dh",
        dims_edges=("decyear" => 1990:(30/365):2026, "height_reference" => 0.:100.:10000.),
        binfunction::T = Altim.binningfun_define(binning_method)
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
    binfunction::T = Altim.binningfun_define(binning_method)
    ) where {T <: Function}

    # bin data by date and elevation
    minmax_date = extrema(df[:,dims_edges[1][1]])
    minmax_height = extrema(df[:, dims_edges[2][1]])

    Δd = dims_edges[1][2][2] - dims_edges[1][2][1]
    Δh = dims_edges[2][2][2] - dims_edges[2][2][1]

    date_ind = (dims_edges[1][2] .>= minmax_date[1] - Δd) .&
            (dims_edges[1][2] .<= (minmax_date[2] + Δd))

    date_ind_center = findall(date_ind)[1:end-1];

    height_ind = (dims_edges[2][2] .>= minmax_height[1] - Δh) .&
                (dims_edges[2][2] .<= (minmax_height[2] + Δh))

    height_ind_center = findall(height_ind)[1:end-1];

    nobs0 = nothing
    var0 = nothing

    # check bounds: binstats will throw an error if no data is passed to median()
    if !any(date_ind)
        return var0, nobs0
    end

    if !Altim.vector_overlap(df[!, dims_edges[1][1]], dims_edges[1][2][date_ind]) ||
    !Altim.vector_overlap(df[!, dims_edges[2][1]], dims_edges[2][2][height_ind])
        
    return var0, nobs0
    end

    df = Altim.binstats(df, [getindex.(dims_edges,1)...], [getindex.(dims_edges,2)...],
        var2bin; col_function=[binfunction], missing_bins=true)

    gdf = DataFrames.groupby(df, dims_edges[1][1])

    dd1 = Dim{Symbol(dims_edges[1][1])}(sort((dims_edges[1][2][1:end-1] .+ dims_edges[1][2][2:end])./2))
    dd2 = Dim{Symbol(dims_edges[2][1])}(sort((dims_edges[2][2][1:end-1] .+ dims_edges[2][2][2:end])./2))

    nobs0 = fill(0, (dd1, dd2))
    var0 = fill(NaN, (dd1, dd2)) # this should really match the type of the input data ... but hey... this is easier right now
        
    # get index into sorted array
    p = sortperm(gdf[1][:,dims_edges[2][1]])

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
    glacier_discharge(; datadir=Altim.pathlocal[:data_dir]) -> DataFrame

Load and combine glacier discharge data from multiple sources.

# Arguments
- `datadir`: Base directory containing glacier data files (default: Altim.pathlocal[:data_dir])

# Returns
- DataFrame containing glacier discharge information with columns:
  - `latitude`: Glacier center latitude
  - `longitude`: Glacier center longitude
  - `discharge_gtyr`: Discharge rate in gigatons per year
  - `discharge_err_gtyr`: Uncertainty in discharge rate
  - `frontal_ablation_gtyr`: Frontal ablation rate in gigatons per year

# Description
Combines discharge data from Kochtitzky (2022) for the Northern Hemisphere and 
Fuerst (2023) for Patagonia. Matches Patagonian glaciers with RGI database entries
to obtain coordinates, and standardizes data format across sources.
"""
function glacier_discharge(; datadir=Altim.pathlocal[:data_dir])
    # Kochtitzky NH discharge and terminus retreate 
    nh_discharge = joinpath(datadir, "GlacierOutlines/GlacierDischarge/Kochtitzky2022/41467_2022_33231_MOESM4_ESM.csv")
    nothern_hemisphere = CSV.read(nh_discharge, DataFrame; header=14, skipto=16)

    fn =  joinpath(datadir, "GlacierOutlines/rgi60/17_rgi60_SouthernAndes/17_rgi60_SouthernAndes.shp")
    df = DataFrame(Shapefile.Table(fn))

    fn =  joinpath(datadir, "GlacierOutlines/GlacierDischarge/Fuerst2023/fuerst_2023_npi_comparison_ice_discharge_v1.0.0.txt")
    patagonia = CSV.read(fn, DataFrame; header=25, skipto=27)
    fn =  joinpath(datadir, "GlacierOutlines/GlacierDischarge/Fuerst2023/fuerst_2023_spi_comparison_ice_discharge_v1.0.0.txt")
    patagonia = vcat(patagonia, CSV.read(fn, DataFrame; header=25, skipto=27))

    # find matching glaciers in RGI 
    patagonia[!, :"Latitude"]  .= NaN;
    patagonia[!, :"Longitude"]  .= NaN;
    df.Name[ismissing.(df.Name)] .= "junk"
    patagonia.Names = replace.(patagonia.Names, "Témpano" => "Tempano")
    patagonia.Names[Base.contains.(patagonia.Names, "Grey")] .= "Grey + Dickson"
    patagonia.Names[Base.contains.(patagonia.Names, "Upsala")] .= "Upsala + Cono"
    patagonia.Names[Base.contains.(patagonia.Names, "O'Higgins")] .= "OHiggins"

    for r in eachrow(patagonia)
        index = findfirst((r.Names .== df.Name) .| (r.Names .==  df.RGIId))
        if isnothing(index)
            println("cold not find match for: $(r.Names)")
            
            continue
        end
        r.Latitude = df.CenLat[index]
        r.Longitude = df.CenLon[index]
    end

    df = DataFrame()
    df[!, :latitude] = vcat(nothern_hemisphere.lat, patagonia.Latitude)
    df[!, :longitude] = vcat(nothern_hemisphere.lon, patagonia.Longitude)
    df[!, :discharge_gtyr] = vcat(nothern_hemisphere[:,"2010_2020_mean_discharge_gt_per_year"], (patagonia[:,"Calving FG"] .* Altim.δice/1000))
    df[!, :discharge_err_gtyr] = vcat(nothern_hemisphere[:,"2010_2020_mean_flux_err_gt"], (patagonia[:,"Unc. Calving FG"].* Altim.δice/1000))
    df[!, :frontal_ablation_gtyr] = vcat(nothern_hemisphere[:,"Frontal_ablation_2010_to_2020_gt_per_yr_mean"], (patagonia[:,"Frontal ablation Minowa (2000-2019)"].* Altim.δice/1000))
    df[!, :frontal_ablation_gtyr] = vcat(nothern_hemisphere[:,"Frontal_ablation_2010_to_2020_gt_per_yr_mean"], (patagonia[:,"Unc. Frontal ablation Minowa (2000-2019)"].* Altim.δice/1000))

    return df
end


"""
    discharge2smb(glaciers; discharge2smb_max_latitude=-60, discharge2smb_equilibrium_period=(Date(1979), Date(2000)))

Calculate discharge from surface mass balance (SMB) trends for glaciers below a specified latitude.

# Arguments
- `glaciers`: DataFrame containing glacier data with columns :CenLat, :CenLon, :smb, and :area_km2
- `discharge2smb_max_latitude`: Maximum latitude threshold for calculating discharge (default: -60)
- `discharge2smb_equilibrium_period`: Time period for equilibrium calculation (default: 1979-2000)

# Returns
- `discharge0`: DataFrame with columns :latitude, :longitude, :discharge_gtyr, :discharge_err_gtyr, :frontal_ablation_gtyr

# Description
This function estimates glacier discharge by calculating SMB trends during an equilibrium period.
For glaciers below the specified latitude threshold, it fits a linear trend to the SMB data
and converts this to discharge in gigatons per year based on glacier area.
"""
function discharge2smb(glaciers; discharge2smb_max_latitude=-60, discharge2smb_equilibrium_period=(Date(1979), Date(2000)))
    index_discharge2smb = glaciers.CenLat .< discharge2smb_max_latitude

    discharge0 = DataFrame(latitude=glaciers[index_discharge2smb, :CenLat], longitude=glaciers[index_discharge2smb, :CenLon], discharge_gtyr=NaN, discharge_err_gtyr=NaN, frontal_ablation_gtyr=NaN)

    ddate = dims(glaciers[1, :smb], :date)
    decyear = Altim.decimalyear.(ddate)
    Δdecyear = decyear .- decyear[1]
    index_date = (ddate .>= discharge2smb_equilibrium_period[1]) .& (ddate .<= discharge2smb_equilibrium_period[2])

    for (i, glacier) in enumerate(eachrow(glaciers[index_discharge2smb, :]))
        y = glacier.smb[index_date]
        fit = curve_fit(Altim.offset_trend, Δdecyear[index_date], y .- mean(y), Altim.offset_trend_p)
        discharge0[i, :discharge_gtyr] = fit.param[2] .* sum(glacier.area_km2) / 1000
    end
    return discharge0
end

"""
    dh2dv(dh, geotiles) -> dv

Convert elevation change (dh) to volume change (dv) using area information from geotiles.

# Arguments
- `dh`: DimArray of elevation changes with dimensions (geotile, date, height)
- `geotiles`: DataFrame containing geotile information with columns :id and :area_km2

# Returns
- `dv`: DimArray of volume changes with dimensions (geotile, date)

# Description
This function calculates volume change by multiplying elevation change by the corresponding 
area for each geotile and height bin, then summing across all height bins. The result is 
converted from km³ to Gt by dividing by 1000.
"""
function dh2dv(dh, geotiles)
    dv = fill(NaN, dims(dh[:, :, 1]))
    index_date = vec(any(.!isnan.(dh[1, :, :]), dims=:height))
    
    for geotile in eachrow(geotiles)
        #geotile = eachrow(geotiles)[560]
        index_height = geotile.area_km2 .> 0
        if !any(index_height)
            dv[At(geotile.id), :] .= 0
        else
            area = (geotile.area_km2[index_height] * ones(1, sum(index_date)))'
            dv[At(geotile.id), index_date] = sum((dh[At(geotile.id), index_date, index_height] .* area ./ 1000); dims=:height)
        end
    end

    return dv
end