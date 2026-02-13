"""
    gemb_read2(gemb_file; vars="all", datebin_edges=nothing)

Read GEMB data from a single file with optional date binning.

# Arguments
- `gemb_file`: Path to a GEMB data file
- `vars`: Either "all" to read all variables, or an array of specific variable names
- `datebin_edges`: Optional array of date edges for binning the data

# Returns
- Dictionary containing the requested variables with optional date binning applied.

# Examples
```julia
julia> gemb = gemb_read2(gemb_file; vars=["smb", "runoff", "acc"])
julia> gemb_all = gemb_read2(gemb_file; datebin_edges=date_edges)
```
"""
function gemb_read2(gemb_file; vars="all", datebin_edges=nothing)

    vars0 = ["latitude", "longitude", "date", "smb", "fac", "ec", "acc", "runoff", "melt", "fac_to_depth", "height", "refreeze", "t_air", "rain"]

    # variables that need to be converted from cumulitive outputs to rates

    var_cumulitive = ["smb", "ec", "acc", "runoff", "melt", "refreeze", "rain"]

    varmap = Dict(
        "latitude" => "lat",
        "longitude" => "lon",
        "date" => "time",
        "smb" => "SMB",
        "fac" => "FAC",
        "ec" => "EC",
        "acc" => "Accumulation",
        "runoff" => "Runoff",
        "melt" => "Melt",
        "fac_to_depth" => "FACtoDepth",
        "height" => "H",
        "refreeze" => "Refreeze",
        "t_air" => "Ta",
        "rain" => "Rain"
    )

    gemb = Dict()
    if vars == "all"
        foo = matread(gemb_file)
        for v in varmap
            f = foo[v[2]]
            f[isnan.(f)] .= 0 #set any random nan to zero
            gemb[v[1]] = f
        end
    else
        if length(intersect(vars0, vars)) != length(vars)
            error("one or more vars not recognized")
        end

        gemb = matopen(gemb_file) do file
            foo = MAT.read.(Ref(file), [varmap[v] for v in vars])
            for (i, v) in enumerate(vars)
                f = foo[i]
                f[isnan.(f)] .= 0 #set any random nan to zero
                gemb[v] = f
            end
            return gemb
        end
    end

    # inputs are in cumulative amounts ... e.g. mm per 5 days
    # convert to cumulitive meters
    for k in keys(gemb)
        if any(var_cumulitive .== k)
            gemb[k] = cumsum(gemb[k] ./ 1000, dims=2)
        end
    end

    if .!isnothing(datebin_edges)
        ind = Vector{Union{UnitRange{Int64},Nothing}}()
        t = vec(gemb["date"])
        for i = eachindex(datebin_edges)[1:end-1]
            a = findfirst((t .>= datebin_edges[i]) .& (t .< datebin_edges[i+1]))
            b = findlast((t .>= datebin_edges[i]) .& (t .< datebin_edges[i+1]))

            if isnothing(a)
                push!(ind, nothing)
            else
                push!(ind, a:b)
            end
        end

        for k in keys(gemb)
            if k == "date"
                continue
            end

            if size(gemb[k], 2) == size(gemb["date"], 2)
                foo = fill(NaN, (size(gemb[k], 1), length(ind)))
                for i in eachindex(ind)
                    if !isnothing(ind[i])
                        foo[:, i] = mean(gemb[k][:, ind[i]], dims=2)
                    end
                end
                gemb[k] = foo
            end
        end
        gemb["date"] = reshape((datebin_edges[1:end-1] .+ datebin_edges[2:end]) ./ 2, (1, length(datebin_edges) - 1))
    end

    return gemb
end


"""
    _matrix_shift_ud!(matrix, shift; exclude_zeros_in_extrapolation=false, npts_linear_extrapolation=Inf)

Shifts a matrix up or down by a specified number of rows, with linear extrapolation at boundaries.

# Arguments
- `matrix`: Input matrix to be shifted in-place
- `shift`: Integer number of rows to shift. Positive shifts down, negative shifts up.
- `exclude_zeros_in_extrapolation`: If true, zero values are excluded when calculating extrapolation
- `npts_linear_extrapolation`: Maximum number of points to use for linear extrapolation

# Returns
- The modified input matrix (in-place).

# Examples
```julia
julia> _matrix_shift_ud!(matrix, 2)
julia> # matrix shifted down by 2 rows with extrapolation at top
```
"""
function _matrix_shift_ud!(matrix, shift; exclude_zeros_in_extrapolation=false, npts_linear_extrapolation=Inf)

    n = size(matrix, 2)
    h0 = (1.:n) ./ n .- 0.5
    M = hcat(ones(n), h0)

    # shift with linear extrapolation
    if shift > 0 # shift values down 
        # shift values down
        matrix[:, 1:(end-shift)] = matrix[:, (shift+1):end]

        # extrapolate upper boundary
        M0 = @view M[1:(end-shift), :]

        for i in axes(matrix, 1)
            y0 = @view matrix[i, 1:(end-shift)]
            if exclude_zeros_in_extrapolation
                valid = y0 .!= 0
                if any(valid)
                    if sum(valid) > npts_linear_extrapolation
                        valid = findall(valid)[(end-npts_linear_extrapolation+1):end]
                    end
                    param1 = M0[valid, :] \ y0[valid]
                else
                    matrix[i, (end-shift+1):end] .= 0
                    continue
                end
            else
                param1 = M0 \ y0
            end

            matrix[i, (end-shift+1):end] = param1[1] .+ param1[2] .* h0[(end-shift+1):end]
        end

    elseif shift < 0 # shift values up

        # shift values up
        matrix[:, (1-shift):end] = matrix[:, 1:(end+shift)]

        # extrapolate lower boundary
        M0 = @view M[(1-shift):end, :]

        for i in axes(matrix, 1)
            y0 = @view matrix[i, (1-shift):end]
            if exclude_zeros_in_extrapolation
                valid = y0 .!= 0
                if any(valid)
                    if sum(valid) > npts_linear_extrapolation
                        valid = findall(valid)[1:npts_linear_extrapolation]
                    end
                    param1 = M0[valid, :] \ y0[valid]
                else
                    matrix[i, 1:(-shift)] .= 0
                    continue
                end
            else
                param1 = M0 \ y0
            end

            matrix[i, 1:(-shift)] = param1[1] .+ param1[2] .* h0[1:(-shift)]
        end
    end
    return matrix
end


"""
    gemb_rate_physical_constraints!(gemb)

Apply physical constraints to GEMB rate variables to ensure consistency.

# Arguments
- `gemb`: Dictionary containing GEMB variables

# Returns
The input dictionary with physically consistent values:
- Negative values for accumulation, rain, melt, and refreeze are set to zero
- Refreeze cannot exceed melt
- Runoff equals melt minus refreeze
- SMB equals accumulation minus runoff minus evaporation/condensation

# Examples
```julia
julia> gemb_rate_physical_constraints!(gemb)
julia> # gemb["runoff"], gemb["smb"], etc. updated in-place
```
"""
function gemb_rate_physical_constraints!(gemb)
    acc = gemb["acc"]
    acc[acc.<0] .= 0

    rain = gemb["rain"]
    rain[rain.<0] .= 0

    melt = gemb["melt"]
    melt[melt.<0] .= 0

    refreeze = gemb["refreeze"]
    refreeze[refreeze.<0] .= 0

    # refreeze can not exceed melt
    refreeze[refreeze.>melt] = melt[refreeze.>melt]

    # runoff must equal melt - refreeze [rain is not included in runoff calculation]
    gemb["runoff"] = melt .- refreeze

    # SMB must equal acc - runoff - ec [rain is not included in SMB calculation]
    gemb["smb"] = acc .- gemb["runoff"] .- gemb["ec"]

    return gemb
end



"""
    gemb_bestfit_grouped(dv_altim, smb, fac, discharge, geotiles; single_geotile_test = nothing)

Calibrate the SMB model to grouped geotiles to handle glaciers that cross multiple geotiles.

# Arguments
- `dv_altim`: Volume change from altimetry data
- `smb`: Surface mass balance data
- `fac`: Firn air content data
- `discharge`: Discharge data
- `geotiles`: DataFrame containing geotile information
- `single_geotile_test`: Optional geotile ID to examine model fits for debugging 

# Returns
- DataFrame with calibrated parameters for each geotile: id, extent, pscale, mscale, and rmse

# Examples
```julia
julia> geotiles_fit = gemb_bestfit_grouped(dv_altim, dv_gemb, geotiles0; seasonality_weight=0.85)
julia> pscale = geotiles_fit.pscale
```
"""
function gemb_bestfit_grouped(dv_altim, dv_gemb, geotiles0;
    distance_from_origin_penalty = 20/100,
    mscale_to_pscale_weight = 1,
    seasonality_weight = 85/100,
    calibrate_to=:all,
    is_scaling_factor=Dict("pscale" => true, "mscale" => true)
)

    if mscale_to_pscale_weight > 1 || mscale_to_pscale_weight < 0
        error("mscale_to_pscale_weight must be between 0 and 1")
    end

    if seasonality_weight > 1 || seasonality_weight < 0
        error("seasonality_weight must be between 0 and 1")
    end
    
    # there are issues with calibrating the dv_gemb model to individual geotiles as glaciers
    # can cross multiple geotiles. To tackle this issue one needs to process mutually 
    # exclusive geotile groups

    # loop for each geotile
    geotiles = DataFrame()
    geotiles[!, :id] = geotiles0[:, :id]
    geotiles[!, :group] = geotiles0[:, :group]
    geotiles[!, :rgi] = geotiles0[:, :rgi]
    geotiles[!, :extent] = geotiles0[:, :extent]
    geotiles[!, :area_km2] = geotiles0[:, :area_km2]
    

    geotiles[!, :pscale] .= 1.0
    geotiles[!, :mscale] .= 1.0
  
    dpscale = dims(dv_gemb, :pscale)
    dmscale = dims(dv_gemb, :mscale)
 
    bounds = boxconstraints(lb=[minimum(dpscale.val), minimum(dmscale.val)].+.01, ub=[maximum(dpscale.val), maximum(dmscale.val)].-.01);


    kwargs2 = (seasonality_weight=seasonality_weight, distance_from_origin_penalty=distance_from_origin_penalty, mscale_to_pscale_weight=mscale_to_pscale_weight, calibrate_to=calibrate_to, is_scaling_factor=is_scaling_factor)

    geotile_groups = copy(geotiles.group);
    geotile_ids = copy(geotiles.id)
    groups_unique = unique(geotile_groups)


     # for full grid search
    begin
        step_size = 0.05
        s0 = -1 / minimum(dpscale.val)
        e0 = maximum(dpscale.val)
        pscale0 = s0:step_size:e0
        dpscale_search = vcat(-1 ./ pscale0[pscale0.<-1], pscale0[pscale0.>=1])

        if all(dmscale .> 0)
            s0 = -1 / minimum(dmscale.val)
            e0 = maximum(dmscale.val)
            mscale0 = s0:step_size:e0
            dmscale_search = vcat(-1 ./ mscale0[mscale0.<-1], mscale0[mscale0.>=1])
        else
            dmscale_search = minimum(dmscale.val):((maximum(dmscale.val)-minimum(dmscale.val))/20):maximum(dmscale.val)
        end

        dpscale_search = Dim{:pscale}(dpscale_search)
        dmscale_search = Dim{:mscale}(dmscale_search)
    end
    # loop for each group
    # NOTE: do not do a groupby on geotiles as you need to keep the original geotiles order

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CAUSE OF DAYS LOST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # This had to be moved out of @threads loop as it is not thread safe - your guess as as good as mine [this block only take 1.5s so it's a non-issue]
    dv_altim0 = [Float32.(dropdims(sum(dv_altim[geotile=At(geotile_ids[groups_unique[i] .== geotile_groups])], dims=:geotile), dims=:geotile)) for i in eachindex(groups_unique)]
    dv_gemb0  = [Float32.(dropdims(sum(dv_gemb[geotile=At(geotile_ids[groups_unique[i]  .== geotile_groups])], dims=:geotile), dims=:geotile)) for i in eachindex(groups_unique)]

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    if length(groups_unique) != 1

        pscale0 = zeros(length(groups_unique));
        mscale0 = zeros(length(groups_unique));
 
        Threads.@threads for i in 1:length(groups_unique)

            #(pscale0[i], mscale0[i]) = gemb_altim_cost_group(dpscale_search, dmscale_search, dv_altim0, dv_gemb0, kwargs2)
           
            # Create closure with thread-local data
            # Each thread gets its own copy of the data, ensuring thread safety
            f2 = Base.Fix{2}(Base.Fix{3}(Base.Fix{4}(gemb_altim_cost, kwargs2), dv_gemb0[i]), dv_altim0[i])

            # Seed the task-local RNG with a deterministic seed unique to this group
            # This ensures each thread has independent RNG state, preventing race conditions
            my_options = Options(f_tol=1e-6, x_tol=1e-10, f_calls_limit=3000, store_convergence=false, seed=i)
            result = optimize(f2, bounds, ECA(; η_max=1.0, K=6, options=my_options))
            (pscale0[i], mscale0[i]) = minimizer(result)
        end

        # modify dataframe outside of loop
        for i in eachindex(groups_unique)
            gindex = groups_unique[i] .== geotile_groups;
            geotiles[gindex, :pscale] .= pscale0[i]
            geotiles[gindex, :mscale] .= mscale0[i]
        end

        return geotiles[:, [:id, :extent, :group, :pscale, :mscale, :rgi]]

    else
        i = 1;
        geotiles_in_group = geotile_ids[groups_unique[i] .== geotile_groups]
        f2 = Base.Fix{2}(Base.Fix{3}(Base.Fix{4}(gemb_altim_cost, kwargs2), dv_gemb0[i]), dv_altim0[i])
       

        time_grid_search = @elapsed begin
            cost = fill(NaN, dpscale_search, dmscale_search)
            for pscale in dpscale_search
                for  mscale in dmscale_search
                    cost[pscale=At(pscale), mscale=At(mscale)] = f2([pscale, mscale])
                end
            end
        end

        show_plot = false
        if show_plot
            center_on_dates = DateTime(2010, 1, 1)..DateTime(2015, 1, 1)
            f = Figure();
            ax = Makie.Axis(f[1, 1]; title="geotile: $(geotiles_in_group[1])");
        
            cost0 = Inf;
            time_grid_search = @elapsed begin
                cost = fill(NaN, dpscale_search, dmscale_search)
                for pscale in dpscale_search
                    for  mscale in dmscale_search
                        cost1 = f2([pscale, mscale])

                        if cost0 > cost1
                            cost0 = cost1
                
                            gemb2X = gemb_dv_sample(pscale, mscale, dv_gemb0[i])
                            gemb2X = gemb2X .- mean(gemb2X[date=center_on_dates])
                            lines!(ax, gemb2X; label = "cost: $(round(cost0, digits=1))")
                        end
                    end
                end
            end

            altim2X = dv_altim0[i] .- mean(dv_altim0[i][date = center_on_dates])

            lines!(ax, altim2X; linewidth=4, color=:black, label="altim")
            Legend(f[1, 2], ax, position=:rt);
            display(f)
        end

    
        index_minimum = argmin(cost)
        pscale_grid_search = DimPoints(cost)[index_minimum][1]
        mscale_grid_search = DimPoints(cost)[index_minimum][2]
        println("grid_search: pscale: $(round(pscale_grid_search, digits=2)), mscale: $(round(mscale_grid_search, digits=2)), time: $(round(time_grid_search, digits=2))")

        # do an optimization search
        time_eca_search = @elapsed begin
            result = optimize(f2, bounds, ECA())
        end
        (pscale_eca, mscale_eca) = minimizer(result)
        println("eca_search : pscale: $(round(pscale_eca, digits=2)), mscale: $(round(mscale_eca, digits=2)), time: $(round(time_eca_search, digits=2))")

        #f = Makie.lines(dv_altim0)
        #Makie.lines!(gemb_dv_sample(pscale_best, mscale_best, dv_gemb0))
        #display(f)

        #println("pscale_best: $(round(pscale_best, digits=1)), mscale_best: $(round(mscale_best, digits=1))")
        return cost, geotiles_in_group
    end
end

"""
    add_pscale_classes(var0, dpscale_new; allow_negative=false)

Interpolate and extrapolate values across a new set of precipitation scale factors.

# Arguments
- `var0`: Input dimensional array with dimensions (:date, :height, :pscale, :mscale)
- `dpscale_new`: New precipitation scale factor dimension to interpolate/extrapolate to
- `allow_negative`: Whether to allow negative values in the output (default: false)

# Returns
A new dimensional array with the same structure as `var0` but with the new pscale dimension.

This function performs linear interpolation for pscale values within the range of the original
data, and linear extrapolation for values outside that range. For extrapolation, it uses either
the first/last three points or all available points depending on data availability.
"""
function add_pscale_classes(var0, dpscale_new; allow_negative=false)

    ddate = dims(var0, :date)
    dheight = dims(var0, :height)
    dpscale = dims(var0, :pscale)
    dmscale = dims(var0, :mscale)

    interpolate_index = (dpscale_new .>= minimum(dpscale)) .& (dpscale_new .<= maximum(dpscale))
    extrapolate_lower_index = .!interpolate_index
    extrapolate_lower_index[findfirst(interpolate_index):end] .= false
    extrapolate_upper_index = .!interpolate_index
    extrapolate_upper_index[1:findlast(interpolate_index)] .= false
    extrapolate_index = extrapolate_lower_index .| extrapolate_upper_index
    x = collect(dpscale)
    x_new = collect(dpscale_new)

    # this can be removed after code rerun
    x_sort_perm = sortperm(x)
    x = x[x_sort_perm]

    M = hcat(ones(length(x)), x)

    gemb_new = fill(NaN, (ddate, dheight, dpscale_new, dmscale))

    valid_range, = validrange(.!isnan.(var0[:, 1, 1, 1]))

    Threads.@threads for date in ddate[valid_range]
        #date = ddate[820]

        y_new = fill(NaN, length(x_new))

        for height in dheight
            #height = 1050

            for mscale in dmscale
                #mscale

                y = var0[At(date), At(height), :, At(mscale)]

                # this can be removed after code rerun
                y[:] = y[x_sort_perm]

                if !allow_negative
                    valid_interp = y .>= 0
                end

                if sum(valid_interp) >= 2

                    # Loess interpolation on so few points == linear interpolation but is much slower
                    #model = Loess.loess(x, y, span=0.3)
                    #y_new[interpolate_index] = Loess.predict(model, x_new[interpolate_index])

                    model = DataInterpolations.LinearInterpolation(y, x)
                    y_new[interpolate_index] = model(x_new[interpolate_index])

                    x0 = x[valid_interp]
                    y0 = y[valid_interp]
                    M0 = M[valid_interp, :]

                    if length(x0) > 3
                        lowest_pts = 1:3
                        highest_pts = length(x0)-2:length(x0)

                        param1 = M0[lowest_pts, :] \ y0[lowest_pts]
                        y_new[extrapolate_lower_index] = param1[1] .+ x_new[extrapolate_lower_index] .* param1[2]

                        param2 = M0[highest_pts, :] \ y0[highest_pts]
                        y_new[extrapolate_upper_index] = param2[1] .+ x_new[extrapolate_upper_index] .* param2[2]
                    else
                        param1 = M0 \ y0
                        y_new[extrapolate_index] = param1[1] .+ x_new[extrapolate_index] .* param1[2]
                    end

                    if !allow_negative
                        y_new[y_new.<0] .= 0
                    end

                    # sanity check
                    # p = lines(x_new, y_new);
                    # CM.plot!(x,collect(y))
                else
                    y_new[:] .= 0
                end

                gemb_new[At(date), At(height), :, At(mscale)] = y_new
            end
        end
    end

    return gemb_new
end


"""
    gemb_calibration(path2runs; surface_mask_default="glacier", gemb_run_id, geotile_grouping_min_feature_area_km2=100, single_geotile_test=nothing, force_remake_before)

Calibrate GEMB (Glacier Energy and Mass Balance) model to altimetry data by finding optimal fits for grouped geotiles.

# Arguments
- `path2runs`: Array of file paths to binned synthesized data files
- `surface_mask_default`: Default surface mask type (default: "glacier")
- `gemb_run_id`: GEMB run identifier
- `geotile_grouping_min_feature_area_km2`: Minimum feature area in km² for geotile grouping (default: 100)
- `single_geotile_test`: Single geotile ID for testing (default: nothing)
- `force_remake_before`: Date to force file regeneration (default: nothing)

# Description
Loads GEMB data and altimetry volume change data, then finds optimal model fits for grouped geotiles.
Saves calibration results as Arrow files with geometry information.
"""
function gemb_calibration(
    path2runs,
    dv_gemb;
    geotile_width=2,
    geotile_grouping_min_feature_area_km2=100,
    single_geotile_test=nothing,
    seasonality_weight=95/100,
    distance_from_origin_penalty=2 / 100,
    mscale_to_pscale_weight=1,
    force_remake_before=nothing,
    calibrate_to = :all
)

    # Load GEMB data
    # Note: GEMB volume change is for a single surface mask - this is a hack since making GEMB dv for multiple surface masks is onerous
    # However, since GEMB data is calibrated to altimetry that uses different surface masks, this is mostly a non-issue

    if !isnothing(single_geotile_test)
        # NOTE: don't prematurely trim data to `single_geotile_test` as fit optimization is done by geotile groups
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)], OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!"
    end

    dgeotile = dims(dv_gemb, :geotile)
    dmscale = dims(dv_gemb, :mscale)
    dpscale = dims(dv_gemb, :pscale)
    params = binned_filled_fileparts.(path2runs)
    surface_masks = unique(getindex.(params, :surface_mask))

    # having issues with pscale and mscale seeming somewhat random... I'm thinking that it's realated to something that is no Thread safe
    dgeotile_test = Dim{:geotile}([geotiles_golden_test[1], "lat[+58+60]lon[-138-136]"])
    pscale_range_test = DimStack(DimArray([1.8, 1.5], dgeotile_test; name="min"), DimArray([3.0, 2.5], dgeotile_test; name="max"))

    # Process each surface mask to load aligned geotiles and calculate hypsometry
    geotiles = Dict()
    area_km2 = Dict()

    if all(dmscale .> 0)
        mscale_is_scaling = true;
    else
        mscale_is_scaling = false;
    end

    if all(dpscale .> 0)
        pscale_is_scaling = true;
    else
        pscale_is_scaling = false;
    end

    is_scaling_factor = Dict("pscale" => pscale_is_scaling, "mscale" => mscale_is_scaling)


    for surface_mask in surface_masks # takes 7 seconds

        # Load and align geotiles with the synthesized data dimensions
        geotiles[surface_mask] = _geotile_load_align(;
            surface_mask,
            geotile_width,
            only_geotiles_w_area_gt_0=true,
            geotile_order=collect(dgeotile)
        )

        area_km2[surface_mask] = _geotile_area_km2(; surface_mask, geotile_width)[geotile=At(val(dgeotile))]

        # Add in groups
        geotiles_groups = geotile_grouping(; surface_mask, min_area_km2=geotile_grouping_min_feature_area_km2, geotile_width, force_remake_before)

        # Add groups to geotiles
        _, grp_index = intersectindices(geotiles[surface_mask].id, geotiles_groups.id)
        geotiles[surface_mask][!, :group] = geotiles_groups[grp_index, :group]
    
        if !isnothing(single_geotile_test)
            index_geotile = findfirst(geotiles[surface_mask].id .== single_geotile_test)
            geotiles[surface_mask] = geotiles[surface_mask][(geotiles[surface_mask].group .== geotiles[surface_mask][index_geotile, :group]), :]
        end

        if false
            @warn "testing subset"
            index_test = falses(nrow(geotiles[surface_mask]))
            for geotile_test in dgeotile_test
                index_geotile_test = findfirst(geotiles[surface_mask].id .== geotile_test)
                index_test .|= (geotiles[surface_mask].group .== geotiles[surface_mask].group[index_geotile_test])
            end
        
            geotiles[surface_mask] = geotiles[surface_mask][index_test, :]
        end

        # add a single rgi column
        add_single_rgi_column!(geotiles[surface_mask])
    end

    # align dv_altim and dv_gemb
    begin
        synthesized_gemb_fit = replace(path2runs[1], ".jld2" => "_gembfit.arrow")
        dh = FileIO.load(first(path2runs), "dh_hyps")
        ddate = dims(dh, :date)
        ddate_gemb = dims(dv_gemb, :date)
        index = .!isnan.(dv_gemb[geotile=1, pscale=1, mscale=1])
        ex_gemb = extrema(ddate_gemb[index])
        index = .!isnan.(dh[geotile=1, height=1])
        ex_dv = extrema(ddate[vec(index)])
        ex = max(ex_gemb[1], ex_dv[1]), min(ex_gemb[2], ex_dv[2])
        index_dv = (ddate .>= ex[1]) .& (ddate .<= ex[2])
        index_gemb = (ddate_gemb .>= ex[1]) .& (ddate_gemb .<= ex[2])
        sum(index_dv) == sum(index_gemb) ? nothing : error("index_dv and index_gemb do not have the same length: $(ex)")
        dates2extract = ddate[index_dv].val

        dv_gemb = dv_gemb[date=At(ddate_gemb[index_gemb].val)]
    end

    # I'm getting an "geotile volume change contains all NaNs" when run using Threads.. no clue why
    @showprogress desc = "Calibrating GEMB model to altimetry data" for binned_synthesized_file in path2runs

        synthesized_gemb_fit = replace(binned_synthesized_file, ".jld2" => "_gembfit.arrow")

        if (isfile(synthesized_gemb_fit) && (isnothing(force_remake_before)) || (Dates.unix2datetime(mtime(synthesized_gemb_fit)) > force_remake_before) && isnothing(single_geotile_test))
            printstyled("    -> Skipping $(synthesized_gemb_fit) because it was created after force_remake_before:$force_remake_before\n"; color=:light_green)
            continue
        else

            run_parameters = binned_filled_fileparts(synthesized_gemb_fit)
            surface_mask = run_parameters.surface_mask

            dh = FileIO.load(binned_synthesized_file, "dh_hyps")

            # Convert elevation change to volume change
            dv_altim = dh2dv_geotile(dh[date=At(dates2extract)], area_km2[surface_mask])

            all_nans = dropdims(all(isnan.(dv_altim), dims=:date), dims=:date)

            if any(all_nans)
                error("geotile volume change contains all NaNs: $(binned_synthesized_file)")
            end

            # Find optimal fit to GEMB data
            # There are issues with calibrating the SMB model to individual geotiles since glaciers 
            # can cross multiple geotiles, therefore we calibrate the model for groups of 
            # distinct geotiles.
            if isnothing(single_geotile_test)
            
                gembfit = gemb_bestfit_grouped(dv_altim, dv_gemb, geotiles[surface_mask]; seasonality_weight, distance_from_origin_penalty, mscale_to_pscale_weight, calibrate_to, is_scaling_factor)

                # having issues with pscale and mscale seeming somewhat random... I'm thinking that it's realated to something that is no Thread safe
                dgeotile_test = Dim{:geotile}([geotiles_golden_test[1], "lat[+58+60]lon[-138-136]"])
                pscale_range_test = DimStack(DimArray([1.8, 1.5], dgeotile_test; name="min"), DimArray([3.0, 2.5], dgeotile_test; name="max"))

                for geotile_test in dgeotile_test
                    geotile_sanity_check = geotile_test
                    gt_index = findfirst(gembfit.id .== geotile_sanity_check)
                    geotile_row = gembfit[gt_index, :]
                    pscale = round(geotile_row["pscale"], digits=2)
                    mscale = round(geotile_row["mscale"], digits=2)
                    dv_altim_mean = dv_altim[geotile=At(geotile_row.id)]
                    dv_altim_mean = mean(dv_altim_mean[.!isnan.(dv_altim_mean)])
                    dv_gemb_mean = dv_gemb[geotile=At(geotile_row.id)]
                    dv_gemb_mean = mean(dv_gemb_mean[.!isnan.(dv_gemb_mean)])

                    println("\n$(geotile_sanity_check): pscale = $pscale, mscale = $mscale, dv_altim_mean = $(round(dv_altim_mean, digits=2)), dv_gemb_mean = $(round(dv_gemb_mean, digits=2)), file = $(binned_synthesized_file)")
                    if !(pscale_range_test[:min][geotile=At(geotile_test)] < pscale < pscale_range_test[:max][geotile=At(geotile_test)])
                        error("pscale for $geotile_sanity_check is out of bounds: $pscale")
                    end
                end

                gembfit[!, :area_km2] = sum.(geotiles[surface_mask].area_km2)
                gembfit.geometry = extent2rectangle.(gembfit.extent)
                gembfit = gembfit[:, Not(:extent)]
                GeoDataFrames.write(synthesized_gemb_fit, gembfit)
            else
                cost, geotiles_in_group = gemb_bestfit_grouped(dv_altim, dv_gemb, geotiles[surface_mask]; seasonality_weight, distance_from_origin_penalty, mscale_to_pscale_weight, calibrate_to, is_scaling_factor)
                return (cost, dv_altim, geotiles_in_group)
            end
        end
    end
end


"""
    gemb2geotile(gemb, geotile_row; geotile_buffer=1000, min_gemb_coverage=0.75, ddate::Dim{:date})

Aggregate GEMB model output to a geotile by averaging over grid points within the geotile extent and height bins, for each precipitation scaling factor.

# Arguments
- `gemb`: Dictionary of GEMB output arrays.
- `geotile_row`: Row from geotile table, containing extent and area information.
- `geotile_buffer`: Buffer distance (meters) to expand geotile extent if coverage is insufficient (default: 1000).
- `min_gemb_coverage`: Minimum fraction of geotile area that must be covered by GEMB data (default: 0.75).
- `ddate`: Date dimension (Dim{:date}).

# Returns
- `Dict` with variables averaged within the geotile, binned by height and precipitation scale.

# Notes
- Buffers geotile extent iteratively until minimum coverage is achieved.
- Returns zeros if geotile area is zero.
"""
function gemb2geotile(gemb, geotile_row; geotile_buffer=1000, min_gemb_coverage=0.75, ddate::Dim{:date})

    if min_gemb_coverage > 1
        error("min_gemb_coverage must be less than 1 but is $min_gemb_coverage")
    end

    unique_precipitation_scale = unique(gemb["precipitation_scale"])
    vars = setdiff(collect(keys(gemb)), ("precipitation_scale", "elevation_delta", "latitude", "longitude", "height", "date", "height_effective"))

    height_range, height_center = project_height_bins()
    dheight = Dim{:height}(height_center)
    dpscale = Dim{:pscale}(sort(unique_precipitation_scale))

    gemb0 = Dict()
    for k in vcat(vars, "nobs")
        gemb0[k] = fill(NaN, (ddate, dheight, dpscale))
    end

    geotile_extent = geotile_row.extent
    latitude_distance, longitude_distance = meters2lonlat_distance(geotile_buffer, mean(geotile_extent.Y))
    has_data = falses(size(height_center))
    in_geotile = falses(size(gemb["longitude"]))

    geotile_total_area = sum(geotile_row.area_km2)
    if geotile_total_area == 0
        for k in vars
            gemb0[k][:] .= 0.
        end
        return gemb0
    end

    cnt = 0
    while (sum(geotile_row.area_km2[has_data]) / geotile_total_area) < min_gemb_coverage
        in_geotile = vec([within(geotile_extent, x, y) for (x, y) in zip(gemb["longitude"], gemb["latitude"])])
        h = gemb["height_effective"][in_geotile]

        for i in 1:(length(height_range)-1)
            has_data[i] = any((h .>= height_range[i]) .& (h .< height_range[i+1]))
        end

        geotile_extent = Extents.buffer(geotile_extent, (X=longitude_distance, Y=latitude_distance))
        cnt += 1

        if cnt > 100
            error("Buffering $(geotile_row.id) extent by $(round(Int,geotile_buffer*cnt/1000)) km, to achieve minimum GEMB coverage of $min_gemb_coverage, has taken too long. Check the code (maybe check coverage of gemb data) for potential issues.")
        end
    end

    if cnt > 10
        @warn "Buffering $(geotile_row.id) extent by $(round(Int,geotile_buffer*cnt/1000)) km, to achieve minimum GEMB coverage of $min_gemb_coverage"
    end

    Threads.@threads for i in 1:(length(height_range)-1)
        index_height = (gemb["height_effective"] .>= height_range[i]) .& (gemb["height_effective"] .< height_range[i+1])
        height0 = height_center[i]

        if !any(index_height)
            continue
        end

        for pscale in dpscale
            index_precipitation = gemb["precipitation_scale"] .== pscale
            index = index_precipitation .& in_geotile .& index_height
            sum_index = sum(index)

            if sum_index == 0
                continue
            elseif sum_index == 1
                for k in vars
                    gemb0[k][:, At(height0), At(pscale)] = gemb[k][index, :]
                end
            else
                for k in vars
                    gemb0[k][:, At(height0), At(pscale)] = mean(gemb[k][index, :], dims=1)
                end
            end

            valid_date = .!isnan.(gemb0[vars[1]][:, At(height0), At(pscale)])
            gemb0["nobs"][valid_date, At(height0), At(pscale)] .= sum(index)
        end
    end

    return gemb0
end

"""
    read_gemb_files(
        gemb_files,
        gembinfo;
        vars2extract = ["acc", "refreeze", "melt", "rain", "ec", "fac", "smb", "dv"],
        date_range,
        date_center,
        path2land_sea_mask = nothing,
        minimum_land_coverage_fraction = 0.70,
    )

Read and combine data from multiple GEMB output files (.mat), extracting requested variables and metadata, and standardizing into a single dictionary suitable for geotile binning.

All mass variables are returned in units of meters ice equivalent, except "fac" which is in meters of air.

# Arguments
- `gemb_files::Vector{String}`: Paths to GEMB output files to read and merge.
- `gembinfo`: Structure or dictionary containing precipitation scale (`precipitation_scale`) and elevation delta (`elevation_delta`) information to assign to outputs.
- `vars2extract::Vector{String}`: List of variables to extract from each file (default includes main GEMB outputs).
- `date_range`: Either two-element tuple or array, used to convert times to decimal years and select matching slices.
- `date_center`: Central value for time dimension for all files (must be provided). !!! Needed to ensure exact date matching with other data sources !!!
- `path2land_sea_mask::Union{String,Nothing}`: Path to land-sea mask raster file to screen data (optional).
- `minimum_land_coverage_fraction::Float64`: Minimum fraction of land required for a point to be included (default 0.70).

# Returns
- `Dict{String, Any}`: Dictionary containing merged variables. Includes:
    - All selected variables as arrays concatenated across all files
    - Coordinates: "latitude", "longitude", "height", "date"
    - "precipitation_scale", "elevation_delta": Array of assigned class values per data point
    - "height_effective": height plus elevation delta (for use in further gridding).
    - If land mask is used: "land_fraction"
"""
function read_gemb_files(gemb_files, gembinfo; vars2extract=["acc", "refreeze", "melt", "rain", "ec", "fac", "smb", "dv"], date_range, date_center, path2land_sea_mask=nothing, minimum_land_coverage_fraction=0.70)
    datebin_edges = decimalyear.(date_range)

    gemb = Vector{Any}(undef, length(gemb_files))
    pscale = gembinfo.precipitation_scale
    Δelevation = gembinfo.elevation_delta
    @showprogress desc = "Reading raw GEMB output into memory, this will take ~2 min [peak memory usage: 150GB]" Threads.@threads for i in eachindex(gemb_files)
        gemb_file = gemb_files[i]
        gemb0 = gemb_read2(gemb_file; vars=vars2extract, datebin_edges)

        # Determine precipitation scale and elevation delta indices
        if (length(gembinfo.elevation_delta) == 1) && (length(gembinfo.precipitation_scale) == 1)
            precip_scale_ind = "p1"
            elev_delta_ind = "t1"
        else
            ind1 = findlast("_p", gemb_file)
            ind2 = findlast("_t", gemb_file)
            ind3 = findlast(".mat", gemb_file)
            precip_scale_ind = gemb_file[ind1[end]:ind2[1]-1]
            elev_delta_ind = gemb_file[ind2[end]:ind3[1]-1]
        end

        gemb0["precipitation_scale"] = pscale[At(precip_scale_ind)]
        gemb0["elevation_delta"] = Δelevation[At(elev_delta_ind)]

        gemb[i] = gemb0
    end

    # If only one class, merge all into a single dictionary
    if (length(gembinfo.elevation_delta) == 1) && (length(gembinfo.precipitation_scale) == 1)
        gemb0 = Dict()
        for k in keys(gemb[1])
            if k !== "date"
                gemb0[k] = reduce(vcat, [g[k] for g in gemb])
            else
                gemb0[k] = gemb[1][k]
            end
        end
        gemb = [gemb0]
    end

    # Merge all simulations into a dictionary of arrays
    latitude = vec(reduce(vcat, getindex.(gemb, Ref("latitude"))))
    longitude = vec(reduce(vcat, getindex.(gemb, Ref("longitude"))))
    height0 = vec(reduce(vcat, getindex.(gemb, Ref("height"))))

    # Ensure longitude is in -180 to 180 range
    index = (longitude .< -180)
    longitude[index] .+= 360
    index = (longitude .> 180)
    longitude[index] .-= 360

    precipitation_scale0 = Float64[]
    elevation_delta0 = Float64[]
    for g in gemb
        precipitation_scale0 = vcat(precipitation_scale0, fill(g["precipitation_scale"], length(g["latitude"])))
        elevation_delta0 = vcat(elevation_delta0, fill(g["elevation_delta"], length(g["latitude"])))
    end

    gemb1 = Dict(
        "date" => date_center,
        "latitude" => latitude,
        "longitude" => longitude,
        "height" => height0,
        "precipitation_scale" => precipitation_scale0,
        "elevation_delta" => elevation_delta0
    )

    foo = Dict()
    for k in setdiff(vars2extract, ["date", "latitude", "longitude", "height"])
        foo[k] = reduce(vcat, getindex.(gemb, Ref(k)))
    end
    gemb = merge(foo, gemb1)

    # Add effective height (height + elevation_delta)
    gemb["height_effective"] = gemb["height"] .+ gemb["elevation_delta"]

    # add land fraction to gemb
    if !isnothing(path2land_sea_mask)
        land_sea_mask = Raster(path2land_sea_mask)[:, :, 1]

        xpt = copy(gemb["longitude"])
        index = xpt .< 0
        xpt[index] .+= 360
        ypt = gemb["latitude"]

        pt = tuple.(X.(Near.(xpt)), Y.(Near.(ypt)))
        gemb["land_fraction"] = land_sea_mask[pt].data

        index_valid = gemb["land_fraction"] .>= minimum_land_coverage_fraction

        n = size(gemb["land_fraction"], 1)
        for k in keys(gemb)
            if size(gemb[k], 1) == n
                if ndims(gemb[k]) == 1
                    gemb[k] = gemb[k][index_valid]
                else
                    gemb[k] = gemb[k][index_valid, :]
                end
            end
        end

    end
    return gemb
end


"""
    dv_gemb(gemb, area_km2, dpscale_new)

Compute volume change (ΔV) for GEMB output, interpolated to new precipitation scaling factors.

# Arguments
- `gemb`: Dictionary of GEMB output arrays, each keyed by variable name.
- `area_km2`: Array of grid cell areas in km².
- `dpscale_new`: Array of new precipitation scaling factors to interpolate to.

# Returns
- `Dict` mapping each variable to a 3D array of volume change (date × pscale × mscale
"""
function dv_gemb(gemb, area_km2, dpscale_new)
    k = first(keys(gemb))
    dpscale = dims(gemb[k], :pscale)
    dmscale = dims(gemb[k], :mscale)
    ddate = dims(gemb[k], :date)

    # Ensure all old pscales are present in new pscales
    if !isempty(setdiff(dpscale, dpscale_new))
        error("old pscales [$(dpscale)] must exist in new pscales [$(dpscale_new)]")
    end

    # Initialize output dictionary
    dv_gemb = Dict()
    for k in keys(gemb)
        dv_gemb[k] = fill(NaN, (ddate, dpscale_new, dmscale))
    end

    Threads.@threads for k in collect(keys(gemb))
        valid = .!isnan.(gemb[k])
        (date_range, height_range, pscale_range, mscale_range) = validrange(valid)
        for date in ddate[date_range]
            for mscale in dmscale[mscale_range]
                dv = @d dropdims(sum(gemb[k][date=At(date), mscale=At(mscale)][height_range, :] .* area_km2[height_range] ./ 1000, dims=:height), dims=:height)
                dv_interp = DataInterpolations.LinearInterpolation(dv, val(dpscale))
                dv_gemb[k][date=At(date), mscale=At(mscale)] = dv_interp(val(dpscale_new))
            end
        end
    end
    return dv_gemb
end


"""
    gemb_ensemble_dv(; gemb_run_id, vars2load=["runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze"])

Load GEMB geotile ensemble data for a given run ID. Returns a dictionary of variables specified in `vars2load` from the corresponding GEMB geotile volume change file. If `vars2load` is `nothing`, loads all variables from the file.

# Arguments
- `gemb_run_id`: Identifier for the GEMB run.
- `dpscale_expanded_increment`: Step size for increasing the resolution or range of `pscale`.
- `dmscale_expanded_increment`: Step size for increasing the resolution or range of `mscale`.
- `vars2load`: Array of variable names to load (default: common mass balance variables). If `nothing`, loads all variables.

# Returns
- `Dict` containing the requested variables (DimStack or similar per variable).

# Examples
```julia
julia> dv_gemb = gemb_ensemble_dv(; gemb_run_id=4)
julia> smb = dv_gemb[:smb]
```
"""
function gemb_ensemble_dv(; gemb_run_id=4)
    gembinfo = gemb_info(; gemb_run_id)
    gemb_geotile_filename_dv = replace(gembinfo.filename_gemb_combined, ".jld2" => "_geotile_dv.jld2")
    dv_gemb = FileIO.load(gemb_geotile_filename_dv, "gemb_dv")
    
    gemb_add_derived_vars!(dv_gemb)

    return dv_gemb

end

"""
    gemb_add_derived_vars!(dv_gemb)

Add derived variables (smb, runoff, dv) to a GEMB volume-change dictionary in-place.

Computes smb = acc - melt + refreeze - ec, runoff = melt - refreeze, dv = smb + fac.

# Arguments
- `dv_gemb`: Dictionary of GEMB variables (must contain acc, melt, refreeze, ec, fac)

# Returns
- The modified dictionary with :smb, :runoff, :dv added.

# Examples
```julia
julia> gemb_add_derived_vars!(dv_gemb)
julia> dv = dv_gemb[:dv]
```
"""
function gemb_add_derived_vars!(dv_gemb)
    dv_gemb = merge(dv_gemb, (; smb = dv_gemb[:acc] .- dv_gemb[:melt] .+ dv_gemb[:refreeze] .- dv_gemb[:ec], runoff = dv_gemb[:melt] .- dv_gemb[:refreeze]))
    dv_gemb = merge(dv_gemb, (; dv = dv_gemb[:smb] .+ dv_gemb[:fac]))

    return dv_gemb
end

"""
    gemb_dv_sample(pscale, mscale, dv_gemb)

Sample GEMB volume change at a single (pscale, mscale) by linear interpolation.

# Arguments
- `pscale`: Precipitation scaling factor
- `mscale`: Elevation offset (m)
- `dv_gemb`: DimArray with dimensions (date, pscale, mscale)

# Returns
- DimArray of volume change with dimension :date only.

# Examples
```julia
julia> dv_sampled = gemb_dv_sample(1.0, 0.0, dv_gemb)
```
"""
function gemb_dv_sample(pscale, mscale, dv_gemb)

    ddate = dims(dv_gemb, :date)
    dpscale = dims(dv_gemb, :pscale)
    dmscale = dims(dv_gemb, :mscale)

    itp = Interpolations.interpolate((eachindex(ddate), dpscale.val, dmscale.val), dv_gemb.data, (NoInterp(), Gridded(Linear()), Gridded(Linear())))
    gemb_dv_out = DimArray(itp.(eachindex(ddate), Ref(pscale), Ref(mscale)),ddate)

    return gemb_dv_out  
end

"""
    gemb_dv_sample!(gemb_dv_out::DimArray, pscale, mscale, dv_gemb::DimArray)

Sample GEMB volume change at (pscale, mscale) into a pre-allocated DimArray (in-place).

# Arguments
- `gemb_dv_out`: Pre-allocated DimArray with :date dimension (modified in-place)
- `pscale`: Precipitation scaling factor
- `mscale`: Elevation offset (m)
- `dv_gemb`: DimArray with dimensions (date, pscale, mscale)

# Returns
- The modified gemb_dv_out.

# Examples
```julia
julia> gemb_dv_sample!(gemb_dv_out, 1.0, 0.0, dv_gemb)
```
"""
function gemb_dv_sample!(gemb_dv_out::DimArray, pscale, mscale, dv_gemb::DimArray)

    ddate = dims(dv_gemb, :date)
    dpscale = dims(dv_gemb, :pscale)
    dmscale = dims(dv_gemb, :mscale)

    itp = Interpolations.interpolate((eachindex(ddate), dpscale.val, dmscale.val), dv_gemb.data, (NoInterp(), Gridded(Linear()), Gridded(Linear())))
    gemb_dv_out[:] .= itp.(eachindex(ddate), Ref(pscale), Ref(mscale))

    return gemb_dv_out
end

"""
    gemb_dv_sample!(gemb_dv_out, pscale_index, mscale_index, dv_gemb)

Sample GEMB volume change at (pscale_index, mscale_index) into a pre-allocated array (in-place).

# Arguments
- `gemb_dv_out`: Pre-allocated output array (modified in-place)
- `pscale_index`: Index into pscale dimension
- `mscale_index`: Index into mscale dimension
- `dv_gemb`: 3D array (date × pscale × mscale)

# Returns
- The modified gemb_dv_out.

# Examples
```julia
julia> gemb_dv_sample!(gemb_dv_out, 1, 1, dv_gemb)
```
"""
function gemb_dv_sample!(gemb_dv_out, pscale_index, mscale_index, dv_gemb)

    itp = Interpolations.interpolate((axes(dv_gemb, 1), axes(dv_gemb, 2), axes(dv_gemb, 3)), dv_gemb, (NoInterp(), Gridded(Linear()), Gridded(Linear())))
    gemb_dv_out[:] .= itp.(axes(dv_gemb, 1), pscale_index, mscale_index)

    return gemb_dv_out
end





"""
    gemb_dv_interpolate(dv_gemb, dpscale_expanded_increment, dmscale_expanded_increment)

Interpolate a GEMB DimStack of volume change (`dv_gemb`) across extended precipitation scaling (`pscale`) and temperature perturbation (`mscale`) classes.

# Arguments
- `dv_gemb`: DimStack or dictionary containing GEMB volume change data, with dimensions including :geotile, :date, :pscale, and :mscale.
- `dpscale_expanded_increment`: Step size for increasing the resolution or range of `pscale`.
- `dmscale_expanded_increment`: Step size for increasing the resolution or range of `mscale`.

# Returns
- `gemb_dv_new`: Interpolated `DimStack` structure with expanded/adjusted :pscale and :mscale dimensions.

# Description
- This function first determines the expanded/desired ranges for `pscale` and `mscale` using the provided increments.
- It then initializes a new DimStack with NaN-filled arrays matching the desired output dimensions.
- For each variable in `dv_gemb`, a linear gridded interpolation is constructed across its original (`pscale`, `mscale`) grid for each geotile and date, and evaluated at the new grid, storing the results in `gemb_dv_new`.

# Examples
```julia
julia> gemb_dv_new = gemb_dv_interpolate(dv_gemb; dpscale_expanded_increment=0.25, dmscale_expanded_increment=0.5)
```
"""
function gemb_dv_interpolate(dv_gemb; dpscale_expanded_increment=0.25, dmscale_expanded_increment=0.5)

    dgeotile = dims(dv_gemb, :geotile)
    ddate = dims(dv_gemb, :date)
    dpscale = dims(dv_gemb, :pscale)
    dmscale = dims(dv_gemb, :mscale)

    pscale_start = minimum(dpscale)
    pscale_start < 1 ? pscale_start = -1 / pscale_start : pscale_start
    pscale_start = ceil(pscale_start, digits=0)
    pscale_end = maximum(dpscale)
    pscale_end < 1 ? pscale_end = -1 / pscale_end : pscale_end
    pscale_end = floor(pscale_end, digits=0)
    pscale_new = pscale_start:dpscale_expanded_increment:pscale_end
    pscale_new = pscale_new[(1 .<= pscale_new).|(pscale_new.<-1)]
    pscale_new[pscale_new.<0] = -1 ./ pscale_new[pscale_new.<1]
    dpscale_new = Dim{:pscale}(pscale_new)

    mscale_start = minimum(dmscale)
    mscale_start = ceil(mscale_start, digits=0)
    mscale_end = maximum(dmscale)
    mscale_end = floor(mscale_end, digits=0)
    mscale_new = mscale_start:dmscale_expanded_increment:mscale_end
    dmscale_new = Dim{:mscale}(mscale_new)

    gemb_dv_new = DimStack([DimArray(fill(NaN, (dgeotile, ddate, dpscale_new, dmscale_new); name=Symbol(k))) for k in keys(dv_gemb)]...)

    index_date_valid = dropdims(any(.!isnan.(dv_gemb[first(keys(dv_gemb))][pscale=At(1), mscale=At(0)]), dims=:geotile), dims=:geotile)
    index_date_range, = validrange(index_date_valid)
    daterange = ddate[first(index_date_range)] .. ddate[last(index_date_range)]

    if false
        k = :melt
        begin
            geotile_id = geotiles_golden_test[1]
            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], title="$geotile_id: $(k)")
            cmap_bar = Makie.resample_cmap(:thermal, (length(dpscale)))
            cnt = 1
            for pscale in dpscale
                lines!(dv_gemb[k][date=daterange, geotile=At(geotile_id), pscale=At(pscale), mscale=At(0)]; color=cmap_bar[cnt], label="pscale: $(pscale)")
                cnt += 1
            end
            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)

            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], title="$geotile_id: $(k)")
            cmap_bar = Makie.resample_cmap(:thermal, (length(dmscale)))
            cnt = 1


            for mscale in dmscale
                lines!(dv_gemb[k][date=daterange, geotile=At(geotile_id), pscale=At(1), mscale=At(mscale)]; color=cmap_bar[cnt], label="mscale: $(mscale)")
                cnt += 1
            end

            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)

            dates2plot = [DateTime(2000, 5, 1), DateTime(2010, 12, 30), DateTime(2015, 9, 25)]
            geotile_id = geotiles_golden_test[1]
            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], title="$geotile_id: $(k) as a function of pscale")
            cmap_bar = Makie.resample_cmap(:thermal, (length(dpscale)))
            cnt = 1
            for date in dates2plot
                lines!(dv_gemb[k][date=Near(date), geotile=At(geotile_id), mscale=At(0)]; color=cmap_bar[cnt], label="date: $(date)")
                cnt += 1
            end
            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)

            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], title="$geotile_id: $(k) as a function of mscale")
            cmap_bar = Makie.resample_cmap(:thermal, (length(dpscale)))
            cnt = 1
            for date in dates2plot
                lines!(dv_gemb[k][date=Near(date), geotile=At(geotile_id), pscale=At(1)]; color=cmap_bar[cnt], label="date: $(date)")
                cnt += 1
            end
            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)
        end
    end

    # linearly interpolate between pscale and mscale... maybe this should be done in fitting function?
    Threads.@threads for k in keys(dv_gemb)
        for date in index_date_range
            for gtidx = eachindex(dgeotile)
                itp = Interpolations.interpolate((dpscale.val, dmscale.val), dv_gemb[k][date=date, geotile=gtidx].data, Gridded(Linear()))
                gemb_dv_new[k][geotile=gtidx, date=date] = itp(dpscale_new.val, dmscale_new.val)
            end
        end
    end

    if false
        k = :melt
        begin
            geotile_id = geotiles_golden_test[1]
            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], title="$geotile_id: $(k)")
            cmap_bar = Makie.resample_cmap(:thermal, (length(dpscale_new)))
            cnt = 1
            for pscale in dpscale_new
                lines!(gemb_dv_new[k][date=daterange, geotile=At(geotile_id), pscale=At(pscale), mscale=At(0)]; color=cmap_bar[cnt], label="pscale: $(pscale)")
                cnt += 1
            end
            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)

            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], title="$geotile_id: $(k)")
            cmap_bar = Makie.resample_cmap(:thermal, (length(dmscale_new)))
            cnt = 1


            for mscale in dmscale_new
                lines!(gemb_dv_new[k][date=daterange, geotile=At(geotile_id), pscale=At(1), mscale=At(mscale)]; color=cmap_bar[cnt], label="mscale: $(mscale)")
                cnt += 1
            end
            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)

            dates2plot = [DateTime(2000, 5, 1), DateTime(2010, 12, 30), DateTime(2015, 9, 25)]
            geotile_id = geotiles_golden_test[1]
            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], title="$geotile_id: $(k) as a function of pscale")
            cmap_bar = Makie.resample_cmap(:thermal, (length(dpscale)))
            cnt = 1
            for date in dates2plot
                lines!(gemb_dv_new[k][date=Near(date), geotile=At(geotile_id), mscale=At(0)]; color=cmap_bar[cnt], label="date: $(date)")
                cnt += 1
            end
            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)

            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], title="$geotile_id: $(k) as a function of mscale")
            cmap_bar = Makie.resample_cmap(:thermal, (length(dpscale)))
            cnt = 1
            for date in dates2plot
                lines!(gemb_dv_new[k][date=Near(date), geotile=At(geotile_id), pscale=At(1)]; color=cmap_bar[cnt], label="date: $(date)")
                cnt += 1
            end
            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)
        end
    end

    return gemb_dv_new
end


"""
    process_gemb_geotiles(gemb, geotiles; maximum_search_iterations=20, height_bins_extrapolate=1, maximum_extrap_fraction=0.1, search_buffer_increment=100_000, show_interp_extrap_plots=false, show_interp_extrap_stats=false)

Process GEMB data and organize it into geotiles with elevation and precipitation classes.

This function processes GEMB model output by:
1. Organizing data into geotiles with consistent spatial and temporal dimensions
2. Applying spatial buffering to ensure minimum coverage requirements
3. Filling data gaps through interpolation and extrapolation across elevation bands
4. Calculating volume changes for each geotile

# Arguments
- `gemb::Dict`: Dictionary containing GEMB model output with keys including variables and metadata
- `geotiles::DataFrame`: DataFrame with geotile definitions, hypsometry, and extent information
- `maximum_search_iterations::Int`: Maximum number of search iterations to populate elevation profile (default: 20)
- `height_bins_extrapolate::Int`: Number of bins above and below valid height range to extrapolate (default: 1)
- `maximum_extrap_fraction::Float64`: Maximum fraction of area below/above valid height range for extrapolation (default: 0.1)
- `search_buffer_increment::Int`: Increment of search buffer in meters with each search iteration (default: 100_000)
- `show_interp_extrap_plots::Bool`: Whether to show interpolation/extrapolation plots (default: false)
- `show_interp_extrap_stats::Bool`: Whether to show interpolation/extrapolation statistics (default: false)


# Returns
- `dv_gemb::DimStack`: DimStack containing processed GEMB data with dimensions (:geotile, :date, :pscale, :mscale)
"""
function process_gemb_geotiles(
    gemb,
    geotiles;
    maximum_search_iterations=20,
    height_bins_extrapolate=1,
    maximum_extrap_fraction=0.1,
    search_buffer_increment=100_000,
    show_interp_extrap_plots=false,
    show_interp_extrap_stats=false,
    single_geotile_test=nothing,
    elevation_classes_method=:none,
)

    unique_precipitation_scale = unique(gemb["precipitation_scale"])
    unique_elevation_delta = unique(gemb["elevation_delta"])
    vars = setdiff(collect(keys(gemb)), ("precipitation_scale", "elevation_delta", "latitude", "longitude", "height", "date", "height_effective", "land_fraction"))

    height_range, height_center = project_height_bins()
    ddate = Dim{:date}(gemb["date"])
    dheight = Dim{:height}(height_center)
    dpscale = Dim{:pscale}(sort(unique_precipitation_scale))
    dΔheight = Dim{:Δheight}(sort(unique_elevation_delta))
    dgeotile = Dim{:geotile}(geotiles.id)

    index_date_valid = .!isnan.(gemb[vars[1]][1, :])
    index_date_range, = validrange(index_date_valid)
    daterange = ddate[first(index_date_range)] .. ddate[last(index_date_range)]

    cmap_bar = Makie.resample_cmap(:dense, (maximum_search_iterations * 3 + 1))

    if !isnothing(single_geotile_test)
        geotiles = geotiles[geotiles.id.==single_geotile_test, :]
        show_interp_extrap_plots = true
    end

    if elevation_classes_method == :none
        gemb_dv0 = DimStack([DimArray(fill(NaN, (dgeotile, ddate, dpscale, dΔheight)); name=k) for k in vars]...)
        gemb_dv0 = DimStack([DimArray(fill(NaN, (dgeotile, ddate, dpscale, dΔheight)); name=k) for k in vars]...)
    else
        mscale_range, mscale_center = project_mscale_bins()
        dmscale = Dim{:mscale}(mscale_center)
        dh1 = dheight.val.step.hi
        max_height_offset = ceil(abs((1000 / -6.5) * maximum(abs.(dmscale.val))) / dh1) * abs(dh1)

        height_center1 = (minimum(dheight.val)+minimum(dΔheight.val)):dh1:((maximum(dheight.val)+maximum(dΔheight.val)))
        height_range1 = (minimum(height_center1)-dh1/2):dh1:(maximum(height_center1)+dh1/2)
        dheight1 = Dim{:height}(height_center1)

        if elevation_classes_method == :Δelevation
            gemb_dv0 = DimStack([DimArray(fill(NaN, (dgeotile, ddate, dpscale, dmscale)); name=k) for k in vars]...)
        elseif elevation_classes_method == :mscale
            mscale_range, mscale_center = project_mscale_bins()
            dmscale = Dim{:mscale}(mscale_center)
            gemb_dv0 = DimStack([DimArray(fill(NaN, (dgeotile, ddate, dpscale, dmscale)); name=k) for k in vars]...)
        end
    end

    @showprogress desc = "Populating geotiles with GEMB data, this will take 2 to 10 min on 128 threads, depending on method" Threads.@threads for geotile_row in eachrow(geotiles)
        #@warn "threads turned off"
        #for geotile_row in eachrow(geotiles)

        geotile_hyps_area_km2 = DimArray(geotile_row.area_km2, dheight; name="area_km2")
        geotile_total_area = sum(geotile_hyps_area_km2)

        if geotile_total_area == 0
            for k in Symbol.(vars)
                gemb_dv0[k][geotile=At(geotile_row.id), date=daterange] .= 0.
            end
            continue
        end

        index_has_ice = geotile_hyps_area_km2 .> 0
        index_has_ice_range, = validrange(index_has_ice)

        search_extent = geotile_row.extent
        cnt = 0
        bar = 1
        search_buffer = 0
        coverage_height_fraction = 0.
        points_in_geotile_count = 0

        if show_interp_extrap_plots
            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], yaxisposition=:right, title="$(geotile_row.id)")
            hidedecorations!(ax)  # hides ticks, grid and labels
            hidespines!(ax)  # hide the frame
        end

        # Populate classes through interpolation and extrapolation... elevation classes are not used to pupulate elevation range
        if elevation_classes_method == :none

            gemb0 = DimStack([DimArray(fill(NaN, (ddate, dheight, dpscale, dΔheight)); name=k) for k in vars]...)

            while isnan(coverage_height_fraction) || (round(coverage_height_fraction, digits=3) < 1.0)
                cnt += 1
                if cnt > 1
                    search_buffer += search_buffer_increment
                    latitude_distance, longitude_distance = meters2lonlat_distance(search_buffer, mean(search_extent.Y))
                    search_extent = Extents.buffer(search_extent, (X=longitude_distance, Y=latitude_distance))
                end

                # points in geotiles 
                index_geotile = within.(Ref(search_extent), gemb["longitude"], gemb["latitude"])

                if !any(index_geotile)
                    bar += 3
                    continue
                end

                index_geotile = findall(index_geotile)

                height0 = gemb["height"][index_geotile]
                elevation_delta0 = gemb["elevation_delta"][index_geotile]
                precipitation_scale0 = gemb["precipitation_scale"][index_geotile]

                index_existing_heights = .!isnan.(gemb0[:fac][pscale=At(1), Δheight=At(0)])
                index_existing_heights = any(index_existing_heights, dims=:date)[:]

                # loop is not thread safe
                for i in eachindex(height_center)[.!index_existing_heights]
                    index_height = height_range[i] .<= height0 .< height_range[i+1]

                    if !any(index_height)
                        continue
                    end

                    for Δheight in dΔheight
                        index_Δheight = Δheight .== elevation_delta0

                        for pscale in dpscale
                            index_precipitation = precipitation_scale0 .== pscale
                            index_valid = index_geotile[index_height.&index_Δheight.&index_precipitation]

                            for k in vars
                                gemb0[Symbol(k)][height=At(height_center[i]), pscale=At(pscale), Δheight=At(Δheight)] = mean(gemb[k][index_valid, :], dims=1)
                            end
                        end
                    end
                end

                valid_ensemble = .!isnan.(gemb0[Symbol(vars[1])][pscale=At(1), Δheight=At(0)])
                index_height_valid = any(valid_ensemble, dims=:date)[:]

                if any(index_height_valid)
                    # Ask questions about the ensemble coverage

                    # plot the ensemble coverage
                    if show_interp_extrap_plots
                        bar_index = (index_height_valid .& .!index_existing_heights)
                        bar_labels = repeat(["raw-$(cnt)"], sum(index_has_ice))
                        bar_labels[.!bar_index[index_has_ice]] .= ""
                        barplot!(ax, dheight[index_has_ice].val, bar_index[index_has_ice]; direction=:x, color=(cmap_bar[bar]), bar_labels=bar_labels, label_position=:center)
                        bar += 1
                    end

                    if false
                        for k in [:fac, :acc, :melt]
                            f = _publication_figure(; columns=1, rows=2)
                            ax = CairoMakie.Axis(f[1, 1]; palette=(:color => Makie.resample_cmap(:thermal, length(findall(index_height_valid)))), title="$(geotile_row.id): $(k)")

                            for height_index in findall(index_height_valid)
                                lines!(ax, gemb0[k][dates=daterange, height=At(dheight[height_index]), pscale=At(1), Δheight=At(0)]; label="$(dheight[height_index])")
                            end

                            f[1, 2] = Legend(f, ax, "height [m]", framevisible=false)
                            display(f)
                        end
                    end

                    # Ask a bunch of questions to determine where extrapolation is appropriate
                    index_height_valid_extrapolate = dilate(index_height_valid, height_bins_extrapolate)
                    index_height_range, = validrange(index_height_valid_extrapolate)

                    if 1 < first(index_height_range) > first(index_has_ice_range)
                        area_below_height_range = sum(geotile_hyps_area_km2[1:(first(index_height_range)-1)])
                        coverage_below_height_range = area_below_height_range ./ geotile_total_area
                    else
                        area_below_height_range = 0
                        coverage_below_height_range = 0
                    end

                    if last(index_height_range) < last(index_has_ice_range)
                        area_above_height_range = sum(geotile_hyps_area_km2[(last(index_height_range)+1):end])
                        coverage_above_height_range = area_above_height_range ./ geotile_total_area
                    else
                        coverage_above_height_range = 0
                    end

                    if (0 < coverage_below_height_range < maximum_extrap_fraction)
                        index_height_range = first(index_has_ice_range):last(index_height_range)
                    end

                    if (0 < coverage_above_height_range < maximum_extrap_fraction)
                        index_height_range = first(index_height_range):last(index_has_ice_range)
                    end

                    index_extrapolate = falses(length(dheight))
                    index_extrapolate[index_height_range] .= true
                    index_interp_range, = validrange(index_height_valid)
                    index_extrapolate[index_interp_range] .= false

                    index_fill = falses(length(dheight))
                    index_fill[index_height_range] .= true
                    index_fill = (index_fill .& .!index_height_valid)

                    # limit size of interpolation gap
                    index_height_range_interpolate, = validrange(index_height_valid)
                    index_height_interpolate = falses(length(dheight))
                    index_height_interpolate[index_height_range_interpolate] .= true
                    index_height_interpolate = (index_height_interpolate .& .!index_height_valid)

                    if any(index_height_interpolate)
                        index_valid_interpolate = .!(true_block_size(index_fill) .> (height_bins_extrapolate * 2 + 1))
                        index_valid_interpolate = dilate(index_valid_interpolate, height_bins_extrapolate)
                        index_fill = index_fill .& index_valid_interpolate
                    end

                    if sum(index_height_valid) > 1
                        for k in Symbol.(vars)
                            itp = extrapolate(Interpolations.interpolate((1:length(ddate), dheight.val[index_height_valid], 1:length(dpscale), 1:length(dΔheight)), gemb0[k][height=findall(index_height_valid)], (NoInterp(), Gridded(Linear()), NoInterp(), NoInterp())), Flat())
                            for date in index_date_range
                                for Δheight in eachindex(dΔheight)
                                    for pscale in eachindex(dpscale)
                                        gemb0[k][date=date, height=index_fill, pscale=pscale, Δheight=Δheight] = itp.(Ref(date), dheight.val[index_fill], Ref(pscale), Ref(Δheight))
                                    end
                                end
                            end
                        end
                    else
                        for k in Symbol.(vars)
                            fill_values = gemb0[k][height=findall(index_height_valid)]

                            for height0 in dheight[index_fill]
                                gemb0[Symbol(k)][height=At(height0)] = fill_values
                            end
                        end
                    end

                    if show_interp_extrap_plots
                        bar_index = index_fill .& .!index_extrapolate
                        bar_labels = repeat(["interp-$(cnt)"], sum(index_has_ice))
                        bar_labels[.!bar_index[index_has_ice]] .= ""
                        barplot!(ax, dheight[index_has_ice].val, bar_index[index_has_ice]; direction=:x, color=(cmap_bar[bar]), bar_labels=bar_labels, label_position=:center)
                        bar += 1

                        bar_index = index_extrapolate
                        bar_labels = repeat(["extrap-$(cnt)"], sum(index_has_ice))
                        bar_labels[.!bar_index[index_has_ice]] .= ""
                        barplot!(ax, dheight[index_has_ice].val, bar_index[index_has_ice]; direction=:x, color=(cmap_bar[bar]), bar_labels=bar_labels, label_position=:center)
                        bar += 1
                    end

                    valid_ensemble = .!isnan.(gemb0[:fac][pscale=At(1), Δheight=At(0)])
                    index_height_valid = any(valid_ensemble, dims=:date)[:]
                    coverage_height_fraction = sum(geotile_hyps_area_km2[index_height_valid]) / geotile_total_area
                end

                if cnt > maximum_search_iterations
                    error("MAXIMUM SEARCH ITERATIONS REACHED: Buffering $(geotile_row.id) search extent by $(round(Int,search_buffer/1000)) km to populate elevation profile")
                end
            end

            # calculate volume change in unit of km3 [GEMB assumes an ice density of 910 kg/m3]
            for k in Symbol.(vars)
                dv = @d gemb0[k][date=daterange] .* geotile_hyps_area_km2 ./ 1000
                dv[isnan.(dv)] .= 0
                dv = sum(dv; dims=:height)
                gemb_dv0[k][geotile=At(geotile_row.id), date=daterange] = dropdims(dv, dims=:height)
            end
        else
            # ensure that height range covers desired mscale classes
            gemb1 = DimStack([DimArray(fill(NaN, (ddate, dheight1, dpscale)); name=k) for k in vars]...)

            index_has_ice1 = (dheight1.val .>= (dheight[first(index_has_ice_range)] - max_height_offset)) .& (dheight1.val .<= (dheight[last(index_has_ice_range)] + max_height_offset))
            index_has_ice_range1, = validrange(index_has_ice1)

            # use simulated elevation classes to populate elevation range
            while isnan(coverage_height_fraction) || (round(coverage_height_fraction, digits=3) < 1.0)

                cnt += 1
                if cnt > 1
                    search_buffer += search_buffer_increment
                    latitude_distance, longitude_distance = meters2lonlat_distance(search_buffer, mean(search_extent.Y))
                    search_extent = Extents.buffer(search_extent, (X=longitude_distance, Y=latitude_distance))
                end

                # points in geotiles 
                index_geotile = within.(Ref(search_extent), gemb["longitude"], gemb["latitude"])
                points_in_geotile_count = sum(index_geotile)
                if points_in_geotile_count == 0
                    bar += 3
                    continue
                end

                index_geotile = findall(index_geotile)

                height0 = gemb["height_effective"][index_geotile]
                elevation_delta0 = gemb["elevation_delta"][index_geotile]
                precipitation_scale0 = gemb["precipitation_scale"][index_geotile]

                index_existing_heights = .!isnan.(gemb1[:fac][pscale=At(1)])
                index_existing_heights = any(index_existing_heights, dims=:date)[:]

                # loop is not thread safe
                for i in eachindex(height_center1)[.!index_existing_heights]
                    index_height = height_range1[i] .<= height0 .< height_range1[i+1]

                    if !any(index_height)
                        continue
                    end

                    for pscale in dpscale
                        index_precipitation = (precipitation_scale0 .== pscale)
                        index_valid = index_geotile[index_height.&index_precipitation]

                        for k in vars
                            gemb1[Symbol(k)][height=At(height_center1[i]), pscale=At(pscale)] = mean(gemb[k][index_valid, :], dims=1)
                        end
                    end
                end

                valid_ensemble = .!isnan.(gemb1[Symbol(vars[1])][pscale=At(1)])
                index_height_valid = any(valid_ensemble, dims=:date)[:]

                if sum(index_height_valid) != 0
                    # Ask questions about the ensemble coverage

                    # plot the ensemble coverage
                    if show_interp_extrap_plots
                        bar_index = (index_height_valid .& .!index_existing_heights)
                        bar_labels = repeat(["raw-$(cnt)"], sum(index_has_ice1))
                        bar_labels[.!bar_index[index_has_ice1]] .= ""
                        barplot!(ax, dheight1[index_has_ice1].val, bar_index[index_has_ice1]; direction=:x, color=(cmap_bar[bar]), bar_labels=bar_labels, label_position=:center)
                        bar += 1
                    end

                    if false
                        for k in [:fac, :acc, :melt]
                            f = _publication_figure(; columns=1, rows=2)
                            ax = CairoMakie.Axis(f[1, 1]; palette=(:color => Makie.resample_cmap(:thermal, length(findall(index_height_valid)))), title="$(geotile_row.id): $(k)")

                            for height_index in findall(index_height_valid)
                                lines!(ax, gemb1[k][dates=daterange, height=At(dheight1[height_index]), pscale=At(1)]; label="$(dheight1[height_index])")
                            end
                            f[1, 2] = Legend(f, ax, "height [m]", framevisible=false)
                            display(f)
                        end
                    end


                    # Ask a bunch of questions to determine where extrapolation is appropriate
                    index_height_valid_extrapolate = dilate(index_height_valid, height_bins_extrapolate)
                    index_height_range, = validrange(index_height_valid_extrapolate)

                    if 1 < first(index_height_range) > first(index_has_ice_range1)
                        area_below_height_range = sum(geotile_hyps_area_km2[height=-Inf .. (dheight1[first(index_height_range)-1] + max_height_offset)])

                        coverage_below_height_range = area_below_height_range ./ geotile_total_area
                    else
                        area_below_height_range = 0
                        coverage_below_height_range = 0
                    end

                    if last(index_height_range) < last(index_has_ice_range1)
                        area_above_height_range = sum(geotile_hyps_area_km2[height=(dheight1[last(index_height_range)+1] - max_height_offset) .. Inf])
                        coverage_above_height_range = area_above_height_range ./ geotile_total_area
                    else
                        coverage_above_height_range = 0
                    end

                    if (0 < coverage_below_height_range < maximum_extrap_fraction)
                        index_height_range = first(index_has_ice_range1):last(index_height_range)
                    end

                    if (0 < coverage_above_height_range < maximum_extrap_fraction)
                        index_height_range = first(index_height_range):last(index_has_ice_range1)
                    end

                    index_extrapolate = falses(length(dheight1))
                    index_extrapolate[index_height_range] .= true

                    index_interp_range, = validrange(index_height_valid)
                    index_extrapolate[index_interp_range] .= false

                    index_fill = falses(length(dheight1))
                    index_fill[index_height_range] .= true

                    # limit size of interpolation gap
                    index_height_range_interpolate, = validrange(index_height_valid)
                    index_height_interpolate = falses(length(dheight1))
                    index_height_interpolate[index_height_range_interpolate] .= true

                    if any(index_height_interpolate)
                        index_fill = index_fill .| index_height_interpolate
                    end

                    if sum(index_height_valid) > 2
                        for k in Symbol.(vars)
                            for date in ddate.val[index_date_range]
                                for pscale in dpscale
                                    itp = DataInterpolations.QuadraticInterpolation(gemb1[k][height=findall(index_height_valid), pscale=At(pscale), date=At(date)], dheight1.val[index_height_valid]; extrapolation=ExtrapolationType.Constant)
                                    gemb1[k][date=At(date), height=At(dheight1.val[index_fill]), pscale=At(pscale)] = itp.(dheight1.val[index_fill])
                                end
                            end
                        end
                    elseif sum(index_height_valid) > 1
                        for k in Symbol.(vars)
                            for date in ddate.val[index_date_range]
                                for pscale in dpscale
                                    itp = DataInterpolations.LinearInterpolation(gemb1[k][height=findall(index_height_valid), pscale=At(pscale), date=At(date)], dheight1.val[index_height_valid]; extrapolation=ExtrapolationType.Constant)
                                    gemb1[k][date=At(date), height=At(dheight1.val[index_fill]), pscale=At(pscale)] = itp.(dheight1.val[index_fill])
                                end
                            end
                        end
                    else
                        for k in Symbol.(vars)
                            fill_values = gemb1[k][height=findall(index_height_valid)]

                            for height0 in dheight1[index_fill]
                                gemb1[Symbol(k)][height=At(height0)] = fill_values
                            end
                        end
                    end

                    if show_interp_extrap_plots
                        bar_index = index_fill .& .!index_extrapolate .& .!index_height_valid
                        bar_labels = repeat(["interp-$(cnt)"], sum(index_has_ice1))
                        bar_labels[.!bar_index[index_has_ice1]] .= ""
                        barplot!(ax, dheight1[index_has_ice1].val, bar_index[index_has_ice1]; direction=:x, color=(cmap_bar[bar]), bar_labels=bar_labels, label_position=:center)
                        bar += 1

                        bar_index = index_extrapolate
                        bar_labels = repeat(["extrap-$(cnt)"], sum(index_has_ice1))
                        bar_labels[.!bar_index[index_has_ice1]] .= ""
                        barplot!(ax, dheight1[index_has_ice1].val, bar_index[index_has_ice1]; direction=:x, color=(cmap_bar[bar]), bar_labels=bar_labels, label_position=:center)
                        bar += 1
                    end

                    valid_ensemble = .!isnan.(gemb1[:fac][pscale=At(1)])
                    index_height_valid = any(valid_ensemble, dims=:date)[:]

                    coverage_height_fraction = 0
                    for mscale in dmscale
                        Δheight = round((mscale * (1000 / -6.5)) / dh1, digits=0) * dh1

                        height_effective = dheight1.val .- Δheight
                        index_height_effective = (height_effective .>= dheight.val[1]) .& (height_effective .<= dheight.val[end])

                        coverage_height_fraction += sum(geotile_hyps_area_km2[height=At(height_effective[index_height_effective])]) / geotile_total_area
                    end
                    coverage_height_fraction = coverage_height_fraction / length(dmscale)
                end

                if cnt > maximum_search_iterations
                    error("MAXIMUM SEARCH ITERATIONS REACHED: Buffering $(geotile_row.id) search extent by $(round(Int,search_buffer/1000)) km to populate elevation profile")
                end
            end

            # calculate volume change in unit of km3 [GEMB assumes an ice density of 910 kg/m3]
            index_height = (dheight1.val .>= dheight.val[1]) .& (dheight1.val .<= dheight.val[end])
            dmscale = dims(gemb_dv0, :mscale)

            for k in [:refreeze, :rain, :acc, :melt]
                # these variables can not be negative change
                v = diff(gemb1[k][date=index_date_range], dims=:date)
                v1 = gemb1[k][date=index_date_range[1]:index_date_range[1]].data
                v = cat(v1, v.data; dims=1)
                v[v.<0] .= 0
                v = cumsum(v; dims=1)

                gemb1[k][date=index_date_range] = v
            end

            for k in Symbol.(vars)

                for pscale in dpscale
                    v = gemb1[k][date=daterange, pscale=At(pscale)]
                    v0 = zeros(ddate, dheight)
                    v0 = v0[date=daterange]

                    for mscale in dmscale
                        if elevation_classes_method == :Δelevation
                            Δheight = round((mscale * (1000 / -6.5)) / dh1, digits=0) * dh1
                            height_effective = dheight1.val .- Δheight
                            index_height_effective = (height_effective .>= dheight.val[1]) .& (height_effective .<= dheight.val[end])
                            v0[:, :] = v[:, index_height_effective]

                        elseif elevation_classes_method == :mscale

                            if k == :melt
                                v0[:, :] = v[:, index_height] .* mscale
                            else
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
        end

        show_interp_extrap_stats && printstyled("$(geotile_row.id): search buffer = $(round(Int,search_buffer/1000)) km, gemb nodes in search extent = $(round(Int,points_in_geotile_count/length(dΔheight)/length(dpscale)))\n"; color=:white)

        if show_interp_extrap_plots
            ax2 = CairoMakie.Axis(f[1, 1], xlabel="area [km²]", ylabel="height [m]")
            lines!(ax2, geotile_hyps_area_km2[index_has_ice].data, dheight[index_has_ice].val)
            display(f)
        end
    end

    # change from Δheight to mscale
    if elevation_classes_method == :none
        # must be sorted in descending order for interpolation to work
        dmscale = Dim{:mscale}(reverse(dΔheight.val .* (-6.5 / 1000)))
        gemb_dv0 = DimStack([DimArray(reverse(gemb_dv0[k].data, dims=4), (dgeotile, ddate, dpscale, dmscale); name=k) for k in Symbol.(vars)]...)
    end

    # these variables can not be negative change
    begin
        index_date_range = findfirst(.!isnan.(gemb_dv0[:melt][1,:, 1, 1])):findlast(.!isnan.(gemb_dv0[:melt][1,:, 1, 1]))
        vx = diff(gemb_dv0[:refreeze][date=index_date_range], dims=:date)
        vx1 = gemb_dv0[:refreeze][date=index_date_range[1]:index_date_range[1]].data
        vx = cat(vx1, vx.data; dims=2)
        
        v_melt = diff(gemb_dv0[:melt][date=index_date_range], dims=:date)
        v1_melt = gemb_dv0[:melt][date=index_date_range[1]:index_date_range[1]].data
        v_melt = cat(v1_melt, v_melt.data; dims=2)

        vx = min.(vx, v_melt)
        vx = cumsum(vx; dims=2)
        gemb_dv0[:refreeze][date=index_date_range] = vx
    end

    gemb_dv0 = gemb_add_derived_vars!(gemb_dv0)

    return gemb_dv0
end

"""
    gemb_altim_cost_group(dpscale_search, dmscale_search, dv_altim, dv_gemb, kwargs2)

Find best (pscale, mscale) for one geotile group by grid search over cost function.

# Arguments
- `dpscale_search`: Iterable of pscale values to try
- `dmscale_search`: Iterable of mscale values to try
- `dv_altim`: Altimetry-derived volume change for the group
- `dv_gemb`: GEMB volume change array for the group
- `kwargs2`: Named tuple of kwargs for model_fit_cost_function

# Returns
- Tuple (pscale_best, mscale_best).

# Examples
```julia
julia> pscale_best, mscale_best = gemb_altim_cost_group(dpscale_search, dmscale_search, dv_altim, dv_gemb, kwargs2)
```
"""
function gemb_altim_cost_group(dpscale_search, dmscale_search, dv_altim, dv_gemb, kwargs2)
    
    f2 = Base.Fix{1}(Base.Fix{2}(Base.Fix{4}(Base.Fix{5}(Base.Fix{6}(gemb_altim_cost!, kwargs2), dv_gemb), dv_altim), deepcopy(dv_altim)), deepcopy(dv_altim))

    # getting very unstable results with ECA, so using grid search instead
    # messing with ECA parameters can lead to getting stuck in a local minimum
    #result = optimize(f2, bounds, ECA())
    #(pscale_best, mscale_best) = minimizer(result)

    #pscale0[i] = pscale_best
    #mscale0[i] = mscale_best

    cost = Inf;
    mscale_best = 1.0;
    pscale_best = 1.0;

    # this can not be threaded
    for pscale in dpscale_search
        for  mscale in dmscale_search
            cost0 = f2([pscale, mscale])
            if cost0 < cost
                cost = cost0
                mscale_best = mscale
                pscale_best = pscale
            end
        end
    end

    return pscale_best, mscale_best
end

"""
    gemb_altim_cost_group!(cost, dv_altim, dv_gemb, kwargs2)

Fill a cost array over (pscale, mscale) for one geotile group (in-place).

# Arguments
- `cost`: DimArray with dimensions (pscale, mscale), filled in-place
- `dv_altim`: Altimetry-derived volume change for the group
- `dv_gemb`: GEMB volume change array for the group
- `kwargs2`: Named tuple of kwargs for model_fit_cost_function

# Returns
- The modified cost array.

# Examples
```julia
julia> gemb_altim_cost_group!(cost, dv_altim, dv_gemb, kwargs2)
```
"""
function gemb_altim_cost_group!(cost, dv_altim, dv_gemb, kwargs2)
    
    f2 = Base.Fix{1}(Base.Fix{2}(Base.Fix{4}(Base.Fix{5}(Base.Fix{6}(gemb_altim_cost!, kwargs2), dv_gemb), dv_altim), deepcopy(dv_altim)), deepcopy(dv_altim))

    # getting very unstable results with ECA, so using grid search instead
    # messing with ECA parameters can lead to getting stuck in a local minimum
    #result = optimize(f2, bounds, ECA())
    #(pscale_best, mscale_best) = minimizer(result)

    #pscale0[i] = pscale_best
    #mscale0[i] = mscale_best


    for pscale in dims(cost, :pscale)
        for  mscale in dims(cost, :mscale)
            cost[pscale=At(pscale), mscale=At(mscale)] = f2([pscale, mscale])
        end
    end

    return cost
end


"""
    dv_adjust4discharge!(gemb, discharge)

Adjust GEMB volume change (dv) for discharge in-place by subtracting cumulative discharge over time.

# Arguments
- `gemb`: Dictionary or DimStack containing :dv (and :date); :dv is modified in-place
- `discharge`: Object with :discharge (rate) and :date; discharge rate is applied over decimal years

# Returns
- The modified gemb (gemb[:dv] updated).

# Examples
```julia
julia> dv_adjust4discharge!(gemb, discharge)
julia> # gemb[:dv] now has discharge removed
```
"""
function dv_adjust4discharge!(gemb, discharge)

    decyear = decimalyear.(dims(gemb, :date)) .- decimalyear.(dims(gemb, :date)[1]);
    gemb[:dv][:] = @d gemb[:dv] .- (decyear .* discharge[:discharge]);
    return gemb
end
