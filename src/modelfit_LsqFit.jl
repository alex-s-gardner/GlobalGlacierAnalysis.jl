# see wave model examples here: https://www.geogebra.org/m/bmyurtag

"""
    ts_fit(
        x, y, t, z, model, coeff0; 
        t_min_range, z_max_outlier, t_bin_edges, t_bin_filter_method, t_bin_filter_threshold, 
    )

Fits a timeseries `model` to `z`. The model is intitialized with `coeff0`` parameters. Returns 
`mean(x)`, `mean(y)`, `count`, optimized `coeff`, and model `rmse`, all calculated after 
removal of outliers that are deteremined using optional inputs of `t_min_range`, 
`z_max_outlier`, `t_bin_edges`, `t_bin_filter_method`, and `t_bin_filter_threshold`.

"""
function ts_fit(
    x::AbstractArray{Float64},                  # x coodinate
    y::AbstractArray{Float64},                  # y coodinate
    t::AbstractArray{Float64},                  # time in decimal years
    z::AbstractArray{Float64},                  # elevation anomaly
    lsq_setup::Altim.LsqSetup;                  # lsq model setup
    count_min::Int = 10,                        # minimum number of valid values for model fit
    t_min_range::Float64 = 0.,                  # only file model it range(t) > t_min_range
    z_max_outlier::Float64 = Inf,               # maximum acceptable absolute deviation from median
    t_bin_edges::StepRangeLen{} = 2002:0.1:2024,# bin edges used for local filtering of outliers
    t_bin_filter_method::F = madnorm,          # method for identifying local outliers
    t_bin_filter_threshold::Number = 10         # threshold for classifying local outliers
) where {F<:Function}

    # initialize default outputs
    x_out = NaN;
    y_out = NaN;
    count = Int32(0);
    coeff = fill(NaN, size(lsq_setup.p))
    rmse = NaN;

    # remove NaNs
    valid = .!isnan.(z)
    x = x[valid];
    y = y[valid];
    t = t[valid];
    z = z[valid];

    if (length(x) >= count_min) && (Altim.range(t) >= t_min_range)
        # global centering 
        z = z .- median(z)

        # remove gross outliesrs
        outlier = (abs.(z) .> z_max_outlier) 

        # using a binned filter, center data and remove outliers
        outlier[.!outlier] = binnedfiltering(t[.!outlier], z[.!outlier], 
            t_bin_edges; method=t_bin_filter_method, threshold=t_bin_filter_threshold)

        # remove any outliers that are > max_outlier m from the centered data
        count = Int32(sum(.!outlier))

        if count < count_min 
            return x_out, y_out, count, coeff, rmse
        end

        x_out = mean(x[.!outlier])
        y_out = mean(y[.!outlier])

        fit = curve_fit(lsq_setup.model, t[.!outlier], z[.!outlier], lsq_setup.p, lower=lsq_setup.p_min, upper=lsq_setup.p_max, autodiff=lsq_setup.autodiff)
        coeff = fit.param
        rmse = sqrt(mean(fit.resid.^2))
    end

    return x_out, y_out, count, coeff, rmse
end


"""
    geotile_ts_fit(geotile, dem, paths, model, coeff0, grid, ts_model_arg; altim_dem_outlier, force_remake)
Fit time series model on a grid
"""
function geotile_ts_fit(
    geotile::Union{DataFrame, DataFrameRow},
    dem, 
    paths::NamedTuple, 
    lsq_setup::Altim.LsqSetup,                               # lsq model setup
    grid::NamedTuple,
    ts_model_arg::NamedTuple;
    altim_dem_outlier = Inf,
    t_minmax = nothing,
    mask = nothing,
    missions = [:gedi, :icesat, :icesat2],
    force_remake = false
) 

    t1 = time()
    outfile = joinpath(paths.height_change, "$(geotile.id).$(dem)")
    
    if isfile(outfile) && !force_remake
        printstyled("    -> model fit to $(geotile.id) $(dem) already exists, skipping\n"; color = :green)
        return
    end

    df = geotile_merge_height(geotile, dem, paths; missions = missions, altim_dem_outlier=altim_dem_outlier, t_minmax=t_minmax , mask=mask)
    if isnothing(df)
        return 
    else
        t2 = time()

        # create a regularly spaced grid using a local low distortion projection
        df, epsg = geotile_utm!(df)

        # this takes 0.3 seconds
        minmax_x = extrema(df.X)
        minmax_y = extrema(df.Y)

        center_extent = (
            x_min=floor(minmax_x[1], digits=-4),
            y_min=floor(minmax_y[1], digits=-4),
            x_max=ceil(minmax_x[2], digits=-4),
            y_max=ceil(minmax_y[2], digits=-4)
        )

        griddef = regular_grid(center_extent, grid.node_spacing, node_width=grid.node_width);

        # sort along y dimension
        p = sortperm(df.Y)
        y = df.Y[p]
        x = df.X[p]
        dh = df.dh[p]
        decyear = df.decyear[p]

        coefficeints = fill(NaN, length(griddef.x_node_center), length(griddef.y_node_center), length(lsq_setup.p))
        x_cord = fill(NaN, length(griddef.x_node_center), length(griddef.y_node_center))
        y_cord = copy(x_cord)
        rmse = copy(x_cord)
        count = zeros(Int32, length(griddef.x_node_center), length(griddef.y_node_center))

        # clear room on memory
        df = [];

        t3 = time()
        Threads.@threads for i in eachindex(griddef.x_node_center)
        #for i in eachindex(griddef.x_node_center)
            x_min = griddef.x_node_center[i] - griddef.node_half_width
            x_max = griddef.x_node_center[i] + griddef.node_half_width
            x_ind = (x .>= x_min) .& (x .<= x_max)

            Int(sum(x_ind)) >= ts_model_arg.count_min || continue

            x0 = x[x_ind]
            y0 = y[x_ind]

            dh0 = dh[x_ind]
            decyear0 = decyear[x_ind]

            # should we continue ?
            Altim.range(decyear0) >= ts_model_arg.t_min_range || continue

            for j in eachindex(griddef.y_node_center)
                y_min = griddef.y_node_center[j] - griddef.node_half_width
                y_max = griddef.y_node_center[j] + griddef.node_half_width

                start = searchsortedfirst(y0, y_min)
                stop = searchsortedlast(y0, y_max)

                (stop-start) >= ts_model_arg.count_min || continue
                Altim.range(@view(decyear0[start:stop])) >= ts_model_arg.t_min_range || continue

                x_cord[i,j], y_cord[i,j], count[i,j], coefficeints[i,j,:], rmse[i,j] = 
                    ts_fit(@view(x0[start:stop]),  @view(y0[start:stop]), @view(decyear0[start:stop]), 
                        dh0[start:stop], lsq_setup; ts_model_arg...)
            end
        end
        t4 = time()
        
        # place into new DataFrame
        begin
            df0 = DataFrame()
            valid = .!isnan.(coefficeints[:, :, 1])
            df0[!, :longitude], df0[!, :latitude] = epsg2epsg(x_cord[valid], y_cord[valid], "EPSG:$epsg", "EPSG:4326", parse_output=true)
            df0[!, :offset] = coefficeints[:, :, 1][valid]
            df0[!, :trend] = coefficeints[:, :, 2][valid]
            df0[!, :amplitude] = (coefficeints[:, :, 3][valid] .^ 2 .+ coefficeints[:, :, 4][valid] .^ 2) .^ 0.5
            df0[!, :phase] = atan.(coefficeints[:, :, 3][valid], coefficeints[:, :, 4][valid]) / (2 * pi)
            df0[!, :count] = count[valid]
            df0[!, :rmse] = rmse[valid]
        end

        tmp = tempname(dirname(outfile))
        Arrow.write(tmp, df0::DataFrame)
        mv(tmp, outfile; force=true)
        t5 = time()

        total_time = round(t5 - t1, digits=1);
        printstyled("    -> $(geotile.id) $dem: height change complete $total_time s (read [$(round(Int16, (t2-t1)/total_time*100))%], model setup [$(round(Int16,(t3-t2)/total_time*100))%], model [$(round(Int16,(t4-t3)/total_time*100))%], save [$(round(Int16,(t5-t4)/total_time*100))%])\n"; color=:light_black)
    end
end


"""
    geotile_merge_height(geotile, dem, paths; missions, altim_dem_outlier, t_minmax, mask)
Return a merged dataframe of heights for given criteria
"""
function geotile_merge_height(
    geotile,
    dem, 
    paths; 
    missions = [:gedi, :icesat, :icesat2],
    altim_dem_outlier=Inf,
    t_minmax=nothing,
    mask=nothing
    )

    # read altimetry dataframes
    df = DataFrame()
    df0 = DataFrame()
    df1 = DataFrame()

    for mission in missions

        if isfile(joinpath(paths[mission].geotile, "$(geotile.id).arrow")) &&
            isfile(joinpath(paths[mission].geotile, "$(geotile.id).$dem"))

            df = append!(df, geotile_aggrigate_reduce(geotile, paths[mission].geotile, 
                ["datetime", "longitude", "latitude", "height"]; extension = ".arrow")); 

            df0 = append!(df0, geotile_aggrigate_reduce(geotile, paths[mission].geotile, 
                ["height"]; extension = ".$dem")) # 1 min for >50N

            df1 = append!(df1, geotile_aggrigate_reduce(geotile, paths[mission].geotile,
                [mask]; extension=".masks")) 
        end
    end

    if isempty(df0) || isempty(df0)
        printstyled("    -> $(geotile.id) no sensor -or- $(dem) data , skipping\n"; color = :light_red)
        return nothing
    end

    df[!, :height_ref] = df0.height;
    df[!, mask] = df1[!, mask];

    # calculate decimal year
    df[!, :decyear] = decimalyear.(df.datetime)

    # only retain valid points
    df[!, :dh] = df.height .- df.height_ref;
    valid = abs.(df.dh) .< altim_dem_outlier;

     # find data within time range
    if !isnothing(t_minmax)
        valid = valid .& (df.decyear .>= t_minmax[1]) .& (df.decyear .<= t_minmax[2])
    end

    if !isnothing(mask)
        valid = valid .& df[!,mask]
    end

    if !any(valid)
        printstyled("    -> no $(geotile.id) $(dem) data within outlier threshold, skipping\n"; color = :light_red)
        return nothing
    end
    
    delete!(df, .!valid)

    return df
end