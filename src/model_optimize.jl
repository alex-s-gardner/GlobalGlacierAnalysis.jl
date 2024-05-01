function elevation_change(p, t, dz; t_intercept=2015)
    p[1] + p[2] * (t - t_intercept) + p[3] * (t - t_intercept)^2 + p[4] * cos(2 * pi * t) + p[5] * sin(2 * pi * t) + p[6] * dz
end


function loss(p, v)
    err = 0.0
    for i in eachindex(v[1])
        pred_i = Altim.elevation_change(p, v[1][i], v[2][i])
        err += v[3][i] * HuberLoss()(pred_i, v[4][i])
    end
    return err
end

""" 
    TSModelParameters
Defines TSModelParameters parameter type
"""
@kwdef struct OptimModelParameters

    # model parameters
    model::Function = elevation_change
    optf::Function = OptimizationFunction(loss, Optimization.AutoForwardDiff())

    coeff_initial::AbstractVector = [0, -1.0, -1.0, 0.5, 0.0, 0.1]
    coeff_name = [:offset, :trend, :accel, :cos365d, :sin365d, :dhdz]
    iterations::Real = 2
    iter_thresh::Real = 7

    # data filtering parameters
    count_min::Int = 3                      # minimum number of valid values for model fit
    t_min_range::Real = 0                   # only include model in range(t) > t_min_range
    dh_max_outlier::Real = Inf               # maximum acceptable absolute deviation from median
    t_bin_edges = 2000:1/6:2025             # bin edges used for local filtering of outliers
    t_bin_filter_method::Function = madnorm  # method for identifying local outliers
    t_bin_filter_threshold = 10             # threshold for classifying local outliers

    # data loading parameters
    altim_dem_outlier::Real = Inf
    t_minmax = nothing
    mask = nothing
    crop2mask::Bool = false
end


# see wave model examples here: https://www.geogebra.org/m/bmyurtag

"""
    ts_fit(
        x, y, h, t, dh; w, p 
    )

Fits a timeseries model, defined by parameters `p`, to `z` with optional weights `w`.

"""
function ts_fit(
    x::AbstractArray,                           # x coodinate
    y::AbstractArray,                           # y coodinate
    h::AbstractArray,                           # y coodinate
    t::AbstractArray,                           # time in decimal years
    dh::AbstractArray;                           # elevation anomaly
    w::Union{AbstractArray, Nothing} = nothing, # weight of elevation anomaly
    p::OptimModelParameters=OptimModelParameters()  # model parameters
)

    # initialize default outputs
    x_out = NaN;
    y_out = NaN;
    count = Int32(0);
    coeff = fill(NaN, size(p.coeff_initial));
    rmse = NaN;

    # remove NaNs
    valid = .!isnan.(dh)
    x = x[valid];
    y = y[valid];
    t = t[valid];
    h = h[valid]
    dh = dh[valid];

    if (length(x) >= p.count_min) && (Altim.range(t) >= p.t_min_range)
        # global centering 
        dh = dh .- median(dh)

        # remove gross outliesrs
        outlier = (abs.(dh) .> p.dh_max_outlier)

        # using a binned filter, center data and remove outliers
        
        # bin filtering turned off
        #outlier[.!outlier] = binnedfiltering(t[.!outlier], dh[.!outlier], 
        #    p.t_bin_edges; method=p.t_bin_filter_method, threshold=p.t_bin_filter_threshold)

        # iterate model fit
        res = Float64[]
        h_out = 0.
        for i in p.iterations
            # remove any outliers that are > max_outlier m from the centered data
            count = Int32(sum(.!outlier))

            if count < p.count_min 
                return x_out, y_out, count, coeff, rmse
            end

            x_out = mean(x[.!outlier])
            y_out = mean(y[.!outlier])
            h_out = mean(h[.!outlier])

            prob = OptimizationProblem(p.optf, p.coeff_initial, [t[.!outlier], h[.!outlier] .- h_out, w[.!outlier], dh[.!outlier]])
            sol = solve(prob, Optim.NelderMead())
           
            pred = elevation_change.(Ref(sol.u), t[.!outlier], h[.!outlier])
            res = dh[.!outlier] .- pred

            if i<p.iterations
                outlier[.!outlier] = outlier[.!outlier] .| (madnorm(res) .> p.iter_thresh)
            else
                coeff = sol.u
            end
        end
        println("test")
        
        rmse = (mean(res .^ 2)) .^ 0.5
    end

    return x_out, y_out, h_out, count, coeff, rmse
end


"""
    geotile_ts_fit(geotile, dem, paths, grid; products, p, force_remake)
Fit time series model on a grid
"""
function geotile_ts_fit(
    geotile::Union{DataFrame, DataFrameRow},
    dem, 
    paths::NamedTuple, 
    grid::NamedTuple;
    products=project_products(; project_id=:v01),
    p::OptimModelParameters=OptimModelParameters(), # model parameters
    force_remake = false
)

    t1 = time()
    outfile = joinpath(paths.height_change, "$(geotile.id).$(dem)")
    
    if isfile(outfile) && !force_remake
        printstyled("    -> model fit to $(geotile.id) $(dem) already exists, skipping\n"; color = :green)
        return
    end

    df = geotile_merge_height(geotile, dem, paths; products=products, altim_dem_outlier=p.altim_dem_outlier, t_minmax=p.t_minmax, mask=p.mask, crop2mask=p.crop2mask)
  
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
        ind = sortperm(df.Y)
        y = df.Y[ind]
        x = df.X[ind]
        h = df.height[ind]
        dh = df.dh[ind]
        w = 1 ./ df.error[ind]
        decyear = df.decyear[ind]

        coefficeints = fill(NaN, length(griddef.x_node_center), 
            length(griddef.y_node_center), length(p.coeff_initial))
        x_cord = fill(NaN, length(griddef.x_node_center), length(griddef.y_node_center))
        y_cord = copy(x_cord)
        height = copy(x_cord)
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

            Int(sum(x_ind)) >= p.count_min || continue

            x0 = x[x_ind]
            y0 = y[x_ind]

            dh0 = dh[x_ind]
            h0 = h[ind]
            w0 = w[x_ind]
            decyear0 = decyear[x_ind]

            # should we continue ?
            Altim.range(decyear0) >= p.t_min_range || continue

            for j in eachindex(griddef.y_node_center)
                y_min = griddef.y_node_center[j] - griddef.node_half_width
                y_max = griddef.y_node_center[j] + griddef.node_half_width

                start = searchsortedfirst(y0, y_min)
                stop = searchsortedlast(y0, y_max)

                (stop-start) >= p.count_min || continue
                Altim.range(@view(decyear0[start:stop])) >= p.t_min_range || continue

                x_cord[i, j], y_cord[i, j], height[i, j], count[i, j], coefficeints[i, j, :], rmse[i, j] =
                    ts_fit(@view(x0[start:stop]),  @view(y0[start:stop]), @view(h0[start:stop]),
                    @view(decyear0[start:stop]), @view(dh0[start:stop]);
                    w = @view(w0[start:stop]), p = p)
            end
        end
        t4 = time()
        
        # place into new DataFrame
        begin
            df0 = DataFrame()
            valid = .!isnan.(coefficeints[:, :, 1])
            df0[!, :longitude], df0[!, :latitude] = epsg2epsg(x_cord[valid], y_cord[valid], epsg, "EPSG:4326", parse_output=true)
            df0[!, :height] = height[valid]
            for i = eachindex(p.coeff_name)
                df0[!, p.coeff_name[i]] = coefficeints[:, :, i][valid]
            end
            #df0[!, :offset] = coefficeints[:, :, 1][valid]
            #df0[!, :trend] = coefficeints[:, :, 2][valid]
            #df0[!, :amplitude] = (coefficeints[:, :, 3][valid] .^ 2 .+ coefficeints[:, :, 4][valid] .^ 2) .^ 0.5
            #df0[!, :phase] = atan.(coefficeints[:, :, 3][valid], coefficeints[:, :, 4][valid]) / (2 * pi)
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
    geotile_merge_height(geotile, dem, paths; products, altim_dem_outlier, t_minmax, mask)
Return a merged dataframe of heights for given criteria
"""
function geotile_merge_height(
    geotile,
    dem, 
    paths; 
    products = project_products(; project_id=:v01),
    altim_dem_outlier=Inf,
    t_minmax = nothing,
    mask = nothing,
    crop2mask = true,
    )

    # read altimetry dataframes
    df = DataFrame()
    df0 = DataFrame()
    df1 = DataFrame()
    df2 = DataFrame()

    cnt = zeros(Int64,length(products))
    for (i, product) in enumerate(products)
        if isfile(joinpath(paths[product.mission].geotile, "$(geotile.id).arrow")) &&
            isfile(joinpath(paths[product.mission].geotile, "$(geotile.id).$dem"))

            if product.apply_quality_filter
                valid = geotile_aggrigate_reduce(geotile, paths[product.mission].geotile,
                    [:quality]; extension=".arrow").quality
            else
                valid = [true]
            end

            if all(valid)
                foo = geotile_aggrigate_reduce(geotile, paths[product.mission].geotile, 
                    [:datetime, :longitude, :latitude, :height, :height_reference]; extension=".arrow")
                df = append!(df, foo)

                df0 = append!(df0, geotile_aggrigate_reduce(geotile, paths[product.mission].geotile, 
                    ["height"]; extension = ".$dem")) # 1 min for >50N

                if !isnothing(mask)
                    df1 = append!(df1, geotile_aggrigate_reduce(geotile, paths[product.mission].geotile,
                        [mask]; extension=".masks"))
                end

                if product.coregister
                    df2 = append!(df2, geotile_aggrigate_reduce(geotile, paths[product.mission].geotile,
                        [:dx, :dy, :dz, :count, :mad_offset, :mad_ref, :ddh]; extension=".offset2$dem"))
                else
                    z =  zeros(nrow(foo));
                    df2 = append!(df2, DataFrame(dx=z, dy=z, dz=z, count=z, mad_offset=z, mad_ref=z, ddh=z))
                end
            else
                foo = geotile_aggrigate_reduce(geotile, paths[product.mission].geotile,
                    [:datetime, :longitude, :latitude, :height, :height_reference]; extension=".arrow")[valid, :]

                df = append!(df, foo)

                df0 = append!(df0, geotile_aggrigate_reduce(geotile, paths[product.mission].geotile, 
                    ["height"]; extension=".$dem")[valid, :]) # 1 min for >50N

                if !isnothing(mask)
                    df1 = append!(df1, geotile_aggrigate_reduce(geotile, paths[product.mission].geotile,
                        [mask]; extension=".masks")[valid, :])
                end

                if product.coregister
                    df2 = append!(df2, geotile_aggrigate_reduce(geotile, paths[product.mission].geotile,
                        [:dx, :dy, :dz, :count, :mad_offset, :mad_ref, :ddh]; extension=".offset2$dem")[valid, :])
                else
                    z = zeros(nrow(foo))
                    df2 = append!(df2, DataFrame(dx=z, dy=z, dz=z, count=z, mad_offset=z, mad_ref=z, ddh=z))
                end

            end
        end
        cnt[i] = nrow(df)
    end

    cnt = vcat(cnt[1], cnt[2:end] .- cnt[1:end-1])

    product_id = vcat([fill(products[i].id, cnt[i]) for i in eachindex(cnt)]...);
    product_error = vcat([fill(products[i].error_sigma, cnt[i]) for i in eachindex(cnt)]...)
    
    if isempty(df) || isempty(df0)
        printstyled("    -> $(geotile.id) no sensor -or- $(dem) data , skipping\n"; color = :light_red)
        return nothing
    end
    
    #@warn("DEM heights are not being used")
    df = rename!(df, :height_reference => :height_reference0)
    df[!, :height_reference] = df0.height
    
    if !isnothing(mask)
        df[!, mask] = df1[!, mask];
    end

    df = hcat(df,df2);

    df[!, :product] = product_id;
    df[!, :error] = product_error;

    # calculate decimal year
    df[!, :decyear] = decimalyear.(df.datetime)

    # only retain valid points
    df[!, :dh] = df.height .- df.height_reference;
    valid = (abs.(df.dh) .< altim_dem_outlier)

     # find data within time range
    if !isnothing(t_minmax)
        valid = valid .& (df.decyear .>= t_minmax[1]) .& (df.decyear .<= t_minmax[2])
    end

    if !isnothing(mask) & crop2mask
        valid = valid .& df[!,mask]
    end

    if !any(valid)
        printstyled("    -> no $(geotile.id) $(dem) data within outlier threshold, skipping\n"; color = :light_red)
        return nothing
    end
    
    df = delete!(df, .!valid)
    return df
end