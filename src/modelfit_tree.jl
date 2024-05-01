
""" 
    TSModelParameters
Defines TSModelParameters parameter type
"""
@kwdef struct TSModelParameters
    #model::Function = model(t, z; t_intercept=2010) = hcat(ones(size(t)), (t .- t_intercept), cos.(2 * pi * t), sin.(2 * pi * t), (z .- mean(z)))
    #coeff_name = [:offset, :trend, :cos365d, :sin365d, :dhdz]
    #model::Function = model(t, z; t_intercept=2010) = hcat(ones(size(t)), (t .- t_intercept), (t .- t_intercept).^2, cos.(2 * pi * t), sin.(2 * pi * t), (z .- mean(z)))
    #coeff_name = [:offset, :trend, :accel, :cos365d, :sin365d, :dhdz]
    model::Function = model(t, z; t_intercept=2010) = hcat(ones(size(t)), (t .- t_intercept), cos.(2 * pi * t), sin.(2 * pi * t))
    coeff_name = [:offset, :trend, :cos365d, :sin365d]
    #model::Function = model(t, z; t_intercept=2010) = hcat(ones(size(t)), (t .- t_intercept))
    #coeff_name = [:offset, :trend]

    regressor = HuberRegression(fit_intercept=false, scale_penalty_with_samples=false)
    #regressor = RobustRegression(fit_intercept=false, scale_penalty_with_samples=false)
    iterations::Real = 1
    iter_thresh::Real = 5
    
    # data filtering parameters
    count_min::Int = 11                     # minimum number of valid values for model fit
    t_min_range::Real = 10                  # only include model in range(t) > t_min_range
    # dh_max_outlier::Real = 120              # maximum acceptable absolute deviation from median
    t_bin_width = 1/365.25                  # bin width used for local filtering of outliers
    t_bin_mad_threshold = 10                # threshold for classifying local outliers

    # data loading parameters
    altim_dem_outlier::Real = Inf
    t_minmax = nothing
    mask = nothing                          # include mask with data [e.g. :landice]
    crop2mask::Bool = false                 # mask observations
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
    dh::AbstractArray;                          # elevation anomaly
    w::Union{AbstractArray, Nothing} = nothing, # weight of elevation anomaly
    p::TSModelParameters = TSModelParameters(), # model parameters
)

    # initialize default outputs
    x_out = NaN;
    y_out = NaN;
    h_out = NaN;
    count = Int32(0);
    θ = fill(NaN, size(p.coeff_name))
    mad_out = NaN;

    # remove NaNs
    valid = .!isnan.(dh)
    x = x[valid];
    y = y[valid];
    t = t[valid];
    h = h[valid]
    dh = dh[valid];
    w = w[valid];

    if (length(x) >= p.count_min) && (Altim.range(t) >= p.t_min_range)
        # global centering 
        #dh = dh .- median(dh)

        # remove gross outliesrs
        #valid = (abs.(dh) .< p.dh_max_outlier)

        #if !any(valid)
        #    return x_out, y_out, h_out, count, θ, mad_out
        #end

        # there is an issue with tracks with lots of data that are huge outliers 
        # (e.g. an icesat2 track with cloud). One way to try and identify these tracks 
        # is to group by time of aquisition, then look for outliers.
        t_bins = floor(minimum(t)):(p.t_bin_width):(ceil(maximum(t) + p.t_bin_width))
        df = binstats(DataFrame(; t, dh), :t, t_bins, :dh, col_function = [median])


        tb = BinStatistics.bincenter.(df.t)
        Xb = hcat(ones(nrow(df)), tb)
        θb = fit(HuberRegression(fit_intercept=false, scale_penalty_with_samples=false), Xb, df.dh_median)
        dh_med = median(abs.(df.dh_median .- (Xb * θb)))

        # only include detrended values with < p.t_bin_mad_threshold
        X = hcat(ones(length(t)), t)
        dh_pred = X * θb
        valid = (abs.((dh .- dh_pred) ./ dh_med) .< p.t_bin_mad_threshold)

        if !any(valid)
            return x_out, y_out, h_out, count, θ, mad_out
        end

        # iterate model fit
        res = Float64[]
    
        for i in p.iterations
            count = Int32(sum(valid))

            if count < p.count_min 
                return x_out, y_out, h_out, count, θ, mad_out
            end

            x_out = mean(x[valid])
            y_out = mean(y[valid])
            h_out = mean(h[valid])

            # create design matrix
            X = p.model(t[valid], h[valid])

            if isnothing(w)
                θ = fit(p.regressor, X, dh[valid])
            else
                θ = fit(p.regressor, X .* w[valid], dh[valid] .* w[valid])
            end

            res = (dh[valid] .- (X * θ))

            if i<p.iterations
                outlier[valid] = outlier[valid] .| (madnorm(res) .> p.iter_thresh)
            end
        end
        count = length(res)
        mad_out = mad(res)
    end
    
    return x_out, y_out, h_out, count, θ, mad_out
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
    p::TSModelParameters = TSModelParameters(), # model parameters
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

        # build tree
        tree = KDTree(vcat(df.X', df.Y'); leafsize=1000)

        coefficeints = fill(NaN, length(griddef.x_node_center), 
            length(griddef.y_node_center), length(p.coeff_name))

        x_cord = fill(NaN, length(griddef.x_node_center), length(griddef.y_node_center))
        y_cord = copy(x_cord)
        height = copy(x_cord)
        mad_out = copy(x_cord)
        count = zeros(Int32, length(griddef.x_node_center), length(griddef.y_node_center))

        t3 = time()
        Threads.@threads for i in eachindex(griddef.x_node_center)
        #ProfileView.@profview  for i in eachindex(griddef.x_node_center)
            idxs = inrange(tree, vcat(griddef.x_node_center[i] * ones(1,length(griddef.y_node_center)), griddef.y_node_center'), griddef.node_half_width)
            !all(isempty.(idxs)) || continue
            for j in eachindex(griddef.y_node_center)

                # should we continue ?
                length(idxs[j]) >= p.count_min || continue

                # should we continue ?
                Altim.range(df.decyear[idxs[j]]) >= p.t_min_range || continue

                @inbounds x_cord[i, j], y_cord[i, j], height[i, j], count[i, j], coefficeints[i, j, :], mad_out[i, j] =
                    ts_fit(df.X[idxs[j]], df.Y[idxs[j]], df.height_reference[idxs[j]], df.decyear[idxs[j]], df.dh[idxs[j]]; w=1 ./ df.error[idxs[j]], p=p)
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
            df0[!, :count] = count[valid]
            df0[!, :mad] = mad_out[valid]
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
        path2h = joinpath(paths[product.mission].geotile, "$(geotile.id).arrow");
        path2dem = joinpath(paths[product.mission].geotile, "$(geotile.id).$dem");

        if !isfile(path2h)
            @warn("no height file, skipping: $path2h")
            if i > 1
                 cnt[i] = cnt[i-1]
            end
            print("$(product.mission):[0 pts] ")
            continue
        end

        if !isfile(path2dem)
            @warn("height file without a dem file, skipping: $path2dem")
            if i > 1
                 cnt[i] = cnt[i-1]
            end
            print("$(product.mission):[0 pts] ")
            continue
        end

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
        cnt[i] = nrow(df)

        print("$(product.mission):[$(nrow(foo)) pts] ")
    end
    print("\n")

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