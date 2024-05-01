""" 
    TSModelParameters
Defines TSModelParameters parameter type
"""
@kwdef struct TSModelParameters
    #NOTE: ---> first column of model must be offset and second column must be trend <---    
    model::Function = model(t; t_intercept=2010) = hcat(ones(size(t)), (t .- t_intercept), cos.(2 * pi * t), sin.(2 * pi * t))
    coeff_name::Vector{Symbol} = [:offset, :trend, :cos365d, :sin365d]

    regressor = HuberRegression(fit_intercept=false, scale_penalty_with_samples=false)
    #regressor = RobustRegression(fit_intercept=false, scale_penalty_with_samples=false)
    
    iterations::Int64 = 1
    iter_thresh::Int64 = 5
    
    # data filtering parameters
    count_min::Int64 = 11                     # minimum number of valid values for model fit
    t_min_range::Float64 = 10                  # only include model in range(t) > t_min_range
    # dh_max_outlier::Real = 120              # maximum acceptable absolute deviation from median
    t_bin_width::Float64 = 1 / 365.25                  # bin width used for local filtering of outliers
    t_bin_mad_threshold::Float64 = 10                # threshold for classifying local outliers

    # data loading parameters
    altim_dem_outlier::Float64 = Inf
    t_minmax::Vector{Float64} = [-Inf, Inf]
    mask::Symbol = nothing                          # include mask with data [e.g. :landice]
    crop2mask::Bool = false                 # mask observations
end


# see wave model examples here: https://www.geogebra.org/m/bmyurtag
"""
    ts_fit(
        x, y, h, t, dh; p 
    )

Fits a timeseries model, defined by parameters `p`, to `z` with optional weights `w`.

"""
function ts_fit(
    x::AbstractArray,                           # x coodinate
    y::AbstractArray,                           # y coodinate
    h::AbstractArray,                           # y coodinate
    t::AbstractArray,                           # time in decimal years
    dh::AbstractArray;                          # elevation anomaly
    p::TSModelParameters = TSModelParameters(), # model parameters
)

    if !(p.coeff_name[1] == :offset) || !(p.coeff_name[2] == :trend)
        # this is to avoid the need to remake the design matrix for outlier filtering
        error("first column of model must be offset and second column must be trend")
    end

    # initialize default outputs
    x_out = NaN;
    y_out = NaN;
    h_out = NaN;
    count = Int32(0);
    θ = fill(NaN, size(p.coeff_name))
    mad_out = NaN;

    if (length(x) >= p.count_min) && (Altim.range(t) >= p.t_min_range)
        
        # data points are grouped by time to remove outliers and improve speed

        # there is an issue with tracks with lots of data that are huge outliers 
        # (e.g. an icesat2 track with cloud). One way to try and identify these tracks 
        # is to group by time of aquisition, then look for outliers.

        halfwidth = (p.t_bin_width / 2)::Float64

        # nearly all of the time is spent here 
        # NOTE: increasing bin size has little impact on speed so no need to filter at 
        # daily then fit model at monthly [tested but made process slower]
        t0, idx = binevenly(t::AbstractVector{Float64}, halfwidth)
       
        t0 = t0[.!isempty.(idx)]
        idx = idx[.!isempty.(idx)]

        dh0 = Float32[]
        for ind in idx
            push!(dh0,median(dh[ind]))
        end

        # create design matrix for full model
        X = p.model(t0)

        # fit simple line + offset to identify gross outliers [subset of full model]
        θ = fit(p.regressor, X[:, 1:2], dh0)
        dh0_pred = (X[:, 1:2] * θ)
        dh0_absres = @. abs(dh0 - dh0_pred)
        dh0_mad = median(dh0_absres)
        valid = @. (dh0_absres / dh0_mad) < p.t_bin_mad_threshold

        if !any(valid)
            return x_out, y_out, h_out, count, θ, mad_out
        end

        # fit full model to cleaned data weighted by number of observations
        w = @. sqrt(length(idx[valid]))

        #θ = fit(LinearRegression(fit_intercept=false), @view(X[valid, :]) .* w, @view(dh0[valid]) .* w)
        θ = fit(p.regressor, X[valid, :] .* w, dh0[valid] .* w)
        #θ = fit(HuberRegression(fit_intercept=false, scale_penalty_with_samples=false), @view(X[valid,:]), @view(dh0[valid]))
        dh0_pred = X[valid, :] * θ
        mad_out = median(abs.(dh0[valid] .- dh0_pred))
        idx = reduce(vcat,idx[valid])

        x_out = mean(x[idx])
        y_out = mean(y[idx])
        h_out = mean(h[idx])

        # number of days with data will be more valuable then total number of data points
        count = length(dh0_pred)
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

        # remove rows with NaNs
        df = deleteat!(df,isnan.(df.dh))


        if isnothing(df)
            return 
        end

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
            # for i in eachindex(griddef.x_node_center)
            #ProfileView.@profview  for i in eachindex(griddef.x_node_center)
            idxs = inrange(tree, vcat(griddef.x_node_center[i] * ones(1,length(griddef.y_node_center)), griddef.y_node_center'), griddef.node_half_width)

            !all(isempty.(idxs)) || continue
            
            for (j, idx) in enumerate(idxs)
                
                # should we continue ?
                length(idx) >= p.count_min || continue

                # should we continue ?
                Altim.range(@view(df.decyear[idx])) >= p.t_min_range || continue

                # sort as a funciton of time to speedup ts_fit
                #k = @view(idx[sortperm(@view(decyear[idx]))])
                
                #k = @view(idx)
                x_cord[i, j], y_cord[i, j], height[i, j], count[i, j], coefficeints[i, j, :], mad_out[i, j] =
                    ts_fit(@view(df.X[idx]), @view(df.Y[idx]), @view(df.height_reference[idx]), @view(df.decyear[idx]), @view(df.dh[idx]); p=p)
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

function binevenly(x, halfwidth; reorder = true)
    # create center nodes
    (xmin, xmax) = extrema(x)
    x0 = (floor(xmin / halfwidth) * halfwidth + halfwidth):halfwidth:(ceil(xmax / halfwidth) * halfwidth + halfwidth)

    # build kdtree
    tree = KDTree(x'; leafsize=10, reorder)

    # find index into ranges
    idx = inrange(tree, x0', halfwidth)

    return x0, idx
end