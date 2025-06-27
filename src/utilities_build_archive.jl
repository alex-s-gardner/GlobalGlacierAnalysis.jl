"""
    geotile_search_granules(geotiles, mission, product, version, outgranulefile; rebuild_dataframe=false, after=nothing, page_size=2000)

Find all granules that intersect with specified geotiles and save results as an Arrow table.

# Arguments
- `geotiles`: DataFrame containing geotile information
- `mission`: Mission symbol (e.g., `:icesat`, `:ICESat2`, `:GEDI`)
- `product`: Product identifier for the mission
- `version`: Version of the product
- `outgranulefile`: Path to save the output Arrow file
- `rebuild_dataframe`: Whether to rebuild the dataframe if it already exists (default: false)
- `after`: Optional date filter to only include granules after this date
- `page_size`: Number of results per page for API requests (default: 2000)

# Returns
- Path to the saved Arrow file
"""
function geotile_search_granules(
    geotiles,
    mission,
    product,
    version,
    outgranulefile;
    rebuild_dataframe=false,
    after=nothing,
    page_size=2000
)

    # ensure SpaceLiDAR capitilization ot mission 
    mission = mission2spacelidar(mission)

    pathfile = splitdir(outgranulefile)
    if !isdir(pathfile[1])
        error("$(pathfile[1]) does not exist")
    end

    printstyled("identifying grannules that intersect each geotile\n"; color=:blue, bold=true)

    # get intersecting granules [run threads]
    n = nrow(geotiles)
    if mission == :ICESat2
        granules = Vector{Vector{ICESat2_Granule{product}}}(undef, n)
    elseif mission == :ICESat
        granules = Vector{Vector{ICESat_Granule{product}}}(undef, n)
    elseif mission == :GEDI
        granules = Vector{Vector{GEDI_Granule{product}}}(undef, n)
    else
        error("mission and/or product not recognized")
    end


    # this returns a non descript "Something went wrong: ["An Internal Error has occurred."]"
    # ------------------- IF YOU GET THIS ERROR, REMOVE THREADS.@THREADS -------------------
    @showprogress dt = 1 desc = "Finding granules $(mission) $(product)..." Threads.@threads for i in 1:n
        granules[i] = search(mission, product; bbox=geotiles[i, :extent], version=version, after=after)

        #printstyled("    -> finding granules in geotile $i of $n\n"; color=:light_black)
        #foo = true
        #while foo
        #    try
        #        foo = true
        #    catch e
        #        throw(e)
        #        println("error returned from search... wait 1 second and try again")
        #        sleep(1)
        #    end
        #end
    end

    # save remote granules before download in case there is a crash
    geotile_granules = hcat(geotiles, DataFrame(granules=granules))

    # check if file already exists 
    if isfile(outgranulefile) && !rebuild_dataframe
        #if file exists read it in and only overwrite updated granules
        geotile_granules = leftmerge(geotile_granules, granules_load(outgranulefile, mission), :id)
    end

    tmp = tempname(dirname(outgranulefile))
    Arrow.write(tmp, geotile_granules::DataFrame)
    mv(tmp, outgranulefile; force=true)
    return outgranulefile
end

"""
    geotile_download_granules!(
        geotile_granules,
        mission,
        savedir,
        outgranulefile;
        threads=true,
        rebuild_dataframe=false,
        aria2c=false,
        downloadstreams=16,
    )

Downloads satellite data granules for each geotile and updates their URLs to local paths.

# Arguments
- `geotile_granules`: DataFrame containing geotiles and their associated granules
- `mission`: Satellite mission identifier (e.g., :ICESat2, :ICESat, :GEDI)
- `savedir`: Directory where downloaded granules will be saved
- `outgranulefile`: Path to save the updated granule information

# Keywords
- `threads=true`: Whether to use multithreading for downloads
- `rebuild_dataframe=false`: Whether to rebuild the dataframe from scratch
- `aria2c=false`: Whether to use aria2c for downloading
- `downloadstreams=16`: Number of download streams when using aria2c

# Returns
- Path to the saved granule file
"""
function geotile_download_granules!(
    geotile_granules,
    mission,
    savedir,
    outgranulefile;
    threads=true,
    rebuild_dataframe=false,
    aria2c=false,
    downloadstreams=16,
)

    printstyled("downloading granules for each geotile\n"; color=:blue, bold=true)

    # ensure SpaceLiDAR capitilization ot mission 
    mission = mission2spacelidar(mission)

    # remove empty granules
    geotile_granules = geotile_granules[.!isempty.(geotile_granules.granules), :]

    filesOnDisk = readdir(savedir)

    n = size(geotile_granules, 1)
    for (i, row) in enumerate(eachrow(geotile_granules))
        t1 = time()
        printstyled("    -> downloading granules $i of $n ... "; color=:light_black)

        # download seems to get killed when it makes too many requests... try just requesting files that actually need downloading. 
        files2download = [g.id for g in row.granules]
        ia, _ = intersectindices(files2download, filesOnDisk; bool=true)
        ia = .!ia

        # download seems to get killed when it makes too many requests... try just requesting files that actually need downloading. 
        if any(ia)
            granules = row.granules[ia]
        else
            print("no new files to download, skipping\n")
            # update granule urls to local paths
            for g in row.granules
                g.url = joinpath(savedir, g.id)
            end
            continue
        end

        # using async is ~4x faster than without
        #TODO: Replace this with Aria2_jll
        if aria2c
            urls = [g.url for g in granules]
            #urls = vcat(urls...)

            fn = tempname()
            url_list = write_urls!(fn, urls)

            cmd = `aria2c --max-tries=10 --retry-wait=1 -x $downloadstreams -k 1M -j 1 -i --max-connection-per-server 15 -c -d $savedir -i $url_list`

            println(cmd)
            run(cmd)

            # update granule urls to local paths
            for g in row.granules
                g.url = joinpath(savedir, g.id)
            end
        else
            if threads
                asyncmap(granules; ntasks=10) do g
                    flag = 0
                    while flag == 0
                        try
                            download!(g, savedir)
                            flag = 1
                        catch e
                            println(e)
                            wait_time = 10
                            printstyled("download hand an issue... will try again in $wait_time s\n"; color=:yellow)
                            sleep(wait_time)
                        end
                    end
                end
            else
                for g in granules
                    flag = 0
                    while flag == 0
                        try
                            download!(g, savedir)
                            flag = 1
                        catch e
                            println(e)
                            wait_time = 10
                            printstyled("download hand an issue... will try again in $wait_time s\n"; color=:yellow)
                            sleep(wait_time)
                        end
                    end
                end
            end
        end
        print("done [$(round((time()-t1)/60, digits=1)) min]\n")
    end

    # check if file already exists 
    if isfile(outgranulefile) && !rebuild_dataframe
        geotile_granules0 = copy(granules_load(outgranulefile, mission))
        geotile_granules = leftmerge(geotile_granules, geotile_granules0, :id)
    end

    # save granules with local paths
    tmp = tempname(dirname(outgranulefile))
    Arrow.write(tmp, geotile_granules::DataFrame)
    mv(tmp, outgranulefile; force=true)
    return outgranulefile
end

"""
    granules_load(filename, mission; geotiles=nothing, rebuild_dataframe=false)

Load satellite data granules from an Arrow file with proper type reconstruction.

# Arguments
- `filename`: Path to the Arrow file containing granule data
- `mission`: Mission symbol (e.g., `:ICESat`, `:ICESat2`, `:GEDI`)
- `geotiles`: Optional DataFrame of geotiles to filter granules by ID
- `rebuild_dataframe`: Whether to ignore geotile filtering (default: false)

# Returns
- DataFrame of granules with proper types restored, or `nothing` if no matching granules found
"""
function granules_load(
    filename,
    mission;
    geotiles::Union{Nothing,DataFrame}=nothing,
    rebuild_dataframe=false,
)

    # ensure SpaceLiDAR capitilization ot mission 
    mission = mission2spacelidar(mission)

    # load all granules
    granules = DataFrame(Arrow.Table(filename))
    granules[!, :extent] = [Extent(ext.bounds) for ext in granules.extent]

    # subset if requested
    if !isnothing(geotiles) && !rebuild_dataframe
        ind = [findfirst(granules[!, :id] .== id) for id in geotiles[!, :id]]
        valid = .!isnothing.(ind)
        if any(valid)
            granules = granules[ind[valid], :]
        else
            granules = nothing
            return granules
        end
    end

    # type is stripped on df save, add back
    polygon = Vector{Vector{Vector{Vector{Float64}}}}[]
    if mission == :ICESat2
        granules[!, :granules] =
            [[ICESat2_Granule{Symbol(g.info.type)}(g.id, g.url, g.info, polygon) for g in grans] for grans in granules[!, :granules]]
    elseif mission == :ICESat
        granules[!, :granules] =
            [[ICESat_Granule{Symbol(g.info.type)}(g.id, g.url, g.info, polygon) for g in grans] for grans in granules[!, :granules]]
    elseif mission == :GEDI
        granules[!, :granules] =
            [[GEDI_Granule{Symbol(g.info.type)}(g.id, g.url, g.info, polygon) for g in grans] for grans in granules[!, :granules]]
    else
        error("need to add mission")
    end

    return granules
end

"""
    pointextract(x, y, point_epsg, ga::GeoArray; kwargs...) -> Union{Vector, Vector{Vector}}

Extract values from a GeoArray at specified coordinates, with optional filtering and derivative calculation.

# Arguments
- `x`: X-coordinates (longitude or easting) of points to extract
- `y`: Y-coordinates (latitude or northing) of points to extract
- `point_epsg`: EPSG code string for the coordinate system of input points
- `ga`: GeoArray containing the data to extract from

# Keywords
- `nodatavalue=0.0`: Value to use for points outside the array or at nodata locations
- `replace_nodatavalue_withnans=false`: Whether to replace nodata values with NaN
- `convert_to_float=false`: Whether to convert integer arrays to floating point
- `filter_kernel=nothing`: Optional kernel for filtering the GeoArray before extraction
- `derivative_kernel=nothing`: Optional kernel(s) for calculating derivatives
- `interp=Constant()`: Interpolation method (Constant(), Linear(), Cubic(), etc.)

# Returns
- If no derivatives requested: Vector of extracted values
- If derivatives requested: Vector of vectors containing [values, derivative1, derivative2, ...]

# Performance Notes
- Constant() and Linear() interpolation are significantly faster than Cubic()
- Most computation time is spent in the interpolation step
"""
function pointextract(
    x,
    y,
    point_epsg,
    ga::GeoArray;
    nodatavalue=0.0,
    replace_nodatavalue_withnans=false,
    convert_to_float=false,
    filter_kernel=nothing,
    derivative_kernel=nothing,
    interp::F=Constant()
) where {F<:Interpolations.Flag}

    #################################### PERFORMANCE NOTE ##################################
    # >80% of time is spent on Interpolations when using ubic(Line(OnGrid()))
    # significant speedups when interp = Constant() or Linear()
    # 
    # example: for ICESat lat[+54+56]lon[+096+098].arrow for 4 derivative variables
    #       1. interp = Constant(): 30s
    #       2. interp = Linear(): 30s 
    #       3. interp = Cubic(Line(OnGrid())):  86s
    #       4. interp = Quadratic(Line(OnGrid())): 85s
    ########################################################################################

    # create projeciton transfomation function
    # NOTE: not worth checking for same projection... faster just to do it this way
    begin # 0.06s

        _, epsg_num = split(point_epsg, ":")
        if !Base.contains(ga.crs.val[end-7:end], "\"EPSG\",\"$epsg_num\"")
            xy = epsg2epsg(x, y, point_epsg, ga.crs.val; parse_output=false)
            x = getindex.(xy, 1)
            y = getindex.(xy, 2)
        end

        # loop for multiple derivative
        if !isnothing(derivative_kernel)
            if derivative_kernel isa Tuple
                n_derivative = length(derivative_kernel)
            else
                n_derivative = 1
                derivative_kernel = (derivative_kernel)
            end
        else
            n_derivative = 0
        end

        # exclude points where x or y is NaN ... also include old x = 0, y = 0 until old files are updated. this was used as a no data value in the datacubes
        ind = .!(isnan.(x) .| isnan.(y))

        if !any(ind)
            # extents do not overlap
            if replace_nodatavalue_withnans
                val = fill(convert(eltype(ga.A), NaN), length(x))
            else
                val = fill(convert(eltype(ga.A), nodatavalue), length(x))
            end

            if n_derivative > 0
                val = [val for i = 1:(n_derivative+1)]
            end
            return val
        end

        # find x y extents
        minmax_x = extrema(x[ind])
        minmax_y = extrema(y[ind])
        dx = ga.f.linear[1, 1]
        dy = ga.f.linear[2, 2]

        # deteremine how much to pad crop by
        if !isnothing(derivative_kernel) .| !isnothing(filter_kernel)
            if isnothing(derivative_kernel)
                px = ceil(size(filter_kernel, 2) / 2)
                py = ceil(size(filter_kernel, 1) / 2)
            elseif isnothing(filter_kernel)
                px = 2
                py = 2
            else
                px = ceil(size(filter_kernel, 2) + 2)
                py = ceil(size(filter_kernel, 1) + 2)
            end
        else
            px = 1
            py = 1
        end

        if typeof(interp) <: Linear
            px += 1
            py += 1
        elseif typeof(interp) <: Cubic
            px += 2
            py += 2
        elseif typeof(interp) <: Quadratic
            px += 3
            py += 3
        end

        extent = (
            min_x=minmax_x[1] - abs(dx) * px,
            min_y=minmax_y[1] - abs(dy) * py,
            max_x=minmax_x[2] + abs(dx) * px,
            max_y=minmax_y[2] + abs(dy) * py
        )
    end

    begin # 4s
        if !bbox_overlap(extent2nt(bbox(ga)), extent)

            # extents do not overlap
            if replace_nodatavalue_withnans
                val = fill(convert(eltype(ga.A), NaN), length(x))
            else
                val = fill(convert(eltype(ga.A), nodatavalue), length(x))
            end

            if n_derivative > 0
                val = [val for i = 1:(n_derivative+1)]
            end

            return val
        else
            # crop and read into  memory 
            #@infiltrate

            ## ~65% of all time is spent here ##
            ga0 = GeoArrays.crop(ga, extent)

            # println("size of cropped array = $(size(ga0)), crop to extent = $(extent)")

            ###################################
            if convert_to_float
                if (eltype(ga0.A) <: Integer)
                    f = ga0.f
                    ga0 = GeoArray(Float32.(ga0.A))
                    ga0.crs = ga.crs
                    ga0.f = f
                end
            end


            if replace_nodatavalue_withnans
                if (eltype(ga0.A) <: Integer)
                    f = ga0.f
                    ga0 = GeoArray(Float32.(ga0.A))
                    ga0.crs = ga.crs
                    ga0.f = f
                end

                for i = 1:size(ga0.A, 3)
                    @view(ga0.A[:, :, i])[ga0.A[:, :, i].==nodatavalue] .= NaN
                end
                nodatavalue = NaN
            end

            # apply filter (typically a smoothing filter)
            if !isnothing(filter_kernel)
                m, nd = size(filter_kernel)
                m = Int((m - 1) / 2)
                nd = Int((nd - 1) / 2)
                stencil = Stencils.Rectangle((-m, m), (-nd, nd))
                k = Stencils.Kernel(stencil, filter_kernel)

                for i = 1:size(ga0.A, 3)
                    A = StencilArray(ga0.A[:, :, i], k)
                    mapstencil!(kernelproduct, @view(ga0.A[:, :, i]), A)
                    #@time imfilter!(@view(ga0.A[:,:,i]), ga0.A[:,:,i], filter_kernel);
                end
            end
        end
    end

    begin #16s for Cubic(Line(OnGrid()))
        x0, y0 = GeoArrays.ranges(ga0)
        gridsize = xy_gridsize(ga0)
        y2x_scale = round(abs((gridsize.y / gridsize.x)) / (dy / dx), digits=5)
        val = Matrix{eltype(ga0.A)}(undef, (length(x), size(ga0.A, 3)))

        for i in eachindex(ga0.A[1, 1, :])
            itp = Interpolations.extrapolate(Interpolations.scale(Interpolations.interpolate(ga0.A[:, :, i], BSpline(interp)), x0, y0 .* y2x_scale), nodatavalue)

            val[:, i] = itp.(x, y .* y2x_scale)
        end
    end

    # check if a derivative has been requested
    if !isnothing(derivative_kernel)
        begin # 11s
            ga1 = [copy(ga0) for i = 1:n_derivative]
            for j = 1:n_derivative
                for i = 1:size(ga1[j].A, 3)
                    mapstencil!(derivative_kernel[j], @view(ga1[j].A[:, :, i]), StencilArray(ga0.A[:, :, i], Window(1)))
                end
            end
        end

        begin # 53s for Cubic(Line(OnGrid()))
            deriv = [Matrix{eltype(ga0.A)}(undef, (length(x), size(ga0.A, 3))) .+ nodatavalue for i = 1:n_derivative]
            for j in 1:n_derivative
                for i = 1:size(ga0.A, 3)
                    itp = Interpolations.extrapolate(Interpolations.scale(Interpolations.interpolate(ga1[j].A[:, :, i], Interpolations.BSpline(interp)), x0, y0 .* y2x_scale), nodatavalue)
                    deriv[j][:, i] = itp.(x, y .* y2x_scale)
                end
            end
        end
        val = vcat([val], deriv) #0.002s
    end
    return val
end

"""
    epsg2epsg(x, y, height, from_epsg, to_epsg; parse_output=true, threaded=true)

Convert coordinates with height between different coordinate reference systems.

# Arguments
- `x`: X-coordinates (longitude or easting) as vector or scalar
- `y`: Y-coordinates (latitude or northing) as vector or scalar
- `height`: Height values as vector or scalar
- `from_epsg`: Source EPSG code as string
- `to_epsg`: Target EPSG code as string

# Keywords
- `parse_output=true`: Whether to return separate x,y,height arrays (true) or tuples (false)
- `threaded=true`: Whether to use multithreading for coordinate transformation

# Returns
- If `parse_output=true`: Tuple of (x, y, height) in target coordinate system
- If `parse_output=false`: Vector of coordinate tuples or single tuple
"""
function epsg2epsg(
    x::Union{AbstractVector{<:Any},Number},
    y::Union{AbstractVector{<:Any},Number},
    height::Union{AbstractVector{<:Number},Number},
    from_epsg::String,
    to_epsg::String;
    parse_output=true,
    threaded=true
)
    # is not stable when using @threads

    # enamble network for geoid conversions
    Proj.enable_network!()

    # build transformation 
    if !threaded
        trans = Proj.Transformation(from_epsg, to_epsg; always_xy=true)

        # project points
        if x isa Vector{}
            xyh = trans.(x, y, height)
        else
            xyh = trans(x, y, height)
        end
    else
        # This will work with threads (and you can add your own Proj context in ctxs), but not on the GPU - that would pretty much require a Julia implementation of Proj.
        ctxs = [Proj.proj_context_clone() for _ in 1:Threads.nthreads()]
        transforms = [Proj.Transformation(from_epsg, to_epsg; always_xy=true, ctx) for ctx in ctxs]

        xyh = Vector{Tuple{Float64,Float64,Float64}}(undef, size(x))
        Threads.@threads for i in eachindex(latitude)
            xyh[i] = transforms[Threads.threadid()](x[i], y[i], height[i])
        end
    end

    # parse from coordinate structure
    if parse_output
        if x isa Vector{}
            x = getindex.(xyh, 1)
            y = getindex.(xyh, 2)
            height = getindex.(xyh, 3)
        else
            x = getindex(xyh, 1)
            y = getindex(xyh, 2)
            height = getindex(xyh, 3)
        end

        return x, y, height
    else
        return xyh
    end
end


"""
    epsg2epsg(x, y, from_epsg, to_epsg; parse_output=true, threaded=true)

Convert coordinates between different coordinate reference systems.

# Arguments
- `x`: X-coordinate(s) (longitude or easting) in source projection
- `y`: Y-coordinate(s) (latitude or northing) in source projection
- `from_epsg`: Source EPSG code as string (e.g., "EPSG:4326")
- `to_epsg`: Target EPSG code as string (e.g., "EPSG:3857")

# Keywords
- `parse_output=true`: If true, returns separate x,y vectors; if false, returns tuples
- `threaded=true`: Whether to use multithreading for batch conversions

# Returns
- If `parse_output=true`: Tuple of (x, y) arrays in target projection
- If `parse_output=false`: Vector of coordinate tuples or single tuple
"""
function epsg2epsg(
    x::Union{AbstractVector{<:Any},Number},
    y::Union{AbstractVector{<:Any},Number},
    from_epsg::String,
    to_epsg::String;
    parse_output=true,
    threaded=true
)

    # build transformation 
    if !threaded
        trans = Proj.Transformation(from_epsg, to_epsg; always_xy=true)

        # project points
        if x isa Vector{}
            xy = trans.(x, y)
        else
            xy = trans(x, y)
        end
    else
        # This will work with threads (and you can add your own Proj context in ctxs), but not on the GPU - that would pretty much require a Julia implementation of Proj.
        ctxs = [Proj.proj_context_clone() for _ in 1:Threads.nthreads()]
        transforms = [Proj.Transformation(from_epsg, to_epsg; always_xy=true, ctx) for ctx in ctxs]

        xy = Vector{Tuple{Float64,Float64}}(undef, size(x))
        Threads.@threads for i in eachindex(xy)
            xy[i] = transforms[Threads.threadid()](x[i], y[i])
        end
    end

    if parse_output
        if x isa Vector{}
            x = getindex.(xy, 1)
            y = getindex.(xy, 2)
        else
            x = getindex(xy, 1)
            y = getindex(xy, 2)
        end
        return x, y
    else
        return xy
    end
end

"""
    dem_height(lon, lat, dem; filter_halfwidth=nothing, filter_kernel=:gaussian, slope=false, curvature=false)

Extract elevation values from a digital elevation model (DEM) at specified coordinates.

# Arguments
- `lon`: Longitude coordinates (degrees)
- `lat`: Latitude coordinates (degrees)
- `dem`: Digital elevation model identifier

# Keywords
- `filter_halfwidth`: Half-width of filter to apply to DEM before sampling (in distance units)
- `filter_kernel`: Type of filter kernel to use (`:gaussian` or `:average`)
- `slope`: Whether to calculate surface slope (default: false)
- `curvature`: Whether to calculate surface curvature (default: false)

# Returns
- If `slope=false` and `curvature=false`: Vector of elevation values in WGS84 datum
- Otherwise: Vector of vectors containing [elevation, derivatives] where derivatives depend on 
  the combination of `slope` and `curvature` parameters
"""
function dem_height(
    lon,
    lat,
    dem;
    filter_halfwidth::Union{Number,Nothing}=nothing,
    filter_kernel=:gaussian,
    slope=false,
    curvature=false,
)

    begin # 0.1s
        # use gaussian sampeling filter only if the center pixel contains < filter_thresh of the signal
        filter_thresh = 0.85
        interp = Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid())) # Constant(), Linear(), Cubic(Line(OnGrid())), Quadratic()
        deminfo, goids_folder = dem_info(dem)

        # exclude points where x or y is NaN ... also include old x = 0, y = 0 until old files are updated. this was used as a no data value in the datacubes
        ind = .!(isnan.(lon) .| isnan.(lat)) 

        # read in dem
        dem_ga = GeoArrays.read(deminfo.fn; masked=false)

        # check if within extents
        if (deminfo.extent.min_x != -180.0) ||
        (deminfo.extent.min_y != -90.0) ||
        (deminfo.extent.max_x != 180.0) ||
        (deminfo.extent.max_y != 90.0)

            ind = within.(Ref(deminfo.extent), lon, lat) .& ind
        end

        if !any(ind)
            fill_val = convert(eltype(dem_ga.A), NaN)
            n = 1 + 2 * slope + 2 * curvature
            m = length(lat)
            wgs84_height = fill(fill_val, m)
            wgs84_height = [wgs84_height for i = 1:n]
            return wgs84_height
        end

        # extract data from DEMs

        # determine approximate x and y grid size
        # check if geographic
        isgeographic = Base.contains(dem_ga.crs.val, "AXIS[\"Latitude\",NORTH]")

        if !isnothing(filter_halfwidth)
            if isgeographic

                (a, b) = extrema(lat[ind])
                (c, d) = extrema(lon[ind])

                cent = (x=(c + d) / 2, y=(a + b) / 2)
                delta = 0.01
                x_lla = LLA(cent.y, cent.x, 0.0)
                y_lla = LLA(cent.y + delta, cent.x, 0.0)
                dist_lat = euclidean_distance(x_lla, y_lla) / delta

                x_lla = LLA(cent.y, cent.x, 0.0)
                y_lla = LLA(cent.y, cent.x + delta, 0.0)
                dist_lon = euclidean_distance(x_lla, y_lla) / delta

                dist = (x=dem_ga.f.linear[1, 1] * dist_lon, y=dem_ga.f.linear[2, 2] * dist_lat)
            else
                dist = (x=dem_ga.f.linear[1, 1], y=dem_ga.f.linear[2, 2])
            end

            if filter_kernel == :gaussian
                # build gaussain sampeling kernel
                k = gaussian([abs(filter_halfwidth / dist.x), abs(filter_halfwidth / dist.y)]) # for GeoArrays x is rows
                k = k[1] * ones(1, length(k[2])) .+ ones(length(k[1]), 1) * k[2]'
                k = k / sum(k)

            elseif filter_kernel == :average
                #round to the nearest odd integer
                dr = round(Int16, abs((filter_halfwidth) / dist.x) * 2 + 1) # for GeoArrays x is rows
                dc = round(Int16, abs((filter_halfwidth) / dist.y) * 2 + 1) # for GeoArrays y is columns
                k = ones(dr, dc) ./ (dr * dc)

            end

            # only apply kernel if needed
            cind = Int.((size(k) .- 1) ./ 2)
            if k[cind[1], cind[2]] > filter_thresh
                k = nothing
            end
        else
            k = nothing
        end
    end

    begin
        if !slope
            if deminfo.hight_datum == "wgs84"
                wgs84_height = pointextract(lon, lat, "EPSG:4326", dem_ga;
                    nodatavalue=deminfo.nodatavalue, filter_kernel=k, interp=interp,
                    replace_nodatavalue_withnans=false, convert_to_float=true)
            else
                # convert from geoid to wgs84 height
                geoid_height = pointextract(lon, lat, "EPSG:4326", dem_ga;
                    nodatavalue=deminfo.nodatavalue, filter_kernel=k, interp=interp,
                    replace_nodatavalue_withnans=false, convert_to_float=true)
                geoid_to_wgs84 = geoid(deminfo.hight_datum; folder=goids_folder)
                geoid_to_wgs84 = pointextract(lon, lat, "EPSG:4326", geoid_to_wgs84;
                    nodatavalue=deminfo.nodatavalue, filter_kernel=k, interp=interp,
                    replace_nodatavalue_withnans=false, convert_to_float=true)
                wgs84_height = geoid_height .+ geoid_to_wgs84

            end
            return wgs84_height
        else
            gsx = dem_ga.f.linear[1, 1] # for GeoArrays x is rows
            gsy = dem_ga.f.linear[2, 2] # for GeoArrays y is columns

            # slope kernels
            if slope
                dx = let X = gsx
                    v -> _dx(v, X)
                end

                dy = let X = gsy
                    v -> _dy(v, X)
                end

            end

            if curvature
                ddx = let X = gsx
                    v -> _ddx(v, X)
                end

                ddy = let X = gsy
                    v -> _ddy(v, X)
                end
            end

            if slope && !curvature
                derivative_kernel = (dx, dy)
            elseif !slope && curvature
                derivative_kernel = (ddx, ddy)
            elseif slope && curvature
                derivative_kernel = (dx, dy, ddx, ddy)
            end

            if deminfo.hight_datum == "wgs84"
                    wgs84_height_slope = pointextract(lon, lat, "EPSG:4326", dem_ga;
                    nodatavalue=deminfo.nodatavalue, filter_kernel=k,
                    derivative_kernel=derivative_kernel, replace_nodatavalue_withnans=false, 
                    interp=interp, convert_to_float=true)
            else
                # convert from geoid to wgs84 height
                geoid_height = pointextract(lon, lat, "EPSG:4326", dem_ga;
                    nodatavalue=deminfo.nodatavalue, filter_kernel=k,
                    derivative_kernel=derivative_kernel, replace_nodatavalue_withnans=false, 
                    interp=interp, convert_to_float=true)

                geoid_to_wgs84 = geoid(deminfo.hight_datum; folder=goids_folder)
                geoid_to_wgs84 = pointextract(lon, lat, "EPSG:4326", geoid_to_wgs84;
                    nodatavalue=deminfo.nodatavalue, filter_kernel=k, interp=interp, convert_to_float=true)

                geoid_height[1] = geoid_height[1] .+ geoid_to_wgs84
                wgs84_height_slope = geoid_height
            end
        end
    end    
    return wgs84_height_slope
end

"""
    geotile_extract_dem(geotile_id, geotile_dir, dem; kwargs...) -> Nothing

Extract DEM (Digital Elevation Model) data for a specific geotile and save as an Arrow file.

# Arguments
- `geotile_id::String`: Identifier for the geotile
- `geotile_dir::String`: Directory containing geotile files
- `dem::Symbol`: DEM source to extract from (e.g., `:copernicus`, `:arcticdem`)

# Keywords
- `filter_halfwidth::Union{Number,Nothing}=nothing`: Half-width of filter to apply to DEM
- `filter_kernel=:gaussian`: Type of filter kernel to use
- `slope=false`: Whether to calculate surface slope
- `curvature=false`: Whether to calculate surface curvature
- `job_id=""`: Identifier for logging purposes
- `xoffset::Real=0`: Offset to add to x-coordinates in local UTM [m]
- `yoffset::Real=0`: Offset to add to y-coordinates in local UTM [m]
- `force_remake=false`: Whether to regenerate existing files

# Returns
Nothing, but creates a DEM file for the geotile with extracted elevation data
"""
function geotile_extract_dem(
    geotile_id::String,
    geotile_dir::String,
    dem::Symbol;
    filter_halfwidth::Union{Number,Nothing}=nothing,
    filter_kernel=:gaussian,
    slope=false,
    curvature=false,
    job_id="",
    xoffset::Real=0, # add offset to altimetry locations in local UTM coordinates [m]
    yoffset::Real=0, # add offset to altimetry locations in local UTM coordinates [m]
    force_remake=false
)

    path2geotile = joinpath(geotile_dir, geotile_id * ".arrow")
    outfile = replace(path2geotile, ".arrow" => ".$dem")

    if isfile(outfile) && !force_remake
        printstyled("\n    ->$job_id $geotile_id $dem already exists, skipping\n"; color=:green)
    elseif isfile(path2geotile)
        t1 = time()

        begin #0.002s
        # load dem info
        deminfo, _ = dem_info(dem)
        
        # check if geotile extents intersect DEM extents 
        if !(Extents.intersects(nt2extent(deminfo.extent), GeoTiles.extent(geotile_id)))
            printstyled("\n    ->$job_id $geotile_id $dem outside of dem limits, skipping\n"; color=:light_red)
            return
        end

        # load lat an lon from geotiles
        df = select!(DataFrame(Arrow.Table(path2geotile)), [:longitude, :latitude])
        end

        #if isempty(df)
        #    printstyled("    ->$job_id $geotile_id had an empty DataFrame, deleting\n"; color=:red, bold=true)
        #    rm(path2geotile)
        #    return
        #end

        # remove empty rows
        # df = df[.!isnan.(df.longitude),:]

        if isempty(df)
            printstyled("\n    ->$job_id $geotile_id had an empty DataFrame, skipping\n"; color=:grey, bold=true)
            return
        end

        begin #0.00002s
            #  stat(infile).size
            if isa(df.longitude[1], AbstractVector)
                lon = reduce(vcat, df[:, :longitude])
                lat = reduce(vcat, df[:, :latitude])
            else
                lon = @view df[:, :longitude];
                lat = @view df[:, :latitude];
            end
        end

        # this code is for  x and y offset testing
        if !isequal(0, xoffset) || !isequal(0, yoffset)
            epsg = utm_epsg(mean(lon), mean(lat); always_xy=true)
            trans = Proj.Transformation("EPSG:4326", epsg, always_xy=true)
            xy = trans.(lon, lat)
            xy = inv(trans).((getindex.(xy, 1) .+ xoffset), (getindex.(xy, 2) .+ yoffset))
            lon = getindex.(xy, 1)
            lat = getindex.(xy, 2)
        end

        if isempty(lon)
            printstyled("\n    ->$job_id $geotile_id is all empty, skipping\n"; color=:light_red)
            return
        end

        # check if within extents
        if (deminfo.extent.min_x != -180.0) || (deminfo.extent.min_y != -90.0) || (deminfo.extent.max_x != 180.0) || (deminfo.extent.max_y != 90.0)
            ind = within.(Ref(deminfo.extent), lon, lat)
            if !any(ind) 
                printstyled("\n    ->$job_id $geotile_id $dem outside of dem limits, skipping\n"; color=:light_red)
                return
            end
        end

        begin # ALL TIME IS SPENT HERE
            if isnothing(filter_halfwidth)
                height = dem_height(lon, lat, dem, slope, curvature)
            else
                height = dem_height(lon, lat, dem; filter_halfwidth=filter_halfwidth, filter_kernel=filter_kernel, slope, curvature)
            end
        end

        begin # 0.02s
            if !slope && !curvature
                if !isa(df.longitude[1], AbstractVector)
                    df[!, :height] = height
                else
                    df[!, :height] = [Array{eltype(height)}(undef, size(row)) for row in df[:, :longitude]]
                    start = 1
                    for row = eachrow(df)
                        stop = length(row.height) + start - 1
                        row.height = height[start:stop]
                        start = stop + 1
                    end
                end
            else
                if slope && !curvature
                    var_names = (:height, :dhdx, :dhdy)
                elseif !slope && curvature
                    var_names = (:height, :dhddx, :dhddy)
                elseif slope && curvature
                    var_names = (:height, :dhdx, :dhdy, :dhddx, :dhddy)
                end

                if !isa(df.longitude[1], AbstractVector)
                    for (k, v) in enumerate(var_names)
                        df[!, v] = vec(height[k])
                    end
                else
                    for v in var_names
                        df[!, v] = [Array{eltype(height[1])}(undef, size(row)) for row in df[:, :longitude]]
                    end

                    start = 1
                    for row = eachrow(df)
                        stop = length(row.height) + start - 1
                        for (k, v) in enumerate(var_names)
                            row[v] = height[k][start:stop]
                        end
                        start = stop + 1
                    end
                end
            end

            tmp = tempname(dirname(outfile))
            Arrow.write(tmp, df::DataFrame)
            mv(tmp, outfile; force=true)

            total_time = round((time() - t1) / 60, digits=2)
            printstyled("\n    ->$job_id $geotile_id $dem extracted: $(total_time) min \n"; color=:light_black)
        end
    else
        printstyled("\n    ->$job_id $geotile_id does not exist\n"; color=:yellow)
        return
    end
end

"""
    geotile_extract_dem(geotiles, geotile_dir, dem; kwargs...)

Extract DEM (Digital Elevation Model) heights for multiple geotiles in parallel.

# Arguments
- `geotiles::DataFrame`: DataFrame containing geotile information
- `geotile_dir::String`: Directory where geotiles are stored
- `dem::Symbol`: Symbol representing the DEM to extract from (e.g., `:MERIT`, `:COP30`)

# Keywords
- `filter_halfwidth::Union{Number,Nothing}=nothing`: Half-width of filter to apply to DEM
- `filter_kernel=:gaussian`: Type of filter kernel to use
- `slope=false`: Whether to calculate slope derivatives
- `curvature=false`: Whether to calculate curvature derivatives
- `job_id=""`: Identifier for the job (used in progress messages)
- `xoffset=0`: Offset to add to X coordinates in local UTM [meters]
- `yoffset=0`: Offset to add to Y coordinates in local UTM [meters]
- `force_remake=false`: Whether to force reprocessing of existing files
"""
function geotile_extract_dem(
    geotiles::DataFrame,
    geotile_dir::String,
    dem::Symbol;
    filter_halfwidth::Union{Number,Nothing}=nothing,
    filter_kernel=:gaussian,
    slope=false,
    curvature=false,
    job_id="",
    xoffset=0, # add offset to altimetry locations in local UTM coordinates [m]
    yoffset=0, # add offset to altimetry locations in local UTM coordinates [m]
    force_remake=false
)


    @showprogress dt = 10 desc = "extracting $dem heights for $job_id geotiles..." Threads.@threads for geotile in eachrow(geotiles)
    #for geotile in eachrow(geotiles)
        geotile_extract_dem(geotile.id, geotile_dir, dem; 
            filter_halfwidth, filter_kernel, slope, curvature,
            job_id, xoffset, yoffset, force_remake)
    end
end

"""
    geotile_extract_dem(products, dems, geotiles, paths; kwargs...)

Extract DEM data for multiple products and DEM sources across all geotiles.

# Arguments
- `products::NamedTuple`: Collection of products with mission, filter halfwidth, and kernel properties
- `dems::Vector{Symbol}`: List of DEM sources to extract from (e.g., [:MERIT, :COP30])
- `geotiles::DataFrame`: DataFrame containing geotile information
- `paths::NamedTuple`: Paths organized by mission containing geotile directories

# Keywords
- `slope::Bool=false`: Whether to calculate slope derivatives
- `curvature::Bool=false`: Whether to calculate curvature derivatives
- `xoffset::Real=0`: Offset to add to X coordinates in local UTM [meters]
- `yoffset::Real=0`: Offset to add to Y coordinates in local UTM [meters]
- `force_remake::Bool=false`: Whether to force reprocessing of existing files

# Returns
Nothing, but creates DEM files for each product-DEM combination across all geotiles
"""
function geotile_extract_dem(
    products::NamedTuple,
    dems::Vector{Symbol},
    geotiles::DataFrame,
    paths::NamedTuple;
    slope=false,
    curvature=false,
    xoffset=0, # add offset to altimetry locations in local UTM coordinates [m]
    yoffset=0, # add offset to altimetry locations in local UTM coordinates [m]
    force_remake=false
)

    for product in products
        for dem in dems
            geotile_extract_dem(geotiles, paths[product.mission].geotile, dem;
                filter_halfwidth=product.halfwidth, filter_kernel=product.kernel, slope=slope, curvature=curvature,
                job_id=product.mission, xoffset=xoffset, yoffset=yoffset, force_remake=force_remake)
        end
    end
end


"""
    geotile_aggrigate_reduce(geotiles::DataFrame, geotile_dir::String, vars::Vector{String}; 
                            extension=".arrow", subset_fraction=1.0) -> DataFrame

Aggregate and combine variables from multiple geotile files into a single DataFrame.

# Arguments
- `geotiles::DataFrame`: DataFrame containing geotile information
- `geotile_dir::String`: Directory containing geotile files
- `vars::Vector{String}`: Variables to extract from geotile files

# Keywords
- `extension::String=".arrow"`: File extension of geotile files
- `subset_fraction=1.0`: Fraction of data to randomly sample (between 0 and 1)

# Returns
- `DataFrame`: Combined data with requested variables and reduced arrays
"""
function geotile_aggrigate_reduce(
    geotiles::DataFrame,
    geotile_dir::String,
    vars::Vector{String};
    extension::String=".arrow",
    subset_fraction=1.0
)

    # printstyled("aggrigating and reducing $extension geotiles\n"; color = :blue, bold = true)
    df = DataFrame()
    #asyncmap(eachrow(geotiles); ntasks = 4) do geotile
    for geotile in eachrow(geotiles)
        path2geotile = joinpath(geotile_dir, "$(geotile.id)$(extension)")
        if isfile(path2geotile)
            append!(df, select!(DataFrame(Arrow.Table(path2geotile)), vars))
        end
    end

    # random subsample requested?
    ind = eachindex(eachrow(df))
    if subset_fraction < 1
        ind = randsubseq(ind, 0.5)
    end

    # rearrange data into vectors 
    df0 = DataFrame()
    for var in vars
        df0[!, var] = reduct(vcat, df[ind, var])
    end

    # remove empty rows
    delete(df0, isempty.(df0[:, first(vars)]))

    return df0
end

"""
    geotile_aggrigate_reduce(geotiles::DataFrameRow, geotile_dir::String, vars::Union{Vector{String},Vector{Symbol}}; kwargs...) -> DataFrame

Aggregate variables from a single geotile file and reduce array data to vectors of elements.

# Arguments
- `geotiles::DataFrameRow`: Single row from a DataFrame containing geotile information
- `geotile_dir::String`: Directory containing geotile files
- `vars::Union{Vector{String},Vector{Symbol}}`: Variables to extract from geotile file

# Keywords
- `extension::String=".arrow"`: File extension of geotile files
- `subset_fraction=1.0`: Fraction of data to randomly sample (between 0 and 1)

# Returns
- `DataFrame`: Data with requested variables and reduced arrays, or empty DataFrame if file not found
"""
function geotile_aggrigate_reduce(
    geotiles::DataFrameRow{DataFrame,DataFrames.Index},
    geotile_dir::String,
    vars::Union{Vector{String},Vector{Symbol}};
    extension::String=".arrow",
    subset_fraction=1.0
)

    # printstyled("aggrigating and reducing $extension geotiles\n"; color = :blue, bold = true)
    path2geotile = joinpath(geotile_dir, "$(geotiles.id)$(extension)")
    if !isfile(path2geotile)
        df0 = DataFrame()
        return df0
    else
        df = select!(DataFrame(Arrow.Table(path2geotile)), vars)
        isarray = [isa(d, AbstractArray) for d in df[1, :]]
        refarray = findfirst(isarray)

        # rearrange data into vectors 
        df0 = DataFrame()

        # random subsample requested?
        ind = eachindex(eachrow(df))
        if subset_fraction < 1
            ind = randsubseq(ind, 0.5)
        end

        for (i, var) in enumerate(vars)
            # in some cases geotile variables are saved as single values if they are constant 
            # for all row values while others are saved as vectors

            # check it it's a vector quantitiy
            if isarray[i]
                df0[!, var] = vcat(df[ind, var]...)
            else
                #if not expand to a vector 
                df0[!, var] = vcat([fill(df[k, var], size(df[k, refarray])) for k = eachindex(df[ind, var])]...)
            end
        end

        # remove empty rows
        delete!(df0, isempty.(df0[:, last(vars)]))
        return df0
    end
end


"""
    itslive_zone(lon, lat; always_xy=true) -> Tuple{Int, Bool}

Determine the UTM zone and hemisphere for ITS_LIVE projection based on coordinates.

# Arguments
- `lon`: Longitude in degrees
- `lat`: Latitude in degrees
- `always_xy=true`: If true, inputs are interpreted as (lon, lat); if false, as (lat, lon)

# Returns
- Tuple containing:
  - UTM zone number (0 for polar stereographic, -1 for invalid coordinates)
  - Boolean indicating northern hemisphere (true) or southern hemisphere (false)
"""
function itslive_zone(lon, lat; always_xy=true)
    if !always_xy
        lat, lon = (lon, lat)
    end

    # check for invalid conditions and return zone = -1
    if isnan(lon) || isnan(lat) 
        return (-1, false)
    end


    if lat > 55
        return (0, true)
    elseif lat < -56
        return (0, false)
    end

    # int versions
    ilon = floor(Int64, Geodesy.bound_thetad(lon))

    # zone
    zone = fld((ilon + 186), 6)

    isnorth = lat >= 0
    return (zone, isnorth)
end

"""
    itslive_paramfiles(lon, lat; gridsize=240, path2param=pathlocal.itslive_parameters, always_xy=true) -> NamedTuple

Get paths to ITS_LIVE parameter files for a given location.

# Arguments
- `lon`: Longitude in degrees
- `lat`: Latitude in degrees

# Keywords
- `gridsize=240`: Grid resolution in meters
- `path2param=pathlocal.itslive_parameters`: Base directory containing parameter files
- `always_xy=true`: If true, inputs are interpreted as (lon, lat); if false, as (lat, lon)

# Returns
- NamedTuple containing paths and binary flags for all ITS_LIVE parameter files
"""
function itslive_paramfiles(
    lon,
    lat;
    gridsize=240,
    path2param=pathlocal.itslive_parameters,
    always_xy=true)

    zone, isnorth = itslive_zone(lon, lat; always_xy=always_xy)

    if isnorth
        if zone == 0
            region = "NPS"
        else
            region = @sprintf("N%02.0f", zone)
        end
    else
        if zone == 0
            region = "SPS"
        else
            region = @sprintf("S%02.0f", zone)
        end
    end

    grid = @sprintf("%04.0fm", gridsize)

    paramfiles = (
        ROI=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "ROI")),
            binary=true
        ), StableSurface=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "StableSurface")),
            binary=true
        ), dhdx=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "dhdx")),
            binary=false
        ), dhdxs=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "dhdxs")),
            binary=false
        ), dhdy=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "dhdy")),
            binary=false
        ), dhdys=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "dhdys")),
            binary=false
        ), floatingice=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "floatingice")),
            binary=true
        ), glacierice=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "glacierice")),
            binary=true
        ), h=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "h")),
            binary=false
        ), inlandwater=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "inlandwater")),
            binary=true
        ), land=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "land")),
            binary=true
        ), landice=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "landice")),
            binary=true
        ), landice_2km_inbuff=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "landice_2km_inbuff")),
            binary=true
        ), ocean=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "ocean")),
            binary=true
        ), region=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "region")),
            binary=false
        ), sp=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "sp")),
            binary=true
        ), thickness=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "thickness")),
            binary=false
        ), vx=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "vx")),
            binary=false
        ), vx0=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "vx0")),
            binary=false
        ), vxSearchRange=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "vxSearchRange")),
            binary=false
        ), vy=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "vy")),
            binary=false
        ), vy0=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "vy0")),
            binary=false
        ), vySearchRange=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "vySearchRange")),
            binary=false
        ), xMaxChipSize=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "xMaxChipSize")),
            binary=false
        ), xMinChipSize=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "xMinChipSize")),
            binary=false
        ), yMaxChipSize=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "yMaxChipSize")),
            binary=false
        ), yMinChipSize=(
            path=joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid, "yMinChipSize")),
            binary=false
        )
    )

    return paramfiles
end

"""
    itslive_extract(lon, lat, vars; gridsize=240, path2param=pathlocal.itslive_parameters, always_xy=true)

Extract point values from ITS_LIVE parameter files for given coordinates.

# Arguments
- `lon::AbstractVector`: Vector of longitude coordinates
- `lat::AbstractVector`: Vector of latitude coordinates
- `vars::Vector{Symbol}`: Variables to extract (e.g., [:vx, :vy, :thickness])

# Keywords
- `gridsize=240`: Grid size in meters
- `path2param=pathlocal.itslive_parameters`: Path to parameter files
- `always_xy=true`: Whether coordinates are in x,y (longitude,latitude) order

# Returns
- `DataFrame`: Contains requested variables at the specified coordinates
"""
function itslive_extract(
    lon::AbstractVector,
    lat::AbstractVector,
    vars::Vector{Symbol};
    gridsize=240,
    path2param=pathlocal.itslive_parameters,
    always_xy=true)

    # determine unique zone
    zone_isnorth = itslive_zone.(lon, lat; always_xy=always_xy)
    uz = unique(zone_isnorth)

    df0 = DataFrame[]
    for z in uz
        ind = zone_isnorth .== Ref(z)
        if z[1] == -1
            df = DataFrame()
            paramfiles = itslive_paramfiles(0., -80.;
                gridsize=gridsize, path2param=path2param, always_xy=always_xy)

            # lat or lon is non-physical -or- nans
            for v in vars
                ga = GeoArrays.read(paramfiles[v].path, masked=false)

                if paramfiles[v].binary
                    df[!, v] = fill(false, sum(ind))
                else
                    df[!, v] = fill(convert(eltype(ga.A), nodatavalue), sum(ind))
                end
            end
            push!(df0, df)
            continue
        end

        paramfiles = itslive_paramfiles(mean(lon[ind]), mean(lat[ind]);
            gridsize=gridsize, path2param=path2param, always_xy=always_xy)

        # get regional parameter files
        for v in vars
            if !haskey(paramfiles, v)
                error(error("$(v) is not a valid variable, must be one of $(keys(paramfiles))"))
            end
        end

        # create projeciton transfomation function
        # NOTE: not worth checking for same projection... faster just to do it this way
        v = vars[1]
        ga = GeoArrays.read(paramfiles[v].path, masked=false)
        if always_xy
            xy = epsg2epsg(lon[ind], lat[ind], "EPSG:4326", ga.crs.val; parse_output=false)
        else
            xy = epsg2epsg(lat[ind], lon[ind], "EPSG:4326", ga.crs.val; parse_output=false)
        end

        # find xy extents
        x = getindex.(xy, 1)
        y = getindex.(xy, 2)
        minmax_x = extrema(x)
        minmax_y = extrema(y)

        extent = (min_x=minmax_x[1], min_y=minmax_y[1], max_x=minmax_x[2] + abs(ga.f.linear[1, 1]), max_y=minmax_y[2] + abs(ga.f.linear[2, 2]))

        if !bbox_overlap(extent2nt(bbox(ga)), extent)
            df = DataFrame()
            # extents do not overlap
            for v in vars
                ga = GeoArrays.read(paramfiles[v].path, masked=false)

                if paramfiles[v].binary
                    df[!, v] = fill(false, length(x))
                else
                    df[!, v] = fill(convert(eltype(ga.A), nodatavalue), length(x))
                end
            end
            push!(df0, df)
            continue
        else
            df = DataFrame()
            v = vars[1]
            ga0 = GeoArrays.crop(ga, extent)
            ind = CartesianIndex.(Tuple.(GeoArrays.indices.(Ref(ga0), xy)))
            
            if paramfiles[v].binary
                df[!, v] = ga0[ind] .== 1
            else
                df[!, v] = ga0[ind]
            end

            for v in vars[2:end]
                ga = GeoArrays.read(paramfiles[v].path, masked=false)
                ga0 = GeoArrays.crop(ga, extent)

                if paramfiles[v].binary
                    df[!, v] = ga0[ind] .== 1
                else
                    df[!, v] = ga0[ind]
                end
            end
            push!(df0, df)
        end
    end
    df0 = vcat(df0...)
    df = copy(df0)

    # this places data from multiple zones back into the original request order
    if length(uz) > 1
        s = 1
        for z in uz
            ind = zone_isnorth .== Ref(z)
            e = s + sum(ind) - 1

            df[ind, :] .= df0[s:e, :]

            s = e + 1
        end
    end

    return df
end

"""
    dem_info(dem::Symbol) -> (NamedTuple, String)

Get information about a Digital Elevation Model (DEM).

# Arguments
- `dem::Symbol`: DEM identifier (`:cop30_v1`, `:cop30_v2`, `:nasadem_v1`, 
                 `:arcticdem_v3_10m`, `:arcticdem_v4_10m`, or `:rema_v2_10m`)

# Returns
- `NamedTuple`: Contains DEM metadata (filename, date, nodatavalue, height datum, 
                EPSG code, resolution, and extent)
- `String`: Path to the geoid folder

# Throws
- `ErrorException`: If the DEM identifier is not recognized
"""
function dem_info(dem::Symbol)
    goids_folder = pathlocal.geoid_dir

    if dem == :cop30_v1
        info = (
            fn=pathlocal.cop30_v1,
            nominal_date=2013.0,
            nodatavalue=0,
            hight_datum="egm2008",
            epsg=4326,
            res=30,
            extent=(min_x=-180.0, min_y=-90.0, max_x=180.0, max_y=90.0)
        )
    elseif dem == :cop30_v2
        info = (
            fn=pathlocal.cop30_v2,
            nominal_date=2013.0,
            nodatavalue=0,
            hight_datum="egm2008",
            epsg=4326,
            res=30,
            extent=(min_x=-180.0, min_y=-90.0, max_x=180.0, max_y=90.0),
        )
    elseif dem == :nasadem_v1
        info = (
            fn=pathlocal.nasadem_v1,
            nominal_date=2000.1328,
            nodatavalue=-32767,
            hight_datum="wgs84",
            epsg=4326,
            res=30,
            extent=(min_x=-180.0, min_y=-56.0, max_x=180.0, max_y=60.0)
        )
    elseif dem == :arcticdem_v3_10m
        info = (
            fn=pathlocal.arcticdem_v3_10m,
            nominal_date=(2013.33 + 2023.92) / 2,
            nodatavalue=-9999.0,
            hight_datum="wgs84",
            epsg=3413,
            res=10,
            extent=(min_x=-180.0, min_y=50.0, max_x=180.0, max_y=90.0)
        )
    elseif dem == :arcticdem_v4_10m
        info = (
            fn=pathlocal.arcticdem_v4_10m,
            nominal_date=(2013.33 + 2023.92) / 2,
            nodatavalue=-9999.0,
            hight_datum="wgs84",
            epsg=3413,
            res=10,
            extent=(min_x=-180.0, min_y=50.0, max_x=180.0, max_y=90.0)
        )
    elseif dem == :rema_v2_10m
        info = (
            fn=pathlocal.rema_v2_10m,
            ominal_date=(2009.67 + 2023.0) / 2,
            nodatavalue=-9999.0,
            hight_datum="wgs84",
            epsg=3031,
            res=10,
            extent=(min_x=-180.0, min_y=-90.0, max_x=180.0, max_y=-61.0)
        )
    else
        error("dem not recognized, valid dem names are :cop30_v1, :nasadem_v1, :rema_v2_10m, and :arcticdem_v3_10m")
    end

    return info, goids_folder
end

"""
    geotile_extract_mask(
        geotile_id::String,
        geotile_dir::String;
        masks::Vector{Symbol}=[:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean],
        masks_highres::Vector{Symbol}=[:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km], # or nothing
        geotile_ext::String=".arrow",
        job_id="",
        force_remake=false
    ) -> Nothing

Extract mask data for a specific geotile and save as an Arrow file.

# Arguments
- `geotile_id::String`: Identifier for the geotile
- `geotile_dir::String`: Directory containing geotile files

# Keywords
- `masks::Vector{Symbol}`: Mask variables to extract (default: [:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean])
- `geotile_ext::String`: File extension of geotile files (default: ".arrow")
- `job_id::String`: Identifier for logging purposes (default: "")
- `force_remake::Bool`: Whether to regenerate existing files (default: false)

# Returns
Nothing, but creates a mask file for the geotile with extracted data
"""
function geotile_extract_mask(
    geotile_id::String,
    geotile_dir::String;
    masks::Vector{Symbol}=[:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean],
    masks_highres::Vector{Symbol}=[:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km], # or nothing
    geotile_ext::String=".arrow",
    job_id="",
    force_remake=false
)

    path2geotile = joinpath(geotile_dir, geotile_id * geotile_ext)
    outfile = replace(path2geotile, ".arrow" => ".masks")

    if isfile(outfile) && !force_remake
        printstyled("\n    ->$job_id $geotile_id masks already exists, skipping\n"; color=:green)
    elseif isfile(path2geotile)
        t1 = time()

        df0 = select!(DataFrame(Arrow.Table(path2geotile)), [:longitude, :latitude])
        #vlength = length.(df0.longitude)

        #lon = vcat(df0.longitude...)
        #lat = vcat(df0.latitude...)

        if isempty(df0.longitude)
            printstyled("\n    ->$job_id $geotile_id is all empty, skipping\n"; color=:light_red)
            return
        end

        mask0 = itslive_extract(df0.longitude, df0.latitude, unique(masks), path2param=pathlocal.itslive_parameters)

        if !isnothing(masks_highres)

            for highres_mask in masks_highres
                # handle to shapefiles 
                # [!!! this needs to be done inside of the thread or you will get a segfault !!!]
                if highres_mask == :land
                    shp = Symbol("$(:water)_shp")
                    fn_shp = pathlocal[shp]
                    feature = Shapefile.Handle(fn_shp)
                    invert = true

                    shp = Symbol("$(:landice)_shp")
                    fn_shp = pathlocal[shp]
                    excludefeature = Shapefile.Handle(fn_shp)
                else
                    shp = Symbol("$(highres_mask)_shp")
                    fn_shp = pathlocal[shp]
                    feature = Shapefile.Handle(fn_shp)
                    invert = false
                    excludefeature = nothing
                end
                mask0 = highres_mask!(mask0, df0, geotile_extent(geotile_id), feature, invert, excludefeature, highres_mask)
            end
        end

        tmp = tempname(dirname(outfile))
        Arrow.write(tmp, mask0::DataFrame)
        mv(tmp, outfile; force=true)

        total_time = round((time() - t1) / 60, digits=2)
        printstyled("\n    ->$job_id $geotile_id masks extracted: $(total_time) min \n"; color=:light_black)
    else
        printstyled("\n    ->$job_id $geotile_id does not exist\n"; color=:yellow)
        return
    end
end

"""
    geotile_extract_mask(geotiles::DataFrame, geotile_dir::String; masks, masks_highres, job_id, force_remake) -> Nothing

Extract mask data for multiple geotiles and save as Arrow files.

# Arguments
- `geotiles::DataFrame`: DataFrame containing geotile information with 'id' column
- `geotile_dir::String`: Directory containing geotile files

# Keywords
- `masks::Vector{Symbol}=[:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean]`: Mask variables to extract
- `masks_highres::Vector{Symbol}=[:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean]`: High-resolution mask variables to extract (or nothing)
- `job_id::String=""`: Identifier for logging purposes
- `force_remake::Bool=false`: Whether to regenerate existing files

# Returns
Nothing, but creates mask files for each geotile
"""

function geotile_extract_mask(
    geotiles::DataFrame,
    geotile_dir::String;
    masks::Vector{Symbol}=[:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean],
    masks_highres::Vector{Symbol}=[:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km], # or nothing
    job_id="",
    force_remake=false
)

    #asyncmap(eachrow(geotiles); ntasks = 4) do geotile
    @showprogress dt = 10 desc = "extracting masks for $job_id geotiles..." Threads.@threads for geotile in eachrow(geotiles)
        geotile_extract_mask(geotile.id, geotile_dir; masks, masks_highres,job_id, force_remake)
    end
end

"""
    pointextract(df0::DataFrame, ga::GeoArray; nodatavalue=0.0, var_name::Symbol=:var, interp=Constant()) -> DataFrame

Extract values from a GeoArray at geographic coordinates specified in a DataFrame.

# Arguments
- `df0::DataFrame`: DataFrame containing longitude and latitude columns
- `ga::GeoArray`: GeoArray containing the data to extract from

# Keywords
- `nodatavalue=0.0`: Value to use for points outside the array or at nodata locations
- `var_name::Symbol=:var`: Name of the column to create in the output DataFrame
- `interp=Constant()`: Interpolation method (Constant(), Linear(), Cubic(), etc.)

# Returns
- DataFrame containing extracted values in a column named by `var_name`
"""
function pointextract(
    df0::DataFrame,
    ga::GeoArray;
    nodatavalue=0.0,
    var_name::Symbol=:var,
    interp::F=Constant()
) where {F<:Interpolations.Flag}

    df = DataFrame()

    if isa(df0.longitude[1], AbstractVector)
        vlength = length.(df0.longitude)
        lon = reduce(vcat, df0.longitude)
        lat = reduce(vcat, df0.latitude)
    else
        lon = df0.longitude;
        lat = df0.latitude;
    end

    if isempty(lon)
        return df
    end

    var = pointextract(lon, lat, "EPSG:4326", ga; nodatavalue, interp)

    if all(var .== nodatavalue) || all(isnan.(var))

    elseif !isa(df0.longitude[1], AbstractVector)
        df[!, var_name] = vec(var)
    else
        df[!, var_name] = [Vector{eltype(var)}(undef, l) for l in vlength]

        start = 1
        for row = eachrow(df)
            stop = length(row[var_name]) + start - 1
            row = var[start:stop]
            start = stop + 1
        end
    end
    return df
end

"""
    geotile_pointextract(geotiles, geotile_dirs, ga; kwargs...)

Extract values from a GeoArray at points defined in geotiles and save results to files.

# Arguments
- `geotiles::DataFrame`: DataFrame containing geotile information
- `geotile_dirs`: Directory or vector of directories containing geotile files
- `ga::GeoArray`: GeoArray containing the data to extract from

# Keywords
- `var_name=:var`: Name of the column to create in the output files
- `nodatavalue=0.0`: Value to use for points outside the array or at nodata locations
- `job_ids=""`: Job identifier(s) for logging purposes
- `force_remake=false`: Whether to overwrite existing output files
- `interp=Constant()`: Interpolation method (Constant(), Linear(), Cubic(), etc.)

# Returns
Nothing, but creates files with extracted values in the specified directories
"""
function geotile_pointextract(
    geotiles::DataFrame,
    geotile_dirs,
    ga::GeoArray;
    var_name=:var,
    nodatavalue=0.0,
    job_ids="",
    force_remake=false,
    interp::F=Constant()
) where {F<:Interpolations.Flag}

    geotile_dirs isa AbstractVector || (geotile_dirs = [geotile_dirs])
    job_ids isa AbstractVector || (job_ids = [job_ids])

    #asyncmap(eachrow(geotiles); ntasks = 4) do geotile
    #@showprogress dt = 10 desc = "extracting $var_name for geotiles..." Threads.@threads for geotile in eachrow(geotiles)
    @showprogress dt = 1 desc = "extracting $var_name for geotiles..." for geotile in eachrow(geotiles)
    #for geotile in eachrow(geotiles)
        t1 = time()

        # check if ga is in EPSG:4326, if so then it's possible to check if lats and lons fall within ga 
        if Base.contains(ga.crs.val, "\"EPSG\",\"4326\"")

            if !bbox_overlap(extent2nt(bbox(ga)), extent2nt(geotile.extent))
                printstyled("\n    ->$(geotile.id) does not overlap $var_name, skipping\n"; color=:light_red)
                continue
            end
        end

        # geotile_pointextract() handles multiple geotile_dirs [multiple missions can be extracted at the same time... this saves greatly on I/O]
        path2geotiles = joinpath.(geotile_dirs, Ref(geotile.id * ".arrow"))
        path2outfiles = replace.(path2geotiles, Ref(".arrow" => ".$var_name"))

        outfiles_exist = isfile.(path2outfiles)
        geotiles_exist = isfile.(path2geotiles)

        extract = geotiles_exist .& (.!outfiles_exist .| force_remake)

        if any(extract)
            path2geotiles = path2geotiles[extract]
            path2outfile = path2outfiles[extract]
            job_id = job_ids[extract]

            df = DataFrame()
            stop = zeros(size(path2geotiles))
            for (i, path2geotile) in enumerate(path2geotiles)
                df = vcat(df, DataFrame(Arrow.Table(path2geotile))[!, [:longitude, :latitude]])
                stop[i] = nrow(df)
            end
            stop = Int.(stop)

            #try
            df = pointextract(df, ga; var_name, nodatavalue, interp)
            #catch
            #    return df 
            #end
            if isempty(df)
                printstyled("\n    ->$(geotile.id) no valid $var_name, skipping\n"; color=:light_red)
            else
                start = Int(1)
                for (i, outifle) = enumerate(path2outfile)
                    tmp = tempname(dirname(outifle))
                    Arrow.write(tmp, df[start:stop[i], :]::DataFrame)
                    mv(tmp, outifle; force=true)
                    start = stop[i] + 1
                end
            end
            total_time = round((time() - t1) / 60, digits=2)
            printstyled("\n    ->$job_id $(geotile.id) $var_name extracted: $(total_time) min \n"; color=:light_black)
            
        elseif all(.!geotiles_exist)

            job_id = job_ids[.!geotiles_exist]
            printstyled("\n    ->$job_id $(geotile.id) does not exist,  skipping\n"; color=:light_red)
        else
            job_id = job_ids[geotiles_exist .& outfiles_exist]
            printstyled("\n    ->$(geotile.id) $var_name already extracted, skipping\n"; color=:light_green)
        end
    end
end

"""
    intersectindices(a, b; bool=false)

Find indices of intersecting elements between collections `a` and `b`.

# Arguments
- `a`: First collection
- `b`: Second collection
- `bool`: If true, returns boolean masks instead of indices (default: false)

# Returns
- If `bool=false`: Tuple of indices `(ia, ib)` where `a[ia]` equals `b[ib]`
- If `bool=true`: Tuple of boolean masks `(ia0, ib0)` marking intersecting elements
"""
function intersectindices(a, b; bool=false)
    ia = findall(in(b), a)
    ib = findall(in(view(a, ia)), b)
    # same order as ia such that a[ia] == b[ib]
    ib = ib[indexin(view(a, ia), view(b, ib))]

    if bool != true
        return ia, ib
    else
        ia0 = falses(size(a))
        ia0[ia] .= true
        ib0 = falses(size(b))
        ib0[ib] .= true

        return ia0, ib0
    end
end

"""
    rightmerge(df_left::DataFrame, df_right::DataFrame, id_unique::Symbol)

Merge two DataFrames based on a unique identifier column, prioritizing values from df_right.

# Arguments
- `df_left`: First DataFrame
- `df_right`: Second DataFrame
- `id_unique`: Symbol representing the column name containing unique identifiers

# Returns
- DataFrame with merged data where df_right values take precedence for matching rows
"""
function rightmerge(df_left::DataFrame, df_right::DataFrame, id_unique::Symbol)
    ileft, iright = intersectindices(df_left[:, id_unique], df_right[:, id_unique], bool=true)
    if any(ileft)
        df_left[ileft, :] = df_right[iright, :] # replace duplicates
    end
    if any(.!iright)
        df_left = vcat(df_left, df_right[.!iright, :]) # add new rows
    end
    return df_left
end

"""
    leftmerge(df_left::DataFrame, df_right::DataFrame, id_unique::Symbol)

Merge two DataFrames based on a unique identifier column, keeping all rows from df_left
and adding only new rows from df_right.

# Arguments
- `df_left`: First DataFrame (takes precedence for duplicate IDs)
- `df_right`: Second DataFrame (only non-duplicate rows are added)
- `id_unique`: Symbol representing the column name containing unique identifiers

# Returns
- DataFrame with merged data where all df_left rows are preserved
"""
function leftmerge(df_left::DataFrame, df_right::DataFrame, id_unique::Symbol)
    _, iright = intersectindices(df_left[:, id_unique], df_right[:, id_unique], bool=true)
    if any(.!iright)
        df_left = vcat(df_left, df_right[.!iright, :]) # add new rows
    end
    return df_left
end

"""
    geotile_build(geotile_granules, geotile_dir; warnings=true, fmt=:arrow, replace_corrupt_h5=true)

Build geotiles from satellite data granules.

# Arguments
- `geotile_granules`: DataFrame containing geotiles and their associated granules
- `geotile_dir`: Directory where geotiles will be saved

# Keywords
- `warnings=true`: Whether to display warning messages
- `fmt=:arrow`: Output file format (`:arrow` or other supported format)
- `replace_corrupt_h5=true`: Whether to attempt recovery of corrupt HDF5 files

# Returns
Nothing, but creates geotile files in the specified directory
"""
function geotile_build(geotile_granules, geotile_dir; warnings=true, fmt=:arrow, replace_corrupt_h5=true)
    printstyled("building geotiles\n"; color=:blue, bold=true)

    # remove empty granules
    geotile_granules = geotile_granules[.!isempty.(geotile_granules.granules), :]

    if !warnings
        Logging.disable_logging(Logging.Warn) # or e.g. Logging.Info
    end

    #asyncmap(eachrow(geotile_granules); ntasks = 10) do row
    @showprogress dt = 1 desc = "Building geotiles for $(geotile_dir)..." for row in eachrow(geotile_granules)

        # tiles
        outfile = joinpath(geotile_dir, row.id * ".$fmt")

        if isempty(row.granules)
            # do nothing
            printstyled("\n    -> $(row.id): has no granules\n"; color=:light_red)
            continue
        else
            if isfile(outfile)
                # printstyled("$(row.id) file exisits, adding to"; color = :light_green)
                # if the GeoTile already exists, check for overlap

                if fmt == :arrow
                    df0 = DataFrame(Arrow.Table(outfile))
                else
                    df0 = FileIO.load(outfile, "df")
                end

                if isempty(df0)
                    ind1 = falses(size(row.granules))
                    ind0 = nothing
                else
                    id0 = unique(df0.id)
                    id1 = [g.id for g in row.granules]
                    ind1, ind0 = intersectindices(id1, id0; bool=true)
                end
            else
                # printstyled("$(row.id) new file "; color = :green)
                ind1 = falses(size(row.granules))
                ind0 = nothing
            end

            if !all(ind1)
                # subset granules to only those that are not already in table df0
                granules = row.granules[.!(ind1)]

                # start read timer
                t1 = time()

                df = reduce(vcat, DataFrame.(SpaceLiDAR.add_id.(getpoints.(granules; extent=row.extent, replace_corrupt_h5=replace_corrupt_h5))))

                # in very rare occations the data can fall outside of the suplied extent as 
                # the subseting is only done on the hdf5 for start and stop read locations
                ind = within.(Ref(row.extent), df.longitude, df.latitude) .| isnan.(df.longitude)
                deleteat!(df, .!ind)

                # add a placeholder for geotiles that were searched but did not return anything
                # this is done to record which granules have been searched... otherwise searched 
                # tiles without data will be seached again

                id0X = unique(df.id)
                id1X = [g.id for g in granules]
                ind1X, ind0X = intersectindices(id1X, id0X; bool=true)

                if any(.!ind1X)
                    er = emptyrow(df)
                    for idX = id1X[.!ind1X]
                        er[end] = idX
                        df = push!(df, er)
                    end
                end

                if !isnothing(ind0) && any(ind0)
                    df = vcat(df0::DataFrame, df::DataFrame)
                end
                read_time = round((time() - t1) / 60, digits=1)

                # start write timer
                t1 = time()
                tmp = tempname(dirname(outfile)) * ".$fmt"
                if fmt == :arrow
                    Arrow.write(tmp, df::DataFrame)
                else
                    save(tmp, Dict("df" => df::DataFrame))
                end

                mv(tmp, outfile; force=true)
                write_time = round((time() - t1) / 60, digits=1)
                printstyled("\n    -> $(row[:id]): generation complete [read: $read_time min, write: $write_time min]\n"; color=:light_black)
            else
                printstyled("\n    -> $(row[:id]): no new granules to add to exisitng GeoTile\n"; color=:light_green)
            end
        end
    end
end


"""
    emptyrow(df)

Create an empty row for a DataFrame with appropriate default values for each column type.

# Arguments
- `df`: DataFrame to create an empty row for

# Returns
- Vector with default values for each column (empty strings for String columns, 
  NaN for floating point columns, and 0 for other numeric types)
"""
function emptyrow(df)
    tps = eltype.(eachcol(df))
    f = []
    for t in tps
        if t <: String
            f = push!(f, "0")
        elseif t <: AbstractFloat
            f = push!(f, t(NaN))
        else
            f = push!(f, t(0))
        end
    end
    return f
end

"""
    mission2spacelidar(mission::Symbol) -> Symbol

Convert mission symbol to standardized SpaceLiDAR capitalization format.

# Arguments
- `mission::Symbol`: Input mission symbol (e.g., `:icesat`, `:ICESat2`, `:GEDI`)

# Returns
- `Symbol`: Standardized mission symbol (`:ICESat`, `:ICESat2`, or `:GEDI`)

# Throws
- `ErrorException`: If the mission is not recognized
"""
function mission2spacelidar(mission::Symbol)
    if !(mission == :ICESat || mission == :ICESat2 || mission == :GEDI)
        if mission == :icesat
            mission = :ICESat
        elseif mission == :icesat2
            mission = :ICESat2
        elseif mission == :gedi
            mission = :GEDI
        else
            error("unrecognized mission")
        end
    end

    return mission
end

"""
    write_urls!(fn::String, urls::Vector{String}) -> String

Write a list of URLs to a file, one URL per line.

# Arguments
- `fn::String`: Path to the output file
- `urls::Vector{String}`: List of URLs to write

# Returns
- Absolute path to the created file
"""
function write_urls!(fn::String, urls::Vector{String})
    open(fn, "w") do f
        for url in urls
            println(f, url)
        end
    end
    abspath(fn)
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
    ras0 = Rasters.read(Rasters.crop(ras, to=geotile.extent))

    # using rasterize is 15 times faster than reading its_live masks !!!!
    mask0 = Rasters.rasterize!(count, zeros(ras0.dims), feature; threaded=false, shape=:polygon, progress=false, verbose=false) .> 0

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
    d = meters2lonlat_distance.(Ref(1), lat)
    a = abs.((1 ./ getindex.(d, 2) * (lat[2] .- lat[1])) .* (1 / d[1][1] * (lon[2] - lon[1])))
    area_m2 = repeat(a', outer=[length(lon), 1])


    feature_area_km2 = mask0 .* area_m2 / (1000 * 1000)

    foo = geotile[var_name]

    df = DataFrame(ras=ras0[mask0], feature_area_km2=feature_area_km2[mask0])

    dfbs = BinStatistics.binstats(df, :ras, bin_edges, :feature_area_km2, col_function=sum, missing_bins=true)

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
    gt_file = joinpath(binned_folder, "geotile_$(surface_mask)_hyps.arrow")
    geotiles = DataFrame(Arrow.Table(gt_file))
    geotiles.extent = Extent.(getindex.(geotiles.extent, 1))
    geotiles = copy(geotiles)


    # As of now I can not figure out how to restore geometry using GI.Polygon
    df2 = project_geotiles(; geotile_width)
    if !any(Base.contains.(names(geotiles), "geometry"))
        geotiles = innerjoin(geotiles, df2[!, [:id, :geometry]], on=[:id])
    else
        geotiles = innerjoin(geotiles[!, Not(:geometry)], df2[!, [:id, :geometry]], on=[:id])
    end

    # println(gt_file)
    return geotiles
end