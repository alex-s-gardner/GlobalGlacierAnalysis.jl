# utilities for working building geogrid point database
include("local_paths.jl")
global pathlocal = setpaths()
const world = Extent(X = (-180, 180), Y = (-90, 90))

"""
    setpaths(geotile_width, mission, product, version)

populate named tuple of paths
"""
function setpaths(geotile_width, mission, product, version)
    geotile_dir = @sprintf "%.0fdeg" geotile_width;
    data_dir = joinpath(pathlocal.data_dir,lowercase(string(mission)), 
        string(product), lpad(string(version), 3, '0'))
    raw_data_dir = joinpath(data_dir, "raw");
    geotile_dir = joinpath(data_dir, "geotile", geotile_dir)
    granules_remote = joinpath(data_dir, "geotile", geotile_dir, "granules.remote")
    granules_local = joinpath(data_dir, "geotile", geotile_dir, "granules.local")
    
    altim_paths = (
        raw_data = raw_data_dir,
        geotile = geotile_dir,
        granules_remote = granules_remote,
        granules_local = granules_local
        )
    return altim_paths
end

"""
    geotile_extent(lat, lon, width)

Return extents of geotile
"""
function geotile_extent(lat, lon, width)
    (min_x = (lon - width/2), min_y = (lat - width/2), 
    max_x = (lon + width/2), max_y = (lat + width/2));
end

"""
    within(extent, x, y)

Determine if a point falls within extents
"""
function within(extent::NamedTuple{(:min_x, :min_y, :max_x, :max_y)}, x, y)
    in = (x >= extent.min_x) .&& (x <= extent.max_x) .&& (y >= extent.min_y) .&& (y <= extent.max_y)
    return in
end

"""
    within(extent, x, y)

Determine if a point falls within Extent
"""
function within(extent::Extent, x, y)
    in = (x >= extent.X[1]) .&& (x <= extent.X[2]) .&& (y >= extent.Y[1]) .&& (y <= extent.Y[2])
    return in
end

"""
    crop!(df::DataFrame, geocell::NamedTuple)

Crop dataframe to only include data that falls within a given geocell
"""
function crop!(df::DataFrame, geocell::NamedTuple)
    for row = eachrow(df)
        ind = ((row.latitude .>= geocell[2]) .& (row.latitude .< geocell[4]) 
            .& (row.longitude .>= geocell[1]) .& (row.longitude .< geocell[3]))
        for j = eachindex(row)
            row[j] = row[j][ind]
            return nothing
        end
        return nothing
    end
end

function crop!(pt::Vector{NamedTuple}, geocell::NamedTuple)
    for row = eachrow(pt)
        ind = ((row.latitude .>= geocell[2]) .& (row.latitude .< geocell[4]) 
        .& (row.longitude .>= geocell[1]) .& (row.longitude .< geocell[3]))
        for j = eachindex(row)
            row[j] = row[j][ind]
            return nothing
        end
        return nothing
    end
end

# seach directory using two keys
searchdir(path,key1,key2) = filter(x->(occursin(key1,x) .& occursin(key2,x)), readdir(path))
searchdir(path,key1) = filter(x->occursin(key1,x), readdir(path))

"""
    points_plus(granule::ICESat2_Granule{}; bbox = (min_x = -Inf, min_y = -Inf, max_x = Inf, max_y = Inf))

returns the ICESat2 granual *WITH* granual infomation for each track
"""
function points_plus(
    granule;
    extent::Extent=world
)

    try
        p = points(granule, bbox=extent)
        ptsplus = merge.(SpaceLiDAR.parent(p), Ref((; granule_info=granule)))
        return ptsplus
    catch ex
        println("-------- error thrown when trying to read: --------")
        println("$(granule.url)")
        println("pease delete file and re-run `geotile_download_granules`")
        println("------------------------------------------------------")
        throw(error("issue reading $(granule.url)"))
    end
end

"""
    points_plus(granule::ICESat_Granule{}; bbox = (min_x = -Inf, min_y = -Inf, max_x = Inf, max_y = Inf))

returns the ICESat granual *WITH* granual infomation for each track
"""
function points_plus(
    granule::ICESat_Granule{};
    extent::Extent=world
)

    try
        p = points(granule, bbox=extent)
        ptsplus = merge(SpaceLiDAR.parent(p), (; granule_info=granule))
        return ptsplus
    catch ex
        println("-------- error thrown when trying to read: --------")
        println("$(granule.url)")
        println("pease delete file and re-run `geotile_download_granules`")
        println("------------------------------------------------------")
        throw(error("issue reading $(granule.url)"))
    end
end

function rm_corrupt_h5(folder; minsize = Inf) 
    # minsize = 1E7 # file size is in bytes
    files = filter!(readdir(folder)) do fname
        fname[end] == '5' .&& filesize(joinpath(folder,fname)) < minsize
    end

    # this checks each file before deleting
    for fname in files
        if fname[end] == '5'
            try 
                h5open(joinpath(folder,fname), "r")
            catch e
                if !any(names(e) .== "msg")
                    error("check that HDF5 library has been imported (e.g. using HDF5)")
                else
                    println(e.msg)
                    println("!!!! deleting file !!!!")
                    rm(joinpath(folder,fname))
                end
            end
        end
    end
end

"""
    geotile_define(geotile_width)

Returns a DataFrame with geotile ids and extents
"""
function geotile_define(geotile_width::Number)

    if mod(180, geotile_width) != 0
        error("a geotile width of $geotile_width does not divide evenly into 180")
    end
    
    # geotile package uses fixed extents for consistancy, once defined a supset of tiles can be selected
    dt = geotile_width/2;
    lat_center = (world.Y[1] + dt):geotile_width:(world.Y[2] - dt)
    lon_center = (world.X[1] + dt):geotile_width:(world.X[2] - dt)

    extent = vec([geotile_extent(lat, lon, geotile_width) for lon in(lon_center), lat in(lat_center)]);
    id = geotile_id.(extent)
    
    DataFrame(id = id, extent = extent)
end

"""
    mission2spacelidar(mission)

Returns mission symbol with spacelidar capitilization
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
geotile_search_granules(geotiles, mission, product, version, outgranulefile; rebuild_dataframe = false )

Find all granuels that intersect geotile and same as an Arrow table
"""
function geotile_search_granules(
    geotiles, 
    mission, 
    product, 
    version, 
    outgranulefile; 
    rebuild_dataframe = false 
)

    # ensure SpaceLiDAR capitilization ot mission 
    mission = mission2spacelidar(mission)

    pathfile = splitdir(outgranulefile)
    if !isdir(pathfile[1])
        error("$(pathfile[1]) does not exist")
    end

    printstyled("identifying grannules that intersect each geotile\n"; color = :blue, bold = true)
    
    # get intersecting granules [run threads]
    n = size(geotiles,1)
    if mission == :ICESat2
        granules = Vector{Vector{ICESat2_Granule{product}}}(undef,n)
    elseif mission == :ICESat
        granules = Vector{Vector{ICESat_Granule{product}}}(undef,n)
    elseif mission == :GEDI
        granules = Vector{Vector{GEDI_Granule{product}}}(undef,n)
    else
        error("mission and/or product not recognized")
    end

    Threads.@threads for i in 1:n
    #for i in 1:n
        printstyled("    -> finding granules in geotile $i of $n\n"; color = :light_black)
        granules[i] = search(mission, product; bbox = convert(Extent, geotiles[i,:extent]), version = version)
    end

    # save remote granules before download in case there is a crash
    geotile_granules = hcat(geotiles, DataFrame(granules = granules))

    # check if file already exists 
    if isfile(outgranulefile) && !rebuild_dataframe
        #if file exists read it in and only overwrite updated granules
        geotile_granules = leftmerge(geotile_granules, granules_load(outgranulefile, mission), :id)
    end

    tmp = tempname(dirname(outgranulefile))
    Arrow.write(tmp,  geotile_granules::DataFrame)
    mv(tmp, outgranulefile; force = true)
    return outgranulefile
end

function geotile_download_granules!(
    geotile_granules,
    mission,
    savedir, 
    outgranulefile;
    threads = true,
    rebuild_dataframe = false,
    aria2c = false, 
    downloadstreams = 16
    )

    printstyled("downloading granules for each geotile\n"; color = :blue, bold = true)

    # ensure SpaceLiDAR capitilization ot mission 
    mission = mission2spacelidar(mission)

    # remove empty granules
    geotile_granules = geotile_granules[.!isempty.(geotile_granules.granules),:]

    filesOnDisk = readdir(savedir)

    n = size(geotile_granules,1)
    for (i, row) in enumerate(eachrow(geotile_granules))
        t1 = time();
        printstyled("    -> downloading granules $i of $n ... "; color = :light_black)
        
         # download seems to get killed when it makes too many requests... try just requesting files that actually need downloading. 
         files2download = [g.id for g in row.granules]
         ia, _ = intersectindices(files2download, filesOnDisk; bool = true)
         ia = .!ia;

        # download seems to get killed when it makes too many requests... try just requesting files that actually need downloading. 
        if any(ia)
            granules = row.granules[ia]
        else
            print("no new files to download, skipping\n")
            # update granule urls to local paths
            for g in row.granules
                g.url = joinpath(savedir,g.id)
            end
            continue
        end

        # using async is ~4x faster than without
        if aria2c
            urls = [g.url for g in granules]
            #urls = vcat(urls...)

            fn  = tempname()
            url_list = write_urls!(fn, urls)
            
            cmd = `aria2c --max-tries=10 --retry-wait=1 -x $downloadstreams -k 1M -j 1 -i --max-connection-per-server 15 -c -d $savedir -i $url_list`

            println(cmd)
            run(cmd)

            # update granule urls to local paths
            for g in row.granules
                g.url = joinpath(savedir,g.id)
            end
        else
            if threads
                asyncmap(granules; ntasks = 10) do g
                    flag = 0
                    while flag == 0
                        try
                            download!(g, savedir)
                            flag = 1
                        catch e
                            println(e)
                            wait_time = 10;
                            printstyled("download hand an issue... will try again in $wait_time s\n"; color = :yellow)
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
                            wait_time = 10;
                            printstyled("download hand an issue... will try again in $wait_time s\n"; color = :yellow)
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
    Arrow.write(tmp, geotile_granules::DataFrame);
    mv(tmp, outgranulefile; force = true)
    return outgranulefile
end

function write_urls!(fn::String, urls::Vector{String})
    open(fn, "w") do f
        for url in urls
            println(f, url)
        end
    end
    abspath(fn)
end

function granules_load(
    filename, 
    mission;
    geotiles::Union{Nothing,DataFrame} = nothing,
    rebuild_dataframe = false
    )
    
    # ensure SpaceLiDAR capitilization ot mission 
    mission = mission2spacelidar(mission)

    # load all granules
    granules = DataFrame(Arrow.Table(filename))
    
    # subset if requested
    if !isnothing(geotiles) && !rebuild_dataframe
        ind = [findfirst(granules[!,:id] .== id) for id in geotiles[!,:id]];
        valid = .!isnothing.(ind);
        if any(valid)
            granules = granules[ind[valid],:]
        else
            granules = nothing
            return granules
        end
    end

    # type is stripped on df save, add back
    polygon = []; # this was added as grannules were extracted using version 2.2 of SpaceLiDAR that didn't have polygons yet
    if mission == :ICESat2
        granules[!,:granules] =
            [[ICESat2_Granule{Symbol(g.info.type)}(g.id, g.url, g.info, polygon) for g in grans]  for grans in granules[!,:granules]]
    elseif mission == :ICESat
        granules[!,:granules] =
            [[ICESat_Granule{Symbol(g.info.type)}(g.id, g.url, g.info, polygon) for g in grans]  for grans in granules[!,:granules]]
    elseif mission == :GEDI
        granules[!,:granules] = 
            [[GEDI_Granule{Symbol(g.info.type)}(g.id, g.url, g.info, polygon) for g in grans]  for grans in granules[!,:granules]]
    else
        error("need to add mission")
    end
    return granules
end

"""
    intersectindices(a,b; bool = false)

Returns the index mapping between intersecting elements in `a` and `b`
"""
function intersectindices(a, b; bool = false)
    ia = findall(in(b), a)
    ib = findall(in(view(a,ia)), b)
    # same order as ia such that a[ia] == b[ib]
    ib = ib[indexin(view(a,ia), view(b,ib))]

    if bool != true
        return ia, ib
    else
        ia0 = falses(size(a))
        ia0[ia] .= true;
        ib0 = falses(size(b))
        ib0[ib] .= true;

        return ia0, ib0
    end
end

"""
    rightmerge(df_left::DataFrame, df_right::DataFrame, id_unique::Symbol)

Merge DataFrames with identical columns based on id_unique
- retains df_left rows if not contianed in df_right
- retains df_right rows if not contianed in df_left
- df_right rows are retained over df_left if rows have the same id_unique
"""
function rightmerge(df_left::DataFrame, df_right::DataFrame, id_unique::Symbol)
    ileft, iright = intersectindices(df_left[:,id_unique], df_right[:,id_unique], bool = true)
    if any(ileft)
        df_left[ileft,:] = df_right[iright,:] # replace duplicates
    end
    if any(.!iright)
        df_left = vcat(df_left, df_right[.!iright,:]) # add new rows
    end
    return df_left
end

"""
    leftmerge(df_left::DataFrame, df_right::DataFrame, id_unique::Symbol)

Merge DataFrames with identical columns based on id_unique
- retains df_left rows if not contianed in df_right
- retains df_right rows if not contianed in df_left
- df_left rows are retained over df_rigth if rows have the same id_unique
"""
function leftmerge(df_left::DataFrame, df_right::DataFrame, id_unique::Symbol)
    _, iright = intersectindices(df_left[:,id_unique], df_right[:,id_unique], bool = true)
    if any(.!iright)
        df_left = vcat(df_left, df_right[.!iright,:]) # add new rows
    end
    return df_left
end

function geotile_build(geotile_granules, geotile_dir, mission; warnings = true)

    # ensure SpaceLiDAR capitilization ot mission 
    mission = mission2spacelidar(mission)

    printstyled("building geotiles\n"; color = :blue, bold = true)

    # remove empty granules
    geotile_granules = geotile_granules[.!isempty.(geotile_granules.granules),:]

    if !warnings
        Logging.disable_logging(Logging.Warn) # or e.g. Logging.Info
    end

    #asyncmap(eachrow(geotile_granules); ntasks = 10) do row
    for row in eachrow(geotile_granules)
        # skip tiles
        outfile = joinpath(geotile_dir, row[:id]*".arrow")

        if isempty(row.granules)
            # do nothing
            printstyled("    -> $(row[:id]): has no granules\n"; color = :light_red)
            continue
        else
            if isfile(outfile)
                # printstyled("$(row.id) file exisits, adding to"; color = :light_green)
                # if the GeoTile already exists, check for overlap
                df0 = DataFrame(Arrow.Table(outfile))
                if isempty(df0)
                    ind1 = falses(size(row.granules))
                    ind0 = nothing
                else
                    id0 = [g.id for g in df0.granule_info]
                    id1 = [g.id for g in row.granules]
                    ind1, ind0 = intersectindices(id1,id0; bool = true)
                end 
            else
                # printstyled("$(row.id) new file "; color = :green)
                ind1 = falses(size(row.granules))
                ind0 = nothing
            end

            if !all(ind1)
                # start read timer
                t1 = time(); 
                if mission == :ICESat 
                    df = DataFrame(points_plus.(row.granules, extent = convert(Extent, row.extent)));    
                elseif mission == :ICESat2 || mission == :GEDI
                    df = reduce(vcat, (DataFrame.(points_plus.(row.granules, extent = convert(Extent, row.extent)))));   
                    # fill arrays [:track, :detector_id] track to a single string
                    df[!,:track] = [!isempty(row.track) ? row.track[1] : row.track for row in eachrow(df)];
                    # replace binary vector with single binary value
                    df[!,:strong_beam] = [!isempty(row.strong_beam) ? row.strong_beam[1] : row.strong_beam for row in eachrow(df)];

                    if mission == :ICESat2
                        df[!,:detector_id] = [!isempty(row.detector_id) ? row.detector_id[1] : row.detector_id for row in eachrow(df)];
                    end

                else
                    error("unrecognized mission")
                end

                if !isnothing(ind0) && any(ind0)
                    df = vcat(df::DataFrame, df0::DataFrame)
                end
                read_time = round((time()-t1)/60, digits=1)

                # start write timer
                t1 = time();
                tmp = tempname(dirname(outfile))
                Arrow.write(tmp, df::DataFrame);
                mv(tmp, outfile; force = true)
                write_time = round((time()-t1)/60, digits=1)
                printstyled("    -> $(row[:id]): generation complete [read: $read_time min, write: $write_time min]\n"; color = :light_black)
            else
                printstyled("    -> $(row[:id]): no new granules to add to exisitng GeoTile\n"; color = :light_green)
            end
        end
    end
end

"""
    geotile_id(extent)

Returns the geotile id given a geotile `extent`
"""
function geotile_id(extent)
    id = @sprintf("lat[%+03.0f%+03.0f]lon[%+04.0f%+04.0f]", extent.min_y, extent.max_y, extent.min_x, extent.max_x)
    return id
end

"""
    geotile_utm!(df::DataFrame{}; height = nothing)

add x and y coodinates for local utm or polar stereo zone to an geotile DataFrame
"""
function geotile_utm!(df::DataFrame; height = nothing)
    # create a regularly spaced grid using a local low distortion projection
    epsg = utm_epsg(mean(df.longitude), mean(df.latitude))
    df[!, :X], df[!, :Y] = epsg2epsg(df.longitude, df.latitude, "EPSG:4326", "$epsg", parse_output=true)

    return df, epsg
end


"""
    pointextract(x::Vector{<:Vector{<:Number}}, latitude::Vector{<:Vector{<:Number}}, 
        ga::GeoArray; nodatavalue = 0.0, replace_nodatavalue_withnan = false,
        filter_kernel::Union{OffsetMatrix{Float64, Matrix{Float64}}, Nothing} = nothing, 
        derivative_kernel::Union{OffsetArray, Nothing, Vector{<:OffsetArray}} = nothing,
        interp = Constant())
"""
function pointextract(
    x::Vector,
    y::Vector,
    point_epsg,
    ga::GeoArray;
    nodatavalue = 0.0,
    replace_nodatavalue_withnans = false,
    filter_kernel::Union{OffsetArray, Nothing,Vector{<:OffsetArray}}=nothing,
    derivative_kernel::Union{OffsetArray, Nothing, Vector{<:OffsetArray}} = nothing,
    interp::F = Constant()
    ) where {F<:Interpolations.Flag}

    # create projeciton transfomation function
    # NOTE: not worth checking for same projection... faster just to do it this way
    
    _, epsg_num = split(point_epsg, ":")
    if !contains(ga.crs.val[end-7:end], "\"EPSG\",\"$epsg_num\"")
        xy = epsg2epsg(x, y, point_epsg, ga.crs.val; parse_output = false)
        x = getindex.(xy, 1);
        y = getindex.(xy, 2);
    end

    # find x y extents
    minmax_x = extrema(x)
    minmax_y = extrema(y)
    dx = ga.f.linear[1,1];
    dy = ga.f.linear[2,2];

    # deteremine how much to pad crop by
    if !isnothing(derivative_kernel) .| !isnothing(filter_kernel)
        if isnothing(derivative_kernel)
            px = ceil(size(filter_kernel, 2) / 2)
            py = ceil(size(filter_kernel, 1) / 2)
        elseif isnothing(filter_kernel)
            px = ceil(size(derivative_kernel, 2) / 2)
            py = ceil(size(derivative_kernel, 1) / 2)
        else
            px = ceil(size(filter_kernel, 2) + size(derivative_kernel,2)/2)
            py = ceil(size(filter_kernel, 1) + size(derivative_kernel, 1)/2)
        end
    else
        px = 1;
        py = 1;
    end

    if interp == Linear()
        px += 1 
        py += 1 
    elseif interp == Cubic()
        px += 2 
        py += 2 
    elseif interp == Quadratic()
        px += 3 
        py += 3
    end 

    extent = (
        min_x = minmax_x[1] - abs(dx)*px, 
        min_y = minmax_y[1] - abs(dy)*py, 
        max_x = minmax_x[2] + abs(dx)*px, 
        max_y = minmax_y[2] + abs(dy)*py
        )

    if !bbox_overlap(bbox(ga), extent)
        # extents do not overlap
        if replace_nodatavalue_withnans
            val = fill(convert(eltype(ga.A), NaN), length(x))
        else
            val = fill(convert(eltype(ga.A), nodatavalue), length(x))
        end
        return val
    else
        # crop and read into  memory 

        ## ~65% of all time is spent here ##
        ga0 = crop(ga, extent)
        ###################################

        if replace_nodatavalue_withnans
            if (eltype(ga0.A) <:Integer)
                f = ga0.f
                ga0 = GeoArray(Float64.(ga0.A))
                ga0.crs = ga.crs
                ga0.f = f
            end

            for i = 1:size(ga0.A,3)
                @view(ga0.A[:,:,i])[ga0.A[:,:,i] .== nodatavalue] .= NaN
            end
            nodatavalue = NaN;            
        end

        # apply filter (typically a smoothing filter)
        if !isnothing(filter_kernel)
            for i = 1:size(ga0.A,3)
                @time imfilter!(@view(ga0.A[:,:,i]), ga0.A[:,:,i], filter_kernel);
            end
        end

        x0, y0 = GeoArrays.ranges(ga0)
        gridsize = xy_gridsize(ga0)
        y2x_scale = round(abs((gridsize.y/gridsize.x)) / (dy/dx), digits = 5);
        val = Matrix{eltype(ga0.A)}(undef, (length(x), size(ga0.A,3)))

        ## ~30% of all time is spent here ##
        for i in eachindex(ga0.A[1,1,:])
            itp = extrapolate(scale(interpolate(ga0.A[:,:,i], BSpline(interp)), x0, y0 .* y2x_scale), nodatavalue)
            val[:,i] = itp.(x, y .* y2x_scale)
        end

        ###################################
 
        # check if a derivative has been requested
        if !isnothing(derivative_kernel)
            # loop for multiple derivative
            if derivative_kernel isa Vector
                n = length(derivative_kernel)
            else
                n = 1
                derivative_kernel = [derivative_kernel]
            end

            ga1 = [copy(ga0) for i = 1:n]
            for j = 1:n
                for i = 1:size(ga1[j].A,3)
                        imfilter!(@view(ga1[j].A[:,:,i]), ga0.A[:,:,i], derivative_kernel[j]);
                end
            end

            deriv = [Matrix{eltype(ga0.A)}(undef,(length(x),size(ga0.A,3))) .+ nodatavalue for i = 1:n]
            for j in 1:n
                for i = 1:size(ga0.A,3)
                    itp = extrapolate(scale(interpolate(ga1[j].A[:,:,i], BSpline(interp)), x0, y0 .* y2x_scale), nodatavalue)
                    deriv[j][:,i] = itp.(x, y .* y2x_scale)
                end
            end
            val = vcat([val], deriv)
        end
        return val
    end
end

"""
    epsg2epsg_nodata(x, y, height, from_epsg, to_epsg, nodatavalue)

A wrapper for `epsg2epsg` to prevent modification of `height` no data values
"""
function epsg2epsg_nodata(x, y, height, from_epsg, to_epsg, nodatavalue)
    valid = (height .!= nodatavalue) .& (.!isnan.(height))
    if any(valid)
        #(x[valid], y[valid], height[valid]) = epsg2epsg(x[valid], y[valid], height[valid], from_epsg, to_epsg; parse_output = false)
        x[valid], y[valid], height[valid] = epsg2epsg(x[valid], y[valid], height[valid], from_epsg, to_epsg; parse_output = true)
    end
    return x, y, height
end


"""
    epsg2epsg(x, y, height, from_epsg, to_epsg)

Returns `x`, `y`,`height` in `to_epsg` projection
"""
function epsg2epsg(
    x::Union{Vector{<:Number}, Number}, 
    y::Union{Vector{<:Number}, Number}, 
    height::Union{Vector{<:Number}, Number},
    from_epsg::String,
    to_epsg::String;
    parse_output = true,
    threaded = true
    )
    # is not stable when using @threads

    # enamble network for geoid conversions
    Proj.enable_network!()

    # build transformation 
    if !threaded
        trans = Proj.Transformation(from_epsg, to_epsg; always_xy=true)

        # project points
        if x isa Vector{}
            xyh = trans.(x,y,height)
        else
            xyh = trans(x,y,height)
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
    epsg2epsg(x, y, from_epsg, to_epsg)

Returns `x`, `y` in `to_epsg` projection
"""
function epsg2epsg(
    x::Union{Vector{<:Number}, Number}, 
    y::Union{Vector{<:Number}, Number}, 
    from_epsg::String,
    to_epsg::String;
    parse_output = true,
    threaded = true
    )

     # build transformation 
    if !threaded
        trans = Proj.Transformation(from_epsg, to_epsg; always_xy=true)

        # project points
        if x isa Vector{}
            xy = trans.(x,y)
        else
            xy = trans(x,y)
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
    decimalyear(datetime::DateTime)

converts DateTime to decimalyear [e.g. decimalyear(DateTime("2018-10-24")) = 2018.813698630137]
"""
function decimalyear(datetime)

    # check if date is a leapyear
    if isleapyear(datetime)
        number_of_days = 366
    else
        number_of_days = 365
    end

    # calculate decimal year
    decyear = year(datetime) + dayofyear(datetime) / number_of_days

    return decyear
end

"""
    regular_grid(center_extent, node_spacing; node_width = node_spacing; node_shape = "square")

Returns a named tuple of node centers (`x_node_center`, `y_node_center`) and `mode_half_width`. 
`node_shape` is also provided
"""
function regular_grid(center_extent, node_spacing; node_width = node_spacing, node_shape = "square")
    x_node_center = center_extent.x_min:node_spacing:center_extent.x_max
    y_node_center = center_extent.y_min:node_spacing:center_extent.y_max
    node_half_width = node_width/2

    node_center = ndgrid(y_node_center, x_node_center)
    (; node_center, x_node_center, y_node_center, node_half_width, node_shape)
end


"""
    regular_grid_extents(x_center, y_center, half_width)

Returns the extents of a box with  width 2 x `half_width` centered at `x_center`, `y_center`. 
"""
function regular_grid_extents(x_center, y_center, half_width)

    (x_min = x_center - half_width, y_min = y_center - half_width, 
    x_max = x_center + half_width, y_max = y_center + half_width)
end


"""
    bin(x::Vector{<:Number}, y::Array{<:Number}, xbin_edges::Union{Vector{<:Number},StepRangeLen{}}; method::F = mean) where {F<:Function}

Fast binning of `y` as a function of `x`. Returns `x_binned`, `y_binned`, `bin_count`. `method` 
specifies approached used for aggregating binned values. NOTE: `x_binned` == 0 and 
`y_binned` == 0 when `bin_count` == 0.

# see https://discourse.julialang.org/t/performance-optimization-of-a-custom-binning-funciton/91616
"""
function bin(
    x::Vector{<:Number}, 
    y::Array{<:Number}, 
    xbin_edges::Union{Vector{<:Number},StepRangeLen{}};
    method::F = mean,
    ) where {F<:Function}

    # find bin breakpoints
    p = sortperm(vcat(x, xbin_edges))
    bins = findall(>(length(x)), p)

    # initialize outputs
    bin_count = Int.(diff(bins).-1)
    x_binned = zeros(eltype(x), length(bin_count))
    y_binned = zeros(eltype(y), length(bin_count), size(y, 2))

    # calculate binned metrics
    for i = findall(bin_count .> 0)
        x_binned[i] = method(@view(x[p[bins[i]+1:bins[i+1]-1]]))
        y_binned[i] = method(@view(y[p[bins[i]+1:bins[i+1]-1]]))
    end
    
    return x_binned, y_binned, bin_count 
end

"""
    bin(
    x::Vector{<:Number}, y::Vector{<:Number}, z::Array{<:Number}, 
    xbin_edges::Union{Vector{<:Number},StepRangeLen{}}, 
    ybin_edges::Union{Vector{<:Number},StepRangeLen{}}; 
    method::F = mean) where {F<:Function}
    )

Fast binning of `z` as a function of `x` and `y`. Returns `x_binned`, `y_binned`, `z_binned`, 
`bin_count`. `method` specifies approached used for aggregating binned values. NOTE: 
`x_binned` == 0, `y_binned` == 0 and z_binned == 0 when `bin_count` == 0.

# see https://discourse.julialang.org/t/performance-optimization-of-a-custom-binning-funciton/91616
"""
function bin(
    x::Vector{<:Number}, 
    y::Vector{<:Number}, 
    z::Vector{<:Number},
    xbin_edges::Union{Vector{<:Number},StepRangeLen{}},
    ybin_edges::Union{Vector{<:Number},StepRangeLen{}};
    method::F = median
    ) where {F<:Function}

    # find bin breakpoints
    py = sortperm(vcat(y, ybin_edges))
    binsy = findall(>(length(y)), py)
    biny_count = Int.(diff(binsy).-1)

    x_binned = zeros(eltype(x), length(ybin_edges)-1, length(xbin_edges)-1)
    y_binned = zeros(eltype(y), length(ybin_edges)-1, length(xbin_edges)-1)
    z_binned = zeros(eltype(z), length(ybin_edges)-1, length(xbin_edges)-1)
    bin_count =  zeros(Int32, length(ybin_edges)-1, length(xbin_edges)-1)

    # loop for each x bin
    Threads.@threads for i = findall(biny_count .> 0)
        y0 = @view y[py[binsy[i]+1:binsy[i+1]-1]]
        x0 = @view x[py[binsy[i]+1:binsy[i+1]-1]]
        z0 = @view z[py[binsy[i]+1:binsy[i+1]-1]]

        # find bin breakpoints
        px = sortperm(vcat(x0, xbin_edges))
        binsx = findall(>(length(x0)), px)

        # initialize outputs
        bin_count[i,:] .= Int.(diff(binsx).-1)

        # calculate binned metrics
        for j = findall(bin_count[i,:] .> 0)
            x_binned[i,j] = method(@view x0[px[binsx[j]+1:binsx[j+1]-1]])
            y_binned[i,j] = method(@view y0[px[binsx[j]+1:binsx[j+1]-1]])
            z_binned[i,j] = method(@view z0[px[binsx[j]+1:binsx[j+1]-1]])
        end
    end

    return x_binned, y_binned, z_binned, bin_count 
end


"""
    binnedfiltering(x::Vector{<:Number}, y::Vector{<:Number}, xbin_edges::Union{Vector{<:Number},StepRangeLen{}}; method:: = :mad)

Fast filtering of data within descrete bins of `y` as a funciton `x`. returns `x_binned`, `y_binned`, `bin_count`. `method` specifies approached used for aggregating binned values. 
"""
function binnedfiltering(
        x::Vector{<:Number}, 
        y::Vector{<:Number}, 
        xbin_edges::Union{Vector{<:Number},StepRangeLen{}};
        method::F = madnorm,
        threshold::Number = 10,
        ) where {F<:Function}

    # find bin breakpoints
    p = sortperm(vcat(x, xbin_edges))
    bins = findall(>(length(x)), p)

    # initialize outputs
    outlier = falses(size(y))
    bin_count = Int.(diff(bins).-1)

    # identify outliers
    for i = findall(bin_count .> 0)
        outlier[p[bins[i]+1:bins[i+1]-1]] = method(@view(y[p[bins[i]+1:bins[i+1]-1]])) .> threshold 
    end

    return outlier
end    


"""
    madnorm(x)

Returns `x_madnorm` the normalized median absolute deviation of `y` from the global median
"""
function madnorm(x)
    x_abs = abs.(x .- median(x))
    x_madnorm = x_abs ./ median(x_abs)
    return x_madnorm
end

"""
    mad(x)

Returns the median absolute deviation from the median
"""
function mad(x)
    if isempty(x)
        mad = NaN
    else
        mad = median(abs.(x .- median(x)))
    end
    return mad
end


"""
    range(x)

Return the range (maximum minus minimum) of x
"""
function range(x)
    x_minmax = extrema(x)
    r = x_minmax[2] - x_minmax[1]
    return r
end


"""
    utm_epsg(lon::Real, lat::Real, always_xy=true)

returns the EPSG code for the intersecting universal transverse Mercator (UTM) zone -OR- 
the relevant polar stereographic projection if outside of UTM limits.

modified from: https://github.com/JuliaGeo/Geodesy.jl/blob/master/src/utm.jl    
"""
function utm_epsg(lon::Real, lat::Real; always_xy=true)

    if !always_xy
        lat, lon = (lon, lat)
    end

    if lat > 84
        # NSIDC Sea Ice Polar Stereographic North
        epsg = 3413
        epsg = "EPSG:$epsg"
    elseif lat < -80
        # Antarctic Polar Stereographic
        epsg = 3031
        epsg = "EPSG:$epsg"
    end

    # make sure lon is from -180 to 180
    lon = lon - floor((lon+180) / (360)) * 360

    # int versions
    ilat = floor(Int64, lat)
    ilon = floor(Int64, lon)

    # get the latitude band
    band = max(-10, min(9,  fld((ilat + 80), 8) - 10))

    # and check for weird ones
    zone = fld((ilon + 186), 6)
    if ((band == 7) && (zone == 31) && (ilon >= 3)) # Norway
        zone = 32
    elseif ((band == 9) && (ilon >= 0) && (ilon < 42)) # Svalbard
        zone = 2 * fld((ilon + 183), 12) + 1
    end

    if lat >= 0
        epsg = 32600 + zone
    else
        epsg = 32700 + zone
    end

    # convert to proj string
    epsg = "EPSG:$epsg"
    return epsg
end


"""
    geoid(geoid::String; folder::String = "~/Documents/geoids")
returns geoid::GeoArray{} height file [relative to WGS 84 (EPSG::4979) ellipsoid]

If geoid file not found in `folder` it will be downloaded from:
https://www.agisoft.com/downloads/geoids/ 

Models Converted from USA NGA data by agisoft under Public Domain license.
"""
function geoid(geoid::String; folder::String="~/Documents/geoids")

    if geoid == "egm84"
        # EGM84 30' geoid model
        # WGS 84 (EPSG::4979) to EGM84 height (EPSG::5798)
        geoid_fn = "us_nga_egm84_30.tif"
    elseif geoid == "egm96"
        # EGM96 15' geoid model
        # WGS 84 (EPSG::4979) to EGM96 height (EPSG::5773
        geoid_fn = "us_nga_egm96_15.tif"
    elseif geoid == "egm2008"
        # EGM2008 1' geoid model
        # WGS 84 (EPSG::4979) to EGM2008 height (EPSG::3855)
        geoid_fn = "us_nga_egm2008_1.tif"
    else
        error("geoid not recognized, valid geoid names are \"egm84\", \"egm96\" and \"egm2008\"")
    end

    path2goid = joinpath(folder, geoid_fn)

    # download if file does not exist
    if !isfile(path2goid)
        if !isdir(folder)
            error("goids folder does not exist: $folder")
        else
            url = joinpath("https://s3-eu-west-1.amazonaws.com/download.agisoft.com/gtg/", geoid_fn)
            printstyled("local copy of $geoid file not found, downloading from: $url \n"; color = :blue, bold = true)
            HTTP.download(url, path2goid)
        end
    end

    # not sure if this Type decleration helps at all, feel free to delete
    return GeoArrays.read(path2goid; masked=false)
end

""" 
    dem_height(lon, lat, dem; filter_halfwidth::Union{Number,Nothing}= nothing)
Returns dem `height` for `lon`, `lat` locations. If filter_halfwidth is supplied, a 2D Gaussain
filter with sigma = `filter_halfwidth` is applied to the `dem` prior to sampeling
"""
function dem_height(
    lon, 
    lat, 
    dem; 
    filter_halfwidth::Union{Number,Nothing}= nothing,
    filter_kernel = :gaussian,
    slope = false
    )
    
    # use gaussian sampeling filter only if the center pixel contains < filter_thresh of the signal
    filter_thresh = 0.85;
    interp = Linear() # Constant(), Linear(), Cubic(), Quadratic()
    deminfo, goids_folder  = dem_info(dem)

    # check if within extents
    if (deminfo.extent.min_x != -180.) || 
        (deminfo.extent.min_y != -90.) || 
        (deminfo.extent.max_x != 180.) || 
        (deminfo.extent.max_y != 90.)

        ind = within.(Ref(deminfo.extent), lon, lat)

        if !any(ind)
            wgs84_height = fill(NaN32, size(lat));
            return wgs84_height
        end
    end
        
    # extract data from DEMs
    dem_ga = GeoArrays.read(deminfo.fn; masked=false)

    # determine approximate x and y grid size
    # check if geographic
    isgeographic = contains(dem_ga.crs.val, "AXIS[\"Latitude\",NORTH]")

    if !isnothing(filter_halfwidth)
        if isgeographic
            (a, b) = extrema(lat)
            (c, d) = extrema(lon)
            cent = (x = (c+d)/2, y = (a+b)/2)

            delta = 0.01
            x_lla = LLA(cent.y, cent.x, 0.0) 
            y_lla = LLA(cent.y+delta, cent.x, 0.0) 
            dist_lat = euclidean_distance(x_lla, y_lla) / delta           

            x_lla = LLA(cent.y, cent.x, 0.0) 
            y_lla = LLA(cent.y, cent.x+delta, 0.0) 
            dist_lon = euclidean_distance(x_lla, y_lla) / delta

            dist = (x = dem_ga.f.linear[1,1] * dist_lon, y = dem_ga.f.linear[2,2] * dist_lat)
        else
            dist = (x = dem_ga.f.linear[1,1], y = dem_ga.f.linear[2,2])
        end

        if filter_kernel == :gaussian
            # build gaussain sampeling kernel
            k = Kernel.gaussian([abs(filter_halfwidth/dist.x), abs(filter_halfwidth/dist.y)]);
        elseif filter_kernel == :average
            #round to the nearest odd integer
            dr = round(Int16, abs((filter_halfwidth)/dist.x)*2 + 1) # for GeoArrays x is rows
            dc = round(Int16, abs((filter_halfwidth)/dist.y)*2 + 1) # for GeoArrays y is columns
            k = centered(ones(dr, dc)./(dr*dc))
        end

        # only apply kernel if needed
        if k[0, 0] > filter_thresh
            k = nothing
        end
    else
        k = nothing
    end

    if !slope
        if deminfo.hight_datum == "wgs84"
            wgs84_height = pointextract(lon, lat, "EPSG:4326", dem_ga; 
            nodatavalue = deminfo.nodatavalue, filter_kernel = k, interp = interp, 
            replace_nodatavalue_withnans = true)
        else
            # convert from geoid to wgs84 height
            geoid_height = pointextract(lon, lat, "EPSG:4326", dem_ga; 
            nodatavalue = deminfo.nodatavalue, filter_kernel = k, interp = interp, 
            replace_nodatavalue_withnans = true)
            geoid_to_wgs84 = geoid(deminfo.hight_datum; folder = goids_folder)
            geoid_to_wgs84 = pointextract(lon, lat, "EPSG:4326", geoid_to_wgs84; 
            nodatavalue = deminfo.nodatavalue, filter_kernel = k, interp = interp, 
            replace_nodatavalue_withnans = true)
            wgs84_height = geoid_height .+ geoid_to_wgs84;

        end
        return wgs84_height
    else
        # form Gaussain sample kernel
        dx = centered([-1.0 -2.0 -1.0; 0.0 0.0 0.0; 1.0 2.0 1.0]) ./ (dem_ga.f.linear[1, 1]* 8) # for GeoArrays x is rows
        dy = centered([-1.0 0.0 1.0; -2.0 0.0 2.0; -1.0 0.0 1.0]) ./ (dem_ga.f.linear[2, 2] * 8) # for GeoArrays y is columns
        dk = [dx, dy]

        if deminfo.hight_datum == "wgs84"
            wgs84_height_slope = pointextract(lon, lat, "EPSG:4326", dem_ga; 
            nodatavalue = deminfo.nodatavalue, filter_kernel = k, derivative_kernel = dk, 
            interp = interp, replace_nodatavalue_withnans = true)
        else
            # convert from geoid to wgs84 height
            geoid_height = pointextract(lon, lat, "EPSG:4326", dem_ga; 
                nodatavalue = deminfo.nodatavalue, filter_kernel = k, 
                derivative_kernel = dk, replace_nodatavalue_withnans = true, interp = interp)
            
            geoid_to_wgs84 = geoid(deminfo.hight_datum; folder = goids_folder)
            geoid_to_wgs84 = pointextract(lon, lat, "EPSG:4326", geoid_to_wgs84; 
            nodatavalue = deminfo.nodatavalue, filter_kernel = k, interp = interp)
            geoid_height[1] = geoid_height[1] .+ geoid_to_wgs84;
            wgs84_height_slope = geoid_height;
        end
        
        return wgs84_height_slope
    end
end

"""
    geotile_extract_dem(geotile, geotile_dir, dem; filter_halfwidth, filter_kernel, slope, job_id, force_remake)
"""
function geotile_extract_dem(
    geotile_id::String,
    geotile_dir::String,
    dem::Symbol;
    filter_halfwidth::Union{Number,Nothing}=nothing,
    filter_kernel = :gaussian,
    slope = false,
    job_id = "",
    xoffset::Real = 0, # add offset to altimetry locations in local UTM coordinates [m]
    yoffset::Real = 0, # add offset to altimetry locations in local UTM coordinates [m]
    force_remake=false
    )

    path2geotile = joinpath(geotile_dir, geotile_id * ".arrow")
    outfile = replace(path2geotile, ".arrow" => ".$dem")

    if isfile(outfile) && !force_remake
        printstyled("    ->$job_id $geotile_id $dem already exists, skipping\n"; color = :green)
    elseif isfile(path2geotile)
        t1 = time(); 

        # use this temporarily to get rid of empty datacubes
        df = DataFrame(Arrow.Table(path2geotile))
        if isempty(df)
            printstyled("    ->$job_id $geotile_id had an empty DataFrame, deleting\n"; color = :red, bold = true)
            rm(path2geotile)
            return
        end
            
        #  stat(infile).size
        df = select!(df, [:longitude, :latitude])
        lon = reduce(vcat, df[:, :longitude]);
        lat = reduce(vcat, df[:, :latitude]);
   
        # this code is for  x and y offset testing
        if !isequal(0,xoffset) || !isequal(0,yoffset)
            epsg = utm_epsg(mean(lon), mean(lat); always_xy=true)
            trans = Proj.Transformation("EPSG:4326", epsg, always_xy=true)
            xy = trans.(lon,lat)
            xy = inv(trans).((getindex.(xy, 1) .+ xoffset), (getindex.(xy, 2) .+ yoffset))
            lon = getindex.(xy, 1);
            lat = getindex.(xy, 2);
        end

        if isempty(lon)
            printstyled("    ->$job_id $geotile_id is all empty, skipping\n"; color = :light_red)
            return
        end
        
        deminfo, _  = dem_info(dem)

        # check if within extents
        if (deminfo.extent.min_x != -180.) || (deminfo.extent.min_y != -90.) || (deminfo.extent.max_x != 180.) || (deminfo.extent.max_y != 90.)
            ind = within.(Ref(deminfo.extent), lon, lat)
            if !any(ind)
                printstyled("    ->$job_id $geotile_id $dem outside of dem limits, skipping\n"; color = :light_red)
                return
            end
        end

        if isnothing(filter_halfwidth)
            height = dem_height(lon, lat, dem, slope)
        else
            height = dem_height(lon, lat, dem; filter_halfwidth=filter_halfwidth, filter_kernel=filter_kernel, slope)
        end

        if !slope
            if all(isnan.(height))
                printstyled("    ->$job_id $geotile_id $dem no valid heights, skipping\n"; color = :light_red)
                return
            else   
                df[!,:height] = [Array{eltype(height)}(undef,size(row)) for row in df[:, :longitude]]   
                start = 1;
                for row = eachrow(df)
                    stop = length(row.height) + start - 1;
                    row.height = height[start:stop]
                    start = stop + 1
                end
            end
        else
            if all(isnan.(height[1]))
                printstyled("    ->$job_id $geotile_id $dem no valid heights, skipping\n"; color = :light_red)
                return
            else
                df[!,:height] = [Array{eltype(height[1])}(undef,size(row)) for row in df[:, :longitude]]
                df[!,:dhdx] = [Array{eltype(height[2])}(undef,size(row)) for row in df[:, :longitude]]
                df[!,:dhdy] = [Array{eltype(height[3])}(undef,size(row)) for row in df[:, :longitude]]   
                start = 1;

                for row = eachrow(df)
                    stop = length(row.height) + start - 1;
                    row.height = height[1][start:stop]
                    row.dhdx = height[2][start:stop]
                    row.dhdy = height[3][start:stop]
                    start = stop + 1
                end
            end

            tmp = tempname(dirname(outfile))
            Arrow.write(tmp, df::DataFrame);
            mv(tmp, outfile; force = true)

            total_time = round((time()-t1)/60, digits = 2);
            printstyled("    ->$job_id $geotile_id $dem extracted: $(total_time) min \n"; color = :light_black)
        end
    else
        println(path2geotile)
        printstyled("    ->$job_id $geotile_id does not exist\n"; color = :yellow)
        return
    end
end

"""
    geotile_extract_dem(geotile, geotile_dir, dem; filter_halfwidth, force_remake)
"""
function geotile_extract_dem(
    geotiles::DataFrame,
    geotile_dir::String,
    dem::Symbol; 
    filter_halfwidth::Union{Number,Nothing} = nothing, 
    filter_kernel=:gaussian,
    slope = false,
    job_id = "",
    xoffset=0, # add offset to altimetry locations in local UTM coordinates [m]
    yoffset=0, # add offset to altimetry locations in local UTM coordinates [m]
    force_remake = false
    )

    printstyled("extracting $dem heights for $job_id geotiles\n"; color = :blue, bold = true)

    #asyncmap(eachrow(geotiles); ntasks = 4) do geotile
    for geotile in eachrow(geotiles)
        geotile_extract_dem(geotile.id, geotile_dir, dem; 
            filter_halfwidth=filter_halfwidth, filter_kernel=filter_kernel, slope=slope,
            job_id = job_id, xoffset=xoffset, yoffset=yoffset, force_remake = force_remake)
    end
end

"""
    geotile_extract_dem(products, dems, geotiles, paths; slope, force_remake)
"""
function geotile_extract_dem(
    products::NamedTuple,
    dems::Vector{Symbol},
    geotiles::DataFrame,
    paths::NamedTuple;
    slope = false,
    xoffset=0, # add offset to altimetry locations in local UTM coordinates [m]
    yoffset=0, # add offset to altimetry locations in local UTM coordinates [m]
    force_remake = false
    )

    for product in products
        for dem in dems
            geotile_extract_dem(geotiles, paths[product.mission].geotile, dem;
                filter_halfwidth=product.halfwidth, filter_kernel=product.kernel, slope=slope,
                job_id=product.mission, xoffset=xoffset, yoffset=yoffset, force_remake = force_remake)
        end
    end
end


"""
    geotile_aggrigate_reduce(geotiles::DataFrame, geotile_dir::String, 
        vars::Union{Vector{String}, String}; extension::String = ".arrow")

Aggrigate variables in `vars` for list of `geotiles`` and reduce to `Vectors` of `elements`.
subset_fraction [0 to 1] allows for random subsampeling of the data  
"""
function geotile_aggrigate_reduce(
    geotiles::DataFrame,
    geotile_dir::String,
    vars::Vector{String};
    extension::String = ".arrow",
    subset_fraction = 1.
    )

    # printstyled("aggrigating and reducing $extension geotiles\n"; color = :blue, bold = true)
    df = DataFrame()
    #asyncmap(eachrow(geotiles); ntasks = 4) do geotile
    for geotile in eachrow(geotiles)
        path2geotile = joinpath(geotile_dir, "$(geotile.id)$(extension)")
        if isfile( path2geotile)
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
        df0[!,var] = vcat(df[ind, var]...)
    end

    # remove empty rows
    delete(df0, isempty.(df0[:, first(vars)]))

    return df0
end

"""
    geotile_aggrigate_reduce(geotiles::DataFrame, geotile_dir::String, 
        vars::Union{Vector{String}, String}; extension::String = ".arrow")

Aggrigate variables in `vars` for list of `geotiles`` and reduce to `Vectors` of `elements`. 
subset_fraction [0 to 1] allows for random subsampeling of the data 
"""
function geotile_aggrigate_reduce(
    geotiles::DataFrameRow{DataFrame, DataFrames.Index},
    geotile_dir::String,
    vars::Union{Vector{String},Vector{Symbol}};
    extension::String = ".arrow",
    subset_fraction = 1.
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
                df0[!,var] = vcat(df[ind, var]...)
            else
                #if not expand to a vector 
                df0[!,var] = vcat([fill(df[k, var], size(df[k, refarray])) for k = eachindex(df[ind, var])]...)
            end
        end

        # remove empty rows
        delete!(df0, isempty.(df0[:, last(vars)]))
        return df0
    end
end

"""
    itslive_proj!(df::DataFrame{}; height = nothing)
add x and y coodinates for local itslive projection to the DataFrame `df`
"""
function itslive_proj!(df; height = nothing)
    epsg = itslive_epsg(mean(df.longitude), mean(df.latitude); always_xy=true)
    
    if isnothing(height)
        df[!, :X], df[!, :Y] = epsg2epsg(df.longitude, df.latitude, "EPSG:4326", epsg, parse_output = true)
    else
        df[!, :X], df[!, :Y], df[!, :H]= epsg2epsg(df.longitude, df.latitude, height, epsg, parse_output = true)
    end
    
    return df, epsg
end

"""
    itslive_epsg(lon, lat)
Return epsg code for the ITS_LIVE projection
"""
function itslive_epsg(longitude, latitude; always_xy = true)

    if !always_xy
        latitude, longitude = (longitude, latitude)
    end

    if latitude > 55
        # NSIDC Sea Ice Polar Stereographic North
        return epsg = 3413
    elseif latitude < -56
        # Antarctic Polar Stereographic
        return epsg = 3031
    end

    # make sure lon is from -180 to 180
    lon = longitude - floor((longitude+180) / (360)) * 360

    # int versions
    ilat = floor(Int64, latitude)
    ilon = floor(Int64, lon)

    # get the latitude band
    band = max(-10, min(9,  fld((ilat + 80), 8) - 10))

    # and check for weird ones
    zone = fld((ilon + 186), 6)
    if ((band == 7) && (zone == 31) && (ilon >= 3)) # Norway
        zone = 32
    elseif ((band == 9) && (ilon >= 0) && (ilon < 42)) # Svalbard
        zone = 2 * fld((ilon + 183), 12) + 1
    end

    if latitude >= 0
        epsg = 32600 + zone
    else
        epsg = 32700 + zone
    end

    # convert to proj string
    epsg = "EPSG:$epsg"
    return epsg
end


"""
    itslive_zone(lon, lat; always_xy = true)
Return the utm `zone` and `isnorth` variables for the ITS_LIVE projection
"""
function itslive_zone(lon, lat; always_xy = true)
    if !always_xy
        lat, lon = (lon, lat)
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
    itslive_paramfiles(lon, lat; gridsize = 240, path2param = "/Users/...", always_xy = true)
Return paths to its_live parameter files
"""
function itslive_paramfiles(
    lon, 
    lat; 
    gridsize = 240, 
    path2param = pathlocal.itslive_parameters, 
    always_xy = true)

    zone, isnorth = itslive_zone(lon, lat; always_xy = always_xy);

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
        ROI = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"ROI")),
            binary = true
            ),

        StableSurface = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"StableSurface")),
            binary = true
            ),

        dhdx = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"dhdx")),
            binary = false
            ),

        dhdxs = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"dhdxs")),
            binary = false
            ),

        dhdy = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"dhdy")),
            binary = false
            ),

        dhdys = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"dhdys")),
            binary = false
            ),

        floatingice = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"floatingice")),
            binary = true
            ),

        glacierice = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"glacierice")),
            binary = true
            ),

        h = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"h")),
            binary = false
            ),

        inlandwater = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"inlandwater")),
            binary = true
            ),

        land = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"land")),
            binary = true
            ),

        landice = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"landice")),
            binary = true
            ),

        landice_2km_inbuff = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"landice_2km_inbuff")),
            binary = true
            ),

        ocean = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"ocean")),
            binary = true
            ),

        region = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"region")),
            binary = false
            ),

        sp = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"sp")),
            binary = true
            ),

        thickness = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"thickness")),
            binary = false
            ),

        vx = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"vx")),
            binary = false
            ),

        vx0 = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"vx0")),
            binary = false
            ),

        vxSearchRange = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"vxSearchRange")),
            binary = false
            ),

        vy = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"vy")),
            binary = false
            ),

        vy0 = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"vy0")),
            binary = false
            ),

        vySearchRange = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"vySearchRange")),
            binary = false
            ),

        xMaxChipSize = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"xMaxChipSize")),
            binary = false
            ),

        xMinChipSize = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"xMinChipSize")),
            binary = false
            ),

        yMaxChipSize = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"yMaxChipSize")),
            binary = false
            ),

        yMinChipSize = (
            path = joinpath(path2param, @sprintf("%s_%s_%s.tif", region, grid,"yMinChipSize")), 
            binary = false
            )
    )
    
    return paramfiles
end

"""
    itslive_extract(lon, lat, varname; gridsize = 240, path2param = "/Users/...", always_xy = true)
Return point values from its_live parameter files
"""
function itslive_extract(
    lon::Vector, 
    lat::Vector,
    vars::Union{Symbol,Vector{Symbol}}; 
    gridsize = 240,
    path2param = pathlocal.itslive_parameters, 
    always_xy = true)

    # get regional parameter files
    paramfiles = itslive_paramfiles(mean(lon), mean(lat); 
        gridsize = gridsize, path2param = path2param, always_xy = always_xy)
    
    for v in vars
        if !haskey(paramfiles, v)
            error(error("$(variable) is not a valid variable, must be one of $(keys(paramfiles))"))
        end
    end

    # create projeciton transfomation function
    # NOTE: not worth checking for same projection... faster just to do it this way
    v = vars[1]
    ga = GeoArrays.read(paramfiles[v].path, masked=false)
    xy = epsg2epsg(lon, lat, "EPSG:4326", ga.crs.val; parse_output = false)

    # find xy extents
    x = getindex.(xy, 1);
    y = getindex.(xy, 2)
    minmax_x = extrema(x)
    minmax_y = extrema(y)

    extent = (min_x = minmax_x[1], min_y = minmax_y[1], max_x = minmax_x[2] + abs(ga.f.linear[1,1]), max_y = minmax_y[2] + abs(ga.f.linear[2,2])) 

    if !bbox_overlap(bbox(ga), extent)
        # extents do not overlap
        for v in vars
            ga = GeoArrays.read(paramfiles[v].path, masked=false)

            if paramfiles[v].binary
                df[!,v] = fill(false, length(x))
            else
                df[!,v] = fill(convert(eltype(ga.A), nodatavalue), length(x))
            end
        end
        return df
    else
        df = DataFrame()
        v = vars[1]
        ga0 = crop(ga, extent)
        ind = CartesianIndex.(Tuple.(indices.(Ref(ga0), xy)))
        df[!,v] = ga0[ind]
        
        for v in vars[2:end]
            ga = GeoArrays.read(paramfiles[v].path, masked=false)
            ga0 = crop(ga, extent)

            if paramfiles[v].binary
                df[!,v] = (ga0[ind] .== 1)
            else
                df[!,v] = ga0[ind]
            end
        end

        return df
    end
end
    
"""
    info, goids_folder  = dem_info()
Return dem information and goids_folder
"""
function dem_info(dem::Symbol)
    goids_folder = pathlocal.geoid_dir

    if dem == :cop30_v1
        info = (
            fn = pathlocal.cop30_v1,
            nodatavalue = 0,
            hight_datum = "egm2008",
            epsg = 4326,
            res = 30,
            extent = (min_x = -180., min_y = -90., max_x = 180., max_y = 90.)
        )
    elseif dem == :cop30_v2
        info = (
            fn = pathlocal.cop30_v2,
            nodatavalue = 0,
            hight_datum = "egm2008",
            epsg = 4326,
            res = 30,
            extent = (min_x = -180., min_y = -90., max_x = 180., max_y = 90.),
        )
    elseif dem == :nasadem_v1
        info = (
            fn = pathlocal.nasadem_v1,
            nodatavalue = -32767,
            hight_datum = "wgs84",
            epsg = 4326,
            res = 30,
            extent = (min_x = -180., min_y = -56., max_x = 180., max_y = 60.)
        )
    elseif dem == :arcticdem_v3_10m
        info = (
            fn = pathlocal.arcticdem_v3_10m,
            nodatavalue = -9999.0,
            hight_datum = "wgs84",
            epsg = 3413,
            res = 10,
            extent = (min_x = -180., min_y = 50., max_x = 180., max_y = 90.)
        )
    elseif dem == :rema_v2_10m
        info = (
            fn = pathlocal.rema_v2_10m,
            nodatavalue = -9999.0,
            hight_datum = "wgs84",
            epsg = 3031,
            res = 10,
            extent = (min_x = -180., min_y = -90., max_x = 180., max_y = -61.)
        )
    else    
        error("dem not recognized, valid dem names are :cop30_v1, :nasadem_v1, :rema_v2_10m, and :arcticdem_v3_10m")
    end

    return info, goids_folder
end

"""
    normalize(data)
Return a nomalized version of data
"""
function normalize(data)
    m = mean(data);
    s = std(data);
    n = (data .- m) ./ s

    datan = (std = s, mean = m, norm = n)
    
    return datan
end


"""
    centroid(ga::GeoArray)
return cetroid of GeoArray
"""
function centroid(ga::GeoArray)
    bbox = GeoArrays.bbox(ga)
    (x = (bbox.min_x + bbox.max_x)/2, y = (bbox.min_y + bbox.max_y)/2)
end


"""
xy_gridsize(ga::GeoArray)
Return x and y grid size in meters
"""
function xy_gridsize(ga::GeoArray)

    # check if geographic
    isgeographic = contains(ga.crs.val, "AXIS[\"Latitude\",NORTH]")
    if isgeographic
        cent = centroid(ga)

        x_lla = LLA(cent.y, cent.x, 0.0) 
        y_lla = LLA(cent.y+1, cent.x, 0.0) 
        dist_lat = euclidean_distance(x_lla, y_lla)                   

        x_lla = LLA(cent.y, cent.x, 0.0) 
        y_lla = LLA(cent.y, cent.x+1, 0.0) 
        dist_lon = euclidean_distance(x_lla, y_lla)        

        gridsize = (x = ga.f.linear[1,1] * dist_lon, y = ga.f.linear[2,2] * dist_lat)
    else
        gridsize = (x = ga.f.linear[1,1], y = ga.f.linear[2,2])
    end
    return gridsize
end

"""
    dist_ll2xy(longitude, latitude)
Return local utm x and y distance in meters per degree latitude and longitude
"""
function dist_ll2xy(longitude, latitude)
    delta = 0.0001;
    lla = LLA.(latitude, longitude, Ref(0));
    lla_dlat = LLA.(latitude .+ delta, longitude, Ref(0));
    lla_dlon = LLA.(latitude, longitude .+ delta, Ref(0));

    zone, north = Geodesy.utm_zone(mean(latitude), mean(longitude))
    utm_from_lla = UTMfromLLA(zone, north, Geodesy.wgs84); 

    lonlat2xy = [Matrix{eltype(latitude[1])}(undef,(2,2)) for i in latitude]

    Threads.@threads for i in eachindex(lla)
        p = utm_from_lla(lla[i])::UTM{Float64}
        p_dlat = utm_from_lla(lla_dlat[i])::UTM{Float64}
        p_dlon = utm_from_lla(lla_dlon[i])::UTM{Float64}

        lonlat2xy[i] = [(p_dlon.x - p.x) (p_dlat.y - p.y);  
                    (p_dlat.x - p.x) (p_dlon.y - p.y)] ./ delta;
    end
    return lonlat2xy
end

"""
    dh_ll2xy(lon, lat, dhdlat, dhdlon)
Return slope in x and y UTM coodinates per meter from WGS84 lat lon slopes per degree
"""
function dh_ll2xy(lon, lat, dhdlon, dhdlat)
    delta = 0.0001;
    lla = LLA.(lat, lon, Ref(0));
    lla_dlat = LLA.(lat .+ delta, lon, Ref(0));
    lla_dlon = LLA.(lat, lon .+ delta, Ref(0));

    zone, north = Geodesy.utm_zone(mean(lat), mean(lon))
    utm_from_lla = UTMfromLLA(zone, north, Geodesy.wgs84); 

    dhdx = similar(lat)
    dhdy = similar(lon)

    Threads.@threads for i in eachindex(lla)
        p = utm_from_lla(lla[i])::UTM{Float64}
        p_dlat = utm_from_lla(lla_dlat[i])::UTM{Float64}
        p_dlon = utm_from_lla(lla_dlon[i])::UTM{Float64}

        a = (p_dlat.x - p.x) / delta
        b = (p_dlat.y - p.y) / delta

        c = (p_dlon.x - p.x) / delta
        d = (p_dlon.y - p.y) / delta

        dhdx[i] = ((dhdlon[i] * a)/sqrt(a.^2+b.^2) + (dhdlat[i] * c)/(c+d))/(a+c);
        dhdy[i] = ((dhdlon[i] * b)/sqrt(a.^2+b.^2) + (dhdlat[i] * d)/(c+d))/(b+d);
    end
    return dhdx, dhdy
end


"""
    dh_ll2aa(x, y, dhdx, dhdy, epsg_pts, epsg_dem)
Return slope in along-track and across-track in DEM coodinates per meter
"""
function dh_xy2aa(x, y, dhdx, dhdy, epsg_pts, epsg_dem)

    if epsg_pts !== epsg_dem
        x, y = epsg2epsg(x, y, "EPSG:4326", "EPSG:$epsg_dem"; parse_output = true)
    end

    dx = x[(begin+1):end] .- x[begin:(end-1)]
    push!(dx, dx[end])
    
    dy = y[(begin+1):end] .- y[begin:(end-1)]
    push!(dy, dy[end])

    theta = atan.(dx,dy)
    sn = sin.(theta);
    cs = cos.(theta);

    along  = dhdx .* sn .+ dhdy .* cs;
    across = -(dhdx .* cs  .- dhdy .* sn); 

    return along, across 
end


"""
    dh_ll2aa(lon, lat, dhdlat, dhdlon)
Return slope in along-track and across-track in UTM coodinates per meter from WGS84 lat lon slopes per degree
"""
function dh_ll2aa(lon, lat, dhdlon, dhdlat)

    if length(lon) == 1
        along = zeros(eltype(lon),1)
        across = zeros(eltype(lon),1)
        return along, across
    else
        delta = 0.0001;
        lla = LLA.(lat, lon, Ref(0));
        lla_dlat = LLA.(lat .+ delta, lon, Ref(0));
        lla_dlon = LLA.(lat, lon .+ delta, Ref(0));

        zone, north = Geodesy.utm_zone(mean(lat), mean(lon))
        utm_from_lla = UTMfromLLA(zone, north, Geodesy.wgs84); 

        dhdx = similar(lat)
        dhdy = similar(lon)
        p = Vector{UTM{Float64}}(undef,length(lat))

        Threads.@threads for i in eachindex(lla)
            p[i] = utm_from_lla(lla[i])::UTM{Float64}
            p_dlat = utm_from_lla(lla_dlat[i])::UTM{Float64}
            p_dlon = utm_from_lla(lla_dlon[i])::UTM{Float64}

            a = (p_dlat.x - p[i].x) / delta
            b = (p_dlat.y - p[i].y) / delta

            c = (p_dlon.x - p[i].x) / delta
            d = (p_dlon.y - p[i].y) / delta

            dhdx[i] = ((dhdlon[i] * a)/(a+b) + (dhdlat[i] * c)/(c+d))/(a+c);
            dhdy[i] = ((dhdlon[i] * b)/(a+b) + (dhdlat[i] * d)/(c+d))/(b+d);
        end

        dx = [k.x for k in  p]
        dx = dx[(begin+1):end] .- dx[begin:(end-1)]
        push!(dx, dx[end])
        
        dy = [k.y for k in  p]
        dy = dy[(begin+1):end] .- dy[begin:(end-1)]
        push!(dy, dy[end])

        theta = atan.(dx,dy)
        sn = sin.(theta);
        cs = cos.(theta);
        across = -(dhdx .* cs  .- dhdy .* sn); 
        along  = dhdx .* sn .+ dhdy .* cs;
        
        return along, across
    end
end

"""
    geotile_extract_masks(geotile, geotile_dir; vars, force_remake)
"""
function geotile_extract_mask(
    geotile_id::String,
    geotile_dir::String;
    vars::Vector{Symbol} = [:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean],
    geotile_ext::String = ".arrow",
    job_id="",
    force_remake = false
    )

    path2geotile = joinpath(geotile_dir, geotile_id * geotile_ext)
    outfile = replace(path2geotile, ".arrow" => ".masks")

    if isfile(outfile) && !force_remake
        printstyled("    ->$job_id $geotile_id masks already exists, skipping\n"; color=:green)
    elseif isfile(path2geotile)
        t1 = time(); 

        df0 = DataFrame(Arrow.Table(path2geotile))[!,[:longitude, :latitude]]
        vlength = length.(df0.longitude)

        lon = vcat(df0.longitude...)
        lat = vcat(df0.latitude...)

        if isempty(lon)
            printstyled("    ->$job_id $geotile_id is all empty, skipping\n"; color=:light_red)
            return
        end

        df0 = itslive_extract(lon, lat, vars, path2param = pathlocal.itslive_parameters)

        df = DataFrame()
        for v in vars
            df[!,v] = [Array{eltype(df0[!,v])}(undef,l) for l in vlength] 
        end

        start = 1
        for row = eachrow(df)
            stop = length(row[vars[1]]) + start - 1;
            for v in vars
                row[v] = df0[!,v][start:stop]
            end
            start = stop + 1
        end

        tmp = tempname(dirname(outfile))
        Arrow.write(tmp, df::DataFrame);
        mv(tmp, outfile; force = true)

        total_time = round((time()-t1)/60, digits = 2);
        printstyled("    ->$job_id $geotile_id masks extracted: $(total_time) min \n"; color=:light_black)
    else
        printstyled("    ->$job_id $geotile_id does not exist\n"; color=:yellow)
        return
    end
end

"""
    geotile_extract_masks(geotile, geotile_dir; vars, force_remake)
"""
function geotile_extract_mask(
    geotiles::DataFrame,
    geotile_dir::String;
    vars::Vector{Symbol} = [:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean],
    job_id="",
    force_remake = false
    )

    printstyled("extracting masks for $job_id geotiles\n"; color=:blue, bold=true)

    #asyncmap(eachrow(geotiles); ntasks = 4) do geotile
    for geotile in eachrow(geotiles)
        geotile_extract_mask(geotile.id, geotile_dir; vars = vars, job_id = job_id, force_remake = force_remake)
    end
end


"""
    _offset_design_matrix(dhdx, dhdy, p)
interanl function for building the design matrix to solve for track offset
"""
_offset_design_matrix(dhdx, dhdy) = hcat(dhdx, dhdy, ones(size(dhdx)))

"""
    track_offset(dhdx, dhdy, dh, p)
Given slope in x [dhdx], slope in y [dhdy], observed  and hieght anomaly [dh], 
solve for the  optimal x offset p[1], y offset p[2], height offset p[3]  
that minmize the distance to dh
"""
function track_offset(
    dhdx, 
    dhdy, 
    dh; 
    weights=nothing,
    regressor=HuberRegression(fit_intercept=false, scale_penalty_with_samples=false),
    interations=1, 
    iter_thresh=5, 
    verbose = true
    )
    
    valid = .!isnan.(dh) 
    θ = []
    mad_offset = 0.;
    mad_ref = 0.
    ddh = [];

    for i = 1:interations
        X = _offset_design_matrix(dhdx[valid], dhdy[valid])

        if isnothing(weights)
            θ = fit(regressor, X, dh[valid])
        else
            θ = fit(regressor, X .* weights[valid], dh[valid] .* weights[valid])
        end

        ddh = X * θ
        delta = dh[valid] .- ddh
        
        if verbose
            println("    [$i]: dx = $(round(θ[1]; digits=2)),  dy = $(round(θ[2]; digits=2)), dz = $(round(θ[3]; digits=2)), mad_offset = $(round(mad(delta); digits=2)) (mad_ref=$(round(mad(dh[valid]); digits=2))) (add to h in dh = h - h0)")
        end

        if i < interations
            valid[valid] = madnorm(delta) .< iter_thresh
        else
            mad_offset = mad(delta)
            mad_ref = mad(dh[valid])
        end
    end

    if verbose
        println("--------------------------------------------------------------------------------------------")
    end

    cnt = sum(valid)
    return θ[1], θ[2], θ[3], cnt, mad_offset, mad_ref
end

"""
    track_offset_dh(dhdx, dhdy, dx, dy, dz,)
Given slope in x [dhdx], slope in y [dhdy], the x offset [dx], 
y offset [dy], height offset [dz], calcualte the heigh anomaly [dh]
"""
function track_offset_dh(dhdx, dhdy, dx, dy, dz)
    dh =
        dx .* dhdx .+ dy .* dhdy .+ dz
    return dh
end

"""
    geotile_track_offset(geotile, geotile_dir, dem; vars, force_remake)
"""
function geotile_track_offset(
    geotile_id::String,
    geotile_dir::String,
    dem::Symbol;
    valid_count_thresh = 100,
    force_remake = false
    )

    path2geotile = joinpath(geotile_dir, geotile_id * ".arrow")
    outfile = replace(path2geotile, ".arrow" => ".$(dem)_offset")

    if isfile(outfile) && !force_remake
        printstyled("    -> $geotile_id $dem track offsets already exists, skipping\n"; color = :green)
    elseif isfile(path2geotile)
        t1 = time(); 

        # check that files exist
        fn_dem = replace(path2geotile, ".arrow" => ".$dem")
        fn_masks = replace(path2geotile, ".arrow" => ".masks")
        fn_canopyh = replace(path2geotile, ".arrow" => ".canopyh")

        if !isfile(fn_dem)
            printstyled("    -> $dem for $geotile_id does not exist\n"; color = :yellow)
            return
        elseif !isfile(fn_masks)
            printstyled("    -> masks for $geotile_id does not exist\n"; color = :yellow)
            return
        elseif !isfile(fn_canopyh)
            printstyled("    -> canopy for $geotile_id does not exist\n"; color = :yellow)
            return 
        end
     
        # read in dem data
        df = DataFrame(Arrow.Table(fn_dem))
        deminfo, goids_folder = Altim.dem_info(:arcticdem_v3_10m)

        # difference altimetry with dem
        df[!, :dh] = Arrow.Table(path2geotile)[:height] .- df.height

        # read in mask
        df[!, :land] = Arrow.Table(fn_masks)[:land]

        df[!, :canopyh] = Arrow.Table(fn_canopyh)[:canopyh]

        # initialize arrays
        vlength = length.(df.longitude)
        df0 = DataFrame()
        vars = [:above, :forward, :right, :count, :weight]
        for v in vars
            df0[!,v] = Vector{Float32}(undef, length(df.longitude))
        end

        vars = [:dh_along, :dh_across]
        for v in vars
            df0[!,v] = [Array{Float32}(undef,l) for l in vlength] 
        end

        # calculate offsets
        for (row, row0) in zip(eachrow(df), eachrow(df0))

            if isempty(row.dh)
                row0.above = NaN32
                row0.forward = NaN32
                row0.right = NaN32
                row0.count = 0
                continue
            elseif length(row.longitude) < 2
                row0.above = NaN32
                row0.forward = NaN32
                row0.right = NaN32
                row0.dh_along .= NaN32
                row0.dh_across .= NaN32
                row0.count = 0
                continue
            end 

            # rotate slope to across and along-track in units of m per m
            if deminfo.epsg == 4326
                row0.dh_along, row0.dh_across = dh_ll2aa(row.longitude, row.latitude, row.dhdx, row.dhdy)
            else
                row0.dh_along, row0.dh_across = dh_xy2aa(row.longitude, row.latitude, row.dhdx, row.dhdy, 4326, deminfo.epsg)
            end

            # filter NaNs and gross outliers and exclude ice and water
            valid = .!isnan.(row.dh) .& (abs.(row.dh) .< 100) .& row.land
            row0.count = sum(valid)

            if (row0.count < valid_count_thresh)
                row0.above = NaN32
                row0.forward = NaN32
                row0.right = NaN32
                row0.weight = NaN32
                continue
            end

            w = 1 ./ max.(row.canopyh, 0.01);
            row0.weight = sum(w[valid])

            # solve for optimal vertical, along and across track offsets 
            row0.right, row0.forward, row0.above, _, _, _ =
                track_offset(row0.dh_along[valid], row0.dh_across[valid], row.dh[valid]; weights = w[valid])
        end
    
        tmp = tempname(dirname(outfile))
        Arrow.write(tmp, df0::DataFrame);
        mv(tmp, outfile; force = true)

        total_time = round((time()-t1)/60, digits = 2);
        printstyled("    -> $geotile_id $dem offsets calculated: $(total_time) min \n"; color = :light_black)
    else
        printstyled("    -> $geotile_id does not exist\n"; color = :yellow)
        return
    end
end

"""
    geotile_track_offset(
        geotile, geotile_dir, dem; valid_count_thresh = 100, 
        force_remake = force_remake
    )
"""
function geotile_track_offset(
    geotiles::DataFrame,
    geotile_dir::String,
    dem::Symbol;
    valid_count_thresh = 100,
    force_remake = false
    )
    
    printstyled("calculating track offsets for geotiles\n"; color = :blue, bold = true)

    #asyncmap(eachrow(geotiles); ntasks = 4) do geotile
    for geotile in eachrow(geotiles)
        geotile_track_offset(geotile.id, geotile_dir, dem; valid_count_thresh = valid_count_thresh, force_remake = force_remake)
    end
end

"""
    pointextract(geotile, geotile_dir; var_name, force_remake)
"""
function pointextract(
    df0::DataFrame,
    ga::GeoArray;
    nodatavalue = 0.0,
    var_name::Symbol = :var
    )

    df = DataFrame()
    vlength = length.(df0.longitude)

    lon = reduce(vcat, df0.longitude)
    lat = reduce(vcat, df0.latitude)
    
    if isempty(lon)
        return df
    end

    var = pointextract(lon, lat, "EPSG:4326", ga; nodatavalue = nodatavalue)

    if all(var .== nodatavalue) || all(isnan.(var))
        return df
    else
        df[!,var_name] = [Vector{eltype(var)}(undef,l) for l in vlength] 

        start = 1
        for row = eachrow(df)
            stop = length(row[var_name]) + start - 1;
            row = var
            start = stop + 1
        end
        return df
    end
end

"""
    geotile_pointextract(geotiles, geotile_dir; var_name, force_remake)
"""
function geotile_pointextract(
    geotiles::DataFrame,
    geotile_dirs,
    ga::GeoArray;
    var_name = :var,
    nodatavalue = 0.0,
    job_ids = "",
    force_remake = false
    )
    
    printstyled("extracting $var_name for geotiles\n"; color = :blue, bold = true)
    geotile_dirs isa AbstractVector || (geotile_dirs = [geotile_dirs])
    job_ids isa AbstractVector || (job_ids = [job_ids])

    #asyncmap(eachrow(geotiles); ntasks = 4) do geotile
    for geotile in eachrow(geotiles)
        t1 = time();

         # check if ga is in EPSG:4326, if so then it's possible to check if lats and lons fall within ga 
        if contains(ga.crs.val, "\"EPSG\",\"4326\"")
            if !bbox_overlap(bbox(ga), geotile.extent)
                printstyled("    ->$(geotile.id) does not overlap $var_name, skipping\n"; color = :light_red)
                continue
            end
        end

        path2geotiles = joinpath.(geotile_dirs, Ref(geotile.id * ".arrow"))
        path2outfiles = replace.(path2geotiles, Ref(".arrow" => ".$var_name"))

        outfiles_exist = isfile.(path2outfiles)
        geotiles_exist = isfile.(path2geotiles) 

        extract = geotiles_exist .& (.!outfiles_exist .| force_remake)

        if any(extract)
            path2geotiles = path2geotiles[extract];
            path2outfile = path2outfiles[extract];
            job_id = job_ids[extract]

            df = DataFrame()
            stop = zeros(size(path2geotiles))
            for (i, path2geotile) in enumerate(path2geotiles)
                df = vcat(df, DataFrame(Arrow.Table(path2geotile))[!,[:longitude, :latitude]])
                stop[i] = nrow(df)
            end
            stop = Int.(stop)

            try
                df = pointextract(df, ga; var_name = var_name, nodatavalue = nodatavalue)
            catch
                return df 
            end
            if isempty(df)
                printstyled("    ->$(geotile.id) no valid $var_name, skipping\n"; color = :light_red)
            else

            start = Int(1);
            for (i, outifle) = enumerate(path2outfile)       
                    tmp = tempname(dirname(outifle))
                    Arrow.write(tmp, df[start:stop[i],:]::DataFrame);
                    mv(tmp, outifle; force = true)
                    start = stop[i]+1;
                end
            end
            total_time = round((time()-t1)/60, digits = 2);
            printstyled("    ->$job_id $(geotile.id) $var_name extracted: $(total_time) min \n"; color = :light_black)
        else
            printstyled("    ->$(geotile.id) does not exist or all $var_name already extracted, skipping\n"; color = :light_red)
        end 
    end
end

"""
    geotile_load(geotile_width, filesuffix, data_dir; extent::Extent = world)

Returns a single dataframe containing all geotiles [within extent]
"""
function geotile_load(geotile_width, filesuffix, data_dir; extent::Extent = world)

    # define geotiles
    geotiles = geotile_define(geotile_width)

    # geotiles within extent
    foo = nt2extent.(geotiles.extent)
    ind = Extents.intersects.(Ref(extent), foo)

    geotiles = geotiles[ind,:]

    df = DataFrame()

    for geotile in eachrow(geotiles)
        filename = joinpath(data_dir, "$(geotile.id).$(filesuffix)")

        if isfile(filename)
            df = append!(df,DataFrame(Arrow.Table(filename)))
        end
    end
    return df
end

"""
    nt2extent(nt::NamedTuple{(:min_x, :min_y, :max_x, :max_y), NTuple{4, <:Number}})

converts a NamedTuple to an Extent
"""
function nt2extent(nt::NamedTuple{(:min_x, :min_y, :max_x, :max_y)})

    Extent(X = (nt.min_x, nt.max_x), Y = (nt.min_y, nt.max_y))
end


"""
    geotile_subset(geotiles, extent::Extent)

returns subset of geotiles that intersect extent
"""
function geotile_subset(geotiles, extent::Extent)
        nt2extent.(geotiles.extent);
        ind = .!isnothing.(Extents.intersect.(Ref(extent), nt2extent.(geotiles.extent)));
        geotiles = geotiles[ind,:];
        return geotiles
end


"""
    geotile_subset!(geotiles, extent::Extent)

returns subset of geotiles that intersect extent
"""
function geotile_subset!(geotiles, extent::Extent)
        nt2extent.(geotiles.extent);
        ind = isnothing.(Extents.intersect.(Ref(extent), nt2extent.(geotiles.extent)));
        geotiles = deleteat!(geotiles, ind);
        return geotiles
end



"""
    region_extent(region)

returns extent of region and epsg of coordinates
"""
function region_extent(region)
    if region == :WNA
        extent = Extent(X=(-131, -109), Y=(42, 55))
        epsg = 4326;
    elseif region == :World
        extent = Extent(X=(-180, +180), Y=(-90, 90))
        epsg = 4326;
    else
        error("unrecognized region")
    end
    
    return extent, epsg
end


"""
    allfiles(rootdir; subfolders=false, fn_startswith=nothing, fn_endswith=nothing, 
        fn_contains=nothing, and_or=&, topdown=true, follow_symlinks=false, onerror=throw)

Returns a vector list of file paths that meet the user defined criteria.\n 
# Arguments
- rootdir: directory from which to begin search
- subfolders [false]: include or exclude subfolders in file search
- fn_startswith: include files that startswith
- fn_endswith: include files that startswith
- fn_contains: include files that contain
- and_or [`&`]: only include if all criteria are met (`&`) or include if any criteria 
  are met (`|`)
- topdown, follow_symlinks, onerror: see walkdir documentation
"""
function allfiles(
    rootdir;
    subfolders=false,
    fn_startswith=nothing,
    fn_endswith=nothing,
    fn_contains=nothing,
    and_or=&,
    topdown=true,
    follow_symlinks=false,
    onerror=throw
)

    filelist = String[]
    for (root, _, files) in walkdir(
        rootdir, 
        topdown=topdown, 
        follow_symlinks=follow_symlinks, 
        onerror=onerror
        )

        for file in files
            endswith_tf = true
            startswith_tf = true
            containsin_tf = true

            if !isnothing(fn_endswith)
                endswith_tf = any(endswith.(file, fn_endswith))
            end

            if !isnothing(fn_startswith)
                startswith_tf = any(startswith.(file, fn_startswith))
            end

            if !isnothing(fn_contains)
                containsin_tf = any(contains.(file, fn_contains))
            end

            tf = and_or(and_or(endswith_tf, startswith_tf), containsin_tf)

            if tf
                push!(filelist, joinpath(root, file))
            end
        end

        if subfolders == false
            return filelist
        end
    end
    return filelist
end


"""
    LsqSetup(model, p, p_min, p_max, autodiff)
Model and parameters to pass to LsqFit
"""
struct LsqSetup
    model::Function         # function used to fit data model(parameters, x)
    p::Vector{Float64}      # starting parameters
    p_min::Vector{Float64}  # minimum parameters
    p_max::Vector{Float64}  # maximum parameters
    autodiff::Symbol        # Eccentricity 
end

"""
    geotile_offset(products, dems, geotiles)
Compute vertical and horizontal offsets between altimetry elevations and DEM
"""
function geotile_offset(
    products,
    dems, 
    geotile, 
    paths; 
    interations = 5,
    iter_thresh = 5,
    use_href = false
    )

    # this if for slope rotation
    interp = Linear() # Constant(), Linear(), Cubic(), Quadratic()

    ## check for offset with 
    for product in products
        for dem in dems
            deminfo, _ =  dem_info(dem)
            
            if !(deminfo.epsg == 4326)
                error("dem must be in geographic (lat, lon) projections")
            end

            fn_h = joinpath(paths[product.mission].geotile, "$(geotile.id).arrow")
            fn_dem = joinpath(paths[product.mission].geotile, "$(geotile.id).$dem")
            fn_masks = joinpath(paths[product.mission].geotile, "$(geotile.id).masks")

            if isfile(fn_h) && isfile(fn_dem) && isfile(fn_masks)
                println("$(dem) $(geotile.id)")

                df_h = DataFrame(Arrow.Table(fn_h))
                df_dem = DataFrame(Arrow.Table(fn_dem))
                df_masks = DataFrame(Arrow.Table(fn_masks))

                if use_href
                    dh = vcat(df_h.height...) .- vcat(df_h.height_reference...)
                else
                    dh = vcat(df_h.height...) .- vcat(df_dem.height...)
                end

                dhdlon = vcat(df_dem.dhdx...)
                dhdlat = vcat(df_dem.dhdy...)
                land = vcat(df_masks.land...)

                valid = .!(isnan.(dhdlon) .| isnan.(dhdlat) .| isnan.(dh) .|
                           (abs.(dh) .> 50) .| (.!land) .| (dhdlon .> 1) .| (dhdlon .> 1))

                if sum(valid) < 1000
                    return (NaN, NaN, NaN)
                end

                dhdlon = dhdlon[valid]
                dhdlat = dhdlat[valid]
                lat = vcat(df_h.latitude...)[valid]
                lon = vcat(df_h.longitude...)[valid]
                dh = dh[valid]

                dx, dy, dz, epsg = 
                    offset_geo_fast(lon, lat, dh, dhdlon, dhdlat;
                        interations=interations, iter_thresh=iter_thresh)

                return (dx, dy, dz, epsg)
            else
                return (NaN, NaN, NaN)
            end
        end
    end
end

"""
    offset_geo_fast(lon, lat, dh, dhdlon, dhdlat; interations=3, iter_thresh=7)

ruturn the x, y, z offsets in meters in local UTM coodinates and the epsg for the utm zone. 
Takes longitude [deg], latitude [deg], dh (elevation anomaly m), slope in lon direction 
(dhdlon, in m per degree), and slope in lat direction (dhdlat, in m per degree). Same as 
offset_geo but uses binning to of values to reduce computations
"""
function offset_geo_fast(lon, lat, dh, dhdlon, dhdlat; interations=3, iter_thresh=7)

    # find local utm epsg and tranform from WGS84
    epsg = utm_epsg(mean(lon), mean(lat); always_xy=true)
    epsg = GFT.EPSG(epsg)
    trans = FastGeoProjections.Transformation(GFT.EPSG("EPSG:4326"), epsg, always_xy=true)

    # estimate on a grid then interpolate, working with the raw point cloud 
    # can be prohibitive 

    # change in x [m] per degree lon
    dlat = 0.05
    lat0 = ((floor(minimum(lat) / dlat))*dlat-2*dlat):dlat:((ceil(maximum(lat) / dlat))*dlat+2*dlat)
    dlon = 0.05
    lon0 = ((floor(minimum(lon) / dlon))*dlon-2*dlon):dlon:((ceil(maximum(lon) / dlon))*dlon+2*dlon)

    dx2lon = fill(NaN, length(lat0), length(lon0))
    dx2lat = fill(NaN, length(lat0), length(lon0))

    dy2lon = fill(NaN, length(lat0), length(lon0))
    dy2lat = fill(NaN, length(lat0), length(lon0))

    for i in eachindex(lat0)
        for j in eachindex(lon0)
            cxy = trans(lon0[j], lat0[i])

            pt = inv(trans).([cxy[1] - 0.5, cxy[1] + 0.5], [cxy[2], cxy[2]])
            dx2lon[i, j] = pt[2][1] - pt[1][1]
            dx2lat[i, j] = pt[2][2] - pt[1][2]

            pt = inv(trans).([cxy[1], cxy[1]], [cxy[2] - 0.5, cxy[2] + 0.5])
            dy2lon[i, j] = pt[2][1] - pt[1][1]
            dy2lat[i, j] = pt[2][2] - pt[1][2]
        end
    end
    interp = Linear()

    itp = scale(interpolate(dx2lon, BSpline(interp)), lat0, lon0)
    dx2lon = itp.(lat, lon)

    itp = scale(interpolate(dx2lat, BSpline(interp)), lat0, lon0)
    dx2lat = itp.(lat, lon)

    itp = scale(interpolate(dy2lon, BSpline(interp)), lat0, lon0)
    dy2lon = itp.(lat, lon)

    itp = scale(interpolate(dy2lat, BSpline(interp)), lat0, lon0)
    dy2lat = itp.(lat, lon)

    # project dhdx and dhdy into local UTM coordinate system

    dhdx, dhdy = (
        (dhdlon .* dx2lon) .+ (dhdlat .* dx2lat),
        (dhdlon .* dy2lon) .+ (dhdlat .* dy2lat)
    )

    df = DataFrame()
    df[!, :dhdx] = vec(dhdx)
    df[!, :dhdy] = vec(dhdy)
    df[!, :dh] = vec(dh)

    # take the median of dh for a range of dhdx, dhdy bins
    df = binstats(df, [:dhdx, :dhdy], [-1:0.001:1, -1:0.001:1], [:dh], col_function=[median])

    # compute h offset to h_reference
    dx, dy, dz, cnt, mad_offset, mad_ref = track_offset(
        BinStatistics.bincenter.(df[:, :dhdx]),
        BinStatistics.bincenter.(df[:, :dhdy]),
        df.dh_median;
        interations=interations,
        iter_thresh=iter_thresh,
        weights=df.nrow
    )

    return dx, dy, dz, cnt, mad_offset, mad_ref, epsg.val
end

"""
    offset_geo(lon, lat, dh, dhdlon, dhdlat; interations=3, iter_thresh=7)

ruturn the x, y, z offsets in meters in local UTM coodinates. Also returns the valid data
count, mad_offset, mad_ref and the epsg for the utm zone and the change in dh (`ddh`). 
Takes longitude [deg], latitude [deg], dh (elevation anomaly m), slope in lon direction 
(dhdlon, in m per degree), and slope in lat direction (dhdlat, in m per degree)
"""
function offset_geo(lon, lat, dh, dhdlon, dhdlat; interations=3, iter_thresh=7, verbose=true)
    # find local utm epsg and tranform from WGS84
    epsg = utm_epsg(mean(lon), mean(lat); always_xy=true)
    epsg = GFT.EPSG(epsg)
    trans = FastGeoProjections.Transformation(GFT.EPSG("EPSG:4326"), epsg, always_xy=true)

    # estimate on a grid then interpolate, working with the raw point cloud 
    # can be prohibitive 

    # change in x [m] per degree lon

    # full solution without using interpolation
    dhdx = similar(dhdlon)
    dhdy = similar(dhdlat)

    for i = eachindex(lon)
        cxy = trans(lon[i], lat[i])
        pt_dx = inv(trans)([cxy[1] - 0.5, cxy[1] + 0.5], [cxy[2], cxy[2]])
        pt_dy = inv(trans)([cxy[1], cxy[1]], [cxy[2] - 0.5, cxy[2] + 0.5])

        dhdx[i], dhdy[i] = (
            (dhdlon[i] * (pt_dx[2][1] - pt_dx[1][1])) + (dhdlat[i] * (pt_dx[2][2] - pt_dx[1][2])),
            (dhdlon[i] * (pt_dy[2][1] - pt_dy[1][1])) + (dhdlat[i] * (pt_dy[2][2] - pt_dy[1][2])))
    end

    # compute h offset to h_reference
    dx, dy, dz, cnt, mad_offset, mad_ref = track_offset(
        dhdx, dhdy, dh;
        interations=interations,
        iter_thresh=iter_thresh, verbose=verbose
    )

    return dx, dy, dz, cnt, mad_offset, mad_ref, epsg.val
end


"""
    geotile_addxy_offsets(geotiles, geotile_dir, xoffset, yoffset)

add x and y offsets in m local UTM to geotile longitude and latitude
"""
function geotile_addxy_offsets(geotiles, geotile_dir, xoffset, yoffset)
    printstyled("add x/y offsets Hugonnet geotiles\n"; color=:blue, bold=true)
    for geotile in eachrow(geotiles)
        outfile = joinpath(geotile_dir, geotile[:id] * ".arrow")
        if isfile(outfile)
            df = DataFrame(Arrow.Table(outfile))
            df.latitude = copy(df.latitude)
            df.longitude = copy(df.longitude)
            epsg = utm_epsg(mean([geotile.extent.min_x, geotile.extent.max_x]), mean([geotile.extent.min_y, geotile.extent.max_y]); always_xy=true)
            trans = Proj.Transformation("EPSG:4326", epsg, always_xy=true)

            for row in eachrow(df)
                xy = trans.(row.longitude,row.latitude)
                x = getindex.(xy,1) .+ xoffset
                y = getindex.(xy,2) .+ yoffset
                lonlat = inv(trans).(x,y)
                row.longitude = getindex.(lonlat,1)
                row.latitude = getindex.(lonlat,2)
            end

            tmp = tempname(dirname(outfile))
            Arrow.write(tmp, df::DataFrame) # do not use compression... greatly slows read time
            mv(tmp, outfile; force=true)

            printstyled("    -> $(geotile.id) x/y offsets added\n"; color=:light_black)
        else
            printstyled("    -> $(geotile.id) does not exist, skipping\n"; color=:light_red)
        end
    end
end


"""
    geotile_coregister(geotile, geotile_dir, dem)
coregister geotile rows to DEM
"""
function geotile_coregister(geotile, geotile_dir, dem)
    t1 = time()

    geotilefn = joinpath(geotile_dir, geotile[:id] * ".arrow");
    outfile = replace(geotilefn, "arrow"=>"offset2$(dem)")
    df = DataFrame(Arrow.Table(geotilefn));
    df_mask = DataFrame(Arrow.Table(replace(geotilefn, "arrow"=>"masks")));
    df_dem = DataFrame(Arrow.Table(replace(geotilefn, "arrow" => dem)));
    df.height_reference = df_dem.height;

    df[!,:land] = df_mask[:,:land];
    df[!,:dhdlon] = df_dem.dhdx;
    df[!,:dhdlat] = df_dem.dhdy;

    df_offset = DataFrame();
    df_offset[!, :dx] = fill(0.0, nrow(df));
    df_offset[!, :dy] = fill(0.0, nrow(df));
    df_offset[!, :dz] = fill(0.0, nrow(df));
    df_offset[!, :count] = fill(0.0, nrow(df))
    df_offset[!, :mad_offset] = fill(0.0, nrow(df))
    df_offset[!, :mad_ref] = fill(0.0, nrow(df))
    df_offset[!, :epsg] = fill(Int64(0), nrow(df));
    df_offset[!, :ddh] = copy(df.height)

    Threads.@threads for i in eachindex(eachrow(df))
        dh = df[i, :height] .- df[i, :height_reference]
        valid = df[i, :land] .& .!isnan.(dh) .& .!isnan.(df[i, :dhdlon])

        if sum(valid) > 100
            df_offset[i, :dx], df_offset[i, :dy], df_offset[i, :dz], df_offset[i, :count], 
            df_offset[i, :mad_offset], df_offset[i, :mad_ref], df_offset[i, :epsg] = 
            offset_geo(df[i, :longitude][valid], df[i, :latitude][valid], dh[valid],
             df[i, :dhdlon][valid], df[i, :dhdlat][valid]; interations=3, iter_thresh=7, verbose=false)
        end

        if df_offset[i, :epsg] == 0 
            df_offset[i, :ddh] .= 0
        else
            df_offset[i, :ddh] = offset_geo(df[i, :longitude], df[i, :latitude], df[i, :dhdlon], df[i, :dhdlat], df_offset[i, :dx], df_offset[i, :dy], df_offset[i, :dz], EPSG(df_offset[i, :epsg]))
        end
    end

    tmp = tempname(dirname(outfile))
    Arrow.write(tmp, df_offset::DataFrame)
    mv(tmp, outfile; force=true)

    total_time = round((time() - t1) / 60, digits=2)
    printstyled("    -> $(geotile.id) $dem offsets calculated: $(total_time) min \n"; color=:light_black)
end


"""
    offset_geo(lon, lat, dh, dhdlon, dhdlat; interations=3, iter_thresh=7)

ruturn the x, y, z offsets in meters in local UTM coodinates. Also returns the valid data
count, mad_offset, mad_ref and the epsg for the utm zone and the change in dh (`ddh`). 
Takes longitude [deg], latitude [deg], dh (elevation anomaly m), slope in lon direction 
(dhdlon, in m per degree), and slope in lat direction (dhdlat, in m per degree)
"""
function offset_geo(lon, lat, dhdlon, dhdlat, dx, dy, dz, epsg)
    # find local utm epsg and tranform from WGS84
    epsg = GFT.EPSG(epsg)
    trans = FastGeoProjections.Transformation(GFT.EPSG("EPSG:4326"), epsg, always_xy=true)

    # estimate on a grid then interpolate, working with the raw point cloud 
    # can be prohibitive 

    # change in x [m] per degree lon

    # full solution without using interpolation
    dhdx = similar(dhdlon)
    dhdy = similar(dhdlat)

    for i = eachindex(lon)
        cxy = trans(lon[i], lat[i])
        pt_dx = inv(trans)([cxy[1] - 0.5, cxy[1] + 0.5], [cxy[2], cxy[2]])
        pt_dy = inv(trans)([cxy[1], cxy[1]], [cxy[2] - 0.5, cxy[2] + 0.5])

        dhdx[i], dhdy[i] = (
            (dhdlon[i] * (pt_dx[2][1] - pt_dx[1][1])) + (dhdlat[i] * (pt_dx[2][2] - pt_dx[1][2])),
            (dhdlon[i] * (pt_dy[2][1] - pt_dy[1][1])) + (dhdlat[i] * (pt_dy[2][2] - pt_dy[1][2])))
    end

    # compute change in h
    dh_offset = track_offset_dh(dhdx, dhdy, dx, dy, dz)
    return dh_offset
end