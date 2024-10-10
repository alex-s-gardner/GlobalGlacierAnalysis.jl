# utilities for working building geogrid point database
include("local_paths.jl")
global pathlocal = setpaths()
const world = Extent(X=(-180, 180), Y=(-90, 90))

struct EpsgPoints
    x::Vector{<:Real}
    y::Vector{<:Real}
    epsg::EPSG
end


"""
    setpaths(geotile_width, mission, product, version)

populate named tuple of paths
"""
function setpaths(geotile_width, mission, product, version)
    geotile_dir = @sprintf "%.0fdeg" geotile_width
    data_dir = joinpath(pathlocal.data_dir, lowercase(string(mission)),
        string(product), lpad(string(version), 3, '0'))
    raw_data_dir = joinpath(data_dir, "raw")
    geotile_dir = joinpath(data_dir, "geotile", geotile_dir)
    granules_remote = joinpath(data_dir, "geotile", geotile_dir, "granules.remote")
    granules_local = joinpath(data_dir, "geotile", geotile_dir, "granules.local")

    altim_paths = (
        raw_data=raw_data_dir,
        geotile=geotile_dir,
        granules_remote=granules_remote,
        granules_local=granules_local
    )
    return altim_paths
end

"""
    geotile_extent(lat, lon, width)

Return extents of geotile
"""
function geotile_extent(lat, lon, width)
    (min_x=(lon - width / 2), min_y=(lat - width / 2),
        max_x=(lon + width / 2), max_y=(lat + width / 2))
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
searchdir(path, key1, key2) = filter(x -> (occursin(key1, x) .& occursin(key2, x)), readdir(path))
searchdir(path, key1) = filter(x -> occursin(key1, x), readdir(path))

"""
    points_plus(granule::ICESat2_Granule{}; bbox = (min_x = -Inf, min_y = -Inf, max_x = Inf, max_y = Inf))

returns the ICESat2 granual *WITH* granual infomation for each track
"""
function points_plus(
    granule;
    extent::Extent=world
)
    try
        p = SpaceLiDAR.points(granule, bbox=extent)
        for i = eachindex(p)
            p[i] = merge(p[i], (; granule_info=Fill(granule, length(p[i].longitude))))
        end

        return p
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
        p = SpaceLiDAR.points(granule, bbox=extent)
        p = merge(p, (; granule_info=Fill(granule, length(p.longitude))))
        return p
    catch ex
        println("-------- error thrown when trying to read: --------")
        println("$(granule.url)")
        println("pease delete file and re-run `geotile_download_granules`")
        println("------------------------------------------------------")
        throw(error("issue reading $(granule.url)"))
    end
end

function rm_corrupt_h5(folder; minsize=Inf)
    # minsize = 1E7 # file size is in bytes
    files = filter!(readdir(folder)) do fname
        fname[end] == '5' .&& filesize(joinpath(folder, fname)) < minsize
    end

    # this checks each file before deleting
    for fname in files
        if fname[end] == '5'
            try
                h5open(joinpath(folder, fname), "r")
            catch e
                if !any(names(e) .== "msg")
                    error("check that HDF5 library has been imported (e.g. using HDF5)")
                else
                    println(e.msg)
                    println("!!!! deleting file !!!!")
                    rm(joinpath(folder, fname))
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
    dt = geotile_width / 2
    lat_center = (world.Y[1]+dt):geotile_width:(world.Y[2]-dt)
    lon_center = (world.X[1]+dt):geotile_width:(world.X[2]-dt)

    extent = vec([geotile_extent(lat, lon, geotile_width) for lon in (lon_center), lat in (lat_center)])
    id = geotile_id.(extent)

    DataFrame(id=id, extent=extent)
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
    rebuild_dataframe=false
)

    # ensure SpaceLiDAR capitilization ot mission 
    mission = mission2spacelidar(mission)

    pathfile = splitdir(outgranulefile)
    if !isdir(pathfile[1])
        error("$(pathfile[1]) does not exist")
    end

    printstyled("identifying grannules that intersect each geotile\n"; color=:blue, bold=true)

    # get intersecting granules [run threads]
    n = size(geotiles, 1)
    if mission == :ICESat2
        granules = Vector{Vector{ICESat2_Granule{product}}}(undef, n)
    elseif mission == :ICESat
        granules = Vector{Vector{ICESat_Granule{product}}}(undef, n)
    elseif mission == :GEDI
        granules = Vector{Vector{GEDI_Granule{product}}}(undef, n)
    else
        error("mission and/or product not recognized")
    end

    Threads.@threads for i in 1:n
        #for i in 1:n
        printstyled("    -> finding granules in geotile $i of $n\n"; color=:light_black)
        granules[i] = search(mission, product; bbox=geotiles[i, :extent], version=version)
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
    intersectindices(a,b; bool = false)

Returns the index mapping between intersecting elements in `a` and `b`
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

Merge DataFrames with identical columns based on id_unique
- retains df_left rows if not contianed in df_right
- retains df_right rows if not contianed in df_left
- df_right rows are retained over df_left if rows have the same id_unique
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

Merge DataFrames with identical columns based on id_unique
- retains df_left rows if not contianed in df_right
- retains df_right rows if not contianed in df_left
- df_left rows are retained over df_rigth if rows have the same id_unique
"""
function leftmerge(df_left::DataFrame, df_right::DataFrame, id_unique::Symbol)
    _, iright = intersectindices(df_left[:, id_unique], df_right[:, id_unique], bool=true)
    if any(.!iright)
        df_left = vcat(df_left, df_right[.!iright, :]) # add new rows
    end
    return df_left
end

function geotile_build(geotile_granules, geotile_dir; warnings=true, fmt=:arrow, replace_corrupt_h5 = true)
    printstyled("building geotiles\n"; color=:blue, bold=true)

    # remove empty granules
    geotile_granules = geotile_granules[.!isempty.(geotile_granules.granules), :]

    if !warnings
        Logging.disable_logging(Logging.Warn) # or e.g. Logging.Info
    end

    #asyncmap(eachrow(geotile_granules); ntasks = 10) do row
    for row in eachrow(geotile_granules)

        # tiles
        outfile = joinpath(geotile_dir, row.id * ".$fmt")

        if isempty(row.granules)
            # do nothing
            printstyled("    -> $(row.id): has no granules\n"; color=:light_red)
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
                printstyled("    -> $(row[:id]): generation complete [read: $read_time min, write: $write_time min]\n"; color=:light_black)
            else
                printstyled("    -> $(row[:id]): no new granules to add to exisitng GeoTile\n"; color=:light_green)
            end
        end
    end
end

function getpoints(granule; extent=nothing, replace_corrupt_h5 = false)

    try
        pts = SpaceLiDAR.points(granule; bbox=extent)
        return pts
    catch e
        if e isa InterruptException
            println("function terminated by user")
            rethrow(e)
        else
            printstyled("issue reading: $(granule.url)\n"; color=:red)
            if replace_corrupt_h5 
                # if a download was killed mid-stream it can lead to corrupt files
                printstyled("               deleting old file and re-downloading from source\n"; color=:red)
                if isfile(granule.url)
                    rm(granule.url)
                end
                (savedir, _) = splitdir(granule.url);
                try
                    granule = search(SpaceLiDAR.mission(granule), granule.info.type; version=granule.info.version, id=granule.id)[1]
                    download!(granule, savedir)
                    pts = SpaceLiDAR.points(granule; bbox=extent)
                    return pts
                catch e
                    throw(e)
                end
            else
                throw(e)
            end
        end
    end
end


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
    geotile_id(extent)

Returns the geotile id given a geotile `extent`
"""
function geotile_id(extent)
    id = @sprintf("lat[%+03.0f%+03.0f]lon[%+04.0f%+04.0f]", extent.min_y, extent.max_y, extent.min_x, extent.max_x)
    return id
end

"""
    geotile_extent(geotile_id)

Returns the geotile extents given a geotile `id`
"""
function geotile_extent(geotile_id)

    min_y = parse(Int, geotile_id[5:7])
    max_y = parse(Int, geotile_id[8:10])

    min_x = parse(Int, geotile_id[16:19])
    max_x = parse(Int, geotile_id[20:23])

    return nt = (min_x=min_x, min_y=min_y, max_x=max_x, max_y=max_y)
end


"""
    geotile_utm!(df::DataFrame{}; height = nothing)

add x and y coodinates for local utm or polar stereo zone to an geotile DataFrame
"""
function geotile_utm!(df::DataFrame; height=nothing)
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
        if !bbox_overlap(bbox(ga), extent)

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
                    itp = Interpolations.extrapolate(scale(Interpolations.interpolate(ga1[j].A[:, :, i], BSpline(interp)), x0, y0 .* y2x_scale), nodatavalue)
                    deriv[j][:, i] = itp.(x, y .* y2x_scale)
                end
            end
        end
        val = vcat([val], deriv) #0.002s
    end
    return val
end


"""
    epsg2epsg_nodata(x, y, height, from_epsg, to_epsg, nodatavalue)

A wrapper for `epsg2epsg` to prevent modification of `height` no data values
"""
function epsg2epsg_nodata(x, y, height, from_epsg, to_epsg, nodatavalue)
    valid = (height .!= nodatavalue) .& (.!isnan.(height))
    if any(valid)
        #(x[valid], y[valid], height[valid]) = epsg2epsg(x[valid], y[valid], height[valid], from_epsg, to_epsg; parse_output = false)
        x[valid], y[valid], height[valid] = epsg2epsg(x[valid], y[valid], height[valid], from_epsg, to_epsg; parse_output=true)
    end
    return x, y, height
end


"""
    epsg2epsg(x, y, height, from_epsg, to_epsg)

Returns `x`, `y`,`height` in `to_epsg` projection
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
    epsg2epsg(x, y, from_epsg, to_epsg)

Returns `x`, `y` in `to_epsg` projection
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
    decimalyear2datetime(decyear)

converts decimalyear to DateTime to  [e.g. decimalyear2datetime(2018.813698630137) = DateTime("2018-10-24")]
"""
function decimalyear2datetime(decyear)

    # check if date is a leapyear
    yr0 = floor(decyear);
    yr = DateTime(yr0)

    if isleapyear(yr)
        number_of_days = 366
    else
        number_of_days = 365
    end

    # integer date
    d = number_of_days * (decyear .- yr0);
    d0 = floor(d)

    # integer hour
    h = (d-d0) * 24
    h0 = floor(h)

    # integer minute
    m = (h-h0)*60
    m0 = floor(m)

    # integer second
    s = (m - m0) * 60
    s0 = floor(s)

    # integer millisecond
    ms = round((s - s0)*100)

    # calculate datetime
    datetime = yr + Day(d0) + Hour(h0) + Minute(m0) + Second(s0) + Millisecond(ms)

    return datetime
end



"""
    regular_grid(center_extent, node_spacing; node_width = node_spacing; node_shape = "square")

Returns a named tuple of node centers (`x_node_center`, `y_node_center`) and `mode_half_width`. 
`node_shape` is also provided
"""
function regular_grid(center_extent, node_spacing; node_width=node_spacing, node_shape="square")
    x_node_center = center_extent.x_min:node_spacing:center_extent.x_max
    y_node_center = center_extent.y_min:node_spacing:center_extent.y_max
    node_half_width = node_width / 2

    node_center = ndgrid(y_node_center, x_node_center)
    (; node_center, x_node_center, y_node_center, node_half_width, node_shape)
end


"""
    regular_grid_extents(x_center, y_center, half_width)

Returns the extents of a box with  width 2 x `half_width` centered at `x_center`, `y_center`. 
"""
function regular_grid_extents(x_center, y_center, half_width)
    (x_min=x_center - half_width, y_min=y_center - half_width,
        x_max=x_center + half_width, y_max=y_center + half_width)
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
    method::F=mean,
) where {F<:Function}

    # find bin breakpoints
    p = sortperm(vcat(x, xbin_edges))
    bins = findall(>(length(x)), p)

    # initialize outputs
    bin_count = Int.(diff(bins) .- 1)
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
    method::F=median
) where {F<:Function}

    # find bin breakpoints
    py = sortperm(vcat(y, ybin_edges))
    binsy = findall(>(length(y)), py)
    biny_count = Int.(diff(binsy) .- 1)

    x_binned = zeros(eltype(x), length(ybin_edges) - 1, length(xbin_edges) - 1)
    y_binned = zeros(eltype(y), length(ybin_edges) - 1, length(xbin_edges) - 1)
    z_binned = zeros(eltype(z), length(ybin_edges) - 1, length(xbin_edges) - 1)
    bin_count = zeros(Int32, length(ybin_edges) - 1, length(xbin_edges) - 1)

    # loop for each x bin
    Threads.@threads for i = findall(biny_count .> 0)
        y0 = @view y[py[binsy[i]+1:binsy[i+1]-1]]
        x0 = @view x[py[binsy[i]+1:binsy[i+1]-1]]
        z0 = @view z[py[binsy[i]+1:binsy[i+1]-1]]

        # find bin breakpoints
        px = sortperm(vcat(x0, xbin_edges))
        binsx = findall(>(length(x0)), px)

        # initialize outputs
        bin_count[i, :] .= Int.(diff(binsx) .- 1)

        # calculate binned metrics
        for j = findall(bin_count[i, :] .> 0)
            x_binned[i, j] = method(@view x0[px[binsx[j]+1:binsx[j+1]-1]])
            y_binned[i, j] = method(@view y0[px[binsx[j]+1:binsx[j+1]-1]])
            z_binned[i, j] = method(@view z0[px[binsx[j]+1:binsx[j+1]-1]])
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
    method::F=madnorm,
    threshold::Number=10,
) where {F<:Function}

    # find bin breakpoints
    p = sortperm(vcat(x, xbin_edges))
    bins = findall(>(length(x)), p)

    # initialize outputs
    outlier = falses(size(y))
    bin_count = Int.(diff(bins) .- 1)

    # identify outliers
    for i = findall(bin_count .> 0)
        outlier[p[bins[i]+1:bins[i+1]-1]] = method(@view(y[p[bins[i]+1:bins[i+1]-1]])) .> threshold
    end

    return outlier
end


"""
    madnorm(x)

Returns `x_madnorm` the normalized median absolute deviation of `y` from the global median. 
Note: consistent_estimator applied
"""
function madnorm(x)
    consistent_estimator = 1.4826 #mad to sigm
    x_abs = abs.(x .- median(x))
    x_madnorm = x_abs ./ (median(x_abs) .* consistent_estimator)
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
        return epsg
    elseif lat < -80
        # Antarctic Polar Stereographic
        epsg = 3031
        epsg = "EPSG:$epsg"
        return epsg
    end

    # make sure lon is from -180 to 180
    lon = lon - floor((lon + 180) / (360)) * 360

    # int versions
    ilat = floor(Int64, lat)
    ilon = floor(Int64, lon)

    # get the latitude band
    band = max(-10, min(9, fld((ilat + 80), 8) - 10))

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
            printstyled("local copy of $geoid file not found, downloading from: $url \n"; color=:blue, bold=true)
            HTTP.download(url, path2goid)
        end
    end

    # not sure if this Type decleration helps at all, feel free to delete
    return GeoArrays.read(path2goid; masked=false)
end

"""
    gaussian(σ::Real, [l]) -> g

Construct a 1d gaussian kernel `g` with standard deviation `σ`, optionally
providing the kernel length `l`. The default is to extend by two `σ`
in each direction from the center. `l` must be odd. 
**Stolen from ImageFiltering.jl**
"""
function gaussian(σ::Real, l::Int=4 * ceil(Int, σ) + 1)
    isodd(l) || throw(ArgumentError("length must be odd"))
    w = l >> 1
    g = σ == 0 ? [exp(0 / (2 * oftype(σ, 1)^2))] : [exp(-x^2 / (2 * σ^2)) for x = -w:w]
    g / sum(g)
end


"""
    gaussian((σ1, σ2, ...), [l]) -> (g1, g2, ...)

Construct a multidimensional gaussian filter as a product of single-dimension
factors, with standard deviation `σd` along dimension `d`. Optionally
provide the kernel length `l`, which must be a tuple of the same
length.
**Stolen from ImageFiltering.jl**
"""
gaussian(σs::NTuple{N,Real}, ls::NTuple{N,Integer}) where {N} =
    map(gaussian, σs, ls)
gaussian(σs::NTuple{N,Real}) where {N} = map(gaussian, σs)
gaussian(σs::AbstractVector, ls::AbstractVector) = gaussian((σs...,), (ls...,))
gaussian(σs::AbstractVector) = gaussian((σs...,))



_dx(v, gsx) = (-v[1] + v[3] - 2 * v[4] + 2 * v[6] - v[7] + v[9]) / (8 * gsx) # for GeoArrays x is rows
_dy(v, gsy) = (-v[1] - 2 * v[2] - v[3] + v[7] + 2 * v[8] + v[9]) / (8 * gsy) # for GeoArrays y is columns
_ddx(v, gsx) = ((v[4] + v[6]) / 2 - v[5]) / (gsx .^ 2) # for GeoArrays x is rows
_ddy(v, gsy) = ((v[2] + v[8]) / 2 - v[5]) / (gsy .^ 2) # for GeoArrays y is columns

""" 
    dem_height(lon, lat, dem; filter_halfwidth::Union{Number,Nothing}= nothing)
Returns dem `height` for `lon`, `lat` locations. If filter_halfwidth is supplied, a 2D Gaussain
filter with sigma = `filter_halfwidth` is applied to the `dem` prior to sampeling
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
        interp = Interpolations.Cubic(Line(OnGrid())); # Constant(), Linear(), Cubic(Line(OnGrid())), Quadratic()
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
    geotile_extract_dem(geotile, geotile_dir, dem; filter_halfwidth, filter_kernel, slope, job_id, force_remake)
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
        printstyled("    ->$job_id $geotile_id $dem already exists, skipping\n"; color=:green)
    elseif isfile(path2geotile)
        t1 = time()

        begin #0.002s
        # load dem info
        deminfo, _ = dem_info(dem)
        
        # check if geotile extents intersect DEM extents 
        if !(Extents.intersects(nt2extent(deminfo.extent), GeoTiles.extent(geotile_id)))
            printstyled("    ->$job_id $geotile_id $dem outside of dem limits, skipping\n"; color=:light_red)
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
            printstyled("    ->$job_id $geotile_id had an empty DataFrame, skipping\n"; color=:grey, bold=true)
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
            printstyled("    ->$job_id $geotile_id is all empty, skipping\n"; color=:light_red)
            return
        end

        # check if within extents
        if (deminfo.extent.min_x != -180.0) || (deminfo.extent.min_y != -90.0) || (deminfo.extent.max_x != 180.0) || (deminfo.extent.max_y != 90.0)
            ind = within.(Ref(deminfo.extent), lon, lat)
            if !any(ind) 
                printstyled("    ->$job_id $geotile_id $dem outside of dem limits, skipping\n"; color=:light_red)
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
            printstyled("    ->$job_id $geotile_id $dem extracted: $(total_time) min \n"; color=:light_black)
        end
    else
        printstyled("    ->$job_id $geotile_id does not exist\n"; color=:yellow)
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
    filter_halfwidth::Union{Number,Nothing}=nothing,
    filter_kernel=:gaussian,
    slope=false,
    curvature=false,
    job_id="",
    xoffset=0, # add offset to altimetry locations in local UTM coordinates [m]
    yoffset=0, # add offset to altimetry locations in local UTM coordinates [m]
    force_remake=false
)

    printstyled("extracting $dem heights for $job_id geotiles\n"; color=:blue, bold=true)

    Threads.@threads for geotile in eachrow(geotiles)
    #for geotile in eachrow(geotiles)
        geotile_extract_dem(geotile.id, geotile_dir, dem; 
            filter_halfwidth, filter_kernel, slope, curvature,
            job_id, xoffset, yoffset, force_remake)
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
    geotile_aggrigate_reduce(geotiles::DataFrame, geotile_dir::String, 
        vars::Union{Vector{String}, String}; extension::String = ".arrow")

Aggrigate variables in `vars` for list of `geotiles`` and reduce to `Vectors` of `elements`.
subset_fraction [0 to 1] allows for random subsampeling of the data  
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
    geotile_aggrigate_reduce(geotiles::DataFrame, geotile_dir::String, 
        vars::Union{Vector{String}, String}; extension::String = ".arrow")

Aggrigate variables in `vars` for list of `geotiles`` and reduce to `Vectors` of `elements`. 
subset_fraction [0 to 1] allows for random subsampeling of the data 
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
    dataframe_reduce(df::DataFrame)

Reduce all variables in `df` to `Vectors` of `elements`.
"""
function dataframe_reduce(df)
    df0 = DataFrame()
    isarray = [isa(d, AbstractArray) for d in df[1, :]]
    refarray = findfirst(isarray)
    refsize = size.(df[:, refarray])

    for var in names(df)
        # in some cases variables are saved as single values if they are constant 
        # for all row values while others are saved as vectors

        # check it it's a vector quantitiy
        if typeof(df[1, var]) <: AbstractVector
            df0[!, var] = reduce(vcat, df[:, var])
        else
            #if not expand to a vector 
            df0[!, var] = reduce(vcat, [fill(df[k, var], refsize[k]) for k = eachindex(df[:, var])])
        end
    end
    return df0
end


"""
    itslive_proj!(df::DataFrame{}; height = nothing)
add x and y coodinates for local itslive projection to the DataFrame `df`
"""
function itslive_proj!(df; height=nothing)
    epsg = itslive_epsg(mean(df.longitude), mean(df.latitude); always_xy=true)

    if isnothing(height)
        df[!, :X], df[!, :Y] = epsg2epsg(df.longitude, df.latitude, "EPSG:4326", epsg, parse_output=true)
    else
        df[!, :X], df[!, :Y], df[!, :H] = epsg2epsg(df.longitude, df.latitude, height, epsg, parse_output=true)
    end

    return df, epsg
end

"""
    itslive_epsg(lon, lat)
Return epsg code for the ITS_LIVE projection
"""
function itslive_epsg(longitude, latitude; always_xy=true)

    if !always_xy
        latitude, longitude = (longitude, latitude)
    end

    if latitude > 55
        # NSIDC Sea Ice Polar Stereographic North
        return epsg = "EPSG:3413"
    elseif latitude < -56
        # Antarctic Polar Stereographic
        return epsg = "EPSG:3031"
    end

    # make sure lon is from -180 to 180
    lon = longitude - floor((longitude + 180) / (360)) * 360

    # int versions
    ilat = floor(Int64, latitude)
    ilon = floor(Int64, lon)

    # get the latitude band
    band = max(-10, min(9, fld((ilat + 80), 8) - 10))

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
    itslive_paramfiles(lon, lat; gridsize = 240, path2param = "/Users/...", always_xy = true)
Return paths to its_live parameter files
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
    itslive_extract(lon, lat, vars; gridsize = 240, path2param = "/Users/...", always_xy = true)
Return point values from its_live parameter files
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

        if !bbox_overlap(bbox(ga), extent)
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
    info, goids_folder  = dem_info()
Return dem information and goids_folder
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
    normalize(data)
Return a nomalized version of data
"""
function normalize(data)
    m = mean(data)
    s = std(data)
    n = (data .- m) ./ s

    datan = (std=s, mean=m, norm=n)

    return datan
end


"""
    centroid(ga::GeoArray)
return cetroid of GeoArray
"""
function centroid(ga::GeoArray)
    bbox = GeoArrays.bbox(ga)
    (x=(bbox.min_x + bbox.max_x) / 2, y=(bbox.min_y + bbox.max_y) / 2)
end


"""
xy_gridsize(ga::GeoArray)
Return x and y grid size in meters
"""
function xy_gridsize(ga::GeoArray)

    # check if geographic
    isgeographic = Base.contains(ga.crs.val, "AXIS[\"Latitude\",NORTH]")
    if isgeographic
        cent = centroid(ga)

        x_lla = LLA(cent.y, cent.x, 0.0)
        y_lla = LLA(cent.y + 1, cent.x, 0.0)
        dist_lat = euclidean_distance(x_lla, y_lla)

        x_lla = LLA(cent.y, cent.x, 0.0)
        y_lla = LLA(cent.y, cent.x + 1, 0.0)
        dist_lon = euclidean_distance(x_lla, y_lla)

        gridsize = (x=ga.f.linear[1, 1] * dist_lon, y=ga.f.linear[2, 2] * dist_lat)
    else
        gridsize = (x=ga.f.linear[1, 1], y=ga.f.linear[2, 2])
    end
    return gridsize
end

"""
    dist_ll2xy(longitude, latitude)
Return local utm x and y distance in meters per degree latitude and longitude
"""
function dist_ll2xy(longitude, latitude)
    delta = 0.0001
    lla = LLA.(latitude, longitude, Ref(0))
    lla_dlat = LLA.(latitude .+ delta, longitude, Ref(0))
    lla_dlon = LLA.(latitude, longitude .+ delta, Ref(0))

    zone, north = Geodesy.utm_zone(mean(latitude), mean(longitude))
    utm_from_lla = UTMfromLLA(zone, north, Geodesy.wgs84)

    lonlat2xy = [Matrix{eltype(latitude[1])}(undef, (2, 2)) for i in latitude]

    Threads.@threads for i in eachindex(lla)
        p = utm_from_lla(lla[i])::UTM{Float64}
        p_dlat = utm_from_lla(lla_dlat[i])::UTM{Float64}
        p_dlon = utm_from_lla(lla_dlon[i])::UTM{Float64}

        lonlat2xy[i] = [(p_dlon.x-p.x) (p_dlat.y-p.y);
            (p_dlat.x-p.x) (p_dlon.y-p.y)] ./ delta
    end
    return lonlat2xy
end

"""
    dh_ll2xy(lon, lat, dhdlat, dhdlon)
Return slope in x and y UTM coodinates per meter from WGS84 lat lon slopes per degree
"""
function dh_ll2xy(lon, lat, dhdlon, dhdlat)
    delta = 0.0001
    lla = LLA.(lat, lon, Ref(0))
    lla_dlat = LLA.(lat .+ delta, lon, Ref(0))
    lla_dlon = LLA.(lat, lon .+ delta, Ref(0))

    zone, north = Geodesy.utm_zone(mean(lat), mean(lon))
    utm_from_lla = UTMfromLLA(zone, north, Geodesy.wgs84)

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

        dhdx[i] = ((dhdlon[i] * a) / sqrt(a .^ 2 + b .^ 2) + (dhdlat[i] * c) / (c + d)) / (a + c)
        dhdy[i] = ((dhdlon[i] * b) / sqrt(a .^ 2 + b .^ 2) + (dhdlat[i] * d) / (c + d)) / (b + d)
    end
    return dhdx, dhdy
end


"""
    dh_ll2aa(x, y, dhdx, dhdy, epsg_pts, epsg_dem)
Return slope in along-track and across-track in DEM coodinates per meter
"""
function dh_xy2aa(x, y, dhdx, dhdy, epsg_pts, epsg_dem)

    if epsg_pts !== epsg_dem
        x, y = epsg2epsg(x, y, "EPSG:4326", "EPSG:$epsg_dem"; parse_output=true)
    end

    dx = x[(begin+1):end] .- x[begin:(end-1)]
    push!(dx, dx[end])

    dy = y[(begin+1):end] .- y[begin:(end-1)]
    push!(dy, dy[end])

    theta = atan.(dx, dy)
    sn = sin.(theta)
    cs = cos.(theta)

    along = dhdx .* sn .+ dhdy .* cs
    across = -(dhdx .* cs .- dhdy .* sn)

    return along, across
end


"""
    dh_ll2aa(lon, lat, dhdlat, dhdlon)
Return slope in along-track and across-track in UTM coodinates per meter from WGS84 lat lon slopes per degree
"""
function dh_ll2aa(lon, lat, dhdlon, dhdlat)

    if length(lon) == 1
        along = zeros(eltype(lon), 1)
        across = zeros(eltype(lon), 1)
        return along, across
    else
        delta = 0.0001
        lla = LLA.(lat, lon, Ref(0))
        lla_dlat = LLA.(lat .+ delta, lon, Ref(0))
        lla_dlon = LLA.(lat, lon .+ delta, Ref(0))

        zone, north = Geodesy.utm_zone(mean(lat), mean(lon))
        utm_from_lla = UTMfromLLA(zone, north, Geodesy.wgs84)

        dhdx = similar(lat)
        dhdy = similar(lon)
        p = Vector{UTM{Float64}}(undef, length(lat))

        Threads.@threads for i in eachindex(lla)
            p[i] = utm_from_lla(lla[i])::UTM{Float64}
            p_dlat = utm_from_lla(lla_dlat[i])::UTM{Float64}
            p_dlon = utm_from_lla(lla_dlon[i])::UTM{Float64}

            a = (p_dlat.x - p[i].x) / delta
            b = (p_dlat.y - p[i].y) / delta

            c = (p_dlon.x - p[i].x) / delta
            d = (p_dlon.y - p[i].y) / delta

            dhdx[i] = ((dhdlon[i] * a) / (a + b) + (dhdlat[i] * c) / (c + d)) / (a + c)
            dhdy[i] = ((dhdlon[i] * b) / (a + b) + (dhdlat[i] * d) / (c + d)) / (b + d)
        end

        dx = [k.x for k in p]
        dx = dx[(begin+1):end] .- dx[begin:(end-1)]
        push!(dx, dx[end])

        dy = [k.y for k in p]
        dy = dy[(begin+1):end] .- dy[begin:(end-1)]
        push!(dy, dy[end])

        theta = atan.(dx, dy)
        sn = sin.(theta)
        cs = cos.(theta)
        across = -(dhdx .* cs .- dhdy .* sn)
        along = dhdx .* sn .+ dhdy .* cs

        return along, across
    end
end

"""
    geotile_extract_masks(geotile, geotile_dir; vars, force_remake)
"""
function geotile_extract_mask(
    geotile_id::String,
    geotile_dir::String;
    vars::Vector{Symbol}=[:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean],
    geotile_ext::String=".arrow",
    job_id="",
    force_remake=false
)

    path2geotile = joinpath(geotile_dir, geotile_id * geotile_ext)
    outfile = replace(path2geotile, ".arrow" => ".masks")

    if isfile(outfile) && !force_remake
        printstyled("    ->$job_id $geotile_id masks already exists, skipping\n"; color=:green)
    elseif isfile(path2geotile)
        t1 = time()

        df0 = DataFrame(Arrow.Table(path2geotile))[!, [:longitude, :latitude]]
        #vlength = length.(df0.longitude)

        #lon = vcat(df0.longitude...)
        #lat = vcat(df0.latitude...)

        if isempty(df0.longitude)
            printstyled("    ->$job_id $geotile_id is all empty, skipping\n"; color=:light_red)
            return
        end

        df = itslive_extract(df0.longitude, df0.latitude, vars, path2param=pathlocal.itslive_parameters)

        #df = DataFrame()
        #for v in vars
        #    df[!, v] = [Array{eltype(df0[!, v])}(undef, l) for l in vlength]
        #end

        # start = 1
        # for row = eachrow(df)
        #     stop = length(row[vars[1]]) + start - 1
        #     for v in vars
        #         row[v] = df0[!, v][start:stop]
        #     end
        #     start = stop + 1
        # end

        tmp = tempname(dirname(outfile))
        Arrow.write(tmp, df::DataFrame)
        mv(tmp, outfile; force=true)

        total_time = round((time() - t1) / 60, digits=2)
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
    vars::Vector{Symbol}=[:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean],
    job_id="",
    force_remake=false
)

    printstyled("extracting masks for $job_id geotiles\n"; color=:blue, bold=true)

    #asyncmap(eachrow(geotiles); ntasks = 4) do geotile
    for geotile in eachrow(geotiles)
        geotile_extract_mask(geotile.id, geotile_dir; vars=vars, job_id=job_id, force_remake=force_remake)
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
    verbose=true
)

    valid = .!isnan.(dh)
    θ = []
    mad_offset = 0.0
    mad_ref = 0.0
    ddh = []

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
    valid_count_thresh=100,
    force_remake=false
)

    path2geotile = joinpath(geotile_dir, geotile_id * ".arrow")
    outfile = replace(path2geotile, ".arrow" => ".$(dem)_offset")

    if isfile(outfile) && !force_remake
        printstyled("    -> $geotile_id $dem track offsets already exists, skipping\n"; color=:green)
    elseif isfile(path2geotile)
        t1 = time()

        # check that files exist
        fn_dem = replace(path2geotile, ".arrow" => ".$dem")
        fn_masks = replace(path2geotile, ".arrow" => ".masks")
        fn_canopyh = replace(path2geotile, ".arrow" => ".canopyh")

        if !isfile(fn_dem)
            printstyled("    -> $dem for $geotile_id does not exist\n"; color=:yellow)
            return
        elseif !isfile(fn_masks)
            printstyled("    -> masks for $geotile_id does not exist\n"; color=:yellow)
            return
        elseif !isfile(fn_canopyh)
            printstyled("    -> canopy for $geotile_id does not exist\n"; color=:yellow)
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
            df0[!, v] = Vector{Float32}(undef, length(df.longitude))
        end

        vars = [:dh_along, :dh_across]
        for v in vars
            df0[!, v] = [Array{Float32}(undef, l) for l in vlength]
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

            w = 1 ./ max.(row.canopyh, 0.01)
            row0.weight = sum(w[valid])

            # solve for optimal vertical, along and across track offsets 
            row0.right, row0.forward, row0.above, _, _, _ =
                track_offset(row0.dh_along[valid], row0.dh_across[valid], row.dh[valid]; weights=w[valid])
        end

        tmp = tempname(dirname(outfile))
        Arrow.write(tmp, df0::DataFrame)
        mv(tmp, outfile; force=true)

        total_time = round((time() - t1) / 60, digits=2)
        printstyled("    -> $geotile_id $dem offsets calculated: $(total_time) min \n"; color=:light_black)
    else
        printstyled("    -> $geotile_id does not exist\n"; color=:yellow)
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
    valid_count_thresh=100,
    force_remake=false
)

    printstyled("calculating track offsets for geotiles\n"; color=:blue, bold=true)

    #asyncmap(eachrow(geotiles); ntasks = 4) do geotile
    for geotile in eachrow(geotiles)
        geotile_track_offset(geotile.id, geotile_dir, dem; valid_count_thresh=valid_count_thresh, force_remake=force_remake)
    end
end

"""
    pointextract(geotile, geotile_dir; var_name, force_remake)
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
    geotile_pointextract(geotiles, geotile_dir; var_name, force_remake, interp)
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

    printstyled("extracting $var_name for geotiles\n"; color=:blue, bold=true)
    geotile_dirs isa AbstractVector || (geotile_dirs = [geotile_dirs])
    job_ids isa AbstractVector || (job_ids = [job_ids])

    #asyncmap(eachrow(geotiles); ntasks = 4) do geotile
    Threads.@threads for geotile in eachrow(geotiles)
    #for geotile in eachrow(geotiles)
        t1 = time()

        # check if ga is in EPSG:4326, if so then it's possible to check if lats and lons fall within ga 
        if contains(ga.crs.val, "\"EPSG\",\"4326\"")

            if !bbox_overlap(bbox(ga), extent2nt(geotile.extent))
                printstyled("    ->$(geotile.id) does not overlap $var_name, skipping\n"; color=:light_red)
                continue
            end
        end

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

                if any((df.longitude .== 0) .& (df.latitude .== 0))
                    @warn "$(path2geotile) contains old longitude = 0, latitude = 0 convention for missing DataFrames"
                end

                stop[i] = nrow(df)
            end
            stop = Int.(stop)

            #try
            df = pointextract(df, ga; var_name, nodatavalue, interp)
            #catch
            #    return df 
            #end
            if isempty(df)
                printstyled("    ->$(geotile.id) no valid $var_name, skipping\n"; color=:light_red)
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
            printstyled("    ->$job_id $(geotile.id) $var_name extracted: $(total_time) min \n"; color=:light_black)
            
        elseif all(.!geotiles_exist)

            job_id = job_ids[.!geotiles_exist]
            printstyled("    ->$job_id $(geotile.id) does not exist,  skipping\n"; color=:light_red)
        else
            job_id = job_ids[geotiles_exist .& outfiles_exist]
            printstyled("    ->$(geotile.id) $var_name already extracted, skipping\n"; color=:light_green)
        end
    end
end

"""
    geotile_load(geotile_width, filesuffix, data_dir; extent::Extent = world)

Returns a single dataframe containing all geotiles [within extent]
"""
function geotile_load(geotile_width, filesuffix, data_dir; extent::Extent=world)

    # define geotiles
    geotiles = geotile_define(geotile_width)

    # geotiles within extent
    foo = nt2extent.(geotiles.extent)
    ind = Extents.intersects.(Ref(extent), foo)

    geotiles = geotiles[ind, :]

    df = DataFrame()

    for geotile in eachrow(geotiles)
        filename = joinpath(data_dir, "$(geotile.id).$(filesuffix)")

        if isfile(filename)
            df = append!(df, DataFrame(Arrow.Table(filename)))
        end
    end
    return df
end

"""
    nt2extent(nt::NamedTuple{(:min_x, :min_y, :max_x, :max_y), NTuple{4, <:Number}})

converts a NamedTuple to an Extent
"""
function nt2extent(nt::NamedTuple{(:min_x, :min_y, :max_x, :max_y)})
    Extent(X=(nt.min_x, nt.max_x), Y=(nt.min_y, nt.max_y))
end

"""
    extent2nt(extent::Extent)

converts a Extent to a GeoArrays NamedTuple
"""
function extent2nt(extent::Extent)
    (min_x=minimum(extent.X), min_y=minimum(extent.Y), max_x=maximum(extent.X), max_y=maximum(extent.Y))
end


"""
    geotile_subset(geotiles, extent::Extent)

returns subset of geotiles that intersect extent
"""
function geotile_subset(geotiles, extent::Extent)
    ind = .!isnothing.(Extents.intersect.(Ref(extent), geotiles.extent))
    geotiles = geotiles[ind, :]
    return geotiles
end


"""
    geotile_subset!(geotiles, extent::Extent)

returns subset of geotiles that intersect extent
"""
function geotile_subset!(geotiles, extent::Extent)
    ind = .!(Extents.intersects.(Ref(extent), geotiles.extent))
    geotiles = deleteat!(geotiles, ind)
    return geotiles
end



"""
    region_extent(region)

returns extent of region and epsg of coordinates
"""
function region_extent(region)
    if string(region)[1:3] == "RGI"
        numid = string(region)[4:5]
        shp = allfiles(pathlocal.RGI_regions; subfolders=true, fn_endswith="o1regions.shp")
        s = DataFrame(Shapefile.Table(shp[1]))
        if numid == "98" # HMA
            ind = (s.o1region .== "13") .| (s.o1region .== "14") .| (s.o1region .== "15")
            ext = GeoInterface.bbox.(s.geometry[ind])
            ext = reduce(Extents.union, ext)
        else
            ind = s.o1region .== numid
            ext = GeoInterface.bbox(s.geometry[ind][1])
        end
        epsg = 4326
    elseif region == :World
        ext = Extent(X=(-180, +180), Y=(-90, 90))
        epsg = 4326
    else
        error("unrecognized region")
    end

    return ext, epsg
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

    if subfolders == false
        files = readdir(rootdir)
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
                push!(filelist, joinpath(rootdir, file))
            end
        end
    else
        # walk dir is very slow so use readdir when you can
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
    interations=5,
    iter_thresh=5,
    use_href=false
)

    ## check for offset with 
    for product in products
        for dem in dems
            deminfo, _ = dem_info(dem)

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
function offset_geo_fast(
    lon,
    lat,
    dh,
    dhdlon, dhdlat; interations=3, iter_thresh=7)

    interp = Linear() # Constant(), Linear(), Cubic(), Quadratic()

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
        df[:, :dhdx],
        df[:, :dhdy],
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
            epsg = utm_epsg(mean(geotile.extent.X), mean(geotile.extent.Y); always_xy=true)
            trans = Proj.Transformation("EPSG:4326", epsg, always_xy=true)

            for row in eachrow(df)
                xy = trans.(row.longitude, row.latitude)
                x = getindex.(xy, 1) .+ xoffset
                y = getindex.(xy, 2) .+ yoffset
                lonlat = inv(trans).(x, y)
                row.longitude = getindex.(lonlat, 1)
                row.latitude = getindex.(lonlat, 2)
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

    geotilefn = joinpath(geotile_dir, geotile[:id] * ".arrow")
    outfile = replace(geotilefn, "arrow" => "offset2$(dem)")
    df = DataFrame(Arrow.Table(geotilefn))
    df_mask = DataFrame(Arrow.Table(replace(geotilefn, "arrow" => "masks")))
    df_dem = DataFrame(Arrow.Table(replace(geotilefn, "arrow" => dem)))
    df.height_reference = df_dem.height

    df[!, :land] = df_mask[:, :land]
    df[!, :dhdlon] = df_dem.dhdx
    df[!, :dhdlat] = df_dem.dhdy

    df_offset = DataFrame()
    df_offset[!, :dx] = fill(0.0, nrow(df))
    df_offset[!, :dy] = fill(0.0, nrow(df))
    df_offset[!, :dz] = fill(0.0, nrow(df))
    df_offset[!, :count] = fill(0.0, nrow(df))
    df_offset[!, :mad_offset] = fill(0.0, nrow(df))
    df_offset[!, :mad_ref] = fill(0.0, nrow(df))
    df_offset[!, :epsg] = fill(Int64(0), nrow(df))
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


"""
    geotile_read(geotile_dir, geotile_id, geotile_ext; buffer = nothing)

load geotile from geotile_dir. If buffer kwarg supplied [in decimal degrees] then all data 
surrounding geotile_id will also be returned 
"""

function geotile_read(geotile_dir, geotile_id, geotile_ext; buffer=nothing)
    if isnothing(buffer)
        fn = joinpath(geotile_dir, "$(geotile_id).$(geotile_ext)")
        df = copy(DataFrame(Arrow.Table(fn)))
        return df
    else
        # find tiles that are within buffer of center tile

        # get geotile extents
        nt = Altim.geotile_extent(geotile_id)

        # buffer extents 
        ext = nt2extent(nt)
        ext_buff = Extents.buffer(ext, buffer)

        geotile_width = nt.max_x - nt.min_x
        bwidth = ceil(maximum(buffer) / geotile_width) * geotile_width

        df = DataFrame()
        for min_x = (nt.min_x-bwidth):geotile_width:(nt.min_x+bwidth)
            for min_y = (nt.min_y-bwidth):geotile_width:(nt.min_y+bwidth)
                nt0 = (min_x=min_x, min_y=min_y, max_x=min_x + geotile_width, max_y=min_y + geotile_width)

                geotile_id0 = Altim.geotile_id(nt0)
                fn = joinpath(geotile_dir, "$(geotile_id0).$(geotile_ext)")

                # only read data if there is a file
                if isfile(fn)
                    df0 = DataFrame(Arrow.Table(fn))
                    ind = within.(Ref(ext_buff), df0.longitude, df0.latitude)
                    df = vcat(df, df0[ind, :])
                end
            end
        end
        return df
    end
end


"""
    points_in_polygons(points::AbstractMultiPoint, polygons::AbstractMultiPolygonTrait)

checks if each point in `points` is within any polygon in `polygons`, returning a Bool of 
length (points).
"""
function points_in_polygons(points, polygons)
    # initialize Bool array
    ind = falses(size(points))

    # check if point and polygon bounds overlap 
    ext0 = Extent(X=extrema(points.is[1]), Y=extrema(points.is[2]))
    ext = GeoInterface.bbox(polygons)

    if !Extents.intersects(ext0, ext)
        # no overlap
        return ind
    else
        # overlap 

        # check points in each polygon
        for p in polygons.shapes
            ext = GeoInterface.bbox(p)
            if Extents.intersects(ext0, ext)
                ind = ind .| LibGEOS.within.(points, Ref(p))
            end
        end
        return ind
    end
end


"""
    meters2lonlat_distance(distance_meters, latitude_degrees)

Returns the decimal degree distance along latitude and longitude lines given a distance in 
meters and a latitude in decimal degrees.

# Example usage:
```julia-repl
julia> distance_meters = 1000.0;  
julia> latitude_degrees = 45.0;  

julia> lat, lon = Altim.meters2lonlat_distance(distance_meters, latitude_degrees)
(0.008997741566866717, 0.012718328120254203)
```
"""
function meters2lonlat_distance(distance_meters, latitude_degrees)
    # Radius of the Earth in meters
    earth_radius = 6371000

    # Calculate the angular distance in radians
    angular_distance = distance_meters / earth_radius

    # Calculate the longitude distance using the Haversine formula
    longitude_distance = angular_distance * (180.0 / π) / cosd(latitude_degrees)

    latitude_distance = distance_meters / 111139

    return latitude_distance, longitude_distance
end

"""

"""
function binstats_tree(z, z0, halfwidth, v, fun::T) where {T<:Function}

    tree = KDTree(z; leafsize=10)
    idx = inrange(tree, z0, halfwidth)

    (n, ~) = size(v)
    (~, m) = size(z0)

    v_out = fill(NaN, [n, m])

    for j in eachindex(idx)
        for i = 1:n
            v_out[i, j] = fun(@view(v[i, idx[j]]))
        end
    end

    return v_out
end


function surface_parameters(z; lx=1, ly=1)

    ## slope
    dx(v) = ((v[7] + 2 * v[8] + v[9]) - (v[1] + 2 * v[2] + v[3])) / 8
    dy(v) = ((v[3] + 2 * v[6] + v[9]) - (v[1] + 2 * v[4] + v[7])) / 8

    # slope_x  = dx/lx
    # slope_y  = dx/ly

    ## planer curvature
    # https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/how-curvature-works.htm
    ddx(v) = ((v[2] + v[8]) / 2 - v[5]) # / lx^2
    ddy(v) = ((v[4] + v[6]) / 2 - v[5]) # / ly^2

    #Profile_Curvature = -2(ddx + ddy) * 100 where lx/y = grid_size in x and y

    stencil = Window(1)
    A = StencilArray(z, stencil)

    dzdx = mapstencil(dx, A)
    dzdy = mapstencil(dy, A)
    dzddx = mapstencil(ddx, A)
    dzddy = mapstencil(ddy, A)

    if lx != 1
        dzdx ./= lx
        dzddx ./= lx .^ 2
    end

    if ly != 1
        dzdy ./= ly
        dzddy ./= ly .^ 2
    end

    return dzdx, dzdy, dzddx, dzddy
end;


"""
    mask_rasterize(feature, extent::Extents, resolution; target_crs=EPSG(4326))
return a rasterized mask for provided features for a given extent and resolution.
"""
function mask_rasterize(feature, extent::Extent, resolution; target_crs=EPSG(4326), source_crs=EPSG(4326), boundary = :center)

    if source_crs !== target_crs
        println("reporjecting vector file... this could take awhile")
        feature = GeometryOps.reproject(feature; target_crs)
    end

    x = X(extent.X[1]:resolution[1]:extent.X[2]);
    y = Y(extent.Y[1]:resolution[2]:extent.Y[2]);
    A = Raster(zeros(UInt8, y, x))
    setcrs(A, target_crs)

    # NOTE: count method is fastest
    mask = Rasters.rasterize!(count, A, feature; threaded=false, shape=:polygon, progress=false, verbos=false, boundary = boundary) .> 0

    return mask
end

"""
    calculate the fraction of glacier, landice, floating ice, and rgi region in each geogrid
"""
function geotiles_w_mask(geotile_width; resolution = (X=.001, Y= .001), remake = false)

    if .!isdir(analysis_paths(; geotile_width).binned)
        mkpath(analysis_paths(; geotile_width).binned)
    end
    
    path2file = joinpath(analysis_paths(; geotile_width).binned, "geotiles_$(geotile_width)deg.arrow")

    if remake || !isfile(path2file)
        geotiles = project_geotiles(; geotile_width);
        masks = [:glacier, :landice, :floating]

        for mask in masks
            sym = Symbol("$(mask)_shp");
            feature = Shapefile.Handle(pathlocal[sym])

            sym = Symbol("$(mask)_frac")
            geotiles[!,sym] = zeros(nrow(geotiles))
            for geotile in eachrow(geotiles)
                mask = Altim.mask_rasterize(feature, geotile.extent, resolution; target_crs=EPSG(4326), source_crs=EPSG(4326), boundary = :center);
                geotile[sym] = sum(mask) / length(mask)
            end
        end

        regions = [:rgi6_regions]
        for region in regions
            sym = Symbol("$(region)_shp")
            feature = Shapefile.Table(pathlocal[sym])

            ids = unique(feature.RGI_CODE)
            rgi_ids =  "rgi" .* string.(ids)

            for rgi_id in rgi_ids
                geotiles[!,rgi_id] .= 0.
            end

            for geotile in eachrow(geotiles)
                for (id, rgi_id) in zip(ids, rgi_ids)
                    idx = feature.RGI_CODE .== id
                    mask = Altim.mask_rasterize(feature.geometry[idx], geotile.extent, (X=.01, Y= .01); target_crs=EPSG(4326), source_crs=EPSG(4326), boundary=:center)

                    if any(mask)
                        geotile[rgi_id] = sum(mask) / length(mask)
                    end
                end
            end
        end

        Arrow.write(path2file, geotiles::DataFrame)
    else
        geotiles = DataFrame(Arrow.Table(path2file))
        geotiles.extent = Extent.(getindex.(geotiles.extent, 1));
    end

    return geotiles
end


"""
    validrange(v)

Finds the enclosing range (start and end indices) of `true` elements along each dimension of a multidimensional boolean array.

**Arguments:**

* `v`: The input multidimensional boolean array.

**Returns:**

* A tuple of ranges, where each range represents the valid range (start and end indices) for a particular dimension of the input array.

**Examples:**

julia> A = [false true true false; false true false false; false false false false]
julia> validrange(A)
```julia
(1:2, 2:3)
```
"""
function validrange(v)
   nd = ndims(v);
   d = collect(1:nd)
   Tuple([findfirst(vec(any(v, dims=Tuple(x for x in d[d.!==n])))):findlast(vec(any(v, dims=Tuple(x for x in d[d.!==n])))) for n in 1:nd])
end


"""
curvature(dhddx, dhddy, epsg)

This function calculates the curvature of a surface at a point, 
accounting for the projection of the geographic coordinates.

**Arguments:**

* `dhddx`: The second derivative of the height field in the x-direction.
* `dhddy`: The second derivative of the height field in the y-direction.
* `epsg`: The EPSG code of the geographic coordinate system. Currently supports
  WGS 84 (EPSG:4326) and non-geographic projections.

**Returns:**

* `curvature`: The curvature of the surface at the point, in units of cm^-2. 
  A negative value indicates a concave surface, while a positive value indicates 
  a convex surface.

**Notes:**

* For geographic projections (EPSG:4326), the function calculates the distance 
  represented by one degree of latitude and longitude in meters at the 
  mean latitude of the geospatial tile. This distance is then used to convert 
  the curvature from units per meter squared to units per centimeter squared.
* For non-geographic projections, the function assumes a constant distance of 
  one meter for both latitude and longitude. This may not be accurate for all 
  projections, and the results should be interpreted with caution.
"""
function curvature(dhddx, dhddy, epsg; lat)

    if epsg == 4326
        ydist, xdist = meters2lonlat_distance.(1, lat) # degree per meter
    else
        ydist = 1
        xdist = 1
    end
    curvature = -2(dhddx * xdist ^ 2 + dhddy * ydist ^ 2) * 100;

    return curvature
end


function slope(dhdx, dhdy, epsg; lat)

    if epsg == 4326
        ydist, xdist = meters2lonlat_distance.(1, lat) # degree per meter
    else
        ydist = 1
        xdist = 1
    end

    dhdx = dhdx * xdist
    dhdy = dhdy * ydist
   
    return dhdx,  dhdy
end


"""
vector_overlap(a, b)

This function checks for overlap between two one-dimensional arrays `a` and `b`. 
Overlap is defined as a condition where the minimum value of `a` is less than 
or equal to the maximum value of `b`, and the maximum value of `a` is greater than or equal 
to the minimum value of `b`. In other words, `a` completely encompasses the range of values in `b`.

Inputs:
  - a: The first one-dimensional array.
  - b: The second one-dimensional array.

Output:
  - tf (Bool): True if there is vectors overlap, False otherwise.
"""
function vector_overlap(a, b)

    ea = extrema(a)
    eb = extrema(b)

    tf = (ea[2] >= eb[1]) && (ea[1] <= eb[2])
    return tf
end


function nanmean(x)
    mean(x[.!isnan.(x)])
end

function nanmedian(x)
    valid = .!isnan.(x)
    if any(valid)
        y = median(x[valid])
    else
        y = NaN
    end
    return y
end


const MATLAB_EPOCH = Dates.DateTime(-1, 12, 31)
function datenum2date(n) 
    d = MATLAB_EPOCH + Dates.Millisecond(round(Int64, n * 1000 * 60 * 60 * 24))
    return d
end



"""
    add_dem_ref!(altim, dem_id, geotile, mission_geotile_folder)

    Add dem info to altim dataframe
"""
function add_dem_ref!(altim, dem_id, geotile, mission_geotile_folder)

    # add dem height and curvature
    if dem_id == :best
        # last dem takes precidence over earlier dems
        dem_id0 = [:cop30_v2, :arcticdem_v4_10m, :rema_v2_10m]
    elseif any([:cop30_v2, :arcticdem_v4_10m, :rema_v2_10m] .== dem_id)
        dem_id0 = [dem_id]
    else
        error("unreconized dem_id: $dem_id")
    end

    altim[!, :height_ref] .= NaN
    altim[!, :curvature] .= NaN

    for dem_id1 in dem_id0
        path2dem = joinpath(mission_geotile_folder, geotile.id * ".$(dem_id1)")
        if isfile(path2dem)
            dem = select!(DataFrame(Arrow.Table(path2dem)), :height, :dhddx, :dhddy)
            # a bit faster to calculate curvature on all data then subset
            curv = Altim.curvature.(dem.dhddx, dem.dhddy, Ref(Altim.dem_info(dem_id1)[1].epsg), lat=mean(geotile.extent.Y))

            ind = .!isnan.(curv) .& (abs.(dem.height) .< 9998)

            altim[ind, :height_ref] = dem[ind, :height]
            altim[ind, :curvature] = curv[ind]
        end
    end

    altim[!, :dh] = altim.height .- altim.height_ref
    return altim
end


function highres_mask(geotile, feature, invert, excludefeature)
    # update mask with high-resolution vector files
    grid_resolution = 0.00027 # ~30m

    x_mask = X(geotile.extent.X[1]:grid_resolution:geotile.extent.X[2],
        sampling=DimensionalData.Intervals(DimensionalData.Start()))
    y_mask = Y(geotile.extent.Y[1]:grid_resolution:geotile.extent.Y[2],
        sampling=DimensionalData.Intervals(DimensionalData.Start()))

    mask1 = Raster(zeros(UInt8, y_mask, x_mask))
    setcrs(mask1, EPSG(4326))

    # NOTE: count method is fastest
    mask1 = Rasters.rasterize!(count, mask1, feature; threaded=false,
        shape=:polygon, progress=false, verbose=false, boundary=:center) .> 0

    if invert
        mask1 = .!(mask1)
    end

    if !isnothing(excludefeature)
        excludemask = Raster(zeros(UInt8, y_mask, x_mask))
        excludemask = Rasters.rasterize!(count, excludemask, excludefeature; threaded=false, shape=:polygon, progress=false, verbose=false) .> 0
        mask1 = mask1 .& .!excludemask
    end

    # calculate area per cell
    lon = lookup(mask1, X)
    lat = lookup(mask1, Y)
    d = Altim.meters2lonlat_distance.(Ref(1), lat)
    a = abs.((1 ./ getindex.(d, 2) * (lat[2] .- lat[1])) .* (1 / d[1][1] * (lon[2] - lon[1])))
    area_m2 = repeat(a', outer=[length(lon), 1])

    return (mask1, area_m2)
end


function highres_mask!(masks0, altim, geotile, feature, invert, excludefeature, surface_mask)
    
    
    (mask1, _) = highres_mask(geotile, feature, invert, excludefeature)
    
    valid = Altim.within.(Ref(geotile.extent), altim.longitude, altim.latitude)
    masks0[!, surface_mask] .= false

    fast_index = true
    if fast_index # fast index is 15x faster than Rasters
        c = floor.(Int64, (altim.longitude[valid] .- first(x_mask)) ./ step(x_mask)) .+ 1
        r = floor.(Int64, (altim.latitude[valid] .- first(y_mask)) ./ step(y_mask)) .+ 1
        pts3 = [CartesianIndex(a, b) for (a, b) in zip(r, c)]
        masks0[valid, surface_mask] .= mask1[pts3]
    else
        #NOTE: 67% of time is taken here for large number of points.
        pts1 = GeoInterface.PointTuple.([((Y=y, X=x)) for (x, y) in
                                         zip(altim.longitude[valid], altim.latitude[valid])])
        pts2 = extract(mask1, pts1, atol=grid_resolution / 2, index=true, geometry=false)
        masks0[:, surface_mask][valid] = getindex.(pts2, 2)
        pts3 = CartesianIndex.(getindex.(pts2, 1))
    end

    return masks0
end


function binningfun_define(binning_method)
    
    if binning_method == "meanmadnorm3"
      x -> mean(x[Altim.madnorm(x).<3])
    elseif binning_method == "meanmadnorm5"
        x -> mean(x[Altim.madnorm(x).<5])
    elseif binning_method == "meanmadnorm10"
        x -> mean(x[Altim.madnorm(x).<10])
    elseif binning_method == "median"
        x -> median(x)
    else
        error("unrecognized binning method")
    end


end



function highres_mask(latitude, longitude, feature; grid_resolution=0.00027)
    # update mask with high-resolution vector files

    ymin = floor(minimum(latitude) ./ grid_resolution) .* grid_resolution
    ymax = ceil(maximum(latitude) ./ grid_resolution) .* grid_resolution
    xmin = floor(minimum(longitude) ./ grid_resolution) .* grid_resolution
    xmax = ceil(maximum(longitude) ./ grid_resolution) .* grid_resolution
    
    x_mask = X(xmin:grid_resolution:xmax, sampling=DimensionalData.Intervals(DimensionalData.Start()))
    y_mask = Y(ymin:grid_resolution:ymax, sampling=DimensionalData.Intervals(DimensionalData.Start()))

    mask1 = Raster(zeros(UInt8, y_mask, x_mask))

    # NOTE: count method is fastest
    mask1 = Rasters.rasterize!(count, mask1, feature; threaded=false, shape=:polygon, progress=false, verbose=false, boundary=:center) .> 0

    fast_index = true
    if fast_index # fast index is 15x faster than Rasters
        c = floor.(Int64, (longitude .- first(x_mask)) ./ step(x_mask)) .+ 1
        r = floor.(Int64, (latitude .- first(y_mask)) ./ step(y_mask)) .+ 1
        pts3 = [CartesianIndex(a, b) for (a, b) in zip(r, c)]
        mask = mask1[pts3]
    else
        #NOTE: 67% of time is taken here for large number of points.
        pts1 = GeoInterface.PointTuple.([((Y=y, X=x)) for (x, y) in
                                         zip(longitude, latitude)])
        pts2 = extract(mask1, pts1, atol=grid_resolution / 2, index=true, geometry=false)
        mask = getindex.(pts2, 2)
    end

    return mask
end


# remove all rows that contian only nans
function remove_allnan(df, column)
    allnan = falses(nrow(df))
    for (i, r) in enumerate(eachrow(df))
        allnan[i] = all(isnan.(r[column])) || isempty(r[column])
    end
    df = df[.!allnan,:]
end


function nanquantile(itr, p;)
    valid = .!isnan.(itr)
    if any(valid)
        q = quantile(itr[valid], p)
    else
        q = fill(NaN, size(p))
    end
end


mission2label = Dict( 
    "hugonnet" => "ASTER",
    "icesat" => "ICESat",
    "icesat2" => "ICESat-2",
    "gedi" => "GEDI",
    "grace" => "GRACE/-FO",
    "synthesis" => "Synthesis",
    "gemb" => "GEMB Model",
)

rgi2label = Dict(
    "rgi1" => "Alaska",
    "rgi2" => "Western Canada",
    "rgi3" => "Arctic Canada N.",
    "rgi4" => "Arctic Canada S.",
    "rgi5" => "Greenland",
    "rgi6" => "Iceland",
    "rgi7" => "Svalbard",
    "rgi8" => "Scandinavia",
    "rgi9" => "Russian Arctic",
    "rgi10" => "North Asia",
    "rgi11" => "Central Europe",
    "rgi12" => "Caucasus",
    "rgi13" => "Central Asia",
    "rgi14" => "South Asia W.",
    "rgi15" => "South Asia E.",
    "rgi16" => "Low Latitudes",
    "rgi17" => "Southern Andes",
    "rgi18" => "New Zealand",
    "rgi19" => "Antarctic",
    "hma" => "High Mtn. Asia",
    "global" => "Global",
    "global_ep" => "Global [excluding periphery]",
    "ghll" => "All Missions: rgi 2, 10, 11, 12, 13, 14, 15, 17, 18",
    "hll" => "ASTER, ICESat/-2: rgi 1, 3, 4, 5, 6, 7, ,8, 9, 19",
    "hll_ep" => "ASTER, ICESat/-2 no periphery: rgi 1, 3, 4, 6, 7, ,8, 9",
)

unit2areaavg = Dict(
    "km⁻³"=>"m", 
    "Gt"=>"m w.e."
)

units2label = Dict(
    "m w.e." => "water equivlent height anomaly", 
    "Gt" => "mass anomaly", 
    "km⁻³" => "volume anomaly",  
    "m" => "height anomaly"
)


var2label = Dict(
    "dm" => "Mass Anomaly",
    "dv" => "Volume Anomaly",
    "fac" => "Firn Air Content Anomaly",
    "smb" => "Total Surface Mass Balance",
    "acc" => "Total Accumulation",
    "runoff" => "Total Runoff",
    "refreeze" => "Total Refreeze",
)

"""
ntpermutations(nt)

given a named tuple [nt] return all a vector of named tuples of all permutiations of nt values, retaining original names

Inputs:
  - nt: a named tuple of iterables
Output:
  - vector of named tuples with all permutations of input named tuple of iterables

  **Examples:**

julia> nt_in = (; category = [1, 3, 4], letter = ["a", "b", "c"], score = [1.1, 4.0], thing = ["hat", "cat"])
julia> nt = ntpermutations(nt_in)
```julia
36-element Vector{NamedTuple}:
 (catagory = 1, letter = "a", score = 1.1, thing = "hat")
 ⋮
 (catagory = 4, letter = "c", score = 4.0, thing = "cat")
```
"""
ntpermutations(nt_in) = [NamedTuple{keys(nt_in)}(reverse(t)) for t in Iterators.product(reverse(nt_in)...)] |> vec



function extent2rectangle(extent)
    xbounds = extent[1]
    ybounds = extent[2]
    rectangle = GeoInterface.Wrappers.Polygon([[(xbounds[1], ybounds[1]), (xbounds[1], ybounds[2]), (xbounds[2], ybounds[2]), (xbounds[2], ybounds[1]), (xbounds[1], ybounds[1])]])
    return rectangle
end
