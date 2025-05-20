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

Constructs file paths for mission data organization.

# Arguments
- `geotile_width`: Geotile width in degrees
- `mission`: Mission name (e.g., "GEDI", "ICESat2")
- `product`: Product identifier
- `version`: Version number

# Returns
Named tuple with data paths (raw_data, geotile, granules_remote, granules_local)
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

Calculate the bounding box of a geotile centered at given coordinates.

# Arguments
- `lat`: Center latitude in degrees
- `lon`: Center longitude in degrees
- `width`: Width of the geotile in degrees

# Returns
NamedTuple with min_x, min_y, max_x, max_y defining the geotile boundaries
"""
function geotile_extent(lat, lon, width)
    (min_x=(lon - width / 2), min_y=(lat - width / 2),
        max_x=(lon + width / 2), max_y=(lat + width / 2))
end

"""
    within(extent, x, y)

Check if point (x,y) falls within the specified extent.

# Arguments
- `extent`: Rectangular boundary (NamedTuple or Extent)
- `x`: X-coordinate (longitude)
- `y`: Y-coordinate (latitude)

# Returns
Boolean indicating whether the point is inside the extent
"""
function within(extent::NamedTuple{(:min_x, :min_y, :max_x, :max_y)}, x, y)
    in = (x >= extent.min_x) .&& (x <= extent.max_x) .&& (y >= extent.min_y) .&& (y <= extent.max_y)
    return in
end

function within(extent::Extent, x, y)
    in = (x >= extent.X[1]) .&& (x <= extent.X[2]) .&& (y >= extent.Y[1]) .&& (y <= extent.Y[2])
    return in
end


# seach directory using two keys
searchdir(path, key1, key2) = filter(x -> (occursin(key1, x) .& occursin(key2, x)), readdir(path))
searchdir(path, key1) = filter(x -> occursin(key1, x), readdir(path))


"""
    geotile_define(geotile_width::Number)

Create a global grid of geotiles with specified width.

# Arguments
- `geotile_width`: Width of each geotile in degrees (must divide evenly into 180)

# Returns
- `DataFrame`: Table with columns for geotile ID and extent (min_x, min_y, max_x, max_y)

# Throws
- `ErrorException`: If the geotile width doesn't divide evenly into 180 degrees
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
function geotile_build(geotile_granules, geotile_dir; warnings=true, fmt=:arrow, replace_corrupt_h5 = true)
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
    geotile_id(extent)

Generate a standardized ID string for a geotile based on its geographic extent.

# Arguments
- `extent`: NamedTuple containing min_x, min_y, max_x, max_y coordinates in degrees

# Returns
- String in format "lat[+YY+YY]lon[+XXX+XXX]" where YY are latitudes and XXX are longitudes
"""
function geotile_id(extent)
    id = @sprintf("lat[%+03.0f%+03.0f]lon[%+04.0f%+04.0f]", extent.min_y, extent.max_y, extent.min_x, extent.max_x)
    return id
end

"""
    geotile_extent(geotile_id::String) -> NamedTuple

Extract geographic extent from a geotile ID string.

# Arguments
- `geotile_id`: String in format "lat[+YY+YY]lon[+XXX+XXX]" where YY are latitudes and XXX are longitudes

# Returns
- NamedTuple with min_x, min_y, max_x, max_y coordinates in degrees
"""
function geotile_extent(geotile_id)
    min_y = parse(Int, geotile_id[5:7])
    max_y = parse(Int, geotile_id[8:10])

    min_x = parse(Int, geotile_id[16:19])
    max_x = parse(Int, geotile_id[20:23])

    return (min_x=min_x, min_y=min_y, max_x=max_x, max_y=max_y)
end


"""
    geotile_utm!(df::DataFrame; height=nothing) -> (DataFrame, String)

Convert geographic coordinates (longitude, latitude) to local UTM or polar stereographic 
coordinates and add them as X and Y columns to the DataFrame.

# Arguments
- `df`: DataFrame containing longitude and latitude columns
- `height`: Optional height parameter (not used in current implementation)

# Returns
- Tuple containing the modified DataFrame and the EPSG code string for the projection
"""
function geotile_utm!(df::DataFrame; height=nothing)
    # create a regularly spaced grid using a local low distortion projection
    epsg = utm_epsg(mean(df.longitude), mean(df.latitude))
    df[!, :X], df[!, :Y] = epsg2epsg(df.longitude, df.latitude, "EPSG:4326", "$epsg", parse_output=true)

    return df, epsg
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
    decimalyear(datetime)

Convert DateTime to decimal year representation.

# Arguments
- `datetime`: DateTime object to convert

# Returns
- Float representing year with decimal fraction (e.g., 2018.813698630137 for 2018-10-24)
"""
function decimalyear(datetime)
    year = Dates.year(datetime)
    day_of_year = Dates.dayofyear(datetime)
    days_in_year = Dates.daysinyear(year)
    decyear = year + (day_of_year / days_in_year)
    return decyear
end


"""
    decimalyear2datetime(decyear)

Convert decimal year to DateTime object.

# Arguments
- `decyear`: Decimal year (e.g., 2018.813698630137)

# Returns
- DateTime object (e.g., DateTime("2018-10-24"))
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
    regular_grid(center_extent, node_spacing; node_width=node_spacing, node_shape="square")

Create a regular grid of nodes within a specified extent.

# Arguments
- `center_extent`: Extent object with x_min, x_max, y_min, y_max properties
- `node_spacing`: Distance between adjacent node centers

# Keywords
- `node_width`: Width of each node (defaults to node_spacing)
- `node_shape`: Shape of the nodes ("square" or other shapes)

# Returns
- NamedTuple containing:
  - `node_center`: Grid of node centers
  - `x_node_center`: Vector of x-coordinates for node centers
  - `y_node_center`: Vector of y-coordinates for node centers
  - `node_half_width`: Half-width of each node
  - `node_shape`: Shape of the nodes
"""
function regular_grid(center_extent, node_spacing; node_width=node_spacing, node_shape="square")
    x_node_center = center_extent.x_min:node_spacing:center_extent.x_max
    y_node_center = center_extent.y_min:node_spacing:center_extent.y_max
    node_half_width = node_width / 2

    node_center = ndgrid(y_node_center, x_node_center)
    (; node_center, x_node_center, y_node_center, node_half_width, node_shape)
end



"""
    madnorm(x)

Calculate the normalized median absolute deviation (MAD) of values in x.

# Arguments
- `x`: Vector of values to analyze

# Returns
- Vector of normalized deviations, where each value represents the number of 
  standard deviations from the median (using MAD scaled by 1.4826)
"""
function madnorm(x)
    consistent_estimator = 1.4826 #mad to sigma conversion factor
    x_abs = abs.(x .- median(x))
    x_madnorm = x_abs ./ (median(x_abs) .* consistent_estimator)
    return x_madnorm
end

"""
    mad(x)

Calculate the median absolute deviation (MAD) of values in x.

# Arguments
- `x`: Vector of values to analyze

# Returns
- The median of absolute deviations from the median, or NaN if x is empty
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

Calculate the range (maximum minus minimum) of a collection of values.

# Arguments
- `x`: Collection of values

# Returns
- The difference between the maximum and minimum values in x
"""
function range(x)
    x_minmax = extrema(x)
    r = x_minmax[2] - x_minmax[1]
    return r
end


"""
    utm_epsg(lon::Real, lat::Real; always_xy=true) -> String

Get the EPSG code for the UTM zone containing the coordinates, or polar stereographic projection if outside UTM limits.

# Arguments
- `lon`: Longitude coordinate
- `lat`: Latitude coordinate
- `always_xy`: If true, interpret inputs as (lon, lat); if false, as (lat, lon)

# Returns
- String containing the EPSG code in the format "EPSG:XXXXX"

# Notes
- Uses EPSG:3413 (NSIDC Sea Ice Polar Stereographic North) for latitudes above 84°N
- Uses EPSG:3031 (Antarctic Polar Stereographic) for latitudes below 80°S
- Handles special cases for Norway and Svalbard
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

Returns a GeoArray containing geoid height data relative to WGS 84 (EPSG:4979) ellipsoid.

Downloads the geoid file from https://www.agisoft.com/downloads/geoids/ if not found locally.

# Arguments
- `geoid::String`: Geoid model name, one of "egm84", "egm96", or "egm2008"
- `folder::String`: Directory where geoid files are stored (default: "~/Documents/geoids")

# Returns
- `GeoArray`: Geoid height data

# Note
Models are converted from USA NGA data by Agisoft under Public Domain license.
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
    normalize(data) -> NamedTuple

Normalize data by subtracting the mean and dividing by the standard deviation.

# Arguments
- `data`: Array of values to normalize

# Returns
- NamedTuple containing:
  - `std`: Standard deviation of the original data
  - `mean`: Mean of the original data
  - `norm`: Normalized data array
"""
function normalize(data)
    m = mean(data)
    s = std(data)
    n = (data .- m) ./ s

    datan = (std=s, mean=m, norm=n)

    return datan
end


"""
    centroid(ga::GeoArray) -> NamedTuple

Calculate the centroid coordinates of a GeoArray.

# Arguments
- `ga::GeoArray`: The GeoArray to find the centroid of

# Returns
- NamedTuple with `x` and `y` coordinates of the centroid
"""
function centroid(ga::GeoArray)
    bbox = extent2nt(GeoArrays.bbox(ga))
    (x=(bbox.min_x + bbox.max_x) / 2, y=(bbox.min_y + bbox.max_y) / 2)
end


"""
    xy_gridsize(ga::GeoArray) -> NamedTuple

Calculate the grid cell size in meters for a GeoArray.

# Arguments
- `ga::GeoArray`: The GeoArray to calculate grid size for

# Returns
- NamedTuple with `x` and `y` grid cell dimensions in meters
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
    dh_xy2aa(x, y, dhdx, dhdy, epsg_pts, epsg_dem) -> (along, across)

Convert slope components from x-y coordinates to along-track and across-track directions.

# Arguments
- `x`: Vector of x-coordinates
- `y`: Vector of y-coordinates
- `dhdx`: Slope in x-direction
- `dhdy`: Slope in y-direction
- `epsg_pts`: EPSG code of input coordinates
- `epsg_dem`: EPSG code of DEM coordinates

# Returns
- Tuple of (along-track, across-track) slope components in meters
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
    dh_ll2aa(lon, lat, dhdlon, dhdlat) -> (along, across)

Convert slope components from longitude-latitude coordinates to along-track and across-track directions.

# Arguments
- `lon`: Vector of longitude coordinates
- `lat`: Vector of latitude coordinates
- `dhdlon`: Slope in longitude direction (per degree)
- `dhdlat`: Slope in latitude direction (per degree)

# Returns
- Tuple of (along-track, across-track) slope components in meters
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
    geotile_extract_mask(
        geotile_id::String,
        geotile_dir::String;
        vars::Vector{Symbol}=[:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean],
        geotile_ext::String=".arrow",
        job_id="",
        force_remake=false
    ) -> Nothing

Extract mask data for a specific geotile and save as an Arrow file.

# Arguments
- `geotile_id::String`: Identifier for the geotile
- `geotile_dir::String`: Directory containing geotile files

# Keywords
- `vars::Vector{Symbol}`: Mask variables to extract (default: [:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean])
- `geotile_ext::String`: File extension of geotile files (default: ".arrow")
- `job_id::String`: Identifier for logging purposes (default: "")
- `force_remake::Bool`: Whether to regenerate existing files (default: false)

# Returns
Nothing, but creates a mask file for the geotile with extracted data
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
        printstyled("\n    ->$job_id $geotile_id masks already exists, skipping\n"; color=:green)
    elseif isfile(path2geotile)
        t1 = time()

        df0 = DataFrame(Arrow.Table(path2geotile))[!, [:longitude, :latitude]]
        #vlength = length.(df0.longitude)

        #lon = vcat(df0.longitude...)
        #lat = vcat(df0.latitude...)

        if isempty(df0.longitude)
            printstyled("\n    ->$job_id $geotile_id is all empty, skipping\n"; color=:light_red)
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
        printstyled("\n    ->$job_id $geotile_id masks extracted: $(total_time) min \n"; color=:light_black)
    else
        printstyled("\n    ->$job_id $geotile_id does not exist\n"; color=:yellow)
        return
    end
end

"""
    geotile_extract_mask(geotiles::DataFrame, geotile_dir::String; vars, job_id, force_remake) -> Nothing

Extract mask data for multiple geotiles and save as Arrow files.

# Arguments
- `geotiles::DataFrame`: DataFrame containing geotile information with 'id' column
- `geotile_dir::String`: Directory containing geotile files

# Keywords
- `vars::Vector{Symbol}=[:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean]`: Mask variables to extract
- `job_id::String=""`: Identifier for logging purposes
- `force_remake::Bool=false`: Whether to regenerate existing files

# Returns
Nothing, but creates mask files for each geotile
"""
function geotile_extract_mask(
    geotiles::DataFrame,
    geotile_dir::String;
    vars::Vector{Symbol}=[:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean],
    job_id="",
    force_remake=false
)

    #asyncmap(eachrow(geotiles); ntasks = 4) do geotile
    @showprogress dt = 10 desc = "extracting masks for $job_id geotiles..." for geotile in eachrow(geotiles)
        geotile_extract_mask(geotile.id, geotile_dir; vars=vars, job_id=job_id, force_remake=force_remake)
    end
end


"""
    _offset_design_matrix(dhdx, dhdy, p)
interanl function for building the design matrix to solve for track offset
"""
_offset_design_matrix(dhdx, dhdy) = hcat(dhdx, dhdy, ones(size(dhdx)))

"""
    track_offset(dhdx, dhdy, dh; kwargs...) -> (dx, dy, dz, count, mad_offset, mad_ref)

Calculate optimal track offsets that minimize height anomalies using surface slopes.

# Arguments
- `dhdx`: Slope in x-direction
- `dhdy`: Slope in y-direction
- `dh`: Observed height anomalies (difference between measured and reference heights)

# Keywords
- `weights=nothing`: Optional weights for regression
- `regressor=HuberRegression(fit_intercept=false, scale_penalty_with_samples=false)`: Regression method
- `interations=1`: Number of iterations for robust fitting
- `iter_thresh=5`: Threshold for outlier removal between iterations (in MAD units)
- `verbose=true`: Whether to print progress information

# Returns
- `dx`: Optimal x-offset [m]
- `dy`: Optimal y-offset [m]
- `dz`: Optimal height offset [m]
- `count`: Number of valid points used in final iteration
- `mad_offset`: Median absolute deviation after correction
- `mad_ref`: Median absolute deviation of original height anomalies
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
    track_offset_dh(dhdx, dhdy, dx, dy, dz) -> dh

Calculate height anomaly based on surface slopes and positional offsets.

# Arguments
- `dhdx`: Slope in x-direction
- `dhdy`: Slope in y-direction
- `dx`: Offset in x-direction
- `dy`: Offset in y-direction
- `dz`: Vertical offset

# Returns
- `dh`: Calculated height anomaly
"""
function track_offset_dh(dhdx, dhdy, dx, dy, dz)
    dh = dx .* dhdx .+ dy .* dhdy .+ dz
    return dh
end

"""
    geotile_track_offset(geotile_id, geotile_dir, dem; valid_count_thresh=100, force_remake=false)

Calculate positional offsets for a geotile by comparing altimetry data with DEM elevations.

# Arguments
- `geotile_id::String`: Identifier for the geotile
- `geotile_dir::String`: Directory containing geotile files
- `dem::Symbol`: DEM source to use (e.g., `:arcticdem`)

# Keywords
- `valid_count_thresh=100`: Minimum number of valid points required for offset calculation
- `force_remake=false`: Whether to regenerate existing offset files

# Returns
Nothing, but creates an offset file for the geotile with calculated positional corrections
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
        deminfo, goids_folder = dem_info(:arcticdem_v3_10m)

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
        geotiles::DataFrame,
        geotile_dir::String,
        dem::Symbol;
        valid_count_thresh=100,
        force_remake=false
    )

Calculate track offsets for multiple geotiles using a digital elevation model.

# Arguments
- `geotiles::DataFrame`: DataFrame containing geotile information with an 'id' column
- `geotile_dir::String`: Directory containing geotile files
- `dem::Symbol`: Digital elevation model to use for offset calculation

# Keywords
- `valid_count_thresh=100`: Minimum number of valid points required for offset calculation
- `force_remake=false`: Whether to recalculate offsets even if they already exist

# Returns
Nothing, but creates offset files for each geotile in the specified directory
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
    nt2extent(nt::NamedTuple{(:min_x, :min_y, :max_x, :max_y)}) -> Extent

Convert a NamedTuple with geographic bounds to an Extent object.

# Arguments
- `nt`: NamedTuple with min_x, min_y, max_x, max_y fields defining a bounding box

# Returns
- `Extent` object with the same geographic bounds
"""
function nt2extent(nt::NamedTuple{(:min_x, :min_y, :max_x, :max_y)})
    Extent(X=(nt.min_x, nt.max_x), Y=(nt.min_y, nt.max_y))
end

"""
    extent2nt(extent::Extent) -> NamedTuple

Convert an Extent object to a NamedTuple with geographic bounds.

# Arguments
- `extent::Extent`: Extent object containing X and Y coordinate ranges

# Returns
- NamedTuple with min_x, min_y, max_x, max_y fields defining a bounding box
"""
function extent2nt(extent::Extent)
    (min_x=minimum(extent.X), min_y=minimum(extent.Y), max_x=maximum(extent.X), max_y=maximum(extent.Y))
end


"""
    geotile_subset!(geotiles, extent::Extent) -> DataFrame

Filter a DataFrame of geotiles to only those that intersect with the given extent.

# Arguments
- `geotiles`: DataFrame containing geotiles with an `extent` column
- `extent::Extent`: Geographic extent to test for intersection

# Returns
- Filtered DataFrame containing only geotiles that intersect with the given extent
"""
function geotile_subset!(geotiles, extent::Extent)
    ind = .!(Extents.intersects.(Ref(extent), geotiles.extent))
    geotiles = deleteat!(geotiles, ind)
    return geotiles
end

"""
    allfiles(rootdir; kwargs...) -> Vector{String}

Return a vector of file paths that meet specified criteria.

# Arguments
- `rootdir`: Directory from which to begin search

# Keywords
- `subfolders=false`: Whether to include subfolders in the search
- `fn_startswith=nothing`: Include files that start with specified string(s)
- `fn_endswith=nothing`: Include files that end with specified string(s)
- `fn_contains=nothing`: Include files that contain specified string(s)
- `and_or=&`: Use `&` to require all criteria or `|` to match any criteria
- `topdown=true`: Whether to traverse directories top-down
- `follow_symlinks=false`: Whether to follow symbolic links
- `onerror=throw`: Function to handle errors during directory traversal

# Returns
- Vector of file paths matching the specified criteria
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
                containsin_tf = any(Base.contains.(file, fn_contains))
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
                    containsin_tf = any(Base.contains.(file, fn_contains))
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
    offset_geo(lon, lat, dh, dhdlon, dhdlat; interations=3, iter_thresh=7, verbose=true)

Calculate geographic offsets in local UTM coordinates based on elevation anomalies and slopes.

# Arguments
- `lon`: Longitude coordinates in degrees
- `lat`: Latitude coordinates in degrees
- `dh`: Elevation anomalies in meters
- `dhdlon`: Slope in longitude direction (m per degree)
- `dhdlat`: Slope in latitude direction (m per degree)

# Keywords
- `interations=3`: Number of iterations for offset calculation
- `iter_thresh=7`: Threshold for iteration convergence
- `verbose=true`: Whether to print progress information

# Returns
- `dx`: X-offset in meters
- `dy`: Y-offset in meters
- `dz`: Z-offset in meters
- `cnt`: Count of valid data points used
- `mad_offset`: Median absolute deviation of offset
- `mad_ref`: Median absolute deviation of reference
- `epsg`: EPSG code for the UTM zone
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
    geotile_coregister(geotile, geotile_dir, dem)

Coregister geotile data to a digital elevation model (DEM).

# Arguments
- `geotile`: Geotile information containing an `id` field
- `geotile_dir`: Directory containing geotile files
- `dem`: Symbol or string identifying the DEM to coregister against

# Returns
Nothing, but creates an offset file with the coregistration parameters
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
    offset_geo(lon, lat, dhdlon, dhdlat, dx, dy, dz, epsg)

Calculate elevation change due to horizontal and vertical offsets.

# Arguments
- `lon`: Longitude coordinates in degrees
- `lat`: Latitude coordinates in degrees
- `dhdlon`: Slope in longitude direction (m per degree)
- `dhdlat`: Slope in latitude direction (m per degree)
- `dx`: X offset in meters (local UTM coordinates)
- `dy`: Y offset in meters (local UTM coordinates)
- `dz`: Z offset in meters
- `epsg`: EPSG code for local UTM projection

# Returns
- Vector of elevation changes (in meters) resulting from the specified offsets
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
    meters2lonlat_distance(distance_meters, latitude_degrees) -> (latitude_distance, longitude_distance)

Convert a distance in meters to equivalent decimal degree distances in latitude and longitude.

# Arguments
- `distance_meters`: Distance in meters
- `latitude_degrees`: Latitude in decimal degrees where the conversion is needed

# Returns
- Tuple containing (latitude_distance, longitude_distance) in decimal degrees
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
    mask_rasterize(feature, extent::Extent, resolution; target_crs=EPSG(4326), source_crs=EPSG(4326), boundary=:center) -> BitMatrix

Convert vector features to a binary raster mask for a given geographic extent and resolution.

# Arguments
- `feature`: Vector geometry feature(s) to rasterize
- `extent::Extent`: Geographic extent of the output raster
- `resolution`: Tuple or NamedTuple specifying pixel size (e.g., (X=0.001, Y=0.001))

# Keywords
- `target_crs=EPSG(4326)`: Coordinate reference system for the output raster
- `source_crs=EPSG(4326)`: Coordinate reference system of the input features
- `boundary=:center`: Method for determining pixel coverage (`:center`, `:touch`, or `:contains`)

# Returns
- `BitMatrix`: Binary mask where `true` indicates presence of the feature
"""
function mask_rasterize(feature, extent::Extent, resolution; target_crs=EPSG(4326), source_crs=EPSG(4326), boundary=:center)

    if source_crs !== target_crs
        println("reporjecting vector file... this could take awhile")
        feature = GeometryOps.reproject(feature; target_crs)
    end

    x = X(extent.X[1]:resolution[1]:extent.X[2]);
    y = Y(extent.Y[1]:resolution[2]:extent.Y[2]);
    A = Raster(zeros(UInt8, y, x))
    setcrs(A, target_crs)

    # NOTE: count method is fastest
    mask = Rasters.rasterize!(count, A, feature; threaded=false, shape=:polygon, progress=false, verbos=false, boundary=boundary) .> 0

    return mask
end

"""
    geotiles_w_mask(geotile_width; resolution=(X=.001, Y=.001), remake=false) -> DataFrame

Calculate the fraction of glacier, land ice, floating ice, and RGI regions in each geotile.

# Arguments
- `geotile_width`: Width of geotiles in degrees

# Keywords
- `resolution=(X=.001, Y=.001)`: Spatial resolution for rasterization
- `remake=false`: Whether to regenerate the file even if it already exists

# Returns
- `DataFrame`: Contains geotile information with added columns for coverage fractions
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
                mask = mask_rasterize(feature, geotile.extent, resolution; target_crs=EPSG(4326), source_crs=EPSG(4326), boundary = :center);
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
                    mask = mask_rasterize(feature.geometry[idx], geotile.extent, (X=.01, Y= .01); target_crs=EPSG(4326), source_crs=EPSG(4326), boundary=:center)

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

Find the enclosing range of `true` elements along each dimension of a boolean array.

# Arguments
- `v`: Multidimensional boolean array

# Returns
- Tuple of ranges, one per dimension, containing the first to last indices where `true` values exist
"""
function validrange(v)
   nd = ndims(v);
   d = collect(1:nd)
   Tuple([findfirst(vec(any(v, dims=Tuple(x for x in d[d.!==n])))):findlast(vec(any(v, dims=Tuple(x for x in d[d.!==n])))) for n in 1:nd])
end


"""
    curvature(dhddx, dhddy, epsg; lat) -> Float64

Calculate surface curvature at a point, accounting for coordinate projection.

# Arguments
- `dhddx`: Second derivative of height in x-direction
- `dhddy`: Second derivative of height in y-direction
- `epsg`: EPSG code of coordinate system
- `lat`: Latitude value(s) for geographic conversion

# Returns
- Surface curvature in cm^-2 (negative = concave, positive = convex)

For EPSG:4326, converts degrees to meters based on latitude.
For other projections, assumes units are already in meters.
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


"""
    slope(dhdx, dhdy, epsg; lat) -> (dhdx_scaled, dhdy_scaled)

Calculate surface slope components, accounting for coordinate projection.

# Arguments
- `dhdx`: Height derivative in x-direction
- `dhdy`: Height derivative in y-direction
- `epsg`: EPSG code of coordinate system
- `lat`: Latitude value(s) for geographic conversion

# Returns
- Tuple of scaled slope components (dhdx_scaled, dhdy_scaled)

For EPSG:4326, converts degrees to meters based on latitude.
For other projections, assumes units are already in meters.
"""
function slope(dhdx, dhdy, epsg; lat)

    if epsg == 4326
        ydist, xdist = meters2lonlat_distance.(1, lat) # degree per meter
    else
        ydist = 1
        xdist = 1
    end

    dhdx = dhdx * xdist
    dhdy = dhdy * ydist
   
    return dhdx, dhdy
end


"""
    vector_overlap(a, b) -> Bool

Check if two vectors have overlapping ranges.

# Arguments
- `a`: First vector
- `b`: Second vector

# Returns
- `Bool`: True if the ranges overlap, false otherwise
"""
function vector_overlap(a, b)
    ea = extrema(a)
    eb = extrema(b)

    tf = (ea[2] >= eb[1]) && (ea[1] <= eb[2])
    return tf
end


"""
    nanmean(x) -> Number

Calculate the mean of non-NaN values in a collection.

# Arguments
- `x`: Collection of values that may contain NaN elements

# Returns
- Mean of all non-NaN values, or NaN if all values are NaN
"""
function nanmean(x)
    mean(x[.!isnan.(x)])
end

"""
    nanmedian(x) -> Number

Calculate the median of non-NaN values in a collection.

# Arguments
- `x`: Collection of values that may contain NaN elements

# Returns
- Median of all non-NaN values, or NaN if all values are NaN
"""
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

"""
    datenum2date(n) -> DateTime

Convert MATLAB datenum format to Julia DateTime.

# Arguments
- `n`: MATLAB datenum value (days since 0000-01-00)

# Returns
- `DateTime`: Equivalent Julia DateTime object
"""
function datenum2date(n) 
    d = MATLAB_EPOCH + Dates.Millisecond(round(Int64, n * 1000 * 60 * 60 * 24))
    return d
end

"""
    add_dem_ref!(altim, dem_id, geotile, mission_geotile_folder) -> DataFrame

Add digital elevation model (DEM) reference information to an altimetry DataFrame.

# Arguments
- `altim::DataFrame`: Altimetry data to which DEM information will be added
- `dem_id::Symbol`: DEM identifier (`:best`, `:cop30_v2`, `:arcticdem_v4_10m`, or `:rema_v2_10m`)
- `geotile`: Geotile object containing extent information
- `mission_geotile_folder::String`: Path to folder containing geotile DEM files

# Returns
- `DataFrame`: The modified altimetry DataFrame with added columns:
  - `height_ref`: Reference height from DEM
  - `curvature`: Surface curvature
  - `dh`: Height difference (altim.height - height_ref)

When `dem_id` is `:best`, tries multiple DEMs in order of preference.
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
            curv = curvature.(dem.dhddx, dem.dhddy, Ref(dem_info(dem_id1)[1].epsg), lat=mean(geotile.extent.Y))

            ind = .!isnan.(curv) .& (abs.(dem.height) .< 9998)

            try
                altim[ind, :height_ref] = dem[ind, :height]
                altim[ind, :curvature] = curv[ind]
            catch e
                @warn "Error adding dem info to altim dataframe, it is most likely that the dem is out of data: $path2dem"
                throw(e)
            end
        end
    end

    altim[!, :dh] = altim.height .- altim.height_ref
    return altim
end


"""
    highres_mask(geotile, feature, invert, excludefeature) -> (mask, area_m2)

Create a high-resolution binary mask from vector features for a geotile.

# Arguments
- `geotile`: Geotile object with extent information
- `feature`: Vector feature to rasterize into a mask
- `invert`: Boolean indicating whether to invert the mask
- `excludefeature`: Optional feature to exclude from the mask (can be nothing)

# Returns
- `mask`: Binary raster mask at ~30m resolution
- `area_m2`: Matrix of cell areas in square meters
"""
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
    d = meters2lonlat_distance.(Ref(1), lat)
    a = abs.((1 ./ getindex.(d, 2) * (lat[2] .- lat[1])) .* (1 / d[1][1] * (lon[2] - lon[1])))
    area_m2 = repeat(a', outer=[length(lon), 1])

    return (mask1, area_m2)
end


"""
    highres_mask!(masks0, altim, geotile, feature, invert, excludefeature, surface_mask) -> DataFrame

Apply a high-resolution binary mask to points within a geotile and update a DataFrame column.

# Arguments
- `masks0::DataFrame`: DataFrame to update with mask values
- `altim`: Object containing longitude and latitude coordinates
- `geotile`: Geotile object with extent information
- `feature`: Vector feature to rasterize into a mask
- `invert::Bool`: Whether to invert the mask
- `excludefeature`: Optional feature to exclude from the mask (can be nothing)
- `surface_mask::Symbol`: Column name in masks0 to update with mask values

# Returns
- Updated DataFrame with the surface_mask column populated
"""
function highres_mask!(masks0, altim, geotile, feature, invert, excludefeature, surface_mask)
    mask1, _ = highres_mask(geotile, feature, invert, excludefeature)
    
    grid_resolution = 0.00027 # ~30m

    x_mask = X(geotile.extent.X[1]:grid_resolution:geotile.extent.X[2],
        sampling=DimensionalData.Intervals(DimensionalData.Start()))
    y_mask = Y(geotile.extent.Y[1]:grid_resolution:geotile.extent.Y[2],
        sampling=DimensionalData.Intervals(DimensionalData.Start()))

    valid = within.(Ref(geotile.extent), altim.longitude, altim.latitude)
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


"""
    binningfun_define(binning_method) -> Function

Create a binning function based on the specified method.

# Arguments
- `binning_method::String`: Method to use for binning data. Options:
  - "meanmadnorm3": Mean of values with MAD normalization < 3
  - "meanmadnorm5": Mean of values with MAD normalization < 5
  - "meanmadnorm10": Mean of values with MAD normalization < 10
  - "median": Median of all values

# Returns
- Function that implements the specified binning method
"""
function binningfun_define(binning_method)
    if binning_method == "meanmadnorm3"
        x -> mean(x[madnorm(x).<3])
    elseif binning_method == "meanmadnorm5"
        x -> mean(x[madnorm(x).<5])
    elseif binning_method == "meanmadnorm10"
        x -> mean(x[madnorm(x).<10])
    elseif binning_method == "median"
        x -> median(x)
    else
        error("unrecognized binning method")
    end
end
"""
    highres_mask(latitude, longitude, feature; grid_resolution=0.00027) -> Vector{Bool}

Create a binary mask for points based on their intersection with a high-resolution vector feature.

# Arguments
- `latitude::Vector{<:Real}`: Vector of latitude coordinates
- `longitude::Vector{<:Real}`: Vector of longitude coordinates
- `feature`: Vector feature (polygon) used to create the mask

# Keywords
- `grid_resolution=0.00027`: Resolution of the rasterized grid in degrees

# Returns
- Vector of boolean values indicating whether each point intersects with the feature
"""
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


"""
    rginum2label(rginum)

Converts an RGI (Randolph Glacier Inventory) region number to its corresponding label.
Takes a numeric RGI region identifier and returns the human-readable label by first
converting the number to a text key using `rginum2txt` and then looking up the label
in the `rgi2label` dictionary.
"""
function rginum2label(rginum) 
    return getindex(rgi2label,getindex(rginum2txt, rginum))
end



rginum2txt = Dict(
    1 => "rgi1" ,
    2 => "rgi2" ,
    3 => "rgi3" ,
    4 => "rgi4" ,
    5 => "rgi5",
    6 => "rgi6" ,
    7 => "rgi7" ,
    8 => "rgi8" ,
    9 => "rgi9" ,
    10 => "rgi10" ,
    11 => "rgi11" ,
    12 => "rgi12" ,
    13 => "rgi13" ,
    14 => "rgi14" ,
    15 => "rgi15" ,
    16 => "rgi16" ,
    17 => "rgi17" ,
    18 => "rgi18" ,
    19 => "rgi19" ,
    98 => "hma" ,
    99 => "global" ,
)

rginum2enclosed_alphanumerics = Dict(
    1 => "①" ,
    2 => "②" ,
    3 => "③" ,
    4 => "④" ,
    5 => "⑤",
    6 => "⑥" ,
    7 => "⑦" ,
    8 => "⑧" ,
    9 => "⑨" ,
    10 => "⑩" ,
    11 => "⑪" ,
    12 => "⑫" ,
    13 => "⑬" ,
    14 => "⑭" ,
    15 => "⑮" ,
    16 => "⑯" ,
    17 => "⑰" ,
    18 => "⑱" ,
    19 => "⑲" ,
    98 => "Ⓗ",
    99 => "Ⓖ",
)


"""
    ntpermutations(nt_in) -> Vector{NamedTuple}

Generate all possible permutations of values from a named tuple of iterables.

# Arguments
- `nt_in`: A named tuple where each field contains an iterable

# Returns
- Vector of named tuples containing all possible combinations of the input values,
  preserving the original field names
"""
ntpermutations(nt_in) = [NamedTuple{keys(nt_in)}(reverse(t)) for t in Iterators.product(reverse(nt_in)...)] |> vec


"""
    extent2rectangle(extent) -> GeoInterface.Wrappers.Polygon

Convert an extent tuple to a rectangular polygon.

# Arguments
- `extent`: Tuple containing x-bounds and y-bounds as ((xmin, xmax), (ymin, ymax))

# Returns
- A GeoInterface Polygon representing the rectangular boundary of the extent
"""
function extent2rectangle(extent)
    xbounds = extent[1]
    ybounds = extent[2]
    rectangle = GeoInterface.Wrappers.Polygon([[(xbounds[1], ybounds[1]), (xbounds[1], ybounds[2]), (xbounds[2], ybounds[2]), (xbounds[2], ybounds[1]), (xbounds[1], ybounds[1])]])
    return rectangle
end




"""
    df_tsfit!(df, tsvars; progress=true, datelimits=nothing) -> DataFrame

Fit a least squares model to time series variables in a DataFrame, calculating offset, trend, amplitude, and phase.

# Arguments
- `df`: DataFrame containing time series data
- `tsvars`: Column names of time series variables to fit

# Keywords
- `progress=true`: Whether to display a progress bar
- `datelimits=nothing`: Optional tuple of (start_date, end_date) to limit the date range for fitting

# Returns
- Modified DataFrame with added columns for offset, trend, amplitude, and phase for each time series variable
"""
function df_tsfit!(df, tsvars; progress=true, datelimits = nothing)
    prog = progress ? ProgressMeter.Progress(length(tsvars); desc="Fitting LSQ model to timeseries...") : nothing

    for tsvar in tsvars
        ddate = colmetadata(df, tsvar, "date")
        x = decimalyear.(ddate)
        x = x .- mean(x)

        out_var_offset = string(tsvar)*"_offset"
        out_var_trend = string(tsvar)*"_trend"
        out_var_amp = string(tsvar)*"_amplitude"
        out_var_phase = string(tsvar)*"_phase"
        df[!, out_var_offset] = zeros(nrow(df))

        df[!, out_var_trend] = zeros(nrow(df))
        df[!, out_var_amp] = zeros(nrow(df))
        df[!, out_var_phase] = zeros(nrow(df))

        if !isnothing(datelimits)
            date_index = (ddate .> datelimits[1]) .& (ddate .< datelimits[2])
        else
            date_index = trues(length(ddate))
        end

        Threads.@threads for g in eachrow(df)
            y = g[tsvar];
            valid = .!isnan.(y) .& date_index

            fit = curve_fit(
                offset_trend_seasonal2, 
                x[valid], 
                y[valid], 
                p_offset_trend_seasonal
                )

            g[out_var_offset] = fit.param[1]
            g[out_var_trend] = fit.param[2]
            g[out_var_amp] = hypot(fit.param[3], fit.param[4])
            g[out_var_phase] = 365.25 * (mod(0.25 - atan(fit.param[4], fit.param[3]) / (2π), 1))
        end

        # Update the progress meter
        progress && update!(prog, prog.counter .+ 1)
    end
    return df
end


"""
    validgaps(valid) -> BitVector

Identify gaps within a valid data range.

# Arguments
- `valid`: BitVector or Boolean array indicating valid data points

# Returns
- BitVector marking gaps (false values) that occur between the first and last valid points

# Details
Returns a BitVector where `true` indicates invalid points (gaps) that occur between the first 
and last valid points in the input array. Points outside this range are always marked as `false`.
"""
function validgaps(valid)
    validgap = falses(length(valid))
    if any(valid) && !all(valid)
        sind = findfirst(valid)
        eind = findlast(valid)
        validgap[sind:eind] = .!valid[sind:eind]
    end
    return validgap
end

"""
    connected_groups(connectivity) -> Vector{Int64}

Group elements into connected components based on their connectivity relationships.

# Arguments
- `connectivity`: Vector of vectors where each element contains indices of connected elements

# Returns
- Vector of integers assigning each element to a group ID

# Details
Identifies connected components in a graph represented by the connectivity list.
Elements with empty connectivity lists are assigned to their own individual groups.
"""
function connected_groups(connectivity)

    # Track which tiles have been assigned to groups
    assigned = falses(length(connectivity))
    group_id = zeros(Int64, length(connectivity))

    # Initialize group counter
    group_id0 = 1

    while !all(assigned)
        # Find first unassigned tile with intersections
        start_idx = findfirst(.!assigned)

        if isempty(connectivity[start_idx])
            assigned[start_idx] = true
            group_id[start_idx] = group_id0
            group_id0 += 1
            continue
        end

        # Build group starting from this tile
        current_group = Int64[start_idx]
        assigned[start_idx] = true
        group_id[start_idx] = group_id0

        # Keep expanding group until no more intersections found
        expanded = true
        while expanded
            expanded = false

            # Check all tiles in current group
            for tile_idx in current_group
                # Look at all intersecting tiles
                for intersect_idx in connectivity[tile_idx]
                    if !assigned[intersect_idx]
                        push!(current_group, intersect_idx)
                        assigned[intersect_idx] = true
                        group_id[intersect_idx] = group_id0
                        expanded = true
                    end
                end
            end
        end

        group_id0 += 1
    end

    return group_id
end




"""
    nc2dd(ds, varname)

Convert a NetCDF variable to a DimensionalData.DimArray.

# Arguments
- `ds`: NCDataset containing the variable
- `varname`: String name of the variable to convert

# Returns
- DimArray with dimensions and data from the NetCDF variable

# Notes
- Preserves dimension names and coordinates from the NetCDF dataset
- Converts dimension names to symbols for DimensionalData compatibility
"""
function nc2dd(cfv)
    dnames = NCDatasets.dimnames(cfv)

    if in("units", keys(cfv.attrib)) 
        unit = cfv.attrib["units"]
        try
            units = Unitful.uparse(unit; unit_context=Unitful.unitmodules)
        catch
            units = unit_convert[unit]
        end
    else
        units = NoUnits
    end

    da = DimArray(collect(cfv[:,:]) * units, tuple([Dim{Symbol(dname)}(cfv[dname]) for dname in dnames]...))
    return da
end


"""
    update_geotile_path(paths; mission = :hugonnet, path_replace ="/2deg" => "/2deg_unfiltered")

Updates the geotile path for a specific mission in the paths NamedTuple.

# Arguments
- `paths`: NamedTuple containing mission paths
- `mission`: Symbol specifying which mission to update (default: `:hugonnet`)
- `path_replace`: Pair of strings for search and replace in the path (default: `"/2deg" => "/2deg_unfiltered"`)

# Returns
- Updated paths NamedTuple with the modified geotile path
"""
function update_geotile_path(paths; mission = :hugonnet, path_replace ="/2deg" => "/2deg_unfiltered")
    k = keys(paths[mission])
    v = collect(values(paths[mission]))

    ind = findfirst(k .== :geotile)
    v[ind] = replace(v[ind], path_replace)

    nt = NamedTuple{k}(v)

    k = keys(paths); 
    v = collect(values(paths))

    ind = findfirst(k .== :hugonnet)
    v[ind] = nt

    paths = NamedTuple{k}(v)

    return paths
end


"""
    wimberly2024(; filename = pathlocal[:wimberly_2024])

Loads and processes glacier model data from Wimberly et al. 2024.

# Arguments
- `filename`: Path to the CSV file containing the data (default: pathlocal[:wimberly_2024])

# Returns
- A DimArray containing glacier model runoff data in units of [km^3] organized by date, 
  GCM, SSP, basin, and model type. The function parses the combined GCM_SSP_Basin column 
  into separate dimensions and structures the data from three glacier models 
  (GloGEM, PyGEM, OGGM) into a DimArray.
"""
function wimberly2024(; filename=setpaths().wimberly_2024)

    df = CSV.read(filename, DataFrame)

    # split the GCM_SSP_Basin column into three columns
    df = transform(df, [:GCM_SSP_Basin] => ByRow(x -> split(x, "_")) => [:GCM, :SSP, :Basin])

    ddate = Ti(sort!(unique(df.Date)))
    dgcm = Dim{:GCM}(unique(df.GCM))
    dssp = Dim{:SSP}(unique(df.SSP))
    dbasin = Dim{:Basin}(unique(df.Basin))
    dmodel = Dim{:Model}(["GloGEM", "PyGEM", "OGGM"])

    gdf = groupby(df, [:GCM, :SSP, :Basin])

    da = zeros(ddate, dgcm, dssp, dbasin, dmodel)

    for gcm in dgcm
        for ssp in dssp
            for basin in dbasin
                for model in dmodel
                    da[At(gdf[(gcm, ssp, basin)][:, :Date]), At(gcm), At(ssp), At(basin), At(model)] = gdf[(gcm, ssp, basin)][:, model]
                end
            end
        end
    end

    return da
end


"""
    read_grace_rgi(; path2file=setpaths()[:grace_rgi])

Read GRACE (Gravity Recovery and Climate Experiment) data for RGI (Randolph Glacier Inventory) regions.

This function loads GRACE data from a MATLAB file and converts region codes from letter-based 
abbreviations to standard RGI numerical identifiers (e.g., "ALA" to "rgi1").

# Arguments
- `path2file`: Path to the GRACE RGI data file (MATLAB format)

# Returns
- A dictionary mapping RGI region identifiers to their corresponding GRACE data

# Note
The function converts region codes from the original letter-based abbreviations (e.g., "ALA" for Alaska)
to the standard RGI numerical format (e.g., "rgi1").
"""
function read_grace_rgi(;path2file=setpaths()[:grace_rgi])

    grace0 = matread(path2file)

    # rgi leter to digit mapping
    old2new = Dict(
        "PAT" => "rgi18", 
        "GRE" => "rgi5", 
        "NAS" => "rgi10", 
        "ICE" => "rgi6", 
        "TRP" => "rgi16",
        "SAW" => "rgi14",
        "SAE" => "rgi15",
        "CDN" => "rgi3",
        "CAS" => "rgi13",
        "AND" => "rgi17",
        "CDS" => "rgi4",
        "ANT" => "rgi19",
        "CEU" => "rgi11",
        "ALA" => "rgi1",
        "SVB" => "rgi7",
        "WNA" => "rgi2", 
        "NEZ" => "rgi18",
        "RAI" => "rgi9",
        "SCA" => "rgi8"
    )


    grace = Dict(any(keys(old2new).== key) ? (old2new[key]) => val : (key) => val for (key, val) in grace0)
    return grace
end

"""
    grace_masschange(; path2file=setpaths()[:grace_rgi])

Load and organize GRACE (Gravity Recovery and Climate Experiment) mass change data for glacier regions.

This function reads GRACE data for RGI (Randolph Glacier Inventory) regions and organizes it into a
dimensional array structure with dimensions for region, date, and error bounds.

# Arguments
- `path2file`: Path to the GRACE RGI data file (MATLAB format). Defaults to the path from `setpaths()[:grace_rgi]`.

# Returns
- A DimensionalData array with dimensions for RGI regions (1-19, plus 98 for HMA and 99 for Global),
  dates, and error bounds (false=value, true=error).

# Note
The error values are multiplied by 2 to represent 95% confidence intervals.
"""
function grace_masschange(; path2file=setpaths()[:grace_rgi])
    grace_raw = read_grace_rgi(; path2file)

    # organize GRACE grace data into a DimArray
    drgi = Dim{:rgi}(collect([1:19..., 98:99...]))
    ddate = Dim{:date}(vec(datenum2date.(grace_raw["rgi1"]["dM_gt_mdl_fill_date"])))
    derror = Dim{:error}([false, true])

    grace = fill(NaN, drgi, ddate, derror)

    for rgi in drgi
        #rgi = drgi[1]
        if rgi == 98
            rgi_id = "HMA"
        elseif rgi == 99
            rgi_id = "Global"
        else
            rgi_id = "rgi$(rgi)"
        end

        if haskey(grace_raw, rgi_id)
            grace[At(rgi), :, At(false)] = vec(grace_raw[rgi_id]["dM_gt_mdl_fill"])

            ## multiply by 2 to get 95% confidence interval
            grace[At(rgi), :, At(true)] = vec(grace_raw[rgi_id]["dM_sigma_gt_mdl_fill"]) * 2
        end
    end

    return grace
end


"""
    check_for_all_nans(d) -> Nothing

Check if any key in a dictionary contains only NaN values.

# Arguments
- `d`: Dictionary to check

# Returns
- Nothing

# Throws
- Error if any key in the dictionary contains only NaN values
"""
function check_for_all_nans(d)
    for k in keys(d)
        if all(isnan.(d[k]))
            error("all NaNs for $k")
        end
    end
end


"""
    resample(colorscheme::Symbol, n) -> ColorSchemes.ColorScheme

Create a new color scheme by resampling an existing one to a specified number of colors.

# Arguments
- `colorscheme::Symbol`: Name of the source color scheme
- `n`: Number of colors in the new scheme

# Returns
- A new ColorScheme with n evenly spaced colors from the original scheme
"""
function resample(colorscheme::Symbol, n)
    newscheme = ColorSchemes.ColorScheme(
        get(ColorSchemes.colorschemes[colorscheme], (0:n-1) / (n-1))
    )
    return newscheme
end


"""
    glambie2024(; path2glambie=paths = setpaths()[:glambie_2024]) -> DimArray

Load and process the GlaMBIE 2024 glacier mass balance dataset.

# Keywords
- `path2glambie`: Path to the GlaMBIE 2024 CSV file (default: from setpaths())

# Returns
- A DimensionalData.DimArray containing cumulative glacier mass balance in Gt
  with dimensions for RGI region (1-19, 98, 99), date (2000-2024), and error flag

# Notes
- Dates are shifted by up to 6 months from original data to align across hemispheres
- Missing values are linearly interpolated
- Region 98 (HMA) is calculated as the sum of regions 13-15
"""
function glambie2024(; path2glambie=paths = setpaths()[:glambie_2024])
    df = CSV.read(path2glambie, DataFrame)

    drgi = Dim{:rgi}(vcat(collect(1:19), 98, 99))
    derror = Dim{:error}([false, true])
    @warn "GlaMBIE dates are shifted by up to 6 months from the original data to align dates across hemispheres"
    ddate = Dim{:date}(collect(decimalyear2datetime.(2000:2024)); name="GlaMBIE 2024 cumulative glacier mass balance [GLAMBIE]")

    glambie_Gt = fill(NaN, drgi, ddate, derror)

    for (i, rgi) in enumerate(drgi)
        if rgi == 98
            continue
        elseif rgi == 99
            i = i-1
        end
        scol = (i - 1) * 3 + 1
        date1 = round.(df[:, scol])
        valid = .!ismissing.(date1)
        date1 = collect(decimalyear2datetime.(date1[valid]))

        glambie_Gt[At(rgi), At(date1), At(false)] = coalesce.(df[valid, scol+1], NaN)
        glambie_Gt[At(rgi), At(date1), At(true)] = coalesce.(df[valid, scol+2], NaN)
    end

 
    glambie_Gt[At(98), :, :] .= 0
    for rgi in 13:15
        glambie_Gt[At(98), :, At(false)] .+= glambie_Gt[At(rgi), :, At(false)]
        glambie_Gt[At(98), :, At(true)] .+= glambie_Gt[At(rgi), :, At(true)].^2
    end

    glambie_Gt[At(98), :, At(true)] = sqrt.(glambie_Gt[At(98), :, At(true)])

    @warn "GlaMBIE data is linearly interpolated over missing date values"
    for rgi in drgi
        for err_flag in derror
            valid = .!isnan.(glambie_Gt[At(rgi), :, At(err_flag)])
            if !all(valid)
                fit = curve_fit(model1_trend, collect(decimalyear.(ddate[valid])), collect(glambie_Gt[At(rgi), :, At(err_flag)][valid]), p_trend)
                glambie_Gt[At(rgi), .!valid, At(err_flag)] = model1_trend(decimalyear.(collect(ddate[.!valid])), fit.param)
            end
        end
    end
    
    return (glambie_Gt)
end
 

"""
    dimarray2netcdf(da, filename; name="var", units=nothing, global_attributes=nothing) -> String

Convert a DimensionalData.DimArray to a NetCDF file.

# Arguments
- `da`: DimensionalData.DimArray to convert
- `filename`: Path where the NetCDF file will be saved

# Keywords
- `name="var"`: Name for the main variable in the NetCDF file
- `units=nothing`: Units for the main variable (uses DimArray units if not specified)
- `global_attributes=nothing`: Dictionary of global attributes to add to the NetCDF file

# Returns
- Path to the created NetCDF file
"""
function dimarray2netcdf(da, filename; name = "var", units = nothing, global_attributes = nothing)
    
    NCDataset(filename, "c") do ds
        # get dimensions
        da_dims = dims(da)

        # step 1: define dimensions
        for dim in da_dims
            dname = string(DimensionalData.name(dim));
            defDim(ds, dname, length(dim))
        end

        # step 2: add dim variables & metadata
        for dim in da_dims
            dname = string(DimensionalData.name(dim));
            d = defVar(ds, dname, val(val(dim)), (dname,))

            if eltype(dim) <: Number
                d.attrib["units"] = string(Unitful.unit(dim[1]))
            end

            # add metadata
            for (k,v) in DD.metadata(dim)
                d.attrib[k] = v
            end
        end

        # step 3: add variable
        v = defVar(ds, name, ustrip.(parent(da)), string.(DimensionalData.name.(da_dims)))
        
        # step 4: add variable metadata
        if units == nothing
            if eltype(da) <: Number
                v.attrib["units"] = string(Unitful.unit(da[1]))
            end
        else
            v.attrib["units"] = string(units)
        end

        for (k,v) in DD.metadata(da)
            d.attrib[k] = v
        end

        # step 5: add global attributes
        if global_attributes != nothing
            for (k,v) in global_attributes
                ds.attrib[k] = v
            end
        end

        return filename
    end;
end


"""
    netcdf2dimarray(filename; varname=nothing) -> DimensionalData.DimArray

Convert a NetCDF file to a DimensionalData.DimArray.

# Arguments
- `filename`: Path to the NetCDF file

# Keywords
- `varname=nothing`: Name of the variable to convert. If not specified, attempts to
  automatically detect the single non-dimension variable in the file.

# Returns
- DimArray with dimensions and data from the NetCDF variable, preserving units and metadata

# Throws
- Error if varname is not specified and there isn't exactly one variable in the file
"""
function netcdf2dimarray(filename; varname=nothing)
    # read the netcdf file
    ds = NCDatasets.Dataset(filename) 

    # get the dimension names
    dnames = NCDatasets.dimnames(ds)

    # get the variable name
    if varname === nothing
        varname = setdiff(keys(ds), dnames)
        if length(varname) != 1
            error("varname must be specified if there is not exactly one variable in the netcdf file")
        end
        varname = first(varname)
    end

    # get the units of variable
    da = _var_from_netcdf(ds, varname, dnames)
    
    return da
end


"""
    _var_from_netcdf(ds, varname, dnames) -> DimensionalData.DimArray

Create a DimArray from a NetCDF variable.

# Arguments
- `ds`: NCDataset containing the variable
- `varname`: Name of the variable to convert
- `dnames`: Names of dimensions in the dataset

# Returns
- DimArray with appropriate dimensions and units (if available)
"""
function _var_from_netcdf(ds, varname, dnames)
    units = _units_from_netcdf(ds, varname)
    if !isnothing(units)
        da = DimArray(ds[varname] * units, tuple(_dim_from_netcdf.(Ref(ds), dnames)...); metadata=ds[varname].attrib)
    else
        da = DimArray(ds[varname], tuple(_dim_from_netcdf.(Ref(ds), dnames)...); metadata=ds[varname].attrib)
    end
    return da
end

"""
    _dim_from_netcdf(ds, dname) -> DimensionalData.Dimension

Create a dimension from a NetCDF dimension variable.

# Arguments
- `ds`: NCDataset containing the dimension
- `dname`: Name of the dimension to convert

# Returns
- Dimension with values and units (if available) and original metadata
"""
function _dim_from_netcdf(ds, dname)
    units = _units_from_netcdf(ds, dname)
    if !isnothing(units)
        dim = Dim{Symbol(dname)}(ds[dname][:]*units; metadata=ds[dname].attrib)
    else
        dim = Dim{Symbol(dname)}(ds[dname][:]; metadata=ds[dname].attrib)
    end
    return dim
end

"""
    _units_from_netcdf(ds, varname) -> Union{Unitful.Units, Nothing}

Extract and parse units from a NetCDF variable.

# Arguments
- `ds`: NCDataset containing the variable
- `varname`: Name of the variable to extract units from

# Returns
- Parsed Unitful units if available and parsable
- `NoUnits` if units attribute exists but is empty
- `nothing` if units are not available or cannot be parsed
"""
function _units_from_netcdf(ds, varname)

    # units can not be applied to non-numeric data
    if !(eltype(ds[varname]) <: Number)
        units = nothing
    else
        if haskey(ds[varname].attrib, "units")
            if isempty(ds[varname].attrib["units"])
                units = NoUnits
            else
                units = ds[varname].attrib["units"];
                if isempty(units)
                    units = nothing
                else
                    try
                        units = uparse(units; unit_context=Unitful.unitmodules)
                    catch
                        @warn "could not parse units for $varname"
                        units = nothing
                    end
                end
            end
        else 
            units = nothing
        end
    end
    return units
end

unit_convert = Dict()
unit_convert["m3 s^-1"] = u"m^3/s"
unit_convert["kg s^-1"] = u"kg/s"