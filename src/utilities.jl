# utilities for working building geogrid point database
pathlocal = setpaths()

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


"""
    geotile_extent(geotile_id::String) -> Extent

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

    return Extent(X=(min_x, max_x), Y=(min_y, max_y))
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
    nmad(x)

Calculate the normalized median absolute deviation (MAD) of values in x.

# Arguments
- `x`: Vector of values to analyze

# Returns
- Vector of normalized deviations, where each value represents the number of 
  standard deviations from the median (using MAD scaled by 1.4826)
"""
function nmad(x)
    consistent_estimator = 1.4826 #mad to sigma conversion factor
    x_abs = abs.(x .- median(x))
    x_nmad = x_abs ./ (median(x_abs) .* consistent_estimator)
    return x_nmad
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

        Threads.@threads for mask in masks
            sym = Symbol("$(mask)_shp");
            feature = Shapefile.Handle(pathlocal[sym])

            sym = Symbol("$(mask)_frac")
            geotiles[!,sym] = zeros(nrow(geotiles))

            # Threads causes issues 
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


"""
    mission_proper_name(mission) -> String

Convert mission identifier to proper display name.

# Arguments
- `mission`: String identifier for the mission (e.g., "icesat2", "gedi")

# Returns
- Properly formatted mission name for display purposes

# Description
Maps internal mission identifiers to their proper display names used in publications
and user interfaces. Returns the original string if no mapping is defined.
"""
function mission_proper_name(mission)
    if mission == "icesat2"
        return "ICESat-2"
    elseif mission == "icesat"
        return "ICESat"
    elseif mission == "gedi"
        return "GEDI"
    elseif mission == "hugonnet"
        return "ASTER/WorldView"
    else
        return mission
    end
end


"""
    overlap_range(x1, x2) -> (range_overlap1, range_overlap2)

Find the overlapping range indices between two arrays.

# Arguments
- `x1`: First array of values
- `x2`: Second array of values

# Returns
- `range_overlap1`: Range of indices in x1 that overlap with x2
- `range_overlap2`: Range of indices in x2 that overlap with x1

# Description
Calculates the minimum and maximum values that exist in both arrays,
then returns the index ranges for each array that fall within this
overlapping range. Useful for finding common data ranges between
two datasets.
"""
function overlap_range(x1, x2)
    x_min = max(minimum(x1), minimum(x2))
    x_max = min(maximum(x1), maximum(x2))

    ind1 = (x1 .>= x_min) .& (x1 .<= x_max)
    range_overlap1 = findfirst(ind1):findlast(ind1)

    ind2 = (x2 .>= x_min) .& (x2 .<= x_max)
    range_overlap2 = findfirst(ind2):findlast(ind2)

    return range_overlap1, range_overlap2
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
    geotile_hypsometry(geotiles, surface_mask; dem_id=:cop30_v2, force_remake_before=nothing)

Calculate glacier hypsometry for geotiles and save to file.

# Arguments
- `geotiles`: DataFrame containing geotile information
- `surface_mask`: Surface mask identifier ("glacier" or "glacier_rgi7")
- `dem_id`: DEM identifier (default: :cop30_v2)
- `force_remake_before`: DateTime threshold for forcing recalculation

# Returns
GeoDataFrame containing glacier hypsometry data

# Description
Loads existing hypsometry file if available and recent enough, otherwise calculates
area-elevation distributions for glaciers within geotiles using DEM data.
"""
function geotile_hypsometry(geotiles, surface_mask; dem_id=:cop30_v2, force_remake_before=nothing)

    path2surface_mask = pathlocal[Symbol("$(surface_mask)_individual")]
    geotile_hyps_file = replace(path2surface_mask, "." * split(path2surface_mask, ".")[end] => "geotile_hyps.jld2")

    # Only process if output file doesn't exist or force_remake is true
    if isfile(geotile_hyps_file) && isnothing(force_remake_before)
        surface_mask_geom = FileIO.load(geotile_hyps_file, surface_mask)
    elseif isfile(geotile_hyps_file) && Dates.unix2datetime(mtime(geotile_hyps_file)) > force_remake_before
        surface_mask_geom = FileIO.load(geotile_hyps_file, surface_mask)
    else

        # Set the glacier ID attribute name based on RGI version
        if Base.contains(path2surface_mask, "rgi6")
            persistent_attribute = :RGIId
        elseif Base.contains(path2surface_mask, "rgi7")
            persistent_attribute = :rgi_id
        else
            error("surface_mask must be either \"glacier\" or \"glacier_rgi7\"")
        end

        # Read glacier outlines
        surface_mask_geom = GeoDataFrames.read(path2surface_mask)

        # Load elevation data and define height bins
        h = Raster(pathlocal[dem_id], lazy=true)
        height_range, height_center = project_height_bins()

        # Calculate hypsometry for each geotile combination
        surface_mask_geom = geotile_zonal_area_hyps(h, height_range, surface_mask_geom, geotiles.id; persistent_attribute)

        # Save results
        FileIO.save(geotile_hyps_file, Dict(surface_mask => surface_mask_geom))
    end
    return surface_mask_geom
end


function _trim_data(data::DimArray, geotile_ids)
    return data[geotile = At(geotile_ids)]
end

function _trim_data(data::DataFrame, geotile_ids)

    if in("id",names(data))
        index, _ = intersectindices(data.id, geotile_ids)   
    elseif in("geotile",names(data))
        index, _ = intersectindices(data.geotile, geotile_ids)
    else
        error("data must have a column named 'id' or 'geotile'")
    end
    return data[index, :]
end

function _trim_data(f::Dict{String, Any}, geotile_ids)
    dgeotile = dims(f["dh_hyps"]["hugonnet"], :geotile)
    if length(dgeotile) != length(geotile_ids)
        for k1 in keys(f)
        #k1 = first(keys(f))
            for k2 in keys(f[k1])
            #k2 = first(keys(f[k1]))
                f[k1][k2] = _trim_data(f[k1][k2], geotile_ids)
            end
        end
    end
    return f
end

function trim_files(path2runs, geotile_ids)
    @showprogress desc="Trimming data to only include geotiles with glacier coverage... !! WARNING DATA WILL BE OVERWRITTEN !!"  for path2file in path2runs
        f = load(path2file)::Dict{String, Any}
        dgeotile = dims(f["dh_hyps"]["hugonnet"], :geotile)

        if length(dgeotile) != length(geotile_ids)
            @warn "Geotile dimensions do not match, removing excess geotiles: $(path2file)"
            f = _trim_data(f, geotile_ids) 

            save(path2file, f)
        end
    end
end


"""
    _geotile_load_align(; surface_mask="glacier", geotile_width=2, geotile_order=nothing)

Load and align geotiles for a given surface mask.

# Arguments
- `surface_mask`: Surface mask type (default: "glacier")
- `geotile_width`: Width of geotiles in degrees (default: 2)
- `geotile_order`: Optional geotile order for alignment (default: nothing)

# Returns
Tuple of (geotiles0, reg) where geotiles0 is the processed geotiles DataFrame and reg contains region information
"""
function _geotile_load_align(; surface_mask="glacier", geotile_width=2, geotile_order=nothing)

    # Generate geotiles for this mask
    geotiles0 = geotiles_mask_hyps(surface_mask, geotile_width)

    # Make geotiles mutually exclusive by RGI region
    geotiles0, reg = geotiles_mutually_exclusive_rgi!(copy(geotiles0))

    # Rename area column to include surface mask
    #rename!(geotiles0,"$(surface_mask)_area_km2" => "area_km2")

    if !isnothing(geotile_order)
        # Align geotile indices with elevation change dimensions
        gt_index, dgt_index = intersectindices(geotiles0.id, geotile_order)

        geotiles0 = geotiles0[gt_index, :]
    end
end

