# Standard libraries
using Random, LinearAlgebra, SparseArrays
# Geospatial libraries
import GeoInterface as GI, GeometryOps as GO, LibGEOS as LG
import Proj, ArchGDAL
import GeoFormatTypes: EPSG
# File reading
using Shapefile, DataFrames, GeoDataFrames
using Rasters # and raster processing
# Spatially binning the geometries
using SortTileRecursiveTree, AbstractTrees


# You can cache a Raster using an LRU cache like this,
# but in this case it doesn't help performance much because the chunks are too small.
# cached_dem = Rasters.rebuild(dem; data = DiskArrays.cache(parent(dem); maxsize = 3000))

"""
    minimum_elevation(masked_dem)

Find the minimum elevation in a masked DEM, and return a 2-tuple of the point at which it occurs `(long, lat)`.
"""
function minimum_elevation(masked_dem)
    # First, check if the masked DEM is empty
    if isempty(masked_dem)
        ext = GI.extent(masked_dem)
        return (ext.X[2] - ext.X[1], ext.Y[2] - ext.Y[1])
    end
    # Otherwise, construct an iterator that skips missing values
    value_iterator = skipmissing(masked_dem)
    if isnothing(iterate(value_iterator)) # check that the iterator is not empty
        ext = GI.extent(masked_dem)
        return (ext.X[2] - ext.X[1], ext.Y[2] - ext.Y[1])
    end
    # finally, find the minimum value and its index
    # this works only because Rasters.jl has a custom skipmissing that allows for this
    minval, minidx = findmin(value_iterator)
    # return the point at which the minimum value occurs
    return DimPoints(masked_dem)[minidx]
end

"""
    is_extent_suitable(raster, node; max_n_chunks = 100^2, chunk_size = 128^2)

Check if the extent of a node is suitable for loading into memory.
"""
function is_extent_suitable(raster, node; max_n_chunks=100^2, chunk_size=128^2)
    # check if the extent of this particular node is less than 100x100 chunks.
    # To do so, we convert the extent to indices within the raster,
    # and then check if the product of the lengths of the indices is less than `max_n_chunks * chunk_size`.
    # TODO: in future, maybe simply `crop` -> `DiskArrays.eachchunk(parent(cropped))` -> count chunks?
    # That seems like a cleaner approach...
    indices = Rasters.DD.dims2indices(raster, GI.extent(node))
    return prod(length.(indices)) < max_n_chunks * chunk_size
end
# NOTE: 
# I tried tuning the max number of chunks here,
# and 100x100 seemed like a sweet spot.
# 200x200 almost doubled the runtime for some reason,
# so I presume there's a balance between memory usage / download speed / decoding speed
# and runtime.
# Plus, since glaciers are nicely clustered in space,
# we can afford to do this relatively easily.

"""
    extract_idxs(node)

Extract the indices of the geometries that an `STRNode` covers.

Returns a vector of integers that are indices of the geometries that made the tree.
"""
function extract_idxs(node)::Vector{Int}
    return mapreduce(Base.Fix2(getproperty, :indices), vcat, AbstractTrees.Leaves(node))
end


function get_minimum_elevations(dem, geometries; progress=true)
    prog = progress ? ProgressMeter.Progress(length(geometries); desc="Getting minimum elevations") : nothing

    # Create a (default) STRtree with the geometries
    # This usually means a leaf node size of around 10 geometries,
    # which is not too onerous.
    tree = STRtree(geometries)

    # Preallocate the minimum elevations array so we can assign to it
    minimum_elevations = fill((NaN, NaN), length(geometries))
    # Now, we construct an iterator over the tree (using AbstractTrees.jl),
    # performing a depth first search over the tree from parent to child.
    # We also use the predicate to stop evaluating once we find a node that is suitable in extent.

    # Note that the predicate causes iteration to assume a leaf node if it's **`false`**,
    # otherwise it will continue iterating over children.
    # So, the predicate is that if the extent is unsuitable, we keep iterating (true),
    # but if it's suitable, we stop iterating over children (false).
    iterator = AbstractTrees.PreOrderDFS(
        Base.Fix1(!is_extent_suitable, dem), # predicate to stop the iterator
        tree.rootnode                        # the root node of the tree
    )

    # Loop over the tree in pre-order depth-first search
    for node::Union{SortTileRecursiveTree.STRNode,SortTileRecursiveTree.STRLeafNode} in iterator
        # We perform rasterization if the node's extent is sufficiently small,
        # or if the node is a leaf node.
        loaded_dem = if is_extent_suitable(dem, node)
            @debug("Loading DEM for node with extent $(GI.extent(node))")
            # Here, we know the node is smaller than the maximum size we're willing to load,
            # so we load the section of the DEM covering the node into memory,
            # and use that for the zonal operation.
            Rasters.read(Rasters.crop(dem; to=GI.extent(node)::GI.Extents.Extent)::Raster)::Raster
        elseif node isa SortTileRecursiveTree.STRLeafNode
            @debug("Using full DEM for node with extent $(GI.extent(node))")
            # Here, we know the node is a leaf node, and that the extent of the node is larger than
            # the maximum size we're willing to load, 
            # so we just use the full (lazy) DEM for the zonal operation.
            dem
        else
            # The node is neither suitable nor a leaf node,
            # so we skip it.
            @debug("Skipping node with extent $(GI.extent(node))")
            continue
        end
        @debug("Zonalling")
        # Extract the indices of the geometries that this node covers
        indices = extract_idxs(node)
        # Perform the zonal operation
        minimum_elevations[indices] .= Rasters.zonal(
            minimum_elevation, loaded_dem;
            of=view(geometries, indices),
            skipmissing=false,  # provide a full masked Raster to the zonal operation
            threaded=true,      # thread over the geometries
            progress=false,     # don't show progress since we are already showing global progress
        )::Vector{Union{Missing,Tuple{Float64,Float64}}}
        # Update the progress meter
        progress && update!(prog, prog.counter + length(indices))
    end

    return minimum_elevations
end


struct UnitCartesianFromGeographic <: CoordinateTransformations.Transformation
end

# median time 19 ns
function (::UnitCartesianFromGeographic)(geographic_point)
    # Longitude is directly translatable to a spherical coordinate
    # θ (azimuth)
    θdeg = GI.x(geographic_point)
    # The polar angle is 90 degrees minus the latitude
    # ϕ (polar angle)
    ϕdeg = 90 - GI.y(geographic_point)
    # Here, we use `sincosd` for two reasons:
    # 1. It's faster than calling `sind` and `cosd` separately
    # 2. `sind` and `cosd` are more numerically stable than 
    # `sin(deg2rad(...))` and `cos(deg2rad(...))` are.
    # If we have coordinates available in radians, using them directly is better.
    sinϕ, cosϕ = sincosd(ϕdeg)
    sinθ, cosθ = sincosd(θdeg)
    # Since this is the unit sphere, the radius is assumed to be 1,
    # and we don't need to multiply by it.
    return GeometryBasics.Point3(
        sinϕ * cosθ,
        sinϕ * sinθ,
        cosϕ
    )
end

struct GeographicFromUnitCartesian <: CoordinateTransformations.Transformation
end

# median time 36 ns
function (::GeographicFromUnitCartesian)(xyz::AbstractVector)
    @assert length(xyz) == 3 "GeographicFromUnitCartesian expects a 3D Cartesian vector"
    sph = SphericalFromCartesian()(xyz)
    return GeometryBasics.Point2(
        rad2deg(sph.θ),
        rad2deg(sph.ϕ)
    )
end

"""
    UnitSphericalCap{T}

A spherical cap, defined by a point on the sphere and a radius.
"""
struct UnitSphericalCap{T}
    point::GeometryBasics.Point3{T}
    radius::T
end

function UnitSphericalCap(point::GeometryBasics.Point3{T1}, radius::T2) where {T1,T2}
    return UnitSphericalCap{promote_type(T1, T2)}(point, radius)
end

UnitSphericalCap(point, radius) = UnitSphericalCap(GI.trait(point), point, radius)

UnitSphericalCap(::GI.PointTrait, point, radius) = UnitSphericalCap(GeometryOpsCore.Spherical(), GI.PointTrait(), point, radius)

function UnitSphericalCap(m::GeometryOpsCore.Spherical{T}, ::GI.PointTrait, point, radius) where {T}
    radius_on_unit_sphere = radius / m.radius
    return UnitSphericalCap(UnitCartesianFromGeographic()(point), radius_on_unit_sphere)
end


function spherical_cap_edge_points(cap::UnitSphericalCap{T}, num_points::Int=100, R::T=one(T)) where {T}
    # The problem with perpendicular-vector-based approaches is that they are vulnerable to oversensitivity, and the normalization actually makes the problem worse.
    # Instead, we construct points on the edge of the cap by rotating the points on the edge of the North Pole,
    # then rotating them to the correct orientation.

    edge_points = Vector{GeometryBasics.Point3{T}}(undef, num_points)
    # NOTE: this is the physical convention of spherical coordinates,
    # where θ is the polar angle (colatitude) and φ is the azimuthal angle.
    center_theta = acos(cap.point[3])
    center_phi = atan(cap.point[2], cap.point[1])
    # see https://math.stackexchange.com/questions/275134/rotation-matrix-sending-the-north-pole-to-a-specified-point-on-the-sphere 
    # for the rotation matrix
    # we directly convert to a static matrix to pre-compute and materialize the matrix,
    rotation = Rotations.RotZX{T}(center_phi + pi / 2, center_theta) |> Rotations.RotMatrix{3}

    sinphi, cosphi = sincos(cap.radius)
    for (i, theta) in enumerate(Base.range(zero(T), stop=2π, length=num_points))
        sintheta, costheta = sincos(theta)
        # Generate a series of Cartesian points on the edge of the North Pole
        cartesian_point = GeometryBasics.Point3{T}(sinphi * costheta, sinphi * sintheta, cosphi)
        # Rotate the points to the correct orientation, such that the "new" north pole is actually the cap's point
        edge_points[i] = LinearAlgebra.normalize(rotation * cartesian_point) * R
    end

    return edge_points
end


function to_latlong_polygon(cap::UnitSphericalCap, n_points::Int)
    points = GeographicFromUnitCartesian().(spherical_cap_edge_points(cap, n_points))
    return GI.Polygon([GI.LinearRing(points)]; crs=GFT.EPSG(4326))
end

import MakieCore
function MakieCore.convert_arguments(::Type{MakieCore.Mesh}, cap::UnitSphericalCap)
    # get the edge points
    points = spherical_cap_edge_points(cap, 50, 1.0)
    # add the center point
    push!(points, cap.point)
    # create the faces
    faces = [GeometryBasics.TriangleFace(i, i + 1, length(points)) for i in 1:(length(points)-2)]
    return (GeometryBasics.normal_mesh(points, faces),)
end

function trace_upstream(COMID, river_traces)
    upstream_trace = [Vector{eltype(COMID)}() for _ in COMID]
    ind = falses(length(river_traces))
    
    for (i, id) in enumerate(COMID)   
        for (j, haystack) in enumerate(river_traces)
            ind[j] = in(id, haystack)
        end
        upstream_trace[i] = copy(COMID[ind])
    end
    return upstream_trace
end


function trace_downstream(start_id::Integer, ids, next_ids; maxiters=length(ids))
    # prog = ProgressUnknown(; dt = 0.1, desc = "Tracing downstream...")
    if start_id == 0
        return Base.nonmissingtype(eltype(ids))[]
    end
    visited_ids = Base.nonmissingtype(eltype(ids))[]
    current_id = start_id
    current_idx = searchsortedfirst(ids, start_id)
    iter_count = 0
    while current_id != 0 && iter_count <= maxiters
        current_idx = searchsortedfirst(ids, current_id)
        push!(visited_ids, current_id)
        current_id = next_ids[current_idx]
        # next!(prog)
        iter_count += 1
    end
    # finish!(prog)
    return visited_ids
end

function trace_downstream_idx(start_id::Integer, ids, next_ids; maxiters=length(ids))
    # prog = ProgressUnknown(; dt = 0.1, desc = "Tracing downstream...")
    if start_id == 0
        return Base.nonmissingtype(eltype(ids))[]
    end
    visited_ids = Base.nonmissingtype(eltype(ids))[]
    current_id = start_id
    current_idx = searchsortedfirst(ids, start_id)
    iter_count = 0
    while current_id != 0 && iter_count <= maxiters
        current_idx = searchsortedfirst(ids, current_id)
        push!(visited_ids, current_idx)
        current_id = next_ids[current_idx]
        # next!(prog)
        iter_count += 1
    end
    # finish!(prog)
    return visited_ids
end



function haversine_distance(x, y)
    return haversine_distance(GOC.Spherical(), x, y)
end

function haversine_distance(manifold::GOC.Spherical, x, y)
    return haversine_distance(manifold, GI.trait(x), x, GI.trait(y), y)
end

haversine_distance(manifold::GOC.Spherical, ::Union{GI.AbstractTrait, Nothing}, x, ::Union{GI.AbstractTrait, Nothing}, y) = error("Not implemented yet for $(GI.trait(x)) and $(GI.trait(y)), please file an issue if you need this.")

haversine_distance(manifold::GOC.Spherical, ::GI.PointTrait, x, ::GI.PointTrait, y) = Distances.Haversine(manifold.radius)((GI.x(x), GI.y(x)), (GI.x(y), GI.y(y)))

haversine_distance(manifold::GOC.Spherical, ::GI.PointTrait, x, ::GI.AbstractCurveTrait, y) = GO.applyreduce(min, GI.PointTrait(), y; init = Inf) do yp
    haversine_distance(manifold, x, yp)
end



"""
    flux_accumulate!(river_inputs, id, nextdown_id, headbasin, majorbasin_id)

Accumulate river fluxes by propagating flows downstream through a river network.

# Arguments
- `river_inputs`: DimensionalArray of basin/river inputs [m s⁻¹] with dimensions [Ti, ID]
- `id`: Vector of river/basin IDs 
- `nextdown_id`: Vector of next downstream river/basin IDs (0 indicates no downstream reach)
- `headbasin`: Boolean vector indicating headwater basins (true = no upstream reaches)
- `majorbasin_id`: Vector of major/main basin IDs used for parallel processing

# Details
For each major basin, this function:
1. Identifies headwater basins and their downstream connections
2. Iteratively propagates flows downstream by adding upstream flows to downstream reaches
3. Continues until all basins in the network have been processed
4. Performs accumulation in-place, modifying the input `river_inputs` array

The function processes each major basin in parallel using multiple threads.
flux_accumulate! is insanely fast.
"""
function flux_accumulate!(river_inputs, id, nextdown_id, headbasin, majorbasin_id; accumulated_length = false)
    # sum river inputs going downstream to get total river flux     

    Threads.@threads for majorbasin in unique(majorbasin_id) #[2 min]
        # majorbasin = unique(majorbasin_id)[1]
        index = majorbasin_id .== majorbasin

        id0 =  id[index]
        nextdown_id0 = nextdown_id[index]
        headbasin0 = headbasin[index]
        
        # identify head basins that do not flow into other basins
        touched = nextdown_id0 .== 0

        # identify head basins that flow into other basins
        flowdown_now = headbasin0 .& (.!touched) # these basins are ready to flow down
        flowdown_later = copy(flowdown_now) # these basins are waiting for other basins catch fill up
        
        # exclude basins that will later be touched by flowdown
        flowdown_now_ids = setdiff(unique(nextdown_id0[flowdown_now]), unique(nextdown_id0[.!(touched .| flowdown_now)]))
        flowdown_now[flowdown_now] = in.(nextdown_id0[flowdown_now], Ref(flowdown_now_ids)) 
        flowdown_later = flowdown_later .& (.!flowdown_now)

        while any(flowdown_now)
            id1 = id0[flowdown_now]
            nextdown_id1 = nextdown_id0[flowdown_now]

            # a loop here seems to be faster than a view
            if accumulated_length
                 for to_id in unique(nextdown_id1)
                    from_id = id1[nextdown_id1.==to_id]
                    river_inputs[At(to_id)] += maximum(river_inputs[At(from_id)])
                end
            else
                for (from_id, to_id) in zip(id1, nextdown_id1) 
                    river_inputs[:, At(to_id)] .+= parent(river_inputs[:, At(from_id)])
                end
            end  

            touched[:] = flowdown_now .| touched

            if all(touched)
                break
            else
                flowdown_now[:] = (in.(id0, Ref(unique(nextdown_id0[flowdown_now]))) .| flowdown_later) .& (.!touched)
                flowdown_later[:] = flowdown_now

                flowdown_now_ids = setdiff(unique(nextdown_id0[flowdown_now]), unique(nextdown_id0[.!(touched .| flowdown_now)]))
                flowdown_now[flowdown_now] = in.(nextdown_id0[flowdown_now], Ref(flowdown_now_ids))
                flowdown_later[:] = flowdown_later .& (.!flowdown_now)
            end
        end

        # this is a sanity check that all basins were touched
        if !all(touched)
            @warn "for major basin $(majorbasin) $(sum(touched)) basins were not touched"
        end
    end
    return river_inputs
end

function linear_reservoir_impulse_response_monthly(Tb)
    #Tb = 45 #days
    Vb = 1 #mm

    Qb =[]
    for i in 1:10000
        if i < 30
            Qsb = 1 #mm/day
        else
            Qsb = 0 #mm/day
        end
        Qb = push!(Qb, Vb/Tb)
        Vb += (Qsb - Qb[i])
    end


    c = 1
    dd = 30
    Qfrac = []
    for _ in 1:10
        push!(Qfrac, round(sum(Qb[c:c+dd-1]) ./ sum(Qb), digits=2))
        c += dd
    end

    Qfrac = Qfrac[1:findfirst(Qfrac .== 0)-1]

    return Qfrac
end


function apply_vector_impulse_resonse!(M, impulse_resonse)
    nir = length(impulse_resonse)
    m, n, k = size(M)
    Q1 = zeros(eltype(M), (k + nir - 1))

    for i in 1:m
        for j in 1:n
            Q1 .= 0
            Q0 = M[i, j, :]
            for r in 1:k
                Q1[r:r+nir-1] .+= Q0[r] .* impulse_resonse
            end
            M[i, j, :] = Q1[1:k]
        end
    end
end



# load in river reaches
function river_reaches(rivers_paths; col_names=nothing)
    rivers = fill(DataFrame(),length(rivers_paths))

    Threads.@threads for i in eachindex(rivers_paths)
        if isnothing(col_names)
            rivers[i] = GeoDataFrames.read(rivers_paths[i])
        else
            rivers[i] = GeoDataFrames.read(rivers_paths[i])[:, col_names]
        end
    end
    
    rivers = reduce(vcat, rivers)
    sort!(rivers, [:COMID])
    return rivers
end


function river_cumulative_lengths(rivers)

    # add basin02 and headbasin flags to glacier rivers
    rivers[!, :basin02] = floor.(Int16, rivers.COMID ./ 1000000)
    rivers[!, :headbasin] = rivers.up1 .== 0

    # if there is no input from upstream then it is a head basin
    # This is needed as only a subset of the full river network is routed
    rivers[!, :headbasin] .|= .!in.(rivers.COMID, Ref(unique(rivers.NextDownID)))

    dcomid = Dim{:COMID}(rivers.COMID)
    river_lengths = DimArray(rivers.lengthkm, dcomid; name = "river lengths [km]")

    river_lengths = Altim.flux_accumulate!(river_lengths, rivers.COMID, rivers.NextDownID, rivers.headbasin, rivers.basin02; accumulated_length=true)

    return river_lengths
end



function compute_population(gpw_ras, rivers, country_polygons, gmax, runoff, buffer_radii; progress=true)
    m2deg = 0.0000089932
    
    countries = unique(country_polygons[!, :country])
    filter!(!=(""), countries)

    dgmax = Dim{:gmax}(gmax; metadata=Dict(:longname => "maximum monthly glacier fraction of total river flux", :units => "fraction"))
    drunoff = Dim{:runoff}(runoff; metadata=Dict(:longname => "minimum runoff flux during gmax", :units => "m³/s"))

    dbuffer = Dim{:buffer}(buffer_radii; metadata=Dict(:longname => "buffer radius", :units => "m"))
    dcountry = Dim{:country}(countries; metadata=Dict(:longname => "country", :units => ""))

    population = zeros(dcountry, dgmax, drunoff,dbuffer; name="population")

    source_crs1 = GFT.EPSG(4326)
    target_crs1 = GFT.EPSG(3857)

    country_polygons[!, :geometry] = _validate_convert.(country_polygons[:, :geometry])

    prog = progress ? ProgressMeter.Progress(length(population); desc="Extracting population for each county...") : nothing
    
    for country in dcountry # having issues with PROJError: pipeline: Pipeline: too deep recursion when using Threads.@threads

        country_polygon = country_polygons[country_polygons.country.==country, :geometry]
        country_polygon = reduce(LibGEOS.union, country_polygon)

        for gmax in dgmax
            for runoff in drunoff

                river_idx = (rivers[:, :gmax_avg] .>= gmax) .& (rivers[:, :country] .== country) .& (rivers[:, :runoff_max_avg] .>= runoff)

                if !any(river_idx)
                    update!(prog, prog.core.counter + length(dbuffer))
                    continue
                end

                geom1 = rivers[river_idx, :geometry]

                # extract points
                geom1 = reduce(vcat, _get_point_subset.(geom1; subsample=10))
    
                # project to web-mercador for equal area projection
                geom1 = GO.reproject(geom1, source_crs1, target_crs1; always_xy=true)

                for buffer_radius in dbuffer
                    # printstyled("\ncountry: $country, gmax: $gmax, buffer_radius: $buffer_radius \n"; color=:red)
                    
                    # buffer
                    geom = GO.buffer(geom1, buffer_radius)

                    # project back to geographic coordinates
                    geom = GO.reproject(geom, target_crs1, source_crs1; always_xy=true)

                    # GO.union is hopeless
                    #GO.union(geom[1], geom[2]; target=GI.PolygonTrait())
                    #GI.Polygon(GI.LinearRing(collect(GI.getpoint(geom[1]))))

                    # convert geometry to LibGEOS to use thier union function
                    geom = _validate_convert.(geom)
                    geom = filter(!isnothing, geom)

                    # union all the polygons
                    try
                        geom = reduce(LibGEOS.union, geom)
                    catch
                        printstyled("\nerror unioning geometry for country: $country, gmax: $gmax, buffer_radius: $buffer_radius\n"; color=:red)
                        progress ? next!(prog) : nothing
                        continue
                    end

                    # trim to country boundaries
                    geom = LibGEOS.intersection(country_polygon, geom)

                    total_pop = Rasters.zonal(sum, gpw_ras; of=geom, skipmissing=true)
                    population[At(country), At(gmax), At(runoff), At(buffer_radius)] += total_pop

                    progress ? next!(prog) : nothing
                end
            end
        end
    end
    finish!(prog)

    return population
end


 function _get_point_subset(geom; subsample=10)

    pts = GI.getpoint(geom)
    if length(pts) < subsample
        pts = GO.centroid(geom)
    else
        pts = collect(pts)
        idx = collect(1:10:length(pts))
        idx[end] = length(pts)
        pts = pts[idx]
    end

    return pts
end 

function _validate_convert(geom)
    try
        geom = GI.convert(LibGEOS, geom)
        if !LibGEOS.isValid(geom)
            geom = nothing
        end
    catch
        geom = nothing
    end
    return geom
end
