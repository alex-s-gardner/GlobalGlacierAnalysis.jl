
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
include("/home/singhvi/.julia/dev/river_mapping/STRTreeAbstractTreesExt.jl")
# Progress meter
using ProgressMeter
#=
include("trace_down_rivers.jl")

df = DataFrame(Shapefile.Table("data/riv_pfaf_81_MERIT_Hydro_v07_Basins_v01_GLDAS_COR.shp"))

maxidx = max(maximum(df.COMID), minimum(df.NextDownID))

sparse_adjmat = sparse(df.COMID, replace(df.NextDownID, 0 => 1), ones(Float32, nrow(df)), maxidx, maxidx)

traces = trace_downstream_idx(df.COMID[33636], df.COMID, df.NextDownID)

river_ls = reverse_point_order(linestring_from_idx(df, traces))
=#


glaciers = GeoDataFrames.read("/mnt/bylot-r3/data/GlacierOutlines/rgi60/rgi60_Global.gpkg")

# We load the DEM lazily as a VRT.  This requires ArchGDAL.jl.
dem = Raster("/mnt/devon-r2/shared_data/copernicus-dem-30m/DEM.vrt"; lazy = true)

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
function is_extent_suitable(raster, node; max_n_chunks = 100^2, chunk_size = 128^2)
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


function get_minimum_elevations(dem, geometries; progress = true)
    prog = progress ? ProgressMeter.Progress(length(geometries); desc = "Getting minimum elevations") : nothing

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
    for node::Union{SortTileRecursiveTree.STRNode, SortTileRecursiveTree.STRLeafNode} in iterator
        # We perform rasterization if the node's extent is sufficiently small,
        # or if the node is a leaf node.
        loaded_dem = if is_extent_suitable(dem, node)
            @debug("Loading DEM for node with extent $(GI.extent(node))")
            # Here, we know the node is smaller than the maximum size we're willing to load,
            # so we load the section of the DEM covering the node into memory,
            # and use that for the zonal operation.
            Rasters.read(Rasters.crop(dem; to = GI.extent(node)::GI.Extents.Extent)::Raster)::Raster
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
            of = view(geometries, indices), 
            skipmissing = false,  # provide a full masked Raster to the zonal operation
            threaded = true,      # thread over the geometries
            progress = false,     # don't show progress since we are already showing global progress
        )::Vector{Union{Missing, Tuple{Float64, Float64}}}
        # Update the progress meter
        progress && update!(prog, prog.counter + length(indices))
    end

    return minimum_elevations
end


