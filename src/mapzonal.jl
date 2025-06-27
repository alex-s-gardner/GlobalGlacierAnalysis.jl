#using Rasters
import Rasters: _zonal, _run, _maybe_skipmissing_call, _get_geometries
import Rasters: RasterStack, crop, mask
import Rasters: Extents, GI
function mapzonal(reducer, operator, data::RasterStack; of=nothing, kw...)
    _mapzonal(reducer, operator, data, of; kw...)
end
function _mapzonal(reducer, operator, x, ext::Extents.Extent; skipmissing=true)
    cropped = Rasters.crop(x; to=ext, touches=true)
    prod(size(cropped)) > 0 || return missing
    # We can't use skipmissing here, since it doesn't work on rasterstacks
    if skipmissing
        return reducer( # reduce the result of the following operations - many reducers don't support iterators.
            map( # apply operator to each value
                operator,
                Iterators.filter( # skip missing values
                    Base.Fix1(any, !ismissing),
                    Iterators.map( # get the values as named tuples from the rasterstack
                        Base.Fix1(getindex, cropped),
                        eachindex(cropped)
                    )
                )
            )
        )
    else
        return reducer( # apply the reducer function
            Iterators.map(eachindex(cropped)) do I
                operator(cropped[I]) # get the value as namedtuples from the rasterstack and apply the operator
            end
        )
    end
end
function _mapzonal(reducer, operator, x, of; kw...)
    # Otherwise of is a geom, table or vector
    _mapzonal(reducer, operator, x, GI.trait(of), of; kw...)
end
function _mapzonal(reducer, operator, x, ::GI.AbstractFeatureCollectionTrait, fc; kw...)
    _mapzonal(reducer, operator, x, nothing, fc; kw...) # treat this as a table of geometries
end
# This handles tables, featurecollections and vectors of geometries.
function _mapzonal(reducer, operator, x, ::Nothing, data; progress=false, threaded=true, geometrycolumn=nothing, kw...)
    geoms = _get_geometries(data, geometrycolumn)
    n = length(geoms)
    n == 0 && return []
    zs, start_index = _alloc_mapzonal(reducer, operator, x, geoms, n; kw...)
    _run(start_index:n, threaded, progress, "Applying $reducer and $operator to each geometry...") do i
        zs[i] = _mapzonal(reducer, operator, x, geoms[i]; kw...)
    end
    return zs
end
function _alloc_mapzonal(reducer, operator, x, geoms, n; kw...)
    # Find first non-missing entry and count number of missing entries
    n_missing::Int = 0
    z1 = _mapzonal(reducer, operator, x, first(geoms); kw...)
    for geom in geoms
        z1 = _mapzonal(reducer, operator, x, geom; kw...)
        if !ismissing(z1)
            break
        end
        n_missing += 1
    end
    zs = Vector{Union{Missing,typeof(z1)}}(undef, n)
    zs[1:n_missing] .= missing
    # Exit early when all elements are missing
    if n_missing == n
        return zs, n_missing + 1
    end
    zs[n_missing+1] = z1
    return zs, n_missing + 1
end
# This handles single features (just decomposes to geometry)
function _mapzonal(reducer, operator, data, ::GI.AbstractFeatureTrait, feature; kw...)
    _mapzonal(reducer, operator, data, GI.geometry(feature); kw...)
end

# Now, we get into the meat of handling actual geometry
function _mapzonal(reducer, operator, st::AbstractRasterStack, ::GI.AbstractGeometryTrait, geom;
    skipmissing=true, kw...
)
    cropped = Rasters.crop(st; to=geom, touches=true)
    prod(size(cropped)) > 0 || return missing # mapzonal should always return ONE value...
    # masked = mask(cropped; with=geom, kw...)
    mask = boolmask(geom; to=cropped, kw...)
    if skipmissing # TODO: decide on whether to use map or Iterators.map.  Iterators version is faster and allocates less, but is less generally applicable.
        # Maybe users can do `sum \circ collect` if they want to get a vector??
        return reducer(map(operator, Iterators.filter(Base.Fix1(any, !ismissing), view(cropped, mask))))
    else
        return reducer(map(view(cropped, mask)) do val
            operator(val)
        end)
    end
end