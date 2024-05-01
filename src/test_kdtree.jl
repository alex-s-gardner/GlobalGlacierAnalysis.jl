using Altim
using Extents
using GeoTiles
using NearestNeighbors
using ProfileView

project_id = :v01;
domain = :landice;

#domain = :all;
extent = Extent(Lon=(-142.0, -141.0), Lat=(58.0, 61.0))

paths = project_paths(project_id=project_id);
df = GeoTiles.readall(paths.icesat2.geotile; suffix=".arrow", extent=extent);


t = Altim.decimalyear.(reduce(vcat, df[2].datetime))

function test(t)
    (tmin, tmax) = extrema(t)
    halfwidth = 1 / 365
    t0 = floor(tmin + halfwidth):halfwidth:ceil(tmax + halfwidth)

    # build kdtree
    tree = KDTree(t'; leafsize=10)
    idx = inrange(tree, collect(t0)', halfwidth)
    return t0, idx
end

 begin
    a, b = test(t)::Tuple{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64},Vector{Vector{Int64}}}
end
