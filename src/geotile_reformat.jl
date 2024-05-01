using Altim
using Extents
using Arrow
using UUIDs

force_remake = true;
project_id = :v01
geotile_width = 2
domain = :landice

paths = project_paths(project_id = project_id)
geotiles = project_geotiles(; geotile_width = geotile_width, domain = domain)
products = project_products(; project_id=project_id)

dems = [:cop30_v2, :nasadem_v1, :rema_v2_10m, :arcticdem_v3_10m]; 
#= --------------------------------------------------------------------------
@warn "--- only processing a subset of geotiles ----"
region = :WNA
extent, epsg = region_extent(region);
geotiles = geotile_subset!(geotiles, extent)
dems = [:cop30_v2]; 
# --------------------------------------------------------------------------
=#


#missions = [:hugonnet];
#sensor_halfwidth = [30];
#sensor_kernel = [:average]
dems = [:cop30_v2]

extent = Extent(X=(-180, 180), Y=(+37, +90))
geotiles = geotile_subset!(geotiles, extent)
products = (hugonnet = products.icesat2,)
for (i, product) in enumerate(products)
    for geotile in eachrow(geotiles)
        fn = joinpath(paths[product.mission].geotile,  geotile.id*".arrow")
        df = DataFrame(Arrow.Table(fn))
    end
end
