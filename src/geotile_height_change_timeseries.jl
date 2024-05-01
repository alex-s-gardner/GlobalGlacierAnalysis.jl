using Altim
using GeoTiles

#force_remake = false;
project_id = :v01;
geotile_width = 2;
domain = :landice;

paths = project_paths(; project_id);
geotiles = project_geotiles(; geotile_width, domain, extent);
products = project_products(; project_id);

folder = paths.icesat2.geotile;
id = geotiles.id[1];

suffix = "arrow"
file = GeoTiles.path2tile(folder, id, suffix)
df = GeoTiles.read(file)

suffix = "cop30_v2"
file = GeoTiles.path2tile(folder, id, suffix)
df = GeoTiles.read(file)

suffix = "canopyh"
file = GeoTiles.path2tile(folder, id, suffix)
df = GeoTiles.read(file)
