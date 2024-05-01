using Altim
using Extents
using Arrow
using DataFrames
using ProfileView
using Plots

geotile_width = 2; #geotile width [degrees]
grid = (node_spacing = 100, node_width = 500); # [(node_spacing = 500, node_width = 500 * 2)]

project_id = :v01;
domain = :landice;
#domain = :all;

paths = project_paths(project_id = project_id);
products = project_products(; project_id=:v01);
geotiles = project_geotiles(; geotile_width = geotile_width, domain = domain);

dems = [:cop30_v2]; 
# missions = [:gedi, :icesat, :icesat2, :hugonnet]
# --------------------------------------------------------------------------
ext = Extent(X=(-126.9, -126.1), Y=(51.1, 51.8));
geotiles = geotile_subset!(geotiles, ext);

geotile = first(geotiles);
dem = first(dems);


for product in products
        Altim.geotile_coregister(geotile, paths[product.mission].geotile, dem)
end

#histogram(df.dx, bins=-20:1:20, label="dx");
#histogram(df.dy, bins=-20:1:20, label="dy");
#histogram(df.dz, bins = -20:1:20, label="dz");
