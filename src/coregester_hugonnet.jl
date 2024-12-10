"""
Coregister altimetry data to reference DEM for a specific geographic region.

This script performs coregistration of various altimetry products to a reference DEM
within a specified geographic region. The coregistration process aligns the altimetry 
data with the DEM to minimize systematic offsets.

Key operations:
1. Sets up project parameters (geotile width, grid spacing, domain)
2. Loads project paths, products and geotiles
3. Defines geographic extent for processing
4. Performs coregistration for each altimetry product

Parameters:
- geotile_width: Width of geotiles in degrees (default: 2)
- grid: Tuple defining node spacing and width for coregistration grid
- project_id: Project identifier (default: :v01) 
- domain: Processing domain (:landice or :all)
- ext: Geographic extent for processing

The script uses the Altim.geotile_coregister function to perform the actual
coregistration for each product.

Dependencies:
Altim, Extents, Arrow, DataFrames, ProfileView, Plots
"""

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
