"""
    coregester_hugonnet.jl

Coregister elevation data products to a reference DEM for a specific geographic region.

This script:
1. Sets up configuration parameters for geotile processing (width, grid spacing)
2. Defines project parameters and domain (landice)
3. Loads project paths, products, and geotiles
4. Subsets geotiles to a specific region in British Columbia, Canada
5. Coregisters each product to the Copernicus 30m DEM (cop30_v2)

The coregistration process aligns elevation data from different sources to minimize
systematic offsets, improving the accuracy of elevation change measurements.
"""

using GlobalGlacierAnalysis
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
        GlobalGlacierAnalysis.geotile_coregister(geotile, paths[product.mission].geotile, dem)
end

