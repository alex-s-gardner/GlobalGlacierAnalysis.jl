# This script extracts surface type masks for geotiles:
#
# Configuration:
# - Uses project ID v01 and 2-degree geotile width
# - Filters to only include geotiles with land ice
# - Extracts masks for: floating ice, glacier ice, inland water, 
#   land, land ice, and ocean
#
# Processing:
# 1. Sets up project paths and loads product definitions
# 2. Filters geotiles to only include those with land ice
# 3. For each product:
#    - Extracts mask data using geotile_extract_mask()
#    - Only processes if force_remake=true or output doesn't exist

# add packages
using Altim

# set paths
force_remake = false

project_id = :v01;
geotile_width = 2;
domain = :landice;

paths = project_paths(; project_id);
products = project_products(; project_id);

geotiles = Altim.geotiles_w_mask(geotile_width);
if domain == :landice
    geotiles = geotiles[geotiles.landice_frac .> 0, :];
end

# --------------------------------------------------------------------------
#begin
# @warn "--- only processing a subset of geotiles ----"
#    region = :RGI02;
#    ext, epsg = region_extent(region);
#    geotiles = geotile_subset!(geotiles, ext);
#    products = (hugonnet=products.hugonnet, )
#end
# --------------------------------------------------------------------------

# geotiles = geotiles[ind,:]
vars = [:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean];

for product in products
    geotile_extract_mask(geotiles, paths[product.mission].geotile; vars=vars, job_id=product.mission, force_remake=force_remake)
end