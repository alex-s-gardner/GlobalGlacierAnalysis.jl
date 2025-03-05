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

# Parameters: user defined 
force_remake = true
project_id = :v01;
geotile_width = 2;
domain = :glacier; # :glacier -or- :landice
missions = (:icesat2,); # (:icesat2, :icesat, :gedi, :hugonnet)

vars2extract = [:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean];

# Initialize: paths, products, geotiles
paths = project_paths(; project_id);
products = project_products(; project_id);
geotiles = Altim.geotiles_w_mask(geotile_width);

# Subset: region & mission 
geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :];
isa(missions, Tuple) || error("missions must be a tuple... maybe you forgot a trailing comma for single-element tuples?")
products = getindex(products, missions)

# Execute: extract masks
for product in products
    geotile_extract_mask(geotiles, paths[product.mission].geotile; vars=vars2extract, job_id=product.mission, force_remake=force_remake)
end