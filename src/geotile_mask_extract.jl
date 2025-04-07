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
force_remake = false
project_id = :v01;
geotile_width = 2;
domain = :landice; # :glacier -or- :landice
missions = (:icesat2, :icesat, :gedi, :hugonnet,); # (:icesat2, :icesat, :gedi, :hugonnet)
hugonnet_unfiltered = true;
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

# include hugonnet unfiltered
if hugonnet_unfiltered && (:hugonnet in missions)
    paths = Altim.update_geotile_path(paths; mission=:hugonnet, path_replace="/2deg" => "/2deg_unfiltered")
    geotile_extract_mask(geotiles, paths[:hugonnet].geotile; vars=vars2extract, job_id=:hugonnet, force_remake=force_remake)
end