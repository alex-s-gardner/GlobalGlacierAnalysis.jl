# This script extracts elevation data and derivatives from global DEMs:
#
# Configuration:
# - Uses project ID v01 and 2-degree geotile width
# - Extracts slope and curvature in addition to elevation
# - Processes 4 DEMs: REMA v2, COP30 v2, ArcticDEM v4, NASADEM v1
#
# Processing:
# 1. Sets up project paths and loads product definitions
# 2. Filters geotiles to only include those with land ice
# 3. Extracts DEM data for each geotile using Altim.geotile_extract_dem()
#    - Includes elevation, slope and curvature
#    - Only processes if force_remake=true or output doesn't exist

using Altim

# Parameters: user defined 
force_remake = false
project_id = :v01;
geotile_width = 2;
domain = :landice; # :glacier -or- :landice
missions = (:icesat2, :icesat, :gedi, :hugonnet,); # (:icesat2, :icesat, :gedi, :hugonnet)
hugonnet_unfiltered = true;

slope = true;
curvature = true;
dems2extract = [:rema_v2_10m, :cop30_v2, :arcticdem_v4_10m, :nasadem_v1]

# Initialize: paths, products, geotiles
paths = project_paths(; project_id);
products = project_products(; project_id);
geotiles = Altim.geotiles_w_mask(geotile_width);

# Subset: region & mission 
geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :];
isa(missions, Tuple) || error("missions must be a tuple... maybe you forgot a trailing comma for single-element tuples?")
products = getindex(products, missions)

# Execute: extract dems
Altim.geotile_extract_dem(products, dems2extract, geotiles, paths; slope, curvature, force_remake)

# include hugonnet unfiltered
if hugonnet_unfiltered && (:hugonnet in missions)
    paths = Altim.update_geotile_path(paths; mission = :hugonnet, path_replace ="/2deg" => "/2deg_unfiltered")
    Altim.geotile_extract_dem(products[(:hugonnet,)], dems2extract, geotiles, paths[(:hugonnet,)]; slope, curvature, force_remake)
end
