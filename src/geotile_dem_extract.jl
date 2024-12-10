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
#using Infiltrator

force_remake = false;
project_id = :v01;
geotile_width = 2;
slope = true;
curvature = true;

extent = nothing
dems = [:rema_v2_10m, :cop30_v2, :arcticdem_v4_10m, :nasadem_v1]

paths = project_paths(; project_id);
geotiles = Altim.geotiles_w_mask(geotile_width);
products = project_products(; project_id);

geotiles = geotiles[geotiles.landice_frac .> 0, :];

# <><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><><><><><><>
# products = (icesat = products.hugonnet,)
#ind = findfirst(geotiles.id .== "lat[+54+56]lon[+096+098]")
#geotiles = geotiles[ind:ind,:]
#dems = [:arcticdem_v4_10m]

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
Altim.geotile_extract_dem(products, dems, geotiles, paths; slope, curvature, force_remake)
