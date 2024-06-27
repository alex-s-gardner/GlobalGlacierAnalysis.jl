"""
    Extract point elevations, slopes and curvature from global DEMS
"""

using Altim
#using Infiltrator

force_remake = false;
project_id = :v01;
geotile_width = 2;
slope = true;
curvature = true;

extent = nothing
dems = [:rema_v2_10m, :cop30_v2, :arcticdem_v4_10m]#[:cop30_v2]; #, :nasadem_v1, :rema_v2_10m, :arcticdem_v4_10m]

paths = project_paths(; project_id);
geotiles = Altim.geotiles_w_mask(geotile_width);
products = project_products(; project_id);

geotiles = geotiles[geotiles.landice_frac .> 0, :];

# <><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><><><><><><>
force_remake = false;
products = (icesat = products.hugonnet,)
#ind = findfirst(geotiles.id .== "lat[+54+56]lon[+096+098]")
#geotiles = geotiles[ind:ind,:]
#dems = [:arcticdem_v4_10m]

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
Altim.geotile_extract_dem(products, dems, geotiles, paths; slope, curvature, force_remake)
