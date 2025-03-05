# This script extracts canopy height data for geotiles:
#
# Configuration:
# - Uses project ID v01 and 2-degree geotile width
# - Processes only geotiles containing land ice
# - Uses ETH Global Canopy Height 10m 2020 dataset
#
# Processing:
# 1. Sets up project paths and loads product definitions
# 2. Filters geotiles to only include those with land ice
# 3. Reads canopy height data from source GeoTIFF
# 4. Extracts canopy height values for each geotile using geotile_pointextract()
#    - Processes all altimetry missions
#    - Only processes if force_remake=true or output doesn't exist
#    - Uses nodata value of 255


using Altim, GeoArrays, Dates, Arrow, DataFrames

## Download canopy height data in teminal
#= 
;cd /mnt/devon-r0/shared_data/canopy_height/
;aria2c -x10 https://share.phys.ethz.ch/~pf/nlangdata/ETH_GlobalCanopyHeight_10m_2020_version1/ETH_GlobalCanopyHeight_10m_2020_version1.tar.gz 
;tar -xvzf ETH_GlobalCanopyHeight_10m_2020_version1.tar.gz
=#

# Parameters: user defined 
force_remake = true
project_id = :v01;
geotile_width = 2;
domain = :glacier; # :glacier -or- :landice
missions = (:icesat2,); # (:icesat2, :icesat, :gedi, :hugonnet)

# Initialize: paths, products, geotiles
paths = project_paths(; project_id);
products = project_products(; project_id);
geotiles = Altim.geotiles_w_mask(geotile_width);

# Subset: region & mission 
geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :];
isa(missions, Tuple) || error("missions must be a tuple... maybe you forgot a trailing comma for single-element tuples?")
products = getindex(products, missions)

# Execute: extract canopy height
ga = GeoArrays.read(setpaths().canopyheight_10m_v1, masked=false);
nodatavalue = 255;
Altim.geotile_pointextract(geotiles, [paths[mission].geotile for mission in missions], ga; var_name = :canopyh, job_ids = [missions...], nodatavalue = nodatavalue, force_remake = force_remake)