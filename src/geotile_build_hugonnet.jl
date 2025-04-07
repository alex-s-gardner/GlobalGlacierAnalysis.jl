# This script processes Hugonnet glacier elevation change data into geotiles
#
# Processing steps:
# 1. Set up configuration:
#    - Define geotile width (2 degrees)
#    - Load project paths and create geotile grid
#    - Optional: Filter geotiles to specific regions/extents
#
# 2. Create output directory for Hugonnet geotiles if needed
#
# 3. Load Hugonnet data:
#    - Get catalogue of Hugonnet elevation change stacks
#    - Detect if using old/new data format
#
# 4. Build geotiles:
#    - Iterate through each geotile
#    - Process Hugonnet data into geotile format
#    - Save geotiles to output directory

using Altim
using Extents


# Parameters: user defined 
force_remake = false
project_id = :v01;
geotile_width = 2;
domain = :landice; # :glacier -or- :landice

# Initialize: paths, products, geotiles
paths = project_paths(; project_id);
products = project_products(; project_id);
geotiles = Altim.geotiles_w_mask(geotile_width);

# Subset: region & mission 
geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :];
isa(missions, Tuple) || error("missions must be a tuple... maybe you forgot a trailing comma for single-element tuples?")
products = getindex(products, missions)


# Execute: build geotiles from Hugonnet data

# make directly if it doesn't exist
if !isdir(paths.hugonnet.geotile)
    mkpath(paths.hugonnet.geotile)
end

# get hstack catalogue
hstacks = hstack_catalogue(paths.hugonnet.raw_data; update_catalogue=force_remake)

if contains(hstacks.path[1], "prefilt")
    old_format = false
else
    old_format = true
end

# build geotiles
begin
    printstyled("building Hugonnet geotiles\n"; color=:blue, bold=true)
    for geotile in eachrow(geotiles)
        geotile_build_hugonnet(geotile, paths.hugonnet.geotile, hstacks; force_remake, old_format)
    end
end