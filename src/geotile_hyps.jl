# This script calculates hypsometric (elevation-binned) area distributions for different masks:
#
# Configuration:
# - Uses 2-degree geotile width
# - Processes COP30 v2 DEM for elevation data
# - Calculates areas for glacier, land, RGI7 glaciers, and glacier buffer zones
#
# Processing steps:
# 1. Set up configuration:
#    - Load required packages and utilities
#    - Define geotile width, DEM source, and masks to process
#    - Set up elevation bin ranges
#
# 2. For each mask type (glacier, land, etc):
#    - Set up output paths and files
#    - Load DEM raster and mask shapefiles
#    - Filter geotiles to those with glacier coverage
#    - Calculate binned areas using parallel processing
#    - Save results to Arrow files
#
# Special handling:
# - Land mask is inverted from water mask
# - Land calculations exclude glacier areas
# - Uses threading for parallel processing
# - Only processes geotiles with glacier coverage

begin
    using Altim
    using Rasters
    using Shapefile
    using DataFrames 
    using Arrow
    using BinStatistics
    include("utilities_hyps.jl")

    geotile_width = 2;
    raster_file = :cop30_v2; #, :nasadem_v1, :rema_v2_10m, :arcticdem_v3_10m]
    subset = :glacier_frac
    masks = [:glacier, :land, :glacier_rgi7, :glacier_b1km, :glacier_b10km]
end

for mask in masks
    
    runid = "geotile_$(mask)_hyps"
    binned_folder = analysis_paths(; geotile_width).binned
    out_file = joinpath(binned_folder, "$(runid).arrow");

    height_range, height_center = Altim.project_height_bins()
    Δh = abs(height_center[2] - height_center[1])

    excludemask_flag = false

    var_name = Symbol("$(mask)_area_km2")
    geotiles = Altim.geotiles_w_mask(geotile_width);
    fn_raster = Altim.pathlocal[raster_file];
    ras = Raster(fn_raster; lazy=true);

    if mask == :land
        shp = Symbol("$(:water)_shp");
        invert = true
    else
        shp = Symbol("$(mask)_shp");
        invert = false
    end

    fn_shp = Altim.pathlocal[shp];
    feature = Shapefile.Handle(fn_shp);

    geotiles[!, var_name] = [zeros(size(height_centers)) for r in 1:nrow(geotiles)];
    mask_frac = Symbol("$(mask)_frac")

    if mask == :land
        shp = Symbol("$(:glacier)_shp")
        fn_shp = Altim.pathlocal[shp]
        excludefeature = Shapefile.Handle(fn_shp)
    else
        excludefeature = nothing
    end

    dfr = eachrow(geotiles)
    Threads.@threads for geotile in dfr[geotiles[:,subset].> 0]
        geotile_binarea!(geotile, ras, feature, height_range; invert, excludefeature, var_name)
    end
    Arrow.write(out_file, geotiles::DataFrame)
end