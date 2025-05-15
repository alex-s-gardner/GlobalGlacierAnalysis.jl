"""
    geotile_hyps.jl

Generate hypsometric area distributions for geotiles.

This script:
1. Calculates area distributions by elevation for different land cover masks
2. Processes data for glacier, land, and buffer zone areas
3. Uses COP30 DEM (or other specified DEM) for elevation binning
4. Saves results as Arrow files in the analysis/binned directory
5. Supports multithreaded processing for efficiency

The output files contain area in km² binned by elevation for each geotile,
which can be used for hypsometric analysis of glacier and land areas.
"""

begin
    using GlobalGlacierAnalysis
    using Rasters
    using Shapefile
    using DataFrames 
    using Arrow
    using BinStatistics
end

# Parameters: user defined 
begin
    force_remake = false
    geotile_width = 2;
    raster_file = :cop30_v2; #, :nasadem_v1, :rema_v2_10m, :arcticdem_v3_10m]
    domain = :landice # :glacier -or- :landice
    masks = [:glacier, :land, :glacier_rgi7, :glacier_b1km, :glacier_b10km] 
end

for mask in masks
    
    runid = "geotile_$(mask)_hyps"
    binned_folder = analysis_paths(; geotile_width).binned
    out_file = joinpath(binned_folder, "$(runid).arrow");

    if !isfile(out_file) || force_remake
        height_range, height_center = GlobalGlacierAnalysis.project_height_bins()
        Δh = abs(height_center[2] - height_center[1])

        excludemask_flag = false

        var_name = Symbol("$(mask)_area_km2")
        geotiles = GlobalGlacierAnalysis.geotiles_w_mask(geotile_width);
        fn_raster = GlobalGlacierAnalysis.pathlocal[raster_file];
        ras = Raster(fn_raster; lazy=true);

        if mask == :land
            shp = Symbol("$(:water)_shp");
            invert = true
        else
            shp = Symbol("$(mask)_shp");
            invert = false
        end

        fn_shp = GlobalGlacierAnalysis.pathlocal[shp];
        feature = Shapefile.Handle(fn_shp);

        geotiles[!, var_name] = [zeros(size(height_centers)) for r in 1:nrow(geotiles)];
        mask_frac = Symbol("$(mask)_frac")

        if mask == :land
            shp = Symbol("$(:glacier)_shp")
            fn_shp = GlobalGlacierAnalysis.pathlocal[shp]
            excludefeature = Shapefile.Handle(fn_shp)
        else
            excludefeature = nothing
        end

        dfr = eachrow(geotiles)
        Threads.@threads for geotile in dfr[geotiles[!, "$(domain)_frac"] .> 0]
            geotile_binarea!(geotile, ras, feature, height_range; invert, excludefeature, var_name)
        end
        Arrow.write(out_file, geotiles::DataFrame)
    end
end