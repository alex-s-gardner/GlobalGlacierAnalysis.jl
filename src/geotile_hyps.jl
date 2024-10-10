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