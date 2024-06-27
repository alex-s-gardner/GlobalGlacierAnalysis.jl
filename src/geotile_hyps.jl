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
    mask = :glacier_b10km #:land #
    subset = :glacier_frac

    runid = "geotile_$(mask)_hyps"
    binned_folder = analysis_paths(; geotile_width).binned
    out_file = joinpath(binned_folder, "$(runid).arrow");

    dbin = 100.;
    bin_centers = dbin/2:dbin:10000-dbin/2;
    bin_edges = 0:dbin:10000;
    excludemask_flag = false

    out_var_name = Symbol("$(mask)_area_km2")
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

    geotiles[!, out_var_name] = [zeros(size(bin_centers)) for r in 1:nrow(geotiles)];
    mask_frac = Symbol("$(mask)_frac")

    if mask == :land
        shp = Symbol("$(:glacier)_shp")
        fn_shp = Altim.pathlocal[shp]
        excludefeature = Shapefile.Handle(fn_shp)
    else
        excludefeature = nothing
    end
end

#Threads.@threads for geotile in eachrow(geotiles)
dfr = eachrow(geotiles)
Threads.@threads for geotile in dfr[geotiles[:,subset].> 0]
    geotile_binarea!(geotile, ras, feature, bin_edges; invert, excludefeature)
end
Arrow.write(out_file, geotiles::DataFrame)