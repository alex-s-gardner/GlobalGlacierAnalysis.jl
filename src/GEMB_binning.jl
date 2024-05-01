# load packages
begin
    using Altim
    using Arrow
    using DataFrames
    using Extents
    using Rasters
    using GeoInterface
    using Shapefile
    using DimensionalData
    using BinStatistics
    using Statistics
    using Dates
    using Plots
    using JLD2
    using FileIO
    using MAT

    # run parameters
    force_remake = false;
    project_id = :v01;
    geotile_width = 2;
    binning_method = "mean"; # "meanmadnorm3"

    # bin method
    runid = "glacier_fac_$(binning_method)_$(project_id)"

    paths = project_paths(; project_id);
    products = project_products(; project_id);
    binned_folder = analysis_paths(; geotile_width).binned

    # load geotile definitions with corresponding hypsometry
    mask = :glacier
    gt_file = joinpath(binned_folder, "geotile_$(mask)_hyps.arrow");
    geotiles = DataFrame(Arrow.Table(gt_file));
    geotiles.extent = Extent.(getindex.(geotiles.extent, 1));

    # filter geotiles
    geotiles = geotiles[geotiles.glacier_frac .> 0.0,:];

    # funtion used for binning data
    if binning_method == "meanmadnorm3"
        binfun::Function = binfun(x) = mean(x[Altim.madnorm(x).<3])
    elseif binning_method == "median"
        binfun::Function = binfun(x) = median(x)
    elseif binning_method == "mean"
        binfun::Function = binfun(x) = mean(x)
    else
        error("unrecognized binning method")
    end

    # open shapefiles
    #glacier_shp = Shapefile.Handle(Altim.pathlocal.glacier_shp);
    glacier_b10km_shp = Shapefile.Handle(Altim.pathlocal.glacier_b10km_shp);
    #glacier_b1km_shp = Shapefile.Handle(Altim.pathlocal.glacier_b1km_shp);

    # define date and hight binning ranges 
    dd =  30;
    date_range = DateTime(1950):Day(dd):now();
    date_center = date_range[1:end-1] .+ Day(dd/2);

    # DO NOT MODIFY: these need to match elevations in hypsometry
    Δh = 100;
    height_range = 0:100:10000;
    height_center = height_range[1:end-1] .+ Δh/2;

    out_file = joinpath(binned_folder, "$(runid).jld2");

    # 3.2 hours for all glaciers, all missions/datasets on 96 threads
    # 45 min for icesat and icesat-2 missions only
    # 2.3 hours for hugonnet only
    showplots = false;

    # initialize dimensional arrays
    ngeotile = nrow(geotiles)
    ndate = length(date_center)
    var_hyps = DimArray(fill(NaN, ngeotile, ndate, length(height_center)), (geotile = geotiles.id, date=date_center, height = height_center));
    nobs_hyps = DimArray(fill(NaN, ngeotile, ndate, length(height_center)), (geotile = geotiles.id, date=date_center, height = height_center));
end;

@time begin
    fac_file = [
        "GEMB_ANT_rv1_0_19500101_20231231.mat",
        "GEMB_ARC_rv1_0_19500101_20231231.mat",
        "GEMB_HMA_rv1_0_19500101_20231231.mat",
        "GEMB_PAT_rv1_0_19500101_20231231.mat"
    ];
    fac_folder = "/home/schlegel/Share/GEMBv1/";

    # this takes 5 m
    fac = Altim.gemb_fac(joinpath.(Ref(fac_folder), fac_file));
end

# 3.66 min for all glacierized geotiles
@time Threads.@threads for geotile in eachrow(geotiles)
    # <><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><><>
    #for geotile in eachrow(geotiles)
    #geotile = geotiles[findfirst((geotiles.rgi3 .> 0.) .& (geotiles.glacier_frac .> 0.3)),:]
    # geotile = geotiles[findfirst(geotiles.id .== "lat[-80-78]lon[+166+168]"), :]
    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    t1 = time();

    # this take 1s
    valid = Altim.within.(Ref(geotile.extent), fac.longitude, fac.latitude)

    if !any(valid)
        continue
    end

    fac0 = fac[valid,:]

    # update glacier mask with high-resolution vector files
    grid_resolution = 0.00027 # ~30m

    x_mask = X(geotile.extent.X[1]:grid_resolution:geotile.extent.X[2], 
        sampling=DimensionalData.Intervals(DimensionalData.Start()));
    y_mask = Y(geotile.extent.Y[1]:grid_resolution:geotile.extent.Y[2], 
        sampling=DimensionalData.Intervals(DimensionalData.Start()));

    mask0 = Raster(zeros(UInt8, y_mask, x_mask))
    setcrs(mask0, EPSG(4326))

    # NOTE: count method is fastest
    mask0 = Rasters.rasterize!(count, mask0, glacier_b10km_shp; threaded=false, 
        shape=:polygon, progress=false, verbose=false, boundary=:center)

    fast_index = true
    var_ind = falses(nrow(fac0))
    if fast_index # fast index is 15x faster than Rasters
        c = round.(Int64, (fac0.longitude .- first(x_mask)) ./ step(x_mask)) .+ 1;
        r = round.(Int64, (fac0.latitude .- first(y_mask)) ./ step(y_mask)) .+ 1;
        pts3 = [CartesianIndex(a, b) for (a, b) in zip(r, c)]
        var_ind = mask0[pts3] .> 0
    else
        #NOTE: 67% of time is taken here for large number of points.
        pts1 = GeoInterface.PointTuple.([((Y=y, X=x)) for (x, y) in 
            zip(fac0.longitude, fac0.latitude)]);
        pts2 = extract(mask0, pts1, atol=grid_resolution/2, index = true, geometry=false);
        var_ind = getindex.(pts2, 2) .> 0
    end;

    if !any(var_ind)
        continue
    end

    # bin data by date and elevation
    minmax_date = extrema(fac0.datetime[var_ind])
    minmax_height = extrema(fac0.height_ref[var_ind])
    date_ind = (date_range .>= minmax_date[1]-Day(dd)) .& (date_range .<= (minmax_date[2]+Day(dd)))
    date_ind_center = findall(date_ind)[1:end-1];
    height_ind = (height_range .>= minmax_height[1]-Δh) .& (height_range .<= (minmax_height[2]+Δh))
    height_ind_center = findall(height_ind)[1:end-1];

    if sum(date_ind) <= 1 || sum(height_ind) <= 1
        continue
    end

    df = binstats(fac0[var_ind, :], [:datetime, :height_ref], 
        [date_range[date_ind], height_range[height_ind]], 
        :fac; col_function=[binfun], missing_bins=true)

    gdf = DataFrames.groupby(df, :datetime)

    # create an array of dh as a function of time and elevation
    obs = fill(NaN, sum(date_ind)-1, sum(height_ind)-1)
    w = fill(Int64(0), sum(date_ind)-1, sum(height_ind)-1)    
    p = sortperm(BinStatistics.bincenter.(gdf[1].height_ref))

    h_center = bincenter.(gdf[1].height_ref)[p];
    t_center = date_center[date_ind_center]

    for (i,df) in enumerate(gdf)
        isval = .!ismissing.(df[p, "fac_binfun"])
        foo1 = @view obs[i,:]
        foo2 = @view w[i, :]
        if any(isval)
            foo1[isval] = df.fac_binfun[p][isval]
            foo2[isval] = df.nrow[p][isval]
        end
    end

    # for troubleshooting 
    showplots && heatmap(h_center, t_center, obs, clim = (-10, 10))

    var_hyps[At(geotile.id), date_ind_center, height_ind_center] = obs
    nobs_hyps[At(geotile.id), date_ind_center, height_ind_center] = w

    t2 = time()
    println("gemb binned: $(geotile.id), time = $(round(t2-t1))s")
end

save(out_file, Dict("fac_hyps" => var_hyps, "nobs_hyps" => nobs_hyps));