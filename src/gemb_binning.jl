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
    runid = "glacier_gemb_$(binning_method)_$(project_id)"

    paths = project_paths(; project_id);
    products = project_products(; project_id);
    binned_folder = analysis_paths(; geotile_width).binned

    # load geotile definitions with corresponding hypsometry
    mask_project  = :glacier
    gt_file = joinpath(binned_folder, "geotile_$(mask_project )_hyps.arrow");
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
end;

file_suffix = "rv1_0_19500101_20231231"
out_df = "/mnt/bylot-r3/data/gemb/raw/$file_suffix.arrow"

# 9 min
if .!(isfile(out_df)) || force_remake   
    gemb_file = [
        "GEMB_ANT_$file_suffix.mat",
        "GEMB_ARC_$file_suffix.mat",
        "GEMB_HMA_$file_suffix.mat",
        "GEMB_PAT_$file_suffix.mat"
    ];

    gemb_folder = "/home/schlegel/Share/GEMBv1/";

    # this takes 4 m
    gemb = Altim.gemb_read(joinpath.(Ref(gemb_folder), gemb_file))


    Arrow.write(out_df, gemb::DataFrame)
else
    gemb = DataFrame(Arrow.Table(out_df))
end

# initialize dimensional arrays
ngeotile = nrow(geotiles)
ndate = length(date_center)

fac_hyps = DimArray(fill(NaN, ngeotile, ndate, length(height_center)), (geotile = geotiles.id, date=date_center, height = height_center));
smb_hyps = DimArray(fill(NaN, ngeotile, ndate, length(height_center)), (geotile=geotiles.id, date=date_center, height=height_center))
nobs_hyps = DimArray(fill(NaN, ngeotile, ndate, length(height_center)), (geotile=geotiles.id, date=date_center, height=height_center));
ec_hyps = DimArray(fill(NaN, ngeotile, ndate, length(height_center)), (geotile=geotiles.id, date=date_center, height=height_center))
acc_hyps = DimArray(fill(NaN, ngeotile, ndate, length(height_center)), (geotile=geotiles.id, date=date_center, height=height_center))
runoff_hyps = DimArray(fill(NaN, ngeotile, ndate, length(height_center)), (geotile=geotiles.id, date=date_center, height=height_center))
melt_hyps = DimArray(fill(NaN, ngeotile, ndate, length(height_center)), (geotile=geotiles.id, date=date_center, height=height_center))
fac_to_depth_hyps = DimArray(fill(NaN, ngeotile, ndate, length(height_center)), (geotile=geotiles.id, date=date_center, height=height_center))
refreeze_hyps = DimArray(fill(NaN, ngeotile, ndate, length(height_center)), (geotile=geotiles.id, date=date_center, height=height_center))
t_air_hyps = DimArray(fill(NaN, ngeotile, ndate, length(height_center)), (geotile=geotiles.id, date=date_center, height=height_center))
rain_hyps = DimArray(fill(NaN, ngeotile, ndate, length(height_center)), (geotile=geotiles.id, date=date_center, height=height_center));

# 2 min for all glacierized geotiles
@time Threads.@threads for geotile in eachrow(geotiles)
    #for geotile in eachrow(geotiles)
    # <><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><><>
    #for geotile in eachrow(geotiles)
    #geotile = geotiles[findfirst((geotiles.rgi3 .> 0.) .& (geotiles.glacier_frac .> 0.3)),:]
    # geotile = geotiles[findfirst(geotiles.id .== "lat[-80-78]lon[+166+168]"), :]
    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    t1 = time();

    # this take 1s
    valid = Altim.within.(Ref(geotile.extent), gemb.longitude, gemb.latitude)

    if !any(valid)
        continue
    end

    gemb0 = gemb[valid,:]

    # update glacier mask with high-resolution vector files
    grid_resolution = 0.001 # ~100m

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
    var_ind = falses(nrow(gemb0))
    if fast_index # fast index is 15x faster than Rasters
        c = round.(Int64, (gemb0.longitude .- first(x_mask)) ./ step(x_mask)) .+ 1;
        r = round.(Int64, (gemb0.latitude .- first(y_mask)) ./ step(y_mask)) .+ 1;
        pts3 = [CartesianIndex(a, b) for (a, b) in zip(r, c)]
        var_ind = mask0[pts3] .> 0
    else
        #NOTE: 67% of time is taken here for large number of points.
        pts1 = GeoInterface.PointTuple.([((Y=y, X=x)) for (x, y) in 
            zip(gemb0.longitude, gemb0.latitude)]);
        pts2 = extract(mask0, pts1, atol=grid_resolution/2, index = true, geometry=false);
        var_ind = getindex.(pts2, 2) .> 0
    end;

    if !any(var_ind)
        continue
    end

    # bin data by date and elevation
    minmax_date = extrema(gemb0.datetime[var_ind])
    minmax_height = extrema(gemb0.height_ref[var_ind])
    date_ind = (date_range .>= minmax_date[1]-Day(dd)) .& (date_range .<= (minmax_date[2]+Day(dd)))
    date_ind_center = findall(date_ind)[1:end-1];
    height_ind = (height_range .>= minmax_height[1]-Δh) .& (height_range .<= (minmax_height[2]+Δh))
    height_ind_center = findall(height_ind)[1:end-1];

    if sum(date_ind) <= 1 || sum(height_ind) <= 1
        continue
    end

    df = binstats(gemb0[var_ind, :], [:datetime, :height_ref], 
        [date_range[date_ind], height_range[height_ind]], 
        [:fac, :smb, :ec, :acc, :runoff, :melt, :fac_to_depth, :refreeze, :t_air, :rain]; col_function=[binfun], missing_bins=true)

    gdf = DataFrames.groupby(df, :datetime)

    # create an array of dh as a function of time and elevation
    fac = fill(NaN, sum(date_ind)-1, sum(height_ind)-1)
    smb = fill(NaN, sum(date_ind) - 1, sum(height_ind) - 1)
    ec = fill(NaN, sum(date_ind) - 1, sum(height_ind) - 1)
    acc = fill(NaN, sum(date_ind) - 1, sum(height_ind) - 1)
    runoff = fill(NaN, sum(date_ind) - 1, sum(height_ind) - 1)
    melt = fill(NaN, sum(date_ind) - 1, sum(height_ind) - 1)
    fac_to_depth = fill(NaN, sum(date_ind) - 1, sum(height_ind) - 1)
    refreeze = fill(NaN, sum(date_ind) - 1, sum(height_ind) - 1)
    t_air = fill(NaN, sum(date_ind) - 1, sum(height_ind) - 1)
    rain = fill(NaN, sum(date_ind) - 1, sum(height_ind) - 1)

    nobs = fill(Int64(0), sum(date_ind)-1, sum(height_ind)-1)    
    p = sortperm(BinStatistics.bincenter.(gdf[1].height_ref))

    h_center = bincenter.(gdf[1].height_ref)[p];
    t_center = date_center[date_ind_center]

    for (i,df) in enumerate(gdf)
        isval = .!ismissing.(df[p, "fac_binfun"])
        foo1 = @view fac[i,:]
        foo2 = @view smb[i, :]
        foo3 = @view nobs[i, :]
        foo4 = @view ec[i, :]
        foo5 = @view acc[i, :]
        foo6 = @view runoff[i, :]
        foo7 = @view melt[i, :]
        foo8 = @view refreeze[i, :]
        foo9 = @view t_air[i, :]
        foo10 = @view rain[i, :]
        foo11 = @view fac_to_depth[i, :]

        if any(isval)
            foo1[isval] = df.fac_binfun[p][isval]
            foo2[isval] = df.smb_binfun[p][isval]
            foo3[isval] = df.nrow[p][isval]
            foo4[isval] = df.ec_binfun[p][isval]
            foo5[isval] = df.acc_binfun[p][isval]
            foo6[isval] = df.runoff_binfun[p][isval]
            foo7[isval] = df.melt_binfun[p][isval]
            foo8[isval] = df.refreeze_binfun[p][isval]
            foo9[isval] = df.t_air_binfun[p][isval]
            foo10[isval] = df.rain_binfun[p][isval]
            foo11[isval] = df.fac_to_depth_binfun[p][isval]
        end
    end

    # for troubleshooting 
    showplots && heatmap(h_center, t_center, fac, clim = (-10, 10))
    fac_hyps[At(geotile.id), date_ind_center, height_ind_center] = fac
    smb_hyps[At(geotile.id), date_ind_center, height_ind_center] = smb
    nobs_hyps[At(geotile.id), date_ind_center, height_ind_center] = nobs
    ec_hyps[At(geotile.id), date_ind_center, height_ind_center] = ec
    acc_hyps[At(geotile.id), date_ind_center, height_ind_center] = acc
    runoff_hyps[At(geotile.id), date_ind_center, height_ind_center] = runoff
    melt_hyps[At(geotile.id), date_ind_center, height_ind_center] = melt
    refreeze_hyps[At(geotile.id), date_ind_center, height_ind_center] = refreeze
    t_air_hyps[At(geotile.id), date_ind_center, height_ind_center] = t_air
    rain_hyps[At(geotile.id), date_ind_center, height_ind_center] = rain
    fac_to_depth_hyps[At(geotile.id), date_ind_center, height_ind_center] = fac_to_depth

    t2 = time()
    println("gemb binned: $(geotile.id), time = $(round(t2-t1))s")
end

out_var = Dict(
    "fac_hyps" => fac_hyps, 
    "smb_hyps" => smb_hyps, 
    "nobs_hyps" => nobs_hyps,
    "ec_hyps" => ec_hyps,
    "acc_hyps" => acc_hyps,
    "runoff_hyps" => runoff_hyps,
    "melt_hyps" =>  melt_hyps,
    "refreeze_hyps" => refreeze_hyps,
    "t_air_hyps" => t_air_hyps,
    "rain_hyps" => rain_hyps,
    "fac_to_depth_hyps" => fac_to_depth_hyps,
)

save(out_file,out_var);