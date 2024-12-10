
using Altim;
using Arrow;
using DataFrames;
using CairoMakie;
using GaussianProcesses;
using Statistics;
using Interpolations;
using BinStatistics;
import FastGeoProjections as FGP;
using Extents
using GeoArrays
using LazyGrids
using Rasters
import GeoFormatTypes as GFT
using Shapefile
using GeoInterface
using GeometryOps

begin
    folder_height_change = "/mnt/bylot-r3/data/height_change/2000_2022";
    pathlocal = setpaths()

    force_remake = true; 
    debug = true
    geotile_width = 2; #geotile width [degrees]
    grid = (node_spacing=50, node_width=100); # [node_spacing=50, node_width=100]
    trend_max_mad = 15;
    length_scale_xy = 500; # 500 works well
    length_scale_z = 100; # 100 works well

    project_id = :v01;
    domain = :landice;
    #domain = :all;

    paths = project_paths(project_id = project_id);
    # add height_change_path
    paths = (height_change = folder_height_change, paths...);
    products = project_products(; project_id=:v01);
    geotiles = project_geotiles(; geotile_width = geotile_width, domain = domain);

    printstyled("calculating height change \n"; color = :blue, bold = true);

    #@time asyncmap(eachrow(geotiles[1:20,:]); ntasks = 100) do geotile 
    dems = ["cop30_v2", "arcticdem_v3_10m", "cop30_v2", "rema_v2_10m"];

    # --------------------------------------------------------------------------
    @warn "--- only processing a subset of geotiles ----"
    region = :RGI01;
    ext, epsg = region_extent(region);

    dems = [:cop30_v2]; 

    #ext = Extent(X=(-126.9, -126.1), Y=(51.1, 51.8));
    geotiles = geotile_subset!(geotiles, ext);
end

t1 = time();

geotiles = geotiles[findfirst("lat[+60+62]lon[-144-142]" .== geotiles.id):end,:]
for geotile in eachrow(geotiles);
#for k = 28:nrow(geotiles)

    #geotile = first(geotiles[geotiles.id .== "lat[+60+62]lon[-142-140]",:])
    #geotile = geotiles[k,:]
    t2 = time();
    dem = first(dems);

    # deteremine buffer size in decimal degrees
    buff_x = (length_scale_xy * 5) / (40075000 * cosd(abs((geotile.extent.min_y .+ geotile.extent.max_y) / 2)) / 360)
    buff_y = (length_scale_xy * 5) / 111111 
    buffer=(X=buff_x, Y=buff_y)
    df = geotile_read(paths.height_change, geotile.id, dem; buffer=buffer)
    
    if isempty(df)
        @warn "$(geotile.id) no heigh change data, skipping"
        continue
    end

    ext = Altim.nt2extent(geotile.extent)
    ext0 = Extents.buffer(ext, buffer)

    df, epsg = itslive_proj!(df);

    # TODO: read data [not sure why height was not properly saved]????
    df.height = vec(dem_height(df.longitude, df.latitude, :cop30_v2));

    ##TODO: save file for working local 
    #Arrow.write("test_data.arrow", df)

    # create grid with heights
    minmax_x = extrema(df.X);
    minmax_y = extrema(df.Y);

    center_extent = (
        x_min=floor(minmax_x[1], digits=-4),
        y_min=floor(minmax_y[1], digits=-4),
        x_max=ceil(minmax_x[2], digits=-4),
        y_max=ceil(minmax_y[2], digits=-4)
    );

    griddef = regular_grid(center_extent, grid.node_spacing, node_width=grid.node_width);

    x = X(griddef.x_node_center);
    y = Y(griddef.y_node_center);

    shp = allfiles(pathlocal.RGI_dissolved; subfolders=true, fn_endswith=".shp")

    A = Raster(zeros(UInt8,x,y); crs = epsg);
    M = copy(A)
    for s in shp
        polygons = Shapefile.Handle(s)
        ext = GeoInterface.bbox(polygons)

        if !Extents.intersects(ext0, ext)
            # no overlap
            continue
        else
            polygons = GeometryOps.reproject(polygons; target_crs=epsg)
            M = M .| (Rasters.rasterize!(count, A, polygons.geom; threaded=false, shape=:polygon) .> 0);
        end
    end

    if debug
        heatmap(M)
    end

    if sum(M) < 10
        @warn "$(geotile.id) does not have ice area, skipping"
        continue
    end

    trans = FGP.Transformation("EPSG:$epsg", "EPSG:4326", always_xy=true);
    lon, lat = trans(collect(vec(griddef.node_center[2])), collect(vec(griddef.node_center[1])));

    # only include lon and lat values within geotile.extent
    ext = Altim.nt2extent(geotile.extent);
    buffX = (grid.node_spacing) / (40075000 * cosd(abs((geotile.extent.min_y .+ geotile.extent.max_y) / 2)) / 360)
    buffY = (grid.node_spacing) / 111139
    ext_buff = Extent((X = ext.X .+ (-buffX, +buffX), Y = ext.Y .+ (-buffY, +buffY)))
    M2 = reshape(Altim.within.(Ref(ext_buff), lon, lat), size(griddef.node_center[2]))
    M0 = (M' .> 0) .& M2

    if debug
        heatmap(M0)
    end

    # transpose M to match griddef definition 
    x = griddef.node_center[2][M0]
    y = griddef.node_center[1][M0]

    # this take about 3.6 seconds
    # transpose M to match griddef definition 
    h = vec(dem_height(lon[M0[:]], lat[M0[:]], :cop30_v2))
    df_grid = DataFrame(; x, y, h);

    ##TODO: save file for working local 
    # Arrow.write("test_grid.arrow", df_grid)

    # bin as a funciton of elevation 
    h_bins = 0:100:9000;
    df0 = binstats(df, :height, h_bins, :trend, col_function=[median]);

    valid = df0.nrow .>= 11;

    # if there are no valid bins then move to next geotile
    if !any(valid)
        @warn "$(geotile.id) does not have enough data, skipping"
        continue
    end

    if debug
        plot(BinStatistics.bincenter.(df0.height[valid]),df0.trend_median[valid])
    end
    
    # detrend data
    nodes = BinStatistics.bincenter.(df0.height[valid]);
    A = df0.trend_median[valid];

    #TODO: This sorting function should not be needed... need to look at binstats.jl
    p = sortperm(nodes);
    nodes = nodes[p];
    A = A[p];

    extrapolation_bc = Linear();

     # if there are no valid bins then move to next geotile
    if length(nodes) < 2
        @warn "$(geotile.id) does not have enough data for hypsometric interpolation, skipping"
        continue
    end

    itp = extrapolate(Interpolations.interpolate((nodes,), A, Gridded(Linear())), Flat());
    df[!, :trend_anom] = df.trend .- itp(Float64.(df.height));

    # m = scatter(df.height, df.trend_anom)

    #Select mean and covariance function
    mZero = MeanZero();

    kern = Matern(5 / 2, [log(length_scale_xy), log(length_scale_xy), log(length_scale_z)], log(0.1)); #Matern 5/2 ARD kernel
    idx = .!isnan.(df.trend_anom) .& (df.mad .< trend_max_mad);

    X0 = vcat(df.X[idx]', df.Y[idx]', df.height[idx]');
    dhdt_anom_foo = df.trend_anom[idx]';

    dhdt_anom_foo = vec(Float64.(dhdt_anom_foo));
    Xp = vcat(df_grid.x', df_grid.y', df_grid.h');
    Yp = zeros(length(Xp[1,:]));
    Sp = zeros(length(Xp[1,:]));

#=
    node_spacing = 2000;
    node_width = node_spacing + (length_scale_xy * 4);

    overlap = (node_width - node_spacing)/2;
    xext = extrema(df_grid.x)
    yext = extrema(df_grid.y)

    center_extent = (
        x_min=(xext[1] + overlap),
        y_min=(yext[1] + overlap),
        x_max=(xext[1] + overlap) + max(0,ceil((xext[2] - xext[1] - 2 * overlap) / node_spacing) * node_spacing),
        y_max=(yext[1] + overlap) + max(0, ceil((yext[2] - yext[1] - 2 * overlap) / node_spacing) * node_spacing),
    )

    chunking = regular_grid(center_extent, node_spacing, node_width=node_width)

    

    node_centers = hcat(chunking.node_center[1][:], chunking.node_center[2][:]);
    @time Threads.@threads for xy in eachrow(node_centers)
    #for xy in zip(chunking.node_center[1], chunking.node_center[2])
        x_min = xy[2] - chunking.node_half_width
        x_max = xy[2] + chunking.node_half_width
        y_min = xy[1] - chunking.node_half_width
        y_max = xy[1] + chunking.node_half_width

        ind1 = (X0[1, :] .> x_min) .& (X0[1, :] .<= x_max) .& (X0[2, :] .> y_min) .& (X0[2, :] .<= y_max)
        ind2 = (Xp[1, :] .> x_min) .& (Xp[1, :] .<= x_max) .& (Xp[2, :] .> y_min) .& (Xp[2, :] .<= y_max)

        if ~any(ind1) || ~any(ind2)
            continue
        end
        
        @time gp = GaussianProcesses.SoR(X0[:,ind1], Xp[:,ind2], dhdt_anom_foo[ind1], mZero, kern, log(0.1))
        #gp = GP(X0[:,ind1], dhdt_anom_foo[ind1], mZero, kern, log(0.24))
        # optimize!(gp)
        Yp[ind2], Sp[ind2] = predict_y(gp, Xp[:,ind2])
    end


    =#
    # this takes 1-hr
    Threads.@threads for i in eachindex(Xp[1,:])
        x_min = Xp[1,i] - length_scale_xy * 3;
        x_max = Xp[1,i] + length_scale_xy * 3;
        y_min = Xp[2,i] - length_scale_xy * 3;
        y_max = Xp[2,i] + length_scale_xy * 3;

        ind1 = (X0[1,:] .> x_min) .& (X0[1,:] .<= x_max) .& (X0[2,:] .> y_min) .& (X0[2,:] .<= y_max);

        if ~any(ind1)
            continue
        end
        
        dhdt_anom_mean = mean(dhdt_anom_foo[ind1]);


        gp = GaussianProcesses.SoR(X0[:, ind1], reshape(Xp[:, i], size(Xp, 1), 1), dhdt_anom_foo[ind1] .- dhdt_anom_mean, mZero, kern, log(0.1));

        #gp = GP(X0[:, ind1], dhdt_anom_foo[ind1] .- dhdt_anom_mean, mZero, kern, log(0.1))
        # optimize!(gp)
        a, b = predict_y(gp, reshape(Xp[:, i], size(Xp, 1), 1));
        Yp[i] = a[1] + dhdt_anom_mean;
        Sp[i] = b[1];
    end

    #Yp[Yp .== 0] .= NaN;
    Yp = (Yp .+ itp(Float64.(vec(df_grid.h))));

    Sp[Sp.==0] .= NaN;
    

    colorrange = (-5, 5);
    cmap = :balance;

    #fig = Figure()
    #sp = scatter(df_grid.x, df_grid.y; color=Yp, colormap=cmap, markersize=1, strokewidth=0, colorrange=(-1, 1))

    #sp = scatter(df_grid.x, df_grid.y; color=Yp, colormap=cmap, markersize=1, strokewidth=0, colorrange=colorrange)

    # save data as a geotiff

    x = X(griddef.x_node_center);
    y = Y(griddef.y_node_center);
    A = Raster(fill(NaN, x, y); crs=epsg, missingval=NaN);

    (A')[M0] = Yp;

    dhdt_file = joinpath(paths.height_change, "$(geotile.id).$(dem)");

    Rasters.write(dhdt_file * ".tif", A; driver="COG", options=Dict("BLOCKSIZE" => "512", "RESAMPLING"=>"NEAREST"), force=true);

    (A')[M0] = Sp;
    Rasters.write(dhdt_file * "_err.tif", A; driver="COG", options=Dict("BLOCKSIZE" => "512", "RESAMPLING" => "NEAREST"), force=true);

    t = round(time() - t2);

    printstyled("    -> $(geotile.id) GP ineterp complete [$t s]\n"; color=:light_black);
end
t = round((time() - t1)/60);
printstyled("    -> GP ineterp complete [$t min]\n"; color=:light_black)