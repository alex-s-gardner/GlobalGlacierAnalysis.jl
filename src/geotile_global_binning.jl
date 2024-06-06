# load packages

begin
    using Arrow
    using Altim
    using DataFrames
    using Extents
    using GeoInterface
    using Shapefile
    using Rasters
    using DimensionalData
    using BinStatistics
    using Statistics
    using Dates
    using CairoMakie
    using JLD2
    using FileIO
    using LsqFit

    # run parameters
    force_remake = false;
    project_id = :v01;
    geotile_width = 2;
    warnings = false
    showplots = false;
    mask = :land

    dem_ids = [:best, :cop30_v2]
    binning_methods = ["median", "meanmadnorm3", "meanmadnorm5"];
    curvature_corrects = [true, false]

    paths = project_paths(; project_id)
    products = project_products(; project_id)
    binned_folder = analysis_paths(; geotile_width).binned
    fig_folder = joinpath(binned_folder, "figures")

    # open shapefiles
    excludemask_flag = false
    if mask == :land
        shp = Symbol("$(:water)_shp");
        fn_shp = Altim.pathlocal[shp];
        feature = Shapefile.Handle(fn_shp);

        invert = true

        shp = Symbol("$(:landice)_shp")
        fn_shp = Altim.pathlocal[shp]
        excludefeature = Shapefile.Handle(fn_shp)
        excludemask_flag = true
    else
        shp = Symbol("$(mask)_shp");
        invert = false
    end


    # define model for curvature correction
    model1::Function = model1(c, p) = p[1] .+ p[2] .* c .+ p[3] .* c .^ 2 .+ p[4] .* c .^ 3 .+ p[5] .* c .^ 4
    p1 = zeros(5)

    #### DON NOT CHANGE THESE PARAMETERS
    max_canopy_height = 1 # do not change

    # filter parameters
    filt = (
        dh_max=200,
    )

    # load geotile definitions with corresponding hypsometry
    gt_file = joinpath(binned_folder, "geotile_$(mask)_hyps.arrow")
    geotiles = DataFrame(Arrow.Table(gt_file))
    geotiles.extent = Extent.(getindex.(geotiles.extent, 1))

    # filter geotiles
    geotiles = geotiles[(geotiles.glacier_frac.>0.0), :]
    #geotiles = geotiles[(geotiles.rgi11.>0.0), :]

    # define date and hight binning ranges 
    dd =  30;
    date_range = DateTime(1990):Day(dd):now();
    date_center = date_range[1:end-1] .+ Day(dd/2);

    Δc = 0.1;
    curvature_range = -1.0:Δc:1.0;
    curvature_center = curvature_range[1:end-1] .+ Δc/2;

    # DO NOT MODIFY: these need to match elevations in hypsometry
    Δh = 100;
    height_range = 0:100:10000;
    height_center = height_range[1:end-1] .+ Δh/2;
end;

#for binning_method = binning_methods
    #for dem_id in dem_ids  # [:rema_v2_10m]#[:cop30_v2]; #, :nasadem_v1, :rema_v2_10m, :arcticdem_v4_10m, :best]
        #for curvature_correct in curvature_corrects
begin
            if  true
                binning_method = "meanmadnorm3"
                dem_id = :best 
                curvature_correct = true
                force_remake = true
            end

            # funtion used for binning data
            if binning_method == "meanmadnorm3"
                binningfun(x) = mean(x[Altim.madnorm(x).<3])
            elseif binning_method == "meanmadnorm5"
                binningfun(x) = mean(x[Altim.madnorm(x).<5])
            elseif binning_method == "median"
                binningfun(x) = median(x)
            else
                error("unrecognized binning method")
            end

            # bin method
            if curvature_correct
                runid = "$(mask)_dh_$(dem_id)_cc_$(binning_method)_$(project_id)"
            else
                runid = "$(mask)_dh_$(dem_id)_$(binning_method)_$(project_id)"
            end

            out_file = joinpath(binned_folder, "$(runid).jld2");

            @time if !isfile(out_file) || force_remake
                # 6.2 hours for all glaciers, all missions/datasets on 96 threads
                # 10hr for land for all glaciers, all missions/datasets on 96 threads

                # initialize dimensional arrays
                dh_hyps = Dict();
                nobs_hyps = Dict();
                curvature_hyps = Dict()
                ngeotile = nrow(geotiles)
                ndate = length(date_center)
                ncurvature = length(curvature_center);
                for product in products
                    mission = String(product.mission)
                    push!(dh_hyps, String(mission) => DimArray(fill(NaN, ngeotile, ndate, length(height_center)), (geotile = geotiles.id, date=date_center, height = height_center)));
                    push!(nobs_hyps, String(mission) => DimArray(fill(0, ngeotile, ndate, length(height_center)), (geotile = geotiles.id, date=date_center, height = height_center)));
                    push!(curvature_hyps, String(mission) => DataFrame(id=geotiles.id, curvature=[curvature_range for i in 1:ngeotile], dh = [fill(NaN,ncurvature) for i in 1:ngeotile], nobs=[zeros(ncurvature) for i in 1:ngeotile], model_coef = [zeros(size(p1)) for i in 1:ngeotile]))
                end

                # 2.6 hours for all 4 missions for all glacierized geotiles
                # 31 min for icesat + iceast2 + gedi for all glacierized geotiles
                @time for mission in keys(dh_hyps)
                    # <><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><><><><>
                    # mission = "icesat"
                    #out_file = "/mnt/devon2-r1/devon0/gardnera/Documents/GitHub/Altim/data/geotiles_glacier_hyps_2deg_dh.arrow";
                    #geotiles  = DataFrame(Arrow.Table(out_file));
                    #geotiles.extent = Extent.(getindex.(geotiles.extent, 1));
                    #geotiles = geotiles[(geotiles.rgi19 .> 0) .& (geotiles.glacier_frac .> 0), :] ;
                    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
                    
                    product = products[Symbol(mission)];

                    Threads.@threads for geotile in eachrow(geotiles)
                        # <><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><><>
                        #for geotile in eachrow(geotiles)
                        # geotile = geotiles[findfirst((geotiles.rgi2 .> 0.) .& (geotiles.glacier_frac .> 0.1)),:];
                        #geotile = geotiles[findfirst(geotiles.id .== "lat[+44+46]lon[+006+008]"), :]
                        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
                        
                        if geotile.glacier_frac == 0.0
                            continue
                        end

                        t1 = time();
                        path2altim = joinpath(paths[Symbol(mission)].geotile, geotile.id * ".arrow");

                        if !isfile(path2altim)
                            continue
                        end

                        altim = select!(DataFrame(Arrow.Table(path2altim)), [:longitude, :latitude, :datetime, :height, :quality]);

                        # add dem height and curvature
                        if dem_id == :best
                            # last dem takes precidence over earlier dems
                            dem_id0 = [:cop30_v2, :arcticdem_v4_10m, :rema_v2_10m]
                        else
                            dem_id0 = [dem_id]
                        end;

                        altim[!, :height_ref] .= NaN;
                        altim[!, :curvature] .= NaN;

                        for dem_id1 in dem_id0           
                            path2dem = joinpath(paths[Symbol(mission)].geotile, geotile.id * ".$(dem_id1)")
                            if isfile(path2dem)
                                dem = select!(DataFrame(Arrow.Table(path2dem)), :height, :dhddx, :dhddy)
                                # a bit faster to calculate curvature on all data then subset
                                curv = Altim.curvature.(dem.dhddx, dem.dhddy, Ref(Altim.dem_info(dem_id1)[1].epsg), lat = mean(geotile.extent.Y))
                                
                                ind = .!isnan.(curv) .& (abs.(dem.height) .< 9998)
                                
                                altim[ind, :height_ref] = dem[ind, :height]
                                altim[ind, :curvature] = curv[ind]
                            end
                        end;

                        altim[!, :dh] = altim.height .- altim.height_ref;

                        # load masks
                        path2masks = joinpath(paths[Symbol(mission)].geotile, geotile.id * ".masks");
                        masks = select!(DataFrame(Arrow.Table(path2masks)), [:inlandwater, :land, :landice, :ocean]);

                        if minimum(geotile.extent.Y) < -54 || maximum(geotile.extent.Y) > 66
                            canopy = DataFrame(canopyh = zeros(size(altim.latitude)))
                        else
                            path2canopy = joinpath(paths[Symbol(mission)].geotile, geotile.id * ".canopyh")
                            canopy = DataFrame(Arrow.Table(path2canopy));
                        end;

                        # update mask with high-resolution vector files
                        grid_resolution = 0.00027; # ~30m

                        x_mask = X(geotile.extent.X[1]:grid_resolution:geotile.extent.X[2], 
                            sampling=DimensionalData.Intervals(DimensionalData.Start()));
                        y_mask = Y(geotile.extent.Y[1]:grid_resolution:geotile.extent.Y[2], 
                            sampling=DimensionalData.Intervals(DimensionalData.Start()));

                        mask0 = Raster(zeros(UInt8, y_mask, x_mask));
                        setcrs(mask0, EPSG(4326));

                        # NOTE: count method is fastest
                        mask0 = Rasters.rasterize!(count, mask0, feature; threaded=false, 
                            shape=:polygon, progress=false, verbose=false, boundary=:center) .> 0;

                        if invert
                            mask0 = .!(mask0)
                        end

                        if excludemask_flag
                            excludemask = Raster(zeros(UInt8, y_mask, x_mask));
                            excludemask = Rasters.rasterize!(count, excludemask, excludefeature; threaded=false, shape=:polygon, progress=false, verbose=false) .> 0
                            mask0 = mask0 .& .!excludemask
                        end

                        isvalid = Altim.within.(Ref(geotile.extent), altim.longitude, altim.latitude);
                        masks[!, mask] .= false;

                        fast_index = true;
                        if fast_index # fast index is 15x faster than Rasters
                            c = floor.(Int64, (altim.longitude[isvalid] .- first(x_mask)) ./ step(x_mask)) .+ 1;
                            r = floor.(Int64, (altim.latitude[isvalid] .- first(y_mask)) ./ step(y_mask)) .+ 1;
                            pts3 = [CartesianIndex(a, b) for (a, b) in zip(r, c)]
                            masks[isvalid,mask] .= mask0[pts3]
                        else
                            #NOTE: 67% of time is taken here for large number of points.
                            pts1 = GeoInterface.PointTuple.([((Y=y, X=x)) for (x, y) in 
                                zip(altim.longitude[isvalid], altim.latitude[isvalid])]);
                            pts2 = extract(mask0, pts1, atol=grid_resolution/2, index = true, geometry=false);
                            masks[:, mask][isvalid] = getindex.(pts2, 2)
                            pts3 = CartesianIndex.(getindex.(pts2, 1))
                        end;

                        if !any(masks[:, mask])
                            continue
                        end

                        ## Check slope relation
                        if false
                            for dem_id1 in dem_id0           
                                path2dem = joinpath(paths[Symbol(mission)].geotile, geotile.id * ".$(dem_id1)")
                                if isfile(path2dem)
                                    dem = select!(DataFrame(Arrow.Table(path2dem)), :height, :dhdx, :dhdy, :dhddx, :dhddy)
                    
                                    dhdx_dhdy = Altim.slope.(dem.dhdx, dem.dhdy, Ref(Altim.dem_info(dem_id1)[1].epsg), lat = mean(geotile.extent.Y));
                                    
                                    var_ind = (.!masks[:, mask]) .& isvalid .& .!isnan.(altim.dh) .& (altim.height .!== 0) .& (altim.height_ref .!== 0) .& .!isnan.(altim.curvature) .& (abs.(altim.dh) .< 10) .& masks.land .& ((canopy.canopyh .<= max_canopy_height) .| (altim.latitude .< -60))
                                    
                                    if sum(var_ind) <= length(p1)
                                        warnings && (@warn ("$(geotile.id): slope offset not calculated, not enought off-ice points"))
                                    else
                                        offset = Altim.track_offset(getindex.(dhdx_dhdy[var_ind], 1), getindex.(dhdx_dhdy[var_ind], 2), altim.dh[var_ind]);
                                    end
                                end
                            end
                        end;

                        # use bin statistics to calculate dh vs. elevation for every 3 months.

                        # identify valid dh data
                        showplots && plot(Altim.decimalyear.(altim.datetime[isvalid]),altim.dh[isvalid]; seriestype=:scatter, label="all")

                        if curvature_correct
                            var_ind = (.!masks[:, :landice]) .& isvalid .& .!isnan.(altim.dh) .& (altim.height .!== 0) .& (altim.height_ref .!== 0) .& .!isnan.(altim.curvature) .& (abs.(altim.dh) .< 10) .& masks.land .& ((canopy.canopyh .<= max_canopy_height) .| (altim.latitude .< -60))

                            if sum(var_ind) <= length(p1)
                                warnings && (@warn ("$(geotile.id): no curvature corecction applied, not enought off-ice points"))
                            else
                                # check bounds: binstats will throw an error if no data is passed to median()
                                if Altim.vector_overlap(altim[var_ind, :curvature], curvature_range)
                                    df = DataFrame();
                                    try
                                        df = binstats(altim[var_ind, :], [:curvature], [curvature_range], :dh; col_function=[median], missing_bins=true);
                                    catch
                                        # having issues with edge cases
                                        println(" ##################### ERRORED ##########################")
                                        println(" ##################### $(mission): $(geotile.id) ########################")
                                        println(" ##################### ERRORED ##########################")
                                    end

                                    if !isempty(df)
                                        bin_center = BinStatistics.bincenter.(df.curvature);
                                        bin_valid = .!(ismissing.(df.dh_median))

                                        if (sum(bin_valid) <= length(p1))
                                            warnings && (@warn ("$(geotile.id): no curvature corecction applied, not enought off-ice points"))
                                        else
                                            fit1 = curve_fit(model1, bin_center[bin_valid], df.dh_median[bin_valid], df.nrow[bin_valid], p1);
                                            #histogram(altim.dh[var_ind], bins= -5:.1:5; label ="raw", xlabel="height change [m]")
                                            altim.dh = altim.dh .- model1(altim.curvature, fit1.param);
                                            #histogram!(altim.dh[var_ind], bins= -5:.1:5; label="curvature corrected")

                                            gt_ind = findfirst(geotile.id .== curvature_hyps[mission].id)
                                            cdh = @view curvature_hyps[mission].dh[gt_ind] 
                                            cdh[1][bin_valid] = df.dh_median[bin_valid];
                                            cnobs = @view curvature_hyps[mission][gt_ind, :nobs]
                                            cnobs[1][bin_valid] = df.nrow[bin_valid]
                                            curvature_hyps[mission][gt_ind, :model_coef] = fit1.param
                                
                                            if showplots
                                                dh_cor = model1(bin_center, fit1.param);
                                                p = plot(bin_center, df.dh_median);

                                                fontsize = 18;

                                                title = "$mission: $(geotile.id)"

                                                title = replace.(title, "hugonnet" => "ASTER")
                                                title = replace.(title, "icesat" => "ICESat")
                                                title = replace.(title, "icesat2" => "ICESat-2")
                                                title = replace.(title, "gedi" => "GEDI")
                                                
                                                f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98), size=(1000, 700), fontsize=fontsize);
                                                ga = f[1:4, 1] = GridLayout()
                                                gb = f[5:6, 1] = GridLayout()
                                
                                                Label(ga[1, 1, Top()],
                                                    title, valign=:bottom,
                                                    font=:bold,
                                                    padding=(0, 0, 5, 0)
                                                )

                                                axmain = Axis(ga[1, 1]; ylabel="height anomaly [m]")
                                                hidexdecorations!(axmain; grid=false, minorgrid=false)
                                                axbottom = Axis(gb[1, 1], ylabel="count [×1000]", xlabel="curvature [cm⁻¹]", )
                                                plot!(axmain, bin_center, df.dh_median; label = "observation");
                                                lines!(axmain, bin_center, dh_cor; label = "model");
                                                plot!(axmain, bin_center, df.dh_median .- dh_cor; label = "corrected");
                                                barplot!(axbottom, bin_center, df.nrow / 1000)
                                                axislegend(axmain, framevisible=false, position = :lt)
                                                

                                                fname = joinpath(fig_folder, "$(mask)_$(mission)_$(geotile.id)_$(dem_id)_curvature.png")
                                                save(fname, f)
                                                display(f)
                                            end
                                        end
                                    end
                                end
                            end
                        end

                        #################################### FILTER 1 ######################################
                        var_ind = masks[:, mask] .& isvalid .& .!isnan.(altim.dh) .& (altim.height .!== 0) .& (altim.height_ref .!== 0)
                        if sum(var_ind) < 100
                            continue
                        end

                        # for troubleshooting
                        showplots && plot!(altim.datetime[var_ind],altim.dh[var_ind]; seriestype=:scatter,label="glacier")

                        var_ind[var_ind] = (abs.(altim.dh[var_ind] .- median(altim.dh[var_ind])) .< filt.dh_max) .| (abs.(altim.dh[var_ind]) .> filt.dh_max*2)

                        if mask == :land
                            var_ind = var_ind .& (canopy.canopyh .<= max_canopy_height)
                        end
                        
                        # for troubleshooting
                        showplots && plot!(altim.datetime[var_ind],altim.dh[var_ind]; seriestype=:scatter, label="glacier-filt", ylims = (-300, +300))
                        
                        if product.apply_quality_filter
                            var_ind = var_ind .& altim.quality

                            # for troubleshooting 
                            showplots && plot!(altim.datetime[var_ind],altim.dh[var_ind]; seriestype=:scatter, ylims = (-100, +100))
                        end
                        ####################################################################################

                        if !any(var_ind)
                            continue
                        end

                        # bin data by date and elevation
                        minmax_date = extrema(altim.datetime[var_ind])
                        minmax_height = extrema(altim.height_ref[var_ind])
                        date_ind = (date_range .>= minmax_date[1]-Day(dd)) .& 
                            (date_range .<= (minmax_date[2]+Day(dd)))
                        date_ind_center = findall(date_ind)[1:end-1];
                        height_ind = (height_range .>= minmax_height[1]-Δh) .& 
                            (height_range .<= (minmax_height[2]+Δh))
                        height_ind_center = findall(height_ind)[1:end-1];
                        
                        # check bounds: binstats will throw an error if no data is passed to median()
                        if !Altim.vector_overlap(altim[var_ind, :datetime], date_range[date_ind]) || 
                            !Altim.vector_overlap(altim[var_ind, :height_ref], height_range[height_ind]) 
                            continue
                        end

                        df = DataFrame()
                        try
                            df = binstats(altim[var_ind, :], [:datetime, :height_ref], 
                                [date_range[date_ind], height_range[height_ind]], 
                                :dh; col_function=[binningfun], missing_bins=true)
                        catch
                            # having issues with edge cases
                            println(" ##################### ERRORED ##########################")
                            println(" ##################### $(geotile.id) ########################")
                            println(" ##################### ERRORED ##########################")
                            continue
                        end

                        gdf = DataFrames.groupby(df, :datetime)

                        # create an array of dh as a function of time and elevation
                        obs1 = fill(NaN, sum(date_ind)-1, sum(height_ind)-1)
                        nobs1 = fill(Int64(0), sum(date_ind)-1, sum(height_ind)-1)    
                        p = sortperm(BinStatistics.bincenter.(gdf[1].height_ref))

                        h_center = bincenter.(gdf[1].height_ref)[p];
                        t_center = date_center[date_ind_center]

                        for (i,df) in enumerate(gdf)
                            isval = .!ismissing.(df[p, "dh_binningfun"])
                            obs2 = @view obs1[i,:]
                            nobs2 = @view nobs1[i, :]
                            if any(isval)
                                obs2[isval] = df.dh_binningfun[p][isval]
                                nobs2[isval] = df.nrow[p][isval]
                            end
                        end

                        # for troubleshooting 
                        showplots && heatmap(h_center, t_center, obs, clim = (-10, 10))

                        dh_hyps[mission][At(geotile.id), date_ind_center, height_ind_center] = obs1
                        nobs_hyps[mission][At(geotile.id), date_ind_center, height_ind_center] = nobs1
                        t2 = time()
                        dt = round(Int16, t2 - t1);
                        println("binned: $(mission) - $(geotile.id): $(dt)s")
                    end
                end
                
                save(out_file, Dict("dh_hyps" => dh_hyps, "nobs_hyps" => nobs_hyps, "curvature_hyps" => curvature_hyps));
            end
        end
    end
end