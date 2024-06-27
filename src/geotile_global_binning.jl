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
    using GeoTiles

    # run parameters
    force_remake = false;
    update_geotile = false; # this will load in prevous results to update select geotiles or missions
    project_id = :v01;
    geotile_width = 2;
    warnings = false
    showplots = false;

    dem_ids = [:best, :cop30_v2]
    binning_methods = ["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"]
    curvature_corrects = [true, false]
    surface_masks = [:glacier] #:glacier_b1km, :land, :glacier_b10km
    
    if true
        dem_ids = [:best]
        binning_methods = ["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"]
        curvature_corrects = [true]
        surface_masks = [:glacier_b1km, :land, :glacier_b10km]
    end

    paths = project_paths(; project_id)
    products = project_products(; project_id)
    binned_folder = analysis_paths(; geotile_width).binned
    fig_folder = joinpath(binned_folder, "figures")
        
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
    geotiles = Altim.geotiles_w_mask(geotile_width)
    #gt_file = joinpath(binned_folder, "geotile_$(:glacier)_hyps.arrow")
    #geotiles = DataFrame(Arrow.Table(gt_file))
    #geotiles.extent = Extent.(getindex.(geotiles.extent, 1))

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

## THIS WILL REMOVE ALL _dh_ files if you need to rebuild binned archive
#files2delete = allfiles("/mnt/bylot-r3/data/binned/2deg"; fn_contains="_dh_")
#rm.(files2delete)
##########

for surface_mask in surface_masks
    #surface_mask = :glacier
    # open shapefiles
    if surface_mask == :land
        shp = Symbol("$(:water)_shp");
        fn_shp = Altim.pathlocal[shp];
        feature = Shapefile.Handle(fn_shp);
        invert = true

        shp = Symbol("$(:landice)_shp")
        fn_shp = Altim.pathlocal[shp]
        excludefeature = Shapefile.Handle(fn_shp)
    else
        shp = Symbol("$(surface_mask)_shp");
        fn_shp = Altim.pathlocal[shp]
        feature = Shapefile.Handle(fn_shp)
        invert = false

        excludefeature = nothing;
    end

    for binning_method = binning_methods
        for dem_id in dem_ids  
            for curvature_correct in curvature_corrects
            
                if  false
                    binning_method = "meanmadnorm10"
                    dem_id = :cop30_v2
                    curvature_correct = true
                    force_remake = true
                end

                # funtion used for binning data
                binningfun = Altim.binningfun_define(binning_method)
            
                # bin method
                if curvature_correct
                    runid = "$(surface_mask)_dh_$(dem_id)_cc_$(binning_method)_$(project_id)"
                else
                    runid = "$(surface_mask)_dh_$(dem_id)_$(binning_method)_$(project_id)"
                end

                out_file = joinpath(binned_folder, "$(runid).jld2");

                # <><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><>
                # binningfun(x) = mean(x[Altim.madnorm(x).<10])
                # <><><><><><><><><><><><><><><><><><><><><><<><><><><><><><><><><><><><><><><><>
                if !isfile(out_file) || force_remake
                    # 6.2 hours for all glaciers, all missions/datasets on 96 threads
                    # 10hr for land for all glaciers, all missions/datasets on 96 threads

                    # initialize dimensional arrays
                    # update_geotile = true
                    if update_geotile
                        # load exisiting
                        dh_hyps = load(out_file, "dh_hyps")
                        nobs_hyps = load(out_file, "nobs_hyps")
                        curvature_hyps = load(out_file, "curvature_hyps")
                    else
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
                    end

                    # 2.6 hours for all 4 missions for all glacierized geotiles
                    # 31 min for icesat + iceast2 + gedi for all glacierized geotiles
                    for mission in keys(dh_hyps)
                        # <><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><><><><>
                        # mission = "hugonnet"
                        #out_file = "/mnt/devon2-r1/devon0/gardnera/Documents/GitHub/Altim/data/geotiles_glacier_hyps_2deg_dh.arrow";
                        #geotiles  = DataFrame(Arrow.Table(out_file));
                        #geotiles.extent = Extent.(getindex.(geotiles.extent, 1));
                        #geotiles = geotiles[(geotiles.rgi2 .> 0) .& (geotiles.glacier_frac .> 0), :] ;
                        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
                        
                        product = products[Symbol(mission)];

                        #idx = (geotiles.rgi9 .> 0) .& (geotiles.glacier_frac .> 0.1)
                        #Threads.@threads 
                        Threads.@threads for geotile in eachrow(geotiles)
                            # <><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><><>
                            #for geotile in eachrow(geotiles[idx,:])
                            # geotile = geotiles[findfirst((geotiles.rgi2 .> 0.) .& (geotiles.glacier_frac .> 0.1)),:];
                            #geotile = geotiles[findfirst(geotiles.id .== "lat[-70-68]lon[+078+080]"), :]
                            #geotile = geotiles[findall(GeoTiles.within.(Ref(79.5), Ref(25.5), geotiles.extent)),:]
                            # Threads.@threads for geotile in eachrow(geotiles[Extents.intersects.(Ref(Extent(X=(23.6, 25.8), Y=(79.2, 79.7))), geotiles.extent), :])

                            #geotile = geotiles[findfirst(geotiles.id .== "lat[-72-70]lon[+160+162]"), :]

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

                            # add height_ref and curvature
                            mission_geotile_folder = paths[Symbol(mission)].geotile
                            altim = Altim.add_dem_ref!(altim, dem_id, geotile, mission_geotile_folder)
                            
                            # load masks
                            path2masks = joinpath(paths[Symbol(mission)].geotile, geotile.id * ".masks");
                            masks0 = select!(DataFrame(Arrow.Table(path2masks)), [:inlandwater, :land, :landice, :ocean]);

                            # add high resolution mask 
                            masks0 = Altim.highres_mask!(masks0, altim, geotile, feature, invert, excludefeature, surface_mask)
                            valid = Altim.within.(Ref(geotile.extent), altim.longitude, altim.latitude)

                            # quick trouble shoot check by plotting time series for each geotile
                            if false
                                valid = masks0.landice .& .!isnan.(altim.dh) 
                                valid[valid] = (abs.(altim.dh[valid] .- median(altim.dh[valid])) .< filt.dh_max) .& (abs.(altim.dh[valid]) .< filt.dh_max*2)
                                #valid[valid] = (abs.(altim.dh[valid] .- median(altim.dh[valid])) .< 2000) .& (altim.dh[valid] .!== 0) .& (altim.height_ref[valid] .> 0)
                                if sum(valid) < 100
                                    #continue
                                end

                                #valid = valid .& altim.quality
                            
                            
                                #binningfun(x) = mean(x[Altim.madnorm(x).<10])

                                #df = binstats(altim[valid, :], :datetime, 2000:(1/12):2025, :dh; col_function=[binningfun]);
                            
                            
                                df = binstats(altim[valid, :], [:datetime, :height_ref], 
                                    [date_range], :dh; col_function=[binningfun], missing_bins=false)

                                x = BinStatistics.bincenter.(df.datetime)
                                
                                p = lines(x,df.dh_binningfun; title = geotile.id, linewidth=2)
                                display(p)
        
                                gdf = DataFrames.groupby(df, :datetime)

                                #continue 
                            end

                            if minimum(geotile.extent.Y) < -54 || maximum(geotile.extent.Y) > 66
                                canopy = DataFrame(canopyh = zeros(size(altim.latitude)))
                            else
                                path2canopy = joinpath(paths[Symbol(mission)].geotile, geotile.id * ".canopyh")
                                canopy = DataFrame(Arrow.Table(path2canopy));
                            end;

                            if !any(masks0[:, surface_mask])
                                continue
                            end

                            ## Check slope relation
                            if false
                                for dem_id1 in dem_id0           
                                    path2dem = joinpath(paths[Symbol(mission)].geotile, geotile.id * ".$(dem_id1)")
                                    if isfile(path2dem)
                                        dem = select!(DataFrame(Arrow.Table(path2dem)), :height, :dhdx, :dhdy, :dhddx, :dhddy)
                        
                                        dhdx_dhdy = Altim.slope.(dem.dhdx, dem.dhdy, Ref(Altim.dem_info(dem_id1)[1].epsg), lat = mean(geotile.extent.Y));
                                        
                                        var_ind = (.!masks0[:, surface_mask]) .& valid .& .!isnan.(altim.dh) .& (altim.height .!== 0) .& (altim.height_ref .!== 0) .& .!isnan.(altim.curvature) .& (abs.(altim.dh) .< 10) .& masks0.land .& ((canopy.canopyh .<= max_canopy_height) .| (altim.latitude .< -60))
                                        
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
                            showplots && plot(Altim.decimalyear.(altim.datetime[valid]),altim.dh[valid]; seriestype=:scatter, label="all")

                            if curvature_correct
                                var_ind = (.!masks0[:, :landice]) .& valid .& .!isnan.(altim.dh) .& (altim.height .!== 0) .& (altim.height_ref .!== 0) .& .!isnan.(altim.curvature) .& (abs.(altim.dh) .< 10) .& masks0.land .& ((canopy.canopyh .<= max_canopy_height) .| (altim.latitude .< -60))

                                if sum(var_ind) <= length(p1)
                                    warnings && (@warn ("$(geotile.id): no curvature corecction applied, not enought off-ice points"))
                                else
                                    # check bounds: binstats will throw an error if no data is passed to median()
                                    if Altim.vector_overlap(altim[var_ind, :curvature], curvature_range)
                                        df = DataFrame();
                                        try
                                            df = binstats(altim[var_ind, :], [:curvature], [curvature_range], :dh; col_function=[median], missing_bins=true);
                                        catch e
                                            # having issues with edge cases
                                            println(" ##################### ERRORED ##########################")
                                            println(" ##################### $(mission): $(geotile.id) ########################")
                                            println(" ##################### ERRORED ##########################")
                                            rethrow(e)
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
                                                    

                                                    fname = joinpath(fig_folder, "$(surface_mask)_$(mission)_$(geotile.id)_$(dem_id)_curvature.png")
                                                    save(fname, f)
                                                    display(f)
                                                end
                                            end
                                        end
                                    end
                                end
                            end

                            #################################### FILTER 1 ######################################
                            var_ind = masks0[:, surface_mask] .& valid .& .!isnan.(altim.dh) .& (altim.height .!== 0) .& (altim.height_ref .!== 0)
                            if sum(var_ind) < 100
                                continue
                            end

                            # for troubleshooting
                            showplots && plot!(altim.datetime[var_ind],altim.dh[var_ind]; seriestype=:scatter,label="glacier")

                            var_ind[var_ind] = (abs.(altim.dh[var_ind] .- median(altim.dh[var_ind])) .< filt.dh_max) .& (abs.(altim.dh[var_ind]) .< (filt.dh_max*2))

                            if surface_mask == :land
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
                            if !any(date_ind)
                                continue
                            end

                            if !Altim.vector_overlap(altim[var_ind, :datetime], date_range[date_ind]) || 
                                !Altim.vector_overlap(altim[var_ind, :height_ref], height_range[height_ind]) 
                                continue
                            end

                            df = DataFrame()
                            try
                                df = binstats(altim[var_ind, :], [:datetime, :height_ref], 
                                    [date_range[date_ind], height_range[height_ind]], 
                                    :dh; col_function=[binningfun], missing_bins=true)
                            catch e
                                # having issues with edge cases
                                println(" ##################### ERRORED ##########################")
                                println(" ##################### $(geotile.id) ########################")
                                println(" ##################### ERRORED ##########################")
                                rethrow(e)
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
                        
                                isval = .!ismissing.(df[p, "dh_function"])
                                obs2 = @view obs1[i,:]
                                nobs2 = @view nobs1[i, :]
                                if any(isval)
                                    obs2[isval] = df[p, "dh_function"][isval]
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
end