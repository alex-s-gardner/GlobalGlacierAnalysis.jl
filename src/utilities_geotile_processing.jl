## DO NOT CHANGE THESE
begin
binned_filling_parameters = Dict()
binned_filling_parameters[1] = (
    bincount_min=Dict("icesat" => 9,
        "icesat2" => 9,
        "gedi" => 9,
        "hugonnet" => 51), smooth_n=Dict("icesat" => 5,
        "icesat2" => 5,
        "gedi" => 5,
        "hugonnet" => 21), 
    smooth_h2t_length_scale=00, # 800 m = 1 year in distance for anomaly from variogram analysis =    <----------------------
    model1_madnorm_max= 25, # = 5  this is a sigma-equivelent threshold
)
binned_filling_parameters[2] = (
    bincount_min=Dict("icesat" => 9,
        "icesat2" => 9,
        "gedi" => 9,
        "hugonnet" => 51,
    ), smooth_n=Dict("icesat" => 5,
        "icesat2" => 5,
        "gedi" => 5,
        "hugonnet" => 21,
    ), smooth_h2t_length_scale=400, # 800 m = 1 year in distance for anomaly from variogram analysis =
    model1_madnorm_max=10, # this is a sigma-equivelent threshold
)
end


function geotile_binning(; 
    project_id = :v01,
    geotile_width = 2,
    warnings = false,
    showplots = false,

    # run parameters
    force_remake = false,
    update_geotile = false, # this will load in prevous results to update select geotiles or missions
    update_geotile_missions = ["icesat2"],

    # run parameters
    all_permutations_for_glacier_only = true,
    surface_masks = [:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km],
    binned_folders=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"),
    dem_ids = [:best, :cop30_v2],
    binning_methods = ["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"],
    curvature_corrects = [true, false],

     #### DON NOT CHANGE THESE PARAMETERS
    max_canopy_height = 1, # do not change

    # filter parameters
    dh_max=400,
    )

    paths = project_paths(; project_id)
    products = project_products(; project_id)

    # define model for curvature correction
    model1::Function = model1(c, p) = p[1] .+ p[2] .* c .+ p[3] .* c .^ 2 .+ p[4] .* c .^ 3 .+ p[5] .* c .^ 4
    p1 = zeros(5)

    # load geotile definitions with corresponding hypsometry
    geotiles = Altim.geotiles_w_mask(geotile_width)
    #gt_file = joinpath(binned_folder, "geotile_$(:glacier)_hyps.arrow")
    #geotiles = DataFrame(Arrow.Table(gt_file))
    #geotiles.extent = Extent.(getindex.(geotiles.extent, 1))

    # filter geotiles
    geotiles = geotiles[(geotiles.glacier_frac.>0.0), :]
    #geotiles = geotiles[(geotiles.rgi11.>0.0), :]

    # define date bin ranges... this should match what's in gemb_classes_binning.jl
    dd =  30;
    date_range = DateTime(1990):Day(dd):DateTime(2026, 1, 1);
    date_center = date_range[1:end-1] .+ Day(dd/2);

    # curvature ranges 
    Δc = 0.1;
    curvature_range = -1.0:Δc:1.0;
    curvature_center = curvature_range[1:end-1] .+ Δc/2;

    # define elevation
    # DO NOT MODIFY: these need to match elevations in hypsometry
    Δh = 100;
    height_range = 0:100:10000;
    height_center = height_range[1:end-1] .+ Δh/2;

    ## THIS WILL REMOVE ALL _dh_ files if you need to rebuild binned archive
    #files2delete = allfiles("/mnt/bylot-r3/data/binned/2deg"; fn_contains="_dh_")
    #rm.(files2delete)
    ##########

    for binned_folder in binned_folders

        fig_folder = joinpath(binned_folder, "figures")
        !isdir(fig_folder) && mkdir(fig_folder)

        for surface_mask in surface_masks
        #surface_mask = first(surface_masks)

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
            # binning_method = first(binning_methods)
                
                for dem_id in dem_ids  
                # dem_id = first(dem_ids)
                    
                    for curvature_correct in curvature_corrects
                    # curvature_correct = first(curvature_corrects)
                    
                        # skip permutations if all_permutations_for_glacier_only = true
                        if all_permutations_for_glacier_only
                            if all_permutations_for_glacier_only &&
                               ((!(surface_mask == :glacier) && contains(binned_folder, "unfiltered")) &&
                                ((dem_id != :best) && (binning_method != "meanmadnorm3") && curvature_correct))

                                continue
                            end
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

                        println("binning:$out_file")
                        
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
                                dh_hyps = FileIO.load(out_file, "dh_hyps")
                                nobs_hyps = FileIO.load(out_file, "nobs_hyps")
                                curvature_hyps = FileIO.load(out_file, "curvature_hyps")
                                missions = update_geotile_missions;
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

                                missions = keys(dh_hyps)
                            end

                            # 2.6 hours for all 4 missions for all glacierized geotiles
                            # 31 min for icesat + iceast2 + gedi for all glacierized geotiles
                            for mission in missions
                                # <><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><><><><>
                                # mission = "icesat2"
                                #out_file = "/mnt/devon2-r1/devon0/gardnera/Documents/GitHub/Altim/data/geotiles_glacier_hyps_2deg_dh.arrow";
                                #geotiles  = DataFrame(Arrow.Table(out_file));
                                #geotiles.extent = Extent.(getindex.(geotiles.extent, 1));
                                #geotiles = geotiles[(geotiles.rgi2 .> 0) .& (geotiles.glacier_frac .> 0), :] ;
                                # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
                                
                                product = products[Symbol(mission)];

                                #idx = (geotiles.rgi9 .> 0) .& (geotiles.glacier_frac .> 0.1)
                                Threads.@threads for geotile in eachrow(geotiles)
                                    # <><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><><>
                                    #for geotile in eachrow(geotiles[idx,:])
                                    # geotile = geotiles[findfirst((geotiles.rgi2 .> 0.) .& (geotiles.glacier_frac .> 0.1)),:];
                                    #geotile = geotiles[findfirst(geotiles.id .== "lat[+28+30]lon[+082+084]"), :]

                                    #geotile = geotiles[findfirst(geotiles.id .== "lat[-72-70]lon[+160+162]"), :]
                                    #geotile = geotiles[findfirst(geotiles.id .== "lat[+44+46]lon[+006+008]"), :]; # used for curvature figure with mission ="icesat"

                                    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

                                    if geotile.glacier_frac == 0.0
                                        continue
                                    end
                    
                                    t1 = time();

                                    mission_geotile_folder = paths[Symbol(mission)].geotile 

                                    # special case for unfiltered hugonnet data

                                    # this is a hack for including 2 folders that contain hugonnet data (filtered and unfiltered)
                                    if contains(binned_folder, "unfiltered") && (mission == "hugonnet")
                                        mission_geotile_folder = replace(mission_geotile_folder, "/2deg" => "/2deg_unfiltered")
                                    end

                                    path2altim = joinpath(mission_geotile_folder, geotile.id * ".arrow");

                                    if !isfile(path2altim)
                                        continue
                                    end

                                    altim = select!(DataFrame(Arrow.Table(path2altim)), [:longitude, :latitude, :datetime, :height, :quality]);

                                    # add height_ref and curvature
                                    altim = Altim.add_dem_ref!(altim, dem_id, geotile, mission_geotile_folder)
                                    
                                    # load masks
                                    path2masks = joinpath(mission_geotile_folder, geotile.id * ".masks");
                                    masks0 = select!(DataFrame(Arrow.Table(path2masks)), [:inlandwater, :land, :landice, :ocean]);

                                    # add high resolution mask 
                                    masks0 = Altim.highres_mask!(masks0, altim, geotile, feature, invert, excludefeature, surface_mask)
                                    valid = Altim.within.(Ref(geotile.extent), altim.longitude, altim.latitude)

                                    # quick trouble shoot check by plotting time series for each geotile
                                    if false
                                        valid = masks0.landice .& .!isnan.(altim.dh) 
                                        valid[valid] = (abs.(altim.dh[valid] .- median(altim.dh[valid])) .< dh_max) .& (abs.(altim.dh[valid]) .< dh_max*2)
                                        #valid[valid] = (abs.(altim.dh[valid] .- median(altim.dh[valid])) .< 2000) .& (altim.dh[valid] .!== 0) .& (altim.height_ref[valid] .> 0)
                                        if sum(valid) < 100
                                            #continue
                                        end
                                    
                                        df = binstats(altim[valid, :], [:datetime, :height_ref], 
                                            [date_range], :dh; col_function=[binningfun], missing_bins=false)

                                        x = df.datetime
                                        
                                        p = lines(x,df.dh_binningfun; title = geotile.id, linewidth=2)
                                        display(p)
                
                                        gdf = DataFrames.groupby(df, :datetime)

                                        #continue 
                                    end

                                    if minimum(geotile.extent.Y) < -54 || maximum(geotile.extent.Y) > 66
                                        canopy = DataFrame(canopyh = zeros(size(altim.latitude))) 
                                    else
                                        path2canopy = joinpath(mission_geotile_folder, geotile.id * ".canopyh")
                                        canopy = DataFrame(Arrow.Table(path2canopy));
                                    end;

                                    if !any(masks0[:, surface_mask])
                                        continue
                                    end

                                    ## Check slope relation
                                    if false
                                        for dem_id1 in dem_id0           
                                            path2dem = joinpath(mission_geotile_folder, geotile.id * ".$(dem_id1)")
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
                                    # plot(Altim.decimalyear.(altim.datetime[valid]),altim.dh[valid]; label="all")

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
                                                    bin_center = df.curvature;
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
                                                            validX = .!ismissing.(df.dh_median)
                                                            dh_medianX = collect(skipmissing(df.dh_median[validX]))
                                                            p = CairoMakie.plot(bin_center[validX], dh_medianX);

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

                                                            
                                                            
                                                            CairoMakie.plot!(axmain, bin_center[validX], dh_medianX; label = "observation");
                                                            lines!(axmain, bin_center[validX], dh_cor[validX]; label = "model");
                                                            CairoMakie.plot!(axmain, bin_center[validX], dh_medianX .- dh_cor[validX]; label = "corrected");
                                                            barplot!(axbottom, bin_center[validX], df.nrow[validX] / 1000)
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
                                    showplots && CairoMakie.plot!(altim.datetime[var_ind],altim.dh[var_ind]; seriestype=:scatter,label="glacier")

                                    var_ind[var_ind] = (abs.(altim.dh[var_ind] .- median(altim.dh[var_ind])) .< dh_max) .& (abs.(altim.dh[var_ind]) .< (dh_max*2))

                                    if surface_mask == :land
                                        var_ind = var_ind .& (canopy.canopyh .<= max_canopy_height)
                                    end
                                    
                                    # for troubleshooting
                                    showplots && CairoMakie.plot!(altim.datetime[var_ind],altim.dh[var_ind]; seriestype=:scatter, label="glacier-filt", ylims = (-300, +300))
                                    
                                    # this is a hack for including 2 folders that contain hugonnet data (filtered and unfiltered)
                                    if product.apply_quality_filter || (contains(binned_folder, "unfiltered") && (mission == "hugonnet"))
                                        var_ind = var_ind .& altim.quality

                                        # for troubleshooting 
                                        showplots && CairoMakie.plot!(altim.datetime[var_ind],altim.dh[var_ind]; seriestype=:scatter, ylims = (-100, +100))
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
                                    #bindstats can only accept numerical indicies for bins
                                    altim[!,:decimalyear] = Altim.decimalyear.(altim.datetime)

                                    try
                                        df = binstats(altim[var_ind, :], [:decimalyear, :height_ref], 
                                            [Altim.decimalyear.(date_range[date_ind]), height_range[height_ind]], 
                                            :dh; col_function=[binningfun], missing_bins=true)
                                    catch e
                                        # having issues with edge cases
                                        println(" ##################### ERRORED ##########################")
                                        println(" ##################### $(geotile.id) ########################")
                                        println(" ##################### ERRORED ##########################")
                                        rethrow(e)
                                        continue
                                    end

                                    gdf = DataFrames.groupby(df, :decimalyear)
                                    for g in gdf
                                        g[!,:datetime] = Altim.decimalyear2datetime.(g.decimalyear)
                                    end

                                    # create an array of dh as a function of time and elevation
                                    obs1 = fill(NaN, sum(date_ind)-1, sum(height_ind)-1)
                                    nobs1 = fill(Int64(0), sum(date_ind)-1, sum(height_ind)-1)    
                                    p = sortperm(gdf[1].height_ref)

                                    h_center = gdf[1].height_ref[p];
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
                                    showplots && Plots.heatmap(h_center, t_center, obs, clim = (-10, 10))

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
    end
end

function geotile_binned_fill(;
    project_id=:v01,
    geotile_width=2,
    force_remake=false,
    update_geotile=false, # this will load in prevous results to update select geotiles or missions
    update_geotile_missions=["icesat2"],
    plot_dh_as_function_of_time_and_elevation=true,
    mission_reference_for_amplitude_normalization="icesat2",
    all_permutations_for_glacier_only=true,
    surface_masks=[:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km],
    binned_folders=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"),
    dem_ids=[:best, :cop30_v2],
    binning_methods=["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"],
    curvature_corrects=[true, false],
    paramater_sets=[1, 2],
    amplitude_corrects=[true, false],
    showplots,
)

    for binned_folder in binned_folders
        
        fig_folder = joinpath(binned_folder, "figures")
        !isdir(fig_folder) && mkdir(fig_folder)

        for surface_mask in surface_masks
            #surface_mask = first(surface_masks)

            for paramater_set in paramater_sets
                #paramater_set = paramater_sets[1]

                param = binned_filling_parameters[paramater_set]
            
                # load geotiles
                geotiles = Altim.geotiles_mask_hyps(surface_mask, geotile_width)

                # make geotile rgi regions mutually exexclusive 
                geotiles, reg = Altim.geotiles_mutually_exclusive_rgi!(geotiles)

                for binning_method = binning_methods
                    # binning_method = first(binning_methods)

                    for dem_id in dem_ids
                        # dem_id = first(dem_ids)

                        for curvature_correct in curvature_corrects
                            # curvature_correct = first(curvature_corrects)


                            for amplitude_correct = amplitude_corrects
                                # amplitude_correct = first(amplitude_corrects)

                                # skip permutations if all_permutations_for_glacier_only = true
                                if all_permutations_for_glacier_only
                                    if all_permutations_for_glacier_only &&
                                       ((!(surface_mask == :glacier) && contains(binned_folder, "unfiltered")) &&
                                        ((dem_id != :best) && (binning_method != "meanmadnorm3") && curvature_correct))

                                        continue
                                    end
                                end

                                # paths to files
                                binned_file = Altim.binned_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct)
                                binned_filled_file, figure_suffix = Altim.binned_filled_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)

                                if (force_remake || !isfile(binned_filled_file)) && isfile(binned_file)

                                    println("filling binned data: surface_mask = $surface_mask; binning_method = $binning_method; dem_id = $dem_id; curvature_correct = $curvature_correct; amplitude_correct = $amplitude_correct")

                                    dh1 = FileIO.load(binned_file, "dh_hyps")

                                    if .!any(.!isnan.(dh1["hugonnet"]))
                                        println("NOT DATA: $binned_file")
                                        continue
                                    end

                                    if update_geotile
                                       
                                        old_keys = setdiff(keys(dh1), update_geotile_missions)
                                        for k in old_keys
                                            delete!(dh1, k)
                                        end
                                    end

                                    nobs1 = FileIO.load(binned_file, "nobs_hyps")

                                    # align geotile dataframe with DimArrays
                                    gt = collect(dims(dh1[first(keys(dh1))], :geotile))
                                    gt_ind = [findfirst(geotiles.id .== g0) for g0 in gt]
                                    geotiles = geotiles[gt_ind, :]

                                    # create a data frame to store model parameters
                                    # initialize dimensional arrays
                                    params = Dict()
                                    for mission in keys(dh1)
                                        n = length(dims(dh1[mission], :geotile))
                                        push!(params, mission => DataFrame(geotile=val(dims(dh1[mission], :geotile)), nobs_raw=zeros(n), nbins_raw=zeros(n), nobs_final=zeros(n), nbins_filt1=zeros(n), param_m1=[fill(NaN, size(Altim.p1)) for i in 1:n], h0=fill(NaN, n), t0=fill(NaN, n), dh0=fill(NaN, n), bin_std=fill(NaN, n), bin_anom_std=fill(NaN, n)))
                                    end

                                    if plot_dh_as_function_of_time_and_elevation
                                        dh_time_elevation_idx = findfirst(geotiles.id .== "lat[+28+30]lon[+082+084]")
                                        Altim.plot_height_time(dh1; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="raw", fig_folder, figure_suffix, mask=surface_mask, showplots)
                                    end

                                    # I think that the model fitting can have a large influence on how off ice dh is treated. When lots of off ice area is included (e.g. surface_mask == glacier_b1km of glacier_b10km) there can be a very sharp transition at lower elevations from high rates to low rates of elevation change. 

                                    Altim.hyps_model_fill!(dh1, nobs1, params; bincount_min=param.bincount_min,
                                        model1_madnorm_max=param.model1_madnorm_max, smooth_n=param.smooth_n,
                                        smooth_h2t_length_scale=param.smooth_h2t_length_scale)

                                    variogram_range_ratio = false

                                    # make plot of the height to time variogram range ratio
                                    if variogram_range_ratio
                                        range_ratio = Altim.hyps_model_fill!(dh1, nobs1, params; bincount_min=param.bincount_min,
                                            model1_madnorm_max=param.model1_madnorm_max, smooth_n=param.smooth_n,
                                            smooth_h2t_length_scale=param.smooth_h2t_length_scale, variogram_range_ratio)

                                        fontsize = 18
                                        f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98), size=(700, 700), fontsize=fontsize)

                                        for (i, mission) in enumerate(keys(dh1))
                                            valid = .!isnan.(range_ratio[mission])
                                            title = replace.(mission, "hugonnet" => "aster")

                                            title = "$title  mean = $(round(Int, mean(range_ratio[mission][valid])))"
                                            Axis(f[i, 1], title=title)
                                            CairoMakie.hist!(collect(range_ratio[mission][valid]); title=mission, bins=0:250:3000)
                                        end

                                        fname = joinpath(fig_folder, "$(figure_suffix)_variogram_range_ratio.png")
                                        save(fname, f)
                                        display(f)
                                    end

                                    # filter and fit model to geotiles [4 min for for all glacierized geotiles and 4 missions]
                                    if plot_dh_as_function_of_time_and_elevation
                                        Altim.plot_height_time(dh1; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="filled", fig_folder, figure_suffix, mask=surface_mask)
                                    end

                                    # apply seasonal amplitude normalization
                                    if amplitude_correct

                                        # check if reference mission is present, if not, add
                                        ref_added = false
                                        if update_geotile && !any(keys(dh1) .== mission_reference_for_amplitude_normalization)
                                            ref_added = true
                                            (dh_hyps, model_param) = FileIO.load(binned_filled_file, ("dh_hyps", "model_param"))

                                            for k in setdiff(keys(dh_hyps), [mission_reference_for_amplitude_normalization])
                                                delete!(dh_hyps, k)
                                                delete!(model_param, k)
                                            end

                                            dh1 = merge(dh1, dh_hyps)
                                            params = merge(params, model_param)
                                        end

                                        Altim.hyps_amplitude_normalize!(dh1, params; mission_reference=mission_reference_for_amplitude_normalization)

                                        # remove mission_reference_for_amplitude_normalization if it was added.
                                        if ref_added
                                            delete!(dh1, mission_reference_for_amplitude_normalization)
                                            delete!(params, mission_reference_for_amplitude_normalization)
                                        end
                                    end

                                    if plot_dh_as_function_of_time_and_elevation
                                        Altim.plot_height_time(dh1; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="filled_ampnorm", fig_folder, figure_suffix, mask=surface_mask, showplots)
                                    end

                                    Altim.hyps_fill_empty!(dh1, params, geotiles; mask=surface_mask)

                                    Altim.hyps_fill_updown!(dh1, geotiles; mask=surface_mask)

                                    if plot_dh_as_function_of_time_and_elevation
                                        Altim.plot_height_time(dh1; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="filled_updown", fig_folder, figure_suffix, mask=surface_mask, showplots)
                                    end

                                    if update_geotile
                                        (dh_hyps, nobs_hyps, model_param) = FileIO.load(binned_filled_file, ("dh_hyps", "nobs_hyps", "model_param"))
                                        
                    
                                        for k in update_geotile_missions
                                            delete!(dh_hyps, k)
                                            delete!(nobs_hyps, k)
                                            delete!(model_param, k)
                                        end

                                        dh1 = merge(dh1, dh_hyps) 
                                        nobs1 = merge(nobs1, nobs_hyps)
                                        params = merge(params, model_param)
                                    end

                                    # save filled geotiles
                                    save(binned_filled_file, Dict("dh_hyps" => dh1, "nobs_hyps" => nobs1, "model_param" => params))
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

function geotile_filled_landfit(; 
    project_id = "v01",
    showplots=false,
    binned_folders=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"),
    geotile_width=2,
    paramater_sets=[1, 2],
    surface_masks=[:land],
    binning_methods=["meanmadnorm3"],
    dem_ids=[:best],
    curvature_corrects=[true],
    amplitude_corrects=[true],
    force_remake_landoffset=true
    )

    df_param = DataFrame()
    for binned_folder in binned_folders
        for paramater_set in paramater_sets
            for surface_mask in surface_masks
                for dem_id in dem_ids
                    for binning_method in binning_methods
                        for curvature_correct in curvature_corrects
                            for amplitude_correct = amplitude_corrects
                                append!(df_param, DataFrame(; binned_folder, paramater_set, surface_mask, dem_id, binning_method, curvature_correct, amplitude_correct))
                            end
                        end
                    end
                end
            end
        end
    end

    Threads.@threads for r in eachrow(df_param)
    #r = first(eachrow(df_param))

        binned_folder = r.binned_folder
        paramater_set = r.paramater_set
        dem_id = r.dem_id
        surface_mask = r.surface_mask
        curvature_correct = r.curvature_correct
        amplitude_correct = r.amplitude_correct
        binning_method = r.binning_method

        geotiles = copy(Altim.geotiles_mask_hyps(surface_mask, geotile_width))

        # make geotile rgi regions mutually exexclusive 
        geotiles, reg = Altim.geotiles_mutually_exclusive_rgi!(geotiles)

        fig_folder = joinpath(binned_folder, "figures")

        # paths to files
        binned_filled_file, figure_suffix = Altim.binned_filled_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)

        landfit_param_file = replace(binned_filled_file, ".jld2" => "_modelfit.jld2")

        if !isfile(landfit_param_file) || force_remake_landoffset

            println("calulating offset of binned data to $surface_mask: surface_mask = $surface_mask; binning_method = $binning_method; dem_id = $dem_id; curvature_correct = $curvature_correct; amplitude_correct = $amplitude_correct")


            if !(isfile(binned_filled_file))
                printstyled("binned file does not exist, skipping: $(binned_filled_file) \n"; color=:yellow)
                continue
            end

            # load data
            dh1 = copy(FileIO.load(binned_filled_file, "dh_hyps"))
            nobs1 = copy(FileIO.load(binned_filled_file, "nobs_hyps"))

            # align geotiles
            gt = collect(dims(dh1[first(keys(dh1))], :geotile))
            gt_ind = [findfirst(geotiles.id .== g0) for g0 in gt]
            geotiles = geotiles[gt_ind, :]

            # create or extract dimensions
            drgi = Dim{:rgi}(reg)
            ddate = dims(dh1[first(keys(dh1))], :date)
            dgeotile = dims(dh1[first(keys(dh1))], :geotile)
            dmission = Dim{:mission}(collect(keys(dh1)))
            dh_reg = fill(NaN, (dmission, drgi, ddate))
            nobs_reg = fill(NaN, (dmission, drgi, ddate))

            # calculate global agreement with DEM
            for mission in dmission
                # mission = "gedi"

                foo = copy(dh1[mission])
                foo[nobs1[mission].==0] .= NaN
                out = copy(foo[:, :, 1])

                for geotile in dgeotile
                    for date in ddate
                        out[At(geotile), At(date)] = Altim.nanmedian(foo[At(geotile), At(date), :])
                    end
                end

                out2 = deepcopy(foo[1, :, 1])
                for date in ddate
                    out2[At(date)] = Altim.nanmedian(out[:, At(date)])
                end

                out2 = dropdims(out2, dims=(:height, :geotile))
                p = Plots.plot(out2; seriestype=:scatter)
                Plots.xlims!(extrema(ddate[.!isnan.(out2)]))
                Plots.ylabel!("height anomaly [m]")
                Plots.title!(mission)

                valid = .!isnan.(out2)
                y = collect(out2[valid])
                x = Altim.decimalyear.(collect(ddate[valid]))

                fit = curve_fit(Altim.model3, x .- mean(x), y, Altim.p3;)

                fit_date = Base.range(extrema(ddate[.!isnan.(out2)])..., step=Day(10))

                fit_y = Altim.model3(Altim.decimalyear.(fit_date) .- mean(x), fit.param)

                Plots.plot!(fit_date, fit_y; label="trend = $(round(fit.param[2], digits = 3)) m yr⁻¹")

                fname = joinpath(fig_folder, "$(figure_suffix)_$(mission)_$("global")_dh.png")
                save(fname, p)
                showplots && display(p)
            end

            for mission in dmission
                #mission = "hugonnet"

                # NOTE: number of observations over are < 1m canopy height !== total land 
                # area hyspmetric distribution. This is why the plots here do not match those 
                # generated by `dv_reg = Altim.hyps_geotile_aggrigate(dv, geotiles, reg; fun=sum)`

                nopbs_sum = sum(nobs1[mission], dims=:height)
                nopbs_sum = dropdims(nopbs_sum, dims=:height)

                foo = dh1[mission] .* nobs1[mission]
                foo[isnan.(foo)] .= 0
                dh_mean = sum(foo, dims=:height) ./ nopbs_sum
                dh_mean = dropdims(dh_mean, dims=:height)

                for rgi in drgi
                    rgi_ind = geotiles[:, rgi] .> 0
                    geotile_id = geotiles[rgi_ind, :id]

                    foo = dh_mean[At(geotile_id), :]
                    foo[isnan.(foo)] .= 0
                    dh_reg[At(mission), At(rgi), :] = sum(foo .* nopbs_sum[At(geotile_id), :], dims=:geotile) ./ sum(nopbs_sum[At(geotile_id), :], dims=:geotile)

                    nobs_reg[At(mission), At(rgi), :] = sum(nopbs_sum[At(geotile_id), :], dims=:geotile)
                end
            end

            # smooth with loess function
            if false
                dh_reg_smooth = deepcopy(dh_reg)
                x = Altim.decimalyear.(collect(ddate))

                for rgi in drgi
                    for mission in dmission
                        y = dh_reg[At(mission), At(rgi), :]
                        valid = .!isnan.(y)
                        if any(valid)
                            # first do a course fit can remove 3 sigma outliers
                            span = 1.0
                            degree = 2
                            rng, = Altim.validrange(valid)
                            model = MLJ.loess(x[valid], y[valid]; span, degree)
                            dh_reg_smooth[At(mission), At(rgi), rng] = MLJ.predict(model, x[rng])

                            # identify outliers
                            valid[valid] = Altim.madnorm(y[valid] .- dh_reg_smooth[At(mission), At(rgi), valid]) .< 3

                            # fit again with with shorter span
                            a, b = extrema(x[valid])
                            span = min(1.0, (b - a) / 4.0)
                            rng, = Altim.validrange(valid)
                            model = MLJ.loess(x[valid], y[valid]; span, degree)
                            dh_reg_smooth[At(mission), At(rgi), rng] = MLJ.predict(model, x[rng])
                        end
                    end

                    # align hugonnet to icesat2
                    dh_reg_smooth
                end

                begin
                    fontsize = 20
                    rgi = "rgi1"
                    title = "Randolph Glacier Inventory: Region $(rgi[4:end])"
                    ylabel = "height anomaly [m]"

                    fig = CairoMakie.Figure(; backgroundcolor=RGBf(0.98, 0.98, 0.98), size=(1000, 700), fontsize)
                    ax = CairoMakie.Axis(fig[1, 1]; ylabel, title)

                    for mission in dmission
                        CairoMakie.plot!(ax, collect(dh_reg[At(mission), At(rgi), :]); label=mission)
                        CairoMakie.lines!(ax, collect(dh_reg_smooth[At(mission), At(rgi), :]); label=mission)
                    end

                    CairoMakie.ylims!(-2, 2)
                    fig
                end
            end

            dmission = dims(dh_reg, :mission)
            dmetric = Dim{:metric}(["mean", "trend", "acceleration", "amplitude", "date_intercept"])
            fit_param = fill(NaN, (dmission, drgi, dmetric))

            for rgi in drgi
                title = "Randolph Glacier Inventory: Region $(rgi[4:end])"
                f, fit_param[:, At(rgi), :] = Altim.plot_dh(dh_reg[:, At(rgi), :], (nobs_reg[:, At(rgi), :]); title, xlims=(DateTime(2000), DateTime(2024)), ylims=(-3, 3))

                fname = joinpath(fig_folder, "$(figure_suffix)_$(rgi)_dh.png")
                save(fname, f)
                showplots && display(f)
            end

            save(landfit_param_file, Dict("landfit" => fit_param))
        end
    end
end


# bin and calibrate firn data
function geotile_bin_fit_fac(;
    gemb_file_binned="/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0.jld2",
    firn_mission_ref = "icesat2",
    project_id=:v01,
    geotile_width=2,
    surface_masks=[:glacier, :glacier_rgi7, :glacier_b1km, :glacier_b10km],
    binned_folders=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"),
    dem_ids=[:best, :cop30_v2],
    binning_methods=["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"],
    curvature_corrects=[true, false],
    paramater_sets=[1, 2],
    amplitude_corrects=[true, false],
    force_remake_fac = false,
)
    # compute regional volume change, firn correction and mass change
    df_param = DataFrame()
    for binned_folder in binned_folders
        for paramater_set in paramater_sets
            for surface_mask in surface_masks
                for dem_id in dem_ids
                    for binning_method in binning_methods
                        for curvature_correct in curvature_corrects
                            for amplitude_correct = amplitude_corrects
                                append!(df_param, DataFrame(; binned_folder, paramater_set, surface_mask, dem_id, binning_method, curvature_correct, amplitude_correct))
                            end
                        end
                    end
                end
            end
        end
    end

    gemb = copy(FileIO.load(gemb_file_binned, "gemb"))
    geotile_buffer = 1

    Threads.@threads for r in eachrow(df_param)
        #r = first(eachrow(df_param))

        binned_folder = r.binned_folder
        paramater_set = r.paramater_set
        surface_mask = r.surface_mask
        dem_id = r.dem_id
        curvature_correct = r.curvature_correct
        amplitude_correct = r.amplitude_correct
        binning_method = r.binning_method

        # paths to files
        binned_filled_file, figure_suffix = Altim.binned_filled_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)
        gemb_fit_file = splitpath(binned_filled_file)
        gemb_fit_file = joinpath(gemb_fit_file[1:end-1]..., "gemb_" * gemb_fit_file[end])

        if (!isfile(gemb_fit_file) || force_remake_fac)
            geotiles = copy(Altim.geotiles_mask_hyps(surface_mask, geotile_width))

            # make geotile rgi regions mutually exexclusive 
            geotiles, reg = Altim.geotiles_mutually_exclusive_rgi!(geotiles)

            # --------------------------------------------------------------------------------------
            if !(isfile(binned_filled_file))
                printstyled("binned file does not exist, skipping: $(binned_filled_file) \n"; color=:yellow)
                continue
            end

            t1 = time()

            dh = copy(FileIO.load(binned_filled_file, "dh_hyps"))

            # find the best matched GEMB data within geotile
            fac, smb = Altim.facfit(; gemb=deepcopy(gemb), dh, geotiles, geotile_buffer, mission_ref_fac=firn_mission_ref)

            save(gemb_fit_file, Dict("fac" => fac, "smb" => smb,))

            t1 = round((time() - t1) / 60, digits=2)
            println("firn model calibrated [$(t1) min]: $(splitpath(gemb_fit_file)[end])")
        end
    end
end



function geotile_regional_dvdm(;
    firn_mission_ref="icesat2",
    fac_scale_apply=true,
    project_id=:v01,
    geotile_width=2,
    surface_masks=[:glacier, :glacier_rgi7, :glacier_b1km, :glacier_b10km],
    binned_folders=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"),
    dem_ids=[:best, :cop30_v2],
    binning_methods=["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"],
    curvature_corrects=[true, false],
    paramater_sets=[1, 2],
    amplitude_corrects=[true, false],
    force_remake_masschange=false,
    )

     # compute regional volume change, firn correction and mass change
    df_param = DataFrame()
    for binned_folder in binned_folders
        for paramater_set in paramater_sets
            for surface_mask in surface_masks
                for dem_id in dem_ids
                    for binning_method in binning_methods
                        for curvature_correct in curvature_corrects
                            for amplitude_correct = amplitude_corrects
                                append!(df_param, DataFrame(; binned_folder, paramater_set, surface_mask, dem_id, binning_method, curvature_correct, amplitude_correct))
                            end
                        end
                    end
                end
            end
        end
    end

    # this is a bit a hack since landfit is overwritten in Altim.regional_dmass()
    landfit_binned_filled_file, figure_suffix = Altim.binned_filled_filepath(; binned_folder=analysis_paths(; geotile_width).binned, surface_mask=:land, dem_id=:best, binning_method="meanmadnorm3", project_id, curvature_correct=true, amplitude_correct=true, paramater_set=1)
    landfit_param_file = replace(landfit_binned_filled_file, ".jld2" => "_modelfit.jld2")
    landfit = FileIO.load(landfit_param_file, "landfit")


    Threads.@threads for r in eachrow(df_param)
    #for (i,r) in enumerate(eachrow(df_param))
    #r = first(eachrow(df_param))

        binned_folder = r.binned_folder
        paramater_set = r.paramater_set
        surface_mask = r.surface_mask
        dem_id = r.dem_id
        curvature_correct = r.curvature_correct
        amplitude_correct = r.amplitude_correct
        binning_method = r.binning_method

        # paths to files
        binned_filled_file, figure_suffix = Altim.binned_filled_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)
        dvdm_reg_file = replace(binned_filled_file, ".jld2" => "_reg.jld2")

        if !isfile(dvdm_reg_file) || force_remake_masschange
            geotiles = copy(Altim.geotiles_mask_hyps(surface_mask, geotile_width))

            # make geotile rgi regions mutually exexclusive 
            geotiles, reg = Altim.geotiles_mutually_exclusive_rgi!(geotiles)
            t1 = time()
            
            # compute regional volume and mass change, scale firn estimate
            dvdm = Altim.regional_dmass(binned_filled_file, geotiles, reg; surface_mask, fac_scale_apply, firn_mission_ref, landfit)

            # save output
            if !isnothing(dvdm)
                save(dvdm_reg_file, dvdm)
                println("regional volume & mass complete [$(round(Int, time()-t1))s]: surface_mask = $surface_mask, binning_method = $binning_method, dem_id = $dem_id, curvature_correct = $curvature_correct, amplitude_correct = $amplitude_correct")
            end

        end
    end
end


function geotile_dvdm_synthesize(;
    # best estimate ,
    surface_mask_best="glacier",
    dem_best="best",
    curvature_correct_best=true,
    amplitude_correct_best=true,
    binning_method_best="meanmadnorm3",
    fill_param_best=1,
    binned_folder_best="/mnt/bylot-r3/data/binned/2deg",

    # to include in uncertainty
    surface_masks=["glacier", "glacier_rgi7", "glacier_b1km"],
    dems=["best", "cop30_v2"],
    curvature_corrects=[false, true],
    amplitude_corrects=[true],
    binning_methods=["median", "meanmadnorm5", "meanmadnorm3", "meanmadnorm10"],# ["median", "meanmadnorm10", "meanmadnorm5", "meanmadnorm3"]
    fill_params=[1, 2],
    binned_folders=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"),

    # manual adjustments
    regions_to_overwrite_hugonnet_data_with_model_fit=["rgi19"],

    ## combine regions as the per-sensor level
    # this gets a bit complicated a only regions with consistnet mission inclusion can be grouped for proper 
    region_col="rgi",
    region_combines=(
        ("rgi13", "rgi14", "rgi15") => "hma",
        ("rgi1", "rgi3", "rgi4", "rgi5", "rgi6", "rgi7", "rgi8", "rgi9", "rgi19") => "hll",
        ("rgi2", "rgi10", "rgi11", "rgi12", "rgi13", "rgi14", "rgi15", "rgi17", "rgi18") => "ghll",
        ("rgi1", "rgi3", "rgi4", "rgi6", "rgi7", "rgi8", "rgi9") => "hll_ep",
    ), combine_vars=["area_km2", "dm_gt", "dv_km3", "nobs", "fac_km3", "smb_km3"],
)

    df = Altim.regionfiles2dataframe()

    missions = unique(df.mission)
    dates = DataFrames.metadata(df, "date")
    decyear = Altim.decimalyear.(dates)
    dems = unique(df.dem)

    df = Altim.region_combine!(df; region_col, region_combines, combine_vars)

    rgis = unique(df.rgi)
    n = nrow(df)

    # compute data density
    df[!, "nobs_km2"] .= 0.0
    for r in eachrow(df)
        hasobs = r.nobs .> 0
        r.nobs_km2 = mean(r.nobs[hasobs] ./ r.area_km2)
    end
    df[isnan.(df.nobs_km2), "nobs_km2"] .= 0.0

    # compute binary indicies
    surface_valid = falses(n)
    for surface_mask in surface_masks
        surface_valid = surface_valid .| (df.surface_mask .== surface_mask)
    end
    surface_valid_best = (df.surface_mask .== surface_mask_best)

    binning_valid = falses(n)
    for binning_method in binning_methods
        binning_valid = binning_valid .| (df.binning_method .== binning_method)
    end
    binning_valid_best = (df.binning_method .== binning_method_best)

    dem_valid = falses(n)
    for dem in dems
        dem_valid = dem_valid .| (df.dem .== dem)
    end
    dem_valid_best = (df.dem .== dem_best)

    curvature_valid = falses(n)
    for curvature_correct in curvature_corrects
        curvature_valid = curvature_valid .| (df.curvature_correct .== curvature_correct)
    end
    curvature_valid_best = (df.curvature_correct .== curvature_correct_best)

    amplitude_valid = falses(n)
    for amplitude_correct in amplitude_corrects
        amplitude_valid = amplitude_valid .| (df.amplitude_correct .== amplitude_correct)
    end
    amplitude_valid_best = (df.amplitude_correct .== amplitude_correct_best)

    fill_valid = falses(n)
    for fill_param in fill_params
        fill_valid = fill_valid .| (df.fill_param .== fill_param)
    end
    fill_valid_best = (df.fill_param .== fill_param_best)

    folder_valid = falses(n)
    for binned_folder in binned_folders
        folder_valid = folder_valid .| (df.binned_folder .== binned_folder)
    end
    folder_valid_best = (df.binned_folder .== binned_folder_best)

    param_valid = surface_valid .& dem_valid .& curvature_valid .& amplitude_valid .& fill_valid .& folder_valid .& binning_valid
    param_valid_best = surface_valid_best .& dem_valid_best .& curvature_valid_best .& amplitude_valid_best .& fill_valid_best .& folder_valid_best .& binning_valid_best

    # create an rgi glacier area data frame
    x = df[(df.mission.=="icesat").&(df.surface_mask.=="glacier"), :]
    idx = unique(z -> x.area_km2[z], 1:length(x.area_km2))

    df_area = DataFrame(rgi=x[idx, :rgi], area_km2=x[idx, :area_km2])

    # find index where each mission has valid data... choose rgi2 since all missions are valid there
    mission_ts_valid = Dict()
    for mission in missions
        valid = .!isnan.(df.dm_gt[findfirst((df.mission .== mission) .& param_valid .& (df.rgi .== "rgi2"))])
        mission_ts_valid[mission] = valid
    end

    # exclude data that does not meet nobs per area criteria [check across rgi and mission]
    nobs_km2_threshold = 0.01

    for rgi in rgis
        #rgi = "rgi1"
        rgi_index = (df.rgi .== rgi)

        for mission in missions
            index = rgi_index .& (df.mission .== mission)
            if mean(df.nobs_km2[index]) < nobs_km2_threshold
                for r in eachrow(@view df[index, :])
                    for var0 in ["dm_gt", "dv_km3"]
                        r[var0] .= NaN
                    end
                    r["nobs"] .= 0.0
                end
            end
        end
    end

    # print excluded regions
    for mission in missions
        b = falses(size(rgis))
        for (i, rgi) in enumerate(rgis)
            rgi_index = (df.rgi .== rgi)

            if mean(df.nobs_km2[(df.mission.==mission).&rgi_index]) >= nobs_km2_threshold
                b[i] = true
            end
        end
        println("$(mission) excluded regions: $(rgis[.!b]))")
    end

    # [1] Set the mean over the icesat2 period to zero by removing the mean offset to the 
    # reference surface and remove that same offset from icesat.
    # [2] Align GEDI and ASTER to the median icesat + icesat2 offset
    mission_index = Dict()
    mission_df = Dict()
    for mission in missions
        mission_index[mission] = df.mission .== mission
        mission_df[mission] = @view df[mission_index[mission], :]
    end

    n = nrow(mission_df["icesat2"])
    for ricesat2 in eachrow(mission_df["icesat2"])
        #ricesat2 = first(df[index_icesat2, :])

        # find matching icesat row
        match_index = Dict()
        for mission in setdiff(missions, ["icesat2"])
            match_index[mission] = trues(n)
            for col in ["rgi", "surface_mask", "dem", "curvature_correct", "binning_method", "project_id", "amplitude_correct", "fill_param", "binned_folder"]
                match_index[mission] = match_index[mission] .& (ricesat2[col] .== mission_df["icesat"][:, col])
            end

            if sum(match_index[mission]) !== 1
                error("non unique match: $mission")
            end

            match_index[mission] = findfirst(match_index[mission])
        end

        ricesat = @view mission_df["icesat"][match_index["icesat"], :]

        # [1] Set the mean over the icesat2 period to zero by removing the mean offset to the 
        # reference surface and remove that same offset from icesat.
        for var0 in ["dm_gt", "dv_km3"]
            offset = Altim.nanmean(ricesat2[var0])

            # remove offset from ICESat2
            ricesat2[var0] = ricesat2[var0] .- offset

            #remove same offset from ICESat
            ricesat[var0] = ricesat[var0] .- offset
        end

        # [2] Align GEDI and ASTER to the median icesat + icesat2 offset
        rhugonnet = @view mission_df["hugonnet"][match_index["hugonnet"], :]
        rgedi = @view mission_df["gedi"][match_index["gedi"], :]

        for var0 in ["dm_gt", "dv_km3"]

            offset = Altim.nanmedian(rgedi[var0] .- ricesat2[var0])
            rgedi[var0] = rgedi[var0] .- offset

            offset = Altim.nanmedian(hcat(rhugonnet[var0] .- ricesat2[var0], rhugonnet[var0] .- ricesat[var0]))
            rhugonnet[var0] = rhugonnet[var0] .- offset
        end
    end


    # replace poor observational records with model fit to better observations
    fit_missions = ["icesat", "icesat2", "gedi"]
    n = nrow(mission_df["hugonnet"])

    if .!isnothing(regions_to_overwrite_hugonnet_data_with_model_fit)

        for rgi in regions_to_overwrite_hugonnet_data_with_model_fit
            #rgi = first(regions_to_overwrite_hugonnet_data_with_model_fit)

            for (k, r) in enumerate(eachrow(mission_df["hugonnet"]))
                #k = 1; r = first(eachrow(mission_df["hugonnet"]))
                if r.rgi !== rgi
                    continue
                end

                match_index = Dict()
                for mission in fit_missions
                    match_index[mission] = trues(n)
                    for col in ["rgi", "surface_mask", "dem", "curvature_correct", "binning_method", "project_id", "amplitude_correct", "fill_param", "binned_folder"]
                        match_index[mission] = match_index[mission] .& (mission_df[mission][:, col] .== r[col])
                    end
                end

                for var0 = ["dm_gt", "dv_km3"]
                    #var0 = "dm_gt"

                    y = fill(NaN, (length(r[var0]), length(fit_missions)))
                    for (i, mission) in enumerate(fit_missions)
                        y[:, i] = mission_df[mission][match_index[mission], var0][1]
                    end

                    y = [Altim.nanmean(foo) for foo in eachrow(y)]

                    notnan = .!isnan.(y)
                    x = decyear .- mean(decyear)
                    foo_fit = Altim.curve_fit(Altim.model3, x[notnan], y[notnan], Altim.p3)

                    notnan = mission_ts_valid["hugonnet"]
                    y .= NaN

                    y[notnan] = Altim.model3(x[notnan], foo_fit.param)

                    mission_df["hugonnet"][k, var0] = y
                end
            end
        end
    end


    # compute central estimate and lower and upper bounds for each mission
    vars = ["dm_gt", "dv_km3", "fac_km3", "smb_km3"]
    quantiles = [0.95]

    df_mission = DataFrame()
    for rgi in rgis
        #rgi = "rgi1"

        for mission in missions
            #mission = "icesat"

            df1 = DataFrame()

            valid = findall((df.rgi .== rgi) .& (mission_index[mission]) .& param_valid)
            valid_best = (df.rgi .== rgi) .& mission_index[mission] .& param_valid_best

            if sum(valid_best) !== 1
                error("non unique best estimate")
            end
            valid_best = findfirst(valid_best)


            for var0 in vars
                #var0 = "dm_gt"

                # central estimate is user selected best estimate
                df1[!, var0] = [df[valid_best, var0]] # this selects the best estimate as the median

                # all estimates
                f = reduce(hcat, df[valid, var0])
                f[f.==0] .= NaN

                # remove central estimate to get anomaly
                f = f .- repeat(df1[!, var0][1], 1, size(f, 2))

                # take absolute as I want a "fair" cemetric error
                f = abs.(f)

                notnan = .!isnan.(f[:, 1])

                q = reduce(hcat, Altim.nanquantile.(eachrow(f), Ref(quantiles)))
                sigma = q[1, :]

                # make error symetric
                df1[!, var0*"_low"] = [df1[!, var0][1] .- sigma]
                df1[!, var0*"_high"] = [df1[!, var0][1] .+ sigma]
            end

            df1[!, "nobs"] .= [df[valid_best, "nobs"]]
            df1[!, "mission"] .= mission
            df1[!, "rgi"] .= rgi

            append!(df_mission, df1)
        end
    end

    # rearange df_mission from columns to rows 
    vars = ["dm_gt", "dv_km3", "fac_km3", "smb_km3"]
    df1 = DataFrame()

    for var0 in vars
        for rgi in rgis
            for mission in missions
                index = (df_mission.mission .== mission) .& (df_mission.rgi .== rgi)

                #println("var0 = $var0; rgi = $rgi, mission = $mission")
                low = df_mission[index, var0*"_low"]
                mid = df_mission[index, var0]
                high = df_mission[index, var0*"_high"]

                area_km2 = df_area[df_area.rgi.==rgi, :area_km2]
                nobs0 = deepcopy(df_mission[index, "nobs"])

                if (var0 == "fac_km3") || (var0 == "smb_km3")
                    nobs0[1] .= NaN
                end

                append!(df1, DataFrame(; var=var0, rgi, mission, area_km2, mid, low, high, nobs=nobs0))
            end
        end
    end
    df_mission = df1

    # store date as metadata
    metadata!(df_mission, "date", dates)


    ## synthesize observations (error weighted average)
    vars = ["dm_gt", "dv_km3", "fac_km3", "smb_km3"]
    df_synth = DataFrame()

    for var0 in vars
        #var0 = "dm_gt"

        var_index = (df_mission.var .== var0)

        for rgi in rgis
            #rgi = "rgi19"

            rgi_index = (df_mission.rgi .== rgi)

            if (var0 == "dm_gt") || (var0 == "dv_km3")
                foo_sum = zeros(size(decyear))
                foo_weight = zeros(size(decyear))
                foo_err = zeros(size(decyear))
                foo_n = zeros(size(decyear))

                for mission in missions
                    mission_index0 = (df_mission.mission .== mission)

                    index = var_index .& rgi_index .& mission_index0

                    low = replace(df_mission[index, "low"][1], NaN => 0)
                    mid = replace(df_mission[index, "mid"][1], NaN => 0)
                    high = replace(df_mission[index, "high"][1], NaN => 0)

                    err = (high .- low) ./ 4
                    foo_err[err.>0] .+= err[err.>0] .^ 2
                    foo_n[err.>0] .+= 1

                    w = replace(1 ./ (err) .^ 2, Inf => 0)
                    foo_sum .+= (mid .* w)
                    foo_weight .+= w
                end

                mid = foo_sum ./ foo_weight
                error0 = sqrt.(foo_err) ./ foo_n .* 2
                notnan = .!isnan.(mid)

                low = mid .- error0
                high = mid .+ error0
            else
                mission = "icesat2"
                mission_index0 = (df_mission.mission .== mission)

                index = var_index .& rgi_index .& mission_index0

                low = replace(df_mission[index, "low"][1], NaN => 0)
                mid = replace(df_mission[index, "mid"][1], NaN => 0)
                high = replace(df_mission[index, "high"][1], NaN => 0)
            end

            index = var_index .& rgi_index
            nobs0 = vec(sum(hcat(df_mission[index, "nobs"]...), dims=2))

            area_km2 = df_area[df_area.rgi.==rgi, :area_km2]

            low = [low]
            high = [high]
            mid = [mid]

            append!(df_synth, DataFrame(; var=var0, rgi, mission="synthesis", area_km2, mid, low, high, nobs=[nobs0]))
        end
    end

    # combine mission and synthesis results
    df_dmass = vcat(df_mission, df_synth)

    # add units as column and remove from var name
    df_dmass[!, "unit"] .= ""
    for (i, r) in enumerate(eachrow(df_dmass))
        if split(r.var, "_")[2] == "gt"
            df_dmass[i, "unit"] = "Gt"
        elseif split(r.var, "_")[2] == "km3"
            df_dmass[i, "unit"] = "km⁻³"
        end

        df_dmass[i, "var"] = split(r.var, "_")[1]
    end

    metadata!(df_dmass, "date", dates)
    return df_dmass
end

# add grace results to dataframe
function geotile_dvdm_addgrace!(df_dmass)

    grace = Altim.read_grace_rgi(; datadir=setpaths()[:grace_rgi])
    rgis = reduce(vcat, (["rgi$i" for i in 1:19], ["HMA"]))
    var0 = "dm"
    mission = "grace"
    unit = "Gt"
    mask_icesheets = true

    dates = DataFrames.metadata(df_dmass, "date")
    decyear = Altim.decimalyear.(dates)

    missions = unique(df_dmass.mission)
    # find index where each mission has valid data... choose rgi2 since all missions are valid there
    mission_ts_valid = Dict()
    for mission in missions
        valid = .!isnan.(df_dmass.mid[findfirst((df_dmass.mission .== mission) .& (df_dmass.var .== "dm") .& (df_dmass.rgi .== "rgi2"))])
        mission_ts_valid[mission] = valid
    end

    x = df_dmass[(df_dmass.mission.=="icesat"), :]
    idx = unique(z -> x.area_km2[z], 1:length(x.area_km2))

    df_area = DataFrame(rgi=x[idx, :rgi], area_km2=x[idx, :area_km2])


    for rgi in rgis
        #rgi = first(rgis)

        if any(keys(grace) .== rgi)
            datenum = grace[rgi]["dM_gt_mdl_fill_date"]
            decyear0 = vec(Altim.decimalyear.(Altim.datenum2date.(datenum)))
            dm_gt = vec(grace[rgi]["dM_gt_mdl_fill"])
            dm_sigma = vec(grace[rgi]["dM_sigma_gt_mdl_fill"])
            notnan = .!isnan.(dm_gt)

            # introplate to df time
            valid = (decyear .<= maximum(decyear0[notnan])) .& (decyear .>= minimum(decyear0[notnan]))

            mid = fill(NaN, length(valid))
            low = fill(NaN, length(valid))
            high = fill(NaN, length(valid))

            # quadratic spline results in some high osilations in some regions (e.g. Svalbard)
            #interp = QuadraticSpline(dm_gt[notnan], decyear0[notnan])
            interp = DataInterpolations.LinearInterpolation(dm_gt[notnan], decyear0[notnan])
            mid[valid] = interp(decyear[valid])

            interp = QuadraticSpline(dm_sigma[notnan], decyear0[notnan])
            dm_sigma = interp(decyear[valid])

            low[valid] = mid[valid] .- dm_sigma
            high[valid] = mid[valid] .+ dm_sigma

            # ceck that interpoalted results are consistent with orginal data
            if false
                p = lines(decyear0, dm_gt)
                lines!(decyear, mid)
            end

            if rgi == "HMA"
                rgi = "hma"
            end
        else
            if rgi == "HMA"
                rgi = "hma"
            end

            index = (df_dmass.var .== "dm") .& (df_dmass.mission .== "synthesis") .& (df_dmass.rgi .== rgi)
            # if region does not exist, substitute with synthesis data
            mid = df_dmass.mid[index][1]
            low = df_dmass.low[index][1]
            high = df_dmass.high[index][1]
        end

        index = (df_dmass.var .== "dm") .& (df_dmass.mission .== "synthesis") .& (df_dmass.rgi .== rgi)

        # align records: subtract the median offset to synthesis over the icesat and icesat-2 period 
        valid = (mission_ts_valid["icesat"] .| mission_ts_valid["icesat2"]) .& .!isnan.(mid)
        mid_synth = df_dmass.mid[index][1]
        offset = Altim.nanmean(mid[valid] .- mid_synth[valid])

        area_km2 = df_area[df_area.rgi.==rgi, :area_km2]
        nobs0 = ones(length(decyear))

        # mask out greenland and antarctic 
        if mask_icesheets && ((rgi == "rgi5") || (rgi == "rgi19"))
            mid .= NaN
            low .= NaN
            high .= NaN
            nobs0 .= 0.0
        end

        mid = [mid .- offset]
        low = [low .- offset]
        high = [high .- offset]
        nobs0 = [nobs0]

        append!(df_dmass, DataFrame(; var=var0, rgi, mission, area_km2, mid, low, high, nobs=nobs0, unit))
    end
    metadata!(df_dmass, "date", dates)
    return df_dmass
end

## combine synthisized regions [can now combine regions covered by different sensors]
function geotile_combine_synth_regions!(df_dmass)
    region_col = "rgi"
    combine_vars = ["area_km2", "mid", "nobs"]
    rss_vars = ["low", "high"]

    missions = unique(df_dmass.mission)

    dates = DataFrames.metadata(df_dmass, "date")

    # to RSS variables you need to frist subtract :mid from :low and high
    for r in eachrow(df_dmass)
        r.low = r.mid .- r.low
        r.high = r.high .- r.mid
    end

    # combine synthesis data for global mass change
    begin
        missions = ["synthesis"]
        combined_regs = ["rgi16", "hll", "ghll"]
        region_combines = (
            tuple(combined_regs...) => "global",
        )

        df_dmass = Altim.region_combine!(df_dmass; region_col, region_combines, combine_vars, rss_vars, missions)
    end

     # combine synthesis data for global mass change excluding periphery
     begin
        missions = ["synthesis"]
        combined_regs = ["rgi16", "hll_ep", "ghll"]
        region_combines = (
            tuple(combined_regs...) => "global_ep",
        )

        df_dmass = Altim.region_combine!(df_dmass; region_col, region_combines, combine_vars, rss_vars, missions)
    end

    # combine grace data for global mass change excluding periphery
    begin
        missions = ["grace"]
        rgis = reduce(vcat, (["rgi$i" for i in 1:4], ["rgi$i" for i in 6:18]))
        region_combines = (
            tuple(rgis...) => "global_ep",
        )
        df_dmass = Altim.region_combine!(df_dmass; region_col, region_combines, combine_vars, rss_vars, missions)
    end

    # add back :mid to :low and :high
    for r in eachrow(df_dmass)
        r.low = r.mid .- r.low
        r.high = r.mid .+ r.high
    end

    #Altim.nanmean(df_dmass[457, :mid] .- df_dmass[457, :low])
    # global mean error = 180 Gt if RSS all rgi region errors together (assumes not correlation between regions)
    # global mean error = 267 Gt if RSS all 3 mission combined region errors (assumes no correlation between mission combined regions)
    # global mean error = 526 Gt if ADD all rgi region errors (assumes 100% correlation between regions)
    # global mean error = 352 Gt if ADD all 3 mission combined region errors (assumes 100% correlation between mission combined regions)
    metadata!(df_dmass, "date", dates)
    return df_dmass
end


## add estimate of trend, acceleration and uncertainty
function geotile_dvdm_add_trend!(df_dmass; iterations = 1000)

    dates = DataFrames.metadata(df_dmass, "date")
    decyear = Altim.decimalyear.(dates)

    vars = ["trend", "trend_err", "acceleration", "acceleration_err", "amplitude", "amplitude_err"]
    for var0 in vars
        df_dmass[!, var0] .= 0.0
    end

    # this takes 80s 
    Threads.@threads for r in eachrow(df_dmass)

        fit0 = Altim.iterative_model2_fit(r.mid, r.low, r.high, decyear; iterations)

        r.trend = fit0.trend
        r.trend_err = fit0.trend_err
        r.acceleration = fit0.acceleration
        r.acceleration_err = fit0.acceleration_err
        r.amplitude = fit0.amplitude
        r.amplitude_err = fit0.amplitude_err
    end
    metadata!(df_dmass, "date", dates)
end


# convert from Gt/km3 to m.w.e/m
function geotile_dvdm_areaaverage(df_dmass)
    dates = DataFrames.metadata(df_dmass, "date")
    df_dheight = deepcopy(df_dmass)
    for (i, r) in enumerate(eachrow(df_dheight))
        for v in ["mid", "low", "high", "trend", "trend_err", "acceleration", "acceleration_err", "amplitude", "amplitude_err"]
            df_dheight[i, v] = r[v] / r.area_km2 * 1000
        end
        df_dheight[i, "unit"] = Altim.unit2areaavg[r.unit]
    end
    metadata!(df_dmass, "date", dates)
end