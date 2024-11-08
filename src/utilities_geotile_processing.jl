## DO NOT CHANGE THESE
begin
binned_filling_parameters = Dict()
binned_filling_parameters[1] = (
    bincount_min=Dict("icesat" => 9,
        "icesat2" => 9,
        "gedi" => 9,
        "hugonnet" => 51), 
    smooth_n=Dict("icesat" => 5,
        "icesat2" => 5,
        "gedi" => 5,
        "hugonnet" => 21), 
    smooth_h2t_length_scale=800, # 800 m = 1 year in distance for anomaly from variogram analysis =    <----------------------
    model1_madnorm_max = 5, # = 5  this is a sigma-equivelent threshold
    adjust2land=false,
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
    adjust2land=false,
)
    binned_filling_parameters[3] = (
        bincount_min=Dict("icesat" => 9,
            "icesat2" => 9,
            "gedi" => 9,
            "hugonnet" => 51), 
        smooth_n=Dict("icesat" => 5,
            "icesat2" => 5,
            "gedi" => 5,
            "hugonnet" => 21),
        smooth_h2t_length_scale=800, # 800 m = 1 year in distance for anomaly from variogram analysis =    <----------------------
        model1_madnorm_max=5, # = 5  this is a sigma-equivelent threshold
        adjust2land=true,
    )
    binned_filling_parameters[4] = (
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
        adjust2land=true,
    )
end

# define model for curvature correction
model_curvature::Function = model_curvature(c, p) = p[1] .+ p[2] .* c .+ p[3] .* c .^ 2 .+ p[4] .* c .^ 3 .+ p[5] .* c .^ 4
p_curvature = zeros(5)

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
    dh_max=200,
    )

    paths = project_paths(; project_id)
    products = project_products(; project_id)

    # load geotile definitions with corresponding hypsometry
    geotiles = Altim.geotiles_w_mask(geotile_width)
    #gt_file = joinpath(binned_folder, "geotile_$(:glacier)_hyps.arrow")
    #geotiles = DataFrame(Arrow.Table(gt_file))
    #geotiles.extent = Extent.(getindex.(geotiles.extent, 1))

    # filter geotiles
    geotiles = geotiles[(geotiles.glacier_frac.>0.0), :]
    #geotiles = geotiles[(geotiles.rgi11.>0.0), :]

    # define date bin ranges... this should match what's in gemb_classes_binning.jl
    # define date and hight binning ranges 
    date_range, date_center = Altim.project_date_bins()
    Δd = abs(date_center[2] - date_center[1])
    height_range, height_center = Altim.project_height_bins()
    Δh = abs(height_center[2] - height_center[1])

    # curvature ranges 
    Δc = 0.1;
    curvature_range = -1.0:Δc:1.0;
    curvature_center = curvature_range[1:end-1] .+ Δc/2;

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
                    
                        binned_file = Altim.binned_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct)

                        println("binning:$binned_file")
                        
                        # <><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><>
                        # binningfun(x) = mean(x[Altim.madnorm(x).<10])
                        # <><><><><><><><><><><><><><><><><><><><><><<><><><><><><><><><><><><><><><><><>
                        if !isfile(binned_file) || force_remake
                            # 6.2 hours for all glaciers, all missions/datasets on 96 threads
                            # 10hr for land for all glaciers, all missions/datasets on 96 threads

                            # initialize dimensional arrays
                            # update_geotile = true
                            if update_geotile
                                # load exisiting
                                dh_hyps = FileIO.load(binned_file, "dh_hyps")
                                nobs_hyps = FileIO.load(binned_file, "nobs_hyps")
                                curvature_hyps = FileIO.load(binned_file, "curvature_hyps")
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
                                    push!(curvature_hyps, String(mission) => DataFrame(id=geotiles.id, curvature=[curvature_range for i in 1:ngeotile], dh = [fill(NaN,ncurvature) for i in 1:ngeotile], nobs=[zeros(ncurvature) for i in 1:ngeotile], model_coef = [zeros(size(p_curvature)) for i in 1:ngeotile]))
                                end

                                missions = keys(dh_hyps)
                            end

                            # 2.6 hours for all 4 missions for all glacierized geotiles
                            # 31 min for icesat + iceast2 + gedi for all glacierized geotiles
                            for mission in missions
                                # <><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><><><><>
                                # mission = "icesat2"
                                #binned_file = "/mnt/devon2-r1/devon0/gardnera/Documents/GitHub/Altim/data/geotiles_glacier_hyps_2deg_dh.arrow";
                                #geotiles  = DataFrame(Arrow.Table(binned_file));
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
                                                
                                                if sum(var_ind) <= length(p_curvature)
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

                                        if sum(var_ind) <= length(p_curvature)
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

                                                    if (sum(bin_valid) <= length(p_curvature))
                                                        warnings && (@warn ("$(geotile.id): no curvature corecction applied, not enought off-ice points"))
                                                    else
                                                        fit1 = curve_fit(model_curvature, bin_center[bin_valid], df.dh_median[bin_valid], df.nrow[bin_valid], p_curvature);
                                                        #histogram(altim.dh[var_ind], bins= -5:.1:5; label ="raw", xlabel="height change [m]")
                                                        altim.dh = altim.dh .- model_curvature(altim.curvature, fit1.param);
                                                        #histogram!(altim.dh[var_ind], bins= -5:.1:5; label="curvature corrected")

                                                        gt_ind = findfirst(geotile.id .== curvature_hyps[mission].id)
                                                        cdh = @view curvature_hyps[mission].dh[gt_ind] 
                                                        cdh[1][bin_valid] = df.dh_median[bin_valid];
                                                        cnobs = @view curvature_hyps[mission][gt_ind, :nobs]
                                                        cnobs[1][bin_valid] = df.nrow[bin_valid]
                                                        curvature_hyps[mission][gt_ind, :model_coef] = fit1.param
                                            
                                                        if showplots
                                                            dh_cor = model_curvature(bin_center, fit1.param);
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

                                    decyear_range = Altim.decimalyear.(date_range)

                                    altim[!, :decimalyear] = Altim.decimalyear.(altim.datetime)
                                    
                                    var0, nobs0 = Altim.geotile_bin2d(
                                        altim[var_ind, :];
                                        var2bin="dh",
                                        dims_edges=("decimalyear" => decyear_range, "height_ref" => height_range),
                                        binfunction=Altim.binningfun_define(binning_method))

                                    if isnothing(var0)
                                        continue
                                    end

                                    # for troubleshooting 
                                    showplots && Plots.heatmap(var0, clim=(-10, 10))

                                    try
                                        dh_hyps[mission][At(geotile.id), :, :] = var0
                                        nobs_hyps[mission][At(geotile.id), :, :] = nobs0
                                    catch e
                                        println((; binned_folderbinned_folder, surface_mask,binning_method,dem_id,curvature_correct, mission, geotile = geotile.id))
                                        println("size var0 = $(size(var0))")
                                        println("size nobs0 = $(size(nobs0))")
                                        throw(e)
                                    end


                                    t2 = time()
                                    dt = round(Int16, t2 - t1);
                                    println("binned: $(mission) - $(geotile.id): $(dt)s")
                                end
                            end
                            save(binned_file, Dict("dh_hyps" => dh_hyps, "nobs_hyps" => nobs_hyps, "curvature_hyps" => curvature_hyps));
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

    params = NamedTuple[]
    for binned_folder in binned_folders
        for paramater_set in paramater_sets
            for surface_mask in surface_masks
                for dem_id in dem_ids
                    for binning_method in binning_methods
                        for curvature_correct in curvature_corrects
                            for amplitude_correct = amplitude_corrects
                                push!(params, (; project_id, binned_folder, paramater_set, surface_mask, dem_id, binning_method, curvature_correct, amplitude_correct))
                            end
                        end
                    end
                end
            end
        end
    end

    Threads.@threads for param in params
        #param = first(params)

        fig_folder = joinpath(param.binned_folder, "figures")
        !isdir(fig_folder) && mkdir(fig_folder)

        param_filling = Altim.binned_filling_parameters[param.paramater_set]

        # load geotiles
        geotiles = Altim.geotiles_mask_hyps(param.surface_mask, geotile_width)

        # make geotile rgi regions mutually exexclusive 
        geotiles, reg = Altim.geotiles_mutually_exclusive_rgi!(geotiles)

        # skip permutations if all_permutations_for_glacier_only = true
        if all_permutations_for_glacier_only
            if all_permutations_for_glacier_only &&
               ((!(param.surface_mask == :glacier) && contains(param.binned_folder, "unfiltered")) &&
                ((param.dem_id != :best) && (param.binning_method != "meanmadnorm3") && param.curvature_correct))

                continue
            end
        end

        # paths to files
        binned_file = Altim.binned_filepath(; param.binned_folder, param.surface_mask, param.dem_id, param.binning_method, project_id, param.curvature_correct)
        binned_file_land = Altim.binned_filepath(; param.binned_folder, surface_mask=:land, param.dem_id, param.binning_method, project_id, param.curvature_correct)
        binned_filled_file, figure_suffix = Altim.binned_filled_filepath(; param.binned_folder, param.surface_mask, param.dem_id, param.binning_method, project_id, param.curvature_correct, param.amplitude_correct, param.paramater_set)

        if (force_remake || !isfile(binned_filled_file)) && isfile(binned_file)

            println("filling binned data: surface_mask = $(param.surface_mask); binning_method = $(param.binning_method); dem_id = $(param.dem_id); curvature_correct = $(param.curvature_correct); amplitude_correct = $(param.amplitude_correct)")

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
            params_fill = Dict()
            for mission in keys(dh1)
                n = length(dims(dh1[mission], :geotile))
                push!(params_fill, mission => DataFrame(geotile=val(dims(dh1[mission], :geotile)), nobs_raw=zeros(n), nbins_raw=zeros(n), nobs_final=zeros(n), nbins_filt1=zeros(n), param_m1=[fill(NaN, size(Altim.p1)) for i in 1:n], h0=fill(NaN, n), t0=fill(NaN, n), dh0=fill(NaN, n), bin_std=fill(NaN, n), bin_anom_std=fill(NaN, n)))
            end

            if plot_dh_as_function_of_time_and_elevation
                dh_time_elevation_idx = findfirst(geotiles.id .== "lat[+28+30]lon[+082+084]")
                Altim.plot_height_time(dh1; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="raw", fig_folder, figure_suffix, mask=param.surface_mask, showplots)
            end

            # I think that the model fitting can have a large influence on how off ice dh is treated. When lots of off ice area is included (e.g. param.surface_mask == glacier_b1km of glacier_b10km) there can be a very sharp transition at lower elevations from high rates to low rates of elevation change. 
            if param_filling.adjust2land && any(keys(dh1) .== "huggonet")
                dh_land = FileIO.load(binned_file_land, "dh_hyps")
                nobs_land = FileIO.load(binned_file_land, "nobs_hyps")

                dh1["hugonnet"] = Altim.geotile_adjust!(
                    dh1["hugonnet"],
                    dh_land["hugonnet"],
                    nobs_land["hugonnet"],
                    dh_land["icesat2"];
                    ref_madnorm_max=10,
                    minnobs=45,
                    ref_minvalid=12,
                )
            end

            Altim.hyps_model_fill!(dh1, nobs1, params_fill; bincount_min=param_filling.bincount_min,
                model1_madnorm_max=param_filling.model1_madnorm_max, smooth_n=param_filling.smooth_n,
                smooth_h2t_length_scale=param_filling.smooth_h2t_length_scale)

            variogram_range_ratio = false

            # make plot of the height to time variogram range ratio
            if variogram_range_ratio
                range_ratio = Altim.hyps_model_fill!(dh1, nobs1, params_fill; bincount_min=param_filling.bincount_min,
                    model1_madnorm_max=param_filling.model1_madnorm_max, smooth_n=param_filling.smooth_n,
                    smooth_h2t_length_scale=param_filling.smooth_h2t_length_scale, variogram_range_ratio)

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
                Altim.plot_height_time(dh1; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="filled", fig_folder, figure_suffix, mask=param.surface_mask)
            end

            # apply seasonal amplitude normalization
            if param.amplitude_correct

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
                    params_fill = merge(params_fill, model_param)
                end

                Altim.hyps_amplitude_normalize!(dh1, params_fill; mission_reference=mission_reference_for_amplitude_normalization)

                # remove mission_reference_for_amplitude_normalization if it was added.
                if ref_added
                    delete!(dh1, mission_reference_for_amplitude_normalization)
                    delete!(params_fill, mission_reference_for_amplitude_normalization)
                end
            end

            if plot_dh_as_function_of_time_and_elevation
                Altim.plot_height_time(dh1; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="filled_ampnorm", fig_folder, figure_suffix, mask=param.surface_mask, showplots)
            end

            dh1 = Altim.hyps_fill_empty!(dh1, params_fill, geotiles; mask=param.surface_mask)

            dh1 = Altim.hyps_fill_updown!(dh1, geotiles; mask=param.surface_mask)

            if plot_dh_as_function_of_time_and_elevation
                Altim.plot_height_time(dh1; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="filled_updown", fig_folder, figure_suffix, mask=param.surface_mask, showplots)
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
                params_fill = merge(params_fill, model_param)
            end

            # save filled geotiles
            save(binned_filled_file, Dict("dh_hyps" => dh1, "nobs_hyps" => nobs1, "model_param" => params_fill))
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

    params = NamedTuple[]
    for binned_folder in binned_folders
        for paramater_set in paramater_sets
            for surface_mask in surface_masks
                for dem_id in dem_ids
                    for binning_method in binning_methods
                        for curvature_correct in curvature_corrects
                            for amplitude_correct = amplitude_corrects
                                push!(params, (; project_id, binned_folder, paramater_set, surface_mask, dem_id, binning_method, curvature_correct, amplitude_correct))
                            end
                        end
                    end
                end
            end
        end
    end

    Threads.@threads for param in params

        geotiles = copy(Altim.geotiles_mask_hyps(param.surface_mask, geotile_width))

        # make geotile rgi regions mutually exexclusive 
        geotiles, reg = Altim.geotiles_mutually_exclusive_rgi!(geotiles)

        fig_folder = joinpath(param.binned_folder, "figures")

        # paths to files
        binned_filled_file, figure_suffix = Altim.binned_filled_filepath(; param...)

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
function old_geotile_bin_fit_fac(;
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


function geotile_regional_dv(;
    project_id = :v01,
    geotile_width = 2,
    surface_masks = [:glacier, :glacier_rgi7, :glacier_b1km, :glacier_b10km],
    binned_folders = ("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"),
    dem_ids = [:best, :cop30_v2],
    binning_methods = ["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"],
    curvature_corrects = [true, false],
    paramater_sets = [1, 2],
    amplitude_corrects = [true, false],
    force_remake = true,
    )

    # compute regional volume change, firn correction and mass change
    params = NamedTuple[]
    for binned_folder in binned_folders
        for paramater_set in paramater_sets
            for surface_mask in surface_masks
                for dem_id in dem_ids
                    for binning_method in binning_methods
                        for curvature_correct in curvature_corrects
                            for amplitude_correct = amplitude_corrects
                                push!(params, (; project_id, binned_folder, paramater_set, surface_mask, dem_id, binning_method, curvature_correct, amplitude_correct))
                            end
                        end
                    end
                end
            end
        end
    end

    Threads.@threads for param in params
    #param = first(params)

        # paths to files
        binned_filled_file, figure_suffix = Altim.binned_filled_filepath(; param...)
        binned_aligned_file = replace(binned_filled_file, ".jld2" => "_aligned.jld2")
        dv_reg_file = replace(binned_aligned_file, ".jld2" => "_reg.jld2")


        if (!isfile(dv_reg_file) || force_remake) && (isfile(binned_aligned_file))
        
            dh = FileIO.load(binned_aligned_file, "dh_hyps")
            nobs0 = FileIO.load(binned_aligned_file, "nobs_hyps")

            geotiles = copy(Altim.geotiles_mask_hyps(param.surface_mask, geotile_width))

            # make geotile rgi regions mutually exexclusive 
            geotiles, reg = Altim.geotiles_mutually_exclusive_rgi!(geotiles)

            if !any(contains.("area_km2", names(geotiles)))
                original_area_name = string(param.surface_mask) * "_area_km2"
                generic_area_name = "area_km2"
                rename!(geotiles, original_area_name => generic_area_name)
            end

            # align geotile dataframe with DimArrays [this is overly cautious]
            gt = collect(dims(dh[first(keys(dh))], :geotile))
            gt_ind = [findfirst(geotiles.id .== g0) for g0 in gt]
            geotiles0 =geotiles[gt_ind, :]

            dv_reg, nobs_reg, area_reg = Altim.geotile_filled_dv_reg(dh, nobs0, geotiles0, reg)

            dv = Dict("dv" => dv_reg, "nobs" => nobs_reg, "area" => area_reg)

            save(dv_reg_file, dv)
            println("regional volume change complete: $(param)")
        end
    end
end

function geotile_dvdm_synthesize_old(;
    # best estimate ,
    project_id = :v01,
    surface_mask_best="glacier",
    dem_best="best",
    curvature_correct_best=true,
    amplitude_correct_best=true,
    binning_method_best="meanmadnorm3",
    fill_param_best=1,
    binned_folder_best="/mnt/bylot-r3/data/binned_unfiltered/2deg",

    # to include in uncertainty
    surface_masks=["glacier" "glacier_rgi7" "glacier_b1km"],
    dems=["best" "cop30_v2"],
    curvature_corrects=[false true],
    amplitude_corrects=[true],
    binning_methods=["median" "meanmadnorm5" "meanmadnorm3" "meanmadnorm10"],# ["median" "meanmadnorm10" "meanmadnorm5" "meanmadnorm3"]
    fill_params=[1 2],
    binned_folders=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"),

    # manual adjustments
    regions_to_overwrite_hugonnet_data_with_model_fit=["rgi19"],

    ## combine regions at the per-sensor level
    # this gets a bit complicated a only regions with consistnet mission inclusion can be grouped for proper 
    region_col="rgi",
    region_combines_before_fac=(
        ("rgi13", "rgi14", "rgi15") => "hma", ), 
    
    region_combines_after_fac=(
    ("rgi1", "rgi3", "rgi4", "rgi5", "rgi6", "rgi7", "rgi8", "rgi9", "rgi19") => "hll",
    ("rgi2", "rgi10", "rgi11", "rgi12", "hma", "rgi17", "rgi18") => "ghll",
    ("rgi1", "rgi3", "rgi4", "rgi6", "rgi7", "rgi8", "rgi9") => "hll_ep",),

    combine_vars=["area_km2", "val", "nobs"],

    missions_to_filter_on_nobs_km2=["icesat", "gedi"],

    path2gemb = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0_geotile_filled_d_reg.jld2",
    )

    pscale_range = nothing #(0.75, 1.25)
    #pscale_range_rgi_exclude = ["rgi12", "rgi16"]

    mission_firn = "gemb"

    df = Altim.regionfiles2dataframe(;
    binned_folders,
    surface_masks,
    dem_ids=dems,
    binning_methods,
    curvature_corrects,
    amplitude_corrects,
    paramater_sets=fill_params,
    project_id,
    )

    # this is a sanity check to make sure that the "best" estimate exists
    # compute binary indecies into best estimate
    param_valid_best = (df.surface_mask .== surface_mask_best) .&
                       (df.binning_method .== binning_method_best) .&
                       (df.dem .== dem_best) .&
                       (df.curvature_correct .== curvature_correct_best) .&
                       (df.amplitude_correct .== amplitude_correct_best) .&
                       (df.fill_param .== fill_param_best) .&
                       (df.binned_folder .== binned_folder_best)

    if !any(param_valid_best)
        error("no `best` estimate found for any rgi, any mission or any variable, are you sure that the `best` estimate that you provided exists? Maybe you have an incorrectly spelled parameter.")
    end

    # load all regional gemb data
    df_gemb = Altim.gemb2dataframe(; path2file=path2gemb)
    
    # densify between gemb classes
    df_gemb = Altim.gemb_classes_densify!(df_gemb; n_densify=4)

    # exclude some pscale values... this is done as an increase in melt can be offset by an 
    # equal increase in prcipitation, creating non-unique solutions that can lead to large 
    # scalings in both precipitation on melt. 
    if .!isnothing(pscale_range)
        index = (df_gemb.pscale .>= minimum(pscale_range)) .& (df_gemb.pscale .<= maximum(pscale_range))
        df_gemb = df_gemb[index, :]
    end

    # load all regional gemb data
    df_gemb = Altim.gemb2dataframe(; path2file=path2gemb)

    # add fac + smb - dischage to get gemb elevation change
    df_gemb = Altim.gem_Δvolume!(df_gemb)

    # df_gemb and df have may have different date ranges... therefore they need some special attention to append. 
    df = Altim.dvdm_append(df, df_gemb)

    missions = unique(df.mission)
    dates = DataFrames.metadata(df, "date")
    decyear = Altim.decimalyear.(dates)
    dems = unique(df.dem)

    rgis = unique(df.rgi)
    n = nrow(df)

    # compute data density
    df[!, "nobs_km2"] .= 0.0
    for r in eachrow(df)
        hasobs = r.nobs .> 0
        r.nobs_km2 = mean(r.nobs[hasobs] ./ r.area_km2)
    end
    df[isnan.(df.nobs_km2), "nobs_km2"] .= 0.0

    # create an rgi surface_mask area data frame
    # round area to the nearest 1/100 km2 to avoid issues with numeric precision
    df.area_km2 = round.(df.area_km2, digits = 2)
    x = df[((df.mission.=="icesat").&(df.surface_mask.==surface_mask_best)), :]
    idx = unique(z -> x.area_km2[z], 1:length(x.area_km2))
    df_area = DataFrame(rgi=x[idx, :rgi], area_km2=x[idx, :area_km2])

    # find index where each mission has valid data... choose rgi2 since all missions are valid there
    mission_ts_valid = Dict()
    for mission in missions
        valid = .!isnan.(df.val[findfirst((df.mission .== mission) .& (df.rgi .== "rgi2"))])
        mission_ts_valid[mission] = valid
    end

    # exclude data that does not meet nobs per area criteria [check across rgi and mission]
    nobs_km2_threshold = 0.01

    for rgi in rgis
        #rgi = "rgi1"
        rgi_index = (df.rgi .== rgi)

        for mission in missions_to_filter_on_nobs_km2
            index = rgi_index .& (df.mission .== mission)
            if mean(df.nobs_km2[index]) < nobs_km2_threshold
                for r in eachrow(@view df[index, :])
                    for col in ["val"]
                        r[col] .= NaN
                    end
                    
                    r["nobs"] .= 0.0
                end
            end
        end
    end

    # print excluded regions
    printstyled("--------------------------------------------------------------------------------------------\n"; color=:blue)
    printstyled("                REGIONS WITH POOR MISSION COVERAGE EXCLUDED FROM ANALYSIS\n"; color=:blue)
    printstyled("--------------------------------------------------------------------------------------------\n"; color=:blue)

    for mission in missions
        b = falses(size(rgis))
        for (i, rgi) in enumerate(rgis)
            rgi_index = (df.rgi .== rgi)

            if mean(df.nobs_km2[(df.mission.==mission).&rgi_index]) >= nobs_km2_threshold
                b[i] = true
            end
        end

        if !any(.!b)
            printstyled("   $(mission): "; color=:blue)
            printstyled("None\n"; color=:green)
        else
            printstyled("   $(mission): "; color=:blue)
            printstyled("$(rgis[.!b])\n"; color=:red)
        end
    end

    printstyled("--------------------------------------------------------------------------------------------\n\n"; color=:blue)

    # [1] Set the mean over the icesat2 period to zero by removing the mean offset to the 
    # reference surface and remove that same offset from icesat.
    # [2] Align GEDI and ASTER to the median icesat + icesat2 offset
    mission_ref1 = "icesat2"
    mission_ref2 = "icesat"
    mission_index = Dict()

    mission_df = Dict()
    for mission in missions
        mission_index[mission] = df.mission .== mission
        mission_df[mission] = @view df[mission_index[mission], :]
    end

    n = nrow(mission_df[mission_ref1])
    
    #=
    for rref1 in eachrow(mission_df[mission_ref1])
        #rref1 = first(eachrow(mission_df[mission_ref1]))

        # find matching icesat2 row
        match_index = Dict()
        for mission in setdiff(missions, [mission_ref1])
            match_index[mission] = trues(n)
            for col in ["var", "rgi", "surface_mask", "dem", "curvature_correct", "binning_method", "project_id", "amplitude_correct", "fill_param", "binned_folder"]
                match_index[mission] = match_index[mission] .& (rref1[col] .== mission_df[mission_ref2][:, col])
            end

            if sum(match_index[mission]) !== 1
                error("non unique match: $mission")
            end

            match_index[mission] = findfirst(match_index[mission])
        end

        rref2 = @view mission_df[mission_ref2][match_index[mission_ref2], :]

        # [1] Set the mean over the icesat2 period to zero by removing the mean offset to the 
        # reference surface and remove that same offset from icesat.
        for col in ["val"]
            offset = Altim.nanmean(rref1[col])

            # remove offset from ICESat2
            rref1[col] = rref1[col] .- offset

            #remove same offset from ICESat
            rref2[col] = rref2[col] .- offset
        end

        # [2] Align other missions to median icesat + icesat2 offset        
        for mission in setdiff(missions, [mission_ref2, mission_ref1])
            dfX = @view mission_df[mission][match_index[mission], :]

            for col in ["val"]
                offset = Altim.nanmedian(hcat(dfX[col] .- rref1[col], dfX[col] .- rref2[col]))
                dfX[col] = dfX[col] .- offset
            end
        end
    end
    =#

    # align firn variables to a common time stamp [alignment to icesat and icesat2 is only done for the common variable "dv_km3"]
    # set mean over the icesat and icesat2 period to zero
    vars = unique(df[df.mission.==mission_firn, :var])

    for r in eachrow(mission_df[mission_firn])

        if !any(vars .== r.var)
            continue
        end

        valid = (mission_ts_valid[mission_ref1] .| mission_ts_valid[mission_ref2])
        r.val = r.val .- mean(r.val[valid])
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

                for col = ["val"]

                    y = fill(NaN, (length(r[col]), length(fit_missions)))
                    for (i, mission) in enumerate(fit_missions)
                        y[:, i] = mission_df[mission][match_index[mission], col][1]
                    end

                    y = [Altim.nanmean(foo) for foo in eachrow(y)]

                    notnan = .!isnan.(y)
                    x = decyear .- mean(decyear)
                    foo_fit = Altim.curve_fit(Altim.model3, x[notnan], y[notnan], Altim.p3)

                    notnan = mission_ts_valid["hugonnet"]
                    y .= NaN

                    y[notnan] = Altim.model3(x[notnan], foo_fit.param)

                    mission_df["hugonnet"][k, col] = y
                end
            end
        end
    end

    # individual HMA regions (rgi13, rgi14, rgi15) do not provide good callibration to smb...
    # to help with this combine hma region before fac callibration ... other regions after
    set2nan = ["nobs_km2"] # need"Δheight", "pscale" to combine gemb variables

    df = Altim.region_combine!(df; region_col, region_combines=region_combines_before_fac, combine_vars, set2nan)

    # recompute mission index
    missions = unique(df.mission)
    mission_index = Dict()
    mission_df = Dict()
    for mission in missions
        mission_index[mission] = df.mission .== mission
        mission_df[mission] = @view df[mission_index[mission], :]
    end

    # calibrate firn model
    mission_df["gemb_calibrate"] = DataFrame()
    mission_df["dm_gt"] = DataFrame()

    match_parameters = ["rgi", "surface_mask", "dem", "curvature_correct", "binning_method", "project_id", "amplitude_correct", "fill_param", "binned_folder"]

    rgis = unique(df.rgi)
    vars = unique(mission_df[mission_firn].var)
    n = nrow(mission_df["icesat2"])
    # iterate through unique model runs (multiple missions)

    df_tmp_plot = DataFrame()
    for r in eachrow(mission_df["icesat2"])
        #r = eachrow(mission_df["icesat2"])[mission_df["icesat2"].rgi .== "rgi2"][10]

        gemb_dv_km3 = mission_df[mission_firn][(mission_df[mission_firn].var.=="dv_km3").&(mission_df[mission_firn].rgi.==r.rgi), :]

        # find matching icesat2 row
        dv = r.val
        match_index = Dict()
        for mission in setdiff(missions, [mission_firn])
            #mission = first(setdiff(missions, [mission_firn]))

            match_index[mission] = trues(n)
            for col in match_parameters
                match_index[mission] = match_index[mission] .& (r[col] .== mission_df["icesat"][:, col])
            end

            if sum(match_index[mission]) !== 1
                error("non unique match: $mission")
            end

            match_index[mission] = findfirst(match_index[mission])
            dv = Altim.nanmean.(eachrow(hcat(dv, mission_df[mission].val[match_index[mission]])))

        end

        se = fill(Inf, nrow(gemb_dv_km3))
        notnan = (.!isnan.(dv)) .& (.!isnan.(gemb_dv_km3.val[1]))
        dv = dv .- mean(dv[notnan])

        for (i, r0) in enumerate(eachrow(gemb_dv_km3))
            #(i, r)  = first(enumerate(eachrow(gemb_dv_km3)))
            gemb_dv0 = (r0.val .- mean(r0.val[notnan]))
            se[i] = std(dv[notnan] .- gemb_dv0[notnan])
        end

        # record best precipitation and delta elevation
        index = findfirst(se .== minimum(se))


        df_tmp_plot = append!(df_tmp_plot, DataFrame(rgi=r.rgi, dv=[dv], gemb=[gemb_dv_km3.val[index]], Δheight=gemb_dv_km3[index, :Δheight][1], pscale=gemb_dv_km3[index, :pscale][1]))


        for mission in setdiff(missions, [mission_firn])
            mission_df[mission][match_index[mission], :Δheight] = gemb_dv_km3[index, :Δheight]
            mission_df[mission][match_index[mission], :pscale] = gemb_dv_km3[index, :pscale]
        end

        for mission in setdiff(missions, [mission_firn])
            #mission = first(setdiff(missions, [mission_firn]))

            foo = DataFrame(mission_df[mission][match_index[mission], :]) #wrapped in DataFrame() to copy row
            fac = mission_df[mission_firn][(mission_df[mission_firn].var.=="fac_km3").&(mission_df[mission_firn].rgi.==r.rgi), :][index, :]
            foo.val[1] = foo.val[1] .- fac.val
            foo.var[1] = "dm_gt"
            mission_df["dm_gt"] = append!(mission_df["dm_gt"], DataFrame(foo))
        end

        # add metadata from altim runs to best fit gemb

        # append other variables 
        for var0 in vars
            #var0 = first(vars)

            foo = DataFrame(mission_df[mission_firn][(mission_df[mission_firn].var.==var0).&(mission_df[mission_firn].rgi.==r.rgi), :][index, :])
            for col in names(foo)
                #col = "dem"

                if (ismissing(foo[1, col]) .| (col == "surface_mask")) #surface_mask need to be overwritten to make smb unique
                    foo[1, col] = r[col]
                end
            end
            mission_df["gemb_calibrate"] = append!(mission_df["gemb_calibrate"], foo)
        end
    end

    # check that everthing is working as expected
    if false
        for rgi in rgis
            index = df_tmp_plot.rgi .== rgi
            foo = df_tmp_plot[index, :]
            for r = eachrow(foo)
                tmp = mission_df[mission_firn][(mission_df[mission_firn].rgi.==rgi).&(mission_df[mission_firn].var.=="dv_km3"), :]
                title = "$rgi: Δheight: $(r.Δheight), pscale: $(r.pscale)"

                p = plot(decyear, tmp.val; size=(1600, 800), legend=false, title)
                plot!(decyear, r.gemb, linewidth=3, color=:black)
                plot!(decyear, r.dv .- Altim.nanmean(r.dv .- r.gemb), linewidth=3, color=:red)

                display(p)
            end
        end
    end

    # only keep best fac estimate for each altim model run (same for each sensor)
    missions = setdiff(keys(mission_df), [mission_firn])
    df0 = DataFrame()
    for mission in missions
        df0 = append!(df0, copy(mission_df[mission]))
    end

    df = df0
    set2nan = ["Δheight", "pscale", "nobs_km2"] # "Δheight", "pscale" must be excluded after best gemb simulations selected
    df = Altim.region_combine!(df; region_col, region_combines=region_combines_after_fac, combine_vars, set2nan)

    # after combining regions, a number of helper variables need to be recomputed
    begin
        # recompute mission index
        missions = unique(df.mission)
        mission_index = Dict()
        mission_df = Dict()
        for mission in missions
            mission_index[mission] = df.mission .== mission
            mission_df[mission] = @view df[mission_index[mission], :]
        end

        # compute binary indecies into best estimate
        param_valid_best = (df.surface_mask .== surface_mask_best) .&
                            (df.binning_method .== binning_method_best) .&
                            (df.dem .== dem_best) .&
                            (df.curvature_correct .== curvature_correct_best) .&
                            (df.amplitude_correct .== amplitude_correct_best) .&
                            (df.fill_param .== fill_param_best) .&
                            (df.binned_folder .== binned_folder_best)

        if !any(param_valid_best)
            error("no `best` estimate found for any rgi, any mission or any variable, are you sure that the `best` estimate that you provided exists?")
        end

        # create an rgi surface_mask area data frame
        x = df[((df.mission.=="icesat").&(df.surface_mask.==surface_mask_best)), :]
        idx = unique(z -> x.area_km2[z], 1:length(x.area_km2))
        df_area = DataFrame(rgi=x[idx, :rgi], area_km2=x[idx, :area_km2])

        # compute data density
        df[!, "nobs_km2"] .= 0.0
        for r in eachrow(df)
            hasobs = r.nobs .> 0
            r.nobs_km2 = mean(r.nobs[hasobs] ./ r.area_km2)
        end

        df[isnan.(df.nobs_km2), "nobs_km2"] .= 0.0
    end

    # compute central estimate and lower and upper bounds for each mission
    vars = unique(df.var)
    rgis = unique(df.rgi)
    quantiles = [0.95]

    df_mission = DataFrame()
    DataFrames.metadata!(df_mission, "date", collect(dates); style=:note)

    for rgi in rgis
    #rgi = "rgi1"

        for mission in missions
            #mission = "hugonnet"
            for var0 in vars
                #var0 = "dm_gt"

                df1 = DataFrame()
                DataFrames.metadata!(df1, "date", collect(dates); style=:note)

                valid = findall((df.rgi .== rgi) .& (mission_index[mission]) .& (df.var .== var0))
                valid_best = (df.rgi .== rgi) .& mission_index[mission] .& param_valid_best .& (df.var .== var0)

                if sum(valid_best) == 0
                    continue
                elseif sum(valid_best) != 1
                    error("non unique best estimate")
                end

                valid_best = findfirst(valid_best)

                # central estimate is user selected best estimate
                df1[!, :mid] = [df[valid_best, :val]] # this selects the best estimate as the median

                # all estimates
                f = reduce(hcat, df[valid, :val])
                f[f.==0] .= NaN

                # remove central estimate to get anomaly
                f = f .- repeat(df1[!, :mid][1], 1, size(f, 2))

                # take absolute as I want a "fair" symmetric error
                f = abs.(f)

                q = reduce(hcat, Altim.nanquantile.(eachrow(f), Ref(quantiles)))
                sigma = q[1, :]

                # make error symetric
                df1[!, :low] = [df1[!, :mid][1] .- sigma]
                df1[!, :high] = [df1[!, :mid][1] .+ sigma]

                df1[!, "nobs"] .= [df[valid_best, "nobs"]]
                df1[!, "mission"] .= mission
                df1[!, "rgi"] .= rgi
                df1[!, "var"] .= var0
                df1[!, "area_km2"] .= df_area[df_area.rgi.==rgi, :area_km2]
                append!(df_mission, df1)
            end
        end 
    end

    ## synthesize observations (error weighted average)
    df_synth = DataFrame()
    DataFrames.metadata!(df_synth, "date", collect(dates); style=:note)
    missions_synthesize = setdiff(unique(df_mission.mission), [mission_firn])
    vars_synthesize = ["dv_km3", "dm_gt"]

    for var0 in vars_synthesize
    #var0 = "dm_gt"

        var_index = (df_mission.var .== var0)

        for rgi in rgis
            #rgi = "rgi19"

            rgi_index = (df_mission.rgi .== rgi)

            foo_sum = zeros(size(decyear))
            foo_weight = zeros(size(decyear))
            foo_err = zeros(size(decyear))
            foo_n = zeros(size(decyear))

            for mission in missions_synthesize
                mission_index0 = (df_mission.mission .== mission)

                index = var_index .& rgi_index .& mission_index0

                if !any(index)
                    continue
                end

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

            low = mid .- error0
            high = mid .+ error0

            index = var_index .& rgi_index
            nobs0 = vec(sum(hcat(df_mission[index, "nobs"]...), dims=2))

            foo = DataFrame(df_mission[findfirst(index), :]) # copy a row

            foo.low[1] = low
            foo.high[1] = high
            foo.mid[1] = mid
            foo.nobs[1] = nobs0
            foo.mission[1] = "synthesis"
            foo[!, :area_km2] .= df_area[df_area.rgi.==rgi, :area_km2]

            append!(df_synth, foo)
        end
    end

    # combine mission and synthesis results
    df_dmass = vcat(df_mission, df_synth)
    DataFrames.metadata!(df_dmass, "date", collect(dates); style=:note)

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

    # align gemb records: subtract the median offset to synthesis over the icesat and icesat-2 period
    valid = (mission_ts_valid[mission_ref1] .| mission_ts_valid[mission_ref2])
    for r in eachrow(df_dmass)

        if r.mission != mission_firn
            continue
        end

        if r.unit == "km⁻³"
            var0 = "dv"
        elseif r.unit == "Gt"
            var0 = "dm"
        end

        index_ref = (df_dmass.mission .== "synthesis") .& (df_dmass.var .== var0) .& (df_dmass.rgi .== r.rgi)
        
        ref0 = df_dmass[index_ref, :mid][1]

        index = .!isnan.(ref0) .& .!isnan.(r.mid) .& valid
        offset = median((r.mid[index] .- ref0[index]))
        r.mid = r.mid .- offset
        r.low = r.low .- offset
        r.high = r.high .- offset
    end

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

    missions = unique(df_dmass[df_dmass.var.=="dm", :mission])
    # find index where each mission has valid data... choose rgi2 since all missions are valid there
    mission_ts_valid = Dict()
    for mission in missions
        valid = .!isnan.(df_dmass.mid[findfirst((df_dmass.mission .== mission) .& (df_dmass.var .== "dm") .& (df_dmass.rgi .== "rgi2"))])
        mission_ts_valid[mission] = valid
    end

    x = df_dmass[(df_dmass.mission.=="icesat2"), :]
    idx = unique(z -> x.area_km2[z], 1:length(x.area_km2))
    df_area = DataFrame(rgi=x[idx, :rgi], area_km2=x[idx, :area_km2])

    # find valid range for ALL grace data (needed when substituting in missing reigons with Synthesis data)
    datenum = grace[first(keys(grace))]["dM_gt_mdl_fill_date"]
    decyear0 = vec(Altim.decimalyear.(Altim.datenum2date.(datenum)))

    dm_gt = vec(grace[first(keys(grace))]["dM_gt_mdl_fill"])
    dm_sigma = vec(grace[first(keys(grace))]["dM_sigma_gt_mdl_fill"])
    notnan = .!isnan.(dm_gt)

    # introplate to df time
    valid = (decyear .<= maximum(decyear0[notnan])) .& (decyear .>= minimum(decyear0[notnan]))

    for rgi in rgis
       

        if any(keys(grace) .== rgi)
            
            dm_gt = vec(grace[rgi]["dM_gt_mdl_fill"])
            dm_sigma = vec(grace[rgi]["dM_sigma_gt_mdl_fill"])
            
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

            # check that interpoalted results are consistent with orginal data
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

            printstyled("--------------------------------------------------------------------------------------------\n"; color=:blue)
            printstyled(" no GRACE solution found for $(rgi): substituting in synthesis record \n"; color=:red)
            printstyled("--------------------------------------------------------------------------------------------\n"; color=:blue)


            index = (df_dmass.var .== "dm") .& (df_dmass.mission .== "synthesis") .& (df_dmass.rgi .== rgi)

            # if region does not exist, substitute with synthesis data
            mid = fill(NaN, length(valid))
            low = fill(NaN, length(valid))
            high = fill(NaN, length(valid))

            mid[valid] = df_dmass.mid[index][1][valid]
            low[valid] = df_dmass.low[index][1][valid]
            high[valid] = df_dmass.high[index][1][valid]

            # extend record for missing years
            ind1 = findfirst(isnan.(mid[valid])) + findfirst(valid) -1
            ind2 = findlast(valid)
            if ind1 < ind2
                mid[ind1:ind2] .= mid[ind1-1]
                low[ind1:ind2] .= low[ind1-1]
                high[ind1:ind2] .= high[ind1-1]
            end
        end

        index = (df_dmass.var .== "dm") .& (df_dmass.mission .== "synthesis") .& (df_dmass.rgi .== rgi)

        # align records: subtract the median offset to synthesis over the icesat and icesat-2 period 
        valid_ref = (mission_ts_valid["icesat"] .| mission_ts_valid["icesat2"]) .& .!isnan.(mid)
        mid_synth = df_dmass.mid[index][1]
        offset = Altim.nanmean(mid[valid_ref] .- mid_synth[valid_ref])

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

    return df_dmass
end

## combine synthisized regions [can now combine regions covered by different sensors]
function geotile_dvdm_global!(df_dmass)
    region_col = "rgi"
    combine_vars = ["area_km2", "mid", "nobs"]
    rss_vars = ["low", "high"]
    set2nan = nothing


    # to RSS variables you need to frist subtract :mid from :low and high
    for r in eachrow(df_dmass)
        r.low = r.mid .- r.low
        r.high = r.high .- r.mid
    end

    # combine synthesis data for global mass change
    begin
        missions = setdiff(unique(df_dmass.mission), ["grace"])
        combined_regs = ["rgi16", "hll", "ghll"]
        region_combines = (
            tuple(combined_regs...) => "global",
        )

        df_dmass = Altim.region_combine!(df_dmass; region_col, region_combines, combine_vars, rss_vars, missions, set2nan)
    end

    # combine synthesis data for global mass change excluding periphery
    begin
        missions = setdiff(unique(df_dmass.mission), ["grace"])
        combined_regs = ["rgi16", "hll_ep", "ghll"]
        region_combines = (
            tuple(combined_regs...) => "global_ep",
        )

        df_dmass = Altim.region_combine!(df_dmass; region_col, region_combines, combine_vars, rss_vars, missions, set2nan)
    end

    # combine grace data for global mass change excluding periphery
    begin
        missions = ["grace"]
        rgis = reduce(vcat, (["rgi$i" for i in 1:4], ["rgi$i" for i in 6:18]))
        region_combines = (
            tuple(rgis...) => "global_ep",
        )
        df_dmass = Altim.region_combine!(df_dmass; region_col, region_combines, combine_vars, rss_vars, missions, set2nan)
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


    df_dmass = Altim.remove_allnan(df_dmass, :mid)

    return df_dmass
end


## add estimate of trend, acceleration and uncertainty
function geotile_dvdm_add_trend!(df_dmass; iterations=1000)

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
    return df_dmass
end


# convert from Gt/km3 to m.w.e/m
function geotile_dvdm_areaaverage(df_dmass)

    df_dheight = deepcopy(df_dmass)
    for (i, r) in enumerate(eachrow(df_dheight))
        for v in ["mid", "low", "high", "trend", "trend_err", "acceleration", "acceleration_err", "amplitude", "amplitude_err"]
            df_dheight[i, v] = r[v] / r.area_km2 * 1000
        end
        df_dheight[i, "unit"] = Altim.unit2areaavg[r.unit]
    end

    return df_dheight

end



# adjust dh estimates according to anomalies measured over land and seasonality from a referecne mission. 
function geotile_adjust!(
    dh,
    dh_land,
    nobs_land,
    dh_land_ref;
    minnobs=45,
    ref_minvalid=12,
    ref_madnorm_max=10,
)

    #dh_land0 = FileIO.load(binned_file_land, "dh_hyps")
    #nobs_land = FileIO.load(binned_file_land, "nobs_hyps")
    #dh = FileIO.load(binned_file, "dh_hyps")
    #dh = dh["hugonnet"],
    #dh_land = dh_land0["hugonnet"],
    #nobs_land = nobs_land["hugonnet"],
    #dh_land_ref = dh_land0["icesat2"]

    dh = copy(dh)

    ddates = dims(dh_land_ref, :date)
    dates = collect(ddates)
    dheight = dims(dh_land_ref, :height)
    dgeotile = dims(dh_land_ref, :geotile)

    Δstd = fill(NaN, dgeotile)
    amplitude = fill(NaN, dgeotile)

    decimalyear = Altim.decimalyear.(dates)
    Δdecimalyear = decimalyear .- mean(decimalyear)

    # loop through each geotile
    Threads.@threads for i in eachindex(dgeotile)
        #i = 1000    
        # time series over land
        ts_ref = Altim.nanmedian.(eachrow(dh_land_ref[i, :, :]))
        ts = Altim.nanmedian.(eachrow(dh_land[i, :, :]))
        ts_nobs = sum(nobs_land[i, :, :], dims=:height)

        ts[vec(ts_nobs).<minnobs] .= NaN

        # fit a seasonal cycle to icesat2
        valid = .!isnan.(ts_ref)
        if (sum(valid) < ref_minvalid) .| all(isnan.(ts))
            continue
        end

        valid[valid] = Altim.madnorm(ts_ref[valid]) .< ref_madnorm_max

        fit1 = curve_fit(Altim.model4, Δdecimalyear[valid], ts_ref[valid], Altim.p4)
        amplitude[i] = fit1.param[2]

        ts_fit = Altim.model4(Δdecimalyear, fit1.param)

        # determine correction
        ts_correction = ts .- ts_fit

        # fill nan's with median correction
        valid = .!isnan.(ts_correction)
        ts_correction[.!valid] .= median(ts_correction[valid]) .- ts_fit[.!valid]  # missing seasonalyity in hugonnet needs to be added back

        ts_unadjusted = Altim.nanmedian.(eachrow(dh[i, :, :]))

        correction = hcat(ts_correction.data) * ones(1, length(dheight))

        dh[i, :, :] -= correction

        ts_adjusted = Altim.nanmedian.(eachrow(dh[i, :, :]))

        Δstd[i] = std(ts_adjusted[.!isnan.(ts_adjusted)]) - std(ts_unadjusted[.!isnan.(ts_unadjusted)])
    end

    return dh

    if false
        # number of tiles adjusted 
        p = Plots.histogram(Δstd; title="Δstd")
        display(p)
        Plots.histogram(abs.(amplitude); title="icesat2 amplitude over land")
        display(p)

        println("numer of tiles adjusted = $(sum(.!isnan.(Δstd))): $((;minnobs, minvalid))")
        println("median Δstd = $(round(Altim.nanmedian(Δstd), digits = 2))")
        println("median amplitude = $(round(Altim.nanmedian(amplitude), digits = 2))")
        println("mean Δstd = $(round(Altim.nanmean(Δstd), digits = 2))")
        println("mean amplitude = $(round(Altim.nanmean(amplitude), digits = 2))")
        println()
        println()
    end

end

offset_trend::Function = offset_trend(t, p) = p[1] .+ p[2] .* t;
offset_trend_p = zeros(2);

function geotile_align!(
    dh;
    mission_ref1="icesat2",
    mission_ref2="icesat",
    min_trend_count=5,
    remove_land_surface_trend=Altim.mission_land_trend(),
    showplots=false,
)

    missions = collect(keys(dh))
    ddata = dims(dh[missions[1]], :date)
    dheight = dims(dh[missions[1]], :height)
    dgeotile = dims(dh[missions[1]], :geotile)

    decyear = Altim.decimalyear.(ddata)
    Δdecyear = decyear .- mean(decyear)

    offset = Dict()
    trend = Dict()
    for mission in missions
        offset[mission] = zeros(dgeotile, dheight)
        trend[mission] = zeros(dgeotile, dheight)
    end

    if showplots
        p = Plots.plot(dh["hugonnet"][700, :, At(4050)]);
        Plots.plot!(dh["icesat"][700, :, At(4050)])
        Plots.plot!(dh["icesat2"][700, :, At(4050)])
    end

    Threads.@threads for geotile in dgeotile
        #geotile = dgeotile[700]
        (ridx, cidx) = Altim.validrange(.!isnan.(dh[mission_ref1][At(geotile), :, :]))
        offset_ref = dropdims(median(dh[mission_ref1][At(geotile), ridx, cidx], dims=:date), dims=:date)
        #fit = curve_fit(Altim.offset_trend, collect(dheight[cidx]), offset_ref[:], Altim.offset_trend_p)
        #offset_ref = DimArray(Altim.offset_trend(dheight, fit.param), dheight)

        for height in dheight[cidx]
            #height = 4050;

            ref1 = dh[mission_ref1][At(geotile), :, At(height)]
            ref2 = dh[mission_ref2][At(geotile), :, At(height)]
            valid_ref1 = .!isnan.(ref1)
            valid_ref2 = .!isnan.(ref2)

            ref = copy(ref2)
            ref[valid_ref1] = ref1[valid_ref1]

            # set reference time period to zero (anomalies are calculated)
            ref .-= offset_ref[At(height)]
            dh[mission_ref1][At(geotile), :, At(height)] .-= offset_ref[At(height)]
            dh[mission_ref2][At(geotile), :, At(height)] .-= offset_ref[At(height)]

            # align overlap in mean and trend
            for mission in setdiff(missions, hcat(mission_ref1, mission_ref2))
                #mission = "hugonnet"

                delta = dh[mission][At(geotile), :, At(height)] .- ref
                valid_delta = .!isnan.(delta)

                if !isnothing(remove_land_surface_trend)
                    trend[mission][At(geotile), At(height)] += remove_land_surface_trend[At(mission)]
                    delta .-= (trend[mission][At(geotile), At(height)] .* Δdecyear)
                end

                if .!any(valid_delta)
                    continue
                end

                offset[mission][At(geotile), At(height)] += median(delta[valid_delta])

                if (sum(valid_ref1 .& valid_delta) >= min_trend_count) && (sum(valid_ref2 .& valid_delta) >= min_trend_count)
                    # correct offset and slope

                    ## we could weight by observation here, but for now just do simply fit ##
                    fit = curve_fit(offset_trend, Δdecyear[valid_delta], (delta[valid_delta] .- offset[mission][At(geotile), At(height)]), offset_trend_p)

                    #offset_trend = delta[valid_delta] \ hcat(ones(sum(valid_delta)))
                    offset[mission][At(geotile), At(height)] += fit.param[1]
                    trend[mission][At(geotile), At(height)] += fit.param[2]
                end

                dh[mission][At(geotile), :, At(height)] .-= offset[mission][At(geotile), At(height)]
                dh[mission][At(geotile), :, At(height)] -= (trend[mission][At(geotile), At(height)] .* Δdecyear)
            end
        end
    end

    if showplots
        Plots.plot!(dh["hugonnet"][700, :, At(4050)]);
        Plots.plot!(dh["icesat"][700, :, At(4050)]);
        Plots.plot!(dh["icesat2"][700, :, At(4050)]);
        display(p)
    end

    return dh
end

function geotile_align_replace(;
    mission_ref1="icesat2",
    mission_ref2="icesat",
    min_trend_count=5,
    remove_land_surface_trend=Altim.mission_land_trend(),
    project_id=:v01,
    surface_masks=[:glacier, :glacier_rgi7, :glacier_b1km, :glacier_b10km],
    binned_folders=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"),
    dem_ids=[:best, :cop30_v2],
    binning_methods=["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"],
    curvature_corrects=[true, false],
    paramater_sets=[1, 2],
    amplitude_corrects=[true, false],
    geotile_width = 2,
    regions2replace_with_model = ["rgi19"],
    mission_replace_with_model = "hugonnet",
    showplots=false,
    force_remake=false,
)

    # compute regional volume change, firn correction and mass change
    params = NamedTuple[]
    for binned_folder in binned_folders
        for paramater_set in paramater_sets
            for surface_mask in surface_masks
                for dem_id in dem_ids
                    for binning_method in binning_methods
                        for curvature_correct in curvature_corrects
                            for amplitude_correct = amplitude_corrects
                                push!(params, (; project_id, binned_folder, paramater_set, surface_mask, dem_id, binning_method, curvature_correct, amplitude_correct))
                            end
                        end
                    end
                end
            end
        end
    end

    for param in params
    #param = first(params)

        # paths to files
        binned_filled_file, figure_suffix = Altim.binned_filled_filepath(; param...)
        binned_aligned_file = replace(binned_filled_file, ".jld2" => "_aligned.jld2")

        if isfile(binned_filled_file) && ((!isfile(binned_aligned_file)) || force_remake)
            t1 = time()

            dh = FileIO.load(binned_filled_file, "dh_hyps")
            nobs_hyps = FileIO.load(binned_filled_file, "nobs_hyps")
            model_param = FileIO.load(binned_filled_file, "model_param")

            dh = geotile_align!(
                dh;
                mission_ref1,
                mission_ref2,
                min_trend_count,
                remove_land_surface_trend,
                showplots,
            )

            # filter geotiles to those that need some replacing with model 
            geotiles = Altim.geotiles_w_mask(geotile_width)
            geotiles, _ = Altim.geotiles_mutually_exclusive_rgi!(copy(geotiles))
            keep = falses(nrow(geotiles))

            for rgi in regions2replace_with_model
                keep = keep .| (geotiles[:,rgi].> 0)
            end
            geotiles2replace = geotiles[keep, :]

            dh = Altim.replace_with_model!(dh, nobs_hyps, geotiles2replace.id; mission_replace=mission_replace_with_model, mission_ref1, mission_ref2)

            save(binned_aligned_file, Dict("dh_hyps" =>
                    dh, "nobs_hyps" => nobs_hyps, "model_param" => model_param, "offset" => offset, "trend" => trend))
            
            println("$binned_aligned_file aligned and replaced: $(round(Int,time() -t1))s")
        end
    end
end


function replace_with_model!(dh, nobs, geotiles2replace::AbstractArray; mission_replace="hugonnet", mission_ref1="icesat2", mission_ref2="icesat")

    dgeotile = dims(dh[first(keys(dh))], :geotile)
    
    # replace data with model fit
    t = Altim.decimalyear.(dims(dh[first(keys(dh))], :date))
    t = repeat(t, 1, length(dims(dh[first(keys(dh))], :height)))

    h = val(dims(dh[first(keys(dh))], :height))'
    h = repeat(h, length(dims(dh[first(keys(dh))], :date)), 1)

    geotiles2replace = intersect(geotiles2replace, dgeotile)

    Threads.@threads for geotile in geotiles2replace
        #geotile = first(geotiles2replace)

        valid0 = .!isnan.(dh[mission_ref1][At(geotile), :, :])
        if .!(any(valid0))
            continue
        end

        _, vheight = Altim.validrange(valid0)
        vdates, _  = Altim.validrange(.!isnan.(dh[mission_replace][At(geotile), :, vheight]))

        dh0 = dh[mission_ref2][At(geotile), :, vheight]
        nobs0 = nobs[mission_ref2][At(geotile), :, vheight]
        dh_ref1 = dh[mission_ref1][At(geotile), :, vheight]
        nobs_ref1 = nobs[mission_ref1][At(geotile), :, vheight]
        valid_ref1 = .!isnan.(dh_ref1)
        dh0[valid_ref1] = dh_ref1[valid_ref1]
        nobs0[valid_ref1] = nobs_ref1[valid_ref1]
        valid0 = .!isnan.(dh0)

        t0 = t[:, vheight] 
        t0 = t0 .- mean(t0)
        h0 = h[:, vheight] 
        h0 = h0 .- mean(h0)

        # some tiles have nobs == 0 where there is not icesat or iceast-2 data and the 
        # geotile has been filled with with neighbor values.. for this reason we need to add
        # one to all nobs
        nobs0 = nobs0 .+ 1;

        # fit global model 
        if length(vheight) == 1
            fit1 = curve_fit(Altim.model3, t0[valid0], dh0[valid0], nobs0[valid0], Altim.p3)
            dh[mission_replace][At(geotile), vdates, vheight] = Altim.model3(t0[vdates, :][:], fit1.param)
        else
            fit1 = curve_fit(Altim.model1, hcat(t0[valid0], h0[valid0]), dh0[valid0], nobs0[valid0], Altim.p1; lower=Altim.lb1, upper=Altim.ub1)
            dh[mission_replace][At(geotile), vdates, vheight] = Altim.model1(hcat(t0[vdates, :][:], h0[vdates, :][:]), fit1.param)
        end

        #if all(fit1.param.==0)
        #    error("model fit to reference mission returned all zero parameter coefficients for geotile: $(geotile)")
        #end
       
        nobs[mission_replace][At(geotile), vdates, vheight] .= 9999;
    end
    return dh
end

function geotile2glacier!(glaciers, var0; varname)

    # initialize fields and group    
    glaciers[!, varname] = [fill(NaN, dims(var0, :date)) for i in 1:nrow(glaciers)]

    glaciers_grp = DataFrames.groupby(glaciers, :geotile)
    geotiles = getindex.(keys(glaciers_grp), "geotile")

    # doing as DimensionalData adds NO increase in processing time
    Threads.@threads for geotile in geotiles
        geotile_glaciers = glaciers_grp[(geotile,)]

        v0 = var0[At(geotile), :, :]
        rrange, _ = Altim.validrange(.!isnan.(v0))
        for glacier0 in eachrow(geotile_glaciers)
        #glacier = first(eachrow(glaciers))
            
            area_index = glacier0.area_km2 .> 0
            if !any(area_index)
                continue
            end
            crange, = Altim.validrange(area_index)
            total_area = sum(glacier0.area_km2[crange])

            out = fill(NaN, length(rrange))
            for (i, v) in enumerate(eachrow(v0[rrange, crange]))
                out[i] = sum(v .* glacier0.area_km2[crange] ./ total_area)
            end
            glacier0[varname][rrange] = out
        end
    end
    return glaciers
end

