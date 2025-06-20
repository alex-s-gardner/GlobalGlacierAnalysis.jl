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
    model1_nmad_max = 5, # = 5  this is a sigma-equivelent threshold
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
    model1_nmad_max=10, # this is a sigma-equivelent threshold
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
        model1_nmad_max=5, # = 5  this is a sigma-equivelent threshold
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
        model1_nmad_max=10, # this is a sigma-equivelent threshold
    )
end

# define model for curvature correction
model_curvature::Function = model_curvature(c, p) = p[1] .+ p[2] .* c .+ p[3] .* c .^ 2 .+ p[4] .* c .^ 3 .+ p[5] .* c .^ 4
p_curvature = zeros(5)

"""
    geotile_binning(; project_id=:v01, geotile_width=2, kwargs...)

Process satellite altimetry data into geotiles by elevation and time.

# Arguments
- `project_id`: Project identifier (default: :v01)
- `geotile_width`: Geotile width in degrees (default: 2)
- `warnings`: Show warnings (default: false)
- `plots_show`: Show plots (default: false)
- `force_remake_before`: Date to force file regeneration (default: nothing)
- `update_geotile`: Update existing geotiles (default: false)
- `update_geotile_missions`: Missions to update (default: ["icesat2"])
- `all_permutations_for_glacier_only`: Process all parameter combinations for glacier mask (default: true)
- `surface_masks`: Surface types to process (default: [:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km])
- `binned_folders`: Output directories (default: ("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"))
- `dem_ids`: DEM sources (default: [:best, :cop30_v2])
- `binning_methods`: Binning methods (default: ["nmad3", "nmad5", "median", "nmad10"])
- `curvature_corrects`: Apply curvature correction (default: [true, false])
- `max_canopy_height`: Max canopy height in meters (default: 1)
- `dh_max`: Max height difference in meters (default: 200)

# Description
Bins satellite altimetry data into geotiles by applying filters and corrections. Processes multiple missions and parameter combinations, organizing results by elevation and time.
"""
function geotile_binning(; 
    project_id = :v01,
    geotile_width = 2,
    warnings = false,
    plots_show = false,
    process_time_show = false,

    # run parameters
    force_remake_before = nothing,
    update_geotile = false, # this will load in prevous results to update select missions
    geotiles2update = nothing,
    update_geotile_missions = ["icesat2"],

    # run parameters
    all_permutations_for_glacier_only = true,
    surface_masks = [:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km],
    binned_folders=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"),
    dem_ids = [:best, :cop30_v2],
    binning_methods = ["nmad5", "nmad3", "median", "mad10"],
    curvature_corrects = [true, false],

     #### DON NOT CHANGE THESE PARAMETERS
    max_canopy_height = 1, # do not change

    # filter parameters
    dh_max=200,

    single_geotile_test=nothing,
    )

    paths = project_paths(; project_id)
    products = project_products(; project_id)

    # load geotile definitions with corresponding hypsometry
    geotiles = geotiles_w_mask(geotile_width)

    # filter geotiles
    if isnothing(single_geotile_test)
        geotiles = geotiles[(geotiles.glacier_frac.>0.0), :]
    else
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)], OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!"
        ind = findfirst(geotiles.id .== single_geotile_test)
        geotiles = geotiles[ind:ind, :]
        plots_show = true
    end

    if !isnothing(geotiles2update)
        index = [findfirst(geotiles.id .== gt) for gt in geotiles2update]
        geotiles = geotiles[index, :]
    end

    # define date bin ranges... this should match what's in gemb_classes_binning.jl
    # define date and hight binning ranges 
    date_range, date_center = project_date_bins()
    height_range, height_center = project_height_bins()

    # curvature ranges 
    Δc = 0.1;
    curvature_range = -1.0:Δc:1.0;
    curvature_center = curvature_range[1:end-1] .+ Δc/2;

    ## THIS WILL REMOVE ALL _dh_ files if you need to rebuild binned archive
    #files2delete = GGA.allfiles("/mnt/bylot-r3/data/binned/2deg"; fn_contains="_dh_")
    #rm.(files2delete)
    ##########

    params = NamedTuple[]
    for binned_folder in binned_folders
        for surface_mask in surface_masks
            for dem_id in dem_ids
                for binning_method in binning_methods
                    for curvature_correct in curvature_corrects
                        push!(params, (; project_id, binned_folder, surface_mask, dem_id, binning_method, curvature_correct))
                    end
                end
            end
        end
    end


    # Threads is throwing errors due to reading of JLD2 files, Threads is implimented at
    # lower level with reasonable performance

    # perfomance could be improved considerably if data was saved per geotile, this would allow 
    # the raw data only be read in once per geotile... this would be complicated to implement
    # as it would require all i
    @showprogress desc = "Binning hypsometric elevation change data ..." for param in params

        fig_folder = pathlocal[:figures]

        # skip permutations if all_permutations_for_glacier_only = true
        if all_permutations_for_glacier_only
            if ((!(param.surface_mask == :glacier) && Base.contains(param.binned_folder, "unfiltered")) &&
                ((param.dem_id != :best) && (param.binning_method != "mad3") && param.curvature_correct))
                continue
            end
        end

        binned_file = binned_filepath(; param.binned_folder, param.surface_mask, param.dem_id, param.binning_method, project_id, param.curvature_correct)

        if !isfile(binned_file) || !isnothing(force_remake_before) || !isnothing(single_geotile_test)

            println("binning:$binned_file")
            t_start = time()

            if isfile(binned_file) && !isnothing(force_remake_before) && isnothing(single_geotile_test)
                if Dates.unix2datetime(mtime(binned_file)) > force_remake_before
                    @warn "Skipping $(binned_file) because it was created after $force_remake_before"
                    continue
                end
            end

            # 6.2 hours for all glaciers, all missions/datasets on 96 threads
            # 10hr for land for all glaciers, all missions/datasets on 96 threads

            # initialize dimensional arrays
            # update_geotile = true
            if isfile(binned_file) && update_geotile
                # load exisiting
                dh_hyps = FileIO.load(binned_file, "dh_hyps")
                nobs_hyps = FileIO.load(binned_file, "nobs_hyps")
                curvature_hyps = FileIO.load(binned_file, "curvature_hyps")
                missions = String.(update_geotile_missions);
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
                #binned_file = "/mnt/devon2-r1/devon0/gardnera/Documents/GitHub/GlobalGlacierAnalysis/data/geotiles_glacier_hyps_2deg_dh.arrow";
                #geotiles  = DataFrame(Arrow.Table(binned_file));
                #geotiles.extent = Extent.(getindex.(geotiles.extent, 1));
                #geotiles = geotiles[(geotiles.rgi2 .> 0) .& (geotiles.glacier_frac .> 0), :] ;
                # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
                
                product = products[Symbol(mission)];
                mission_geotile_folder = paths[Symbol(mission)].geotile

                # special case for unfiltered hugonnet data
                # this is a hack for including 2 folders that contain hugonnet data (filtered and unfiltered)
                if Base.contains(param.binned_folder, "unfiltered") && (mission == "hugonnet")
                    mission_geotile_folder = replace(mission_geotile_folder, "/2deg" => "/2deg_unfiltered")
                    mission_suffix = "-unfiltered"
                else
                    mission_suffix = ""
                end
                
                #idx = (geotiles.rgi9 .> 0) .& (geotiles.glacier_frac .> 0.1)
                Threads.@threads for geotile in eachrow(geotiles)
                #for geotile in eachrow(geotiles)
                    t0 = time()

                    if geotile.glacier_frac == 0.0
                        continue
                    end

                    path2altim = joinpath(mission_geotile_folder, geotile.id * ".arrow");

                    if !isfile(path2altim)
                        continue
                    end

                    altim = select!(DataFrame(Arrow.Table(path2altim)), [:longitude, :latitude, :datetime, :height, :quality]);

                    if process_time_show 
                        t1 = time()
                        dt = round(Int16, t1 - t0)
                        println("loaded: $(mission_proper_name(mission)) - $(geotile.id): $(dt)s")
                    end

                    # add height_ref and curvature
                    altim = add_dem_ref!(altim, param.dem_id, geotile, mission_geotile_folder)

                    if process_time_show 
                        t2 = time()
                        dt = round(Int16, t2 - t1);
                        println("dem ref: $(mission_proper_name(mission)) - $(geotile.id): $(dt)s")
                    end

                    # load masks
                    path2masks = joinpath(mission_geotile_folder, geotile.id * ".masks");
                    masks0 = DataFrame(Arrow.Table(path2masks));

                    if process_time_show 
                        t3 = time()
                        dt = round(Int16, t3 - t2);
                        println("masks: $(mission_proper_name(mission)) - $(geotile.id): $(dt)s")
                    end

                    if nrow(masks0) != nrow(altim)
                        error("masks and altim do not have the same number of rows, try deleting masks and rerunning geotile_mask_extract(): $(path2masks)")
                    end
                    
                    # add high resolution mask
                    valid = within.(Ref(geotile.extent), altim.longitude, altim.latitude) .& .!isnan.(altim.dh)

                    if process_time_show
                        t4 = time()
                        dt = round(Int16, t4 - t3);
                        println("highres mask: $(mission_proper_name(mission)) - $(geotile.id): $(dt)s")
                    end

                    if plots_show
                        title = "raw $(mission_proper_name(mission))$(mission_suffix) minus $(param.dem_id) for $(geotile.id) over $(param.surface_mask): median = $(round(median(altim.dh[valid]), digits=1))"
                        f = GGA.plot_unbinned_height_anomalies(altim.datetime[valid], altim.dh[valid]; title)
                        display(f)
                    end

                    # canopy height
                    if minimum(geotile.extent.Y) < -54 || maximum(geotile.extent.Y) > 66
                        canopy = DataFrame(canopyh = zeros(size(altim.latitude))) 
                    else
                        path2canopy = joinpath(mission_geotile_folder, geotile.id * ".canopyh")
                        if isfile(path2canopy)
                            canopy = DataFrame(Arrow.Table(path2canopy));
                        else
                            @warn "Setting conopy height to zero as canopy height file does not exist: $(path2canopy)"
                            canopy = DataFrame(canopyh=zeros(size(altim.latitude)))
                        end
                    end;

                    if process_time_show
                        t5 = time()
                        dt = round(Int16, t5 - t4);
                        println("canopy: $(mission_proper_name(mission)) - $(geotile.id): $(dt)s")
                    end

                    if !any(masks0[:, param.surface_mask])
                        continue
                    end

                    ## Check slope relation
                    if false
                        for dem_id1 in dem_id0           
                            path2dem = joinpath(mission_geotile_folder, geotile.id * ".$(dem_id1)")
                            if isfile(path2dem)
                                dem = select!(DataFrame(Arrow.Table(path2dem)), :height, :dhdx, :dhdy, :dhddx, :dhddy)
                
                                dhdx_dhdy = slope.(dem.dhdx, dem.dhdy, Ref(dem_info(dem_id1)[1].epsg), lat = mean(geotile.extent.Y));
                                
                                var_ind = (.!masks0[:, param.surface_mask]) .& valid .& .!isnan.(altim.dh) .& (altim.height .!== 0) .& (altim.height_ref .!== 0) .& .!isnan.(altim.curvature) .& (abs.(altim.dh) .< 10) .& masks0.land .& ((canopy.canopyh .<= max_canopy_height) .| (altim.latitude .< -60))
                                
                                if sum(var_ind) <= length(p_curvature)
                                    warnings && (@warn ("$(geotile.id): slope offset not calculated, not enought off-ice points"))
                                else
                                    offset = track_offset(getindex.(dhdx_dhdy[var_ind], 1), getindex.(dhdx_dhdy[var_ind], 2), altim.dh[var_ind]);
                                end
                            end
                        end
                    end;

                    # use bin statistics to calculate dh vs. elevation for every 3 months.

                    # identify valid dh data
                    # plot(decimalyear.(altim.datetime[valid]),altim.dh[valid]; label="all")

                    if param.curvature_correct
                        var_ind = (.!masks0[:, :landice]) .& valid .& .!isnan.(altim.dh) .& (altim.height .!== 0) .& (altim.height_ref .!== 0) .& .!isnan.(altim.curvature) .& (abs.(altim.dh) .< 10) .& masks0.land .& ((canopy.canopyh .<= max_canopy_height) .| (altim.latitude .< -60))

                        if sum(var_ind) <= length(p_curvature)
                            warnings && (@warn ("$(geotile.id): no curvature corecction applied, not enought off-ice points"))
                        else
                            # check bounds: binstats will throw an error if no data is passed to median()
                            if vector_overlap(altim[var_ind, :curvature], curvature_range)
                                df = DataFrame();
                                try
                                    df = binstats(altim[var_ind, :], [:curvature], [curvature_range], :dh; col_function=[median], missing_bins=true);
                                catch e
                                    # having issues with edge cases
                                    println(" ##################### ERRORED ##########################")
                                    println(" ##################### $(mission_proper_name(mission)): $(geotile.id) ########################")
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
                            
                                        if plots_show
                                            dh_cor = model_curvature(bin_center, fit1.param)
                                            title = "curvature corrected $(mission_proper_name(mission))$(mission_suffix) minus DEM[$(param.dem_id)]: $(geotile.id)"
                                            f = plot_curvature(bin_center, df.dh_median, dh_cor, df.nrow; title)

                                            mission_fn = replace(mission_proper_name(mission), "R/W" => "R-W")
                                            fname = joinpath(fig_folder, "$(mission_fn)$(mission_suffix)_$(geotile.id)_$(param.dem_id)_curvature.png")
                                            CairoMakie.save(fname, f)
                                            display(f)
                                      
                                            title = "curvature corrected $(mission_proper_name(mission))$(mission_suffix) minus $(param.dem_id) for $(geotile.id): median = $(round(median(altim.dh[valid]), digits=1))"
                                            f = GGA.plot_unbinned_height_anomalies(altim.datetime[valid], altim.dh[valid]; title)
                                            display(f)
                       
                                        end
                                    end
                                end
                            end
                        end
                    end

                    if process_time_show
                        t6 = time()
                        dt = round(Int16, t6 - t5);
                        println("curvature: $(mission_proper_name(mission)) - $(geotile.id): $(dt)s")
                    end
                    
                    #################################### FILTER 1 ######################################
                    var_ind = masks0[:, param.surface_mask] .& valid .& .!isnan.(altim.dh) .& (altim.height .!== 0) .& (altim.height_ref .!== 0)
                    if sum(var_ind) < 100
                        continue
                    end

                    var_ind[var_ind] = (abs.(altim.dh[var_ind] .- median(altim.dh[var_ind])) .< dh_max) .& (abs.(altim.dh[var_ind]) .< (dh_max*2))

                    if param.surface_mask == :land
                        var_ind = var_ind .& (canopy.canopyh .<= max_canopy_height)
                    end
                    
                    if plots_show
                        title = "filtered [#1] corrected $(mission_proper_name(mission)) minus $(param.dem_id) for $(geotile.id) over $(param.surface_mask): median = $(round(median(altim.dh[var_ind]), digits=1))"
                        f = GGA.plot_unbinned_height_anomalies(altim.datetime[var_ind], altim.dh[var_ind]; title)
                        display(f)
                    end
                    
                    if process_time_show
                        t7 = time()
                        dt = round(Int16, t7 - t6);
                        println("filter 1: $(mission_proper_name(mission)) - $(geotile.id): $(dt)s")
                    end

                    # this is a hack for including 2 folders that contain hugonnet data (filtered and unfiltered)
                    if product.apply_quality_filter || (Base.contains(param.binned_folder, "unfiltered") && (mission == "hugonnet"))
                        var_ind = var_ind .& altim.quality

                        # for troubleshooting 
                        if plots_show
                            title = "quality > 70 $(mission_proper_name(mission)) minus $(param.dem_id) for $(geotile.id) over $(param.surface_mask): median = $(round(median(altim.dh[valid]), digits=1))"
                            f = GGA.plot_unbinned_height_anomalies(altim.datetime[valid], altim.dh[valid]; title)
                            display(f)
                        end

                    end
                    ####################################################################################

                    if !any(var_ind)
                        continue
                    end

                    decyear_range = decimalyear.(date_range)

                    altim[!, :decimalyear] = decimalyear.(altim.datetime)
                    
                    var0, nobs0 = geotile_bin2d(
                        altim[var_ind, :];
                        var2bin="dh",
                        dims_edges=("decimalyear" => decyear_range, "height_ref" => height_range),
                        binfunction=binningfun_define(param.binning_method))

                    if process_time_show
                        t8 = time()
                        dt = round(Int16, t8 - t7);
                        println("binned: $(mission_proper_name(mission)) - $(geotile.id): $(dt)s")
                    end

                    if isnothing(var0)
                        continue
                    end

                    try
                        dh_hyps[mission][At(geotile.id), :, :] = var0
                        nobs_hyps[mission][At(geotile.id), :, :] = nobs0
                    catch e
                        println((; param.binned_folder, param.surface_mask,param.binning_method,param.dem_id,param.curvature_correct, mission, geotile = geotile.id))
                        println("size var0 = $(size(var0))")
                        println("size nobs0 = $(size(nobs0))")
                        throw(e)
                    end

                    if plots_show
                        f = plot_elevation_time(dh_hyps[mission][At(geotile.id), :, :]; colorrange=(-20, 20))
                        f.content[1].title = "$(mission_proper_name(mission)) minus DEM[$(param.dem_id)] over $(param.surface_mask): $(geotile.id)"
                        Colorbar(f[1, 2], f.content[1].scene.plots[1], vertical=true, label="height anomaly [m]")
                        display(f)
                    end

                    if process_time_show
                        t9 = time()
                        dt = round(Int16, t9 - t8);
                        println("saved: $(mission_proper_name(mission)) - $(geotile.id): $(dt)s")
                    end

                    dt = round(Int16, time() - t0)
                    println("Total time: $(mission_proper_name(mission)) - $(geotile.id): $(dt)s")
                end
            end
            if isnothing(single_geotile_test)
                save(binned_file, Dict("dh_hyps" => dh_hyps, "nobs_hyps" => nobs_hyps, "curvature_hyps" => curvature_hyps));
                println("binning complete: $binned_file: $(round(Int16, (time() - t_start)/60))m")
            end
        end
    end
end


"""
    geotile_binned_fill(;
        project_id=:v01,
        geotile_width=2,
        force_remake_before=nothing,
        update_geotile=false,
        update_geotile_missions=["icesat2"],
        mission_reference_for_amplitude_normalization="icesat2",
        all_permutations_for_glacier_only=true,
        surface_masks=[:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km],
        binned_folders=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"),
        dem_ids=[:best, :cop30_v2],
        binning_methods=["nmad3", "nmad5", "median", "nmad10"],
        curvature_corrects=[true, false],
        paramater_sets=[1, 2],
        amplitude_corrects=[true, false],
        remove_land_surface_trend=GGA.mission_land_trend(),
        regions2replace_with_model=["rgi19"],
        mission_replace_with_model="hugonnet",
        missions2align2=["icesat2", "icesat"],
        plots_show=false,
        plots_save=false,
        plot_save_format=".png",
        process_time_show=false,
        geotiles2plot="lat[+28+30]lon[+082+084]",
        single_geotile_test=nothing,
    )

Process and fill gaps in binned altimetry data.

# Arguments
- `project_id`: Project identifier
- `geotile_width`: Width of geotiles in degrees
- `force_remake_before`: Skip files created after this date
- `update_geotile`: Update specific geotiles/missions instead of recomputing all
- `update_geotile_missions`: Missions to update when `update_geotile` is true
- `mission_reference_for_amplitude_normalization`: Reference mission for amplitude normalization
- `all_permutations_for_glacier_only`: Process all parameter combinations for glacier surfaces
- `surface_masks`: Surface types to process
- `binned_folders`: Input/output data folders
- `dem_ids`: Digital elevation models to use
- `binning_methods`: Data binning methods
- `curvature_corrects`: Whether to apply curvature corrections
- `paramater_sets`: Parameter sets for filling algorithms
- `amplitude_corrects`: Whether to apply amplitude corrections
- `remove_land_surface_trend`: Land surface trend to remove
- `regions2replace_with_model`: Regions to replace with model data
- `mission_replace_with_model`: Mission to replace with model data
- `missions2align2`: Missions to align to reference
- `plots_show`: Display diagnostic plots
- `plots_save`: Save diagnostic plots
- `plot_save_format`: Format for saved plots
- `process_time_show`: Display processing times
- `geotiles2plot`: Geotiles to include in plots
- `single_geotile_test`: Test single geotile processing

# Description
Fills temporal/spatial gaps, normalizes between missions, and applies amplitude corrections
to previously binned altimetry data. Can update specific geotiles/missions or process all data.
"""
function geotile_binned_fill(;
    project_id=:v01,
    geotile_width=2,
    force_remake_before=nothing,
    update_geotile=false, # this will load in prevous results to update select geotiles or missions
    update_geotile_missions=["icesat2"],
    mission_reference_for_amplitude_normalization="icesat2",
    all_permutations_for_glacier_only=true,
    surface_masks=[:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km],
    binned_folders=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"),
    dem_ids=[:best, :cop30_v2],
    binning_methods=["nmad3", "nmad5", "median", "nmad10"],
    curvature_corrects=[true, false],
    paramater_sets=[1, 2],
    amplitude_corrects=[true, false],
    remove_land_surface_trend=GGA.mission_land_trend(),
    regions2replace_with_model=["rgi19"],
    mission_replace_with_model="hugonnet",
    missions2align2=["icesat2", "icesat"],
    plots_show=false,
    plots_save=false,
    plot_save_format = ".png",
    process_time_show=false,
    geotiles2plot="lat[+28+30]lon[+082+084]",
    single_geotile_test = nothing,
)

    # do a little input checking
    if isa(binned_folders, String)
        binned_folders = (binned_folders,)
    end

    params = NamedTuple[]
    for binned_folder in binned_folders
        for surface_mask in surface_masks
            for dem_id in dem_ids
                for binning_method in binning_methods
                    for curvature_correct in curvature_corrects
                        push!(params, (; project_id, binned_folder, surface_mask, dem_id, binning_method, curvature_correct))
                    end
                end
            end
        end
    end

    if any(mission_reference_for_amplitude_normalization .== update_geotile_missions) && any(amplitude_corrects)
        error("can not update mission_reference_for_amplitude_normalization [$(mission_reference_for_amplitude_normalization)] without updating all other missions")
    end

    if .!isnothing(single_geotile_test)
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)], OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!"
        geotiles2plot = [single_geotile_test]
    end
    

    # load geotiles
    geotiles0 = Dict()
    for surface_mask in surface_masks
        geotiles0[surface_mask] = geotiles_mask_hyps(surface_mask, geotile_width)

        # make geotile rgi regions mutually exexclusive 
        geotiles0[surface_mask], _ = geotiles_mutually_exclusive_rgi!(geotiles0[surface_mask])
    end

    # usings threads here cuases the memory usage to explode, Threads is implimented at
    # lower level with reasonable performance
    @showprogress desc = "Filling hypsometric elevation change data ..." for param in params
        process_time_show ? t1 = time() : nothing
        fig_folder = pathlocal[:figures]
 
        if occursin("binned_unfiltered", param.binned_folder)
            fig_folder = joinpath(fig_folder, "binned_unfiltered")
        else
            fig_folder = joinpath(fig_folder, "binned")
        end

        !isdir(fig_folder) && mkdir(fig_folder)

        # skip permutations if all_permutations_for_glacier_only = true
        if all_permutations_for_glacier_only
            if ((!(param.surface_mask == :glacier) && Base.contains(param.binned_folder, "unfiltered")) &&
                ((param.dem_id != :best) && (param.binning_method != "nmad3") && param.curvature_correct))
                continue
            end
        end

        binned_file = binned_filepath(; param.binned_folder, param.surface_mask, param.dem_id, param.binning_method, project_id, param.curvature_correct)

        if !isfile(binned_file)
            @warn "binned_file does not exist, skipping: $binned_file"
            continue
        end

        # skip block of files if all have been created after force_remake_before
        skip = falses(length(paramater_sets)*length(amplitude_corrects))
        i = 1
        for paramater_set in paramater_sets
             for amplitude_correct = amplitude_corrects

                # paths to files
                binned_filled_file, figure_suffix = binned_filled_filepath(; param.binned_folder, param.surface_mask, param.dem_id, param.binning_method, project_id, param.curvature_correct, amplitude_correct, paramater_set)

                if isfile(binned_filled_file) && !isnothing(force_remake_before)
                    if Dates.unix2datetime(mtime(binned_file)) > force_remake_before
                        #try
                            #foo = FileIO.load(binned_filled_file)
                            @warn "Skipping $(binned_filled_file) because it was created after $force_remake_before"
                            skip[i] = true
                        #catch e
                        #    printstyled("!!! removing $(binned_filled_file) because it is corrupted, then recomupting !!!\n", color=:red)
                        #    rm(binned_filled_file)
                        #end
                    end
                end
                i += 1
             end
        end


        # load binned data that is the same for all paramater sets
        dh11 = FileIO.load(binned_file, "dh_hyps")

        process_time_show ? t2 = time() : nothing
        process_time_show && printstyled("dh_hyps loaded: $(round(Int16, t2 - t1))s\n", color=:green)

        if .!any(.!isnan.(dh11["hugonnet"]))
            println("NOT DATA: $binned_file")
            continue
        end

        nobs11 = FileIO.load(binned_file, "nobs_hyps")

        process_time_show ? t3 = time() : nothing
        process_time_show && printstyled("nobs_hyps loaded: $(round(Int16, t3 - t2))s\n", color=:green)


        #Threads.@threads for paramater_set in paramater_sets
        # getting some issues with variable scope, so using a for loop instead
        for paramater_set in paramater_sets
             for amplitude_correct = amplitude_corrects

                process_time_show ? t4 = time() : nothing

                dh1 = deepcopy(dh11)
                nobs1 = deepcopy(nobs11)

                if .!isnothing(single_geotile_test)
                    for k in keys(dh1)
                        ind = findfirst(dims(dh1[k], :geotile) .== single_geotile_test)
                        dh1[k] = dh1[k][ind:ind, :, :]
                        nobs1[k] = nobs1[k][ind:ind, :, :]
                    end

                    # remove empty missions
                    for k in keys(dh1)
                        if all(isnan.(dh1[k]))
                            delete!(dh1, k)
                            delete!(nobs1, k)
                        end
                    end

                end

                process_time_show ? t5 = time() : nothing
                process_time_show && printstyled("deepcopy of dh1 and nobs1: $(round(Int16, t5 - t4))s\n", color=:green)

                param_filling = binned_filling_parameters[paramater_set]

                # paths to files

                binned_filled_file, figure_suffix = binned_filled_filepath(; param.binned_folder, param.surface_mask, param.dem_id, param.binning_method, project_id, param.curvature_correct, amplitude_correct, paramater_set)

                if !isnothing(force_remake_before) || !isfile(binned_filled_file) || !isnothing(single_geotile_test)

                    if isfile(binned_filled_file) && !isnothing(force_remake_before) && isnothing(single_geotile_test)
                        if Dates.unix2datetime(mtime(binned_file)) > force_remake_before
                            continue
                        end
                    end

                    println("\n filling binned data: surface_mask = $(param.surface_mask); binning_method = $(param.binning_method); dem_id = $(param.dem_id); curvature_correct = $(param.curvature_correct); amplitude_correct = $(amplitude_correct)")

                    if update_geotile && isfile(binned_filled_file)
                        old_keys = setdiff(keys(dh1), update_geotile_missions)
                        for k in old_keys
                            delete!(dh1, k)
                        end
                    end

                    # align geotile dataframe with DimArrays
                    geotiles = deepcopy(geotiles0[param.surface_mask])
                    gt = collect(dims(dh1[first(keys(dh1))], :geotile))
                    gt_ind = [findfirst(geotiles.id .== g0) for g0 in gt]
                    geotiles = geotiles[gt_ind, :]
                    area_km2 = DimArray((reduce(hcat, geotiles[:, "$(param.surface_mask)_area_km2"])'), (dims(dh1[first(keys(dh1))], :geotile), dims(dh1[first(keys(dh1))], :height)))

                    process_time_show ? t6 = time() : nothing
                    process_time_show && printstyled("geotiles copied and aligned: $(round(Int16, t6 - t5))s\n", color=:green)

                    # create a data frame to store model parameters
                    params_fill = Dict()
                    for mission in keys(dh1)
                        n = length(dims(dh1[mission], :geotile))

                        params_fill[mission] = DataFrame(
                            geotile=val(dims(dh1[mission], :geotile)), 
                            nobs_raw=zeros(n), nbins_raw=zeros(n), 
                            nobs_final=zeros(n), nbins_filt1=zeros(n), 
                            param_m1=[fill(NaN, size(p1)) for i in 1:n], 
                            h0=fill(NaN, n), 
                            t0=fill(NaN, n), 
                            dh0=fill(NaN, n), 
                            bin_std=fill(NaN, n), 
                            bin_anom_std=fill(NaN, n),
                            )

                        for mission_ref in missions2align2
                            params_fill[mission][!, "offset"] = zeros(n)
                            params_fill[mission][!, "offset_$mission_ref"] = fill(NaN, n)
                            params_fill[mission][!, "offset_nmad_$mission_ref"] = fill(NaN, n)
                            params_fill[mission][!, "offset_nobs_$mission_ref"] = zeros(Int64,n)
                        end
                    end

                    process_time_show ? t7 = time() : nothing
                    process_time_show && printstyled("params_fill initialized: $(round(Int16, t7 - t6))s\n", color=:green)

                    # plot raw binned height anomalies
                    if plots_show || plots_save
                        colorbar_label = "binned height anomalies [m]"
                        plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")

                        plot_elevation_time_multimission_geotiles(
                            dh1,
                            geotiles2plot;
                            area_km2,
                            colorrange=(-20, 20),
                            colorbar_label,
                            hypsometry=true,
                            area_averaged=true,
                            plots_show,
                            plots_save,
                            plot_save_path_prefix,
                            plot_save_format,
                        )
                    end

                    # interpolate height anomalies
                    begin
                        hyps_model_fill!(dh1, nobs1, params_fill; bincount_min=param_filling.bincount_min,
                            model1_nmad_max=param_filling.model1_nmad_max, smooth_n=param_filling.smooth_n,
                            smooth_h2t_length_scale=param_filling.smooth_h2t_length_scale, show_times = false)

                        if plots_show || plots_save
                            colorbar_label = "interpolated height anomalies [m]"
                            plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")

                            plot_elevation_time_multimission_geotiles(
                                dh1,
                                geotiles2plot;
                                area_km2,
                                colorrange=(-20, 20),
                                colorbar_label,
                                hypsometry=true,
                                area_averaged=true,
                                plots_show,
                                plots_save,
                                plot_save_path_prefix,
                                plot_save_format,
                            )
                        end


                        process_time_show ? t8 = time() : nothing
                        process_time_show && printstyled("model_fill: $(round(Int16, t8 - t7))s\n", color=:green)
                    end

                    # make plot of the height to time variogram range ratio
                    variogram_range_ratio = false

                    # make plot of the height to time variogram range ratio
                    if variogram_range_ratio
                        range_ratio = hyps_model_fill!(dh1, nobs1, params_fill; bincount_min=param_filling.bincount_min,
                            model1_nmad_max=param_filling.model1_nmad_max, smooth_n=param_filling.smooth_n,
                            smooth_h2t_length_scale=param_filling.smooth_h2t_length_scale, variogram_range_ratio)

                        fontsize = 18
                        f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98), size=(700, 700), fontsize=fontsize)

                        for (i, mission) in enumerate(keys(dh1))
                            valid = .!isnan.(range_ratio[mission])
                            title = replace.(mission, "hugonnet" => "aster")

                            title = "$title  mean = $(round(Int, mean(range_ratio[mission][valid])))"
                            CairoMakie.Axis(f[i, 1:2], title=title)
                            CairoMakie.hist!(collect(range_ratio[mission][valid]); title=mission, bins=0:250:3000)
                        end

                        fname = joinpath(fig_folder, "$(figure_suffix)_variogram_range_ratio.png")
                        save(fname, f)
                        display(f)
                    end

                    # apply seasonal amplitude normalization
                    if amplitude_correct

                        # check if reference mission is present, if not, add
                        ref_added = false
                        if update_geotile && (!any(keys(dh1) .== mission_reference_for_amplitude_normalization)) && isfile(binned_filled_file)
                            ref_added = true
                            (dh_hyps, model_param) = FileIO.load(binned_filled_file, ("dh_hyps", "model_param"))

                            for k in setdiff(keys(dh_hyps), [mission_reference_for_amplitude_normalization])
                                delete!(dh_hyps, k)
                                delete!(model_param, k)
                            end

                            dh1 = merge(dh1, dh_hyps)
                            params_fill = merge(params_fill, model_param)
                        end

                        hyps_amplitude_normalize!(dh1, params_fill; mission_reference=mission_reference_for_amplitude_normalization, plots_show=false)

                        # remove mission_reference_for_amplitude_normalization if it was added.
                        if ref_added
                            delete!(dh1, mission_reference_for_amplitude_normalization)
                            delete!(params_fill, mission_reference_for_amplitude_normalization)
                        end

                        if plots_show || plots_save
                            colorbar_label = "normalized height anomalies [m]"
                            plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")
                            plot_elevation_time_multimission_geotiles(
                                dh1,
                                geotiles2plot;
                                area_km2,
                                colorrange=(-20, 20),
                                colorbar_label,
                                hypsometry=true,
                                area_averaged=true,
                                plots_show,
                                plots_save,
                                plot_save_path_prefix,
                                plot_save_format,
                            )

                        end
                    end

                    process_time_show ? t9 = time() : nothing
                    process_time_show && printstyled("amplitude_normalize: $(round(Int16, t9 - t8))s\n", color=:green)

                    # fill geotiles
                    begin
                        # hyps_fill_empty! can add mission data to geotiles that would 
                        # otherwise be empty. an example of this is lat[+60+62]lon[-142-140] 
                        # which has not GEDI data but GEDI data is added after hyps_fill_empty! 
                        # becuase at least on of its 5 closest neighbors have GEDI data
                        dh1 = hyps_fill_empty!(dh1, params_fill, geotiles; mask=param.surface_mask)

                        process_time_show ? t10 = time() : nothing
                        process_time_show && printstyled("hyps_fill_empty: $(round(Int16, t10 - t9))s\n", color=:green)

                        dh1 = hyps_fill_updown!(dh1, geotiles; mask=param.surface_mask)

                        process_time_show ? t11 = time() : nothing
                        process_time_show && printstyled("hyps_fill_updown: $(round(Int16, t11 - t10))s\n", color=:green)

                        if plots_show || plots_save
                            colorbar_label = "extrapolated height anomalies [m]"
                            plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")
                            plot_elevation_time_multimission_geotiles(
                                dh1,
                                geotiles2plot;
                                area_km2,
                                colorrange=(-20, 20),
                                colorbar_label,
                                hypsometry=true,
                                area_averaged=true,
                                plots_show,
                                plots_save,
                                plot_save_path_prefix,
                                plot_save_format,
                            )
                        end
                    end
                    

                    # correct for any erronious trends found over land
                    begin
                        dgeotile = dims(dh1[first(keys(dh1))], :geotile)
                        dheight = dims(dh1[first(keys(dh1))], :height)
                        
                        for mission in keys(dh1)
                            if !isnothing(remove_land_surface_trend) && (remove_land_surface_trend[At(mission)] != 0)
                                ddate = dims(dh1[mission], :date)
                                decyear = decimalyear.(ddate)

                                # center date arround mission center date
                                _, date_range, _ = GGA.validrange(.!isnan.(dh1[mission]))
                                mid_date = mean(decyear[date_range])
                                delta = (decyear .- mid_date) .* remove_land_surface_trend[At(mission)]

                                for geotile in dgeotile
                                    for height in dheight
                                        dh1[mission][geotile=At(geotile), height=At(height)] .-= delta
                                    end
                                end
                            end
                        end

                        if plots_show || plots_save
                            colorbar_label = "land surface trend corrected height anomalies [m]"
                            plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")

                            plot_elevation_time_multimission_geotiles(
                                dh1,
                                geotiles2plot;
                                area_km2,
                                colorrange=(-20, 20),
                                colorbar_label,
                                hypsometry=true,
                                area_averaged=true,
                                plots_show,
                                plots_save,
                                plot_save_path_prefix,
                                plot_save_format,
                            )
                        end
                    end


                    # align height nomalies to reference missions
                    begin
                        
                        dh1, params_fill = hyps_align_dh!(dh1, nobs1, params_fill, area_km2; missions2align2)

                        if plots_show || plots_save
                            colorbar_label = "adjusted height anomalies [m]"
                            plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")

                            plot_elevation_time_multimission_geotiles(
                                dh1,
                                geotiles2plot;
                                area_km2,
                                colorrange=(-20, 20),
                                colorbar_label,
                                hypsometry=true,
                                area_averaged=true,
                                plots_show,
                                plots_save,
                                plot_save_path_prefix,
                                plot_save_format,
                            )
                        end
                    end

                    # fill geotiles with model
                    begin
                        # filter geotiles to those that need some replacing with model
                        keep = falses(nrow(geotiles))

                        for rgi in regions2replace_with_model
                            keep = keep .| (geotiles[:, rgi] .> 0)
                        end
                        geotiles2replace = geotiles[keep, :]

                        dh1, nobs1 = replace_with_model!(dh1, nobs1, geotiles2replace.id; mission_replace=mission_replace_with_model, missions2align2)

                        if plots_show || plots_save
                            colorbar_label = "model-filled height anomalies [m]"
                            plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")

                            plot_elevation_time_multimission_geotiles(
                                dh1,
                                geotiles2plot;
                                area_km2,
                                colorrange=(-20, 20),
                                colorbar_label,
                                hypsometry=true,
                                area_averaged=true,
                                plots_show,
                                plots_save,
                                plot_save_path_prefix,
                                plot_save_format,
                            )
                        end
                    end

                    # update dh1 and nobs1 with new missions
                    if update_geotile && isfile(binned_filled_file)
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
                    if isnothing(single_geotile_test)
                        save(binned_filled_file, Dict("dh_hyps" => dh1, "nobs_hyps" => nobs1, "model_param" => params_fill))

                        process_time_show ? t12 = time() : nothing
                        process_time_show && printstyled("saving binned_filled_file: $(round(Int16, t12 - t11))s\n", color=:green)
                    end
                end
            end
        end
    end
end



"""
    replace_with_model!(dh, nobs, geotiles2replace::AbstractArray; mission_replace="hugonnet", missions2align2)

Replace elevation change data for specified geotiles with model-fitted values.

# Arguments
- `dh`: Dictionary of DimensionalArrays containing elevation change data by mission
- `nobs`: Dictionary of DimensionalArrays containing observation counts by mission
- `geotiles2replace`: Array of geotile IDs to replace with model-fitted values
- `mission_replace`: Mission whose data will be replaced (default: "hugonnet")
-  `missions2align2`: List of missions to align to the reference mission

# Returns
- Modified `dh` dictionary with replaced values for specified geotiles

# Description
Fits elevation change models to reference mission data and uses these models to replace
values in the target mission for specified geotiles. Uses a combination of reference
missions to ensure data coverage, with the primary reference taking precedence.
"""
function replace_with_model!(dh, nobs, geotiles2replace; mission_replace="hugonnet", missions2align2 = ["icesat2", "icesat"])

    dgeotile = dims(dh[first(keys(dh))], :geotile)
    dheight = dims(dh[first(keys(dh))], :height)
    ddate = dims(dh[first(keys(dh))], :date)
    dmissions = Dim{:mission}(missions2align2)
    
    geotiles2replace = intersect(geotiles2replace, dgeotile)

    Threads.@threads for geotile in geotiles2replace

        valid0 = falses(dmissions, ddate, dheight)

        # in order for the data to be replaced with model, all missions2align2 must have data
        for mission in missions2align2
            valid1 = .!isnan.(dh[mission][geotile = At(geotile)])
            if !any(valid1)
                @warn "No data for $(mission_proper_name(mission)) for geotile: $(geotile), observations not replaced with model"
                return dh
            end
            valid0[At(mission), :, :] = valid1;
        end

        # first missmissions2align2 is the perfered mission and will overwrite any overlaping data
        _, _, vheight = validrange(valid0)
        vdates, _,  = validrange(.!isnan.(dh[mission_replace][At(geotile), :, vheight]))
        
        dh0 = fill(NaN, ddate, dheight)
        nobs0 = fill(NaN, ddate, dheight)
        for mission in missions2align2 # first mission is the perfered mission
            not_valid = isnan.(dh0)
            dh0[not_valid] = dh[mission][At(geotile), :, :][not_valid]
            nobs0[not_valid] = nobs[mission][At(geotile), :, :][not_valid]
        end
        valid0 = .!isnan.(dh0)

        # replace data with model fit
        t0 = decimalyear.(ddate)
        t0 = repeat(t0, 1, length(dheight))

        h0 = val(dheight)'
        h0 = repeat(h0, length(ddate), 1)

        t0 .-= mean(decimalyear.(parent(ddate)[vdates]))
        h0 .-= mean(val(dheight[vheight]))

        # some tiles have nobs == 0 where there is not icesat or iceast-2 data and the 
        # geotile has been filled with with neighbor values.. for this reason we need to add
        # one to all nobs
        nobs0 .+= 1;

        # fit global model 
        if length(vheight) == 1
            fit1 = curve_fit(model3, t0[valid0], dh0[valid0], nobs0[valid0], p3)
            dh[mission_replace][At(geotile), vdates, vheight] = model3(vec(t0[vdates, vheight]), fit1.param)
        else
            fit1 = curve_fit(model1, hcat(t0[valid0], h0[valid0]), dh0[valid0], nobs0[valid0], p1; lower=lb1, upper=ub1)
            dh[mission_replace][At(geotile), vdates, vheight] = model1(hcat(vec(t0[vdates, vheight]), vec(h0[vdates, vheight])), fit1.param)
        end

        nobs[mission_replace][At(geotile), vdates, vheight] .= 9999;
    end
    return dh, nobs
end


"""
    geotile2glacier!(glaciers, var0; varname)

Aggregates geotile data to glacier-level time series using area-weighted averaging.

# Arguments
- `glaciers`: DataFrame with glacier metadata including geotile assignments and area_km2
- `var0`: DimArray of geotile data with dimensions (geotile, date, height)
- `varname`: Column name for storing the resulting glacier time series

# Returns
Modified `glaciers` DataFrame with new `varname` column containing glacier-level time series
"""
function geotile2glacier!(glaciers, var0; varname)
    # Pre-allocate output arrays
    ddate = dims(var0, :date)
    glaciers[!, varname] = [fill(NaN, ddate) for _ in 1:nrow(glaciers)]
    
    # Group glaciers by geotile for faster lookup
    glaciers_by_geotile = DataFrames.groupby(glaciers, :geotile)
    geotiles = getindex.(keys(glaciers_by_geotile), "geotile")
    
    # Process each geotile in parallel
    @showprogress dt=1 desc = "Geotile \"variable\" to glacier..." Threads.@threads for geotile in geotiles
        # Get glaciers in this geotile
        geotile_glaciers = glaciers_by_geotile[(geotile,)]
        
        # Get valid data range for this geotile
        v0 = var0[At(geotile), :, :]
        valid_rows, _ = validrange(.!isnan.(v0))
        
        # Skip if no valid data
        isempty(valid_rows) && continue
        
        # Pre-allocate temporary array for this geotile's calculations
        out = Vector{Float64}(undef, length(valid_rows))
        
        # Process each glacier in this geotile
        for glacier in eachrow(geotile_glaciers)
            # Get valid area indices
            area_index = glacier.area_km2 .> 0
            any(area_index) || continue
            
            # Get valid column range and total area
            valid_cols, = validrange(area_index)
            total_area = sum(view(glacier.area_km2, valid_cols))
            
            # Calculate area weights once
            area_weights = view(glacier.area_km2, valid_cols) ./ total_area
            
            # Calculate weighted averages for each row
            @views for (i, row) in enumerate(eachrow(v0[valid_rows, valid_cols]))
                out[i] = sum(row .* area_weights)
            end
            
            # Assign results to glacier
            glacier[varname][valid_rows] = out
        end
    end
    
    return glaciers
end