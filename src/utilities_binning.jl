## DO NOT CHANGE THESE
begin
binned_filling_parameters = Dict()
binned_filling_parameters[1] = (
    bincount_min=Dict("icesat" => 9,
        "icesat2" => 9,
        "gedi" => 9,
        "hugonnet" => 51), 
    smooth_n=Dict("icesat" => 5, # smooth n controls the number of points used to smooth the data
        "icesat2" => 5,
        "gedi" => 5,
        "hugonnet" => 21), 
    smooth_h2t_length_scale= 800, # 800 m = 1 year in distance for anomaly from variogram analysis =    <----------------------
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
model_curvature(c, p) = p[1] .+ p[2] .* c .+ p[3] .* c .^ 2 .+ p[4] .* c .^ 3 .+ p[5] .* c .^ 4
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
- `missions2update`: Missions to update (default: ["icesat2"])
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

    # run parameters
    force_remake_before = nothing,
    geotiles2update=nothing,# this will load in prevous results to update select missions
    missions2update = ["icesat2"],

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
    #files2delete = allfiles("/mnt/bylot-r3/data/binned/2deg"; fn_contains="_dh_")
    #rm.(files2delete)
    ##########

    param_nt = (; project_id = [project_id], surface_mask=surface_masks, dem_id=dem_ids, binning_method=binning_methods, curvature_correct=curvature_corrects, binned_folder=binned_folders)
    params = ntpermutations(param_nt)

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
                    printstyled("\n    -> Skipping $(binned_file) as it was created after the force_remake_before date: $force_remake_before \n"; color=:light_green)
                    continue
                end
            end

            # 6.2 hours for all glaciers, all missions/datasets on 96 threads
            # 10hr for land for all glaciers, all missions/datasets on 96 threads

            # initialize dimensional arrays
            if isfile(binned_file) && !(isnothing(missions2update))
                # load exisiting
                dh_hyps = FileIO.load(binned_file, "dh_hyps")
                nobs_hyps = FileIO.load(binned_file, "nobs_hyps")
                curvature_hyps = FileIO.load(binned_file, "curvature_hyps")
                missions = String.(missions2update);
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

                    # add height_ref and curvature
                    altim = add_dem_ref!(altim, param.dem_id, geotile, mission_geotile_folder)

                    # load masks
                    path2masks = joinpath(mission_geotile_folder, geotile.id * ".masks");
                    masks0 = DataFrame(Arrow.Table(path2masks));

                    if nrow(masks0) != nrow(altim)
                        error("masks and altim do not have the same number of rows, try deleting masks and rerunning geotile_mask_extract(): $(path2masks)")
                    end
                    
                    # add high resolution mask
                    valid = within.(Ref(geotile.extent), altim.longitude, altim.latitude) .& .!isnan.(altim.dh)

                    if plots_show
                        title = "raw $(mission_proper_name(mission))$(mission_suffix) minus $(param.dem_id) for $(geotile.id) over $(param.surface_mask): median = $(round(median(altim.dh[valid]), digits=1))"
                        f = plot_unbinned_height_anomalies(altim.datetime[valid], altim.dh[valid]; title)
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
                                            f = plot_unbinned_height_anomalies(altim.datetime[valid], altim.dh[valid]; title)
                                            display(f)
                       
                                        end
                                    end
                                end
                            end
                        end
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
                        f = plot_unbinned_height_anomalies(altim.datetime[var_ind], altim.dh[var_ind]; title)
                        display(f)
                    end

                    # this is a hack for including 2 folders that contain hugonnet data (filtered and unfiltered)
                    if product.apply_quality_filter || (Base.contains(param.binned_folder, "unfiltered") && (mission == "hugonnet"))
                        var_ind = var_ind .& altim.quality

                        # for troubleshooting 
                        if plots_show
                            title = "quality > 70 $(mission_proper_name(mission)) minus $(param.dem_id) for $(geotile.id) over $(param.surface_mask): median = $(round(median(altim.dh[valid]), digits=1))"
                            f = plot_unbinned_height_anomalies(altim.datetime[valid], altim.dh[valid]; title)
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
        missions2update=["icesat2"],
        mission_reference_for_amplitude_normalization="icesat2",
        all_permutations_for_glacier_only=true,
        surface_masks=[:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km],
        binned_folders=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"),
        dem_ids=[:best, :cop30_v2],
        binning_methods=["nmad3", "nmad5", "median", "nmad10"],
        curvature_corrects=[true, false],
        fill_params=[1, 2],
        amplitude_corrects=[true, false],
        remove_land_surface_trend=mission_land_trend(),
        regions2replace_with_model=["rgi19"],
        missions2replace_with_model="hugonnet",
        missions2align2=["icesat2", "icesat"],
        plots_show=false,
        plots_save=false,
        plot_save_format=".png",
        geotiles2plot="lat[+28+30]lon[+082+084]",
        single_geotile_test=nothing,
    )

Process and fill gaps in binned altimetry data.

# Arguments
- `project_id`: Project identifier
- `geotile_width`: Width of geotiles in degrees
- `force_remake_before`: Skip files created after this date
- `missions2update`: Missions to update 
- `mission_reference_for_amplitude_normalization`: Reference mission for amplitude normalization
- `all_permutations_for_glacier_only`: Process all parameter combinations for glacier surfaces
- `surface_masks`: Surface types to process
- `binned_folders`: Input/output data folders
- `dem_ids`: Digital elevation models to use
- `binning_methods`: Data binning methods
- `curvature_corrects`: Whether to apply curvature corrections
- `fill_params`: Parameter sets for filling algorithms
- `amplitude_corrects`: Whether to apply amplitude corrections
- `remove_land_surface_trend`: Land surface trend to remove
- `regions2replace_with_model`: Regions to replace with model data
- `missions2replace_with_model`: Mission to replace with model data
- `missions2align2`: Missions to align to reference
- `plots_show`: Display diagnostic plots
- `plots_save`: Save diagnostic plots
- `plot_save_format`: Format for saved plots
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
    missions2update=nothing, #["icesat2"],
    mission_reference_for_amplitude_normalization="icesat2",
    all_permutations_for_glacier_only=true,
    surface_masks=[:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km],
    binned_folders=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"),
    dem_ids=[:best, :cop30_v2],
    binning_methods=["nmad3", "nmad5", "median", "nmad10"],
    curvature_corrects=[true, false],
    fill_params=[1, 2],
    amplitude_corrects=[true, false],
    remove_land_surface_trend=mission_land_trend(),
    regions2replace_with_model=["rgi19"],
    missions2replace_with_model="hugonnet",
    missions2align2=["icesat2", "icesat"],
    plots_show=false,
    plots_save=false,
    plot_save_format = ".png",
    geotiles2plot=["lat[+28+30]lon[+082+084]"],
    single_geotile_test = nothing,
    subsample_fraction = nothing, # used to investigate solution sensitivity to subsampling (e.g. 0.1 will subsample 10% of the data)
)

    # do a little input checking
    if isa(binned_folders, String)
        binned_folders = (binned_folders,)
    end

    param_nt = (; project_id=[project_id], surface_mask=surface_masks, dem_id=dem_ids, binning_method=binning_methods, curvature_correct=curvature_corrects, binned_folder=binned_folders)
    params = ntpermutations(param_nt)


    if any(mission_reference_for_amplitude_normalization .== missions2update) && any(amplitude_corrects)
        error("can not update mission_reference_for_amplitude_normalization [$(mission_reference_for_amplitude_normalization)] without updating all other missions")
    end

    if .!isnothing(subsample_fraction) && !isnothing(missions2update)
        error("subsampling should not be used with missions2update: set `missions2update = nothing`")
    end

    if .!isnothing(subsample_fraction) && isnothing(single_geotile_test)
        @warn "!!!!!!!!!!!!!! SUBSAMLING ALL GEOTILES, OUTPUT WILL NOT BE SAVED TO FILE !!!!!!!!!!!!!!"
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

        fig_folder = joinpath(pathlocal[:figures], splitpath(param.binned_folder)[end-1])
        figure_suffix = replace(splitpath(binned_file)[end], ".jld2" => "")

        
        # load binned data that is the same for all paramater sets
        dh11 = FileIO.load(binned_file, "dh_hyps")

        if .!any(.!isnan.(dh11["hugonnet"]))
            println("!!!!!! NO HUGONNET DATA - skipping: $binned_file !!!!!!!")
            continue
        end

        nobs11 = FileIO.load(binned_file, "nobs_hyps")

        #Threads.@threads for fill_param in fill_params
        # getting some issues with variable scope, so using a for loop instead
        for fill_param in fill_params
             for amplitude_correct = amplitude_corrects

                # paths to files
                binned_filled_file = binned_filled_filepath(; param.binned_folder, param.surface_mask, param.dem_id, param.binning_method, project_id, param.curvature_correct, amplitude_correct, fill_param)

                if isfile(binned_filled_file) && !isnothing(force_remake_before) && isnothing(single_geotile_test)
                    if Dates.unix2datetime(mtime(binned_filled_file)) > force_remake_before
                        printstyled("\n    -> Skipping $(binned_filled_file) as it was created after the force_remake_before date: $force_remake_before \n"; color=:light_green)
                        continue
                    end
                elseif isfile(binned_filled_file) && isnothing(single_geotile_test)
                    printstyled("\n    -> Skipping $(binned_filled_file) as it already exists \n"; color=:light_green)
                    continue
                end

                if !isnothing(missions2update)
                    printstyled("\n   -> Filling and aligning binned data for select missions $(missions2update): $(binned_filled_file)"; color=:light_gray)
                else
                    printstyled("\n   -> Filling and aligning binned data for all missions: $(binned_filled_file)"; color=:light_gray)
                end

                dh1 = deepcopy(dh11)
                nobs1 = deepcopy(nobs11)

                param_filling = binned_filling_parameters[fill_param]
                
                if !isnothing(missions2update) && isfile(binned_filled_file)
                    missions2update = missions2update
                else
                    missions2update = keys(dh1)
                end

                # align geotile dataframe with DimArrays
                geotiles = deepcopy(geotiles0[param.surface_mask])
                
                area_km2 = _geotile_area_km2(;surface_mask=param.surface_mask, geotile_width)
                geotile_extent = _geotile_extent(;surface_mask=param.surface_mask, geotile_width)

                # filter geotiles to those that need some replacing with model
                keep = falses(nrow(geotiles))

                for rgi in regions2replace_with_model
                    keep = keep .| (geotiles[:, rgi] .> 0)
                end
                geotiles2replace = geotiles.id[keep]

                # create a data frame to store model parameters
                params_fill = Dict()
                for mission in keys(dh1)
                    n = length(dims(dh1[mission], :geotile))

                    params_fill[mission] = DataFrame(
                        geotile=collect(dims(dh1[mission], :geotile)), 
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

                # replace non-updated missions with previous results
                if !isnothing(missions2update)
                    (dh_hyps, nobs_hyps, model_param) = FileIO.load(binned_filled_file, ("dh_hyps", "nobs_hyps", "model_param"))

                    for k in setdiff(keys(dh1), missions2update)
                        dh1[k] = dh_hyps[k]
                        nobs1[k] = nobs_hyps[k]
                        params_fill[k] = model_param[k]
                    end
                end

                if .!isnothing(single_geotile_test)
                    for k in keys(dh1)
                        ind = findfirst(dims(dh1[k], :geotile) .== single_geotile_test)
                        dh1[k] = dh1[k][ind:ind, :, :]
                        nobs1[k] = nobs1[k][ind:ind, :, :]
                    end
                end

                
                if .!isnothing(subsample_fraction)
                    for k in keys(dh1)
                        rand_index = rand(dims(dh1[k])) .> subsample_fraction
                        dh1[k][rand_index] .= NaN
                        nobs1[k][rand_index] .= 0
                    end
                end

                # plot raw binned height anomalies
                if plots_show || plots_save
                    colorbar_label = "binned height anomalies"
                    plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")

                    println(plot_save_path_prefix)
                    plot_elevation_time_multimission_geotiles(
                        dh1;
                        geotiles2plot,
                        area_km2,
                        colorrange=(-20, 20),
                        colorbar_label,
                        hypsometry=true,
                        area_averaged=true,
                        plots_show,
                        plots_save,
                        plot_save_path_prefix,
                        plot_save_format,
                        mission_order=plot_order["missions"],
                    )
                end

                # correct for any erronious trends found over land
                begin
                    hyps_remove_land_surface_trend!(dh1; missions2update, remove_land_surface_trend)

                    if plots_show || plots_save
                        colorbar_label = "land surface trend corrected height anomalies"
                        plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")

                        plot_elevation_time_multimission_geotiles(
                            dh1;
                            geotiles2plot,
                            area_km2,
                            colorrange=(-20, 20),
                            colorbar_label,
                            hypsometry=true,
                            area_averaged=true,
                            plots_show,
                            plots_save,
                            plot_save_path_prefix,
                            plot_save_format,
                            mission_order=plot_order["missions"],
                        )
                    end
                end

                # interpolate height anomalies
                begin
                    hyps_model_fill!(dh1, nobs1, params_fill; missions2update, bincount_min=param_filling.bincount_min,
                        model1_nmad_max=param_filling.model1_nmad_max, smooth_n=param_filling.smooth_n,
                        smooth_h2t_length_scale=param_filling.smooth_h2t_length_scale, show_times=false, )

                    if plots_show || plots_save
                        colorbar_label = "interpolated height anomalies"
                        plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")

                        plot_elevation_time_multimission_geotiles(
                            dh1;
                            geotiles2plot,
                            area_km2,
                            colorrange=(-20, 20),
                            colorbar_label,
                            hypsometry=true,
                            area_averaged=true,
                            plots_show,
                            plots_save,
                            plot_save_path_prefix,
                            plot_save_format,
                            mission_order=plot_order["missions"],
                        )
                    end
                end

                # apply seasonal amplitude normalization
                if amplitude_correct

                    for mission in setdiff(missions2update, [mission_reference_for_amplitude_normalization])
                        hyps_amplitude_normalize!(dh1[mission], params_fill[mission], params_fill[mission_reference_for_amplitude_normalization])
                    end

                    if plots_show || plots_save
                        colorbar_label = "normalized height anomalies"
                        plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")
                        plot_elevation_time_multimission_geotiles(
                            dh1;
                            geotiles2plot,
                            area_km2,
                            colorrange=(-20, 20),
                            colorbar_label,
                            hypsometry=true,
                            area_averaged=true,
                            plots_show,
                            plots_save,
                            plot_save_path_prefix,
                            plot_save_format,
                            mission_order=plot_order["missions"],
                        )

                    end
                end

                # fill geotiles
                begin
                    # hyps_fill_empty! can add mission data to geotiles that would 
                    # otherwise be empty. an example of this is lat[+60+62]lon[-142-140] 
                    # which has not GEDI data but GEDI data is added after hyps_fill_empty! 
                    # becuase at least on of its 5 closest neighbors have GEDI data..  
                    # to limit the degree of extrapoaltion mission latitudinal limits are used
                    
                    # NOTE: if valid data extends beyond elevation range of surface_mask then extents of valid output data can differ.. this is not a problem
                    dh1 = hyps_fill_empty!(dh1, params_fill, geotile_extent, area_km2; missions2update)

                    # NOTE: extraploation of data is done after `amplitude_correct` as you 
                    # can get exteem values in poorly measured geotiles if the model is used 
                    # to extraplate the data (partifularly dh as a function of elevation)... 
                    # this should be investigated in future versions... possibly checking if 
                    #a model is appropiate for extrapolation 
                    dh1 = hyps_fill_updown!(dh1, area_km2; missions2update)

                    if plots_show || plots_save
                        colorbar_label = "extrapolated height anomalies"
                        plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")
                        plot_elevation_time_multimission_geotiles(
                            dh1;
                            geotiles2plot,
                            area_km2,
                            colorrange=(-20, 20),
                            colorbar_label,
                            hypsometry=true,
                            area_averaged=true,
                            plots_show,
                            plots_save,
                            plot_save_path_prefix,
                            plot_save_format,
                            mission_order=plot_order["missions"],
                        )
                    end
                end
                
                # align height nomalies to reference missions
                begin
                    dh1, params_fill = hyps_align_dh!(dh1, nobs1, params_fill, area_km2; missions2align2, missions2update)

                    if plots_show || plots_save
                        colorbar_label = "adjusted height anomalies"
                        plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")

                        plot_elevation_time_multimission_geotiles(
                            dh1;
                            geotiles2plot,
                            area_km2,
                            colorrange=(-20, 20),
                            colorbar_label,
                            hypsometry=true,
                            area_averaged=true,
                            plots_show,
                            plots_save,
                            plot_save_path_prefix,
                            plot_save_format,
                            mission_order=plot_order["missions"],
                        )
                    end
                end

                # fill geotiles with model
                begin
                    dh1, nobs1 = replace_with_model!(dh1, nobs1, geotiles2replace; missions2replace=intersect(missions2replace_with_model, missions2update), missions2align2)

                    if plots_show || plots_save
                        colorbar_label = "model-filled height anomalies"
                        plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")

                        plot_elevation_time_multimission_geotiles(
                            dh1;
                            geotiles2plot,
                            area_km2,
                            colorrange=(-20, 20),
                            colorbar_label,
                            hypsometry=true,
                            area_averaged=true,
                            plots_show,
                            plots_save,
                            plot_save_path_prefix,
                            plot_save_format,
                            mission_order=plot_order["missions"],
                        )
                    end
                end

                # save filled geotiles
                if isnothing(single_geotile_test) && isnothing(subsample_fraction)
                    save(binned_filled_file, Dict("dh_hyps" => dh1, "nobs_hyps" => nobs1, "model_param" => params_fill))
                else
                    return dh1, nobs1, params_fill
                end
            end
        end
    end
end

"""
    add_dem_ref!(altim, dem_id, geotile, mission_geotile_folder) -> DataFrame

Add digital elevation model (DEM) reference information to an altimetry DataFrame.

# Arguments
- `altim::DataFrame`: Altimetry data to which DEM information will be added
- `dem_id::Symbol`: DEM identifier (`:best`, `:cop30_v2`, `:arcticdem_v4_10m`, or `:rema_v2_10m`)
- `geotile`: Geotile object containing extent information
- `mission_geotile_folder::String`: Path to folder containing geotile DEM files

# Returns
- `DataFrame`: The modified altimetry DataFrame with added columns:
  - `height_ref`: Reference height from DEM
  - `curvature`: Surface curvature
  - `dh`: Height difference (altim.height - height_ref)

When `dem_id` is `:best`, tries multiple DEMs in order of preference.
"""
function add_dem_ref!(altim, dem_id, geotile, mission_geotile_folder)

    # add dem height and curvature
    if dem_id == :best
        # last dem takes precidence over earlier dems
        dem_id0 = [:cop30_v2, :arcticdem_v4_10m, :rema_v2_10m]
    elseif any([:cop30_v2, :arcticdem_v4_10m, :rema_v2_10m] .== dem_id)
        dem_id0 = [dem_id]
    else
        error("unreconized dem_id: $dem_id")
    end

    altim[!, :height_ref] .= NaN
    altim[!, :curvature] .= NaN

    for dem_id1 in dem_id0
        path2dem = joinpath(mission_geotile_folder, geotile.id * ".$(dem_id1)")

        if isfile(path2dem)
            dem = select!(DataFrame(Arrow.Table(path2dem)), :height, :dhddx, :dhddy)

            # a bit faster to calculate curvature on all data then subset
            curv = curvature.(dem.dhddx, dem.dhddy, Ref(dem_info(dem_id1)[1].epsg), lat=mean(geotile.extent.Y))

            ind = .!isnan.(curv) .& (abs.(dem.height) .< 9998)

            try
                altim[ind, :height_ref] = dem[ind, :height]
                altim[ind, :curvature] = curv[ind]
            catch e
                @warn "Error adding dem info to altim dataframe, it is most likely that the dem is out of data: $path2dem"
                throw(e)
            end
        end
    end

    altim[!, :dh] = altim.height .- altim.height_ref
    return altim
end


"""
    highres_mask(extent, feature, invert, excludefeature) -> (mask, area_m2)

Create a high-resolution binary mask from vector features for a geotile.

# Arguments
- `extent`: Extent object with extent information
- `feature`: Vector feature to rasterize into a mask
- `invert`: Boolean indicating whether to invert the mask
- `excludefeature`: Optional feature to exclude from the mask (can be nothing)

# Returns
- `mask`: Binary raster mask at ~30m resolution
- `area_m2`: Matrix of cell areas in square meters
"""
function highres_mask(extent, feature, invert, excludefeature)
    # update mask with high-resolution vector files
    grid_resolution = 0.00027 # ~30m

    x_mask = X(extent.X[1]:grid_resolution:extent.X[2],
        sampling=DimensionalData.Intervals(DimensionalData.Start()))
    y_mask = Y(extent.Y[1]:grid_resolution:extent.Y[2],
        sampling=DimensionalData.Intervals(DimensionalData.Start()))

    mask1 = Raster(zeros(UInt8, y_mask, x_mask))
    setcrs(mask1, EPSG(4326))

    # NOTE: count method is fastest
    mask1 = Rasters.rasterize!(count, mask1, feature; threaded=false,
        shape=:polygon, progress=false, verbose=false, boundary=:center) .> 0

    if invert
        mask1 = .!(mask1)
    end

    if !isnothing(excludefeature)
        excludemask = Raster(zeros(UInt8, y_mask, x_mask))
        excludemask = Rasters.rasterize!(count, excludemask, excludefeature; threaded=false, shape=:polygon, progress=false, verbose=false) .> 0
        mask1 = mask1 .& .!excludemask
    end

    # calculate area per cell
    lon = lookup(mask1, X)
    lat = lookup(mask1, Y)
    d = meters2lonlat_distance.(Ref(1), lat)
    a = abs.((1 ./ getindex.(d, 2) * (lat[2] .- lat[1])) .* (1 / d[1][1] * (lon[2] - lon[1])))
    area_m2 = repeat(a', outer=[length(lon), 1])

    return (mask1, area_m2)
end


"""
    highres_mask!(masks0, altim, extent, feature, invert, excludefeature, surface_mask) -> DataFrame

Apply a high-resolution binary mask to points within a geotile and update a DataFrame column.

# Arguments
- `masks0::DataFrame`: DataFrame to update with mask values
- `altim`: Object containing longitude and latitude coordinates
- `extent`: Extent object with extent information
- `feature`: Vector feature to rasterize into a mask
- `invert::Bool`: Whether to invert the mask
- `excludefeature`: Optional feature to exclude from the mask (can be nothing)
- `surface_mask::Symbol`: Column name in masks0 to update with mask values

# Returns
- Updated DataFrame with the surface_mask column populated
"""
function highres_mask!(masks0, altim, extent, feature, invert, excludefeature, surface_mask)
    mask1, _ = highres_mask(extent, feature, invert, excludefeature)

    valid = within.(Ref(extent), altim.longitude, altim.latitude)
    masks0[!, surface_mask] .= false

    fast_index = true
    if fast_index # fast index is 6x faster than Rasters [v0.14.2] with identical output
        grid_resolution = val(dims(mask1, :X).val.span) # ~30m

        x_mask = X(extent.X[1]:grid_resolution:extent.X[2],
            sampling=DimensionalData.Intervals(DimensionalData.Start()))
        y_mask = Y(extent.Y[1]:grid_resolution:extent.Y[2],
            sampling=DimensionalData.Intervals(DimensionalData.Start()))

        c = floor.(Int64, (altim.longitude[valid] .- first(x_mask)) ./ step(x_mask)) .+ 1
        r = floor.(Int64, (altim.latitude[valid] .- first(y_mask)) ./ step(y_mask)) .+ 1
        pts3 = [CartesianIndex(a, b) for (a, b) in zip(r, c)]
        masks0[valid, surface_mask] .= mask1[pts3]
    else
        #NOTE: 67% of time is taken here for large number of points.
        pts1 = GeoInterface.Point.(tuple.(altim.longitude[valid], altim.latitude[valid]))
        pts2 = Rasters.extract(mask1, pts1, atol=grid_resolution / 2, index=true, geometry=false)
        masks0[:, surface_mask][valid] = coalesce.(getindex.(pts2, 2), false)
    end

    return masks0
end


"""
    binningfun_define(binning_method) -> Function

Create a binning function based on the specified method.

# Arguments
- `binning_method::String`: Method to use for binning data. Options:
  - "nmad3": Mean of values with MAD normalization < 3
  - "nmad5": Mean of values with MAD normalization < 5
  - "nmad10": Mean of values with MAD normalization < 10
  - "median": Median of all values

# Returns
- Function that implements the specified binning method
"""
function binningfun_define(binning_method)
    if binning_method == "nmad3"
        x -> mean(x[nmad(x).<3])
    elseif binning_method == "nmad5"
        x -> mean(x[nmad(x).<5])
    elseif binning_method == "nmad10"
        x -> mean(x[nmad(x).<10])
    elseif binning_method == "median"
        x -> median(x)
    elseif binning_method == "mean"
        x -> mean(x)
    else
        error("unrecognized binning method")
    end
end
"""
    highres_mask(latitude, longitude, feature; grid_resolution=0.00027) -> Vector{Bool}

Create a binary mask for points based on their intersection with a high-resolution vector feature.

# Arguments
- `latitude::Vector{<:Real}`: Vector of latitude coordinates
- `longitude::Vector{<:Real}`: Vector of longitude coordinates
- `feature`: Vector feature (polygon) used to create the mask

# Keywords
- `grid_resolution=0.00027`: Resolution of the rasterized grid in degrees

# Returns
- Vector of boolean values indicating whether each point intersects with the feature
"""
function highres_mask(latitude, longitude, feature; grid_resolution=0.00027)
    # update mask with high-resolution vector files

    ymin = floor(minimum(latitude) ./ grid_resolution) .* grid_resolution
    ymax = ceil(maximum(latitude) ./ grid_resolution) .* grid_resolution
    xmin = floor(minimum(longitude) ./ grid_resolution) .* grid_resolution
    xmax = ceil(maximum(longitude) ./ grid_resolution) .* grid_resolution

    x_mask = X(xmin:grid_resolution:xmax, sampling=DimensionalData.Intervals(DimensionalData.Start()))
    y_mask = Y(ymin:grid_resolution:ymax, sampling=DimensionalData.Intervals(DimensionalData.Start()))

    mask1 = Raster(zeros(UInt8, y_mask, x_mask))

    # NOTE: count method is fastest
    mask1 = Rasters.rasterize!(count, mask1, feature; threaded=false, shape=:polygon, progress=false, verbose=false, boundary=:center) .> 0

    fast_index = true
    if fast_index # fast index is 15x faster than Rasters
        c = floor.(Int64, (longitude .- first(x_mask)) ./ step(x_mask)) .+ 1
        r = floor.(Int64, (latitude .- first(y_mask)) ./ step(y_mask)) .+ 1
        pts3 = [CartesianIndex(a, b) for (a, b) in zip(r, c)]
        mask = mask1[pts3]
    else
        #NOTE: 67% of time is taken here for large number of points.
        pts1 = GeoInterface.PointTuple.([((Y=y, X=x)) for (x, y) in
                                         zip(longitude, latitude)])
        pts2 = extract(mask1, pts1, atol=grid_resolution / 2, index=true, geometry=false)
        mask = getindex.(pts2, 2)
    end

    return mask
end

"""
    validgaps(valid) -> BitVector

Identify gaps within a valid data range.

# Arguments
- `valid`: BitVector or Boolean array indicating valid data points

# Returns
- BitVector marking gaps (false values) that occur between the first and last valid points

# Details
Returns a BitVector where `true` indicates invalid points (gaps) that occur between the first 
and last valid points in the input array. Points outside this range are always marked as `false`.
"""
function validgaps(valid)
    validgap = falses(length(valid))
    if any(valid) && !all(valid)
        sind = findfirst(valid)
        eind = findlast(valid)
        validgap[sind:eind] = .!valid[sind:eind]
    end
    return validgap
end
