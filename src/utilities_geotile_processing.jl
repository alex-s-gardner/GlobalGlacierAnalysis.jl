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

"""
    geotile_binning(; project_id=:v01, geotile_width=2, kwargs...)

Bin altimetry data into geotiles by elevation and time.

# Arguments
- `project_id`: Project identifier (default: :v01)
- `geotile_width`: Width of geotiles in degrees (default: 2)
- `warnings`: Whether to show warnings (default: false)
- `showplots`: Whether to show plots (default: false)
- `force_remake_before`: Only remake files created before this date (default: nothing)
- `update_geotile`: Whether to update existing geotiles (default: false)
- `update_geotile_missions`: Missions to update if updating geotiles (default: ["icesat2"])
- `all_permutations_for_glacier_only`: Whether to run all permutations only for glacier mask (default: true)
- `surface_masks`: Surface masks to use (default: [:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km])
- `binned_folders`: Folders to store binned data (default: ("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"))
- `dem_ids`: DEM identifiers to use (default: [:best, :cop30_v2])
- `binning_methods`: Methods for binning (default: ["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"])
- `curvature_corrects`: Whether to apply curvature correction (default: [true, false])
- `max_canopy_height`: Maximum canopy height to consider (default: 1)
- `dh_max`: Maximum height difference to consider (default: 200)

# Description
Bins altimetry data from various missions into geotiles, applying filtering and corrections.
The function processes data for each combination of parameters, binning by elevation and time.
"""
function geotile_binning(; 
    project_id = :v01,
    geotile_width = 2,
    warnings = false,
    showplots = false,

    # run parameters
    force_remake_before = nothing,
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
    height_range, height_center = Altim.project_height_bins()

    # curvature ranges 
    Δc = 0.1;
    curvature_range = -1.0:Δc:1.0;
    curvature_center = curvature_range[1:end-1] .+ Δc/2;

    ## THIS WILL REMOVE ALL _dh_ files if you need to rebuild binned archive
    #files2delete = allfiles("/mnt/bylot-r3/data/binned/2deg"; fn_contains="_dh_")
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
    @showprogress desc = "Binning hypsometric elevation change data ..." for param in params

        fig_folder = joinpath(param.binned_folder, "figures")
        !isdir(fig_folder) && mkdir(fig_folder)

        # skip permutations if all_permutations_for_glacier_only = true
        if all_permutations_for_glacier_only
            if ((!(param.surface_mask == :glacier) && Base.contains(param.binned_folder, "unfiltered")) &&
                ((param.dem_id != :best) && (param.binning_method != "meanmadnorm3") && param.curvature_correct))
                continue
            end
        end

        # funtion used for binning data
        binningfun = Altim.binningfun_define(param.binning_method)
    
        binned_file = Altim.binned_filepath(; param.binned_folder, param.surface_mask, param.dem_id, param.binning_method, project_id, param.curvature_correct)

        # <><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><>
        # binningfun(x) = mean(x[Altim.madnorm(x).<10])
        # <><><><><><><><><><><><><><><><><><><><><><<><><><><><><><><><><><><><><><><><>
        if !isfile(binned_file) || !isnothing(force_remake_before)

            println("binning:$binned_file")

            if isfile(binned_file) && !isnothing(force_remake_before)
                if Dates.unix2datetime(mtime(binned_file)) > force_remake_before
                    @warn "Skipping $(binned_file) because it was created after $force_remake_before"
                    continue
                end
            end

            # open shapefiles
            if param.surface_mask == :land
                shp = Symbol("$(:water)_shp")
                fn_shp = Altim.pathlocal[shp]
                feature = Shapefile.Handle(fn_shp)
                invert = true

                shp = Symbol("$(:landice)_shp")
                fn_shp = Altim.pathlocal[shp]
                excludefeature = Shapefile.Handle(fn_shp)
            else
                shp = Symbol("$(param.surface_mask)_shp")
                fn_shp = Altim.pathlocal[shp]
                feature = Shapefile.Handle(fn_shp)
                invert = false

                excludefeature = nothing
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
                    if Base.contains(param.binned_folder, "unfiltered") && (mission == "hugonnet")
                        mission_geotile_folder = replace(mission_geotile_folder, "/2deg" => "/2deg_unfiltered")
                    end

                    path2altim = joinpath(mission_geotile_folder, geotile.id * ".arrow");

                    if !isfile(path2altim)
                        continue
                    end

                    altim = select!(DataFrame(Arrow.Table(path2altim)), [:longitude, :latitude, :datetime, :height, :quality]);

                    # add height_ref and curvature
                    altim = Altim.add_dem_ref!(altim, param.dem_id, geotile, mission_geotile_folder)
                    
                    # load masks
                    path2masks = joinpath(mission_geotile_folder, geotile.id * ".masks");
                    masks0 = select!(DataFrame(Arrow.Table(path2masks)), [:inlandwater, :land, :landice, :ocean]);

                    if nrow(masks0) != nrow(altim)
                        error("masks and altim do not have the same number of rows, try deleting masks and rerunning geotile_mask_extract(): $(path2masks)")
                    end
                    
                    # add high resolution mask 
                    masks0 = Altim.highres_mask!(masks0, altim, geotile, feature, invert, excludefeature, param.surface_mask)
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
                
                                dhdx_dhdy = Altim.slope.(dem.dhdx, dem.dhdy, Ref(Altim.dem_info(dem_id1)[1].epsg), lat = mean(geotile.extent.Y));
                                
                                var_ind = (.!masks0[:, param.surface_mask]) .& valid .& .!isnan.(altim.dh) .& (altim.height .!== 0) .& (altim.height_ref .!== 0) .& .!isnan.(altim.curvature) .& (abs.(altim.dh) .< 10) .& masks0.land .& ((canopy.canopyh .<= max_canopy_height) .| (altim.latitude .< -60))
                                
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

                    if param.curvature_correct
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
                                            
                                            fname = joinpath(fig_folder, "$(param.surface_mask)_$(mission)_$(geotile.id)_$(param.dem_id)_curvature.png")
                                            save(fname, f)
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

                    # for troubleshooting
                    showplots && CairoMakie.plot!(altim.datetime[var_ind],altim.dh[var_ind]; seriestype=:scatter,label="glacier")

                    var_ind[var_ind] = (abs.(altim.dh[var_ind] .- median(altim.dh[var_ind])) .< dh_max) .& (abs.(altim.dh[var_ind]) .< (dh_max*2))

                    if param.surface_mask == :land
                        var_ind = var_ind .& (canopy.canopyh .<= max_canopy_height)
                    end
                    
                    # for troubleshooting
                    showplots && CairoMakie.plot!(altim.datetime[var_ind],altim.dh[var_ind]; seriestype=:scatter, label="glacier-filt", ylims = (-300, +300))
                    
                    # this is a hack for including 2 folders that contain hugonnet data (filtered and unfiltered)
                    if product.apply_quality_filter || (Base.contains(param.binned_folder, "unfiltered") && (mission == "hugonnet"))
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
                        binfunction=Altim.binningfun_define(param.binning_method))

                    if isnothing(var0)
                        continue
                    end

                    # for troubleshooting 
                    showplots && Plots.heatmap(var0, clim=(-10, 10))

                    try
                        dh_hyps[mission][At(geotile.id), :, :] = var0
                        nobs_hyps[mission][At(geotile.id), :, :] = nobs0
                    catch e
                        println((; param.binned_folder, param.surface_mask,param.binning_method,param.dem_id,param.curvature_correct, mission, geotile = geotile.id))
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


"""
    geotile_binned_fill(;
        project_id=:v01,
        geotile_width=2,
        force_remake_before=nothing,
        update_geotile=false,
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
        showplots=false,
        show_times=false,
    )

Processes and fills binned geotile elevation change data using various models and corrections.

# Arguments
- `project_id`: Identifier for the project version
- `geotile_width`: Width of geotiles in degrees
- `force_remake_before`: Skip files created after this datetime
- `update_geotile`: Whether to update specific geotiles rather than recomputing all
- `update_geotile_missions`: List of missions to update when `update_geotile` is true
- `plot_dh_as_function_of_time_and_elevation`: Generate diagnostic plots of elevation change
- `mission_reference_for_amplitude_normalization`: Reference mission for amplitude normalization
- `all_permutations_for_glacier_only`: Only process all parameter permutations for glacier mask
- `surface_masks`: List of surface masks to process
- `binned_folders`: Folders containing binned data
- `dem_ids`: Digital elevation models to use
- `binning_methods`: Statistical methods for binning data
- `curvature_corrects`: Whether to apply curvature correction
- `paramater_sets`: Parameter sets for model filling
- `amplitude_corrects`: Whether to apply amplitude correction
- `showplots`: Display plots during processing
- `show_times`: Print timing information

# Description
This function processes binned geotile data by applying various corrections and filling 
missing data using models. It can update specific geotiles/missions or process all data.
The function handles land adjustments, model fitting, amplitude normalization, and filling
of empty or incomplete data.
"""
function geotile_binned_fill(;
    project_id=:v01,
    geotile_width=2,
    force_remake_before=nothing,
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
    showplots=false,
    show_times=false,
)

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

    # load geotiles
    geotiles0 = Dict()
    for surface_mask in surface_masks
        geotiles0[surface_mask] = Altim.geotiles_mask_hyps(surface_mask, geotile_width)

        # make geotile rgi regions mutually exexclusive 
        geotiles0[surface_mask], _ = Altim.geotiles_mutually_exclusive_rgi!(geotiles0[surface_mask])
    end

    # usings threads here cuases the memory usage to explode, Threads is implimented at
    # lower level with reasonable performance
    @showprogress desc = "Filling hypsometric elevation change data ..." for param in params
        show_times ? t1 = time() : nothing
        fig_folder = joinpath(param.binned_folder, "figures")
        !isdir(fig_folder) && mkdir(fig_folder)

        # skip permutations if all_permutations_for_glacier_only = true
        if all_permutations_for_glacier_only
            if ((!(param.surface_mask == :glacier) && Base.contains(param.binned_folder, "unfiltered")) &&
                ((param.dem_id != :best) && (param.binning_method != "meanmadnorm3") && param.curvature_correct))
                continue
            end
        end

        binned_file = Altim.binned_filepath(; param.binned_folder, param.surface_mask, param.dem_id, param.binning_method, project_id, param.curvature_correct)
        binned_file_land = Altim.binned_filepath(; param.binned_folder, surface_mask=:land, param.dem_id, param.binning_method, project_id, param.curvature_correct)

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
                binned_filled_file, figure_suffix = Altim.binned_filled_filepath(; param.binned_folder, param.surface_mask, param.dem_id, param.binning_method, project_id, param.curvature_correct, amplitude_correct, paramater_set)

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

        if all(skip)
            continue
        end

        # load binned data that is the same for all paramater sets
        dh11 = FileIO.load(binned_file, "dh_hyps")

        show_times ? t2 = time() : nothing
        show_times && printstyled("dh_hyps loaded: $(round(Int16, t2 - t1))s\n", color=:green)

        if .!any(.!isnan.(dh11["hugonnet"]))
            println("NOT DATA: $binned_file")
            continue
        end

        nobs11 = FileIO.load(binned_file, "nobs_hyps")

        show_times ? t3 = time() : nothing
        show_times && printstyled("nobs_hyps loaded: $(round(Int16, t3 - t2))s\n", color=:green)

        Threads.@threads for paramater_set in paramater_sets
             for amplitude_correct = amplitude_corrects

                show_times ? t4 = time() : nothing
                

                dh1 = deepcopy(dh11)
                nobs1 = deepcopy(nobs11)

                show_times ? t5 = time() : nothing
                show_times && printstyled("deepcopy of dh1 and nobs1: $(round(Int16, t5 - t4))s\n", color=:green)

                param_filling = Altim.binned_filling_parameters[paramater_set]

                # paths to files

                binned_filled_file, figure_suffix = Altim.binned_filled_filepath(; param.binned_folder, param.surface_mask, param.dem_id, param.binning_method, project_id, param.curvature_correct, amplitude_correct, paramater_set)

                if !isnothing(force_remake_before) || !isfile(binned_filled_file)

                    if isfile(binned_filled_file) && !isnothing(force_remake_before)
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

                    show_times ? t6 = time() : nothing
                    show_times && printstyled("geotiles copied and aligned: $(round(Int16, t6 - t5))s\n", color=:green)

                    # create a data frame to store model parameters
                    # initialize dimensional arrays
                    params_fill = Dict()
                    for mission in keys(dh1)
                        n = length(dims(dh1[mission], :geotile))
                        params_fill[mission] = DataFrame(geotile=val(dims(dh1[mission], :geotile)), nobs_raw=zeros(n), nbins_raw=zeros(n), nobs_final=zeros(n), nbins_filt1=zeros(n), param_m1=[fill(NaN, size(Altim.p1)) for i in 1:n], h0=fill(NaN, n), t0=fill(NaN, n), dh0=fill(NaN, n), bin_std=fill(NaN, n), bin_anom_std=fill(NaN, n))
                    end

                    show_times ? t7 = time() : nothing
                    show_times && printstyled("params_fill initialized: $(round(Int16, t7 - t6))s\n", color=:green)

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
                        smooth_h2t_length_scale=param_filling.smooth_h2t_length_scale, show_times = false)

                    show_times ? t8 = time() : nothing
                    show_times && printstyled("model_fill: $(round(Int16, t8 - t7))s\n", color=:green)

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

                        Altim.hyps_amplitude_normalize!(dh1, params_fill; mission_reference=mission_reference_for_amplitude_normalization)

                        # remove mission_reference_for_amplitude_normalization if it was added.
                        if ref_added
                            delete!(dh1, mission_reference_for_amplitude_normalization)
                            delete!(params_fill, mission_reference_for_amplitude_normalization)
                        end
                    end

                    show_times ? t9 = time() : nothing
                    show_times && printstyled("amplitude_normalize: $(round(Int16, t9 - t8))s\n", color=:green)

                    if plot_dh_as_function_of_time_and_elevation
                        Altim.plot_height_time(dh1; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="filled_ampnorm", fig_folder, figure_suffix, mask=param.surface_mask, showplots)
                    end

                    dh1 = Altim.hyps_fill_empty!(dh1, params_fill, geotiles; mask=param.surface_mask)

                    show_times ? t10 = time() : nothing
                    show_times && printstyled("hyps_fill_empty: $(round(Int16, t10 - t9))s\n", color=:green)

                    dh1 = Altim.hyps_fill_updown!(dh1, geotiles; mask=param.surface_mask)

                    show_times ? t11 = time() : nothing
                    show_times && printstyled("hyps_fill_updown: $(round(Int16, t11 - t10))s\n", color=:green)

                    if plot_dh_as_function_of_time_and_elevation
                        Altim.plot_height_time(dh1; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="filled_updown", fig_folder, figure_suffix, mask=param.surface_mask, showplots)
                    end

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
                    save(binned_filled_file, Dict("dh_hyps" => dh1, "nobs_hyps" => nobs1, "model_param" => params_fill))
                    
                    show_times ? t12 = time() : nothing
                    show_times && printstyled("saving binned_filled_file: $(round(Int16, t12 - t11))s\n", color=:green)
                end
            end
        end
    end
end

"""
    geotile_filled_landfit(;
        project_id="v01",
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

Calculates and fits land surface elevation change models for geotiles.

# Arguments
- `project_id`: Identifier for the project version
- `showplots`: Whether to display plots during processing
- `binned_folders`: Folders containing binned data
- `geotile_width`: Width of geotiles in degrees
- `paramater_sets`: Parameter sets for model fitting
- `surface_masks`: Surface masks to process (typically land)
- `binning_methods`: Statistical methods for binning data
- `dem_ids`: Digital elevation models to use
- `curvature_corrects`: Whether to apply curvature correction
- `amplitude_corrects`: Whether to apply amplitude correction
- `force_remake_landoffset`: Force recalculation of land offset parameters

# Description
This function processes binned geotile data for land surfaces, calculating regional 
elevation changes and fitting temporal models to the data. It generates diagnostic plots
showing elevation trends and saves the fitted parameters for later use in glacier
elevation change corrections.
"""
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


"""
    geotile_adjust!(dh, dh_land, nobs_land, dh_land_ref; kwargs...)

Adjust elevation change estimates based on land anomalies and reference mission seasonality.

# Arguments
- `dh`: Elevation change DimArray to adjust
- `dh_land`: Land elevation change DimArray for the same mission
- `nobs_land`: Number of observations in land measurements
- `dh_land_ref`: Reference mission land elevation change DimArray (typically icesat2)
- `minnobs`: Minimum number of observations required for valid adjustment (default: 45)
- `ref_minvalid`: Minimum number of valid reference points required (default: 12)
- `ref_madnorm_max`: Maximum MAD normalization threshold for reference data (default: 10)

# Returns
- Adjusted elevation change DimArray with land anomalies and seasonality corrections applied

# Description
This function corrects elevation change estimates by removing anomalies measured over land
and applying seasonality patterns from a reference mission. It fits a seasonal cycle to the
reference data and applies corrections to each geotile independently.
"""
function geotile_adjust!(
    dh,
    dh_land,
    nobs_land,
    dh_land_ref;
    minnobs=45,
    ref_minvalid=12,
    ref_madnorm_max=10,
)

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
end

offset_trend::Function = offset_trend(t, p) = p[1] .+ p[2] .* t;
offset_trend_p = zeros(2);

"""
    geotile_align!(
        dh;
        mission_ref1="icesat2",
        mission_ref2="icesat",
        min_trend_count=5,
        remove_land_surface_trend=Altim.mission_land_trend(),
        showplots=false,
    )

Aligns elevation change data from multiple missions to reference missions.

# Arguments
- `dh`: Dictionary of elevation change DimArrays by mission
- `mission_ref1`: Primary reference mission (default: "icesat2")
- `mission_ref2`: Secondary reference mission (default: "icesat")
- `min_trend_count`: Minimum number of valid points required for trend fitting (default: 5)
- `remove_land_surface_trend`: Dictionary of land surface trends by mission to remove (default: from Altim.mission_land_trend())
- `showplots`: Whether to display diagnostic plots (default: false)

# Description
This function aligns elevation change data from multiple missions by:
1. Setting a common reference level using the primary and secondary reference missions
2. Removing offsets between missions
3. Correcting for trend differences between missions
4. Applying land surface trend corrections if specified

Returns the aligned elevation change dictionary.
"""
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

"""
    geotile_align_replace(;
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
        force_remake_before=nothing,
    )

Aligns elevation change data across missions and replaces specific regional data with model-derived values.

# Arguments
- `mission_ref1`: Primary reference mission for alignment (default: "icesat2")
- `mission_ref2`: Secondary reference mission for alignment (default: "icesat")
- `min_trend_count`: Minimum number of points required for trend calculation (default: 5)
- `remove_land_surface_trend`: Land surface trends to remove from data (default: from Altim.mission_land_trend())
- `project_id`: Project identifier (default: :v01)
- `surface_masks`: Surface masks to process (default: [:glacier, :glacier_rgi7, :glacier_b1km, :glacier_b10km])
- `binned_folders`: Folders containing binned data (default: ("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"))
- `dem_ids`: Digital elevation models to use (default: [:best, :cop30_v2])
- `binning_methods`: Statistical methods for binning data (default: ["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"])
- `curvature_corrects`: Whether to apply curvature correction (default: [true, false])
- `paramater_sets`: Parameter sets for model filling (default: [1, 2])
- `amplitude_corrects`: Whether to apply amplitude correction (default: [true, false])
- `geotile_width`: Width of geotiles in degrees (default: 2)
- `regions2replace_with_model`: RGI regions to replace with model-derived data (default: ["rgi19"])
- `mission_replace_with_model`: Mission whose data will be replaced with model (default: "hugonnet")
- `showplots`: Whether to display plots during processing (default: false)
- `force_remake_before`: Only remake files created before this date (default: nothing)

# Description
This function processes binned geotile data by first aligning elevation change measurements 
across different missions using reference missions, then replacing specific regional data 
with model-derived values. It handles trend corrections, offset adjustments, and saves the 
aligned and replaced data to new files.
"""
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
    force_remake_before=nothing,
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

    @showprogress dt = 1 desc = "Aligning geotile data with $(mission_ref1) and $(mission_ref2)..." for param in params
    #param = first(params)

        # paths to files
        binned_filled_file, figure_suffix = Altim.binned_filled_filepath(; param...)
        binned_aligned_file = replace(binned_filled_file, ".jld2" => "_aligned.jld2")

        if isfile(binned_filled_file) && ((!isfile(binned_aligned_file)) || !isnothing(force_remake_before))

            if isfile(binned_aligned_file) && !isnothing(force_remake_before)
                if Dates.unix2datetime(mtime(binned_aligned_file)) > force_remake_before
                    @warn "Skipping $(binned_aligned_file) because it was created after $force_remake_before"
                    continue
                end
            end
            
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
                    dh, "nobs_hyps" => nobs_hyps, "model_param" => model_param))

            #save(binned_aligned_file, Dict("dh_hyps" =>
            #        dh, "nobs_hyps" => nobs_hyps, "model_param" => model_param, "offset" => offset, "trend" => trend))
            
            println("$binned_aligned_file aligned and replaced: $(round(Int,time() -t1))s")
        end
    end
end


"""
    replace_with_model!(dh, nobs, geotiles2replace::AbstractArray; mission_replace="hugonnet", mission_ref1="icesat2", mission_ref2="icesat")

Replace elevation change data for specified geotiles with model-fitted values.

# Arguments
- `dh`: Dictionary of DimensionalArrays containing elevation change data by mission
- `nobs`: Dictionary of DimensionalArrays containing observation counts by mission
- `geotiles2replace`: Array of geotile IDs to replace with model-fitted values
- `mission_replace`: Mission whose data will be replaced (default: "hugonnet")
- `mission_ref1`: Primary reference mission for model fitting (default: "icesat2")
- `mission_ref2`: Secondary reference mission for model fitting (default: "icesat")

# Returns
- Modified `dh` dictionary with replaced values for specified geotiles

# Description
Fits elevation change models to reference mission data and uses these models to replace
values in the target mission for specified geotiles. Uses a combination of reference
missions to ensure data coverage, with the primary reference taking precedence.
"""
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


"""
    geotile2glacier!(glaciers, var0; varname)

Transfer geotile-based data to glacier-level data using area-weighted averaging.

# Arguments
- `glaciers`: DataFrame containing glacier information with columns for geotile and area_km2
- `var0`: DimArray containing geotile-based data with dimensions for geotile, date, and height
- `varname`: Symbol or String specifying the column name to store results in the glaciers DataFrame

# Description
This function calculates area-weighted averages of geotile data for each glacier and stores
the results in the glaciers DataFrame. It processes geotiles in parallel for efficiency and
handles missing data appropriately.

# Returns
- Modified `glaciers` DataFrame with added `varname` column containing time series data
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
        valid_rows, _ = Altim.validrange(.!isnan.(v0))
        
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
            valid_cols, = Altim.validrange(area_index)
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