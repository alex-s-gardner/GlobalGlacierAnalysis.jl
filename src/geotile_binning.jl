
#TODO: Pass arguments as an object using Base.@kwarg

#@time begin
# load packages
#using Altim

    project_id = :v01
    geotile_width = 2

    binned_folder_filtered = analysis_paths(; geotile_width).binned
    binned_folder_unfiltered = replace(binned_folder_filtered, "binned" => "binned_unfiltered")

    warnings = false 
    showplots = false
    
    # run parameters
    force_remake_binning = false;
    update_geotile = false; # this will load in prevous results to update select geotiles or missions
    update_geotile_missions = ["icesat2"]

    # run parameters
    all_permutations_for_glacier_only = true
    surface_masks = [:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km]
    binned_folders= (binned_folder_filtered, binned_folder_unfiltered)
    dem_ids = [:best, :cop30_v2]
    binning_methods = ["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"]
    curvature_corrects = [true, false]
    max_canopy_height = 1
    dh_max=200

    # filling only parameters
    filling_paramater_sets = [1, 2]
    amplitude_corrects = [true, false]
    force_remake_fill  = true
    plot_dh_as_function_of_time_and_elevation = true;
    mission_reference_for_amplitude_normalization = "icesat2"
    
    # ------------------------------------------------------------------------------------------
    if true
        force_remake_binning = false
        force_remake_fill = false

        update_geotile = true # this will load in prevous results to update select geotiles or missions [only if force_remake = true]
        update_geotile_missions = ["gedi"]

        all_permutations_for_glacier_only = true
        surface_masks = [:glacier]
        binned_folders = (binned_folder_filtered, )
        dem_ids = [:best]
        binning_methods = ["meanmadnorm3"]
        curvature_corrects = [true]
        dh_max=200

        filling_paramater_sets = [1]
        amplitude_corrects = [true]
        plot_dh_as_function_of_time_and_elevation = true;
        mission_reference_for_amplitude_normalization = "icesat2"

    end


    Altim.geotile_binning(; 
        project_id,
        geotile_width,
        warnings,
        showplots,

        # run parameters
        force_remake = force_remake_binning,
        update_geotile, # this will load in prevous results to update select geotiles or missions
        update_geotile_missions,

        # run parameters
        all_permutations_for_glacier_only,
        surface_masks,
        binned_folders,
        dem_ids,
        binning_methods,
        curvature_corrects,

        #### DON NOT CHANGE THESE PARAMETERS
        max_canopy_height, # do not change

        # filter parameters
        dh_max,
)


Altim.geotile_binned_fill(; 
    project_id,
    geotile_width,
    force_remake = force_remake_fill,
    update_geotile, # this will load in prevous results to update select geotiles or missions
    update_geotile_missions,
    plot_dh_as_function_of_time_and_elevation,
    mission_reference_for_amplitude_normalization,
    all_permutations_for_glacier_only,
    surface_masks,
    binned_folders,
    dem_ids,
    binning_methods,
    curvature_corrects,
    paramater_sets = filling_paramater_sets,
    amplitude_corrects,
    showplots,
    )

Altim.geotile_filled_landfit(; 
    project_id = "v01",
    showplots = false,
    binned_folders=(binned_folder_filtered, binned_folder_unfiltered),
    geotile_width=2, 
    paramater_sets=filling_paramater_sets,
    surface_masks=[:land], 
    binning_methods=["meanmadnorm3"], 
    dem_ids=[:best], 
    curvature_corrects=[true], 
    amplitude_corrects=[true], 
    force_remake_landoffset=false
    )

Altim.geotile_bin_fit_fac(;
        gemb_file_binned = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0.jld2",
        firn_mission_ref = "icesat2",
        project_id=:v01,
        geotile_width=2,
        surface_masks=[:glacier, :glacier_rgi7, :glacier_b1km, :glacier_b10km],
        binned_folders=(binned_folder_filtered, binned_folder_unfiltered),
        dem_ids=[:best, :cop30_v2],
        binning_methods=["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"],
        curvature_corrects=[true, false],
        paramater_sets=filling_paramater_sets,
        amplitude_corrects=[true, false],
        force_remake_fac = false,
    )

Altim.geotile_regional_dvdm(;
        firn_mission_ref="icesat2",
        fac_scale_apply=true,
        project_id,
        geotile_width,
        surface_masks,
        binned_folders,
        dem_ids,
        binning_methods,
        curvature_corrects,
        paramater_sets=filling_paramater_sets,
        amplitude_corrects,
        force_remake_masschange = true,
    )


df = Altim.geotile_dvdm_synthesize(;
        # best estimate ,
        surface_mask_best = "glacier",
        dem_best = "best",
        curvature_correct_best = true,
        amplitude_correct_best = true,
        binning_method_best = "meanmadnorm3",
        fill_param_best = 1,
        binned_folder_best = binned_folder_filtered,

        # to include in uncertainty
        surface_masks = ["glacier", "glacier_rgi7", "glacier_b1km"],
        dems = ["best", "cop30_v2"],
        curvature_corrects = [false, true],
        amplitude_corrects = [true],
        binning_methods = ["median", "meanmadnorm5", "meanmadnorm3", "meanmadnorm10"],# ["median", "meanmadnorm10", "meanmadnorm5", "meanmadnorm3"]
        fill_params=filling_paramater_sets,
        binned_folders=(binned_folder_filtered, binned_folder_unfiltered),

        # manual adjustments
        regions_to_overwrite_hugonnet_data_with_model_fit = ["rgi19"],

        ## combine regions as the per-sensor level
        # this gets a bit complicated a only regions with consistnet mission inclusion can be grouped for proper 
        region_col = "rgi",
        region_combines = (
            ("rgi13", "rgi14", "rgi15") => "hma", 
            ("rgi1", "rgi3", "rgi4", "rgi5", "rgi6", "rgi7", "rgi8", "rgi9", "rgi19") => "hll",
            ("rgi2", "rgi10", "rgi11", "rgi12", "rgi13", "rgi14", "rgi15", "rgi17", "rgi18") => "ghll",
            ("rgi1", "rgi3", "rgi4", "rgi6", "rgi7", "rgi8", "rgi9") => "hll_ep",
        ),

        combine_vars =  ["area_km2", "dm_gt", "dv_km3", "nobs", "fac_km3", "smb_km3"],
    )

    df = Altim.geotile_dvdm_addgrace!(df)

    df = Altim.geotile_combine_synth_regions!(df)

    df_dm = Altim.geotile_dvdm_add_trend!(df; iterations = 1000)

    df_dh = Altim.geotile_dvdm_areaaverage(df)

    fig = Altim.plot_regional_dvdm(df_dm;
        rgi="global_ep",
        variable="dm",
        featured_mission="synthesis",
        showmissions=true,
        )


Altim.plot_multiregion_dvdm(df_dm;
    variable = "dm",
    featured_mission = "synthesis",
    regions = reduce(vcat, (["rgi$i" for i in vcat(1:12, 16:19)], ["hma"])),
    showlines = false,
    showmissions = false,
    fontsize = 15,
    cmap=:Dark2_4
    )

end