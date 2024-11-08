
#TODO: Pass arguments as an object using Base.@kwarg

#@time begin

# load packages
using Altim

    project_id = :v01
    geotile_width = 2

    force_remake_binning = false
    force_remake_fill = false

    gemb_file = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0_geotile_filled_d_reg.jld2";

    binned_folder_filtered = analysis_paths(; geotile_width).binned
    binned_folder_unfiltered = replace(binned_folder_filtered, "binned" => "binned_unfiltered")

    warnings = false 
    showplots = false
    
    # run parameters
    update_geotile = false; # this will load in prevous results to update select geotiles or missions
    update_geotile_missions = ["icesat2"]

    # run parameters
    all_permutations_for_glacier_only = false
    surface_masks = [:glacier, :land, :glacier_rgi7, :glacier_b1km, :glacier_b10km]
    surface_masks = [:glacier, :land, :glacier_rgi7]; #, :glacier_b1km, :glacier_b10km]
    binned_folders= (binned_folder_unfiltered, binned_folder_filtered)
    dem_ids = [:best, :cop30_v2]
    binning_methods = ["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"]
    curvature_corrects = [true, false]
    max_canopy_height = 1
    dh_max=200

    # filling only parameters
    filling_paramater_sets = [1, 2, 3, 4]
    amplitude_corrects = [true, false]
    force_remake_fill  = false
    plot_dh_as_function_of_time_and_elevation = true;
    mission_reference_for_amplitude_normalization = "icesat2"
    

    # ------------------------------------------------------------------------------------------
    if false
        force_remake_binning = false
        force_remake_fill = false

        update_geotile = false # this will load in prevous results to update select geotiles or missions [only if force_remake = true]
        update_geotile_missions = nothing

        all_permutations_for_glacier_only = true
        surface_masks = [:glacier]
        binned_folders = (binned_folder_filtered, binned_folder_unfiltered)
        dem_ids = [:best]
        binning_methods = ["meanmadnorm3"]
        curvature_corrects = [true]


        filling_paramater_sets = [1,2,3,4]
        amplitude_corrects = [true]



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

    # align missions 
    @time Altim.geotile_align_replace(;
        mission_ref1 = "icesat2",
        mission_ref2 = "icesat",
        min_trend_count = 5,
        remove_land_surface_trend=Altim.mission_land_trend(),
        project_id,
        surface_masks,
        binned_folders,
        dem_ids,
        binning_methods,
        curvature_corrects,
        paramater_sets = filling_paramater_sets,
        amplitude_corrects,
        geotile_width = 2,
        regions2replace_with_model = ["rgi19"],
        mission_replace_with_model="hugonnet",
        showplots = false,
        force_remake = false,
    )
    
    # calculate regional volume change [~3 min]
     @time Altim.geotile_regional_dv(;
        project_id,
        geotile_width,
        surface_masks,
        binned_folders,
        dem_ids,
        binning_methods,
        curvature_corrects,
        paramater_sets=filling_paramater_sets,
        amplitude_corrects,
        force_remake=false,
    )
end