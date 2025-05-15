"""
    geotile_binning.jl

Process altimetry data into binned geotiles for glacier elevation change analysis.

This script:
1. Bins raw altimetry measurements into geotiles with consistent spatial and temporal dimensions
2. Fills data gaps through interpolation and extrapolation
3. Aligns and normalizes data from different altimetry missions
4. Replaces problematic data with model-based estimates where necessary

The script processes multiple altimetry missions (ICESat-2, ICESat, GEDI, Hugonnet) with
various filtering and processing options to create a consistent global dataset.
"""
begin
    using Altim
    using Dates
        
    # remake files that were created before this date
    force_remake_before = Date(2025, 3, 15) # set == nothing to disable

    project_id = :v01
    geotile_width = 2

    gemb_file = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0_geotile_filled_d_reg.jld2";

    binned_folder_filtered = Altim.analysis_paths(; geotile_width).binned
    binned_folder_unfiltered = replace(binned_folder_filtered, "binned" => "binned_unfiltered")

    warnings = false 
    showplots = false
    
    # run parameters
    update_geotile = true; # this will load in prevous results to update select geotiles or missions [only touched if force_remake_binning = true]
    update_geotile_missions = ["icesat2"]

    # run parameters
    all_permutations_for_glacier_only = true
    surface_masks = [:glacier, :glacier_rgi7]; #[:glacier, :land, :glacier_rgi7, :glacier_b1km, :glacier_b10km]
    binned_folders = (binned_folder_unfiltered, binned_folder_filtered)
    dem_ids = [:best, :cop30_v2]
    binning_methods = ["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"]
    curvature_corrects = [false,true]
    max_canopy_height = 1
    dh_max=200

    # filling only parameters
    filling_paramater_sets = [1, 2, 3, 4]
    amplitude_corrects = [true]
    plot_dh_as_function_of_time_and_elevation = false;
    mission_reference_for_amplitude_normalization = "icesat2"
end
                                


"""
    geotile_binning(; project_id, geotile_width, warnings, showplots, force_remake_before, 
                     update_geotile, update_geotile_missions, all_permutations_for_glacier_only, 
                     surface_masks, binned_folders, dem_ids, binning_methods, curvature_corrects, 
                     max_canopy_height, dh_max)

Bin raw altimetry measurements into geotiles with consistent spatial and temporal dimensions.

This function processes altimetry data into standardized geotiles, applying various filtering 
and processing options. It currently saves output as JLD2 files, which should eventually be 
changed to NetCDF format for better compatibility.

# Arguments
- `project_id`: Identifier for the project
- `geotile_width`: Width of geotiles in degrees
- `warnings`: Whether to display warnings
- `showplots`: Whether to display plots
- `force_remake_before`: Remake files created before this date
- `update_geotile`: Whether to update existing geotiles
- `update_geotile_missions`: List of missions to update
- `all_permutations_for_glacier_only`: Whether to process all parameter combinations for glacier surfaces
- `surface_masks`: List of surface types to process
- `binned_folders`: Output folders for binned data
- `dem_ids`: Digital elevation models to use
- `binning_methods`: Methods for binning data
- `curvature_corrects`: Whether to apply curvature corrections
- `max_canopy_height`: Maximum canopy height to consider (do not change)
- `dh_max`: Maximum elevation difference filter
"""
Altim.geotile_binning(; 
    project_id,
    geotile_width,
    warnings,
    showplots,

    # run parameters
    force_remake_before,
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


"""
    geotile_binned_fill(; project_id, geotile_width, force_remake_before, update_geotile=false, 
                         update_geotile_missions=nothing, plot_dh_as_function_of_time_and_elevation, 
                         mission_reference_for_amplitude_normalization, all_permutations_for_glacier_only, 
                         surface_masks, binned_folders, dem_ids, binning_methods, curvature_corrects, 
                         paramater_sets, amplitude_corrects, showplots, show_times=false)

Fill data gaps in binned geotiles through interpolation and extrapolation.

This function processes previously binned altimetry data to fill temporal and spatial gaps,
normalize data between different missions, and apply amplitude corrections where needed.

# Arguments
- `project_id`: Identifier for the project
- `geotile_width`: Width of geotiles in degrees
- `force_remake_before`: Remake files created before this date
- `update_geotile`: Whether to update existing geotiles
- `update_geotile_missions`: List of missions to update (all missions must be updated if ICESat-2 is updated)
- `plot_dh_as_function_of_time_and_elevation`: Whether to generate diagnostic plots
- `mission_reference_for_amplitude_normalization`: Reference mission for normalizing amplitude
- `all_permutations_for_glacier_only`: Whether to process all parameter combinations for glacier surfaces
- `surface_masks`: List of surface types to process
- `binned_folders`: Input/output folders for binned data
- `dem_ids`: Digital elevation models to use
- `binning_methods`: Methods for binning data
- `curvature_corrects`: Whether to apply curvature corrections
- `paramater_sets`: Parameter sets for filling algorithms
- `amplitude_corrects`: Whether to apply amplitude corrections
- `showplots`: Whether to display plots
- `show_times`: Whether to display processing times
"""
Altim.geotile_binned_fill(; project_id,
    geotile_width,
    force_remake_before,
    update_geotile = false, # this will load in prevous results to update select geotiles or missions
    update_geotile_missions = nothing, # if icesat2 is updated, all other missions must be updated
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
    show_times=false,
)


"""
    geotile_align_replace(; mission_ref1, mission_ref2, min_trend_count, remove_land_surface_trend,
                          project_id, surface_masks, binned_folders, dem_ids, binning_methods, 
                          curvature_corrects, paramater_sets, amplitude_corrects, geotile_width,
                          regions2replace_with_model, mission_replace_with_model, showplots,
                          force_remake_before)

Align altimetry data between missions and replace problematic data with model-based estimates.

This function performs two key operations:
1. Aligns data between different altimetry missions (e.g., ICESat-2 and ICESat) to ensure consistency
2. Replaces problematic altimetry data in specific regions with model-based estimates

# Arguments
- `mission_ref1`: Primary reference mission (e.g., "icesat2")
- `mission_ref2`: Secondary reference mission (e.g., "icesat")
- `min_trend_count`: Minimum number of observations required for trend calculation
- `remove_land_surface_trend`: Whether to remove land surface trends
- `project_id`: Identifier for the project
- `surface_masks`: List of surface types to process
- `binned_folders`: Input/output folders for binned data
- `dem_ids`: Digital elevation models to use
- `binning_methods`: Methods for binning data
- `curvature_corrects`: Whether to apply curvature corrections
- `paramater_sets`: Parameter sets for filling algorithms
- `amplitude_corrects`: Whether to apply amplitude corrections
- `geotile_width`: Width of geotiles in degrees
- `regions2replace_with_model`: Regions where altimetry data should be replaced with model data
- `mission_replace_with_model`: Mission data to use for model-based replacement
- `showplots`: Whether to display diagnostic plots
- `force_remake_before`: Remake files created before this date
"""
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
    force_remake_before,
)