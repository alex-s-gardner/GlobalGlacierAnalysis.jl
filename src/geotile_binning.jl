# Process altimetry data into binned geotiles for glacier elevation change analysis.
# This script:
# 1. Bins raw altimetry measurements into geotiles with consistent spatial and temporal dimensions
# 2. Fills data gaps through interpolation and extrapolation
# 3. Aligns and normalizes data from different altimetry missions
# 4. Replaces problematic data with model-based estimates where necessary
#
# The script processes multiple altimetry missions (ICESat-2, ICESat, GEDI, Hugonnet) with
# various filtering and processing options to create a consistent global dataset.

begin
    import GlobalGlacierAnalysis as GGA
    using Dates
        
    # remake files that were created before this date
    force_remake_before = Date(2025, 3, 15) # set == nothing to disable

    project_id = :v01
    geotile_width = 2

    gemb_file = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0_geotile_filled_d_reg.jld2";

    binned_folder_filtered = GGA.analysis_paths(; geotile_width).binned
    binned_folder_unfiltered = replace(binned_folder_filtered, "binned" => "binned_unfiltered")

    warnings = false 
    plots_show = false
    
    # run parameters
    update_geotile = false; # this will load in prevous results to update select geotiles or missions [only touched if force_remake_binning = true]
    update_geotile_missions = ["icesat2"]

    # run parameters
    all_permutations_for_glacier_only = true
    surface_masks = [:glacier, :glacier_rgi7]; #[:glacier, :land, :glacier_rgi7, :glacier_b1km, :glacier_b10km]
    binned_folders = [binned_folder_unfiltered, binned_folder_filtered]
    dem_ids = [:best, :cop30_v2]
    binning_methods = ["nmad3", "nmad5", "median", "nmad10"]
    curvature_corrects = [false,true]
    max_canopy_height = 1
    dh_max = 200

    # filling only parameters
    filling_paramater_sets = [1, 2, 3, 4]
    amplitude_corrects = [true]
    plot_dh_as_function_of_time_and_elevation = false;
    mission_reference_for_amplitude_normalization = "icesat2"
end
        

if false
    @warn "Test parameter set being used"

    # run parameters
    all_permutations_for_glacier_only = true
    surface_masks = [:glacier]; #[:glacier, :land, :glacier_rgi7, :glacier_b1km, :glacier_b10km]
    binned_folders = [binned_folder_filtered]
    dem_ids = [:best]
    binning_methods = ["nmad5"]
    curvature_corrects = [true]
    max_canopy_height = 1
    dh_max = 200

    # filling only parameters
    filling_paramater_sets = [2]
    amplitude_corrects = [true]
    plot_dh_as_function_of_time_and_elevation = false;
    mission_reference_for_amplitude_normalization = "icesat2"
end


# this is to make the curvature correction figure included in the methods section of the paper
if false
    GGA.geotile_binning(; 
        project_id = :v01,
        geotile_width = 2,
        warnings = false,
        plots_show = true,

        # run parameters
        force_remake_before = Date(2026, 6, 3),
        update_geotile = true, # this will load in prevous results to update select geotiles or missions]
        geotiles2update = GGA.geotiles_golden_test[1:1],
        update_geotile_missions = ["icesat"],

        # run parameters
        all_permutations_for_glacier_only = true,
        surface_masks = [:glacier],
        binned_folders=("/mnt/bylot-r3/data/binned/2deg",),
        dem_ids = [:best],
        binning_methods = ["nmad5"],
        curvature_corrects = [true],

        #### DON NOT CHANGE THESE PARAMETERS
        max_canopy_height = 1, # do not change

        # filter parameters
        dh_max=200,
        )
end


# Bin raw altimetry measurements into geotiles with consistent spatial and temporal dimensions.
#
# This function processes altimetry data into standardized geotiles, applying various filtering 
# and processing options. It currently saves output as JLD2 files, which should eventually be 
# changed to NetCDF format for better compatibility.
#
# Arguments
# - `project_id`: Identifier for the project
# - `geotile_width`: Width of geotiles in degrees
# - `warnings`: Whether to display warnings
# - `plots_show`: Whether to display plots
# - `force_remake_before`: Remake files created before this date
# - `update_geotile`: Whether to update existing geotiles
# - `update_geotile_missions`: List of missions to update
# - `all_permutations_for_glacier_only`: Whether to process all parameter combinations for glacier surfaces
# - `surface_masks`: List of surface types to process
# - `binned_folders`: Output folders for binned data
# - `dem_ids`: Digital elevation models to use
# - `binning_methods`: Methods for binning data
# - `curvature_corrects`: Whether to apply curvature corrections
# - `max_canopy_height`: Maximum canopy height to consider (do not change)
# - `dh_max`: Maximum elevation difference filter

GGA.geotile_binning(; 
    project_id,
    geotile_width,
    warnings,
    plots_show,

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



# Process and fill gaps in binned altimetry data
#
# Fills temporal/spatial gaps, normalizes between missions, and applies amplitude corrections
# to previously binned altimetry data.
#
# Required Arguments:
# - `project_id`: Project identifier
# - `geotile_width`: Geotile width in degrees
# - `mission_reference_for_amplitude_normalization`: Reference mission for amplitude normalization
# - `surface_masks`: Surface types to process
# - `binned_folders`: Input/output data folders
# - `dem_ids`: Digital elevation models to use
# - `binning_methods`: Data binning methods
# - `curvature_corrects`: Whether to apply curvature corrections
# - `amplitude_corrects`: Whether to apply amplitude corrections
#
# Optional Arguments:
# - `force_remake_before`: Remake files before this date (default: 2026-06-05)
# - `update_geotile`: Update existing geotiles (default: false)
# - `update_geotile_missions`: Missions to update (default: nothing)
# - `all_permutations_for_glacier_only`: Process all parameter combinations for glaciers (default: true)
# - `paramater_sets`: Parameter sets for filling algorithms (default: filling_paramater_sets)
# - `plots_show`: Display diagnostic plots (default: false)
# - `show_times`: Display processing times (default: false)
GGA.geotile_binned_fill(;
    project_id,
    geotile_width,
    force_remake_before = Date(2026, 6, 5),
    update_geotile = false, # Update select missions from previous results
    update_geotile_missions = nothing, # All missions must update if ICESat-2 updates
    mission_reference_for_amplitude_normalization,
    all_permutations_for_glacier_only,
    surface_masks,
    binned_folders,
    binning_methods,
    dem_ids,
    curvature_corrects,
    paramater_sets = filling_paramater_sets,
    amplitude_corrects,

    remove_land_surface_trend = GGA.mission_land_trend(),
    regions2replace_with_model = ["rgi19"],
    mission_replace_with_model = "hugonnet",
    missions2align2 = ["icesat2", "icesat"],
    show_times = false,

    plots_show = false,
    plots_save = false,
    plot_save_format = ".png",
    geotiles2plot = GGA.geotiles_golden_test[2:2],
    single_geotile_test=nothing, #GGA.geotiles_golden_test[2],
)