# This script performs elevation data processing and analysis for glacier studies.
# Key processing steps:
#
# 1. Configuration (~1 min):
#    - Sets up project parameters (ID, geotile width, paths)
#    - Configures processing flags and thresholds
#    - Defines surface masks, DEM sources, and statistical methods
#    - Sets parameters for gap filling and mission alignment
#
# 2. Geotile Binning (~30 min):
#    - Bins elevation data into geotiles using Altim.geotile_binning()
#    - Processes multiple surface types (glacier, land etc)
#    - Applies different binning methods and corrections
#    - Handles filtered and unfiltered data streams
#
# 3. Gap Filling (~2 hrs):
#    - Fills data gaps using Altim.geotile_binned_fill()
#    - Applies multiple parameter sets for robust filling
#    - Handles amplitude corrections and mission normalization
#    - Generates diagnostic plots if requested
#
# 4. Mission Alignment (~1 hr):
#    - Aligns data between ICESat-2 and ICESat missions
#    - Corrects for land surface trends
#    - Replaces data with models in specific regions
#    - Ensures consistent measurements across missions
#
# 5. Regional Aggregation (~3 min):
#    - Calculates regional volume changes
#    - Aggregates geotile data to larger regions
#    - Processes multiple surface types and parameter sets
#
# Total runtime is approximately 3.5 hours, with most time spent on
# gap filling and mission alignment steps.

@time begin
    # set force_remake == true to redo all steps from scratch
    force_remake = false;

    # This section initializes the GEMB (Glacier Energy Mass Balance) model:
    #
    # Key steps:
    # 1. Loads required packages for data processing, statistics, and file I/O
    # 2. Sets model parameters:
    #    - Project ID and version
    #    - Geotile width and buffer distances
    #    - Date ranges from 1940 onwards
    #    - Height bins and elevation deltas
    #    - Precipitation scaling factors
    begin
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
    end

    # Calls Altim.geotile_binning() to bin elevation data into geotiles
    # Key parameters:
    # - project_id: Version identifier
    # - geotile_width: Width of geotiles in degrees
    # - force_remake: Whether to regenerate existing results
    # - update_geotile: Whether to update specific geotiles/missions
    # - surface_masks: Different terrain types to process
    # - binned_folders: Output folders for filtered/unfiltered data
    # - dem_ids: DEM sources to use
    # - binning_methods: Statistical methods for binning
    # - curvature_corrects: Whether to apply curvature correction
    # - max_canopy_height: Maximum canopy height (fixed)
    # - dh_max: Maximum elevation change threshold
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

    # Calls Altim.geotile_binned_fill() to fill gaps in binned elevation data
    # Key parameters:
    # - project_id: Version identifier
    # - geotile_width: Width of geotiles in degrees
    # - force_remake: Whether to regenerate existing results
    # - update_geotile: Whether to update specific geotiles/missions
    # - plot_dh_as_function_of_time_and_elevation: Whether to plot elevation changes
    # - mission_reference_for_amplitude_normalization: Reference mission for normalization
    # - surface_masks: Different terrain types to process
    # - binned_folders: Output folders for filtered/unfiltered data
    # - dem_ids: DEM sources to use
    # - binning_methods: Statistical methods for binning
    # - curvature_corrects: Whether to apply curvature correction
    # - paramater_sets: Different parameter sets for filling
    # - amplitude_corrects: Whether to apply amplitude correction
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

    #Altim.geotile_align_replace() performs mission alignment and data replacement for 
    # elevation data across different satellite missions. Here's a breakdown of its main 
    # purposes:
    # Mission Alignment:
    # Aligns data between two reference missions (ICESat-2 and ICESat)
    # Uses ICESat-2 as the primary reference mission (mission_ref1)
    # Requires a minimum number of trends (min_trend_count = 5) for alignment
    #
    # Land Surface Trend Correction:
    # Can remove land surface trends using Altim.mission_land_trend()
    #
    # Model Replacement:
    # Can replace data in specific regions (here "rgi19") with model data
    #
    # Replaces "hugonnet" data with model for specified regions
    #  This is likely for filling gaps or replacing less reliable measurements
    #
    # Processing Parameters:
    # Works with multiple surface masks (glacier, land, etc.)
    # Processes data across different DEM sources
    # Applies various binning methods and corrections (curvature, amplitude)
    # Organizes data in geotiles (2-degree width)
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
    
    # [This might not be valid anymore... this is now handled in synthesis.jl]
    # Calculates regional volume change by aggregating geotile data [~3 min]
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