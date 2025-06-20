# master_run.jl
#
# Main workflow script for the GlobalGlacierAnalysis glacier and river analysis pipeline.
#
# Orchestrates the entire data processing workflow by executing component scripts in sequence:
# 1. Builds satellite altimetry archives (GEDI, ICESat-2, ICESat)
# 2. Processes Hugonnet glacier elevation change data
# 3. Validates existing ancillary data (DEMs, masks) for geotiles
# 4. Extracts DEM data for each geotile
# 5. Extracts mask data (glacier, land ice) for each geotile
# 6. Extracts canopy height data for each geotile
# 7. Performs hypsometric analysis of elevation data
# 8. Creates glacier model (GEMB) classes for each geotile
# 9. Performs statistical binning of elevation data
# 10. Synthesizes processed geotile data into global/regional datasets
# 11. Generates summary files with key statistics and metrics
# 12. Creates visualization plots and GIS outputs
# 13. Routes land surface model runoff through river networks
# 14. Routes glacier runoff through river networks
# 15. Calculates gmax (maximum glacier contribution to river flux)
# 16. Analyzes population affected by glacier-fed river changes
# 17. Generates point-based figures for gmax visualization
# 18. Produces regional results for analysis and sharing
#
# Implements checkpoint logic to avoid redundant processing, allowing efficient restarts.
#
# Processing times (on old RAID @ 100 MB/s):
# - GEDI: ~4 days
# - ICESat-2: ~1 week  
# - ICESat: ~3 hours

use_test_parameters = false
begin
    import GlobalGlacierAnalysis as GGA
    using Dates

    force_remake=false # it is best practice to delete files than to use force_remake
    project_id=:v01
    geotile_width=2
    domain=:glacier # :glacier -or- :landice
    missions= (:icesat2, :icesat, :gedi, :hugonnet,) # (:icesat2, :icesat, :gedi, :hugonnet)
    hugonnet_datasets = ["unfiltered", "filtered"]
    hugonnet_unfiltered=true
    slope=true
    curvature=true
    dems2extract=[:rema_v2_10m, :cop30_v2, :arcticdem_v4_10m, :nasadem_v1]
    masks2extract = [:floatingice, :inlandwater, :landice, :ocean]
    masks2extract_highres = [:glacier, :glacier_rgi7, :land]
    single_geotile_test = nothing

    if use_test_parameters
        @warn "WARNING: TEST PARAMETERS BEING USED"
        force_remake = false
        project_id = :v01
        geotile_width = 2
        domain = :glacier # :glacier -or- :landice
        missions = (:hugonnet,) # (:icesat2, :icesat, :gedi, :hugonnet,) # (:icesat2, :icesat, :gedi, :hugonnet)
        hugonnet_datasets = ["unfiltered"]
        hugonnet_unfiltered = true
        slope = true
        curvature = true
        dems2extract = [:rema_v2_10m, :cop30_v2, :arcticdem_v4_10m, :nasadem_v1]
        masks2extract = [:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean]
        single_geotile_test = GGA.geotiles_golden_test[1]
    end
end


# [1] Build archives from satellite altimetry data (GEDI, ICESat-2, ICESat)
GGA.geotile_build_archive(;
    force_remake,
    project_id,
    geotile_width,
    domain, # :glacier -or- :landice
    missions, # (:icesat2, :icesat, :gedi, :hugonnet)
    single_geotile_test,
)

# [2] Process Hugonnet glacier elevation change data
GGA.geotile_build_hugonnet(;
    # Parameters: user defined 
    hugonnet_datasets,
    force_remake,
    project_id,
    geotile_width,
    domain, # :glacier -or- :landice
    single_geotile_test,
)

# [3] Validate ancillary data that might already exist (DEMs, masks, etc.) for geotiles
GGA.geotile_ancillary_check(;
    project_id,
    missions,
)

# [4] Extract DEM data for each geotile
GGA.geotile_dem_extract(;
    # Parameters: user defined 
    force_remake,
    project_id,
    geotile_width,
    domain, # :glacier -or- :landice
    missions, # (:icesat2, :icesat, :gedi, :hugonnet)
    hugonnet_unfiltered, 
    slope,
    curvature,
    dems2extract,
    single_geotile_test,
)

# [5] Extract mask data (glacier, land ice, etc.) for each geotile
GGA.geotile_mask_extract(;
    # Parameters: user defined 
    force_remake,
    project_id,
    geotile_width,
    domain, # :glacier -or- :landice
    missions, # (:icesat2, :icesat, :gedi, :hugonnet)
    hugonnet_unfiltered,
    masks2extract,
    masks2extract_highres,
    single_geotile_test,
)

# [6] Extract canopy height data for each geotile
GGA.geotile_canopyh_extract(;
    # Parameters: user defined 
    force_remake,
    project_id,
    geotile_width,
    domain, # :glacier -or- :landice
    missions, # (:icesat2, :icesat, :gedi, :hugonnet)
    hugonnet_unfiltered,
    single_geotile_test,
)

# [7] Perform hypsometric analysis of elevation data for each geotile
masks = [:glacier, :land, :glacier_rgi7, :glacier_b1km, :glacier_b10km]
GGA.geotile_hyps_extract(;
    force_remake,
    geotile_width,
    raster_file = :cop30_v2, #, :nasadem_v1, :rema_v2_10m, :arcticdem_v3_10m]
    domain, # :glacier -or- :landice
    masks
)

# [8] Create glacier model (GEMB) classes for each geotile
include("gemb_classes_binning.jl")


# binning parameters
begin
    force_remake_before = DateTime(2025, 6, 20, 9, 0, 0)
    update_geotile = true # Update select missions from previous results
    update_geotile_missions = ["hugonnet"]
    warnings = false
    max_canopy_height = 1
    dh_max = 200

    mission_reference_for_amplitude_normalization = "icesat2"
    all_permutations_for_glacier_only = true
    surface_masks = [:glacier, :glacier_rgi7] #[:glacier, :land, :glacier_rgi7, :glacier_b1km, :glacier_b10km]
    binned_folders = [replace(GGA.analysis_paths(; geotile_width).binned, "binned" => "binned_unfiltered")] # [GGA.analysis_paths(; geotile_width).binned, replace(GGA.analysis_paths(; geotile_width).binned, "binned" => "binned_unfiltered")]
    binning_methods = ["nmad3", "nmad5", "median"]
    dem_ids = [:best, :cop30_v2]
    curvature_corrects = [true, false]

    paramater_sets = [1, 2, 3, 4]
    amplitude_corrects = [true]
    remove_land_surface_trend = GGA.mission_land_trend()
    regions2replace_with_model = ["rgi19"]
    mission_replace_with_model = "hugonnet"
    missions2align2 = ["icesat2", "icesat"]

    process_time_show = false
    plots_show = false
    plots_save = false
    plot_save_format = ".png"
    geotiles2plot = GGA.geotiles_golden_test
    single_geotile_test = nothing #GGA.geotiles_golden_test[1] # nothing 
    
    if use_test_parameters
        @warn "WARNING: TEST PARAMETERS BEING USED"
        force_remake_before = nothing
        update_geotile = false # Update select missions from previous results
        update_geotile_missions = nothing
        warnings = false
        max_canopy_height = 1
        dh_max = 200

        mission_reference_for_amplitude_normalization = "icesat2"
        all_permutations_for_glacier_only = true
        surface_masks = [:glacier] #[:glacier, :land, :glacier_rgi7, :glacier_b1km, :glacier_b10km]
        binned_folders = [GGA.analysis_paths(; geotile_width).binned, replace(GGA.analysis_paths(; geotile_width).binned, "binned" => "binned_unfiltered")]
        binning_methods = ["nmad5"]
        dem_ids = [:best]
        curvature_corrects = [true]

        paramater_sets = [2]
        amplitude_corrects = [true]
        remove_land_surface_trend = GGA.mission_land_trend()
        regions2replace_with_model = ["rgi19"]
        mission_replace_with_model = "hugonnet"
        missions2align2 = ["icesat2", "icesat"]

        process_time_show = false
        plots_show = false
        plots_save = false
        plot_save_format = ".png"
        geotiles2plot = GGA.geotiles_golden_test
        single_geotile_test = GGA.geotiles_golden_test[1] # nothing 
    end
end



# [9] Perform statistical binning of elevation data for each geotile
GGA.geotile_binning(;
    project_id,
    geotile_width,
    warnings,
    plots_show,
    process_time_show=false,

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

    max_canopy_height, # do not change

    # filter parameters
    dh_max,
    single_geotile_test, #GGA.geotiles_golden_test[2],
)

#[10] Fill, extrapolate, and adjust binned data
GGA.geotile_binned_fill(;
    project_id,
    geotile_width,
    force_remake_before,
    update_geotile, # Update select missions from previous results
    update_geotile_missions, # All missions must update if ICESat-2 updates

    mission_reference_for_amplitude_normalization,
    all_permutations_for_glacier_only, 
    surface_masks,
    binned_folders,
    binning_methods,
    dem_ids,
    curvature_corrects,

    paramater_sets,
    amplitude_corrects, remove_land_surface_trend,
    regions2replace_with_model,
    mission_replace_with_model,
    missions2align2,
    process_time_show, 
    plots_show,
    plots_save,
    plot_save_format,
    geotiles2plot,
    single_geotile_test #GGA.geotiles_golden_test[1], #GGA.geotiles_golden_test[2],
)


# [10] Synthesize processed geotile data into global and regional datasets
include("geotile_synthesis.jl")

# [11] Generate summary files with key statistics and metrics for further analysis and sharing
include("summary_files.jl")

# [12] Generate visualization plots and GIS outputs for data analysis and presentation
include("synthesis_plots_gis.jl")

# [13] Route land surface model runoff through river networks for river flux calculation
include("land_surface_model_routing.jl")

# [14] Route glacier runoff through river networks for glacier flux calculation
include("glacier_routing.jl")

# [15] Calculate gmax (maximum glacier contribution to river flux)
include("gmax_global.jl")

# [16] Analyze population affected by glacier-fed river changes using buffer analysis
include("river_buffer_population.jl")

# [17] Generate point-based figures for gmax visualization
include("gmax_point_figure.jl")

# [18] Generate regional results for further analysis and sharing
include("regional_results.jl")
