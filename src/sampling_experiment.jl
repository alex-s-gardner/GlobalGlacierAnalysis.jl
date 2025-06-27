# binning parameters
begin
    import GlobalGlacierAnalysis as GGA
    using Dates

    force_remake = false # it is best practice to delete files than to use force_remake
    project_id = :v01
    geotile_width = 2

    force_remake_before = nothing
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
  
end


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
    curvature_corrects, paramater_sets,
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