begin
    using Altim
    using FileIO
    using DimensionalData
    using DataFrames
    using Plots
   
    project_id = :v01;
    geotile_width = 2;
  
    force_remake = false
    plot_dh_as_function_of_time_and_elevation = true;
    mission_reference_for_amplitude_normalization = "icesat2"


    all_permutations_for_glacier_only = true
    surface_masks = [:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km]
    hugonnet_filtered_flags = [true, false]
    dem_ids = [:best, :cop30_v2]
    binning_methods = ["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"]
    curvature_corrects = [true, false]
    paramater_sets = [1, 2]
    amplitude_corrects = [true, false]


    # ------------------------------------------------------------------------------------------
    if true
        force_remake = true
        all_permutations_for_glacier_only = true
        surface_masks = [:glacier]
        hugonnet_filtered_flags = [true]
        dem_ids = [:best]
        binning_methods = ["meanmadnorm3"]
        curvature_corrects = [true]
        paramater_sets = [1]
        amplitude_corrects = [true]
    end
    # ------------------------------------------------------------------------------------------

    if plot_dh_as_function_of_time_and_elevation
        showplots = true;
        using Plots
    end
end

geotile_binned_fill(; 
    project_id,
    geotile_width,
    force_remake,
    plot_dh_as_function_of_time_and_elevation,
    mission_reference_for_amplitude_normalization,
    all_permutations_for_glacier_only,
    surface_masks,
    hugonnet_filtered_flags,
    dem_ids,
    binning_methods,
    curvature_corrects,
    paramater_sets,
    amplitude_corrects,
    )

