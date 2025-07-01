import GlobalGlacierAnalysis as GGA
using DimensionalData
using Dates
using ProgressMeter

begin
    force_remake = false # it is best practice to delete files than to use force_remake
    project_id = :v01
    geotile_width = 2
    domain = :glacier # :glacier -or- :landice
    missions = (:icesat2, :icesat, :gedi, :hugonnet,) # (:icesat2, :icesat, :gedi, :hugonnet)
    force_remake_before = nothing
   
    mission_reference_for_amplitude_normalization = "icesat2"
    all_permutations_for_glacier_only = true
    surface_masks = [:glacier] #[:glacier, :land, :glacier_rgi7, :glacier_b1km, :glacier_b10km]
    #surface_masks = [:glacier_rgi7]
    binned_folders = [GGA.analysis_paths(; geotile_width).binned] # [GGA.analysis_paths(; geotile_width).binned, replace(GGA.analysis_paths(; geotile_width).binned, "binned" => "binned_unfiltered")]
    binning_methods = ["nmad5"]
    #binning_methods = ["nmad3"]
    dem_ids = [:best]
    curvature_corrects = [true]
    fill_params = [2]
    #fill_params = [1]
    amplitude_corrects = [true]
    remove_land_surface_trend = GGA.mission_land_trend()
    regions2replace_with_model = ["rgi19"]
    mission_replace_with_model = "hugonnet"
    missions2align2 = ["icesat2", "icesat"]
    update_missions = nothing
    process_time_show = false
    plots_show = true
    plots_save = true
    #plots_save = false
    plot_save_format = ".png"
    geotiles2plot = GGA.geotiles_golden_test
    single_geotile_test = GGA.geotiles_golden_test[1] # nothing
    #single_geotile_test =  "lat[-32-30]lon[-072-070]"

    # for sub-sampling experiment
    nsamples = 100;
    subsample_fraction = 0.70;


     missions2include=["hugonnet", "gedi", "icesat", "icesat2"]
    downscale_to_glacier_method = "area"
    gemb_run_id = 4
    binned_file_for_over_land_mission_error = "/mnt/bylot-r3/data/binned_unfiltered/2deg/land_dh_best_cc_nmad5_v01.jld2"

    # file for synthesis error
    path2runs_filled, params = GGA.binned_filled_filepaths(; 
        project_id, surface_masks=["glacier", "glacier_rgi7"], 
        dem_ids=["best", "cop30_v2"], 
        curvature_corrects=[false, true], 
        amplitude_corrects=[true], 
        binning_methods=["median", "nmad3", "nmad5"], 
        fill_params=[1, 2, 3, 4], 
        binned_folders=[GGA.analysis_paths(; geotile_width).binned, replace(GGA.analysis_paths(; geotile_width).binned, "binned" => "binned_unfiltered")], 
        include_existing_files_only=true
        )

end


# this is to make the curvature correction figure included in the methods section of the paper
if false # !!! THIS CURRENTLY DOES NOT WORK BECAUSE ICESAT-2 RAW DATA WAS DELETED BY MISTAKE !!!
    GGA.geotile_binning(;
        project_id,
        geotile_width,
        warnings,
        plots_show,

        # run parameters
        force_remake_before=Date(2026, 6, 3),
        geotiles2update=GGA.geotiles_golden_test[1:1],
        update_missions=nothing, #["icesat"],

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
end


GGA.geotile_binned_fill(;
    project_id,
    geotile_width,
    force_remake_before,
    update_missions, # All missions must update if ICESat-2 updates
    mission_reference_for_amplitude_normalization,
    all_permutations_for_glacier_only,
    surface_masks,
    binned_folders,
    binning_methods,
    dem_ids,
    curvature_corrects, fill_params,
    amplitude_corrects, remove_land_surface_trend,
    regions2replace_with_model,
    mission_replace_with_model,
    missions2align2,
    process_time_show,
    plots_show,
    plots_save,
    plot_save_format,
    geotiles2plot,
    single_geotile_test#GGA.geotiles_golden_test[1], #GGA.geotiles_golden_test[2],
)


mission_error = GGA.binned_mad_mission(binned_file_for_over_land_mission_error)


(dh_err, files_included) = GGA.geotile_synthesis_error(;
    path2runs = path2runs_filled,
    outfile="/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2",
    mission_error,
    update_missions=nothing,
    single_geotile_test, #GGA.geotiles_golden_test[1], #GGA.geotiles_golden_test[2],
    force_remake_before=nothing,
);


begin
    

    path2runs_filled0, params = GGA.binned_filled_filepaths(; project_id, surface_masks, dem_ids, curvature_corrects, amplitude_corrects, binning_methods, fill_params, binned_folders, include_existing_files_only=true)
    path2runs_synthesized = replace.(path2runs_filled0, "aligned.jld2" => "synthesized.jld2")
end


GGA.geotile_synthesize_runs(;
    path2runs=path2runs_filled0,
    dh_err,
    missions2include,
    geotiles2plot=nothing,
    single_geotile_test, #GGA.geotiles_golden_test[1], #GGA.geotiles_golden_test[2],
    force_remake_before=nothing
);

# run sub-sampling experiment
#begin
    dh = GGA._simrun_init(;nsamples, missions2include, single_geotile_test)

    @showprogress desc = "Running sampling experiment ..." for i in 1:nsamples
    
        dh0, _, _ = GGA.geotile_binned_fill(;
            project_id,
            geotile_width,
            force_remake_before = nothing,
            update_missions = nothing, # All missions must update if ICESat-2 updates
            mission_reference_for_amplitude_normalization,
            all_permutations_for_glacier_only,
            surface_masks,
            binned_folders,
            binning_methods,
            dem_ids,
            curvature_corrects, 
            fill_params,
            amplitude_corrects, 
            remove_land_surface_trend,
            regions2replace_with_model,
            mission_replace_with_model,
            missions2align2,
            process_time_show,
            plots_show = false,
            plots_save = false,
            plot_save_format = nothing,
            geotiles2plot = nothing,
            single_geotile_test,
            subsample_fraction#GGA.geotiles_golden_test[1], #GGA.geotiles_golden_test[2],
        );

        dh0["Synthesis"], dh_synth_err =  GGA.geotile_synthesize_runs(;
            path2runs = nothing,
            dh_err,
            missions2include,
            geotiles2plot=nothing,
            single_geotile_test=nothing, #GGA.geotiles_golden_test[1], #GGA.geotiles_golden_test[2],
            dh_override=deepcopy(dh0), # used to override the dh_hyps data with a different set of data
            force_remake_before=nothing,
        );

        for mission in vcat(missions2include, "Synthesis")
            dh[mission][run=At(i), geotile=At(single_geotile_test)] = dh0[mission]
        end
    end

    dh_area_average_median, dh_area_average_error = GGA._simrun2areaaverage(deepcopy(dh); surface_mask=:glacier, geotile_width=2, geotile2extract=single_geotile_test);
    f = GGA.plot_area_average_height_anomaly_with_error(dh_area_average_median, dh_area_average_error; mission_order=mission_order =["Synthesis"], median_in_label=false);
    display(f)
    
    fig_folder = joinpath(GGA.pathlocal[:figures], splitpath(binned_folders[1])[end-1])
    data_dir = joinpath(GGA.pathlocal[:project_dir], splitpath(binned_folders[1])[end-1])
    fname = replace(splitpath(path2runs_filled[1])[end], "aligned.jld2" => "synthesis_$(single_geotile_test)_simulation_sf$(subsample_fraction)_n$(nsamples).png")
    save(joinpath(fig_folder, fname), f)
    
    save(joinpath(data_dir, replace(fname, ".png" => ".jld2")), Dict("dh_area_average_median" => dh_area_average_median, "dh_area_average_error" => dh_area_average_error))
end