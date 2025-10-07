# =============================================================================
# This script generates the Extended Data Figures for the manuscript.
#
# It runs the full analysis pipeline for glacier height change, uncertainty,
# discharge, and model calibration, using the GlobalGlacierAnalysis (GGA) package.
#
# The script is organized by Extended Data Figure number, with each section
# corresponding to a figure or analysis step in the manuscript.
#
# Key steps include:
#   - Setting up analysis parameters and file paths
#   - Running geotile binning, filling, and synthesis
#   - Calculating global glacier discharge
#   - Calibrating the GEMB model
#   - Running ensemble and sub-sampling experiments
#
# NOTE: Some sections (e.g., Extended Data Figure 3) may not run due to missing data.
# =============================================================================

begin
    import GlobalGlacierAnalysis as GGA
    using CairoMakie
    using GeoDataFrames
    using FileIO
    using DimensionalData
    using Dates
    using ProgressMeter

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
    missions2replace_with_model = ["hugonnet"]
    missions2align2 = ["icesat2", "icesat"]
    missions2update = nothing
    plots_show = true
    #plots_save = true
    plots_save = false
    plot_save_format = ".png"
    geotiles2plot = GGA.geotiles_golden_test
    single_geotile_test = GGA.geotiles_golden_test[1] # nothing
    #single_geotile_test = "lat[+78+80]lon[+010+012]"

    # for sub-sampling experiment
    nsamples = 100;
    subsample_fraction = 0.70;

    missions2include=["hugonnet", "gedi", "icesat", "icesat2"]
    downscale_to_glacier_method = "area"
    gemb_run_id = 4
    binned_file_for_over_land_mission_error = "/mnt/bylot-r3/data/binned_unfiltered/2deg/land_dh_best_cc_nmad5_v01.jld2"

    path2runs_filled, params = GGA.binned_filled_filepaths(; project_id, surface_masks, dem_ids, curvature_corrects, amplitude_corrects, binning_methods, fill_params, binned_folders, include_existing_files_only=true)
    path2runs_synthesized = replace.(path2runs_filled, "aligned.jld2" => "synthesized.jld2")

    # file for synthesis error
    path2runs_filled_all_ensembles, params = GGA.binned_filled_filepaths(; 
        project_id, surface_masks=["glacier", "glacier_rgi7"], 
        dem_ids=["best", "cop30_v2"], 
        curvature_corrects=[false, true], 
        amplitude_corrects=[true], 
        binning_methods=["median", "nmad3", "nmad5"], 
        fill_params=[1, 2, 3, 4], 
        binned_folders=[GGA.analysis_paths(; geotile_width).binned, replace(GGA.analysis_paths(; geotile_width).binned, "binned" => "binned_unfiltered")], 
        include_existing_files_only=true
        )

    path2runs_synthesized_all_ensembles = replace.(path2runs_filled_all_ensembles, "aligned.jld2" => "synthesized.jld2")
end;

#[Extended Data Figure 3]
# This is to make the curvature correction figure included in the methods section of the paper
if false # !!! THIS CURRENTLY DOES NOT WORK BECAUSE ICESAT-2 RAW DATA WAS DELETED BY MISTAKE !!!
    GGA.geotile_binning(;
        project_id,
        geotile_width,
        warnings,
        plots_show,

        # run parameters
        force_remake_before=Date(2026, 6, 3),
        geotiles2update=GGA.geotiles_golden_test[1:1],
        missions2update=nothing, #["icesat"],

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

#[Extended Data Figure 4/5/6]
GGA.geotile_binned_fill(;
    project_id,
    geotile_width,
    force_remake_before,
    missions2update, # All missions must update if ICESat-2 updates
    mission_reference_for_amplitude_normalization,
    all_permutations_for_glacier_only,
    surface_masks,
    binned_folders,
    binning_methods,
    dem_ids,
    curvature_corrects, fill_params,
    amplitude_corrects, remove_land_surface_trend,
    regions2replace_with_model,
    missions2replace_with_model,
    missions2align2,
    plots_show,
    plots_save,
    plot_save_format,
    geotiles2plot,
    single_geotile_test#GGA.geotiles_golden_test[1], #GGA.geotiles_golden_test[2],
)

#[Extended Data Figure 7]
begin
    mission_error = GGA.binned_mad_mission(binned_file_for_over_land_mission_error)

    (dh_err, files_included) = GGA.geotile_synthesis_error(;
        path2runs = path2runs_filled_all_ensembles,
        outfile="/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2",
        mission_error,
        missions2update=nothing,
        single_geotile_test, #GGA.geotiles_golden_test[1], #GGA.geotiles_golden_test[2],
        force_remake_before=nothing,
    );

    # example synthesis file
    GGA.geotile_synthesize_runs(;
        path2runs=path2runs_filled,
        dh_err,
        missions2include,
        geotiles2plot=nothing,
        single_geotile_test, #GGA.geotiles_golden_test[1], #GGA.geotiles_golden_test[2],
        plots_show=true,
        plots_save=true,
        force_remake_before=nothing
    );
end

# Calculate global glacier discharge, filled Antarctic data with modeled SMB
discharge = GGA.global_discharge_filled(;
    surface_mask="glacier",
    discharge_global_fn=GGA.pathlocal[:discharge_global],
    gemb_run_id,
    discharge2smb_max_latitude=-60,
    discharge2smb_equilibrium_period=(Date(1979), Date(2000)),
    pscale=1,
    Δheight=0,
    geotile_width=2,
    force_remake_before=nothing
);

vars2load = ["runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze"]
gemb = GGA.gemb_ensemble_dv(; gemb_run_id, vars2load)

# [Extended Data Figure 8]
GGA.gemb_calibration(
    path2runs_synthesized,
    discharge,
    gemb;
    geotile_width,
    geotile_grouping_min_feature_area_km2=100,
    single_geotile_test,
    seasonality_weight=GGA.seasonality_weight,
    distance_from_origin_penalty=GGA.distance_from_origin_penalty,
    plots_show=true,
    plots_save=true,
    force_remake_before=nothing
)

# [Extended Data Figure 9 & 10]
begin
    outfile_prefix = "Gardner2025_geotiles_rates"

    geotiles0 = GeoDataFrames.read(joinpath(GGA.pathlocal.data_dir, "project_data", "$(outfile_prefix)_m3yr.gpkg"))

    # plot a histogram of the rates
    f = GGA.plot_hist_gemb_altim_trend_amplitude(geotiles0)
    display(f)

    outfile = joinpath(GGA.pathlocal[:figures], "hist_gemb_altim_trend_amplitude.png")
    CairoMakie.save(outfile, f)

    f = GGA.plot_hist_pscale_Δheight(geotiles0)
    display(f)

    outfile = joinpath(GGA.pathlocal[:figures], "hist_pscale_dheight.png")
    CairoMakie.save(outfile, f)
end

# [Extended Data Figure 9]
begin
    geotiles2extract = [single_geotile_test]

    dh_area_averaged = GGA.ensemble_area_average_height_anomalies(path2runs_synthesized_all_ensembles; geotiles2extract=[single_geotile_test])

    f = GGA.plot_area_average_height_anomaly_with_ensemble_spread(dh_area_averaged; path2runs_synthesized_reference=path2runs_synthesized, geotiles2plot=[single_geotile_test], ref_period=(Date(2000, 1, 1), Date(2001, 12, 31)), xtickspacing=5, p=0.95)

    plots_show && display(f[single_geotile_test])
    plots_save && GGA.CairoMakie.save(joinpath(GGA.pathlocal[:figures], "dh_area_averaged_ensemble_spread_$(single_geotile_test).png"), f[single_geotile_test])
end

# [Extended Data Figure 10]
begin
    dh = GGA._simrun_init(;nsamples, missions2include, single_geotile_test)

    @showprogress desc = "Running sampling experiment ..." for i in 1:nsamples
    
        dh0, _, _ = GGA.geotile_binned_fill(;
            project_id,
            geotile_width,
            force_remake_before = nothing,
            missions2update = nothing, # All missions must update if ICESat-2 updates
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
            missions2replace_with_model,
            missions2align2,
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
    GGA.CairoMakie.save(joinpath(fig_folder, fname), f);
    
    GGA.FileIO.save(joinpath(data_dir, replace(fname, ".png" => ".jld2")), Dict("dh_area_average_median" => dh_area_average_median, "dh_area_average_error" => dh_area_average_error))
end