# run_all.jl
#
# Main workflow script for the GlobalGlacierAnalysis glacier and river analysis pipeline.
#
# This script orchestrates the entire data processing workflow by executing component scripts in sequence.
# Each major step is numbered below for clarity and reproducibility.
#
# Steps:
#  1.  Build satellite altimetry archives (GEDI, ICESat-2, ICESat)
#  2.  Process Hugonnet glacier elevation change data
#  3.  Validate existing ancillary data (DEMs, masks) for geotiles
#  4.  Extract DEM data for each geotile
#  5.  Extract mask data (glacier, land ice) for each geotile
#  6.  Extract canopy height data for each geotile
#  7.  Perform hypsometric analysis of elevation data
#  8.  Create glacier model (GEMB) classes for each geotile
#  9.  Perform statistical binning of elevation data
# 10.  Fill, extrapolate, and adjust binned data
# 11.  Synthesize processed geotile data into global/regional datasets
# 12.  Calculate global glacier discharge, including filled Antarctic data with modeled SMB
# 13.  Calibrate GEMB model to altimetry data for grouped geotiles
# 14.  Calibrate GEMB model to altimetry data for each synthesized geotile dataset
# 15.  Generate glacier-level summary files with key statistics and metrics
# 16.  Export geotile-level glacier change trends and amplitudes to GIS-compatible files
# 17.  Route land surface model runoff through river networks
# 18.  Route glacier runoff through river networks
# 19.  Calculate gmax (maximum glacier contribution to river flux)
# 20.  Analyze population affected by glacier-fed river changes
# 21.  Generate point-based figures for gmax visualization
# 22.  Produce regional results for analysis and sharing
# 23.  Generate extended data figures for manuscript
#
# Implements checkpoint logic to avoid redundant processing, allowing efficient restarts.
#
# Approximate processing times (on old RAID @ 100 MB/s):
#   - GEDI: ~4 days
#   - ICESat-2: ~1 week
#   - ICESat: ~3 hours
begin
    import GlobalGlacierAnalysis as GGA
    using Unitful
    Unitful.register(GGA.MyUnits)    
    using Dates
    using DimensionalData

    # Downloading and preprocessing elevation data
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

    # Binning parameters
    force_remake_before = DateTime(2025, 7, 1, 0, 0, 0)
    update_geotile = true # Update select missions from previous results
    missions2update = nothing #["hugonnet"]
    warnings = false
    max_canopy_height = 1
    dh_max = 200

    # Filling parameters
    mission_reference_for_amplitude_normalization = "icesat2"
    all_permutations_for_glacier_only = true
    surface_masks = [:glacier, :glacier_rgi7] #[:glacier, :land, :glacier_rgi7, :glacier_b1km, :glacier_b10km]
    binned_folders = [replace(GGA.analysis_paths(; geotile_width).binned, "binned" => "binned_unfiltered"), GGA.analysis_paths(; geotile_width).binned, replace(GGA.analysis_paths(; geotile_width).binned, "binned" => "binned_unfiltered")]
    binning_methods = ["nmad3", "nmad5", "median"]
    dem_ids = [:best, :cop30_v2]
    curvature_corrects = [true, false]
    fill_params = [1, 2, 3, 4]
    amplitude_corrects = [true]
    remove_land_surface_trend = GGA.mission_land_trend()
    regions2replace_with_model = ["rgi19"]
    missions2replace_with_model = ["hugonnet"]
    missions2align2 = ["icesat2", "icesat"]
    plots_show = false
    plots_save = false
    plot_save_format = ".png"
    geotiles2plot = nothing; #GGA.geotiles_golden_test
    single_geotile_test = nothing #GGA.geotiles_golden_test[1] # nothing 

    # Synthesis parameters
    gemb_run_id = 5
    path2runs_filled, params = GGA.binned_filled_filepaths(;
        project_id,
        surface_masks=["glacier", "glacier_rgi7"],
        dem_ids=["best", "cop30_v2"],
        curvature_corrects=[false, true],
        amplitude_corrects=[true],
        binning_methods=["median", "nmad3", "nmad5"],
        fill_params=[1, 2, 3, 4],
        binned_folders=[GGA.analysis_paths(; geotile_width).binned, replace(GGA.analysis_paths(; geotile_width).binned, "binned" => "binned_unfiltered")],
        include_existing_files_only=true
    )

    path2runs_synthesized = replace.(path2runs_filled, "aligned.jld2" => "synthesized.jld2")

    binned_file_for_over_land_mission_error = "/mnt/bylot-r3/data/binned_unfiltered/2deg/land_dh_best_cc_nmad5_v01.jld2"

    gembfit_outfile_suffix = ""

    binned_synthesized_dv_files = replace.(path2runs_synthesized, ".jld2" => "_gembfit_dv$(gembfit_outfile_suffix).jld2")
    binned_synthesized_dv_file_ref = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_rgi7_dh_best_cc_nmad5_v01_filled_ac_p1_synthesized_gembfit_dv$(gembfit_outfile_suffix).jld2"
    #GGA.match_rgi6_rgi7_attribute_names(GGA.pathlocal.glacier_rgi7_individual)
end

# 1. Build archives from satellite altimetry data (GEDI, ICESat-2, ICESat)
GGA.geotile_build_archive(;
    force_remake,
    project_id,
    geotile_width,
    domain, # :glacier -or- :landice
    missions, # (:icesat2, :icesat, :gedi, :hugonnet)
    single_geotile_test,
)

# 2. Process Hugonnet glacier elevation change data
GGA.geotile_build_hugonnet(;
    # Parameters: user defined 
    hugonnet_datasets,
    force_remake,
    project_id,
    geotile_width,
    domain, # :glacier -or- :landice
    single_geotile_test,
)

# 3. Validate ancillary data that might already exist (DEMs, masks, etc.) for geotiles
GGA.geotile_ancillary_check(;
    project_id,
    missions,
)

# 4. Extract DEM data for each geotile
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

# 5. Extract mask data (glacier, land ice, etc.) for each geotile
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

# 6. Extract canopy height data for each geotile
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

# 7. Perform hypsometric analysis of elevation data for each geotile
GGA.geotile_hyps_extract(;
    force_remake = true,
    geotile_width,
    raster_file = :cop30_v2, #, :nasadem_v1, :rema_v2_10m, :arcticdem_v3_10m]
    domain, # :glacier -or- :landice
    masks=[:glacier_rgi7] #[:glacier, :land, :glacier_rgi7]
)

# 9. Perform statistical binning of elevation data for each geotile
GGA.geotile_binning(;
    project_id,
    geotile_width,
    warnings,
    plots_show,

    # run parameters
    force_remake_before,
    missions2update,# this will load in prevous results to update select geotiles or missions

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

# 10. Fill, extrapolate, and adjust binned data
GGA.geotile_binned_fill(; project_id,
    geotile_width,
    missions2update, # All missions must update if ICESat-2 updates

    mission_reference_for_amplitude_normalization,
    all_permutations_for_glacier_only, 
    surface_masks,
    binned_folders,
    binning_methods,
    dem_ids,
    curvature_corrects,

    fill_params,
    amplitude_corrects, remove_land_surface_trend,
    regions2replace_with_model,
    missions2replace_with_model,
    missions2align2,
    plots_show,
    plots_save,
    plot_save_format,
    geotiles2plot,
    
    single_geotile_test,

    force_remake_before=DateTime("2025-07-14T01:00:00") + GGA.local2utc
)

# 11. Synthesize geotile data by combining multiple altimetry missions and applying error corrections
# ~36 min for 192 runs 
GGA.geotile_synthesize(path2runs_filled;
    error_file="/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2",
    mission_error=GGA.binned_mad_mission(binned_file_for_over_land_mission_error),
    missions2update=nothing,
    force_remake_before=DateTime("2025-07-16T8:00:00") + GGA.local2utc,
)

# 8. Create glacier model (GEMB) classes for each geotile
include("gemb_classes_binning.jl")

# parms_ref = GGA.binned_filled_fileparts(binned_synthesized_dv_file_ref)


force_remake_before_gemb = DateTime("2026-01-26T22:30") + GGA.local2utc

gemb = GGA.gemb_ensemble_dv(; gemb_run_id);

# 12. Calculate global glacier solid ice discharge, filled Antarctic data with modeled SMB
discharge = GGA.global_discharge_filled(;
    surface_mask="glacier", # = parms_ref.surface_mask, # having issues changing this to :glacier_rgi7
    discharge_global_fn=GGA.pathlocal[:discharge_global],
    gemb,
    discharge2smb_max_latitude = -60,
    discharge2smb_equilibrium_period=(Date(1979), Date(2000)),
    pscale=1,
    ΔT=1,
    geotile_width=2,
    force_remake_before=force_remake_before_gemb,
    force_remake_before_hypsometry=nothing
);


discharge0 = GGA.discharge2geotile(discharge, dims(gemb,:geotile); mass2volume = true);
gemb[:dv][:] = @d gemb[:dv] .- discharge0[:discharge];


# 13. Calibrate GEMB (Glacier Energy and Mass Balance) model to altimetry data by finding optimal fits for grouped geotiles 
# [5 hrs on 128 cores, 150GB of memory]
@time GGA.gemb_calibration(
    path2runs_synthesized, 
    gemb;
    geotile_width,
    geotile_grouping_min_feature_area_km2=100, 
    single_geotile_test = nothing,
    seasonality_weight=GGA.seasonality_weight,
    distance_from_origin_penalty=GGA.distance_from_origin_penalty,
    ΔT_to_pscale_weight = GGA.ΔT_to_pscale_weight,
    force_remake_before = force_remake_before_gemb,
    calibrate_to=:all
)


# 14. Calibrate GEMB model to altimetry data for each synthesized geotile dataset [3.5 min, 150GB of memory]
GGA.geotile_synthesis_gembfit_dv(
    path2runs_synthesized, 
    discharge0,
    gemb;
    geotile_grouping_min_feature_area_km2 = 100, 
    geotile_width, 
    force_remake_before=DateTime(2027,1,1) #force_remake_before_gemb,
)


# 15. Generate glacier level summary nc file with key statistics and metrics for further analysis and sharing  [8 min]
glacier_summary_file = GGA.glacier_summary_file(
    binned_synthesized_dv_files,
    binned_synthesized_dv_file_ref;
    error_quantile=0.95,
    geotile_width,
    error_scaling=1.5, 
    reference_period=(DateTime(2000, 4, 1), DateTime(2024, 12, 31)), 
    surface_mask="glacier",  # = parms_ref.surface_mask, # having issues changing this to :glacier_rgi7
    force_remake_before=DateTime(2027, 1, 1) #force_remake_before_gemb,
)
    

# 16. Export geotile-level glacier change trends and amplitudes to GIS-compatible files [13 s]
GGA.gembfit_dv2gpkg(binned_synthesized_dv_file_ref; 
    outfile_prefix="Gardner2025_geotiles_rates", 
    datelimits=(DateTime(2000, 3, 1), DateTime(2025, 1, 1))
)

# 17. Route land surface model runoff through river networks for river flux calculation
include("land_surface_model_routing.jl")

# 18. Route glacier runoff through river networks for glacier flux calculation [5 min]
GGA.glacier_routing(; surface_mask = "glacier"); # !!! having issues with reading DEM.. leave as rgi6 for now !!!!

# 19. Calculate gmax (maximum glacier contribution to river flux) [3.4 min]
include("gmax_global.jl")

# 20. Generate point-based figures for gmax visualization [1 min]
include("gmax_point_figure.jl")

# 21. Generate regional results for further analysis and sharing [7 min]
include("regional_results.jl")

# 22. Generate extended data figures for manuscript
include("manuscript_extended_data_figures.jl")

# 23. Analyze population affected by glacier-fed river changes using buffer analysis [~5hrs]
include("river_buffer_population.jl")

# 24. Create output files for the manuscript
include("create_output_files.jl")

# run on local machine to copy over figures to my computer
#=
rsync -rav devon:/mnt/bylot-r3/altim_figs/ /Users/gardnera/Research/24_01_GlobalGlacier/figures/
rsync -rav devon:/mnt/bylot-r3/data/project_data/ ~/data/GlobalGlacier/project_data/
=#