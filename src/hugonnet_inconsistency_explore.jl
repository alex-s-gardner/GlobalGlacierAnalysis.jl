begin
    using FileIO
    import GlobalGlacierAnalysis as GGA
    import GeometryOps as GO
    using DataFrames
    using CairoMakie
    using Statistics
    using Dates
    using NCDatasets
    using Rasters

    # stacks to compare
    stack_id = "N30E078"
    h1 = "/mnt/bylot-r3/data/hugonnet/HSTACK/001/raw/13_14_15_rgi60/stacks/44N/$(stack_id).nc"
    h2 = replace(h1, "/raw/" => "/raw_filtered/")
    h2 = replace(h2,  ".nc" => "_prefilt.nc")
end;

# compare using Rasters
r1 = Raster(h1, name=:z, crs=EPSG(32644), missingval=NaN)
r2 = Raster(h2, name=:z, crs=EPSG(32644), missingval=NaN)

r1 = resample(r1, to = r2);


dem_names1 = NCDataset(h1)["dem_names"][:];
dem_names2 = NCDataset(h2)["dem_names"][:];

dem_index = parent([findfirst(dem_names1 .== dem) for dem in dem_names2]);

r1 = r1[:, :, dem_index];
dem_names1 = dem_names1[dem_index];

time_indexes = [1, 50, 100, 200];

for time_index in time_indexes
    rx1 = r1[:, :, time_index];
    rx2 = r2[:, :, time_index];
    valid = .!isnan.(rx1) .&& .!isnan.(rx2);

    f = Figure(size=(1200, 500));

    ax1 = Axis(f[1, 1], title="$(stack_id) raw: $(dem_names1[time_index])");
    ax2 = Axis(f[1, 2], title="$(stack_id) prefilt: $(dem_names2[time_index])");
    ax3 = Axis(f[1, 3], title="unfiltered vs filtered");
    heatmap!(ax1, rx1);
    heatmap!(ax2, rx2);
    scatter!(ax3, rx1[valid], rx2[valid]);
    display(f)
end

reference_file1 = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_dh_best_nmad5_v01_filled_ac_p2_aligned.jld2"
model_param1 = load(reference_file1, "model_param")["hugonnet"]

reference_file2 = replace(reference_file1, "_unfiltered" => "")
model_param2 = load(reference_file2, "model_param")["hugonnet"]

var_names = ["offset", "offset_icesat2", "offset_icesat"]

for var_name in var_names
    f = Figure(size=(1200, 500));
    ax1 = Axis(f[1, 1], title="unfiltered $(var_name)");
    ax2 = Axis(f[1, 2], title="filtered $(var_name)");
    ax3 = Axis(f[1, 3], title="unfiltered vs filtered $(var_name)");

    hist!(ax1, model_param1[:,var_name]);
    hist!(ax2, model_param2[:,var_name]);
    scatter!(ax3, model_param1[:,var_name], model_param2[:,var_name]);
    display(f)
end

# explore model fit for a single geotile
geotile_index1 = findfirst(model_param1.geotile .== GGA.geotiles_golden_test[1])
geotile_index2 = findfirst(model_param2.geotile .== GGA.geotiles_golden_test[1])

model_param1[geotile_index1, :param_m1]
model_param2[geotile_index2, :param_m1]

geotile_file = "/mnt/bylot-r3/data/binned_unfiltered/2deg/geotile_$(geotile_id).jld2"
geotile_data = load(geotile_file)
geotile_data["geotile_$(geotile_id)"]["dh_best_nmad5_v01_filled_ac_p2_aligned"][:, :, 1]


begin
    project_id = :v01
    geotile_width = 2
    force_remake_before = nothing
    update_geotile = false
    update_geotile_missions = nothing
    mission_reference_for_amplitude_normalization = "icesat2"
    all_permutations_for_glacier_only = true
    surface_masks = [:glacier]
    binned_folders = ["/mnt/bylot-r3/data/binned_unfiltered/2deg", "/mnt/bylot-r3/data/binned/2deg"]
    binning_methods = ["nmad5"]
    dem_ids = [:best]
    curvature_corrects = [false]
    paramater_sets = [2]
    amplitude_corrects = [true]
    remove_land_surface_trend = GGA.mission_land_trend()
    regions2replace_with_model = ["rgi19"]
    mission_replace_with_model = "hugonnet"
    missions2align2 = ["icesat2", "icesat"]
    process_time_show = false
    plots_show = true
    plots_save = true
    plot_save_format = ".png"
    geotiles2plot = GGA.geotiles_golden_test[2:2]
    single_geotile_test = GGA.geotiles_golden_test[2]
end


GGA.geotile_binning(;
    project_id,
    geotile_width,
    warnings,
    plots_show,
    process_time_show,

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
    single_geotile_test,
)



GGA.geotile_binned_fill(;
    project_id=:v01,
    geotile_width=2,
    force_remake_before=nothing,
    update_geotile=false, # Update select missions from previous results
    update_geotile_missions=nothing, # All missions must update if ICESat-2 updates
    mission_reference_for_amplitude_normalization="icesat2",
    all_permutations_for_glacier_only=true,
    surface_masks=[:glacier],
    binned_folders=["/mnt/bylot-r3/data/binned_unfiltered/2deg", "/mnt/bylot-r3/data/binned/2deg"],
    binning_methods=["nmad5"],
    dem_ids=[:best],
    curvature_corrects=[false],
    paramater_sets=[2], #filling_paramater_sets,
    amplitude_corrects=[true], 
    remove_land_surface_trend=GGA.mission_land_trend(),
    regions2replace_with_model=["rgi19"],
    mission_replace_with_model="hugonnet",
    missions2align2=["icesat2", "icesat"],
    process_time_show=false,
    plots_show=true,
    plots_save=true,
    plot_save_format=".png",
    geotiles2plot=GGA.geotiles_golden_test[2:2],
    single_geotile_test=GGA.geotiles_golden_test[1], #GGA.geotiles_golden_test[2],
)





