using Altim
using FileIO

project_id = :v01
geotile_width = 2

gemb_file = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0_geotile_filled_d_reg.jld2";

binned_folder_filtered = Altim.analysis_paths(; geotile_width).binned
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
binned_folders = (binned_folder_unfiltered, binned_folder_filtered)
dem_ids = [:best, :cop30_v2]
binning_methods = ["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"]
curvature_corrects = [true, false]
max_canopy_height = 1
dh_max = 200

# filling only parameters
filling_paramater_sets = [1, 2, 3, 4]
amplitude_corrects = [true, false]
plot_dh_as_function_of_time_and_elevation = true;
mission_reference_for_amplitude_normalization = "icesat2"

binned_file = Altim.binned_filepath(; binned_folder = binned_folder_filtered, surface_mask = :glacier, dem_id = :best, binning_method = "meanmadnorm3", project_id, curvature_correct = true)






