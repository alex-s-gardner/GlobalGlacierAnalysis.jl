begin
    using Altim
    using Arrow
    using DataFrames
    using NearestNeighbors
    using Plots
    using Statistics
    using DimensionalData
    using Dates
    using LsqFit
    using FileIO

    project_id = :v01
    geotile_width = 2;
    paths = project_paths(; project_id)

    binning_method = "meanmadnorm3"

    binned_folder = "/mnt/bylot-r3/data/binned/2deg"
    dem_id = :best
    binning_methods = "meanmadnorm3"
    curvature_correct = true
    paramater_set = 1
    amplitude_correct = true
end

#@time begin



   
        # load land surface data 
        surface_mask = :glacier

        binned_file_land = Altim.binned_filepath(; binned_folder, surface_mask = :land, dem_id, binning_method, project_id, curvature_correct)
        
        binned_file = Altim.binned_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct)

        binned_file_corrected = replace(binned_file, ".jld2" => "_adjusted.jld2")

        if !isfile(binned_file_land) || !isfile(binned_file)
            continue
        end

        dh_land = FileIO.load(binned_file_land, "dh_hyps")
        nobs_land = FileIO.load(binned_file_land, "nobs_hyps")
        dh = FileIO.load(binned_file, "dh_hyps")

        @time Altim.geotile_adjust!(
            dh["hugonnet"],
            dh_land["hugonnet"],
            nobs_land["hugonnet"],
            dh_land["icesat2"];
            ref_madnorm_max = 10,
            minnobs = 45,
            ref_minvalid = 12,
            )