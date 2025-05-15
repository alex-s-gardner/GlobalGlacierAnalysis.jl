"""
    geotile_ancillary_check.jl

Validates and cleans ancillary data files associated with altimetry measurements.

This script:
1. Checks that all altimetry, mask, and DEM files exist for each geotile
2. Verifies that ancillary files have the same number of rows as their corresponding altimetry files
3. Deletes any ancillary files with mismatched row counts to prevent data inconsistencies

The script processes multiple altimetry missions (ICESat-2, ICESat, GEDI, Hugonnet) and
various ancillary data types (masks, DEMs, canopy height).
"""

using Arrow
using Altim
using ProgressMeter

# Parameters: user defined 
missions = (:icesat2, :icesat, :gedi, :hugonnet)
suffix2check = ["masks", "cop30_v2", "canopyh", "rema_v2_10m", "arcticdem_v4_10m", "nasadem_v1"]

# Initialize: paths, products, geotiles
paths = project_paths(; project_id);
mission_geotile_folders = [paths[mission].geotile for mission in missions]

# add hugonnet unfiltered to the list
if :hugonnet in missions
    mission_geotile_folders = vcat(mission_geotile_folders, replace(paths[:hugonnet].geotile, "/2deg" => "/2deg_unfiltered"))
end

# Execute: check length of all files
for mission_geotile_folder in mission_geotile_folders

    mission_dir_files = Altim.allfiles(mission_geotile_folder)
    paths2altim = filter(x -> occursin(".arrow", x), mission_dir_files)

    if isempty(paths2altim)
        @warn "no altimetry files found in: $mission_geotile_folder"
        continue
    end

    printstyled("$(mission_geotile_folder)\n", color=:blue)
    @showprogress desc = "checking that number of rows in altimetry == ancillary:$(mission_geotile_folder)..." Threads.@threads for path2altim in paths2altim
     
        geotile_id = replace(splitpath(path2altim)[end], ".arrow" => "")

        paths2ancillary = filter(x -> (occursin(geotile_id, x) && !occursin(".arrow", x)), mission_dir_files)

        n_rows = length(Arrow.Table(path2altim)[1])

        for path2ancillary in paths2ancillary
            n_rows_ancillary = length(Arrow.Table(path2ancillary)[1])
            if n_rows_ancillary != n_rows
                @warn "number of rows in altimetry  != ancillary, deleting: $path2ancillary"
                rm(path2ancillary)
            end
        end
    end
end