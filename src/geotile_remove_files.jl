"""
    geotile_remove_files.jl

Utility script for removing specific file types from geotile directories.

This script:
1. Deletes specified file types (e.g., .arrow, .masks, DEM files) for selected missions
2. Optionally removes temporary files (prefixed with 'jl_')
3. Includes safety measures to prevent accidental deletion of all files
4. Requires explicit configuration of missions and file suffixes

The script is designed for maintenance and cleanup of geotile data directories
when specific derived products need to be regenerated.
"""
begin
    import GlobalGlacierAnalysis as GGA
    project_id = :v01;
    paths = GGA.project_paths(; project_id)
    delete_temp_files = true

    #mission = ["icesat2", "icesat",  "gedi",  "hugonnet"]
    #filesufx = [".arrow", ".masks", ".arcticdem_v4_10m", ".cop30_v2", ".rema_v2_10m", ".canopyh", ".offset2cop30_v2", "dhdt_ref"]

    # these are set to nothing as a safety measure to prevent accidental deletion of all files.. 
    # you can set them to the desired values to delete specific files
    missions = nothing; #["hugonnet"]
    filesufxs = nothing; # [".masks", ".arcticdem_v4_10m", ".cop30_v2", ".rema_v2_10m", ".canopyh"]


    begin
        if isnothing(missions) || isnothing(filesufxs)
            printstyled("!! You are about to delete ALL files ending in $(filesufx) for $(mission) missions !!\n\n"; color=:red, bold = true)
            error("You must set either missions and filesufxs to delete specific files")
        else
            printstyled("-------------------------------------------------------------------------------\n", ; color=:red)
            printstyled("-------------------------------------------------------------------------------\n\n", ; color=:red)
            printstyled("!! You are about to delete ALL files ending in $(filesufx) for $(mission) missions !!\n\n"; color=:red, bold = true)
            printstyled("-------------------------------------------------------------------------------\n", ; color=:red)
            printstyled("-------------------------------------------------------------------------------\n", ; color=:red)
        end
    end

    for mission in missions
        for filesufx in filesufxs
            files2delete = GGA.allfiles(paths[Symbol(mission)].geotile; fn_endswith = filesufx)

            if isempty(files2delete)
                continue
            else
                rm.(files2delete)
            end
        end

        if delete_temp_files
            files2delete = GGA.allfiles(paths[Symbol(mission)].geotile; fn_startswith="jl_")

            if isempty(files2delete)
                continue
            else
                rm.(files2delete)
            end
        end
    end
end