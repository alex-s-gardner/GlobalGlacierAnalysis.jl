# this script removes specified file from geotile folders. This is useful when recreating 
# geotile files... this is better than using force_remake within each script as it will 
# allow failed jobs to skip already created files

project_id = :v01;
paths = project_paths(; project_id);
delete_temp_files = true

#mission = ["icesat2", "icesat",  "gedi",  "hugonnet"]
#filesufx = [".arrow", ".masks", ".arcticdem_v4_10m", ".cop30_v2", ".rema_v2_10m", ".canopyh", ".offset2cop30_v2", "dhdt_ref"]

missions = ["hugonnet"]
filesufxs = [".masks", ".arcticdem_v4_10m", ".cop30_v2", ".rema_v2_10m", ".canopyh"]


begin
printstyled("-------------------------------------------------------------------------------\n", ; color=:red)
printstyled("-------------------------------------------------------------------------------\n\n", ; color=:red)
printstyled("!! You are about to delete ALL files ending in $(filesufx) for $(mission) missions !!\n\n"; color=:red, bold = true)
printstyled("-------------------------------------------------------------------------------\n", ; color=:red)
printstyled("-------------------------------------------------------------------------------\n", ; color=:red)
end

for mission in missions
    for filesufx in filesufxs
        files2delete = allfiles(paths[Symbol(mission)].geotile; fn_endswith = filesufx)

        if isempty(files2delete)
            continue
        else
            rm.(files2delete)
        end
    end

    if delete_temp_files
        files2delete = allfiles(paths[Symbol(mission)].geotile; fn_startswith="jl_")

        if isempty(files2delete)
            continue
        else
            rm.(files2delete)
        end
    end
end
       


