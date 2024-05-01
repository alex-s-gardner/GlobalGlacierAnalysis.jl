using Altim

force_remake = true;

project_id = :v01
geotile_width = 2
domain = :all

paths = project_paths(project_id = project_id)
geotiles = project_geotiles(; geotile_width = geotile_width, domain = domain)

# geotiles = geotiles[ind,:]
dems = [:cop30_v2, :nasadem_v1, :rema_v2_10m, :arcticdem_v3_10m];
missions = [:icesat, :icesat2, :gedi]

valid_count_thresh = 50;

#
#ind = geotiles.id .== "lat[-74-72]lon[-114-112]";
#geotiles = geotiles[ind,:];

for mission in missions
    for dem in [dems[end]]
        geotile_track_offset(geotiles, paths[missions].geotile, dem; valid_count_thresh = valid_count_thresh, force_remake = force_remake)
    end
end

