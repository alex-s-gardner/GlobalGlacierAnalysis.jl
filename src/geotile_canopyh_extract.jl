using Altim, GeoArrays, Dates, Arrow, DataFrames

# canopy_height

## Download source data at teminal
#= 
;cd /mnt/devon-r0/shared_data/canopy_height/

;aria2c -x10 https://share.phys.ethz.ch/~pf/nlangdata/ETH_GlobalCanopyHeight_10m_2020_version1/ETH_GlobalCanopyHeight_10m_2020_version1.tar.gz 

;tar -xvzf ETH_GlobalCanopyHeight_10m_2020_version1.tar.gz
=#

force_remake = false;

project_id = :v01;
geotile_width = 2;
#domain = :landice;

paths = project_paths(project_id = project_id);
geotiles = project_geotiles(; geotile_width = geotile_width);
pathlocal = setpaths();

binned_folder = analysis_paths(; geotile_width).binned
mask = :glacier
gt_file = joinpath(binned_folder, "geotile_$(mask)_hyps.arrow");
geotiles = DataFrame(Arrow.Table(gt_file));
geotiles.extent = Extent.(getindex.(geotiles.extent, 1));

# filter geotiles
geotiles = geotiles[geotiles.glacier_frac .> 0.0,:];

#=
mission = :icesat
flist = allfiles(paths[mission].geotile, fn_endswith=".canopyh");
last_modified = Dates.unix2datetime.(mtime.(joinpath.(paths[mission].geotile, flist)));
ind = findall(maximum(last_modified) .== last_modified);
println("last file written = $(flist[ind])")

ext = Extent(X=(-180, 180), Y=(71, 90));
geotiles = geotile_subset!(geotiles, ext);
=#

missions = [keys(paths)...];
ga = GeoArrays.read(pathlocal.canopyheight_10m_v1, masked=false);
nodatavalue = 255;


# remove all older files
#for mission in missions
#    rm.(allfiles(paths[mission].geotile, fn_endswith=".canopyh"))
#end

Altim.geotile_pointextract(geotiles, [paths[mission].geotile for mission in missions], ga; var_name = :canopyh, job_ids = missions, nodatavalue = nodatavalue, force_remake = force_remake)