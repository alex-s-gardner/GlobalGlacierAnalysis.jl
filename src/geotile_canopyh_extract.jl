# This script extracts canopy height data for geotiles:
#
# Configuration:
# - Uses project ID v01 and 2-degree geotile width
# - Processes only geotiles containing land ice
# - Uses ETH Global Canopy Height 10m 2020 dataset
#
# Processing:
# 1. Sets up project paths and loads product definitions
# 2. Filters geotiles to only include those with land ice
# 3. Reads canopy height data from source GeoTIFF
# 4. Extracts canopy height values for each geotile using geotile_pointextract()
#    - Processes all altimetry missions
#    - Only processes if force_remake=true or output doesn't exist
#    - Uses nodata value of 255

begin
    using Altim, GeoArrays, Dates, Arrow, DataFrames

    # canopy_height

    ## Download source data in teminal
    #= 
    ;cd /mnt/devon-r0/shared_data/canopy_height/

    ;aria2c -x10 https://share.phys.ethz.ch/~pf/nlangdata/ETH_GlobalCanopyHeight_10m_2020_version1/ETH_GlobalCanopyHeight_10m_2020_version1.tar.gz 

    ;tar -xvzf ETH_GlobalCanopyHeight_10m_2020_version1.tar.gz
    =#

    force_remake = false;
    project_id = :v01;
    geotile_width = 2;

    paths = project_paths(; project_id);
    pathlocal = setpaths();

    geotiles = Altim.geotiles_w_mask(geotile_width);
    geotiles = geotiles[geotiles.landice_frac .> 0, :];
    products = project_products(; project_id);

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

    Altim.geotile_pointextract(geotiles, [paths[mission].geotile for mission in missions], ga; var_name = :canopyh, job_ids = missions, nodatavalue = nodatavalue, force_remake = force_remake)
end