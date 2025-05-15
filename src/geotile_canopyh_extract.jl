"""
    geotile_canopyh_extract.jl

Extract canopy height data for altimetry measurements in geotiles.

This script:
1. Loads global canopy height data at 10m resolution
2. Extracts canopy height values at altimetry measurement locations
3. Saves the extracted values as ancillary data files for each geotile
4. Processes multiple altimetry missions (ICESat-2, ICESat, GEDI, Hugonnet)
5. Optionally processes unfiltered Hugonnet data separately

The canopy height data is from ETH's Global Canopy Height 10m 2020 dataset.
"""

begin
    using Altim, GeoArrays, Dates, Arrow, DataFrames

    ## Download canopy height data in teminal
    #= 
    ;cd /mnt/devon-r0/shared_data/canopy_height/
    ;aria2c -x10 https://share.phys.ethz.ch/~pf/nlangdata/ETH_GlobalCanopyHeight_10m_2020_version1/ETH_GlobalCanopyHeight_10m_2020_version1.tar.gz 
    ;tar -xvzf ETH_GlobalCanopyHeight_10m_2020_version1.tar.gz
    =#

    # Parameters: user defined 
    force_remake = false
    project_id = :v01;
    geotile_width = 2;
    domain = :landice; # :glacier -or- :landice
    missions = (:icesat2, :icesat, :gedi, :hugonnet,); # (:icesat2, :icesat, :gedi, :hugonnet)
    hugonnet_unfiltered = true;

    # Initialize: paths, products, geotiles
    paths = project_paths(; project_id);
    products = project_products(; project_id);
    geotiles = Altim.geotiles_w_mask(geotile_width);

    # Subset: region & mission 
    geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :];
    isa(missions, Tuple) || error("missions must be a tuple... maybe you forgot a trailing comma for single-element tuples?")
    products = getindex(products, missions)

    # Execute: extract canopy height
    ga = GeoArrays.read(setpaths().canopyheight_10m_v1, masked=false);
    nodatavalue = 255;
    Altim.geotile_pointextract(geotiles, [paths[mission].geotile for mission in missions], ga; var_name = :canopyh, job_ids = [missions...], nodatavalue = nodatavalue, force_remake = force_remake)

    # include hugonnet unfiltered
    if hugonnet_unfiltered && (:hugonnet in missions)
        paths = Altim.update_geotile_path(paths; mission=:hugonnet, path_replace="/2deg" => "/2deg_unfiltered")
        Altim.geotile_pointextract(geotiles, paths[:hugonnet].geotile, ga; var_name=:canopyh, job_ids=[:hugonnet,], nodatavalue=nodatavalue, force_remake=force_remake)
    end
end
