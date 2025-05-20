"""
    geotile_mask_extract.jl

Extract mask data for altimetry measurements in geotiles.

This script:
1. Extracts various mask types (floatingice, glacierice, inlandwater, etc.) for altimetry points
2. Processes data for multiple altimetry missions (ICESat-2, ICESat, GEDI, Hugonnet)
3. Optionally processes unfiltered Hugonnet data separately
4. Organizes extracted mask data into consistent geotile structure

The script supports both landice and glacier domains and can be configured to
force remake existing extractions if needed.
"""

begin
    import GlobalGlacierAnalysis as GGA

    # Parameters: user defined 
    force_remake = false
    project_id = :v01;
    geotile_width = 2;
    domain = :landice; # :glacier -or- :landice
    missions = (:icesat2, :icesat, :gedi, :hugonnet,); # (:icesat2, :icesat, :gedi, :hugonnet)
    hugonnet_unfiltered = true;
    vars2extract = [:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean];

    # Initialize: paths, products, geotiles
    paths = GGA.project_paths(; project_id)
    products = GGA.project_products(; project_id)
    geotiles = GGA.geotiles_w_mask(geotile_width);

    # Subset: region & mission 
    geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :];
    isa(missions, Tuple) || error("missions must be a tuple... maybe you forgot a trailing comma for single-element tuples?")
    products = getindex(products, missions)

    # Execute: extract masks
    for product in products
        GGA.geotile_extract_mask(geotiles, paths[product.mission].geotile; vars=vars2extract, job_id=product.mission, force_remake=force_remake)
    end

    # include hugonnet unfiltered
    if hugonnet_unfiltered && (:hugonnet in missions)
        paths = GGA.update_geotile_path(paths; mission=:hugonnet, path_replace="/2deg" => "/2deg_unfiltered")
        GGA.geotile_extract_mask(geotiles, paths[:hugonnet].geotile; vars=vars2extract, job_id=:hugonnet, force_remake=force_remake)
    end
end