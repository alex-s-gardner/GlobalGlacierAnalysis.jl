"""
    geotile_dem_extract.jl

Extract DEM data for altimetry measurements in geotiles.

This script:
1. Extracts elevation data from multiple DEMs (REMA, COP30, ArcticDEM, NASADEM)
2. Calculates slope and curvature metrics at altimetry measurement locations
3. Processes data for multiple altimetry missions (ICESat-2, ICESat, GEDI, Hugonnet)
4. Optionally processes unfiltered Hugonnet data separately
5. Organizes extracted data into consistent geotile structure

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

    slope = true;
    curvature = true;
    dems2extract = [:rema_v2_10m, :cop30_v2, :arcticdem_v4_10m, :nasadem_v1]

    # Initialize: paths, products, geotiles
    paths = GGA.project_paths(; project_id);
    products = GGA.project_products(; project_id)
    geotiles = GGA.geotiles_w_mask(geotile_width);

    # Subset: region & mission 
    geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :];
    isa(missions, Tuple) || error("missions must be a tuple... maybe you forgot a trailing comma for single-element tuples?")
    products = getindex(products, missions)

    # Execute: extract dems
    GGA.geotile_extract_dem(products, dems2extract, geotiles, paths; slope, curvature, force_remake)

    # include hugonnet unfiltered
    if hugonnet_unfiltered && (:hugonnet in missions)
        paths = GGA.update_geotile_path(paths; mission = :hugonnet, path_replace ="/2deg" => "/2deg_unfiltered")
        GGA.geotile_extract_dem(products[(:hugonnet,)], dems2extract, geotiles, paths[(:hugonnet,)]; slope, curvature, force_remake)
    end
end
