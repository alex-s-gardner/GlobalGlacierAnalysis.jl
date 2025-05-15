"""
    geotile_build_hugonnet.jl

Process and organize Hugonnet glacier elevation change data into geotiles.

This script:
1. Loads Hugonnet glacier elevation change data from raw files
2. Organizes data into standardized geotiles with consistent spatial dimensions
3. Handles both old and new Hugonnet data formats
4. Creates output directories if they don't exist
5. Processes each geotile individually for efficient parallel processing

The script is configured to work with landice or glacier domains and can be
forced to rebuild existing geotiles if needed.
"""

begin
    using Altim
    using Extents


    # Parameters: user defined 
    force_remake = false
    project_id = :v01;
    geotile_width = 2;
    domain = :landice; # :glacier -or- :landice

    # Initialize: paths, products, geotiles
    paths = project_paths(; project_id);
    products = project_products(; project_id);
    geotiles = Altim.geotiles_w_mask(geotile_width);

    # Subset: region & mission 
    geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :];
    isa(missions, Tuple) || error("missions must be a tuple... maybe you forgot a trailing comma for single-element tuples?")
    products = getindex(products, missions)


    # Execute: build geotiles from Hugonnet data

    # make directly if it doesn't exist
    if !isdir(paths.hugonnet.geotile)
        mkpath(paths.hugonnet.geotile)
    end

    # get hstack catalogue
    hstacks = hstack_catalogue(paths.hugonnet.raw_data; update_catalogue=force_remake)

    if contains(hstacks.path[1], "prefilt")
        old_format = false
    else
        old_format = true
    end

    # build geotiles
    begin
        printstyled("building Hugonnet geotiles\n"; color=:blue, bold=true)
        for geotile in eachrow(geotiles)
            geotile_build_hugonnet(geotile, paths.hugonnet.geotile, hstacks; force_remake, old_format)
        end
    end
end