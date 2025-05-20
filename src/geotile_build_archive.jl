"""
    geotile_build_archive.jl

Process and organize altimetry data into geotiles for glacier analysis.

This script:
1. Searches for altimetry granules from specified missions (ICESat-2, ICESat, GEDI, Hugonnet)
2. Downloads granules from remote repositories
3. Organizes data into geotiles with consistent spatial dimensions
4. Optimizes processing by sorting granules by location

The script handles permission issues with .netrc files and includes special handling for
polar regions which require more processing time.
"""

# IF IT'S NOT WORKING IT MIGHT BE THE PERMISSIONS ON THE .NETRC, fix it with 
# ;chmod 600 ~/.netrc

begin
    # add packages
    import GlobalGlacierAnalysis as GGA
    using Statistics
    using Dates


    # Parameters: user defined 
    project_id = :v01;
    geotile_width = 2;
    domain = :glacier; # :glacier -or- :landice
    missions = (:icesat2,); # (:icesat2, :icesat, :gedi, :hugonnet)

    rebuild_geotiles_dataframe = false;

    # Initialize: paths, products, geotiles
    paths = GGA.project_paths(; project_id);
    products = GGA.project_products(; project_id);
    geotiles = GGA.geotiles_w_mask(geotile_width);

    # Subset: region & mission 
    geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :];
    isa(missions, Tuple) || error("missions must be a tuple... maybe you forgot a trailing comma for single-element tuples?")
    products = getindex(products, missions)

    # Execute: find granules, download granules, build geotiles
    for product in products
    #product = first(products)
        if false
            # make directly if it doesn't exist
            if !isdir(paths[product.mission].raw_data)
                mkpath(paths[product.mission].raw_data)
            end

            if !isdir(paths[product.mission].geotile)
                mkpath(paths[product.mission].geotile)
            end

            ## 'find' can run in parallel.. therfore run first [GLAH06 = 30 min from scratch]
            # do not use kward `after` as it will cuase downstream issues as earlier data will be excluded
            GGA.geotile_search_granules(geotiles, product.mission, product.name, product.version, paths[product.mission].granules_remote; rebuild_dataframe=rebuild_geotiles_dataframe)

            # load remote granule list
            geotile_granules = GGA.granules_load(paths[product.mission].granules_remote, product.mission; geotiles = geotiles)

            # download granules from list [ATL06 17 min, no data]

            # download can sometimes be unhappy when run with more than one downloadstream
            # NOTE: I'm not sure why download keeps failing ... wrapping in a while loop as a bandaid
            foo = true
            while foo
                try
                    GGA.geotile_download_granules!(geotile_granules, product.mission, paths[product.mission].raw_data, 
                        paths[product.mission].granules_local; threads = false, aria2c = true, 
                        rebuild_dataframe = rebuild_geotiles_dataframe, downloadstreams = 8)

                    foo = false;
                catch e
                    println(e)
                end
            end
        end

        # load local granule list [most load time used adding back granule_type]
        geotile_granules = GGA.granules_load(paths[product.mission].granules_local, product.mission; geotiles = geotiles)


        # sort by longitude for imporved speed (less cashing as granules run north-south)
        geotile_granules[!, :longitude] = mean.(getindex.(geotile_granules.extent, :X))
        geotile_granules[!, :latitude] = mean.(getindex.(geotile_granules.extent, :Y))    
        sort!(geotile_granules, [:longitude, :latitude])
        # split at equator as this is where most granules are custom
        ind = geotile_granules[:,:latitude] .>0
        geotile_granules = vcat(geotile_granules[ind, :], geotile_granules[.!ind, :])

        # extract data and save locally 
        # build geotiles for each mission 

        # place < -65 last as they take FOREVER
        ind = (geotile_granules[:, :latitude] .> -91) .& (geotile_granules[:, :latitude] .< -84.9)
        geotile_granules = vcat(geotile_granules[ind, :], geotile_granules[.!ind, :])

        GGA.geotile_build(geotile_granules, paths[product.mission].geotile; warnings=false)
    end
end