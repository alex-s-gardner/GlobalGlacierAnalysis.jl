# This script processes altimetry data from ICESat, ICESat-2, and GEDI missions into geotiles.
#
# The script:
# 1. Sets up project configuration:
#    - Defines project ID and geotile width (2 degrees)
#    - Loads project paths and product definitions
#    - Filters geotiles to only include those with land ice
#
# 2. For each altimetry product (ICESat/ICESat-2/GEDI):
#    - Creates required directories for raw data and geotiles
#    - Searches for granules matching geotile extents using parallel 'find'
#    - Downloads granules from remote source using aria2c with retry logic
#    - Loads and processes granule data
#    - Sorts granules by longitude/latitude for improved processing speed
#    - Builds geotiles from the granule data
#
# Processing times (on old RAID @ 100 MB/s):
# - GEDI: ~4 days
# - ICESat-2: ~1 week  
# - ICESat: ~3 hours


# user credentials for downloading files
# using SpaceLiDAR
# netrc_credentials = (username = "alex.s.gardner", password = "0c!Y7Pfcg35a")  # replace with your credentials
# SpaceLiDAR.netrc!(netrc_credentials...)

# IF IT'S NOT WORKING IT MIGHT BE THE PERMISSIONS ON THE .NETRC, fix it with 
# ;chmod 600 ~/.netrc


# add packages
using Altim
using Statistics
using Dates


# Parameters: user defined 
force_remake = true
project_id = :v01;
geotile_width = 2;
domain = :glacier; # :glacier -or- :landice
missions = (:icesat2,); # (:icesat2, :icesat, :gedi, :hugonnet)

after = nothing # only search for after this date
rebuild_geotiles_dataframe = true;

# Initialize: paths, products, geotiles
paths = project_paths(; project_id);
products = project_products(; project_id);
geotiles = Altim.geotiles_w_mask(geotile_width);

# Subset: region & mission 
geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :];
isa(missions, Tuple) || error("missions must be a tuple... maybe you forgot a trailing comma for single-element tuples?")
products = getindex(products, missions)

# Execute: find granules, download granules, build geotiles
for product in products
#product = first(products)
    
    # make directly if it doesn't exist
    if !isdir(paths[product.mission].raw_data)
        mkpath(paths[product.mission].raw_data)
    end

    if !isdir(paths[product.mission].geotile)
        mkpath(paths[product.mission].geotile)
    end

    ## 'find' can run in parallel.. therfore run first [GLAH06 = 30 min from scratch]
    geotile_search_granules(geotiles, product.mission, product.name, product.version, paths[product.mission].granules_remote; rebuild_dataframe=rebuild_geotiles_dataframe, after =after)

    # load remote granule list
    geotile_granules = granules_load(paths[product.mission].granules_remote, product.mission; geotiles = geotiles)

    # download granules from list [ATL06 17 min, no data]

    # download can sometimes be unhappy when run with more than one downloadstream
    # NOTE: I'm not sure why download keeps failing ... wrapping in a while loop as a bandaid
    foo = true
    while foo
        try
            geotile_download_granules!(geotile_granules, product.mission, paths[product.mission].raw_data, 
                paths[product.mission].granules_local; threads = false, aria2c = true, 
                rebuild_dataframe = rebuild_geotiles_dataframe, downloadstreams = 8)

            foo = false;
        catch e
            println(e)
        end
    end

    # load local granule list [most load time used adding back granule_type]
    geotile_granules = granules_load(paths[product.mission].granules_local, product.mission; geotiles = geotiles)


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

    geotile_build(geotile_granules, paths[product.mission].geotile; warnings=false)
end