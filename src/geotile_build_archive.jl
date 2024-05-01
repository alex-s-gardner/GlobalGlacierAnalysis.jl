# user credentials for downloading files
# using SpaceLiDAR
# netrc_credentials = (username = "alex.s.gardner", password = "0c!Y7Pfcg35a")  # replace with your credentials
# SpaceLiDAR.netrc!(netrc_credentials...)

# IF IT'S NOT WORKING IT MIGHT BE THE PERMISSIONS ON THE .NETRC, fix it with 
# ;chmod 600 ~/.netrc

# add packages
using Altim
using Statistics

# set paths
rebuild_geotiles_dataframe = false;

project_id = :v01;
geotile_width = 2;

paths = project_paths(project_id = project_id);
products = project_products(project_id = project_id);
geotiles = Altim.geotiles_w_mask(geotile_width);

geotiles = geotiles[geotiles.landice_frac .> 0, :];

## -------------------------------------
#products = (icesat=products.icesat,);
products = (icesat2=products.icesat2,);
#ext = Extent(X=(-126.9, -126.1), Y=(51.1, 51.8));
#geotiles = project_geotiles(; geotile_width=geotile_width, domain=domain, extent=ext);
## -------------------------------------

# on old RAID (100 MB/s)
# 4 days for GEDI
# 1-week day for ICESat2
# 3 hours for ICESat

for product in products
    
    #=
    # make directly if it doesn't exist
    if !isdir(paths[product.mission].raw_data)
        mkpath(paths[product.mission].raw_data)
    end

    if !isdir(paths[product.mission].geotile)
        mkpath(paths[product.mission].geotile)
    end

    ## 'find' can run in parallel.. therfore run first [GLAH06 = 30 min from scratch]
    geotile_search_granules(geotiles, product.mission, product.name, product.version, paths[product.mission].granules_remote; rebuild_dataframe=rebuild_geotiles_dataframe)

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
    =#

    # load local granule list [most load time used adding back granule_type]
    geotile_granules = granules_load(paths[product.mission].granules_local, product.mission; geotiles = geotiles)

    #=
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
    =#

    geotile_build(geotile_granules, paths[product.mission].geotile; warnings=false)
end