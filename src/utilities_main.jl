"""
    geotile_build_archive(; force_remake=false, project_id=:v01, geotile_width=2, 
                         domain=:glacier, missions=(:icesat2,), single_geotile_test=nothing)

Build satellite altimetry archives organized into geotiles for glacier analysis.

This function processes altimetry data from multiple missions (ICESat-2, ICESat, GEDI, Hugonnet)
by searching for granules, downloading them from remote repositories, and organizing them
into geotiles with consistent spatial dimensions. It optimizes processing by sorting granules
by location and includes special handling for polar regions.

# Arguments
- `force_remake`: Force recreation of existing archives (default: false)
- `project_id`: Project identifier (default: :v01)
- `geotile_width`: Width of geotiles in degrees (default: 2)
- `domain`: Processing domain, either :glacier or :landice (default: :glacier)
- `missions`: Tuple of altimetry missions to process (default: (:icesat2,))
- `single_geotile_test`: Single geotile ID for testing (default: nothing)

# Returns
Nothing. Processed altimetry data is organized into geotiles and saved locally.

# Notes
- Downloads can fail with multiple download streams; wrapped in retry loop
- Polar regions (< -65° latitude) are processed last as they require more time
- Granules are sorted by longitude for improved processing speed
- When `single_geotile_test` is specified, only processes that specific geotile for testing
"""

function geotile_build_archive(;
    force_remake = false,
    project_id = :v01,
    geotile_width = 2,
    domain = :glacier, # :glacier -or- :landice
    missions = (:icesat2,), # (:icesat2, :icesat, :gedi, :hugonnet)
    single_geotile_test = nothing,
    )

    rebuild_geotiles_dataframe = force_remake

    # Initialize: paths, products, geotiles
    paths = project_paths(; project_id);
    products = project_products(; project_id);
    geotiles = geotiles_w_mask(geotile_width; remake=false);

    # Subset: region & mission 
    if isnothing(single_geotile_test)
        geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :]
    else
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)] !!!!!!!!!!!!!!"
        geotiles = geotiles[geotiles.id .== single_geotile_test, :]
    end

    isa(missions, Tuple) || error("missions must be a tuple... maybe you forgot a trailing comma for single-element tuples?")
    products = getindex(products, missions)

    # Execute: find granules, download granules, build geotiles
    for product in products
        if product.mission == :hugonnet
            continue
        end
        #product = first(products)

        # make directly if it doesn't exist
        if !isdir(paths[product.mission].raw_data)
            mkpath(paths[product.mission].raw_data)
        end

        if !isdir(paths[product.mission].geotile)
            mkpath(paths[product.mission].geotile)
        end

        ## 'find' can run in parallel.. therfore run first [GLAH06 = 30 min from scratch]
        # do not use kward `after` as it will cuase downstream issues as earlier data will be excluded
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
end


"""
    geotile_build_hugonnet(; hugonnet_datasets=["unfiltered", "filtered"], force_remake=false, 
                          project_id=:v01, geotile_width=2, domain=:glacier, single_geotile_test=nothing)

Build geotiles from Hugonnet glacier elevation change data.

# Arguments
- `hugonnet_datasets`: Dataset types to process ("filtered" or "unfiltered")
- `force_remake`: Force recreation of existing files
- `project_id`: Project identifier
- `geotile_width`: Width of geotiles in degrees
- `domain`: Processing domain (:glacier or :landice)
- `single_geotile_test`: Single geotile ID for testing (default: nothing)

# Description
Processes Hugonnet glacier elevation change data into geotiles. Handles both filtered and unfiltered datasets,
creating separate output directories for each. Filters geotiles based on domain coverage and builds individual
geotile files from the Hugonnet data stacks.
"""
function geotile_build_hugonnet(;
    # Parameters: user defined 
    hugonnet_datasets = ["unfiltered", "filtered"], # "filtered" or "unfiltered"
    force_remake = false,
    project_id = :v01,
    geotile_width = 2,
    domain = :glacier, # :glacier -or- :landice
    single_geotile_test = nothing,
    )

    # Initialize: paths, products, geotiles
    paths = project_paths(; project_id);
    geotiles = geotiles_w_mask(geotile_width);

    # Subset: region & mission 
    if isnothing(single_geotile_test)
        geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :];
    else
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)] !!!!!!!!!!!!!!"
        geotiles = geotiles[geotiles.id.==single_geotile_test, :]
    end

    for hugonnet_dataset in hugonnet_datasets
        # Execute: build geotiles from Hugonnet data
        if hugonnet_dataset .== "filtered"
            hugonnet_dataset = "filtered"
            raw_data = replace(paths.hugonnet.raw_data, "raw" => "raw_filtered")
            geotile_dir = paths.hugonnet.geotile
        else
            hugonnet_dataset = "unfiltered"
            raw_data = paths.hugonnet.raw_data
            geotile_dir = replace(paths.hugonnet.geotile, "2deg" => "2deg_unfiltered")
        end

        # get hstack catalogue
        hstacks = hstack_catalogue(raw_data; force_remake)

        # make directly if it doesn't exist
        if !isdir(geotile_dir)
            mkpath(geotile_dir)
        end

        # build geotiles
        printstyled("building Hugonnet geotiles\n"; color=:blue, bold=true)
        @warn "!!! hugonnet quality flag currently excludes all ArcticDEM data !!!"

        @showprogress dt = 10 desc = "building hugonnet [$(hugonnet_dataset)] geotiles..." Threads.@threads for geotile in eachrow(geotiles)
            geotile_build_hugonnet(geotile, geotile_dir, hstacks; force_remake)
        end
    end
end

"""
    geotile_ancillary_check(; missions=(:icesat2, :icesat, :gedi, :hugonnet))

Validate and clean ancillary data files associated with altimetry measurements.

This function checks that all altimetry, mask, and DEM files exist for each geotile and verifies
that ancillary files have the same number of rows as their corresponding altimetry files.
Any ancillary files with mismatched row counts are deleted to prevent data inconsistencies.

# Arguments
- `missions`: Tuple of altimetry missions to process (default: (:icesat2, :icesat, :gedi, :hugonnet))

# Returns
Nothing. Invalid ancillary files are deleted from the filesystem.

# Description
Processes multiple altimetry missions and various ancillary data types (masks, DEMs, canopy height).
For Hugonnet data, both filtered and unfiltered versions are checked. Files with row count
mismatches are automatically removed to maintain data integrity.
"""
function geotile_ancillary_check(;
    project_id = :v01,
    missions = (:icesat2, :icesat, :gedi, :hugonnet),
    )

    # Initialize: paths, products, geotiles
    paths = project_paths(; project_id);
    mission_geotile_folders = [paths[mission].geotile for mission in missions]

    # add hugonnet unfiltered to the list
    if :hugonnet in missions
        mission_geotile_folders = vcat(mission_geotile_folders, replace(paths[:hugonnet].geotile, "/2deg" => "/2deg_unfiltered"))
    end

    # Execute: check length of all files
    for mission_geotile_folder in mission_geotile_folders

        mission_dir_files = allfiles(mission_geotile_folder)
        paths2altim = filter(x -> occursin(".arrow", x), mission_dir_files)

        if isempty(paths2altim)
            @warn "no altimetry files found in: $mission_geotile_folder"
            continue
        end

        printstyled("$(mission_geotile_folder)\n", color=:blue)
        @showprogress desc = "checking that number of rows in altimetry == ancillary:$(mission_geotile_folder)..." Threads.@threads for path2altim in paths2altim
        
            geotile_id = replace(splitpath(path2altim)[end], ".arrow" => "")

            paths2ancillary = filter(x -> (occursin(geotile_id, x) && !occursin(".arrow", x)), mission_dir_files)

            n_rows = length(Arrow.Table(path2altim)[1])

            for path2ancillary in paths2ancillary
                n_rows_ancillary = length(Arrow.Table(path2ancillary)[1])
                if n_rows_ancillary != n_rows
                    @warn "number of rows in altimetry  != ancillary, deleting: $path2ancillary"
                    rm(path2ancillary)
                end
            end
        end
    end
end


"""
    geotile_dem_extract(; force_remake=false, project_id=:v01, geotile_width=2, 
                       domain=:landice, missions=(:icesat2, :icesat, :gedi, :hugonnet,), 
                       hugonnet_unfiltered=true, slope=true, curvature=true, 
                       dems2extract=[:rema_v2_10m, :cop30_v2, :arcticdem_v4_10m, :nasadem_v1])

Extract DEM data for altimetry measurements in geotiles.

This function processes multiple altimetry missions and extracts elevation data from various DEMs,
calculating slope and curvature metrics at measurement locations. It supports both glacier and 
land ice domains and can optionally process unfiltered Hugonnet data separately.

# Arguments
- `force_remake`: Force recreation of existing extractions (default: false)
- `project_id`: Project identifier (default: :v01)
- `geotile_width`: Width of geotiles in degrees (default: 2)
- `domain`: Processing domain, either :glacier or :landice (default: :landice)
- `missions`: Tuple of altimetry missions to process (default: (:icesat2, :icesat, :gedi, :hugonnet,))
- `hugonnet_unfiltered`: Process unfiltered Hugonnet data separately (default: true)
- `slope`: Calculate slope metrics (default: true)
- `curvature`: Calculate curvature metrics (default: true)
- `dems2extract`: Array of DEMs to extract (default: [:rema_v2_10m, :cop30_v2, :arcticdem_v4_10m, :nasadem_v1])
- `single_geotile_test`: Test with single geotile (default: nothing)

# Returns
Nothing. Extracted DEM data is saved to geotile structure.
"""
function geotile_dem_extract(;
    # Parameters: user defined 
    force_remake=false,
    project_id=:v01,
    geotile_width=2,
    domain=:landice, # :glacier -or- :landice
    missions=(:icesat2, :icesat, :gedi, :hugonnet,), # (:icesat2, :icesat, :gedi, :hugonnet)
    hugonnet_unfiltered=true, 
    slope=true,
    curvature=true,
    dems2extract=[:rema_v2_10m, :cop30_v2, :arcticdem_v4_10m, :nasadem_v1],
    single_geotile_test = nothing,
    )

    # Initialize: paths, products, geotiles
    paths = project_paths(; project_id)
    products = project_products(; project_id)
    geotiles = geotiles_w_mask(geotile_width)

    # Subset: region & mission 
    if isnothing(single_geotile_test)
        geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :]
    else
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)] !!!!!!!!!!!!!!"
        geotiles = geotiles[geotiles.id.==single_geotile_test, :]
    end

    isa(missions, Tuple) || error("missions must be a tuple... maybe you forgot a trailing comma for single-element tuples?")
    products = getindex(products, missions)

    # Execute: extract dems
    geotile_extract_dem(products, dems2extract, geotiles, paths; slope, curvature, force_remake)

    # include hugonnet unfiltered
    if hugonnet_unfiltered && (:hugonnet in missions)
        paths = update_geotile_path(paths; mission=:hugonnet, path_replace="/2deg" => "/2deg_unfiltered")
        geotile_extract_dem(products[(:hugonnet,)], dems2extract, geotiles, paths[(:hugonnet,)]; slope, curvature, force_remake)
    end
end

"""
    geotile_mask_extract(; force_remake=false, project_id=:v01, geotile_width=2, 
                        domain=:glacier, missions=(:icesat2, :icesat, :gedi, :hugonnet,),
                        hugonnet_unfiltered=true, 
                        masks2extract=[:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean],
                        masks2extract_highres=[:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km],
                        single_geotile_test=nothing)

Extract mask data for altimetry measurements in geotiles.

This function processes various mask types (floatingice, glacierice, inlandwater, etc.) for altimetry points
across multiple missions and can optionally process unfiltered Hugonnet data separately.

# Arguments
- `force_remake`: Force recreation of existing extractions (default: false)
- `project_id`: Project identifier (default: :v01)
- `geotile_width`: Width of geotiles in degrees (default: 2)
- `domain`: Processing domain, either :glacier or :landice (default: :glacier)
- `missions`: Tuple of altimetry missions to process (default: (:icesat2, :icesat, :gedi, :hugonnet,))
- `hugonnet_unfiltered`: Process unfiltered Hugonnet data separately (default: true)
- `masks2extract`: Array of mask variables to extract (default: [:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean])
- `masks2extract_highres`: Array of high-resolution mask variables to extract (default: [:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km])
- `single_geotile_test`: Single geotile ID for testing (default: nothing)

# Returns
Nothing. Extracted mask data is saved to geotile structure.
"""
function geotile_mask_extract(;
    # Parameters: user defined 
    force_remake = false,
    project_id = :v01,
    geotile_width = 2,
    domain = :glacier, # :glacier -or- :landice
    missions = (:icesat2, :icesat, :gedi, :hugonnet,), # (:icesat2, :icesat, :gedi, :hugonnet)
    hugonnet_unfiltered = true,
    masks2extract=[:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean],
    masks2extract_highres = [:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km],
    single_geotile_test = nothing,
    )

    # Initialize: paths, products, geotiles
    paths = project_paths(; project_id)
    products = project_products(; project_id)
    geotiles = geotiles_w_mask(geotile_width);

    # Subset: region & mission 
    if isnothing(single_geotile_test)
        geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :];
    else
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)] !!!!!!!!!!!!!!"
        geotiles = geotiles[geotiles.id.==single_geotile_test, :]
    end

    isa(missions, Tuple) || error("missions must be a tuple... maybe you forgot a trailing comma for single-element tuples?")
    products = getindex(products, missions)

    # Execute: extract masks
    for product in products
        geotile_extract_mask(geotiles, paths[product.mission].geotile; masks=masks2extract, masks_highres=masks2extract_highres, job_id=product.mission, force_remake)
    end

    # include hugonnet unfiltered
    if hugonnet_unfiltered && (:hugonnet in missions)
        paths = update_geotile_path(paths; mission=:hugonnet, path_replace="/2deg" => "/2deg_unfiltered")
        geotile_extract_mask(geotiles, paths[:hugonnet].geotile; masks=masks2extract, masks_highres=masks2extract_highres, job_id=:hugonnet, force_remake)
    end
end

"""
    geotile_canopyh_extract(; force_remake=false, project_id=:v01, geotile_width=2, 
                           domain=:landice, missions=(:icesat2, :icesat, :gedi, :hugonnet,), 
                           hugonnet_unfiltered=true, single_geotile_test=nothing)

Extract canopy height data for altimetry measurements in geotiles.

This function processes canopy height data from ETH's Global Canopy Height 10m 2020 dataset for altimetry points
across multiple missions and can optionally process unfiltered Hugonnet data separately.

# Arguments
- `force_remake`: Force recreation of existing extractions (default: false)
- `project_id`: Project identifier (default: :v01)
- `geotile_width`: Width of geotiles in degrees (default: 2)
- `domain`: Processing domain, either :glacier or :landice (default: :landice)
- `missions`: Tuple of altimetry missions to process (default: (:icesat2, :icesat, :gedi, :hugonnet,))
- `hugonnet_unfiltered`: Process unfiltered Hugonnet data separately (default: true)
- `single_geotile_test`: Single geotile ID for testing (default: nothing)

# Returns
Nothing. Extracted canopy height data is saved to geotile structure.

# Notes
- Requires ETH Global Canopy Height 10m 2020 dataset to be downloaded and extracted
- Uses nodatavalue of 255 for canopy height data
- Processes both filtered and unfiltered Hugonnet data when specified
"""
function geotile_canopyh_extract(;
    # Parameters: user defined 
    force_remake = false,
    project_id = :v01,
    geotile_width = 2,
    domain = :landice, # :glacier -or- :landice
    missions = (:icesat2, :icesat, :gedi, :hugonnet,), # (:icesat2, :icesat, :gedi, :hugonnet)
    hugonnet_unfiltered = true,
    single_geotile_test = nothing,
    )

    ## Download canopy height data in teminal
    #= 
    ;cd /mnt/devon-r0/shared_data/canopy_height/
    ;aria2c -x10 https://share.phys.ethz.ch/~pf/nlangdata/ETH_GlobalCanopyHeight_10m_2020_version1/ETH_GlobalCanopyHeight_10m_2020_version1.tar.gz 
    ;tar -xvzf ETH_GlobalCanopyHeight_10m_2020_version1.tar.gz
    =#

    # Initialize: paths, products, geotiles
    paths = project_paths(; project_id);
    products = project_products(; project_id);
    geotiles = geotiles_w_mask(geotile_width);

    if isnothing(single_geotile_test)
        # Subset: region & mission 
        geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :]
    else
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)] !!!!!!!!!!!!!!"
        geotiles = geotiles[geotiles.id.==single_geotile_test, :]
    end

    isa(missions, Tuple) || error("missions must be a tuple... maybe you forgot a trailing comma for single-element tuples?")
    products = getindex(products, missions)

    # Execute: extract canopy height
    ga = GeoArrays.read(setpaths().canopyheight_10m_v1, masked=false)
    nodatavalue = 255;
    geotile_pointextract(geotiles, [paths[mission].geotile for mission in missions], ga; var_name = :canopyh, job_ids = [missions...], nodatavalue = nodatavalue, force_remake = force_remake)

    # include hugonnet unfiltered
    if hugonnet_unfiltered && (:hugonnet in missions)
        paths = update_geotile_path(paths; mission=:hugonnet, path_replace="/2deg" => "/2deg_unfiltered")
        geotile_pointextract(geotiles, paths[:hugonnet].geotile, ga; var_name=:canopyh, job_ids=[:hugonnet,], nodatavalue=nodatavalue, force_remake=force_remake)
    end
end


"""
    geotile_hyps_extract(; force_remake=false, geotile_width=2, raster_file=:cop30_v2, 
                        domain=:landice, masks=[:glacier, :land, :glacier_rgi7, :glacier_b1km, :glacier_b10km])

Generate hypsometric area distributions for geotiles by calculating area distributions by elevation for different land cover masks.

This function processes elevation data to create area distributions binned by elevation for specified land cover types.
It uses a DEM raster (default: COP30) and shapefile masks to calculate area in km² for each elevation bin.
Results are saved as Arrow files in the analysis/binned directory for further hypsometric analysis.

# Arguments
- `force_remake`: Force recreation of existing files (default: false)
- `geotile_width`: Width of geotiles in degrees (default: 2)
- `raster_file`: DEM raster to use for elevation data (default: :cop30_v2)
- `domain`: Processing domain, either :glacier or :landice (default: :landice)
- `masks`: Array of land cover mask types to process (default: [:glacier, :land, :glacier_rgi7, :glacier_b1km, :glacier_b10km])

# Returns
Nothing. Processed hypsometric data is saved as Arrow files in the analysis/binned directory.

# Notes
- Uses project height bins for consistent elevation binning across analyses
- Processes geotiles with domain coverage > 0
- Supports multithreaded processing for efficiency
- Special handling for land mask (uses water shapefile with invert=true)
- Excludes glacier areas when processing land mask
"""
function geotile_hyps_extract(;
    force_remake = false,
    geotile_width = 2,
    raster_file = :cop30_v2, #, :nasadem_v1, :rema_v2_10m, :arcticdem_v3_10m]
    domain = :landice, # :glacier -or- :landice
    masks = [:glacier, :land, :glacier_rgi7, :glacier_b1km, :glacier_b10km],
    )

    Threads.@threads for mask in masks

        runid = "geotile_$(mask)_hyps"
        binned_folder = analysis_paths(; geotile_width).binned
        out_file = joinpath(binned_folder, "$(runid).arrow")

        if !isfile(out_file) || force_remake
            height_range, height_center = project_height_bins()
            Δh = abs(height_center[2] - height_center[1])

            excludemask_flag = false

            var_name = Symbol("$(mask)_area_km2")
            geotiles = geotiles_w_mask(geotile_width)
            fn_raster = pathlocal[raster_file]
            ras = Raster(fn_raster; lazy=true)

            if mask == :land
                shp = Symbol("$(:water)_shp")
                invert = true
            else
                shp = Symbol("$(mask)_shp")
                invert = false
            end

            geotiles[!, var_name] = [zeros(size(height_center)) for r in 1:nrow(geotiles)]

            fn_shp = GGA.pathlocal[shp]

            if mask == :land
                shp = Symbol("$(:glacier)_shp")
                fn_shp_ex = pathlocal[shp]
            else
                fn_shp_ex = nothing
                excludefeature = nothing
            end

            feature = Shapefile.Handle(fn_shp)
            if !isnothing(fn_shp_ex)
                excludefeature = Shapefile.Handle(fn_shp_ex)
            else
                excludefeature = nothing
            end

            # using Threads here does not improve performance
            for geotile in eachrow(geotiles)[geotiles[!, "$(domain)_frac"].>0]

                GGA.geotile_binarea!(geotile, ras, feature, height_range; invert, excludefeature, var_name)
            end

            Arrow.write(out_file, select!(geotiles, Not(:geometry))::DataFrame)

            # file can not be written with geometry column
            Arrow.write(out_file, select!(geotiles, Not(:geometry))::DataFrame)
        end
    end
end