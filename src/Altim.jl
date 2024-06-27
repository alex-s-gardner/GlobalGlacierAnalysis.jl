module Altim

    using Proj
    using Arrow
    using DataFrames
    using GeoArrays
    using SpaceLiDAR
    using Dates
    using Statistics
    using Geodesy
    using ImageFiltering
    using OffsetArrays
    using HTTP
    using Logging
    using Printf
    using Extents
    using Interpolations
    using StaticArrays
    import GeoFormatTypes as GFT
    using NetCDF
    using BinStatistics
    using FastGeoProjections
    using GeoTiles
    using NearestNeighbors
    using Stencils  
    using Infiltrator
    using Rasters
    using MAT

    #using Optimization
    #using LossFunctions
    #using OptimizationOptimJL
    using MLJ
    using MLJLinearModels
    using LazyGrids
    #using GeometryOps
    using LibGEOS
    using GeoInterface
    using Shapefile
    using FileIO
    #using ProfileView

    #include("utilities.jl")
    include("utilities.jl")
    include("utilities_project.jl")
    include("utilities_hugonnet.jl")
    include("utilities_gemb.jl")
    include("modelfit_tree_fast.jl")
    #include("model_optimize.jl")

    export Extent
    export bin, binnedfiltering, crop!, decimalyear, download!, epsg2epsg, epsg2epsg_nodata
    export geotile_build, geotile_define, geotile_download_granules!, geotile_extent
    export geotile_search_granules, geotile_id, granules_load, pointextract, points_plus
    export range, regular_grid, regular_grid_extents, searchdir, setpaths, utm_epsg, within
    export madnorm, ts_fit, geotile_utm!, geoid, geotile_extract_dem, dem_height
    export geotile_aggrigate_reduce, geotile_ts_fit, itslive_zone, itslive_epsg, itslive_proj!
    export normalize, centroid, dist_ll2xy, dh_ll2xy, dh_ll2aa, itslive_paramfiles 
    export itslive_extract, geotile_extract_mask, geotile_track_offset, track_offset_dh
    export project_paths, project_geotiles, geotile_pointextract, project_products, nt2extent
    export region_extent, geotile_subset,  geotile_subset!, geotile_merge_height, allfiles
    export geotile_build_hugonnet, hstack_catalogue, geotile_offset, geotile_read, EpsgPoints
    export geotiles_w_mask, analysis_paths

end # module