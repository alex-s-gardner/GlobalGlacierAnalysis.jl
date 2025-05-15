module GlobalGlacierAnalysis

    module MyUnits
       using Unitful; @unit Gt "Gt" M 1u"Pg" false
    end

    function __init__()
       Unitful.register(MyUnits)
    end

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
    using NetCDF
    using BinStatistics
    using FastGeoProjections
    using GeoTiles
    using NearestNeighbors
    using Stencils  
    using NonlinearSolve
  
    using Rasters
    using MAT
    using ScatteredInterpolation
    using DimensionalData
    using LsqFit
    using CSV
    using NCDatasets

    using Distances
    using DataInterpolations
    using ColorSchemes
    using CairoMakie
    using JLD2
    using Distributions
    using Images
    using FlexiJoins
    using ProgressMeter
    using Loess
    using Unitful
    Unitful.register(MyUnits)

    using LinearAlgebra
    using Rotations
    using CoordinateTransformations
    using GeometryBasics 
    using SparseArrays
    using SortTileRecursiveTree
    using AbstractTrees
    using ProgressMeter
    using Random

    using RangeExtractor
    import DimensionalData as DD
    using LazyGrids

    using LibGEOS
    using GeoInterface
    using Shapefile
    using FileIO

    import GeometryOps as GO
    import GeoInterface as GI
    import GeoFormatTypes as GFT
    import GeometryOpsCore
    import LibGEOS as LG
    import GeometryOpsCore as GOC

    #include("utilities.jl")
    include("utilities.jl")
    include("mapzonal.jl")
    include("utilities_project.jl")
    include("utilities_hugonnet.jl")
    include("utilities_gemb.jl")
    include("modelfit_tree_fast.jl")
    include("utilities_hyps.jl")
    include("utilities_plot.jl")
    include("utilities_geotile_processing.jl")
    include("utilities_routing.jl")
    include("utilities_postprocessing.jl")
    # include("/home/gardnera/Documents/GitHub/GlobalGlacierAnalysis.jl/src/STRTreesAbstractTreesExt.jl")

    export Extent
    export bin, binnedfiltering, crop!, decimalyear, download!, epsg2epsg, epsg2epsg_nodata
    export geotile_build, geotile_define, geotile_download_granules!, geotile_extent
    export geotile_search_granules, geotile_id, granules_load, pointextract
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
