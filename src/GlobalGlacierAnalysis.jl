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
end 
