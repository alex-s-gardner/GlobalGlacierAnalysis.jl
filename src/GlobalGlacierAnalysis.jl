module GlobalGlacierAnalysis

   module MyUnits
      using Unitful; @unit Gt "Gt" M 1u"Pg" false
   end

   using Unitful
   Unitful.register(MyUnits)
   using GlobalGlacierAnalysis.MyUnits
   
   # import geographic packages
   using Proj
   using GeoArrays
   using SpaceLiDAR
   using Geodesy
   using FastGeoProjections
   using Rasters
   using Shapefile
   using LibGEOS
   using GeoInterface

   # file IO
   using FileIO
   using Arrow
   using HTTP
   using MAT
   using CSV
   using NCDatasets
   using JLD2
   

   # standard packages
   using DataFrames
   using Dates
   using Statistics

   # plotting packages
   using ColorSchemes
   using CairoMakie
   
   # other packages
   using OffsetArrays
   using ImageFiltering
   using Extents
   using Interpolations
   using StaticArrays
   using Logging
   using Printf
   using NearestNeighbors
   using Stencils  
   using NonlinearSolve
   using ScatteredInterpolation
   using LsqFit
   using Distances
   using DataInterpolations
   using Distributions
   using Images
   using FlexiJoins
   using Loess
   using LinearAlgebra
   using Rotations
   using CoordinateTransformations
   using SparseArrays
   using SortTileRecursiveTree
   using AbstractTrees
   using ProgressMeter
   using Random
   using LazyGrids
   using GeometryBasics
   

   # custom packages create as part of this project
   using RangeExtractor
   using GeoTiles
   using BinStatistics

   # import packages
   import DimensionalData as DD
   import GeometryOps as GO
   import GeoInterface as GI
   import GeoFormatTypes as GFT
   import GeometryOpsCore as GOC


   # this is where the local paths are defined
   include("local_paths.jl")

   # add utilities
   include("utilities_project.jl")
   include("utilities_build_archive.jl")
   include("utilities_hugonnet.jl")
   include("utilities_gemb.jl")
   include("utilities_main.jl")
   include("utilities_binning.jl")
   include("utilities_binning_lowlevel.jl")
   include("utilities_routing.jl")
   include("utilities_synthesis.jl")
   include("utilities_postprocessing.jl")
   include("utilities.jl")
   include("utilities_plotting.jl")
   include("utilities_readers.jl")
   include("mapzonal.jl")
end
