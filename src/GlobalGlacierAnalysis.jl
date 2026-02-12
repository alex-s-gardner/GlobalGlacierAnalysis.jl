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
   #using SpaceLiDAR
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
   using Metaheuristics
   using Extents
   using Interpolations
   using StaticArrays
   using Logging
   using Printf
   using NearestNeighbors
   using StatsBase
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
   import SortTileRecursiveTree


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
   include("utilities_manuscript.jl")
   include("mapzonal.jl")

   function Makie._register_argument_conversions!(::Type{P}, attr::Makie.ComputeGraph, user_kw) where {P}
      dim_converts = to_value(get!(() -> Makie.DimConversions(), user_kw, :dim_conversions))
      args = attr.args[]
      Makie.add_convert_kwargs!(attr, user_kw, P, args)
      kw = attr.convert_kwargs[]
      args_converted = Makie.convert_arguments(P, args...; kw...)
      status = Makie.got_converted(P, Makie.conversion_trait(P, args...), args_converted)
      force_dimconverts = Makie.needs_dimconvert(dim_converts)
      if force_dimconverts && status == true
         Makie.add_dim_converts!(attr, dim_converts, args)
      elseif (status === true || status === Makie.SpecApi)
         # Nothing needs to be done, since we can just use convert_arguments without dim_converts
         # And just pass the arguments through
         map!(attr, :args, :dim_converted) do args
               return Ref{Any}(args)
         end
      elseif isnothing(status) || status == true # we don't know (e.g. recipes)
         Makie.add_dim_converts!(attr, dim_converts, args)
      elseif status === false
         if args_converted !== args
               # Not at target conversion, but something got converted
               # This means we need to convert the args before doing a dim conversion
               map!(attr, :args, :recursive_convert) do args
                  return Makie.convert_arguments(P, args...)
               end
               Makie.add_dim_converts!(attr, dim_converts, args_converted, :recursive_convert)
         else
               Makie.add_dim_converts!(attr, dim_converts, args)
         end
      end
      #  backwards compatibility for plot.converted (and not only compatibility, but it's just convenient to have)

      map!(attr, [:dim_converted, :convert_kwargs], :converted) do dim_converted, convert_kwargs
         x = Makie.convert_arguments(P, dim_converted...; convert_kwargs...)
         if x isa Tuple
               return x
         elseif x isa Union{Makie.PlotSpec, AbstractVector{Makie.PlotSpec}, Makie.GridLayoutSpec}
               return (x,)
         else
               error("Result needs to be Tuple or SpecApi")
         end
      end
      converted = attr[:converted][]
      n_args = length(converted)
      map!(attr, :converted, [Makie.argument_names(P, n_args)...]) do converted
         return converted # destructure
      end

      Makie.add_input!((k, v) -> Ref{Any}(v), attr, :transform_func, identity)

      Makie.add_input!(attr, :f32c, :uninitialized)

      return
   end
end