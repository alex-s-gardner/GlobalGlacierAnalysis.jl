# ============================================================================
# Glacier Routing Module
#
# This module handles the routing of glacier meltwater through river networks.
#
# Key Functionality:
# - Imports necessary packages for glacier routing analysis
# - Sets up paths and constants for data processing
# - Finds lowest elevation points for each glacier to determine meltwater entry points
# - Maps glaciers to drainage basins and river networks
# - Traces downstream flow paths for glacier meltwater
# - Calculates and accumulates glacier runoff through river networks
# - Exports results to NetCDF format for further analysis
#
# Data Sources:
# - Glacier outlines and attributes
# - Digital elevation models (30m resolution)
# - River network data from MERIT Hydro
# - GLDAS flux data (corrected and uncorrected versions)
#
# Workflow:
# 1. Identify lowest points in glacier boundaries
# 2. Map glaciers to drainage basins and nearest rivers
# 3. Trace downstream flow paths
# 4. Calculate glacier runoff inputs to river segments
# 5. Accumulate runoff through river networks
# 6. Export results to NetCDF format
# ============================================================================

# Import packages needed for glacier:
begin
    import GlobalGlacierAnalysis as GGA
    using GeoDataFrames
    using Rasters
    using Statistics
    import GeometryOps as GO
    import GeoInterface as GI
    import GeoFormatTypes as GFT
    using CairoMakie
    using ProgressMeter
    using SortTileRecursiveTree
    using Extents
    using FileIO
    using DataFrames
    using NCDatasets
    using Dates
    using LsqFit
    using Interpolations
    using Arrow
    using LsqFit
    using SortTileRecursiveTree
    using Unitful
    Unitful.register(GGA.MyUnits)
end


# Sets up paths and constants needed for glacier routing analysis.
#
# Configuration includes:
# - File paths for input/output data
# - Physical constants for unit conversions
# - River network data paths (corrected and uncorrected versions)
# - GLDAS flux data configuration
#
# This section establishes all necessary file paths and constants required for the
# glacier routing workflow, including paths to glacier outlines, river networks,
# and output locations for processed data.
begin
    glacier_summary_riverflux_file = replace(GGA.pathlocal[:glacier_summary], ".nc" => "_riverflux.nc")
    route_rgi2ocean = ["19"] 

    # Time conversion constant
    SECONDS_PER_YEAR = 365.25 * 24 * 60 * 60

    #NOTE: glacier variables are in unit of m i.e.
    volume2mass = GGA.δice / 1000

    # Load local configuration paths
    paths = GGA.pathlocal

    # Derived data paths
    glacier_lowest_point_path = replace(paths[:glacier_individual], ".gpkg" => "_lowest_point.gpkg")
    glacier_routing_path = replace(paths[:glacier_individual], ".gpkg" => "_routing.arrow")

    # River network paths
    # Note: Two options for river data - using corrected or uncorrected data

    # Uncorrected river paths (currently used)
    rivers_paths = GGA.allfiles(paths[:river]; fn_endswith="MERIT_Hydro_v07_Basins_v01.shp", fn_startswith="riv_pfaf")
    glacier_rivers_path = joinpath(paths[:river], "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.arrow")

    # Output paths for processed data
    glacier_rivers_runoff_path = replace(glacier_rivers_path, ".arrow" => "_runoff.arrow")
    glacier_rivers_runoff_qout_path = replace(glacier_rivers_runoff_path, ".arrow" => "_qout.nc")
    glacier_sinks_path = joinpath(paths[:data_dir], "glacier_sinks.fgb")
    glacier_sinks_grouped_path = joinpath(paths[:data_dir], "glacier_sinks_grouped.fgb")


    # TODO: I need to update to use the new glacier_flux path
    # glacier_flux = joinpath(paths[:project_dir], "gardner2025_glacier_summary.nc")
    
    glacier_vars_fns = reduce(vcat,GGA.allfiles.(["/mnt/bylot-r3/data/binned_unfiltered/2deg/", "/mnt/bylot-r3/data/binned/2deg/"]; fn_endswith="synthesized_perglacier.jld2"))
    # glacier_vars_fns = ["/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_dh_best_cc_nmad3_v01_filled_ac_p1_synthesized_perglacier.jld2"]

    # GLDAS flux data configuration
    use_corrected_flux = false
    if use_corrected_flux 
        # Paths for corrected flux data using in-situ observations
        river_Q_path = joinpath(paths.data_dir, "rivers/Collins2024/Qout_pfaf_ii_GLDAS_COR_M_1980-01_2009-12_utc/")
        river_Q_files = readdir(river_Q_path; join=true)
        lev2_basins = [splitpath(fn)[end][11:12] for fn in river_Q_files]

        # Backfill missing data with ensemble simulations
        river_Q_path_ens = joinpath(paths.data_dir, "rivers/Collins2024/Qout_pfaf_ii_GLDAS_ENS_M_1980-01_2009-12_utc/")
        river_Q_files_ens = readdir(river_Q_path_ens; join=true)

        for fn in river_Q_files_ens
            if !any(splitpath(fn)[end][11:12] .== lev2_basins)
                river_Q_files = push!(river_Q_files, fn)
            end
        end
    else
        # Use uncorrected ensemble simulations only
        river_Q_path_ens = joinpath(paths.data_dir, "rivers/Collins2024/Qout_pfaf_ii_GLDAS_ENS_M_1980-01_2009-12_utc/")
        river_Q_files = readdir(river_Q_path_ens; join=true)
    end
end

# Find the lowest elevation point for each glacier
#
# Identifies the lowest elevation point within each glacier boundary by overlaying
# glacier geometries on a 30m Copernicus DEM. These points represent where glacier
# meltwater likely enters the river network.
#
# The function:
# 1. Reads glacier outline geometries from GeoPackage
# 2. Overlays geometries on the 30m DEM to find minimum elevation points
# 3. Extracts longitude and latitude coordinates of these points
# 4. Saves the augmented glacier dataset with minimum elevation coordinates
#
# Only runs if the output file doesn't already exist to avoid redundant processing.
if !isfile(glacier_lowest_point_path)
    glaciers = GeoDataFrames.read(paths[:glacier_individual])
    dem = Raster("/mnt/devon-r2/shared_data/copernicus-dem-30m/DEM.vrt"; lazy=true)
    mes = get_minimum_elevations(dem, glaciers.geom)
    glaciers.min_elev_long = first.(mes)
    glaciers.min_elev_lat = last.(mes)
    GeoDataFrames.write(glacier_lowest_point_path, glaciers)
end


# Find drainage basins and rivers for each glacier
#
# Maps glaciers to their drainage basins and river networks by:
# 1. Identifying which drainage basin contains each glacier's lowest point
# 2. Finding the closest river to each glacier's lowest point
# 3. Tracing the downstream flow path for glacier meltwater
# 4. Mapping glaciers to all downstream river segments
# 5. Identifying which rivers terminate at the ocean
# 6. Assigning continent and country information to each river segment
#
# The function uses spatial indexing for efficient querying and processes data in parallel
# where possible. Results are saved to disk to avoid redundant processing in future runs.
#
# Only executes if the output file doesn't already exist.
if !isfile(glacier_rivers_path)
#begin
    glaciers = GeoDataFrames.read(glacier_lowest_point_path)

    # read the file containing the basins
    basins = GeoDataFrames.read(paths[:river_basins])

    # construct points from minimum latitude and longitude
    pts = GI.Point.(glaciers.min_elev_long, glaciers.min_elev_lat)

    # Convert the `basin` geometries to Julia geometries,
    # so that access is super fast
    basin_geometries = GO.tuples(basins.geometry)

    # Build a spatial index for the basins
    tree = SortTileRecursiveTree.STRtree(basins.geometry; nodecapacity=3)

    # Iterate over the points and find the first basin that intersects with each point
    # Prefilter by the basins that could possibly intersect with each point
    basin_indices = zeros(Union{Int64,Missing}, size(pts))
    @showprogress dt=1 desc="Locating drainage basins..." Threads.@threads :greedy for (idx, pt) in enumerate(pts)
        # Query the tree for the point
        potential_basin_idxs = SortTileRecursiveTree.query(tree, pt)
        # Find the first basin that intersects with the point
        first_intersecting_basin = findfirst(Base.Fix1(GO.intersects, pt), view(basin_geometries, potential_basin_idxs))
        # Store the basin index, or missing if no intersection was found
        basin_indices[idx] = isnothing(first_intersecting_basin) ? missing : potential_basin_idxs[first_intersecting_basin]
    end

    valid = .!ismissing.(basin_indices)
    glaciers[!, :BasinID] = Vector{Union{Missing,Int64}}(fill(missing, length(pts)))
    glaciers[valid, :BasinID] = basin_indices[valid]

    # identify the river that each glacier is closest to
    # load in river reaches
    col_names = [ "geometry",  "COMID", "NextDownID", "up1", "lengthkm"]
    riv = GeoDataFrames.read.(rivers_paths)
    rivers = reduce(vcat, r[:,col_names] for r in riv)
    tree = SortTileRecursiveTree.STRtree(rivers.geometry; nodecapacity=3)

    river_indices = zeros(Union{Int64,Missing}, size(pts))
    @showprogress dt=1 desc="Locating closest river..." Threads.@threads :greedy for (idx, pt) in enumerate(pts)
    #idx = 3;
    #pt = pts[idx]
        # Query the tree for the point
        potential_river_idxs = Int64[];
        buff = 0;
        while isempty(potential_river_idxs)
            buff += 1_000 # distance in meters
            potential_river_idxs = SortTileRecursiveTree.query(tree, GGA.to_latlong_polygon(GGA.UnitSphericalCap(pt, buff), 16))
        end

        if length(potential_river_idxs) == 1
            river_indices[idx] = potential_river_idxs[1]
        else
            river_distances = [GGA.haversine_distance(pt, rivers.geometry[river_idx]) for river_idx in potential_river_idxs]
            river_indices[idx] = first(potential_river_idxs[river_distances .== minimum(river_distances)])
        end
    end

    glaciers[!, :RiverID] = rivers.COMID[river_indices]

    # set the river id to 0 for glaciers that flow directly to the ocean
    for reg in route_rgi2ocean
        idx = glaciers.O1Region .== reg
        if any(idx)
            glaciers[idx, :RiverID] .= 0
        end
    end

    glaciers[!, :RiverIDTrace] .= [Int64[]]
    @showprogress dt=1 desc="Tracing downstream..." for r in eachrow(glaciers)
    #r = eachrow(glaciers)[1]
        r.RiverIDTrace = GGA.trace_downstream(r.RiverID, rivers.COMID, rivers.NextDownID; maxiters=length(rivers.COMID))
    end

    # save the glaciers with the routing information
    GeoDataFrames.write(glacier_routing_path, glaciers)

    # identify just those rivers that recieve glacier meltwater
    glacier_melt_rivers = unique(reduce(vcat, glaciers.RiverIDTrace))
    rivers = rivers[in.(rivers.COMID, Ref(glacier_melt_rivers)), :]

    # map glaciers to each downstream rivers
    rivers[!,  :glacier_table_row] .= [Vector{Int64}()]

    @showprogress dt=1 desc="Mapping glaciers to downstream rivers..." Threads.@threads :greedy for r in eachrow(rivers)
        r.glacier_table_row =  findall(in.(Ref(r.COMID), glaciers.RiverIDTrace))
    end

    rivers[!, :RGIId] = [glaciers.RGIId[idx] for idx in rivers.glacier_table_row]

    # determine which rivers flow to the ocean
    river2ocean_path = joinpath(paths.data_dir, "rivers/Collins2024/riv_pfaf_ii_MERIT_Hydro_v07_Basins_v01_coast")
    river2ocean_files = readdir(river2ocean_path; join=true)
    river2ocean_files = filter(fn -> Base.contains(fn, ".shp"), river2ocean_files)
    rivers[!, :ocean_terminating] .= false

    river2ocean = GeoDataFrames.read(river2ocean_files[1])
    for fn in river2ocean_files[2:end]
        try
            river2ocean = reduce(vcat, (river2ocean , GeoDataFrames.read(fn)))
        catch
        end
    end

    idx = [any(rid .== river2ocean.COMID) for rid in rivers.COMID]
    rivers[idx, :ocean_terminating] .= true

    rivers = GeoDataFrames.read(glacier_rivers_path)
    rivers[!, :continent] .= ""
    rivers[!, :country] .= ""
    rivers[!, :centroid] = GO.centroid.(rivers.geometry)
    tree = SortTileRecursiveTree.STRtree(rivers[!, :centroid])

    # Assign continent to each river segment
    continents = GeoDataFrames.read(paths.continents)

    @time Threads.@threads for r in eachrow(continents) # [4.5 min]
        query_result = SortTileRecursiveTree.query(tree, GI.extent(r.geometry))
        index = GO.intersects.(rivers.centroid[query_result], Ref(r.geometry))
        rivers[query_result[index], :continent] .= r.CONTINENT
    end

    # Assign country to each river segment
    countries = GeoDataFrames.read(paths[:countries])
    @time Threads.@threads for r in eachrow(countries) # [2.5 min]
        query_result = SortTileRecursiveTree.query(tree, GI.extent(r.geometry))
        index = GO.intersects.(rivers.centroid[query_result], Ref(r.geometry))
        rivers[query_result[index], :country] .= r.SOVEREIGNT
    end

    GeoDataFrames.write(glacier_rivers_path, rivers[:, Not(:centroid)])
end


# Route glacier runoff through river networks.
#
# This function processes glacier runoff data and routes it through downstream river networks:
# 1. Loads glacier runoff data from a NetCDF file
# 2. Converts runoff from Gt to kg/s
# 3. Interpolates values to monthly time steps
# 4. Maps glacier runoff to river segments
# 5. Accumulates runoff downstream through the river network
# 6. Saves the resulting river flux data to a NetCDF file
#
# Arguments
# - `glacier_routing_path`: Path to glacier routing information
# - `GGA.pathlocal[:glacier_summary]`: Path to input glacier runoff data
# - `glacier_summary_riverflux_file`: Path for output river flux data
# - `glacier_rivers_path`: Path to river network data
# - `paths`: Dictionary of project paths

begin # [4 minutes]
    vars2route=["runoff"]

    glaciers = copy(GeoDataFrames.read(glacier_routing_path)) # getting a strange issue with scope later on so using copy

    glacier_flux_nc = NCDataset(GGA.pathlocal[:glacier_summary])
        
    glacier_flux0 = Dict()

    for tsvar in vars2route
    #tsvar = "runoff"
        var0 = GGA.nc2dd(glacier_flux_nc[tsvar])

        # convert from Gt to kg³/s
        dTi = dims(var0, :Ti)
        Δdate = Second.(unique(Day.(val(dTi)[2:end] .- val(dTi)[1:end-1])))

        # Gt per month
        glacier_flux0[tsvar] = var0[:, 2:end]
        glacier_flux0[tsvar][:,:] =  diff(parent(var0), dims=2)
        glacier_flux0[tsvar] = uconvert.(u"kg", glacier_flux0[tsvar]) ./ Δdate
    end

    # interpolate to the 1st of every month
    glacier_flux = Dict()
    for tsvar in vars2route # [30 seconds]
        var0 = glacier_flux0[tsvar]
        dTi = dims(var0, :Ti)
        drgi = dims(var0, :rgiid)
        
        yrange = extrema(Dates.year.(dTi))
        t0 = [Dates.DateTime(y,m,1) for m in 1:12, y in yrange[1]:yrange[2]][:];
        t0 = t0[(t0 .>= dTi[1]) .& (t0 .<= dTi[end])]
        dTi0 = Ti(t0)
        glacier_flux[tsvar] = zeros(drgi,dTi0)* u"kg/s"

        tms = Dates.datetime2epochms.(val(dTi));
        dt = unique(tms[2:end] .- tms[1:end-1]);
        @assert length(dt) == 1

        t0ms = Dates.datetime2epochms.(t0)

        #scale for efficiency of interpolation
        t0ms = (t0ms .- tms[1]) ./  dt[1] .+ 1

        Threads.@threads for rgi in drgi
            A = var0[At(rgi), :]
            valid = .!isnan.(A)
            if any(valid)
                itp = extrapolate(interpolate(A[valid], BSpline(Quadratic())), NaN * u"kg/s")
                foo = itp.(t0ms .- findfirst(valid) .+ 1)
                (tsvar == "runoff") ? (foo[foo .< 0* u"kg/s"] .= 0 * u"kg/s") : nothing
                glacier_flux[tsvar][At(rgi), :] = foo
            else
                glacier_flux[tsvar][At(rgi), :] = fill(NaN, length(t0ms))* u"kg/s"
            end
        end
    end

    # this takes 1 min to load
    rivers = GeoDataFrames.read(glacier_rivers_path)
    rivers[!, :basin02] = floor.(Int16, rivers.COMID./1000000)
    rivers[!, :headbasin] = rivers.up1 .== 0
    rivers[!, :longitude] = getindex.(GO.centroid.(rivers.geometry), 1)
    rivers[!, :latitude] = getindex.(GO.centroid.(rivers.geometry), 2)

    # if there is no input from upstream then it is a head basin
    # This is needed as only a subset of the full river network is routed
    rivers[!, :headbasin] .|= .!in.(rivers.COMID, Ref(unique(rivers.NextDownID)))

    dcomid = Dim{:COMID}(rivers.COMID)
    dTi = dims(glacier_flux[first(vars2route)], :Ti)
    river_inputs = Dict()
    river_flux = Dict()

    for tsvar in vars2route
        tsvar = "runoff"
        river_inputs[tsvar] = zeros(eltype(glacier_flux[tsvar]), (dTi, dcomid))

        # need to sum as multiple glaciers can flow into the same river
        for r in eachrow(glaciers[glaciers.RiverID .!= 0, :])
            # r = eachrow(glaciers[glaciers.RiverID .!= 0, :])[1]
            river_inputs[tsvar][:, At(r.RiverID)] .+= glacier_flux[tsvar][At(r.RGIId), :]
        end

        # accululate inputs downstream
        river_flux[tsvar] = GGA.flux_accumulate!(river_inputs[tsvar], rivers.COMID, rivers.NextDownID, rivers.headbasin, rivers.basin02) #[< 1 second]

        # save data as netcdf
        # data needs to be Float32 (not Real) for saving as netcdf
        river_flux[tsvar] = Float32.(river_flux[tsvar])
    end

    dstack = DimStack(getindex.(Ref(river_flux), vars2route); name=Symbol.(vars2route))
    
    # save to netcdf
    NCDataset(glacier_summary_riverflux_file, "c") do ds

        data_dims = dims(dstack)

        # First all dimensions need to be defined
        for dim in data_dims
            dname = string(DimensionalData.name(dim));
            defDim(ds, dname, length(dim))
        end

        # now add the variables
        for dim in data_dims
            dname = string(DimensionalData.name(dim));
            d = defVar(ds, dname, val(val(dim)), (dname,))
            if DateTime <: eltype(dim)
                d.attrib["cf_role"] = "timeseries_id";
            end
        end

        for vaname in keys(dstack)
            v = defVar(ds, "$vaname", ustrip.(parent(dstack[vaname])), string.(DimensionalData.name.(data_dims)))
            v.attrib["units"] = string(Unitful.unit(dstack[vaname][1]))
        end

        # add latitude and longitude [This is a hack until I can get geometry into a NetCDF]
        defVar(ds, "latitude", rivers.latitude, string.(DimensionalData.name.(data_dims[2:2])))
        defVar(ds, "longitude", rivers.longitude, string.(DimensionalData.name.(data_dims[2:2])))

        # add global attributes
        ds.attrib["title"] = "river flux from glacier model"
        ds.attrib["version"] = "beta - " * Dates.format(now(), "yyyy-mm-dd")
        ds.attrib["featureType"] = "timeSeries";

        println("glacier output saved to $glacier_summary_riverflux_file")
    end;
end
