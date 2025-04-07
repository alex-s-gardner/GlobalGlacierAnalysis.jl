# Glacier Routing Analysis Script
#
# This script analyzes glacier meltwater routing and discharge by:
#
# 1. Finding Lowest Points:
#    - Identifies lowest elevation point within each glacier boundary
#    - Uses 30m DEM to determine where meltwater enters river network
#
# 2. Mapping Drainage Networks:
#    - Maps glaciers to drainage basins and closest river
#    - Traces downstream river paths for each glacier
#    - Identifies rivers receiving glacier meltwater
#    - Determines ocean vs inland termination points
#
# 3. Processing Glacier Variables:
#    - Converts per-glacier geotile data to single values
#    - Interpolates variables (elevation change, mass balance, etc) to monthly values
#    - Converts units from m/year to m³/s where needed
#
# 4. Routing Glacier Runoff:
#    - Maps glacier runoff to river segments
#    - Calculates monthly mean runoff rates
#    - Computes fraction of river flow from glacier melt
#
# 5. Analyzing Discharge Points:
#    - Identifies terminal points where rivers/glaciers end
#    - Maps direct ocean drainage points
#    - Incorporates discharge measurement points
#    - Aggregates data into gridded format for visualization
#
# 6. Calculating Statistics:
#    - Converts runoff to water equivalent units
#    - Fits trends and seasonal cycles to time series
#    - Calculates spatial averages on lat/lon grid
#    - Produces final datasets for analysis and visualization
#
# The script processes global glacier and river datasets to understand:
# - How glacier meltwater enters and moves through river networks
# - Relative contributions of glacier melt to river discharge
# - Spatial and temporal patterns in glacier changes and runoff


#begin

# Import packages needed for glacier:
begin
    using Altim
    using GeoDataFrames
    using Rasters
    using Statistics
    import GeometryOps as GO
    import GeoInterface as GI
    import GeoFormatTypes as GFT
    using CairoMakie
    using ProgressMeter
    using SortTileRecursiveTree
    using FileIO
    using DataFrames
    using NCDatasets
    using Dates
    using NetCDF
    using LsqFit
    using Interpolations
    using Arrow
    using LsqFit
end


# This section sets up paths and constants needed for glacier routing analysis
begin
    route_rgi2ocean = ["19"] 

    # Time conversion constant
    SECONDS_PER_YEAR = 365.25 * 24 * 60 * 60

    #NOTE: glacier variables are in unit of m i.e.
    volume2mass = Altim.δice / 1000

    # Load local configuration paths
    paths = Altim.pathlocal

    # Derived data paths
    glacier_lowest_point_path = replace(paths[:glacier_individual_outlines], ".gpkg" => "_lowest_point.gpkg")
    glacier_routing_path = replace(paths[:glacier_individual_outlines], ".gpkg" => "_routing.arrow")

    # River network paths
    # Note: Two options for river data - using corrected or uncorrected data

    # Uncorrected river paths (currently used)
    rivers_paths = Altim.allfiles(paths[:river]; fn_endswith= "MERIT_Hydro_v07_Basins_v01.shp")
    glacier_rivers_path = joinpath(paths[:river], "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.arrow")


    # Output paths for processed data
    glacier_rivers_runoff_path = replace(glacier_rivers_path, ".arrow" => "_runoff.arrow")
    glacier_rivers_runoff_qout_path = replace(glacier_rivers_runoff_path, ".arrow" => "_qout.nc")
    glacier_sinks_path = joinpath(paths[:data_dir], "glacier_sinks.fgb")
    glacier_sinks_grouped_path = joinpath(paths[:data_dir], "glacier_sinks_grouped.fgb")


    glacier_vars_fns = reduce(vcat,Altim.allfiles.(["/mnt/bylot-r3/data/binned_unfiltered/2deg/", "/mnt/bylot-r3/data/binned/2deg/"]; fn_endswith="synthesized_perglacier.jld2"))
    # glacier_vars_fns = ["/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_dh_best_cc_meanmadnorm3_v01_filled_ac_p1_synthesized_perglacier.jld2"]

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

## Find the lowest elevation point for each glacier
# This section reads glacier outlines and finds the lowest elevation point within each glacier boundary
# by overlaying the glacier geometries on a 30m digital elevation model (DEM).
# The longitude and latitude coordinates of these lowest points are stored in the glacier dataset.
# This is important for determining where glacier meltwater enters the river network.
if !isfile(glacier_lowest_point_path)
    glaciers = GeoDataFrames.read(paths[:glacier_individual_outlines])
    dem = Raster("/mnt/devon-r2/shared_data/copernicus-dem-30m/DEM.vrt"; lazy=true)
    mes = get_minimum_elevations(dem, glaciers.geom)
    glaciers.min_elev_long = first.(mes)
    glaciers.min_elev_lat = last.(mes)
    GeoDataFrames.write(glacier_lowest_point_path, glaciers)
end


## Find drainage basins and rivers for each glacier
# This section maps glaciers to their drainage basins and closest rivers:
#
# 1. Reads glacier lowest point data and basin boundaries
# 2. Uses spatial indexing to efficiently find which basin each glacier point falls within
# 3. Finds the closest river to each glacier using spatial search with expanding buffers
# 4. Traces downstream river paths for each glacier to determine full routing
# 5. Identifies rivers that receive glacier meltwater and maps glaciers to downstream rivers
# 6. Determines which rivers terminate at the ocean vs inland
# 7. Saves the processed routing information to files
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
    tree = STRtree(basins.geometry; nodecapacity=3)

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
    col_names = [ "geometry",  "COMID", "NextDownID", "up1"]
    riv = GeoDataFrames.read.(rivers_paths)
    rivers = reduce(vcat, r[:,col_names] for r in riv)
    tree = STRtree(rivers.geometry; nodecapacity=3)

    river_indices = zeros(Union{Int64,Missing}, size(pts))
    @showprogress dt=1 desc="Locating closest river..." Threads.@threads :greedy for (idx, pt) in enumerate(pts)
    #idx = 3;
    #pt = pts[idx]
        # Query the tree for the point
        potential_river_idxs = Int64[];
        buff = 0;
        while isempty(potential_river_idxs)
            buff += 1_000 # distance in meters
            potential_river_idxs = SortTileRecursiveTree.query(tree, Altim.to_latlong_polygon(Altim.UnitSphericalCap(pt, buff), 16))
        end

        if length(potential_river_idxs) == 1
            river_indices[idx] = potential_river_idxs[1]
        else
            river_distances = [Altim.haversine_distance(pt, rivers.geometry[river_idx]) for river_idx in potential_river_idxs]
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
        r.RiverIDTrace = Altim.trace_downstream(r.RiverID, rivers.COMID, rivers.NextDownID; maxiters=length(rivers.COMID))
    end

    # save the glaciers with the routing information
    GeoDataFrames.write(glacier_routing_path, glaciers)

    # identify just those rivers that recieve glacier meltwater
    glacier_melt_rivers = unique(reduce(vcat, glaciers.RiverIDTrace))
    rivers_glacier_melt = rivers[in.(rivers.COMID, Ref(glacier_melt_rivers)), :]

    # map glaciers to each downstream rivers
    rivers_glacier_melt[!,  :glacier_table_row] .= [Vector{Int64}()]

    @showprogress dt=1 desc="Mapping glaciers to downstream rivers..." Threads.@threads :greedy for r in eachrow(rivers_glacier_melt)
        r.glacier_table_row =  findall(in.(Ref(r.COMID), glaciers.RiverIDTrace))
    end

    rivers_glacier_melt[!, :RGIId] = [glaciers.RGIId[idx] for idx in rivers_glacier_melt.glacier_table_row]

    # determine which rivers flow to the ocean
    river2ocean_path = joinpath(paths.data_dir, "rivers/Collins2024/riv_pfaf_ii_MERIT_Hydro_v07_Basins_v01_coast")
    river2ocean_files = readdir(river2ocean_path; join=true)
    river2ocean_files = filter(fn -> Base.contains(fn, ".shp"), river2ocean_files)
    rivers_glacier_melt[!, :ocean_terminating] .= false

    river2ocean = GeoDataFrames.read(river2ocean_files[1])
    for fn in river2ocean_files[2:end]
        try
            river2ocean = reduce(vcat, (river2ocean , GeoDataFrames.read(fn)))
        catch
        end
    end

    idx = [any(rid .== river2ocean.COMID) for rid in rivers_glacier_melt.COMID]
    rivers_glacier_melt[idx, :ocean_terminating] .= true

    GeoDataFrames.write(glacier_rivers_path, rivers_glacier_melt)
end


# This section processes glacier variables, interpolates them to a consistent time series, routes runoff to rivers and saves the results:
#
# 1. Reads glacier river network and per-glacier geotile variables
# 2. Converts per-glacier geotile data (where glaciers can span multiple geotiles) into single values per glacier
# 3. Verifies glacier IDs match between datasets and copies over area and time series variables
# 4. Converts units from m/year to m³/s for relevant variables
# 5. Interpolates glacier variables (runoff) to the 1st of every month:
#    - Uses quadratic B-spline interpolation
#    - Handles missing data (NaN values)
#    - Ensures runoff remains non-negative
# 6. Routes glacier runoff by:
#    - Copying runoff values to corresponding river segments
#    - Summing runoff from multiple glaciers per segment
#    - Accumulating flow downstream through river network
# 7. Saves results as NetCDF with:
#    - River segment IDs
#    - Monthly glacier runoff flux (m3/s)
#    - Metadata and attributes
#
# The output NetCDF contains:
# - Dimensions: time and river segment ID (COMID)
# - Variables: glacier runoff flux in m3/s
# - Attributes: units, descriptions, etc.


begin # [4 minutes]
    @showprogress dt=1 desc="routing glacier variables..." for glacier_vars_fn in glacier_vars_fns
        glaciers = copy(GeoDataFrames.read(glacier_routing_path)) # getting a strange issue with scope later on so using copy
        geotile_perglacier = load(glacier_vars_fn, "glaciers")

        tsvars=["runoff"]

         #[1 min]

        # ensure that glaciers0 RGIId order matches glaciers order
        @assert glaciers0.RGIId == glaciers.RGIId
        glaciers[!, :area_km2] = glaciers0[:, :area_km2]

        for tsvar in tsvars
            glaciers[!, tsvar] = glaciers0[:, tsvar]
            glaciers = colmetadata!(glaciers, tsvar, "date", colmetadata(glaciers0, tsvar, "date"), style=:note)
        end

        glaciers = Altim.mie2cubicms!(glaciers; tsvars)

        # interpolate glacier runoff to the 1st of every month
        for tsvar in tsvars # [30 seconds]
        #tsvar = "runoff"
            t = colmetadata(glaciers, tsvar, "date")
            yrange = extrema(Dates.year.(t))
            t0 = [Dates.DateTime(y,m,1) for m in 1:12, y in yrange[1]:yrange[2]][:];

            tms = Dates.datetime2epochms.(t)
            dt = unique(tms[2:end] .- tms[1:end-1])
            @assert length(dt) == 1

            t0ms = Dates.datetime2epochms.(t0)

            #scale for efficiency of interpolation
            t0ms = (t0ms .- tms[1]) ./  dt[1] .+ 1

            Threads.@threads for g in eachrow(glaciers) #[30 seconds]
                A = g[tsvar];
                valid = .!isnan.(A)
                if any(valid)
                    itp = extrapolate(interpolate(A[valid], BSpline(Quadratic())), NaN)
                    foo = itp(t0ms .- findfirst(valid) .+ 1)
                    (tsvar == "runoff") ? (foo[foo .< 0] .= 0) : nothing
                    g[tsvar] = foo
                else
                    g[tsvar] = fill(NaN, length(t0ms))
                end
            end
            
            #update column metadata
            colmetadata!(glaciers, tsvar, "date", t0; style=:note)
        end

        # this takes 1 min to load
        rivers = GeoDataFrames.read(glacier_rivers_path)

        rivers[!, :basin02] = floor.(Int16, rivers.COMID./1000000)
        rivers[!, :headbasin] = rivers.up1 .== 0

        # if there is no input from upstream then it is a head basin
        # This is needed as only a subset of the full river network is routed
        rivers[!, :headbasin] .|= .!in.(rivers.COMID, Ref(unique(rivers.NextDownID)))

        # add glacier runoff to the rivers
        for tsvar in tsvars
            rivers[!,  tsvar] = [zeros(length(glaciers[1,tsvar])) for r in 1:nrow(rivers)]
            colmetadata!(rivers, tsvar, "date", colmetadata(glaciers, tsvar, "date"), style=:note)
        end
        
        # glacier inputs to rivers
        valid_temporal_range = falses(size(glaciers[1,first(tsvars)]))
        map(tsvar -> valid_temporal_range .|= (.!isnan.(glaciers[1, tsvar])), tsvars)

        dcomid = Dim{:COMID}(rivers.COMID)
        dTi = Ti(colmetadata(glaciers, first(tsvars), "date")[valid_temporal_range])
        river_inputs = Dict()
        river_flux = Dict()
        
        for tsvar in tsvars
            river_inputs[tsvar] = zeros(eltype(glaciers[1,tsvar]), (dTi, dcomid))

            # need to sum as multiple glaciers can flow into the same river
            for r in eachrow(glaciers[glaciers.RiverID .!= 0, :])
                river_inputs[tsvar][:, At(r.RiverID)] .+= r[tsvar][valid_temporal_range]
            end

            # accululate inputs downstream
            river_flux[tsvar] = Altim.flux_accumulate!(river_inputs[tsvar], rivers.COMID, rivers.NextDownID, rivers.headbasin,  rivers.basin02) #[< 1 second]

            # save data as netcdf
            # data needs to be Float32 (not Real) for saving as netcdf
            river_flux[tsvar] = Float32.(river_flux[tsvar])
        end

        # save as netcdf
        glacier_rivers_runoff_qout_path = replace(glacier_vars_fn, ".jld2" => "_rivers.nc")
        NCDataset(glacier_rivers_runoff_qout_path,"c") do ds
            data_dims = dims(river_flux[first(tsvars)])

            # First all dimensions need to be defined
            for dim in data_dims
                dname = string(DimensionalData.name(dim));
                defDim(ds, dname, length(dim))
            end

            # now add the variables
            for dim in data_dims
                dname = string(DimensionalData.name(dim));
                defVar(ds, dname, val(val(dim)), (dname,))
            end

            for tsvar in tsvars
                v = defVar(ds, "$tsvar", parent(river_flux[tsvar]), string.(DimensionalData.name.(data_dims)))
                v.attrib["units"] = "m3 s-1"
                v.attrib["long_name"] = "glacier $tsvar river flux"

                # add global attributes
                ds.attrib["title"] = "river flux from modeled glacier $tsvar"
            end
        end
    end
end