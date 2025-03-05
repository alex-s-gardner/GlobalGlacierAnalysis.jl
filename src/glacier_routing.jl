# Glacier Routing Analysis Script
#
# This script analyzes glacier meltwater routing and discharge by:
#
# 1. Finding Lowest Points:
#    - Identifies lowest elevation point within each glacier boundary
#    - Uses 30m DEM to determine where meltwater enters river network
#
# 2. Mapping Drainage Networks:
#    - Maps glaciers to drainage basins and closest rivers
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

    # Core input data paths
    glacier_outlines_path = "/mnt/bylot-r3/data/GlacierOutlines/rgi60/rgi60_Global.gpkg"
    copernicus_dem_path = "/mnt/baffin-r1/shared_data/copernicus-dem-30m/DEM.vrt"
    basins_path = "/mnt/bylot-r3/data/rivers/BasinATLAS_Data_v10.gdb/BasinATLAS_v10_lev02.geojson"

    # Derived data paths
    glacier_lowest_point_path = replace(glacier_outlines_path, ".gpkg" => "_lowest_point.gpkg")
    glacier_routing_path = replace(glacier_outlines_path, ".gpkg" => "_routing.arrow")

    # River network paths
    # Note: Two options for river data - using corrected or uncorrected data
    if false
        # Corrected river paths (not currently used)
        rivers_path1 = joinpath(paths.data_dir, "rivers/Collins2024/riv_pfaf_ii_MERIT_Hydro_v07_Basins_v01_GLDAS_COR")
        rivers_path2 = joinpath(paths.data_dir, "rivers/Collins2024/riv_pfaf_ii_MERIT_Hydro_v07_Basins_v01_NOTIN_GLDAS_COR/")
        rivers_paths = readdir(rivers_path1; join=true)
        rivers_paths = reduce(vcat, (rivers_paths, readdir(rivers_path2; join=true)))
        glacier_rivers_path = joinpath(rivers_path1, "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.arrow")

        # Filter river paths to only include shapefiles
        rivers_paths = rivers_paths[Base.contains.(rivers_paths, ".shp")]
    else
        # Uncorrected river paths (currently used)
        rivers_path = joinpath(paths.data_dir, "rivers/MERIT_Hydro_v07_Basins_v01")
        rivers_paths = Altim.allfiles(rivers_path; fn_endswith= "MERIT_Hydro_v07_Basins_v01.shp")
        glacier_rivers_path = joinpath(rivers_path, "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.arrow")
    end

    # Output paths for processed data
    glacier_rivers_runoff_path = replace(glacier_rivers_path, ".arrow" => "_runoff.arrow")
    glacier_rivers_runoff_qout_path = replace(glacier_rivers_runoff_path, ".arrow" => "_qout.nc")
    glacier_sinks_path = joinpath(paths.data_dir, "glacier_sinks.fgb")
    glacier_sinks_grouped_path = joinpath(paths.data_dir, "glacier_sinks_grouped.fgb")

    path2runs_override = ["/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_dh_best_cc_meanmadnorm3_v01_filled_ac_p1_synthesized.jld2"]
    binned_synthesized_file = path2runs_override[1]

    # for sanity checking in QGIS
    glacier_vars_fn = replace(binned_synthesized_file, ".jld2" => "_perglacier.jld2")

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
    glaciers = GeoDataFrames.read(glacier_outlines_path)
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
    basins = GeoDataFrames.read(basins_path)

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
    glaciers = GeoDataFrames.read(glacier_routing_path)
    geotile_perglacier = load(glacier_vars_fn, "glaciers")

    tsvars=["runoff", "dm"]

    # convert from perglacier geotile (glacier can call in multiple geotiles) to perglacier
    glaciers0 = Altim.perglacier_geotile2glacier(geotile_perglacier; tsvars) #[1 min]

    # ensure that glaciers0 RGIId order matches glaciers order
    @assert glaciers0.RGIId == glaciers.RGIId
    glaciers[!, :area_km2] = glaciers0[:, :area_km2]

    for tsvar in tsvars
        glaciers[!, tsvar] = glaciers0[:, tsvar]
        glaciers = colmetadata!(glaciers, tsvar, "date", colmetadata(glaciers0, tsvar, "date"), style=:note)
    end

    glaciers = Altim.mie2cubicms!(glaciers; tsvars)

    # interpolate glacier runoff to the 1st of every month
    @time for tsvar in tsvars # [30 seconds]
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

        @showprogress dt=1 desc="Interpolating glacier $(tsvar)..." Threads.@threads for g in eachrow(glaciers) #[30 seconds]
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


#=
# ----------------------------- THIS IS OLD CODE -------------------------------------------

# This section of code processes river discharge data and calculates fracitonal glacier melt contribution:
#
# 1. Loads river network data with glacier runoff from previous step
# 2. Loads river discharge data from NetCDF files and maps it to glacier river reaches:
#    - Reads discharge time series for each river segment
#    - Interpolates to monthly values on the 15th of each month
#    - Handles missing data and ensures discharge is non-negative
# 3. Calculates monthly mean discharge rates for each river segment
# 4. Computes the maximum fraction of river flow from glacier melt:
#    - For each river segment, calculates glacier runoff / total discharge ratio
#    - Caps the fraction at 1.0 to handle cases where estimated river flow is less than glacier discharge
#    - Handles edge cases where either value is zero or NaN
# 5. Saves the processed data in multiple formats:
#    - Full dataset in Arrow and JLD2 formats with all variables
#    - Simplified version in FlatGeobuf format with key variables for visualization:
#      * River geometry
#      * River segment ID (COMID) 
#      * Whether river terminates at ocean
#      * Maximum glacier melt fraction
begin
    in_fn = glacier_rivers_runoff_path
    rivers = GeoDataFrames.read(in_fn)

    # metadatagymnastics
    metadata = FileIO.load(replace(in_fn, ".arrow" => "metadata.jld2"), "metadata")
    # transfer the metadata
    for k in keys(metadata)
        colmetadata!(rivers, k, "date", metadata[k], style=:note)
    end

    # load in river discharge data and map to the glacier river reaches
    # get dimensions
    #seconds since 1970-01-01 00:00:00 +00:00
    fn = river_Q_files[1]

    # interpolate to the 15th of every month
    t = Dates.DateTime("1970-01-01") .+ Dates.Second.(ncread(fn, "time"))
    yrange = extrema(Dates.year.(t))
    t0 = [Dates.DateTime(y,m,15) for m in 1:12, y in yrange[1]:yrange[2]][:];

    tms = Dates.datetime2epochms.(t)
    dt = unique(tms[2:end] .- tms[1:end-1])
    @assert length(dt) == 1

    t0ms = Dates.datetime2epochms.(t0)

    #scale for efficiency of interpolation
    t0ms = (t0ms .- tms[1]) ./  dt[1] .+ 1

    # get river ids
    rivid = ncread(fn, "rivid")

    #for fn in river_Q_files [m3 s-1]
    ddates = Dim{:date}(DateTime.(t0))
    drivid = Dim{:rivid}(rivers.COMID)
    riverQ = zeros(Float32, (drivid, ddates))

    # using Treads causes crash
    # 1 min on single thread
    @showprogress dt=1 desc="Mapping river discharge..." for fn in river_Q_files
    #fn = river_Q_files[4]
        rividsX = ncread(fn, "rivid")
        matching_rivids = intersect(drivid, rividsX)
        if !isempty(matching_rivids)
            QX = ncread(fn, "Qout")
            for rividX in matching_rivids
            #rividX = first(matching_rivids)
                if any(drivid .== rividX)
                    A = QX[findfirst(rividsX .== rividX), :]
                    valid = .!isnan.(A)
                    if any(valid)
                        itp = extrapolate(interpolate(A[valid], BSpline(Quadratic())), NaN)
                        foo = itp(t0ms .- findfirst(valid) .+ 1)
                        foo[foo .< 0] .= 0
                        riverQ[At(rividX), :] = foo
                    end
                end
            end
        end
    end
 
    # Convert river discharge array to vector of vectors, one per river
    rivers[!, :Qout] = [zeros(Float64, length(ddates)) for r in 1:nrow(rivers)]
    for (i, r) in enumerate(eachrow(riverQ))
        rivers[i, :Qout] = r
    end
    
    # Add metadata to track the dates associated with discharge values
    colmetadata!(rivers, :Qout, "date", collect(ddates), style=:note)
    
    # Initialize monthly discharge array for each river (12 months)
    rivers[!, :Qout_monthly] = [zeros(12) for r in 1:nrow(rivers)]

    for r in eachrow(rivers)
        r.Qout_monthly = Altim.nanmean.(groupby(DimArray(r.Qout, ddates), :date => Bins(month, 1:12)))
    end

    rivers[!,:glacfrac_max] .= 0.
    for (i, r) in enumerate(eachrow(rivers))
        # 0 / 0 = NaN
        fac = r.glacier_runoff_monthly ./ (r.Qout_monthly .+ r.glacier_runoff_monthly);
        valid = .!isnan.(fac)

        # for lots of mountian ranges the estimated river flow is << the glacier dischage
        if any(valid)
            r.glacfrac_max = min(maximum(fac[valid]), 1.0)
        else
            r.glacfrac_max = 0
        end
    end

    out_fn = glacier_melt_rivers_runoff_qout_path
    GeoDataFrames.write(out_fn, rivers)

    # metadatagymnastics
    metadata = Dict()
    for n in colmetadatakeys(rivers)
        metadata[n[1]] = colmetadata(rivers, n[1], collect(n[2])[1])
    end
    FileIO.save(replace(out_fn, ".arrow" => "metadata.jld2"), Dict("metadata" => metadata))

    # for plotting
    GeoDataFrames.write(glacier_melt_rivers_runoff_qout_path_fgb, rivers[:, [:geometry, :COMID, :ocean_terminating, :glacfrac_max]])
end




# This section determines glacier sinks (where glacier meltwater ultimately flows)
# It combines three types of sinks:
# 1. Terminal river points (where rivers end, either at ocean or endorheic basins)
# 2. Glaciers that drain directly to the ocean (not connected to river network)
# 3. Discharge measurement points from observations
begin
    # Load river and glacier data with runoff values
    in_fn = glacier_melt_rivers_runoff_qout_path
    rivers = GeoDataFrames.read(in_fn)

    # metadatagymnastics
    metadata = FileIO.load(replace(in_fn, ".arrow" => "metadata.jld2"), "metadata")
    # transfer the metadata
    for k in keys(metadata)
        colmetadata!(rivers, k, "date", metadata[k], style=:note)
    end
    
    # Load glacier data and transfer metadata similarly
    in_fn = glacier_vars_interpolated_fn
    glaciers = GeoDataFrames.read(in_fn)
    # metadatagymnastics
    metadata = FileIO.load(replace(in_fn, ".arrow" => "metadata.jld2"), "metadata")
    # transfer the metadata
    for k in keys(metadata)
        colmetadata!(glaciers, k, "date", metadata[k], style=:note)
    end

    # Remove a specific river in West Africa that doesn't receive glacier melt
    COMID_exclude = 14088718
    rivers = rivers[rivers.COMID .!= COMID_exclude, :]

    # Load global discharge data
    globaldischage_fn = joinpath(paths.data_dir, "GlacierOutlines/GlacierDischarge/global_glacier_dischage.jld2")
    discharge = load(globaldischage_fn, "discharge")

    # Process per-glacier geotile data
    begin
        fn = "/mnt/bylot-r3/data/binned/2deg/glacier_dh_best_cc_meanmadnorm3_v01_filled_ac_p2_synthesized_perglacier.jld2"
        geotile_perglacier = load(fn, "glaciers")

        # Convert from per-glacier geotile to per-glacier values
        tsvars=["dh", "smb", "fac", "runoff"]
        glaciers0 = Altim.perglacier_geotile2glacier(geotile_perglacier; tsvars)

        # Verify glacier order matches and transfer data
        @assert glaciers0.RGIId == glaciers.RGIId
        glaciers[!, :area_km2] = glaciers0[:, :area_km2]
        for tsvar in tsvars
            glaciers[!, tsvar] = glaciers0[:, tsvar]
            colmetadata!(glaciers, tsvar, "date", colmetadata(glaciers0, tsvar, "date"), style=:note)
        end

        # Convert units from mie/year to m³/s
        glaciers = Altim.mie2cubicms!(glaciers; tsvars)
    end

    # Flag terminal points (where rivers/glaciers end)
    rivers[!, :terminal] = rivers.NextDownID .== 0
    glaciers[!, :terminal] = glaciers.RiverID .== 0

    source_crs1 = GFT.EPSG(4326)

    # Create sinks dataframe and populate with river terminators
    sinks = DataFrame()
    foo = rivers[rivers.terminal, :]
    endpts = [GI.Point(first(GI.coordinates(x)), crs=source_crs1) for x in foo.geometry]
    sinks[!, :geometry] = endpts
    sinks[!, :river] .= true
    sinks[!, :endorheic] = .!(foo.ocean_terminating)
    sinks[!, :frontal] .= false
    sinks[!, :glacier_runoff] .= foo.glacier_runoff
    colmetadata!(sinks, :glacier_runoff, "date", colmetadata(glaciers, :runoff, "date"), style=:note)

    # Add glaciers that drain directly to ocean
    foo = glaciers[glaciers.RiverID .== 0, :]
    endpts = [GI.Point(f.min_elev_long, f.min_elev_lat, crs=source_crs1) for f in eachrow(foo)]
    sinks = vcat(sinks, DataFrame(geometry=endpts, river=false, endorheic=false, frontal=false, glacier_runoff=foo.runoff))

    # Add discharge measurement points
    discharge[!, :discharge_m3s] = max.(Ref(0.0), (discharge.discharge_gtyr * 1E3^3) / SECONDS_PER_YEAR)
    discharge_ts = [d*ones(length(colmetadata(glaciers, :runoff, "date"))) for d in discharge.discharge_m3s]
    endpts = [GI.Point(f.longitude, f.latitude, crs=source_crs1) for f in eachrow(discharge)]
    sinks = vcat(sinks, DataFrame(geometry=endpts, river=false, endorheic=false, frontal=true, glacier_runoff=discharge_ts))

    # Calculate mean runoff and add colors for visualization
    sinks[!, :glacier_runoff_mean] = Altim.nanmean.(sinks.glacier_runoff)
    ocean = "rgba(183, 72, 75, 0.60)"
    endorheic = "rgba(69, 116, 40, 0.60)"
    sinks[!, :color] .= ocean
    sinks[sinks.endorheic, :color] .= endorheic

    # Save results
    out_fn = glacier_sinks_path
    GeoDataFrames.write(out_fn, sinks[:, Not(:glacier_runoff)]; crs=GFT.EPSG(4326))

    # metadatagymnastics
    metadata = Dict()
    for n in colmetadatakeys(sinks)
        metadata[n[1]] = colmetadata(sinks, n[1], collect(n[2])[1])
    end
    FileIO.save(replace(out_fn, ".fgb" => "metadata.jld2"), Dict("metadata" => metadata))
end


# This section processes and aggregates glacier sink data into a gridded format
# for visualization purposes:

# 1. Reads glacier sink data from file
# 2. Creates a 5-degree lat/lon grid to aggregate sinks into
# 3. For each grid cell:
#    - Identifies sinks that fall within the cell bounds
#    - Calculates weighted average location based on runoff values
#    - Sums total runoff for the cell
# 4. Filters out grid cells with no runoff
# 5. Saves the aggregated results to file for plotting
begin
    sinks = GeoDataFrames.read(glacier_sinks_path)
    
    Δdeg = 5.0
    lat = (-90+Δdeg/2):Δdeg:(90-Δdeg/2)
    lon = (-180+Δdeg/2):Δdeg:(180-Δdeg/2)
    source_crs1 = GFT.EPSG(4326)

    sinks_grouped = DataFrame(extent=[extent = Extent(X=(x - Δdeg / 2, x + Δdeg / 2), Y=(y - Δdeg / 2, y + Δdeg / 2)) for y in lat, x in lon][:])
    sinks_grouped[!, :geometry] .= [GI.Point(0., 0., crs=source_crs1) for _ in 1:nrow(sinks_grouped)]

    # group for each lat/lon
    sinks_coords = GI.coordinates.(sinks.geometry)
    sinks_x = getindex.(sinks_coords,1)
    sinks_y= getindex.(sinks_coords,2)

    sinks_grouped[!, :glacier_runoff_mean] .= 0.0
    for sg in eachrow(sinks_grouped)
    #sg = eachrow(sinks_grouped)[1]
        index = Altim.within.(Ref(sg.extent), sinks_x, sinks_y) .& (.!isnan.(sinks.glacier_runoff_mean))

        if any(index)
            # weight location by dischage/glacier_runoff_monthly
            x0 = sum(sinks_x[index] .* sinks[index, :glacier_runoff_mean]) / sum(sinks[index, :glacier_runoff_mean])
            y0 = sum(sinks_y[index] .* sinks[index, :glacier_runoff_mean]) / sum(sinks[index, :glacier_runoff_mean])

            sg.geometry = GI.Point(x0, y0, crs=source_crs1)
            sg.glacier_runoff_mean = sum(sinks[index, :glacier_runoff_mean])
        end
    end

    sinks_grouped = sinks_grouped[(sinks_grouped.glacier_runoff_mean .!= 0.0) .& (.!isnan.(sinks_grouped.glacier_runoff_mean)), :]

    GeoDataFrames.write(glacier_sinks_grouped_path, sinks_grouped[:, [:geometry, :glacier_runoff_mean]])
    # rsync -ar devon:/mnt/bylot-r3/data/glacier_sinks_grouped.fgb ~/data/Altim/
end


# This section processes glacier time series data and calculates gridded statistics:
#
# 1. Sets up paths for input/output files
# 2. Loads glacier data from file and restores metadata
# 3. Converts glacier runoff values from m³/s to meters water equivalent (mwe):
#    - Integrates flow rates over time periods
#    - Normalizes by glacier area
#    - Converts units from km³ to mwe
# 4. Calculates temporal statistics for each glacier:
#    - Fits trends and seasonal cycles to dh, smb, fac, runoff time series
#    - Extracts offset, trend, amplitude and phase parameters
# 5. Creates 1-degree lat/lon grid for spatial averaging:
#    - Groups glaciers by grid cell
#    - Calculates area-weighted averages of statistics
#    - Filters out empty grid cells
# 6. Saves results to files:
#    - Individual glacier statistics
#    - Gridded averages
#
# The output files can be used to analyze and visualize spatial patterns
# in glacier changes and runoff across different regions.

#begin
    glacier_rates_path = joinpath(paths.data_dir, "glacier_rates.fgb")
    glacier_rates_grouped_path = joinpath(paths.data_dir, "glacier_rates_grouped.fgb")

    # 
    glaciers = GeoDataFrames.read(glacier_vars_interpolated_fn)


    ## this is a hack until arrow can support metadata
    tsvars= ["dh",  "smb", "fac", "runoff"]

    for tsvar in tsvars
        geotile_perglacier = load(glacier_vars_fn, "glaciers")
        t = collect(dims(geotile_perglacier[1, tsvar], :date))
        yrange = extrema(Dates.year.(t))
        t0 = [Dates.DateTime(y,m,15) for m in 1:12, y in yrange[1]:yrange[2]][:];

        # add back  metadata
        colmetadata!(glaciers, tsvar, "date", t0, style=:note)
    end
    ##

    #convert from m^3 s-1 to mwe
    @showprogress dt=1 desc="Converting to mwe..." Threads.@threads for glacier in eachrow(glaciers)
        for tsvar in tsvars
            #glacier = first(eachrow(glaciers))
            #tsvar = first(tsvars)
            v = copy(glacier[tsvar])
            valid = .!isnan.(glacier[tsvar])
            t0 = colmetadata(glaciers, tsvar, "date")
            t0 = vcat(Second(0), Second.(diff(t0))) # elapsed time in seconds
            t0 = [t.value for t in t0]
            v = v .* t0 # convert to m^3
            glacier[tsvar][valid] = cumsum(vcat(0, v[valid]))[2:end] / 1E3^3 / glacier.area_km2 * 1000
        end
    end

    # sanity check
    #glaciers1 = FileIO.load(glacier_vars_fn, "glaciers")
    #id = "RGI60-01.01000"
    #lines(glaciers[findfirst(glaciers.RGIId .== id), :runoff])
    #lines(glaciers1[findfirst(glaciers1.RGIId .== id), :runoff])

    # this takes ~ 10 minutes
    @time Altim.df_tsfit!(glaciers, [:dh, :smb, :fac, :runoff]; datelimits = (DateTime(2000,1,1), DateTime(2023,1,1)))
    
    outvars = ["geom", "RGIId"]
    for tsvar in tsvars
        push!(outvars, "$(tsvar)_offset")
        push!(outvars, "$(tsvar)_trend")
        push!(outvars, "$(tsvar)_amplitude")
        push!(outvars, "$(tsvar)_phase")
    end

    # find glacier centroids
    glaciers[!, :geom] = GI.Point.(GO.centroid.(glaciers.geom))
    valid = .!isnan.(getindex.(GI.coordinates.(glaciers.geom),1))
    GeoDataFrames.write(glacier_rates_path, glaciers[valid, outvars]; crs=GFT.EPSG(4326))
    # rsync -ar devon:/mnt/bylot-r3/data/glacier_rates.fgb ~/data/Altim/

    # now average rates by 1 deg lat x 1 deg long
    outvars = setdiff(outvars, ["geom", "RGIId"])
    Δdeg = 1.0
    lat = (-90+Δdeg/2):Δdeg:(90-Δdeg/2)
    lon = (-180+Δdeg/2):Δdeg:(180-Δdeg/2)
    source_crs1 = GFT.EPSG(4326)

    glaciers_grouped = DataFrame(extent=[extent = Extent(X=(x - Δdeg / 2, x + Δdeg / 2), Y=(y - Δdeg / 2, y + Δdeg / 2)) for y in lat, x in lon][:])
    glaciers_grouped[!, :geom] .= [GI.Point(0., 0.) for _ in 1:nrow(glaciers_grouped)]

    # group for each lat/lon
    glaciers_coords = GI.coordinates.(glaciers.geom)
    glaciers_x = getindex.(glaciers_coords,1)
    glaciers_y= getindex.(glaciers_coords,2)

    for vname in outvars
        glaciers_grouped[!, vname] .= NaN
    end

    @showprogress dt=1 desc="Grouping glaciers..." Threads.@threads for gg in eachrow(glaciers_grouped)
        #gg = eachrow(glaciers_grouped)[1]
        index = Altim.within.(Ref(gg.extent), glaciers_x, glaciers_y)
        gg.geom = GI.Point(mean(gg.extent[1]), mean(gg.extent[2]))
        if any(index)
            
            for vname in outvars

                # weight location by dischage/glacier_runoff_monthly
                #x0 = sum(glaciers_x[index] .* glaciers[index, vname]) / sum(glaciers[index, vname])
                #y0 = sum(glaciers_y[index] .* glaciers[index, vname]) / sum(glaciers[index, vname])

                #gg.geometry = GI.Point(x0, y0, crs=source_crs1)
                
                glacier_area = glaciers[index, :area_km2]
                gg[vname] = sum(glaciers[index, vname] .* glacier_area) ./ sum(glacier_area)
            end
        end
    end

    glaciers_grouped = glaciers_grouped[:, Not(:extent)]
    valid = .!isnan.(glaciers_grouped.dh_offset)
    rename!(glaciers_grouped, :geom => :geometry)
    GeoDataFrames.write(glacier_rates_grouped_path, glaciers_grouped[valid, :]; crs=GFT.EPSG(4326))

    # rsync -ar devon:/mnt/bylot-r3/data/glacier_rates_grouped.arrow ~/data/Altim/
end
=#
