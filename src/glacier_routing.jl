begin
    using Altim
    using GeoDataFrames
    using Rasters
    using Statistics
    import GeometryOps as GO
    import GeoInterface as GI
    using CairoMakie
    using ProgressMeter
    using SortTileRecursiveTree
    using FileIO
    using DataFrames
    using NCDatasets
    using Dates
end


## Populate Paths
begin
    route_rgi2ocean = ["19"] # list of RGIIds that flow to the ocean.. MUST BE A STRING

    # load in local paths
    paths = Altim.pathlocal

    #geomfile_rgi6 = joinpath(paths.data_dir, "GlacierOutlines/rgi60/rgi60_Global.gpkg")
    # paths are hardcoded for now to make it easier collaborate
    glacier_outlines_path = "/mnt/bylot-r3/data/GlacierOutlines/rgi60/rgi60_Global.gpkg"
    copernicus_dem_path = "/mnt/baffin-r1/shared_data/copernicus-dem-30m/DEM.vrt"
    basins_path = "/mnt/bylot-r3/data/rivers/BasinATLAS_Data_v10.gdb/BasinATLAS_v10_lev02.geojson"

    # path to the file containing the glaciers with their minimum elevations
    # this file was created by the script "get_glacier_min_elev.jl"
    glacier_lowest_point_path = replace(glacier_outlines_path, ".gpkg" => "_lowest_point.gpkg")
    glacier_routing_path= replace(glacier_outlines_path, ".gpkg" => "_routing.jld2")
    basins_path = "/mnt/bylot-r3/data/rivers/BasinATLAS_Data_v10.gdb/BasinATLAS_v10_lev02.geojson"


    # corrected rivers are contianed in two separate folders 
    if false
        rivers_path1 = joinpath(paths.data_dir, "rivers/Collins2024/riv_pfaf_ii_MERIT_Hydro_v07_Basins_v01_GLDAS_COR")
        rivers_path2 = joinpath(paths.data_dir, "rivers/Collins2024/riv_pfaf_ii_MERIT_Hydro_v07_Basins_v01_NOTIN_GLDAS_COR/")
        rivers_paths = readdir(rivers_path1; join=true);
        rivers_paths =reduce(vcat, (rivers_paths, readdir(rivers_path2; join=true)));
        glacier_melt_rivers_path = joinpath(rivers_path1, "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.jld2")

    else
        #rivers_path = joinpath(paths.data_dir, "rivers/MERIT_Hydro_v07_Basins_v01_bugfix1")
        rivers_path = joinpath(paths.data_dir, "rivers/MERIT_Hydro_v07_Basins_v01")
        rivers_paths =readdir.(rivers_path; join=true);
        glacier_melt_rivers_path = joinpath(rivers_path, "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.jld2")
    end

    # path to the file containing the basins
    #basins_path = joinpath(paths.data_dir, "rivers/BasinATLAS_Data_v10.gdb/BasinATLAS_v10_lev02.geojson")
    rivers_paths = rivers_paths[Base.contains.(rivers_paths, ".shp")];
end


## find the lowest point on the glacier
if !isfile(glacier_lowest_point_path)
    glaciers = GeoDataFrames.read(glacier_outlines_path)
    dem = Raster("/mnt/devon-r2/shared_data/copernicus-dem-30m/DEM.vrt"; lazy=true)
    mes = get_minimum_elevations(dem, glaciers.geom)
    glaciers.min_elev_long = first.(mes)
    glaciers.min_elev_lat = last.(mes)
    GeoDataFrames.write(glacier_lowest_point_path, glaciers)
end


## find hydrology level 2 basins
if !isfile(glacier_melt_rivers_path)
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
    @showprogress Threads.@threads :greedy for (idx, pt) in enumerate(pts)
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
    col_names = [ "geometry",  "COMID", "NextDownID"]
    riv = GeoDataFrames.read.(rivers_paths)
    rivers = reduce(vcat, r[:,col_names] for r in riv)
    tree = STRtree(rivers.geometry; nodecapacity=3)

    river_indices = zeros(Union{Int64,Missing}, size(pts))
    @showprogress Threads.@threads :greedy for (idx, pt) in enumerate(pts)
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

    glaciers[!, :RiverIDTrace] .= [Int64[]]
    @showprogress  for r in eachrow(glaciers)
        r.RiverIDTrace = Altim.trace_downstream(r.RiverID, rivers.COMID, rivers.NextDownID; maxiters=length(rivers.COMID))
    end

    # RiverIDTrace is a vector of integers, and therefor must be saved as a JLD2 file
    FileIO.save(glacier_routing_path, Dict("glaciers" => glaciers))

    # identify just those rivers that recieve glacier meltwater
    glacier_melt_rivers = unique(reduce(vcat, glaciers.RiverIDTrace))
    rivers_glacier_melt = rivers[in.(rivers.COMID, Ref(glacier_melt_rivers)), :]

    FileIO.save(glacier_melt_rivers_path, Dict("rivers" => rivers_glacier_melt))
end

if !isfile(glacier_routing_path)
#begin
    glaciers = load(glacier_routing_path, "glaciers")
    rivers = load(glacier_melt_rivers_path, "rivers")
    #rivers = GeoDataFrames.read(glacier_melt_rivers_path)

    # load in per-glacier data
    # p = Altim.allfiles("/mnt/bylot-r3/data/binned/2deg/"; fn_endswith = "_perglacier.jld2")
    fn = "/mnt/bylot-r3/data/binned/2deg/glacier_dh_best_cc_meanmadnorm3_v01_filled_ac_p2_synthesized_perglacier.jld2"
    glacier_runoff = load(fn, "glaciers")

    # loop for each glacier
    glaciers[!, :glacier_runoff] = [zeros(dims(glacier_runoff[1,:runoff])) for r in 1:nrow(glaciers)]

    # NOTE: This can not be parallelized because the glaciers.max_runoff variable is modified
    @showprogress for g in eachrow(glacier_runoff)
        valid = .!isnan.(g.runoff) 
        if any(valid)
            idx = findfirst(g.RGIId .== glaciers.RGIId)
            # needs to be += because glaciers can span multiple geotiles
            glaciers[idx, :glacier_runoff] .+= g.runoff
        end
    end

    rivers[!,  :glacier_runoff] = [zeros(dims(glacier_runoff[1,:runoff])) for r in 1:nrow(rivers)]
    # NOTE: This can not be parallelized because the glaciers.glacier_runoff variable is modified
    @showprogress for g in eachrow(glaciers)
        if g["O1Region"] in route_rgi2ocean
            continue
        end

        idx = [findfirst(rivers.COMID .== id) for id in g.RiverIDTrace]
        for i in idx
            rivers[i, :glacier_runoff] .+= g.glacier_runoff
        end
    end

    # runoff units are in native area averaged units (often m w.e.) multiplied by km^2  
    rivers[!, :glacier_runoff_monthly] = [zeros(12) for r in 1:nrow(rivers)]
    for r in eachrow(rivers)
        #r = first(eachrow(rivers))
        # divide by 1000 to convert to m w.e. to km w.e.
        r.glacier_runoff[2:end] = (r.glacier_runoff[2:end] .- collect(r.glacier_runoff[1:end-1]))/1000
        r.glacier_runoff_monthly = Altim.nanmean.(groupby(r.glacier_runoff, :date => Bins(month, 1:12)))
    end

    # RiverIDTrace is a vector of integers, and therefor must be saved as a JLD2 file
    FileIO.save(glacier_melt_rivers_path, Dict("rivers" => rivers))
end



rivers = load(glacier_melt_rivers_path, "rivers")
lines(rivers[1500, :glacier_runoff])
lines(rivers[1500,:glacier_runoff_monthly])


river_Q_path = joinpath(paths.data_dir, "rivers/Collins2024/Qout_pfaf_ii_GLDAS_COR_M_1980-01_2009-12_utc/")
river_Q_files = readdir(river_Q_path; join=true)

#for fn in river_Q_files
fn = river_Q_files[1]
ds = NCDatasets.Dataset(fn)
ddates = Dim{:date}(DateTime.(ds[:time]))
drivid = Dim{:rivid}(rivers.COMID)

for (i, rivid) in enumerate(ds[:rivid])
#i = 10; rivid = ds[:rivid][i]
    index = rivers.COMID .== rivid
    if any(index)
        v = collect(skipmissing(ds[:Qout][i,:]))
        rivers[findfirst(index), :Q][:] = v;
    end
end


Q = ds[:Qout][:,:]
Q = DimArray(ds[:Qout], (rivid=ds[:rivid], time=ds[:time]))
T = ds[:time][:]


f = groupby(r.glacier_runoff, :date => Bins(month, 12))

#rsync -r devon:/mnt/bylot-r3/data/rivers/MERIT_Hydro_v07_Basins_v01_bugfix1 ~/data/rivers/