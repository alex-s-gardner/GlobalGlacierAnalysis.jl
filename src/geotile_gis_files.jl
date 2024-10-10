using Altim
using GeoDataFrames
import GeoInterface as GI
using DataFrames


paths = Altim.pathlocal
geotile_width = 2;

geotiles_file = "geotiles.gpkg"
if !isfile(geotiles_file)
    geomfile_rgi6 = joinpath(paths.data_dir, "GlacierOutlines/rgi60/rgi60_Global.gpkg")

    glacier_geom = GeoDataFrames.read(geomfile_rgi6)
    geotiles = Altim.geotiles_w_mask(geotile_width)
    geotiles[!, :geometry] = Altim.GeoTiles.define(2)[:, :geometry]
    geotiles = geotiles[:, Not([:extent])]
    geotiles = geotiles[geotiles.glacier_frac.>0, :]

    geotiles[!, :glacierized] .= false
    Threads.@threads for geotile in eachrow(geotiles)
        #geotile = first(geotiles)
        geotile.glacierized = any(GeometryOps.intersects.(glacier_geom[:, 1], Ref(geotile.geometry)))
    end
    geotiles = geotiles[geotiles.glacierized, :]
    GeoDataFrames.write("geotiles.gpkg", geotiles)
else
    geotiles = GeoDataFrames.read("geotiles.gpkg")
end



mission0 = ["GEDI", "ASTER", "ICESat", "ICESat-2", "WorldView"]
lat_limit0 = [51.6, 83.0, 86.0, 88., 90.]

lat_limit = Vector{GI.LineString{false,false,Vector{GI.Point{false,false,Tuple{Int64,Float64},Nothing}},Nothing,Nothing}}()
mission = Vector{String}()

for i in eachindex(mission0)
    for s in [1, -1]
        p1 = GI.Point.(-180, lat_limit0[i]*s);
        p2 = GI.Point.(180, lat_limit0[i]*s);
        lat_limit =  push!(lat_limit, GI.LineString([p1, p2]))
        mission = push!(mission,mission0[i])
    end
end

mission_limits = DataFrame(geometry=lat_limit, mission = mission)
GeoDataFrames.write("mission_limits.gpkg", mission_limits)