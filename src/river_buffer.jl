import GeoFormatTypes as GFT
using Proj
import GeometryOps as GO
import GeoInterface as GI
using CairoMakie
using LibGEOS
using Rasters, ArchGDAL
#using Altim
using ArchGDAL
using GeoDataFrames
using DataFrames
using RangeExtractor
using ProgressMeter
using BinStatistics

glacier_melt_rivers_runoff_qout_path = "/mnt/bylot-r3/data/rivers/MERIT_Hydro_v07_Basins_v01/riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier_runoff_qout.arrow"

# build a vrt of global high resolution population data
begin
    hrsl_dir = "/mnt/bylot-r3/data/dataforgood-fb-data/hrsl-cogs/hrsl_general/v1.5/"
    foo = splitpath(hrsl_dir)
    hrsl_vrt = joinpath(foo[1:end-1]..., foo[end] * ".vrt")
    if .!isfile(hrsl_vrt)
        hrsl_files = Altim.allfiles(hrsl_dir; fn_endswith = ".tif")

        hrsl = ArchGDAL.read.(hrsl_files)
        ArchGDAL.gdalbuildvrt(hrsl; dest=hrsl_vrt) do vrt
            ArchGDAL.write(vrt, hrsl_vrt)
        end
    end
end

# load rivers
rivers = GeoDataFrames.read(glacier_melt_rivers_runoff_qout_path)

# lets make life easier for now
rivers = rivers[vcat([1:20, 80:89, 5220:5230]...), [:geometry]]
GeoDataFrames.write("test_rivers.geojson", rivers)

# buffer rivers
# project to ease grid
source_crs1 = GFT.EPSG(4326)
target_crs1 = GFT.EPSG(6933)
trans = Proj.Transformation(source_crs1, target_crs1; always_xy=true)
buffer_radius = 10_000 #meters


#Threads does not help here
geometry_buffered = []
@showprogress dt = 1 desc = "buffering rivers..." for r in eachrow(rivers)
    #r = eachrow(rivers)[2]
    foo = GO.transform(trans, r.geometry)
    foo_buff = GO.buffer(foo, buffer_radius)
    foo_buff = GO.transform(inv(trans), foo_buff)
    push!(geometry_buffered, foo_buff)
end

# sanity check
#p = plot(geometry_buffered[1])
#plot!.(geometry_buffered[2:end])
#p
#

#add to dataframe
rivers[!, :geometry_buffered] = geometry_buffered

# first geometry

g = rivers.geometry_buffered[1]
g0 =  [GI.Polygon(GI.LinearRing(GI.getgeom(g.geom[1])))]
geometry_buffered_cropped = [];
#p = plot(g0[1])

@showprogress dt = 1 desc = "buffering rivers..." for (i, g) in enumerate(rivers.geometry_buffered[2:end])

    #g = rivers.geometry_buffered[40]
    #i = 36
#println(i)
#g = rivers.geometry_buffered[i+1]

    # find intersection of geometry with all previous geometries
    g = GI.Polygon(GI.LinearRing(GI.getgeom(g.geom[1])))
    gint = GO.intersection.(g0, Ref(g); target=GI.PolygonTrait())

    if isempty(gint[1])
        gdiff = [g]
    #elseif gint[1] isa Array{}
    #    gint = GI.Polygon(GI.LineString(GI.getgeom(gint[1])))
    #    gdiff = GO.difference.(Ref(g), gint; target=GI.PolygonTrait())
    else
        # find difference
        gx = gint[1][1]
        gx = GI.Polygon(GI.LinearRing(GI.getgeom(gx.geom[1])))
        gdiff = GO.difference(g, gx; target=GI.PolygonTrait())
    end
    
   #plot!(gdiff[1])
    push!(geometry_buffered_cropped, gdiff)

    # union with all previous geometries
    g0 = GO.union(g0[1], g; target=GI.PolygonTrait())
end
p