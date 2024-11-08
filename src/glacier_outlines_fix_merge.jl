 
using ArchGDAL
using GeoTiles
using GeoDataFrames

datadir = "/Users/gardnera/data";

geomfile_rgi7_gpkg = joinpath(datadir, "GlacierOutlines/RGI2000-v7.0-G-global-fix/rgi70_Global.gpkg")
folder = joinpath(datadir, "GlacierOutlines/RGI2000-v7.0-G-global")

cp(folder, folder * "-fix", force=true)

shp = GeoTiles.allfiles(folder; subfolders=true, fn_endswith=".shp")

files2merge = String[]
for s in shp
    d = ArchGDAL.read(s)
    fn_out = replace(s, folder => (folder * "-fix"))
    push!(files2merge, fn_out)
    ArchGDAL.unsafe_gdalvectortranslate([d], ["-nlt" "MULTIPOLYGON" "-dim" "XY"], dest=fn_out)
end

gdf = GeoDataFrames.read.(files2merge)
gdf = reduce(vcat, gdf)
GeoDataFrames.write(geomfile_rgi7_gpkg, gdf)