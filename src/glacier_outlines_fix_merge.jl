 
"""
    glacier_outlines_fix_merge.jl

Process and merge RGI 7.0 glacier outlines into a single GeoPackage file.

This script:
1. Copies the original RGI 7.0 glacier outlines to a new directory
2. Converts all shapefiles to MULTIPOLYGON geometry with XY dimensions
3. Reads all converted shapefiles into GeoDataFrames
4. Merges all GeoDataFrames into a single dataframe
5. Writes the merged data to a GeoPackage file

The script ensures consistent geometry types across all glacier outlines
for improved compatibility with GIS applications.
"""

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