"""
    rgi_fix.jl

Fix RGI (Randolph Glacier Inventory) shapefile geometries.

This script processes RGI glacier outline shapefiles by:
1. Creating a copy of the original directory with "-fix" suffix
2. Converting all polygon geometries to MULTIPOLYGON type
3. Ensuring all geometries use XY dimensions only

The conversion is necessary for compatibility with certain GIS operations
that require consistent geometry types across all features.
"""
using ArchGDAL

folder = "/Users/gardnera/data/GlacierOutlines/RGI2000-v7.0-C-global"

cp(folder, folder * "-fix", force=true)

shp = allfiles(folder; subfolders=true, fn_endswith=".shp")

for s in shp
    d = ArchGDAL.read(s)
    fn_out = replace(s, folder => (folder * "-fix"))
    ArchGDAL.unsafe_gdalvectortranslate([d], ["-nlt" "MULTIPOLYGON" "-dim" "XY"], dest=fn_out)
end