"""
    cop30_vrt_build.jl

Build virtual raster (VRT) files for Copernicus DEM 30m data.

This script creates consolidated VRT files for each type of Copernicus DEM 30m data product.
It recursively searches through the specified directory for all TIF files and groups them
by suffix (DEM, WBM, EDM, FLM, HEM) to create a single VRT file for each product type.

The Copernicus DEM uses WGS84-G1150 (EPSG:4326) horizontal coordinate system and 
EGM2008 (EPSG:3855) vertical datum.

Available product types:
- DEM.tif: Digital Elevation Model
- WBM.tif: Water Body Mask (8-bit unsigned integer)
- EDM.tif: Editing Mask (8-bit unsigned integer)
- FLM.tif: Filling Mask (8-bit unsigned integer)
- HEM.tif: Height Error Mask (32-bit floating point)
"""

using ArchGDAL, GlobalGlacierAnalysis

# info page
# https://spacedata.copernicus.eu/web/guest/collections/copernicus-digital-elevation-model?p_l_back_url=%2Fweb%2Fguest%2Fsearch%3Fq%3DDEM

# 2021 release available on aws:
# https://registry.opendata.aws/copernicus-dem/
# download all dems using AWS CLI:
# aws s3 cp s3://copernicus-dem-30m/ /mnt/devon-r2/shared_data/copernicus-dem-30m/ --recursive --no-sign-request
path2folder = "/mnt/devon-r2/shared_data/copernicus-dem-30m/"

# 2019 release was downloaded with ancillary files: /mnt/devon-r2/shared_data/COP-DEM_GLO-30-DGED_PUBLIC/
# path2folder = "/mnt/devon-r2/shared_data/COP-DEM_GLO-30-DGED_PUBLIC/";

# Coordinate Reference System: Horizontal WGS84-G1150 (EPSG 4326) (DGED & DTED format), Vertical EGM2008 (EPSG 3855)
suffix = ["WBM.tif", "DEM.tif", "DEM_int16.tif", "EDM.tif", "FLM.tif", "HEM.tif"];

# EDM.tif = Editing Mask: 8 Bit unsigned integer, GeoTIFF
# FLM.tif = Filling Mask: 8 Bit unsigned integer, GeoTIFF
# HEM.tif = Height Error Mask: 32 Bit floating point, GeoTIFF
# WBM.tif = Water Body Mask: 8 Bit unsigned integer, GeoTIFF

# this takes a long time as it maps 440972 files
items = [item for item in walkdir(path2folder)]
files = [[joinpath(root, file) for file in files] for (root, dirs, files) in items]
files = reduce(vcat,files)

# build single .vrt for all .tif files with suffix
for suff in suffix 
    out_vrt = joinpath(path2folder, replace(suff, ".tif" => ".vrt"))
    ind = endswith.(files,suff)
    in_tifs = files[ind]
    in_tifs = ArchGDAL.read.(in_tifs)
    ArchGDAL.gdalbuildvrt(in_tifs; dest=out_vrt) do vrt
        ArchGDAL.write(vrt, out_vrt)
    end
end