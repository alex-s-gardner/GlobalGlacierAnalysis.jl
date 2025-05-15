"""
    nasadem_vrt_build.jl

Build Virtual Raster Tiles (VRTs) for NASADEM datasets.

This script creates VRT files that virtually combine multiple NASADEM GeoTIFF files
into single, seamless datasets. It processes four NASADEM data types:
- wgs84_hgt.tif: Elevation data in WGS84 reference frame
- egm08_hgt.tif: Elevation data with EGM2008 geoid
- num.tif: Number of observations per pixel
- swb.tif: Surface water body mask

The VRTs enable efficient access to the entire NASADEM dataset without
requiring all component files to be loaded into memory simultaneously.

For downloading the NASADEM data see:
https://portal.opentopography.org/datasetMetadata?otCollectionID=OT.032021.4326.2 
"""
using ArchGDAL, GlobalGlacierAnalysis

path2folder = "/mnt/devon-r2/shared_data/NASADEM/mosaic/";
suffix = ["wgs84_hgt.tif", "egm08_hgt.tif", "num.tif", "swb.tif"];

# build single .vrt for all .tif files with suffix
for suff in suffix
    suff = first(suffix)
    out_vrt = joinpath(path2folder, replace(suff, ".tif" => ".vrt"))
    in_tifs =  searchdir(path2folder, suff);
    in_tifs = joinpath.(path2folder, in_tifs)
    in_tifs = ArchGDAL.read.(in_tifs)

    ArchGDAL.gdalbuildvrt(in_tifs; dest=out_vrt) do vrt
        ArchGDAL.write(vrt, out_vrt)
    end
end
