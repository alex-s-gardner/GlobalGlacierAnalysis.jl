"""
    pgc_dems_get.jl

Download and process Polar Geospatial Center (PGC) Digital Elevation Models (DEMs).

This script downloads ArcticDEM or REMA data from PGC's public S3 bucket, 
processes the individual tiles, and creates a virtual raster (VRT) file 
that combines all tiles into a single seamless dataset.

The script:
1. Constructs the appropriate URLs based on the selected DEM, version, and resolution
2. Downloads all available tiles using aria2 for parallel downloading
3. Creates a VRT file that virtually combines all downloaded GeoTIFF files

Requirements:
- aria2 command-line download utility must be installed
- Sufficient disk space for the downloaded DEM tiles
"""

using HTTP, JSON, Printf, ArchGDAL

localdatafolder = "/mnt/devon-r2/shared_data/"
dem = "arcticdem" # "arcticdem" -or- "rema" 
version = 4.1 # 2 or 3 
resolution = 10; # 2, 10, 32
downloadstreams = 16

# JSON files found here: 
# https://polargeospatialcenter.github.io/stac-browser/#/external/pgc-opendata-dems.s3.us-west-2.amazonaws.com/pgc-data-stac.json
# navigate products to determine `dataset` and `suffex` 

# this script requires `aria2` to be insalled in the shell: https://github.com/q3aql/aria2-static-builds


# build dataset specific strings
dataset = @sprintf("%s/mosaics/v%.1f/%.0fm", dem, version, resolution)
if dem == "arcticdem"
    if version == 3
        suffex = @sprintf("_%.0fm_v%.1f_reg_dem.tif", resolution, version)
    else 
         suffex = @sprintf("_%.0fm_v%.1f_dem.tif", resolution, version)
    end
else
    suffex = @sprintf("_%.0fm_v%.1f_dem.tif", resolution, version)
end   

# folder where data should be written 
outfolder = joinpath(localdatafolder, dataset)

# read in json stack and parse urls
#https://pgc-opendata-dems.s3.us-west-2.amazonaws.com/arcticdem/mosaics/v4.1/10m/07_40/07_40_10m_v4.1.json
#https://pgc-opendata-dems.s3.us-west-2.amazonaws.com/arcticdem/mosaics/v3.0/10m/07_40/07_40_10m_v3.0.json


resp = HTTP.get("https://pgc-opendata-dems.s3.us-west-2.amazonaws.com/$dataset.json")
str = String(resp.body)
jobj = JSON.Parser.parse(str)
urls = [link["href"] for link in jobj["links"]]
# remove extra urls
urls = urls[contains.(urls,"$(last(split(dataset,"/")))/")]

# make a function to write urls to a txt file
function write_urls!(fn::String, urls::Vector{String})
    open(fn, "w") do f
        for url in urls
            url = replace(url, ".json"=>"")
            url = joinpath(url,"$(last(split(url,"/")))$suffex")
            println(f, url)
        end
    end
    abspath(fn)
end

# write urls to file
fn  = tempname()
url_list = write_urls!(fn, urls)

# make directly if it doesn't exist
if !isdir(outfolder)
    mkpath(outfolder)
end

# download files
cmd = `aria2c -x$downloadstreams -c -d $outfolder -i $url_list`
run(cmd)

# delete url list
rm(url_list)

# build single .vrt for all .tif files
out_vrt = joinpath(outfolder, @sprintf("%s_v%.1f_%.0fm.vrt", dem, version, resolution))
in_tifs =  readdir(outfolder);
in_tifs = in_tifs[endswith.(in_tifs, ".tif")];
in_tifs = joinpath.(outfolder, in_tifs)
in_tifs = ArchGDAL.read.(in_tifs)

ArchGDAL.gdalbuildvrt(in_tifs; dest=out_vrt) do vrt
    ArchGDAL.write(vrt,out_vrt)
end
