using Rasters # https://github.com/rafaqz/Rasters.jl/tree/use_BitArray
using Shapefile
using ArchGDAL
using GeoFormatTypes
using Downloads
using ZipFile
using GeometryOps
using GeoInterface

# the shapefiles can be downloaded from NSIDC: https://nsidc.org/data/nsidc-0770/versions/7 
# but you will require an Earth Data Login

# for convenience I've temporarily put the files here
fn = "RGI2000-v7.0-C-02_western_canada_usa.zip"
url2shp =  "https://its-live-data.s3.amazonaws.com/test/$fn.zip"

# download and unzip
# Downloads.download(url2shp, "$fn.zip")
# run(`unzip $fn`)

# Load glacier polygons
glacier = Shapefile.Handle("$fn/$fn.shp")

# reporject files [doesn't work but let's keeep going]
epsg = GeoFormatTypes.EPSG("EPSG:32609")
glacier_rep = GeometryOps.reproject(glacier; target_crs=epsg);

# reutrns: ERROR: StackOverflowError:

# define a Raster to wich the shapefile will be rasterized
gridsize = 25;
min_x = 273500;
max_x = 2507100;
min_y = 4077900;
max_y = 7234200;

x = X(min_x:gridsize:max_x);
y = Y(min_y:gridsize:max_y);
A = Raster(zeros(UInt8,y,x))
setcrs(A, epsg)

# rasterize polygons
@time foo = Rasters.rasterize!(count, A, glacier_rep.geom; threaded=true, shape=:polygon)

GeoArrays.write("junk.tif", A; tiled=true, compress="ZSTD")

# at the command line this takes about 10s
run(`ogr2ogr -t_srs EPSG:3857 $(pwd())/$(fn)/$(fn)_32609.shp $(pwd())/$fn/$fn.shp`)# This takse a few seconds
run(`gdal_rasterize -burn 1.0 -tr 25 25 -a_nodata 0.0 -ot Byte  $(pwd())/$(fn)/$(fn)_32609.shp $(pwd())/$(fn)/$(fn)_32609.tif`)



using Rasters
using GeoFormatTypes

epsg = GeoFormatTypes.EPSG("EPSG:32609")
gridsize = 250;
min_x = 273500;
max_x = 2507100;
min_y = 4077900;
max_y = 7234200;

x = X(min_x:gridsize:max_x);
y = Y(min_y:gridsize:max_y);
A = Raster(zeros(UInt8,y,x))
setcrs(A, epsg)

Rasters.write( "out.tif", A)

Rasters.write("out.tif", A, tiled = true, compress="ZSTD")
