# validate dem extration

using Altim
using DataFrames
using Arrow
using Plots
using Statistics
using GeoArrays
using ArchGDAL
using GDAL



test_data = "/mnt/bylot-r2/shared_data/altim/h_validate.arrow";
df = copy(DataFrame(Arrow.Table(test_data)));

# extract data from a dem
dem = :cop30_v2;

h= vec(dem_height(df.lon[1:100], df.lat[1:100], dem; filter_sigma=nothing, slope=false));
d = df.lidar[1:100] - h

std(d[.!isnan.(d)])

plot(df.hgt[1:100] .- df.dh[1:100], h; seriestype=:scatter)

mean(d)


using GeoArrays
deminfo, goids_folder = Altim.dem_info(dem)
ga = GeoArrays.read(deminfo.fn; masked=false)

x = df.lon[1:100]
y =  df.lat[1:100]
point_epsg = "EPSG:4326"
_, epsg_num = split(point_epsg, ":")
if !contains(ga.crs.val[end-7:end], "\"EPSG\",\"$epsg_num\"")
    xy = epsg2epsg(x, y, point_epsg, ga.crs.val; parse_output=false)
    x = getindex.(xy, 1)
    y = getindex.(xy, 2)
end

# find x y extents
minmax_x = extrema(x)
minmax_y = extrema(y)
dx = ga.f.linear[1, 1];
dy = ga.f.linear[2, 2];

extent = (
    min_x=minmax_x[1] - abs(dx),
    min_y=minmax_y[1] - abs(dy),
    max_x=minmax_x[2] + abs(dx),
    max_y=minmax_y[2] + abs(dy)
)

ga0 = crop(ga, extent)
x0, y0 = GeoArrays.ranges(ga0)
length(x0)
length(y0)

k = nothing

using HTTP

url = "https://s3-eu-west-1.amazonaws.com/download.agisoft.com/gtg/us_nga_egm2008_1.tif"
HTTP.download(url, abspath("junk.tif"))
ga = GeoArrays.read(abspath("junk.tif"); masked=false)

gdl = ArchGDAL.read("junk.tif")
ArchGDAL.gdalinfo(gdl)

b = ArchGDAL.geomdim(gdl)


dataset = ArchGDAL.getgeotransform(ArchGDAL.read("junk.tif"))


b = ArchGDAL.getcoorddim(gdl)
?ArchGDAL.getcoorddim
function ranges(ga::GeoArray, strategy::AbstractStrategy=Center())
    extra = typeof(strategy) == Center ? 0 : 1
    lx, ly = coords(ga, (1, 1), strategy)
    hx, hy = coords(ga, size(ga)[1:2] .+ extra, strategy)
    dx = ga.f.linear[1, 1]
    dy = ga.f.linear[2, 2]
    lx:dx:hx, ly:dy:hy
end






x0, y0 = GeoArrays.ranges(ga)
length(x0)
length(y0)
size(ga)

ga = pointextract(x, y, "EPSG:4326", geoid_to_wgs84)




minmax_x = extrema(x)
minmax_y = extrema(y)
dx = ga.f.linear[1, 1];
dy = ga.f.linear[2, 2];

extent = (
    min_x=minmax_x[1] - abs(dx),
    min_y=minmax_y[1] - abs(dy),
    max_x=minmax_x[2] + abs(dx),
    max_y=minmax_y[2] + abs(dy)
)

ga0 = crop(geoid_to_wgs84, extent)
x0, y0 = GeoArrays.ranges(ga0)
length(x0)
length(y0)
