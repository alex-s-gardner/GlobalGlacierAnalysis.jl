using Shapefile
using DataFrames
using Altim
using GeoTiles
using Extents 
using Statistics

test_pts_shpfile = "data/test_modelfit.shp"
folder_height_change = "/mnt/bylot-r3/data/height_change/2000_2022"
project_id = :v01;
domain = :landice;
#domain = :all;
geotile_width = 2

paths = project_paths(project_id=project_id);
# add height_change_path
paths = (height_change=folder_height_change, paths...);
products = project_products(; project_id=:v01);
geotiles = project_geotiles(; geotile_width=geotile_width, domain=domain);
pts = DataFrame(Shapefile.Table(test_pts_shpfile))


pts_lon = [p.x for p in pts.geometry];
pts_lat = [p.y for p in pts.geometry];

# buffer extent
i = 1;
pts_extent = Extent(Lat=extrema(pts_lat[i]), Lon=extrema(pts_lon[i]))
pts_extent_buff = Extents.buffer(pts_extent, (Lat=lat_deg, Lon=lon_deg))
buff_meters = 1000;
lat_deg, lon_deg = Altim.meters2lonlat_distance(buff_meters, mean(pts_lat))

gt = GeoTiles.readall(folder_height_change; suffix=".cop30_v2", extent=pts_extent_buff)


# load in ttest points
test_pts = DataFrame(Shapefile.Table(test_pts_shpfile))

# find geotile for each point
