#rsync -ra bylot:/mnt/bylot-r3/data/height_change/test/ ~/data/height_change/test 
# add https://github.com/JuliaGeo/TileProviders.jl https://github.com/JuliaGeo/MapTiles.jl https://github.com/MakieOrg/Tyler.jl.git

using Altim, DataFrames, Arrow, GLMakie, Statistics, Tyler

tyler = Tyler.Map(Rect2f(-0.0921, 51.5, 0.04, 0.025))

geotile_width = 2 #geotile width [degrees] 
epsgs = (WGS84="EPSG:4326", EGM2008="EPSG:3855", NPS="EPSG:3413")
dem = "arcticdem_v3_10m"; #["cop30_v2", "nasadem_v1", "arcticdem_v3_10m","rema_v2_10m"]; # "cop30_v1", "rema_v2_10m", 
dem = "rema_v2_10m"
hostname = gethostname()
if hostname == "MT-304223"
    data_dir = "/Users/gardnera/data/"
elseif (hostname == "bylot.jpl.nasa.gov") || (hostname == "devon.jpl.nasa.gov")
    data_dir = "/mnt/bylot-r3/data/"
end

# set paths
paths = (; height_change=joinpath(data_dir, "height_change", "test"))

dem = "cop30_v2"
#dem = "rema_v2_10m"
# list of geotile files
geotiles = searchdir(paths.height_change, ".$dem")

extent = (min_x = 102., min_y = -70., max_x = 104., max_y = -64.)


ind = findfirst(geotiles .== "lat[+50+52]lon[+156+158].$dem")
geotiles = geotiles[ind]

df = DataFrame()
@time for geotile in geotiles
    infile = joinpath(paths.height_change, geotile)
    df = append!(df, DataFrame(Arrow.Table(infile)))
end

# find data that falls within the Greenland extents
ind = within.(Ref(extent), df.longitude, df.latitude)
df = df[ind,:];

df, zone, north = itslive_proj!(df, height = 0.)

# bin data
gridsize = 1E3
variable = "trend"; # "amplitude"
begin
    minmax_x = extrema(df.x)
    minmax_y = extrema(df.y)
    xbin_edges= (floor(minmax_x[1]/gridsize)*gridsize - (gridsize/2)):gridsize:(ceil(minmax_x[2]/gridsize)*gridsize + (gridsize/2))
    ybin_edges = (floor(minmax_y[1]/gridsize)*gridsize - (gridsize/2)):gridsize:(ceil(minmax_y[2]/gridsize)*gridsize + (gridsize/2))
    x_binned, y_binned, z_binned, bin_count = bin(df.x, df.y, df[:,"$variable"], xbin_edges, ybin_edges; method = median)
end

# select a map provider
provider = TileProviders.Esri(:WorldImagery)
#provider = TileProviders.Google(:roadmap)

geotile_width = 2 #geotile width [degrees] 
epsgs = (WGS84="EPSG:4326", EGM2008="EPSG:3855", NPS="EPSG:3413")
dem = "arcticdem_v3_10m"; #["cop30_v2", "nasadem_v1", "arcticdem_v3_10m","rema_v2_10m"]; # "cop30_v1", "rema_v2_10m", 
dem = "rema_v2_10m"
hostname = gethostname()
if hostname == "MT-304223"
    data_dir = "/Users/gardnera/data/"
elseif (hostname == "bylot.jpl.nasa.gov") || (hostname == "devon.jpl.nasa.gov")
    data_dir = "/mnt/bylot-r3/data/"
end

# set paths
paths = (; height_change=joinpath(data_dir, "height_change", "test"))

dem = "cop30_v2"
#dem = "rema_v2_10m"
# list of geotile files
geotiles = searchdir(paths.height_change, ".$dem")

# Greenland
extent = (min_x = -70., min_y = 60., max_x = -13., max_y = 83.)

x,y = epsg2epsg(df.longitude, df.latitude, "EPSG:4326", "EPSG:3857")

frame = Rect2f(extent.min_x, extent.min_y, extent.max_x - extent.min_x, extent.max_y - extent.min_y)

# show map
m = Tyler.Map(frame;
    provider, min_tiles=8, max_tiles=16, figure=Figure(resolution=(1000, 600)))
    
# choose color map [https://docs.juliahub.com/AbstractPlotting/6fydZ/0.12.10/generated/colors.html]
cmap = reverse(ColorSchemes.balance.colors)
n = length(cmap);
alpha = ones(n)
nmod = 50;
mod = -nmod:1:nmod
alpha[mod .+ round(Int64,n/2)] = abs.(mod ./ nmod)
cmap = RGBA.(cmap, alpha)

# set color axis limits which by default is automatically equal to the extrema of the color values
colorrange = [-.25, .25];
objscatter = scatter!(m.axis, x, y; color = df.trend, colormap = cmap, colorrange = colorrange, markersize = 20)

# hide ticks, grid and lables
hidedecorations!(m.axis) 

# hide frames
hidespines!(m.axis)