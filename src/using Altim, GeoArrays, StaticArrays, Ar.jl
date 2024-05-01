using Altim, GeoArrays, StaticArrays, Arrow, DataFrames

geotile_width = 2 #geotile width [degrees] 
force_remake = false;

hostname = gethostname()
if hostname == "MT-304223"
    data_dir = "/Users/gardnera/data/"
elseif (hostname == "bylot.jpl.nasa.gov") || (hostname == "devon.jpl.nasa.gov")
    data_dir = "/mnt/bylot-r3/data/"
end

paths = (
    icesat2=setpaths(data_dir, geotile_width, :ICESat2, "ATL06", "005"),
    icesat=setpaths(data_dir, geotile_width, :ICESat, "GLAH06", "034"),
    gedi=setpaths(data_dir, geotile_width, :GEDI, "GEDI02_A", "002")
)

geotiles = geotile_define(geotile_width)

#geotiles = DataFrame(geotiles[findfirst(contains.(geotiles.id, "lat[+76+78]lon[-082-080]")), :]) # Devon Island
#geotiles = DataFrame(geotiles[findfirst(contains.(geotiles.id, "lat[+30+32]lon[+078+080]")), :]) # HMA
geotiles = DataFrame(geotiles[findfirst(contains.(geotiles.id, "lat[+44+46]lon[+006+008]")), :]) # Alps


sensor = :icesat2

fn = joinpath(paths[sensor].geotile, geotiles.id[1] * ".arrow")
df0 = DataFrame(Arrow.Table(fn))[!,[:longitude, :latitude]]

vars = [:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean]

lon = vcat(df0.longitude...)
lat = vcat(df0.latitude...)

@time df0 = itslive_extract(lon, lat, vars, path2param = "/mnt/devon-r2/data/its-live-data/autorift_parameters/v001/")

df = DataFrame()
vlength = length.(df0.longitude)
for v in vars
    df[!,v] = [Array{Bool}(undef,l) for l in vlength] 
end  

start = 1
for row = eachrow(df)
    stop = length(row[vars[1]]) + start - 1;
    for v in vars
        row[v] = df0[!,v][start:stop]
    end
    start = stop + 1
end
