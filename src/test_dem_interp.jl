using Altim, Statistics, GeoArrays, DataFrames, Arrow, CairoMakie

geotile_width = 2 #geotile width [degrees]
force_remake = true;
sensor = :GEDI

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

# geotiles = DataFrame(geotiles[findfirst(contains.(geotiles.id, "lat[+76+78]lon[-082-080]")), :]) # Devon Island
geotiles = DataFrame(geotiles[findfirst(contains.(geotiles.id, "lat[+30+32]lon[+078+080]")), :]) # HMA

# geotiles = geotiles[ind,:]

#dems = ["cop30_v2", "nasadem_v1","rema_v2_10m","arcticdem_v3_10m"]; 
dems = ["nasadem_v1"];
slope = true
# ICESat1
begin
    sensor_footrpint = 35; # in meters
    for dem in dems
        geotile_extract_dem(geotiles, paths.icesat.geotile, dem; sensor_footrpint = sensor_footrpint, slope = slope, force_remake = force_remake);
    end

    fn = joinpath(paths.icesat.geotile, geotiles.id[1]*".arrow");
    df = DataFrame(Arrow.Table(fn));
    h1 = vcat(df[:, :height]...);

    fn = joinpath(paths.icesat.geotile, geotiles.id[1]*"."*dems[1]);
    df = DataFrame(Arrow.Table(fn));
    
    h2 = vcat(df[:, :height]...);

    dh = h1 .- h2;
    valid = abs.(dh) .< 300;
    std(h1[valid] .- h2[valid])

    # HMA: lat[+30+32]lon[+078+080]
    # cop30_v2 - w filter
    # 18.03 = nearest
    # 15.52 = linear
    # 15.46 = Quadratic

    # nasadem_v1 - w/o filter
    # 16.82 = nearest
    # 15.86 = linear
    # 15.98 = Quadratic

    # nasadem_v1 - w filter
    # 16.81 = nearest
    # 15.86 = linear
    # 15.96 = Quadratic

    # Ellesmere Island: lat[+76+78]lon[-082-080]
    # arcticdem_v3 - w filter
    # 6.26 = nearest
    # 6.82 = linear
    # 5.77 = Quadratic

    # arcticdem_v3 - w/o filter
    # 4.67 = nearest 
    # 4.70 = linear
    # 7.07 = Quadratic

    # cop30_v2 - w filter
    # 5.687 = nearest
    # 5.37 = linear
    # NaN = linear

    # cop30_v2 - w/o filter
    # 5.71 = nearest
    # 5.40 = linear
    # NaN = linear
end

# ICESat2
begin
    sensor_footrpint = 11; # in meters
    for dem in dems
        geotile_extract_dem(geotiles, paths.icesat2.geotile, dem; sensor_footrpint = sensor_footrpint, slope = slope, force_remake = force_remake)
    end


    fn = joinpath(paths.icesat2.geotile, geotiles.id[1]*".arrow");
    df = DataFrame(Arrow.Table(fn));
    h1 = vcat(df[:, :height]...);

    fn = joinpath(paths.icesat2.geotile, geotiles.id[1]*"."*dems[1]);
    df = DataFrame(Arrow.Table(fn));
    h2 = vcat(df[:, :height]...);


    dh = h1 .- h2;
    valid = abs.(dh) .< 300;
    std(h1[valid] .- h2[valid])

    # HMA: lat[+30+32]lon[+078+080]
    # 24.62 = nearest
    # 23.95 = linear
    # 24.05 = cubic
    # 24.05 = Quadratic
end

# GEDI
begin
    sensor_footrpint = 22; # in meters
    for dem in dems
        geotile_extract_dem(geotiles, paths.gedi.geotile, dem; sensor_footrpint = sensor_footrpint, slope = slope, force_remake = force_remake)
    end

    fn = joinpath(paths.gedi.geotile, geotiles.id[1]*".arrow");
    df = DataFrame(Arrow.Table(fn));
    h1 = vcat(df[:, :height]...);

    fn = joinpath(paths.gedi.geotile, geotiles.id[1]*"."*dems[1]);
    df = DataFrame(Arrow.Table(fn));
    h2 = vcat(df[:, :height]...);


    dh = h1 .- h2;
    valid = abs.(dh) .< 300;
    std(h1[valid] .- h2[valid])

    # cop30_v2 - w filter
    # 20.56 = nearest
    # 18.65 = linear
    # 18.63 = Quadratic

    # nasadem_v1 - w filter
    # 19.66 = nearest
    # 18.96 = linear
    # 19.06 = Cubic 
    # 19.05 = Quadratic

    # with sensor filter
end


