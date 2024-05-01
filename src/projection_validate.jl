using Altim

lats = -90:1+90;
lons = -180:1:180;

for lat in lats
    for lon in lons
        epsg = itslive_epsg(lon, lat, alwaysxy = true)
        zone, north = itslive_zone(lon, lat; alwaysxy = true)

        x, y = AltimPlots.epsg2epsg(lon, lat, "EPSG:4326", ga.crs.val; parse_output=true)

        # having issues getting this fuction to work
        # x,y = lonlat2utm.(lon, lat, getindex(geotile.zone, 1), getindex(geotile.zone, 2))

        #points_lla = LLA.(lat, lon, Ref(0))
        #f = UTMfromLLA(getindex(geotile.zone, 1), getindex(geotile.zone, 2), Geodesy.wgs84)
        #utm = f.(points_lla)
        #x = [a.x for a in utm]
        #y = [a.y for a in utm]
        
        x, y = AltimPlots.epsg2epsg(lon, lat, "EPSG:4326", ga.crs.val; parse_output=true)

        # find xy extents
        minmax_x = extrema(x)
        minmax_y = extrema(y)

        extent = (min_x = minmax_x[1], min_y = minmax_y[1], max_x = minmax_x[2] + abs(ga.f.linear[1,1]), max_y = minmax_y[2] + abs(ga.f.linear[2,2])) 

        ga0 = crop(ga, extent)

        foo = GeoArrays.coords(ga0);
        xy = foo.f.(collect(foo.iter))

        lon, lat = AltimPlots.epsg2epsg(getindex.(xy, 1), getindex.(xy, 2), ga.crs.val, "EPSG:4326"; parse_output=true)

    #lon1, lat1 = utm2lonlat(getindex.(xy, 1), getindex.(xy, 2), getindex(geotile.zone, 1), getindex(geotile.zone, 2))

