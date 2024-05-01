# fast index = true
#0.996989 seconds (29.39 k allocations: 678.275 MiB, 38.24% gc time, 3.98% compilation time)

# Rasters v0.10.1
#9.457382 seconds (115.50 M allocations: 5.703 GiB, 7.77% gc time, 0.24% compilation time)

lon, lat = X(25:.1:30), Y(25:.1:30);
mask0 = Raster(falses(lon, lat));

lon = rand(25:0.1:30, 100);
lat = rand(25:0.1:30, 100);

pts1 = ([((y, x)) for (x, y) in zip(lon, lat)]);
pts2 = extract(mask0, pts1, atol=0.0001, index=true, geometry=false);
