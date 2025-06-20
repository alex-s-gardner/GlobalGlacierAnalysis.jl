import GlobalGlacierAnalysis as GGA
using DimensionalData
using CairoMakie


geotiles2plot = ["lat[+28+30]lon[+082+084]"]
dheight = Dim{:height}(50:100:10000)
date_range, date_center = GGA.project_date_bins()
ddate = Dim{:date}(date_center)
dgeotile = Dim{:geotile}(geotiles2plot)
mission = ["icesat2", "icesat", "hugonnet", "gedi"]
dh = Dict(m => rand(dgeotile, ddate, dheight) for m in mission)
area_km2 = rand(dheight)


dh0 = Dict(mission => dh[mission][geotile=At("lat[+28+30]lon[+082+084]")] for mission in keys(dh))


f = GGA.plot_elevation_time(dh0["icesat2"]; colorrange=(-20, 20))

#GGA.plot_elevation_time_multimission(dh0; colorrange=(-20, 20), linkaxes=true, colorbar_label="height anomaly [m]")






