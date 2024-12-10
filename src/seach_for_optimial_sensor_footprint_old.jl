using Altim
using Extents

force_remake = true;
project_id = :v01
geotile_width = 2
domain = :landice

paths = project_paths(project_id = project_id)
geotiles = project_geotiles(; geotile_width = geotile_width, domain = domain)

dems = [:cop30_v2, :nasadem_v1, :rema_v2_10m, :arcticdem_v3_10m]; 
# --------------------------------------------------------------------------
@warn "--- only processing a subset of geotiles ----"
region = :WNA
extent, epsg = region_extent(region);
geotiles = geotile_subset!(geotiles, extent)
dems = [:cop30_v2];
# --------------------------------------------------------------------------
#

slope = true
missions = [:icesat, :icesat2, :gedi, :hugonnet];

# -> sensor_halfwidt: 1-sigma for gaussian and 1/2 width for average
sensor_halfwidth == [35/4, 11/4, 22/4, 100/2]; 
sensor_kernel = [:gaussian, :gaussian, :gaussian, :average]

missions = [:hugonnet];
sensor_halfwidth = [100/2];
sensor_kernel = [:average]
dems = [:cop30_v2]

using Arrow, DataFrames, Arrow, Statistics, Plots
dem = first(dems)
mission = first(missions)
geotile = first(geotiles)

sensor_halfwidth = 50;
sensor_kernel = :average;
#sensor_kernel = :gaussian;
geotiles0 = DataFrame(geotiles[4,:];)


for sensor_halfwidth in [30, 100, ]/2; #, 100, 120, 140, 150, 160, 180, 200]/4

    geotile_extract_dem(geotiles0, paths[mission].geotile, dem;
        filter_halfwidth=sensor_halfwidth, filter_kernel=sensor_kernel, slope=slope,
        job_id=mission, force_remake=force_remake)

    path2geotile = joinpath(paths[mission].geotile, geotiles0[1, :id] * ".arrow")

    df1 = DataFrame(Arrow.Table(path2geotile));
    df2 = DataFrame(Arrow.Table(replace(path2geotile, ".arrow" => ".$dem")));
    df3 = DataFrame(Arrow.Table(replace(path2geotile, ".arrow" => ".masks")));

    # plot off ice histoframs
    h1 = vcat(df1.height...)
    h2 = vcat(df1.height_reference...)
    h3 = vcat(df2.height...)
    valid = .!(vcat(df3.landice...))

    dh1 = h1[valid] .- h2[valid];
    dh2 = h1[valid] .- h3[valid];

    valid2 = (abs.(dh1) .< 200) .& (abs.(dh2) .< 200)
    dh1 = dh1[valid2]
    dh2 = dh2[valid2]
    println("$(sensor_halfwidth) m, STD = $(round(std(dh2)/std(dh1)*100))%")
end

#histogram(dh1, bins = -100:10:100)
#histogram!(dh2, bins=-100:10:100)