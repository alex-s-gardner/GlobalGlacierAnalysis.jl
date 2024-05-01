# this is 30 times slower than simple sign fit. 


# geotile_play.jl
using Altim, DataFrames, Statistics, GeoArrays, Arrow, Extents, EasyModelAnalysis

folder_height_change = "/mnt/bylot-r3/data/height_change/WNA/2018_2022";

force_remake = true; 
geotile_width = 2; #geotile width [degrees]
altim_dem_outlier = 200;
t_minmax = (2018.,  2022.);                  # only incude data aquired between min & max date
grid = (node_spacing = 250, node_width = 250 * 2); # [(node_spacing = 500, node_width = 500 * 2)]


const y_intercept = 2020.0; # [2015]

# see https://math.stackexchange.com/questions/2430564/equation-of-a-tilted-sine
function model(t, p)
    amp_scale = 1 / ((1 - (sqrt(1.0 - p[5]^2) + 2 * sqrt(1.0 - p[5]^4)) / 3) * (pi / 2 - 1) + 1)
    t_cycle = 2.0 * pi * (t .+ p[4])
    y = p[1] .+ p[2] * (t .- y_intercept) .+ p[3] * (amp_scale / p[5]) * atan.(p[5] * sin.(t_cycle) ./ (1 .- p[5] * cos.(t_cycle)))
    return y
end

p = [0.0, 0.0, 100, 0.5, 0.5]
p_min = [-20000., -2000., 0., 0., 0.001]
p_max = [20000., 2000., 5000., 1., 0.999]

lsq_setup = Altim.LsqSetup(model, p, p_min, p_max, :forwarddiff)

ts_model_arg = (
    count_min = 5,                 # minimum number of valid values for model fit [15]
    t_min_range = 3.,               # only file model in range(t) > t_min_range [10]
    z_max_outlier = 150.,           # maximum acceptable absolute deviation from median
    t_bin_edges = 2002:1/6:2023,    # bin edges used for local filtering of outliers
    t_bin_filter_method = madnorm,  # method for identifying local outliers
    t_bin_filter_threshold = 5,     # threshold for classifying local outliers
);

project_id = :v01;
geotile_width = 2;
domain = :landice;

paths = project_paths(project_id = project_id);
# add height_change_path
paths = (height_change = folder_height_change, paths...);
geotiles = project_geotiles(; geotile_width = geotile_width, domain = domain);

# make directly if it doesn't exist
if !isdir(paths.height_change)
    mkpath(paths.height_change);
end

# loop for each geotile

#point_data_file = joinpath(paths.height_change, "$(geotile_id(extent)).$(dem)")
#if isfile(point_data_file) && !force_remake
#    df = DataFrame(Arrow.Table(point_data_file)) # 20s to load lazilly
#else

printstyled("calculating height change \n"; color = :blue, bold = true);

#@time asyncmap(eachrow(geotiles[1:20,:]); ntasks = 100) do geotile 
dems = [:cop30_v2] #, :nasadem_v1, :rema_v2_10m, :arcticdem_v3_10m];

# --------------------------------------------------------------------------
# @warn "--- only processing a subset of geotiles ----"
# region = :WNA;
# extent, epsg = region_extent(region);
# geotiles = geotile_subset!(geotiles, extent);
# dems = [:cop30_v2];
# delete old files
#run(`/bin/bash rm -r /mnt/bylot-r3/data/height_change/WNA/2018_2022/*`)

# --------------------------------------------------------------------------

for dem in dems
    for geotile in eachrow(geotiles) #208s, 97s with 8 @Threads, 218s asyncmap, 69s asyncmap & @Threads
        Altim.geotile_ts_fit(geotile, dem, paths, lsq_setup, grid, ts_model_arg;
            altim_dem_outlier=altim_dem_outlier, t_minmax=t_minmax, mask=:landice, force_remake = force_remake)
    end
end

vars = [:h, :inlandwater, :landice, :floatingice, :land, :ocean, :thickness, :region, :vx0, :vy0, :dhdxs, :dhdys];
for dem in dems
    for geotile in eachrow(reverse(geotiles))
        t1 = time()
        fname = joinpath(paths.height_change, "$(geotile.id).$(dem)")
        if isfile(fname)
            df = Arrow.Table(fname)
            if !isempty(df.latitude)
                outifle = joinpath(paths.height_change, "$(geotile.id).$(dem)+")
                df0 = itslive_extract(copy(df.longitude), copy(df.latitude), vars)

                tmp = tempname(dirname(outifle))
                Arrow.write(tmp, df0::DataFrame);
                mv(tmp, outifle; force = true)

                total_time = round((time()-t1)/60, digits = 3);
                printstyled("    -> $(geotile.id) $dem height change ancillary data extracted: $(total_time) min \n"; color = :light_black)
            end
        end
    end
end