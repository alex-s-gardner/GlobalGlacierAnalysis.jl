

using FileIO
import GlobalGlacierAnalysis as GGA
import GeometryOps as GO
using DataFrames
using CairoMakie
using Statistics
using Dates
using NCDatasets


h1 = "/mnt/bylot-r3/data/hugonnet/HSTACK/001/raw/13_14_15_rgi60/stacks/44N/N30E078.nc"
h2 = "/mnt/bylot-r3/data/hugonnet/HSTACK/001/raw_filtered/13_14_15_rgi60/stacks/44N/N30E078_prefilt.nc"

h1_mtime = Dates.unix2datetime(mtime(h1))
h2_mtime = Dates.unix2datetime(mtime(h2))

h1 = NCDataset(h1);
h2 = NCDataset(h2);

h1["z"]
h2["z"]

h1["time"][1]
h2["time"][1]

z1 = h1["z"]
z2 = h2["z"]

function range_overlap(x1, x2)
    x_min = max(minimum(x1), minimum(x2))
    x_max = min(maximum(x1), maximum(x2))

    range_overlap1 = (x1 .>= x_min) .& (x1 .<= x_max)
    range_overlap2 = (x2 .>= x_min) .& (x2 .<= x_max)

    return range_overlap1, range_overlap2
end

ind_x1, ind_x2 = range_overlap(h1["x"], h2["x"]);
ind_y1, ind_y2 = range_overlap(h1["y"], h2["y"]);

z1 = h1["z"][ind_x1, ind_y1, 1];
z2 = h2["z"][ind_x2, ind_y2, 1];

scatter(vec(z1), vec(z2))





GGA.geotiles_golden_test[i]

reference_run = "binned/2deg/glacier_dh_best_nmad5_v01_filled_ac_p2_aligned.jld2"
path2reference = joinpath(GGA.pathlocal[:data_dir], reference_run)

foo = load(path2reference)

geotiles = GGA.geotiles_mask_hyps(:glacier, 2)

i = 1;
area = round(Int, sum(geotiles.glacier_area_km2[findfirst(geotiles.id .== GGA.geotiles_golden_test[i])]))

println("$(GGA.geotiles_golden_test[i]) contains $(area) km2 of glaciers")

i = 2;
sum(geotiles.glacier_area_km2[findfirst(geotiles.id .== GGA.geotiles_golden_test[i])])
println("$(GGA.geotiles_golden_test[i]) contains $(sum(geotiles.glacier_area_km2[findfirst(geotiles.id .== GGA.geotiles_golden_test[i])])) km2 of glaciers")


geotiles = geotiles[sum.(geotiles.glacier_area_km2).>0.0, :]

glacier_dh_best_cc_nmad5_v01.jld2
path2reference = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_dh_best_nmad5_v01_filled_ac_p1_aligned.jld2"
foo = load(path2reference)

Dates.unix2datetime(mtime(path2reference))

df = foo["model_param"]["hugonnet"];
i = 1;
offset = round(df.offset[findfirst(df.geotile .== GGA.geotiles_golden_test[i])], digits=1)
println("$(GGA.geotiles_golden_test[i]) has an ASTER offset of $(offset) m")

df[findall(abs.(df.offset) .> 20), :geotile]


index = [findfirst(geotile .== df.geotile) for geotile in geotiles.id]
hist(df.offset)

median(df.offset)
std(df.offset)

# make lat lon scatter plot
pts = GO.centroid.(GGA.geotile2rectangle.(df.geotile))

# set color axis limits which by default is automatically equal to the extrema of the color values
colorrange = [-15, 15];

# choose color map [https://docs.juliahub.com/AbstractPlotting/6fydZ/0.12.10/generated/colors.html]
cmap = :balance;

# make colored scatter plot
f = Figure();
scatter(f[1, 1], getindex.(pts, 1), getindex.(pts, 2); color=df.offset_icesat, colormap=cmap, markersize=10, strokewidth=0, colorrange=colorrange)

Colorbar(f[1, 2], limits=colorrange, colormap=cmap, flipaxis=false, label="ASTER minus ICESat (m)")

f


egm2008 = GGA.geoid("egm2008"; folder=GGA.pathlocal.geoid_dir)
fig, ax, hm = heatmap(-egm2008[1:10:end, 1:10:end,1], colormap=cmap, colorrange=colorrange)
scatter!(ax, getindex.(pts, 1), getindex.(pts, 2); color=df.offset_icesat, colormap=cmap, markersize=10, strokewidth=0, colorrange=colorrange)
Colorbar(fig[:, end+1], hm)
fig