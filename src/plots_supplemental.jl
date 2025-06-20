

using FileIO
import GlobalGlacierAnalysis as GGA
import GeometryOps as GO
using DataFrames
using CairoMakie
using Statistics

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


geotiles = geotiles[sum.(geotiles.glacier_area_km2).>1, :]

path2reference = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_dh_best_nmad5_v01_filled_ac_p1_aligned.jld2"
foo = load(path2reference)

mission = "hugonnet"
df = foo["model_param"][mission];

index = [findfirst(geotile .== df.geotile) for geotile in geotiles.id]
df = df[index,:]
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