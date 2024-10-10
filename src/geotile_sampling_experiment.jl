# makes figure incuded in 2024 ROSES Cryopshere propsal
using Altim
using Rasters
using CairoMakie
using GeoDataFrames
import GeometryOps as GO
using DataFrames
using BinStatistics
using Statistics
using DataInterpolations
using FileIO
using DimensionalData

path2geotile_dhdt = "/mnt/bylot-r3/data/hugonnet/2000-01-01_2020-01-01/test/hugonnet_dhdt.tif"
paths = Altim.pathlocal

geotiles = Altim.GeoTiles.define(2)

index = findfirst(Altim.within.(geotiles.extent, Ref(-125), Ref(51)))
geotile = geotiles[index,:]

dhdt = Raster(path2geotile_dhdt)
h = Raster(paths.cop30_v2, lazy=true)

h = Rasters.crop(h; to=geotile.geometry)

h0 = Rasters.resample(h, to=dhdt)
heatmap(h0)


shp = GeoDataFrames.read(paths.glacier_shp)

index = GO.intersects.(shp.geometry, Ref(geotile.geometry))
shp = shp[index,:]

source_crs
geom = GO.reproject(shp[:, :geometry]; source_crs=crs(shp), target_crs=crs(h0))

ice_mask = rasterize!(count, zeros(Int16, dims(h0)), geom) .>0
heatmap(ice_mask)


dims(h0, X)[2] - dims(h0, X)[1]
valid = ice_mask .& (h0 .> 0) .& (dhdt .!= 0)

df = DataFrame(dhdt=dhdt[valid], h=h0[valid], area_km2=0.010)

h_bins = 200:100:3300

df1 = binstats(df, :h, h_bins, [:dhdt, :area_km2], col_function=[median, sum], missing_bins=true)

total_area = df1.area_km2_sum;
total_dv = sum(total_area .* df1.dhdt_median) /1000

bar(df1.h, df1.area_km2_sum)


sample_size = Int.(nrow(df)/2:-1000:10)
iter = 100

dv = fill(NaN, iter)
df2 = DataFrame(;sample_size)
df2[!,:dv] = [fill(NaN, iter) for r in 1:nrow(df2)]

n = nrow(df);
for r in eachrow(df2)
#r = first(eachrow(df2))
    for i = 1:iter
        rand_index = round.(Int,rand(r.sample_size)*(n-1).+1)
        df1 = binstats(df[rand_index,:], :h, h_bins, :dhdt, col_function=[median], missing_bins=true)

        dhdt0 = df1.dhdt_median
        valid = .!ismissing.(dhdt0)
        if .!all(valid)
            interp = DataInterpolations.LinearInterpolation(df1.h[valid], dhdt0[valid]; extrapolate=true)
            dhdt0[.!valid] .= interp(df1.h[.!valid])
        end

        r.dv[i] = sum(dhdt0 .*total_area ./1000)
    end
end

r = first(eachrow(df2))
p = plot(Float64.(r.sample_size * ones(length(r.dv))), r.dv; seriestype=:scatter, legend=false)
for r in eachrow(df2)
    plot!(Float64.(r.sample_size * ones(length(r.dv))), r.dv; seriestype=:scatter)
end

p
ylims!(-2.5, -2.)


# cnvert to % error in retrieved volume change
for r in eachrow(df2)
    r.dv = (r.dv .- total_dv) ./ total_dv
end
df2[!,:dv_median] = median.(df2.dv)*100
df2[!,:dv_std] = Altim.mad.(df2.dv)*100


f = Figure(fontsize = 20)
ax = Axis(f[1, 1]; ylabel="% error in volume change", xlabel="% of data sampled")

μ = df2.dv_median;
t = df2.sample_size ./ nrow(df) * 100
σ = df2.dv_std
lines!(t, μ)              # plot mean line
band!(t, μ + σ, μ - σ)   # plot stddev band
Makie.ylims!(ax,-1,1)

f


# what is the arverage number of ICESat points in a geotile
project_id = "v01"
paths = Altim.project_paths(; project_id)
products = project_products(; project_id)



binned_folder = "/mnt/bylot-r3/data/binned/2deg"
surface_mask = :glacier
binning_method = "meanmadnorm3"
dem_id = :best
curvature_correct = true


binned_file = Altim.binned_filepath(; binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct)

nobs0 = load(binned_file, "nobs_hyps")
nobs_count = sum(nobs0["icesat"], dims=(:date, :height))

d = readdir(binned_folder);
index = findfirst(contains.(d, "glacier_dh_best_cc_meanmadnorm3_v01_filled_ac_p1_synthesized.jld2"))
d = d[index]

dh0 = load(joinpath(binned_folder, d), "dh_hyps")


geotiles = Altim.geotiles_mask_hyps(:glacier, 2)
geotile = "lat[+28+30]lon[+082+084]"

index = findfirst(geotiles.id .== geotile)
dh0 = dh0[At(geotile),:,:]


area = geotiles[index,:].glacier_area_km2

# average surface height
sum(area .* dims(dh0, :height)) ./ sum(area)

dh0[isnan.(dh0)] .= 0
f = sum((hcat(ones(size(dh0, 1))) * area') .* dh0; dims = :height)/1000