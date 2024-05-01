var = :glacier
fileIn = "/Users/gardnera/data/GlacierOutlines/WorldMask_20190513/$(var)_0.0083.tif"
fileOut = joinpath(pwd(), "data", "$(var)_0.0083_cog.tif")

# add GeoArrays#master
using GeoArrays

ga0 = GeoArrays.read(fileIn);

# convert to a UInt8
ga = GeoArray(UInt8.(ga0.A))
ga.crs = ga0.crs
ga.f = ga0.f


# plot data
using GLMakie
fig = Figure()
ax = Axis(fig[1, 1])
hm = heatmap!(ax, A[1:100:end,1:100:end,1])
Colorbar(fig[1, 2], hm)
fig

# write as a cog
GeoArrays.write(fileOut, ga, shortname="COG");