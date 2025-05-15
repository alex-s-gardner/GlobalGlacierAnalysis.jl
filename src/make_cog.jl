"""
    make_cog.jl

Convert a GeoTIFF file to a Cloud Optimized GeoTIFF (COG).

This script:
1. Reads a glacier mask GeoTIFF file
2. Converts the data to UInt8 format
3. Preserves the coordinate reference system and transformation
4. Visualizes a downsampled version of the data
5. Writes the result as a Cloud Optimized GeoTIFF

COGs allow for more efficient access to geospatial data, especially when
stored in cloud environments, by enabling access to specific portions of
the data without downloading the entire file.
"""

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