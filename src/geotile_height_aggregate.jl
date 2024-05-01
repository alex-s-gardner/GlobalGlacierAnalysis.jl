```
aggregate elevation data from multiple altimetry missions and save to file
```
using Altim
using DataFrames
using Arrow

geotile_width = 2; #geotile width [degrees]
altim_dem_outlier = 200;
t_minmax = nothing;

project_id = :v01;
geotile_width = 2;
domain = :landice;

paths = project_paths(project_id = project_id);
geotiles = project_geotiles(; geotile_width = geotile_width, domain = domain);

# make directly if it doesn't exist
# if !isdir(paths.height_change)
#     mkpath(paths.height_change);
# end

#@time asyncmap(eachrow(geotiles[1:20,:]); ntasks = 100) do geotile 
dem = [:cop30_v2, :nasadem_v1, :rema_v2_10m, :arcticdem_v3_10m];
# missions = [:gedi, :icesat, :icesat2, :hugonnet]
products = project_products(; project_id=:v01)

# --------------------------------------------------------------------------
@warn "--- only processing a subset of geotiles ----"
region = :WNA;
extent, epsg = region_extent(region);
geotiles = geotile_subset!(geotiles, extent);
dem = :cop30_v2
mask = :landice
crop2mask = false

ext = Extent(X=(-128, -124), Y=(50, 52))
geotiles = geotile_subset!(geotiles, ext);


# delete old files
# run(`/bin/bash rm -r /mnt/bylot-r3/data/height_change/WNA/2018_2022/*`)

outfile = "/mnt/bylot-r2/shared_data/altim/$(region)_$(dem)_$(mask).arrow"

geotiles = geotiles[end:end,:]
# --------------------------------------------------------------------------
df = DataFrame();
for geotile in eachrow(geotiles)
    append!(df, geotile_merge_height(geotile, dem, paths; products, altim_dem_outlier=altim_dem_outlier, t_minmax=t_minmax, mask=mask, crop2mask=crop2mask))
end

# save as an arrow file
Arrow.write(outfile, df::DataFrame)

using Plots
valid = .!(df.landice) .& (df.product .== "I206") #.& .!(df.height .== df.height_reference0)
histogram(df.height[valid] .- df.height_reference[valid], linecolor = nothing)
histogram!(df.height[valid] .- df.height_reference0[valid], linecolor=nothing)