begin
using Altim
using NCDatasets
using Rasters
using Statistics
using GeoTiles
using Shapefile
using DataFrames
using Extents
using Plots
using CSV
using NearestNeighbors

include("utilities_hyps.jl")

project_id = :v01;
geotile_width = 2;

g = 9.80665;

gemb_folder = "/home/schlegel/Share/GEMBv1/"

surface_mask = :glacier

pathlocal = setpaths()
fn_geopotential = "/mnt/bylot-r3/data/era5/Geopotential.nc"
fn_precip = "/mnt/bylot-r3/data/era5/total_precip_monthly.nc"
fn_landsea = "/mnt/bylot-r3/data/era5/LandSeaMask.nc"

# DO NOT MODIFY: these need to match elevations in hypsometry
Δh = 100.;
height_range = 0.:Δh:10000.;
height_center = height_range[1:end-1] .+ Δh / 2;

precip = Raster(fn_precip; name = :tp)
geopotential = Raster(fn_geopotential; name=:z)
landseamask = Raster(fn_landsea; name=:lsm)

# expver = 1 is ERA5 and expver = 5 is ERA5T (low latency data that will be updated)
idx = ismissing.(precip[1, 1, At(1), :]);
precip[:, :, At(1), idx] = precip[:, :, At(5), idx]
precip = precip[:, :, At(1), :]
precip = mean(precip, dims=3)

precip = mean(precip, dims=3) .* (365.25)# in units of meters per day
precip = dropdims(precip, dims=:Ti)

geopotential = dropdims(geopotential, dims=:Ti)
height = geopotential ./ g

landseamask = dropdims(landseamask, dims=:Ti)

# find points near ice
lon = val(lookup(precip, :X))
lat = val(lookup(precip, :Y))
lon[lon.>180] = lon[lon.>180] .- 360

n = length(lat)
m = length(lon)
lon = lon * ones(1, n)
lat =  ones(m, 1) * lat'

mask = falses(size(lon))

foo = Altim.itslive_extract(lon[:], lat[:], [:glacierice])
mask[:] = foo.glacierice

surface_mask0 = "glacier_b10km"
shp = Symbol("$(surface_mask0)_shp");
fn_shp = Altim.pathlocal[shp]
feature = Shapefile.Handle(fn_shp)

Δlat = abs(lon[1] - lon[2])
mask[:] = Altim.highres_mask(lat[:], lon[:], feature; grid_resolution=Δlat)

# remove points where points are over ocean
mask[landseamask.data .== 0] .= false

df = DataFrame(latitude=lat[mask], longitude=lon[mask], era5_height=height[mask], era5_precip=precip[mask])

region = :rgi6_regions
shp = Symbol("$(region)_shp")
fn_shp = pathlocal[shp]
feature = Shapefile.Table(fn_shp)

ids = unique(feature.RGI_CODE)
rgi_ids = "rgi" .* string.(ids)

for (id, rgi_id) in zip(ids, rgi_ids)
    idx = feature.RGI_CODE .== id
    df[!, rgi_id] = Altim.highres_mask(df.latitude, df.longitude, feature.geometry[idx]; grid_resolution=Δlat)
end


# check that everyingthing looks correct
# heatmap(height)
# heatmap(precip, clims=(0, 15))
end

# find subset of tiles with glacier ice
geotile_min_area = 0
geotiles = geotiles_mask_hyps(surface_mask, geotile_width)
geotiles = geotiles[sum.(geotiles.glacier_area_km2).>0, :]
ntile_total = nrow(geotiles)
area_total = sum(sum(geotiles.glacier_area_km2))
geotiles = geotiles[sum.(geotiles.glacier_area_km2).>geotile_min_area, :]
ntile_subset = nrow(geotiles)
area_subset = sum(sum(geotiles.glacier_area_km2))
geotiles, reg = geotiles_mutually_exclusive_rgi!(geotiles)

println("$ntile_subset of $ntile_total geotiles with > $geotile_min_area km2 include $(round(Int64,(area_subset / area_total)*100)) % of total ice coverage")

# internal function
_inbin(x, bins) = (x .> bins[1:end-1]) .& (x .<= bins[2:end])

df[!, "geotile"] .= ""
for r in eachrow(geotiles)
    # move max value slightly down to avoid overlap between geotiles
    ext = r.extent
    delta = 0.00001
    x = ext.X
    y = ext.Y

    ext = Extent((X=(x[1], x[2] - delta), Y=(y[1], y[2] - delta)))
    in_geotile = Altim.within.(Ref(ext), df.longitude, df.latitude)
    df[in_geotile, "geotile"] .= r.id
end

# remove points that do not fall within a geotile (probably edge points)
df = df[.!(df[!, "geotile"].==""), :]

begin
    height_delta0 = [
        [0],
        [-200, 0, 200],
        [-200, 0, 200, 500, 1000]
        ]

    min_area_improve = 1.
    df_vec = [];

    buffer = 1;

    for height_delta in height_delta0
        geotiles[!, :era5_nodes] .= 0
        geotiles[!, :era5_height_extrema] .= Ref((NaN,NaN))
        geotiles[!, :era5_precip_extrema] .= Ref((NaN, NaN))
        geotiles[!, :era5_height_bins] .= Ref(falses(size(height_center)))

        for r in eachrow(geotiles)

            # move max value slightly down to avoid overlap between geotiles
            ext = Extents.buffer(r.extent, (X=buffer, Y=buffer))
            delta = 0.00001
            x = ext.X
            y = ext.Y

            ext = Extent((X=(x[1], x[2] - delta), Y=(y[1], y[2] - delta)))

            in_geotile = Altim.within.(Ref(ext), df.longitude, df.latitude)

            if !any(in_geotile)
                continue
            end

            h = reduce(vcat, [df[in_geotile, :era5_height] .+ hd for hd in height_delta])

            r.era5_nodes = sum(in_geotile)
            r.era5_height_extrema = extrema(h)
            r.era5_precip_extrema = extrema(df[in_geotile, :era5_precip])

            r.era5_height_bins = reduce(.|,   _inbin.(h, Ref(height_range)))
        end

        minimum_glacier_area_per_elevation = 0.01
        extrema_exclude(x; minarea=0.01, hrange=height_center) = any(x .> minarea) ? extrema(hrange[x.>minarea]) : (0, 0)
        geotiles[!, :glacier_height_range] = extrema_exclude.(geotiles.glacier_area_km2)

        geotiles[!, :missing_low_elevation] = max.(getindex.(geotiles.era5_height_extrema, 1) .- getindex.(geotiles.glacier_height_range, 1), Ref(0))
        geotiles[!, :missing_high_elevation] = max.(getindex.(geotiles.glacier_height_range, 2) .- getindex.(geotiles.era5_height_extrema, 2), Ref(0))

        plot(geotiles[!, :missing_high_elevation]; label="missing_high")
        plot!(geotiles[!, :missing_low_elevation]; label="missing_low")

        # sampled area
        geotiles[!, :era5_sampled_area_km2]= sum.(map(.*, geotiles.era5_height_bins, geotiles.glacier_area_km2))

        # unsampled glacier area
        geotiles[!, :era5_unsampled_area_km2] = sum.(map(.*, map(.!, geotiles.era5_height_bins), geotiles.glacier_area_km2))

        # % unsampled
        sa = sum(geotiles[!, :era5_sampled_area_km2])
        ua = sum(geotiles[!, :era5_unsampled_area_km2])
        pu = round(Int64, ua /(sa +ua)*100)
        println("% area without FAC = $pu%")

        push!(df_vec, copy(geotiles))
    end

    for i = 1:length(height_delta0)-1
        idx = (df_vec[i+1].era5_sampled_area_km2 .- df_vec[i].era5_sampled_area_km2) .> min_area_improve
        println("geotiles with more than $(min_area_improve) km² covered with height_delta: $(sum(idx))")
    end
end


#df = select!(df, Not(["rgi_id"]))
rgi_names = names(df)
rgi_names = rgi_names[contains.(rgi_names, "rgi")]
rgi_id = reduce(hcat, eachcol(df[:, rgi_names]))
rgi_id = findfirst.(eachrow(rgi_id))
df[!, "rgi_id"] = rgi_id
# remove old ids
df = select!(df, Not(rgi_names))

#tree = KDTree(hcat(df.longitude, df.latitude); leafsize=10)
valid = .!isnothing.(df.rgi_id);
tree = KDTree(hcat(df.longitude[valid], df.latitude[valid])'; leafsize=10)
idx, dis = nn(tree, hcat(df.longitude[.!valid], df.latitude[.!valid])')
df[.!valid, "rgi_id"] = df[idx, "rgi_id"]
df.rgi_id = Int64.(df.rgi_id)

# groupings
run_groups = [[1,2], [3], [4], [5], [6, 7, 9], [8, 10, 11, 12], [13], [14,15], [16, 17, 18], [19]]
df[!,"run_group"] .= 0;
for i in eachindex(run_groups)
    rg = run_groups[i]
    idx = falses(size(df.rgi_id))
    for g in rg
        idx = idx .| (df.rgi_id .== g)
    end
    
    df[idx,"run_group"] .= i
end

df[!,:northern_hemisphere] .= df.latitude .>= 0

gdf = DataFrames.groupby(df, :run_group)
DataFrames.combine(gdf, nrow)

output = joinpath(gemb_folder, "global_glacier_nodes_v1.csv")
CSV.write(output, df)

plot(df.longitude, df.latitude; seriestype=scatter)
plot(df.longitude, df.latitude; seriestype=scatter)
