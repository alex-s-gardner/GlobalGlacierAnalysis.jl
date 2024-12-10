# extract raw data for a given geotile and dump statistics


begin
mission = "hugonnet"
surface_mask = :landice
project_id = :v01
geotile_width = 2
dh_max = 200;

dem_id = :cop30_v2
using Altim
using Arrow
using DataFrames
using Statistics

geotiles = Altim.geotiles_w_mask(geotile_width)
paths = Altim.project_paths(; project_id)
products = project_products(; project_id)
binned_folder = analysis_paths(; geotile_width).binned

product = products[Symbol(mission)];

#idx = (geotiles.rgi9 .> 0) .& (geotiles.glacier_frac .> 0.1)
#Threads.@threads for geotile in eachrow(geotiles)


# <><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><><>
#for geotile in eachrow(geotiles[idx,:])
#geotile = geotiles[findfirst((geotiles.rgi2 .> 0.) .& (geotiles.glacier_frac .> 0.1)),:];
#geotile = geotiles[findfirst(geotiles.id .== "lat[+28+30]lon[+082+084]"), :]

#geotile = geotiles[findfirst(geotiles.id .== "lat[-72-70]lon[+160+162]"), :]
geotile = geotiles[findfirst(geotiles.id .== "lat[+44+46]lon[+006+008]"), :]; # used for curvature figure with mission ="icesat"

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

#if geotile.glacier_frac == 0.0
#    continue
#end

t1 = time();

mission_geotile_folder = paths[Symbol(mission)].geotile 

# special case for unfiltered hugonnet data

# this is a hack for including 2 folders that contain hugonnet data (filtered and unfiltered)
if contains(binned_folder, "unfiltered") && (mission == "hugonnet")
    mission_geotile_folder = replace(mission_geotile_folder, "/2deg" => "/2deg_unfiltered")
end

path2altim = joinpath(mission_geotile_folder, geotile.id * ".arrow");

#if !isfile(path2altim)
#    continue
#end

altim = select!(DataFrame(Arrow.Table(path2altim)), [:longitude, :latitude, :datetime, :height, :quality]);

altim = Altim.add_dem_ref!(altim, dem_id, geotile, mission_geotile_folder)

path2masks = joinpath(mission_geotile_folder, geotile.id * ".masks");
masks0 = select!(DataFrame(Arrow.Table(path2masks)), [:inlandwater, :land, :landice, :ocean]);

valid = (masks0[:,surface_mask]) .& (.!isnan.(altim.dh))
valid[valid] = (abs.(altim.dh[valid] .- median(altim.dh[valid])) .< dh_max) .& (abs.(altim.dh[valid]) .< dh_max * 2)

Altim.std(altim.dh[valid])

println("$(mission) $(geotile.id) $(surface_mask): 
    median = $(round(median(altim.dh[valid]), digits = 2))
    mean = $(round(mean(altim.dh[valid]), digits = 2))
    mad = $(round(Altim.mad(altim.dh[valid]), digits = 2))
    std = $(round(std(altim.dh[valid]), digits = 2))
       ")
end