"""
Calculate and analyze elevation offsets between altimetry data and reference DEMs.

This script calculates horizontal and vertical offsets between altimetry data 
(primarily Hugonnet and ICESat data) and reference DEMs (primarily Copernicus 30m DEM).
It processes data by geotiles and can handle different domains and regions.

Key operations:
1. Sets up project parameters (geotile width, domain, products)
2. Optionally filters geotiles to specific geographic regions
3. Extracts DEM data and masks for each geotile
4. Calculates offsets between altimetry and DEM data
5. Visualizes offset distributions through histograms

Parameters:
- project_id: Project identifier (default: :v01)
- geotile_width: Width of geotiles in degrees (default: 2)
- domain: Processing domain (:landice)
- force_remake: Whether to force remaking of extracted data
- slope: Whether to include slope calculations
- ext: Geographic extent for processing

The script outputs:
- Horizontal (dx, dy) and vertical (dz) offsets
- Histograms comparing offsets against both COP30 and reference heights
- Statistics on offset distributions

Dependencies:
Altim, Extents, Arrow, DataFrames, GeoFormatTypes, Proj, Statistics, Plots
"""

using Altim
using Extents
using Arrow
using DataFrames
using GeoFormatTypes
using Proj 
using Statistics

force_remake = true;
project_id = :v01
geotile_width = 2
domain = :landice

paths = project_paths(project_id = project_id)
geotiles = project_geotiles(; geotile_width = geotile_width, domain = domain)
products = project_products(; project_id=project_id)

#= --------------------------------------------------------------------------
@warn "--- only processing a subset of geotiles ----"
region = :WNA
extent, epsg = region_extent(region);
geotiles = geotile_subset!(geotiles, extent)
dems = [:cop30_v2]; 
# --------------------------------------------------------------------------
=#

slope = true
#missions = [:icesat, :icesat2, :gedi, :hugonnet];

# -> sensor_halfwidt: 1-sigma for gaussian and 1/2 width for average
#sensor_halfwidth = [35 / 4, 11 / 4, 22 / 4, 150 / 2]; #150/2 was selected for hugonnet based on tests
#sensor_kernel = [:gaussian, :gaussian, :gaussian, :average]

#missions = [:hugonnet];
#sensor_halfwidth = [30];
#sensor_kernel = [:average]
dems = [:cop30_v2]

ext = Extent(X=(-126.9, -126.1), Y=(51.1, 51.8))
#ext = Extent(X=(-180, 180), Y=(80, 81))
geotiles = geotile_subset!(geotiles, ext);

vars = [:floatingice, :glacierice, :inlandwater, :land, :landice, :ocean];

df_offsets = DataFrame()
df_offsets[!, :dx] = zeros(size(geotiles.id))
df_offsets[!, :dy] .= 0.
df_offsets[!, :dz] .= 0.

df_offsets[!, :dx0] .= 0.
df_offsets[!, :dy0] .= 0.
df_offsets[!, :dz0] .= 0.

products = project_products(; project_id=project_id);
products = (hugonnet = products.hugonnet,);
#products = (icesat2=products.icesat2,)

@time for (i, geotile) in enumerate(eachrow(geotiles))
    println("---------------------------------------------------------")
    products = project_products(; project_id=project_id);
    products = (hugonnet=products.hugonnet,)
    #products = (icesat2=products.icesat2,);
    
    println(products)

    #hstacks = hstack_catalogue(paths.hugonnet.raw_data; update_catalogue=force_remake)
    #printstyled("building Hugonnet geotiles\n"; color=:blue, bold=true)
    #for geotile in eachrow(geotiles)
    #    geotile_build_hugonnet(geotile, paths.hugonnet.geotile, hstacks; force_remake=force_remake)
    #end

    
    geotile_extract_dem(products, dems, DataFrame(geotile), paths; slope=slope, force_remake=force_remake, xoffset=0, yoffset=0) # TRYING OPOSITE OF RETRIEVED 

    

    for product in products
        geotile_extract_mask(DataFrame(geotile), paths[product.mission].geotile; vars=vars, job_id=product.mission,force_remake=force_remake)
    end
    #=

    df_offsets[i, :dx], df_offsets[i, :dy], df_offsets[i, :dz] = 
    geotile_offset(products, dems, geotile, paths; interations=3, iter_thresh=5)

    df_offsets[i, :dx0], df_offsets[i, :dy0], df_offsets[i, :dz0] = 
        geotile_offset(products, dems, geotile, paths; interations=3, iter_thresh=5, use_href=true)
    println("---------------------------------------------------------")

    =#
    
    #=
    geotile_extract_dem(products, dems, DataFrame(geotile), paths; slope=slope, force_remake=force_remake, xoffset=+50, yoffset=-50) # TRYING OPOSITE OF RETRIEVED 


    for product in products
        geotile_extract_mask(DataFrame(geotile), paths[product.mission].geotile; vars=vars, job_id=product.mission, force_remake=force_remake)
    end

    df_offsets[i, :dx], df_offsets[i, :dy], df_offsets[i, :dz] =
        geotile_offset(products, dems, geotile, paths; interations=3, iter_thresh=5)       
    =#
end

#Arrow.write("/mnt/bylot-r3/data/hugonnet/v001_cop30_v2_offsets.arrow", df_offsets)

using Plots
valid = .!isnan.(df_offsets.dx)
histogram(
    df_offsets.dx[valid]; 
    bins=-130:5:50, 
    label="h - cop30: dx", 
    alpha=0.5, 
    subplot=1, 
    layout=4,
    legend=:topright,
    )

histogram!(df_offsets.dy[valid]; bins=-130:5:50, label="h - cop30: dy", alpha=0.5, subplot=2, legend = :topleft)
histogram!(df_offsets.dx0[valid]; bins=-130:5:50, label="h - h_ref: dx", alpha=0.5, subplot=3, legend=:topleft)
histogram!(df_offsets.dy0[valid]; bins=-130:5:50, label="h - h_ref: dy", alpha=0.5, subplot=4, legend=:topleft)

#=cop30_v2 lat[+50+52]lon[-128-126]: above = -1.09 m, north = -6.0 m, east = 4.0 m
#                                 : std(dh)= 14.03 m, std(dh_shifted)= 15.34 m
# cop30_v2 lat[+50+52]lon[-128-126]: above = -0.32 m, north = -1.0 m, east = -4.0 m
                                 : std(dh)= 12.42 m, std(dh_shifted)= 12.37 m
geotile_offset(products, dems, geotiles, paths)
=#

products = project_products(; project_id=project_id)
products = (hugonnet=products.icesat,)
for geotile in eachrow(geotiles)
    println("-----------------------")
    geotile = DataFrame(geotile)

    geotile_offset(products, dems, geotile, paths; interations=3, iter_thresh = 7)
    geotile_offset(products, dems, geotile, paths; interations=3, iter_thresh=7, use_href=true)
    println("-----------------------")
end
