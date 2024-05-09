using Altim
using Rasters
using Shapefile
using DataFrames 
using Arrow

geotile_width = 2;
raster_file = :cop30_v2; #, :nasadem_v1, :rema_v2_10m, :arcticdem_v3_10m]
mask = :glacier; #, :nasadem_v1, :rema_v2_10m, :arcticdem_v3_10m]

runid = "geotile_$(mask)_hyps"
binned_folder = analysis_paths(; geotile_width).binned
out_file = joinpath(binned_folder, "$(runid).arrow");

dbin = 100.;
bin_centers = dbin/2:dbin:10000-dbin/2;

out_var_name = Symbol("$(mask)_area_km2")
shp = Symbol("$(mask)_shp");
geotiles = Altim.geotiles_w_mask(geotile_width);

fn_raster = Altim.pathlocal[raster_file];
fn_shp = Altim.pathlocal[shp];

feature = Shapefile.Handle(fn_shp);
ras = Raster(fn_raster; lazy=true);

geotiles[!, out_var_name] = [zeros(size(bin_centers)) for r in 1:nrow(geotiles)];

mask_frac = Symbol("$(mask)_frac")
Threads.@threads for geotile in eachrow(geotiles)

    if (geotile[mask_frac] > 0.0)

        t1 = time();
        ras0 = read(Rasters.crop(ras, to=geotile.extent));

        mask0 = Rasters.rasterize!(count, zeros(ras0.dims), feature; threaded=false, shape=:polygon, progress=false, verbos=false) .> 0;

        # calculate area per cell
        lon = lookup(ras0, X)
        lat = lookup(ras0, Y)
        d = Altim.meters2lonlat_distance.(Ref(1), lat)
        a = abs.((1 ./ getindex.(d, 2) * (lat[2] .- lat[1])) .* (1 / d[1][1] * (lon[2] - lon[1])))
        area_m2 = repeat(a', outer = [length(lon), 1])

        feature_area_km2 = mask0 .* area_m2 / (1000 * 1000)

        foo = geotile[out_var_name]
        for (i, cntr) in enumerate(bin_centers)
            bs = cntr - dbin / 2
            be = cntr + dbin / 2
            foo[i] = sum(feature_area_km2[(ras0.>=bs).&(ras0.<be)]) 
        end
        
        t = round(time() - t1)
        printstyled("    -> $(geotile.id) hypsometry calculated: $(t)s\n"; color=:light_black)
    end
end

Arrow.write(out_file, geotiles::DataFrame)