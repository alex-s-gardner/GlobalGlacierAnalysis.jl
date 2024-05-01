# geotile_play.jl
using Altim, DataFrames, Statistics, GeoArrays, Arrow, MLJLinearModels
#using ProfileView

folder_height_change = "/mnt/bylot-r3/data/height_change/2003_2022"

#foreach(rm, readdir(folder_height_change))
force_remake = true; 
geotile_width = 2; #geotile width [degrees]
grid = (node_spacing = 50., node_width = 100.); # [(node_spacing = 500, node_width = 500 * 2)]

project_id = :v01;
domain = :landice;
#domain = :all;

paths = project_paths(project_id = project_id);
# add height_change_path
paths = (height_change = folder_height_change, paths...);
products = project_products(; project_id=:v01);

geotiles = project_geotiles(; geotile_width = geotile_width, domain = domain);

parameters = Altim.TSModelParameters(; mask=:landice, crop2mask=true)

# make directly if it doesn't exist
if !isdir(paths.height_change)
    mkpath(paths.height_change);
end

printstyled("calculating height change \n"; color = :blue, bold = true);

#@time asyncmap(eachrow(geotiles[1:20,:]); ntasks = 100) do geotile 
dems = ["cop30_v2", "arcticdem_v3_10m", "cop30_v2", "rema_v2_10m"];

# --------------------------------------------------------------------------
@warn "--- only processing a subset of geotiles ----"
begin
    region = :RGI01;
    ext, epsg = region_extent(region);

    dems = [:cop30_v2]; 
    #products = (; products.icesat, products.icesat2)
    products = (; products.hugonnet)

    #ext = Extent(X=(-126.9, -126.1), Y=(51.1, 51.8));
    #ext = Extent(X=(-141.0, -141.5), Y=(59.0, 59.5))

    geotiles = geotile_subset!(geotiles, ext);

    #geotiles = geotiles[geotiles.id.=="lat[+60+62]lon[-142-140]", :]

end
# --------------------------------------------------------------------------
#using ProfileView
@time begin
    for dem in dems
        for geotile in eachrow(geotiles) #208s, 97s with 8 @Threads, 218s asyncmap, 69s asyncmap & @Threads
            Altim.geotile_ts_fit(geotile, dem, paths, grid; products=products, p=parameters, force_remake=force_remake)
        end
    end
end

begin
    vars = [:h, :inlandwater, :landice, :floatingice, :land, :ocean, :thickness, :region, :vx0, :vy0, :dhdxs, :dhdys];
    for dem in dems
        for geotile in eachrow(reverse(geotiles))
            t1 = time()
            fname = joinpath(paths.height_change, "$(geotile.id).$(dem)")
            if isfile(fname)
                df = Arrow.Table(fname)
                if !isempty(df.latitude)
                    outifle = joinpath(paths.height_change, "$(geotile.id).$(dem)+")
                    df0 = itslive_extract(copy(df.longitude), copy(df.latitude), vars)

                    tmp = tempname(dirname(outifle))
                    Arrow.write(tmp, df0::DataFrame);
                    mv(tmp, outifle; force = true)

                    total_time = round((time()-t1)/60, digits = 3);
                    printstyled("    -> $(geotile.id) $dem height change ancillary data extracted: $(total_time) min \n"; color = :light_black)
                end
            end
        end
    end
end