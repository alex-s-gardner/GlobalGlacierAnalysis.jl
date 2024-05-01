using Altim
using GeoTiles
using Plots
using BinStatistics
using DataFrames
using Statistics
#Pkg.add("OnlineStats")

project_id = :v01;
paths = project_paths(project_id = project_id);
products = project_products(; project_id=project_id);

(extent, epsg) = Altim.region_extent(:RGI02)

p1 = plot();
xlabel!("curvature");
ylabel!("elevation anomaly [altim - dem]");

bins = -1:0.05:1
for mission in (:icesat, :icesat2, :gedi, :hugonnet)
    folder = paths[mission].geotile

    # find geotiles with a region region
    suffix = ".cop30_v2";
    dem = GeoTiles.listtiles(folder; suffix, extent);
    suffix = ".arrow";
    altim = GeoTiles.listtiles(folder; suffix, extent);
    suffix = ".masks";
    masks = GeoTiles.listtiles(folder; suffix, extent);

    # find geotiles that exisit for all datasets 
    ids = intersect(dem.id, altim.id, masks.id);
    dem = dem[findall(in(ids),dem.id), :];
    altim = altim[findall(in(ids),altim.id), :];
    masks = masks[findall(in(ids),masks.id), :];

    # read matching geotiles
    altim = reduce(vcat, GeoTiles.read.(altim.path2file));
    dem = reduce(vcat, GeoTiles.read.(dem.path2file));
    masks = reduce(vcat, GeoTiles.read.(masks.path2file));

    dh = reduce(vcat, altim.height .- dem.height);
    #dh = reduce(vcat, is2.height .- is2.height_reference)
    # plot histogram
    ice = reduce(vcat,masks.landice);

    foo = Altim.meters2latlon_distance.(Ref(1), reduce(vcat, altim.latitude)); # degree per meter
    lat_dist = getindex.(foo, 1);
    lon_dist = getindex.(foo, 2);

    curvature = -2(reduce(vcat, dem.dhddx) .* lon_dist.^2 + reduce(vcat, dem.dhddy) .* lat_dist.^2) * 100;
    height = reduce(vcat, dem.height);

    #meters2latlon_distance(distance_meters, latitude_degrees)
    df = DataFrame(; dh, ice, curvature, height);

    idx = .!isnan.(dh) .& .!isnan.(curvature) .& .!ice
    df = df[idx,:]
    
    bs = binstats(df, :curvature, bins, :dh; col_function=[mean, median, std])
    scatter!(BinStatistics.bincenter.(bs.curvature), bs.dh_median, label= "$(mission)")
end

#=
p2 = plot();
xlabel!("curvature");
ylabel!("frequency");

p2 = histogram!(curvature[.!ice];  normalize=:pdf, color=:gray);
xlims!(-1,1)

plot(p, p2, layout=(2, 1), legend=false, size=(500, 500))
=#
