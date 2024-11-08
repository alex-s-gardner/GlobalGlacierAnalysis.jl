## this script calculate ice dischage automatically... it was perty much done until I moved to using published values... it works well so I may reserect at some point


#begin
# load packages
begin
    using Downloads
    using GDAL
    using Rasters
    import ArchGDAL
    using Shapefile
    using DataFrames
    using Plots
    using Images
    import GeometryOps as GO
    import GeoInterface as GI
    import GeoFormatTypes as GFT
    import GeoDataFrames as GDF
    import DimensionalData as DD
    import GeoDataFrames as GDF
    using GeoJSON
    using Proj
    using Statistics
    using StatsBase
    using FlexiJoins
    using Altim
    using CSV
    
end

#--------------- set local paths --------------- 
paths = Altim.pathlocal
datadir = paths.data_dir;
outshapepath = joinpath(datadir, "global_vector_mask")

# path2shapefile must be manually downloaded... google drive does not provide direct access... place files in outshapepath 
path2shapefile = "https://drive.google.com/file/d/1NeKbo_tiaiyAnd15QLN2mqzt-oPuWan3/view?usp=sharing"
#------------------------------------------------


rgi = 1

varnames = ["not_ocean", "glacier"]

#--------------- set local paths --------------- 
paths = Altim.pathlocal
datadir = paths.data_dir;
outshapepath = joinpath(datadir, "global_vector_mask")

# path2shapefile must be manually downloaded... google drive does not provide direct access... place files in outshapepath 
path2shapefile = "https://drive.google.com/file/d/1NeKbo_tiaiyAnd15QLN2mqzt-oPuWan3/view?usp=sharing"
#------------------------------------------------


rgi = 1
varnames = ["not_ocean", "glacier"]

# get paths to ITS_LIVE component velocities (download locally if it does not exist)
begin
    rasterpaths = Dict()
    rasterpaths[rgi] = Dict()

    for var0  in ["vx", "vy"]
        path2rasterfile = "https://its-live-data.s3.amazonaws.com/velocity_mosaic/v2/static/cog/ITS_LIVE_velocity_120m_RGI$(lpad(string(rgi), 2, '0'))A_0000_v02_$(var0).tif"
        rasterpaths[rgi][var0] = joinpath(datadir, "its-live-data", splitpath(path2rasterfile)[3:end]...)

        if !isdir(joinpath(splitpath(rasterpaths[rgi][var0])[1:end-1]))
            mkpath(joinpath(splitpath(rasterpaths[rgi][var0])[1:end-1]))
        end

        if !isfile(rasterpaths[rgi][var0])
            f = Downloads.download(path2rasterfile, rasterpaths[rgi][var0])
        end
    end

    rasterpaths[rgi]["h"] = paths.cop30_v2
end

# build a vrt of millan thickness data + errors
begin
    thickness_dir = joinpath(datadir, "Millan2021/thickness/RGI-$(rgi)/")
    d = readdir(thickness_dir)
    d = d[startswith.(d, "THICKNESS_RGI").& endswith.(d, ".tif")]
    d = joinpath.(Ref(thickness_dir), d)

    thickness_vrt = joinpath(thickness_dir, "THICKNESS_$(splitpath(thickness_dir)[end]).vrt")

    d = ArchGDAL.read.(d)
    ArchGDAL.gdalbuildvrt(d; dest=thickness_vrt) do vrt
        ArchGDAL.write(vrt, thickness_vrt)
    end
end

# define some helper functions
begin
    function clip_or_empty(polygon, clipper)
        if GO.contains(polygon, clipper) # early return if polygon is completely inside clipper
            return GO.tuples(polygon)::Union{GI.Polygon,GI.MultiPolygon}
        end

        if GO.disjoint(polygon, clipper) # early return if polygon is completely outside clipper
            null_point = GI.is3d(polygon) ? (GI.ismeasured(polygon) ? (NaN, NaN, NaN, NaN) : (NaN, NaN, NaN)) : (NaN, NaN)
            contents = [GI.LinearRing([null_point, null_point, null_point])]
            return GO.tuples(GI.Polygon{GI.is3d(polygon),GI.ismeasured(polygon),typeof(contents),Nothing, typeof(GI.crs(polygon))}(contents, nothing, GI.crs(polygon)))
        end

        result = GO.intersection(polygon, clipper; target = GI.PolygonTrait())
        if isempty(result)
            null_point = GI.is3d(polygon) ? (GI.ismeasured(polygon) ? (NaN, NaN, NaN, NaN) : (NaN, NaN, NaN)) : (NaN, NaN)
            contents = [GI.LinearRing([null_point, null_point, null_point])]
            return GO.tuples(GI.Polygon{GI.is3d(polygon),GI.ismeasured(polygon),typeof(contents),Nothing, typeof(GI.crs(polygon))}(contents, nothing, GI.crs(polygon)))
        else
            return (GI.MultiPolygon(result; crs = GI.crs(polygon)))
        end
        # return GO.intersection(GO.GEOS(), polygon, clipper)
    end

    function bounds2rectangle(xbounds, ybounds)
        rectangle = GI.Wrappers.Polygon([[(xbounds[1], ybounds[1]), (xbounds[1], ybounds[2]), (xbounds[2], ybounds[2]), (xbounds[2], ybounds[1]), (xbounds[1], ybounds[1])]])
        return rectangle
    end

    function area_of_intersection(polya, polyb)
        if GO.intersects(polya, polyb)
            return GO.area(GO.intersection(GO.GEOS(), polya, polyb))
        else
            return 0.0
        end
    end
end

# read in raster data
@time begin
    vx = Raster(rasterpaths[rgi]["vx"])
    
    # find bounds of raster in geographic coodinates
    (xbounds, ybounds) = Rasters.bounds(vx)
    bounds_geom = bounds2rectangle(xbounds, ybounds)
    bounds_geom = GO.segmentize(bounds_geom, max_distance= 1000)
    bounds_geom = GO.reproject(bounds_geom, source_crs = GI.crs(vx), target_crs = GFT.EPSG(4236))
    xbounds_geo = extrema(getindex.(collect(GI.getpoint(bounds_geom)),1))
    ybounds_geo = extrema(getindex.(collect(GI.getpoint(bounds_geom)),2))

    # build back inot a polygon
    bounds_geom = bounds2rectangle(xbounds_geo, ybounds_geo)

    # buffer by 1 degree to make sure nothing is missed
    buff = 1;
    xmin = xbounds_geo[1] .- buff;
    ymin = xbounds_geo[2] .+ buff;
    xmax = ybounds_geo[1] .- buff;
    ymax = ybounds_geo[2] .+ buff;

    vx[abs.(vx).>20000] .= 0
    vy = Raster(rasterpaths[rgi]["vy"])
    vy[abs.(vy).>20000] .= 0

    h = Raster(rasterpaths[rgi]["h"], lazy=true)
    h = Rasters.crop(h, to = bounds_geom)
    h = resample(h; to=vx)
    h[abs.(h).>20000] .= 0

    T = Raster(thickness_vrt, lazy=true)
    T = resample(T; to=vx)

    # extract dimensions of reference dataset
    dx = dims(vx, :X)
    dy = dims(vx, :Y)
end

# rasterize shapefile data
#begin
    # create masks
    mask0 = Dict()
    # USING OGR2OGR [10s]
    @time for varname = varnames
        shapefile =getindex(paths,Symbol("$(varname)_shp"))
        outfile = tempname()*".shp"
        
        GDAL.ogr2ogr_path() do ogr2ogr
            run(`$ogr2ogr -clipsrc $xmin $ymin $xmax $ymax -t_srs EPSG:3413 $outfile $shapefile`)
        end

        feature = DataFrame(Shapefile.Table(outfile))
        mask0[varname] = rasterize!(count, zeros(dx, dy), feature; threaded=false, verbose=false) .> 0
        p = heatmap(mask0[varname]; title=varname)
        display(p)
    end

    # USING GeometryOps
    #=
    @time for varname = varnames
        shapefile = joinpath(outshapepath, "$varname.shp")
        outfile = tempname()*".shp"
        
        feature = DataFrame(Shapefile.Table(shapefile));
        clipper = bounds2rectangle((xmin, xmax), (ymin, ymax))

        glaicer_geom.geometry = clip_or_empty.(glaicer_geom.geometry, (clipper,));
        feature.geometry = GO.reproject.(feature.geometry; source_crs = GFT.EPSG(4326), target_crs = GFT.EPSG(3031), calc_extent=true);

        mask0[varname] = rasterize!(count, zeros(dx, dy),feature; threaded=false, verbose=false) .> 0

        p = heatmap(mask0[varname]; title=varname)
        display(p)
    end

=#

    # encode lake terminating and ocean terminating
    #TODO: Need to automatically locate correct file... maybe just combine all rgi6 files into one or build a vrt?

    geomfile_rgi6 = joinpath(datadir, "GlacierOutlines/rgi60/rgi60_Global.gpkg")
    glaicer_geom = GDF.read(geomfile_rgi6)
    
    # find intersecting geometries
    index = GO.intersects.(Ref(bounds_geom), foo.geom)
    glaicer_geom = glaicer_geom[index,:]

    glaicer_geom.geom = GO.reproject(glaicer_geom.geom; source_crs = GFT.EPSG(4326), target_crs=Rasters.crs(vx), calc_extent=true)

    ocean_terminating = glaicer_geom[glaicer_geom.TermType .== 1 .| glaicer_geom.TermType .== 9, :]
    lake_terminating = glaicer_geom[glaicer_geom.TermType .== 2, :]

    mask0["id"] = rasterize!(last, zeros(dx, dy), feature[(glaicer_geom.TermType.==1).|(glaicer_geom.TermType.==2), :]; threaded=false, verbose=false, fill=:CenLon)
    heatmap(mask0["id"])

    mask0["ocean_terminating"] = rasterize!(count, zeros(dx, dy), ocean_terminating; threaded=false, verbose=false) .> 0
    heatmap(mask0["ocean_terminating"])

    mask0["lake_terminating"] = rasterize!(count, zeros(dx, dy), lake_terminating; threaded=false, verbose=false) .> 0
    heatmap(mask0["lake_terminating"])
end

# define glacier terminus
begin
    # set all velocities outside of lake and ocean terminating glaciers to zero
    vx[.!(mask0["ocean_terminating"] .| mask0["lake_terminating"])] .= 0
    vy[.!(mask0["ocean_terminating"] .| mask0["lake_terminating"])] .= 0

    # identify terminus as lowest x% of ice
    terminus2glaicer_frac = 0.02;
    mask0["lake_terminus"] = falses(dx,dy);
    mask0["ocean_terminus"]= falses(dx,dy);

    valid_term = falses(dx, dy);

    Threads.@threads for id = unique(mask0["id"])
    #id = unique(mask0["id"])[2]

        if id == 0
            continue
        else
            valid_id = id .== mask0["id"];
            term_h = quantile(h[valid_id], terminus2glaicer_frac )
            valid_term[valid_id] .|= (h[valid_id] .<= term_h)
        end
    end

    # buffer slightly to get zero velocity outside of glacier domain
    buffer_pixels = 500 / 120;
    valid_term = Images.distance_transform(feature_transform(valid_term)) .<= buffer_pixels

    y = lookup(vx, Y)[1]:(lookup(vx, Y)[2]-lookup(vx, Y)[1]):lookup(vx, Y)[end];
    x = lookup(vx, X)[1]:(lookup(vx, X)[2]-lookup(vx, X)[1]):lookup(vx, X)[end];

    pol = GO.polygonize(x, y, valid_term.data);

    # polygonize returns a multipart polygon for each unique value in the raster... so in 
    # this case pol is a single multipart polygon with each polygon being a connected 
    # cluster of pixels. I want the polygons as sperated geometires.. so us GI.getgeom
    pol = DataFrame(geometry=GI.getgeom(pol))
end

# migrate attibutes and glacier plygons to terminus plygons
begin
    # do a left join to on the polygonize features 
    df = leftjoin((pol, feature), by_pred(:geometry, GO.intersects, :geometry));
    df[:,:delete] .= false;
    rename!(df, :geometry_1 => :geometry_glacier)

    # some left polygons intersect more than one right polygon...
    gdf = groupby(df, :geometry)

    # loop though each group, if a group has more than one feature, select the feature with the 
    # largest overlap
    for df in gdf
        #df = gdf[4]
        if nrow(df) > 1
            area = zeros(nrow(df))
            for (i, r) in enumerate(eachrow(df))
                area[i] = area_of_intersection(r.geometry, r.geometry_glacier)
            end
            df.delete = area .< maximum(area)
        end
    end

    df = reduce(vcat,gdf)
    df = df[.!df.delete,:]

    df[!, :geometry] = [GI.MultiPolygon([r.geometry]) for r in eachrow(df)]

    # check if a glacier has more than one plygon
    gdf = groupby(df, :RGIId)
    for df in gdf
    #df = gdf[7]
        if nrow(df) > 1
            df[1,:geometry] = GI.MultiPolygon(reduce(vcat,GI.getgeom.(df.geometry)))
            df.delete[2:end] .= true
        end
    end
    df = reduce(vcat, gdf)
    df = df[.!df.delete, :]
    select!(df, Not(:delete))
end

# calculate flux
begin
    df[:, :glacier_flux] .= 0.;
  
    for r in eachrow(df)
        pts = GI.Point.(GI.getpoint(r.geometry))
        foo = collect(Rasters.extract(vx, pts))
        vx0 = getindex.(foo, :var"")
        foo = collect(Rasters.extract(vy, pts))
        vy0 = getindex.(foo, :var"")

        foo = collect(Rasters.extract(T, pts))
        T0 = getindex.(foo, :var"")

        pts = [[p.geom[1], p.geom[2]] for p in pts]
        pts = vcat(pts, pts[1,:])
        Δx = getindex.(pts[2:end,:], 1) .- getindex.(pts[1:end-1,:], 1)
        Δy = getindex.(pts[2:end,:], 2) .- getindex.(pts[1:end-1], 2)

        # this is flux into a polygon that has points defined with clockwise ordering
        fx = vx0 .* Δy .* T0;
        fy = vy0 .* -Δx .* T0;

        r.glacier_flux = (sum(fx) .+ sum(fy)) ./ 1E9
    end
end

# reporject to 4326
begin
    df.geometry = GI.MultiPolygon.(df.geometry, crs=GFT.EPSG(4326))
    df.geometry_glacier = GI.MultiPolygon.(df.geometry_glacier, crs=GFT.EPSG(4326))
    df.geometry = GI.MultiPolygon.(df.geometry, crs=GFT.EPSG[3413])
    df.geometry_glacier = GI.MultiPolygon.(df.geometry_glacier, crs=GFT.EPSG[3413])
end

# write to file
begin
    rename!(df, :glacier_flux => :FLUX)
    filename = replace(splitpath(shapefile_rgi6)[end], ".shp" => "_flux.shp")
    filename = joinpath("/Users/gardnera/Downloads/", filename)
    Shapefile.write(filename, df; force=true)
end
end