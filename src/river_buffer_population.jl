begin
    using ProgressMeter
    using Altim
    import GeoFormatTypes as GFT
    using Proj
    import GeometryOps as GO
    import GeoInterface as GI
    using CairoMakie
    using LibGEOS
    using Rasters, ArchGDAL
    using JLD2
    using FileIO
    #using Altim
    using ArchGDAL
    using GeoDataFrames
    using DataFrames
    #using RangeExtractor
    using GeometryBasics
    using DimensionalData.Lookups

    using BinStatistics
    using LibGEOS # was not able to use GeometryOps

    # Load local configuration paths
    paths = Altim.pathlocal

    # exclude rivers at high latitudes and near the dateline
    latlim = 65.0

    # path to flux data
    rivers_path = joinpath(paths.data_dir, "rivers/MERIT_Hydro_v07_Basins_v01")
    glacier_rivers_path = joinpath(rivers_path, "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.arrow")
    glacier_rivers_population_path = joinpath(rivers_path, "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier_population.arrow")


    glacier_rivers_runoff_path = replace(glacier_rivers_path, ".arrow" => "_runoff.arrow")
    glacier_melt_rivers_runoff_qout_path = replace(glacier_rivers_runoff_path, ".arrow" => "_qout.arrow")

end

# define a setcrs function
begin
    setcrs(x, crs) = setcrs(GI.trait(x), x, crs)

    function setcrs(trait::Nothing, x, crs)
        error("Trait not found for $(typeof(x))")
        # or, use GeometryOpsCore at the top level
    end

    function setcrs(trait::GI.AbstractGeometryTrait, x, crs)
        return GI.Wrappers.geointerface_geomtype(trait)(x; crs)
    end


    function setcrs(geom::Missing, crs)
        return missing
    end
end

# build a vrt of global high resolution population data
begin
    hrsl_dir = "/mnt/bylot-r3/data/dataforgood-fb-data/hrsl-cogs/hrsl_general/v1.5/"
    foo = splitpath(hrsl_dir)
    hrsl_vrt = joinpath(foo[1:end-1]..., foo[end] * ".vrt")
    if .!isfile(hrsl_vrt)
        hrsl_files = Altim.allfiles(hrsl_dir; fn_endswith = ".tif")

        hrsl = ArchGDAL.read.(hrsl_files)
        ArchGDAL.gdalbuildvrt(hrsl; dest=hrsl_vrt) do vrt
            ArchGDAL.write(vrt, hrsl_vrt)
        end
    end
    ras1 = Raster(hrsl_vrt; lazy=true);
end

# load in rivers, trim to lat limits, group by major basin
begin
    # load rivers
    rivers = GeoDataFrames.read(glacier_rivers_path)[:, ["geometry", "COMID"]]

    # exclude rivers at high latitudes and near the dateline
    latlim = 65.
    midlatitude_poly = GI.Polygon(GI.LinearRing(
        [
            (-180., -latlim),
            (-180., latlim),
            (180., latlim),
            (180., -latlim),
            (-180., -latlim)
        ]
    ))

    rivers[!, :midlatitude] = GO.intersects.(getindex.(collect.(GI.getpoint.(rivers.geometry)),1), Ref(midlatitude_poly))

    # exclude rivers at high latitudes and near the dateline
    rivers = rivers[rivers[:, :midlatitude], :]

    # loop over major basins
    rivers[!, :majorbasin_id] = floor.(Int16, rivers.COMID ./ 1e6)

    # group by major basin
    gdf = groupby(rivers, :majorbasin_id)

    # sanity check
    #=
    p = plot(rivers[(majorbasin_id .== uid[1]) .& (rivers[!, :midlatitude]), :geometry])
    for u in uid[2:end]
        index =(majorbasin_id .== u) .& (rivers[!, :midlatitude]);
        if any(index)
            plot!(rivers[index, :geometry])
        end
    end
    p
    =#
end


buffer_radii = [1_000, 5_000, 25_000, 50_000] #meters

for buffer_radius in buffer_radii
    pop_var = "population_$(round(Int16, buffer_radius/1000))km"
    println(pop_var)
    rivers[!, pop_var] = Array{Union{Float64,Missing}}(missing, nrow(rivers))
end


# This code processes river data in batches by major basin:
# 1. For each group of rivers:
#    - Projects river geometries from WGS84 to EASE-Grid 2.0 (EPSG:6933)
#    - Creates 5km buffers around each river segment
#    - Projects buffers back to WGS84
#    - Converts geometries to LibGEOS format for spatial operations
#    - Handles overlapping buffers by:
#      * Finding intersections between current and previous geometries
#      * Taking differences to avoid double-counting areas
#      * Maintaining a union of all processed geometries
#    - Sets CRS and handles invalid/empty geometries
#    - Extracts population data from high-res raster within buffers
#    - Stores population sums for each river segment

if !isfile(glacier_rivers_population_path)
    for (k, rivers) in enumerate(gdf)
        println("#### EXTRACTING POPULATION ARROUND RIVER BUFFER: $k of $(length(gdf))")

        if k < 20
            continue
        end

        for buffer_radius in buffer_radii

            if (k == 20) && (buffer_radius < 25_000)
                continue
            end

            pop_var = "population_$(round(Int16, buffer_radius/1000))km"
            #rivers = gdf[4]

            #GeoDataFrames.write("test_rivers.geojson", rivers)

            # buffer rivers
            # project to ease grid
            source_crs1 = GFT.EPSG(4326)
            target_crs1 = GFT.EPSG(6933)
            trans = Proj.Transformation(source_crs1, target_crs1; always_xy=true)

            #Threads does not help here
            geometry_buffered = []
            @showprogress dt = 1 desc = "buffering rivers..." for r in eachrow(rivers)
                #r = eachrow(rivers)[2]
                foo = GO.transform(trans, r.geometry)
                foo_buff = GO.buffer(foo, buffer_radius)
                foo_buff = GO.transform(inv(trans), foo_buff)
                push!(geometry_buffered, foo_buff)
            end

            # sanity check
            #p = plot(geometry_buffered[1])
            #plot!.(geometry_buffered[2:end])
            #p
            #

            # convert geometry to LibGEOS
            geom = []
            for (i, g) in enumerate(geometry_buffered)
                try
                    geom = push!(geom, GI.convert(LibGEOS, g))
                catch e
                    geom = push!(geom, missing)
                end
            end

            # remove invalid geometries
            for (i, g) in enumerate(geom)
                if !ismissing(g)
                    if !LibGEOS.isValid(g)
                        geom[i] = missing
                    end
                end
            end

            #An array of all previous geometries, as polygons
            g0 =  geom[1]
            geometry_buffered_cropped = [];
            push!(geometry_buffered_cropped, g0)
            #p = plot(g0[1])

            gdiff = g0;
            gint = g0;

            #@showprogress dt = 1 desc = "buffering rivers..." for (i, g) in enumerate(geom[2:end])
            @showprogress dt = 1 desc = "find valid buffered geometries..." for (i, g) in enumerate(geom[2:end])

            # i = 1
            # g = geom[i+1]

                # println(i)

                # find intersection of geometry with all previous geometries
                # there are edge cases where the geometry is not valid
                try
                    gint = LibGEOS.intersection(g0, g)
                catch e
                    println(e)
                    gint = missing
                end

                # there are edge cases where the difference is not valid
                try
                    gdiff = LibGEOS.difference(g, gint)
                catch e
                    println(e)
                    gdiff = missing
                end

                #plot!(gdiff[1])
                push!(geometry_buffered_cropped, gdiff)

                # union with all previous geometries
                # there are edge cases at dateline and poles
                try
                    g0 = LibGEOS.union(g0, g)
                catch e
                    println(e)
                end
            end

            geom = geometry_buffered_cropped

            # remove invalid geometries
            for (i, g) in enumerate(geom)
                if !ismissing(g)
                    if !LibGEOS.isValid(g)
                        geom[i] = missing
                    end
                end
            end

            geom = setcrs.(geom, Ref(GFT.EPSG(4326)));
            geom = convert(Vector{Union{Missing, eltype(geom)}}, geom)

            # find empty geometries and replace with missing [empty geometries excite all sorts of cryptic errors]
            for (i, g) in enumerate(geom)
                if !ismissing(g)    
                    if GI.isempty(g.geom)
                        geom[i] = missing
                    end
                end
            end

            # Load raster dataset
            ras = read(Rasters.crop(ras1; to=skipmissing(geom)))

            # Define tiling scheme
            rivers[.!ismissing.(geom), pop_var] = Rasters.zonal(sum, ras; of=skipmissing(geom), threaded=false)
        end
    end
    GeoDataFrames.write(glacier_rivers_population_path, rivers)
else
    rivers = GeoDataFrames.read(glacier_rivers_population_path)
end

# load in river glacier runoff data
begin
    glacier_flux = GeoDataFrames.read(glacier_melt_rivers_runoff_qout_path)
    glacier_flux = leftjoin(glacier_flux, rivers[!, Not(:geometry)], on=:COMID)
end


# calculate population as a function of maximum glacier fraction and buffer radius
begin
    buffer_radii = Sampled([1, 5, 25, 50]; order=ForwardOrdered(), span=Irregular(0, 50), sampling=Intervals(End()))

    glacier_fractions = Sampled(0.0:0.1:0.8; order=ForwardOrdered(), span=Irregular(0., 1.), sampling=Intervals(Start()))


    glacier_fractions = Dim{Symbol("gmax")}(glacier_fractions)
    buffer_radii = Dim{Symbol("distance from river (km)")}(buffer_radii) #kilometers
    population = zeros(glacier_fractions, buffer_radii)
    for gf in glacier_fractions
        for br in buffer_radii
            pop_var = "population_$(br)km"

            valid = .!ismissing.(glacier_flux[!, pop_var]) .& .!ismissing.(glacier_flux[!, :glacfrac_max]) .& (glacier_flux[!, :glacfrac_max] .>= gf)
            population[At(gf), At(br)] = sum(glacier_flux[valid, pop_var])
        end
    end


    fig, ax, hm = CairoMakie.heatmap(population; colorscale=log10, colormap=:dense)
    CairoMakie.Colorbar(fig[1, 2], hm; label= "population")
    fig
end



gmax = 0
buff = 50
println("buffer = $(buff)km, gmax = $(gmax): population = $(round(Int, population[At(gmax), At(buff)] / 1E6)) million")


gmax = 0.3
buff = 50
println("buffer = $(buff)km, gmax = $(gmax): population = $(round(Int, population[At(gmax), At(buff)] / 1E6)) million")