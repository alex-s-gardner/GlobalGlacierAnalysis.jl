function hstack_catalogue(hstack_parent_dir; update_catalogue = false)
    outfile = joinpath(hstack_parent_dir,"hstack_catalogue.arrow")

    if !isfile(outfile) || update_catalogue

        # find files
        files = allfiles(hstack_parent_dir; subfolders=true, fn_endswith=".nc")

        # build a geospatial catalogue of
        hstacks = DataFrame(path=files);
        hstacks[!, :ProjString] .= Ref(GFT.ProjString("+proj="));
        hstacks[!, :x] .= Ref(0:1:0);
        hstacks[!, :y] .= Ref(0:1:0);
        hstacks[!, :time_range] .= Ref((Dates.DateTime(0), Dates.DateTime(0)));
        hstacks[!, :nlayers] .= 0;
        hstacks[!, :keys] .= Ref(["time", "x", "y", "crs", "dem_names", "z", "uncert", "ref_z", "corr"])

        # this take 12s
        for (i, file) in enumerate(hstacks.path)
        
            hstacks[i, :ProjString] = convert(GFT.ProjString, GFT.WellKnownText(ncgetatt(file, "crs", "spatial_ref")))

            ##----------------------------------------
            # edgest => centers [this was an issue with v1 hstacks]
            if !contains(hstack_parent_dir, "HSTACK/001")
                error("hugonnet verson 1 data had a 1/2 pixel offset that is corrected for here, check different versions before useing")
            end

            x = ncread(file, "x")
            hstacks[i, :x] = (x[1]:(x[2]-x[1]):(x[1]+(x[2]-x[1])*(length(x)-1))) .+ 50

            y = ncread(file, "y")
            hstacks[i, :y] = (y[1]:(y[2]-y[1]):(y[1]+(y[2]-y[1])*(length(y)-1))) .- 50
            ##----------------------------------------

            t = Date(1900, 1, 1) .+ Day.(ncread(file, "time"))
            t_ext = extrema(t)
            if any(t_ext .< Date(2000, 1, 1)) || any(t_ext .> Date(2025, 1, 1))
                error("data is being read incorrectly")
            end
            hstacks[i, :time_range] = t_ext;

            hstacks[i, :nlayers] = length(t)

            # neet to fix this once I figure out how.
            ### stacks[i, :keys] = = NetCDF.ncinfo(file)
        end

        # find bounding extents in lat lon
        # this take 15s

        hstacks[!, :GeoExtent] .= Ref(Extent(X=(0.0, 0.0), Y=(0.0, 0.0)));

        for dfr = eachrow(hstacks)
            tran = Proj.Transformation(dfr.ProjString.val, "EPSG:4326")

            # create a point bounding box
            a = tran.(collect(dfr.x), Ref(maximum(dfr.y)))
            b = tran.(collect(dfr.x), Ref(minimum(dfr.y)))
            c = tran.(Ref(maximum(dfr.x)), collect(dfr.y))
            d = tran.(Ref(minimum(dfr.x)), collect(dfr.y))
            a = vcat(a,b,c,d)
            #dfr.proj2geo = Proj.Transformation( dfr.ProjString.val,"EPSG:4326")
            dfr.GeoExtent = Extent(X = (minimum([x[2] for x in a]), maximum([x[2] for x in a])), Y = (minimum([x[1] for x in a]), maximum([x[1] for x in a])))
        end

        tmp = tempname(dirname(outfile))
        Arrow.write(tmp, hstacks::DataFrame)
        mv(tmp, outfile; force = true)
    end

    hstacks = DataFrame(Arrow.Table(outfile))

    # add back types
    ps = [h.ProjString.val for h in eachrow(hstacks)]
    hstacks[!, :ProjString] = [GFT.ProjString(p) for p in ps]

    ex = [h.GeoExtent.bounds for h in eachrow(hstacks)]
    hstacks[!, :GeoExtent] = [Extent(x) for x in ex]
    #hstacks[!, :proj2geo] .= [Proj.Transformation(hstack.ProjString.val, "EPSG:4326") for hstack in eachrow(hstacks)]

    return hstacks
end

function hstacks2geotile(geotile, hstacks; old_format=false)
    # find intersecting tiles
    stacksind = findall(Extents.intersects.(hstacks.GeoExtent, Ref(geotile.extent)))

    gt = DataFrame()

    for sind in stacksind

        # extract data
        x = collect(hstacks[sind, :x]) * ones(1, length(hstacks[sind, :y])) 
        y = ones(length(hstacks[sind, :x])) * collect(hstacks[sind, :y])'

        proj2geo = Proj.Transformation(hstacks[sind, :ProjString].val, "EPSG:4326")
        lat_lon = proj2geo.(x, y)

        isin = map(x -> within(geotile.extent, x[2], x[1]), lat_lon)
        
        if !any(isin)
            continue
        else
            ll = lat_lon[isin]
            ind = findall(isin)

            # build geotile [98% of all time is spent reading in data in this block]
            begin
                if old_format
                    # old unfilered data 
                    t =  Date(2000, 1, 1) .+ Day.(ncread(hstacks[sind, :path], "time"))
                    z = NetCDF.open(hstacks[sind, :path], "z") do v v[ind, :]::Matrix{Float32} end
                    ref_z = NetCDF.open(hstacks[sind, :path], "ref_z") do v v[ind]::Vector{Float32} end
                    uncert = ncread(hstacks[sind, :path], "uncert")::Vector{Float32}
                    dem_names = ncread(hstacks[sind, :path], "dem_names")::Vector{String}
                    corr = NetCDF.open(hstacks[sind, :path], "corr") do v v[ind, :]::Matrix{Int8} end
                else
                    # new unfilered data 
                    t =  Date(1900, 1, 1) .+ Day.(ncread(hstacks[sind, :path], "time"))
                    z = ncread(hstacks[sind, :path], "z")[ind, :]::Matrix{Float32}
                    ref_z = ncread(hstacks[sind, :path], "ref_z")[ind]::Vector{Float32}
                    uncert = ncread(hstacks[sind, :path], "uncert")::Vector{Float32}
                    dem_names = ncread(hstacks[sind, :path], "dem_names")::Vector{String}
                    corr = ncread(hstacks[sind, :path], "corr")[ind, :]::Matrix{UInt8}
                end

            end

            begin
                gt0 = DataFrame()

                gt0[!, :datetime] = t
                gt0[!, :granule_info] = [(dem_name=dem_names[i], uncert=uncert[i]) for i = 1:nrow(gt0)]

                f = (x, valid) -> vec(x[valid])
                latitude = Float32.([l[1] for l in ll])
                longitude = Float32.([l[2] for l in ll])

                gt0[!, :height] .= [f(z[:, i], .!isnan.(z[:, i])) for i = 1:nrow(gt0)]
                gt0[!, :latitude] .= [f(latitude, .!isnan.(z[:, i])) for i = 1:nrow(gt0)]
                gt0[!, :longitude] .= [f(longitude, .!isnan.(z[:, i])) for i = 1:nrow(gt0)]

                @warn(" hugonnet quality flag currently excludes all ArcticDEM data")
                # isaster .= getindex.(altim.id, 1) .== 'A'
                # isarcticdem .= getindex.(altim.id, 1) .== 'S'
                gt0[!, :quality] .= [f(corr[:, i] .>= 70, .!isnan.(z[:, i])) for i = 1:nrow(gt0)]
                gt0[!, :height_error] .= uncert
                gt0[!, :height_reference] .= [f(ref_z[:], .!isnan.(z[:, i])) for i = 1:nrow(gt0)]

                append!(gt, gt0)
            end
        end
    end

    # delete empty rows
    if !isempty(gt)
        delete!(gt, isempty.(gt[!, :height]))
    end

    # convert from older rows of vectors format to rows of elemets [takes about 1 second]
    df0 = DataFrame()
    for r in eachrow(gt)
        n = length(r.latitude)
        df1 = DataFrame()

        df1[!, :latitude] = r.latitude
        df1[!, :longitude] = r.longitude
        df1[!, :datetime] = fill(r.datetime, n)
        df1[!, :height] = r.height
        df1[!, :quality] = r.quality
        df1[!, :height_error] = fill(r.height_error, n)
        df1[!, :id] = fill(r.granule_info[1], n)
        df1[!, :height_reference] = r.height_reference

        append!(df0, df1)
    end

    return df0
end


function geotile_build_hugonnet(geotile, geotile_dir, hstacks; force_remake=false, old_format=false)

    t1 = time()
    outfile = joinpath(geotile_dir, geotile[:id] * ".arrow")

    if !isfile(outfile) || force_remake

        gt = hstacks2geotile(geotile, hstacks; old_format)
        t2 = time()
        if !isempty(gt)
            tmp = tempname(dirname(outfile))
            Arrow.write(tmp, gt::DataFrame) # do not use compression... greatly slows read time
            mv(tmp, outfile; force=true)

            total_time = round((time()-t1)/60, digits = 2);
            save_time = round((time() - t2) / (time() - t1) * 100, digits=0)
            printstyled("    -> Hugonnet $(geotile.id) created: $(total_time) min [$(save_time)% for save] \n"; color=:light_black)
        end
    end
end

