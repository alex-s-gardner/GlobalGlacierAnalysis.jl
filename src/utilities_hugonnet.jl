"""
    hstack_catalogue(hstack_parent_dir; force_remake = false)

Create or load a geospatial catalogue of Hugonnet elevation change data stacks.

# Arguments
- `hstack_parent_dir`: Directory containing Hugonnet data stack files (.nc)
- `force_remake`: Whether to force rebuilding the catalogue even if it exists (default: false)

# Returns
- DataFrame containing metadata for each stack file, including spatial extents, 
  projection information, time ranges, and file paths

# Description
This function scans the provided directory for Hugonnet elevation change data files,
extracts their metadata (projection, spatial extent, time range, etc.), and compiles
it into a searchable catalogue. The catalogue is saved as an Arrow file for faster
future access. For version 1 data, it applies a half-pixel offset correction.
"""
function hstack_catalogue(hstack_parent_dir; force_remake = false)
    outfile = joinpath(hstack_parent_dir,"hstack_catalogue.arrow")

    if !isfile(outfile) || force_remake

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
        
            ##----------------------------------------
            # edgest => centers [this was an issue with v1 hstacks]
            if !Base.contains(hstack_parent_dir, "HSTACK/001")
                error("hugonnet verson 1 data had a 1/2 pixel offset that is corrected for here, check different versions before useing")
            end

            NCDataset(file) do ds

                hstacks[i, :ProjString] = convert(GFT.ProjString, GFT.WellKnownText(ds["crs"].attrib["spatial_ref"]))

                x = ds["x"]
                hstacks[i, :x] = (x[1]:(x[2]-x[1]):(x[1]+(x[2]-x[1])*(length(x)-1))) .+ 50

                y = ds["y"]
                hstacks[i, :y] = (y[1]:(y[2]-y[1]):(y[1]+(y[2]-y[1])*(length(y)-1))) .- 50
                ##----------------------------------------

                hstacks[i, :time_range] = extrema(ds["time"])
                hstacks[i, :nlayers] = length(ds["time"])
            end

            # neet to fix this once I figure out how.
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

"""
    hstacks2geotile(geotile, hstacks)

Extract Hugonnet elevation change data for a specific geotile from a collection of data stacks.

# Arguments
- `geotile`: The geotile object containing spatial extent information
- `hstacks`: DataFrame of Hugonnet data stacks with metadata

# Returns
- DataFrame containing extracted elevation data points within the geotile, with columns for
  latitude, longitude, datetime, height, quality, height_error, id, and height_reference

# Description
This function identifies which Hugonnet data stacks intersect with the provided geotile,
extracts the relevant data points, and transforms them into a standardized DataFrame format.
The function handles coordinate transformations, filters points to those within the geotile
boundary, and restructures the data from a row-of-vectors format to a standard tabular format.
"""
function hstacks2geotile(geotile, hstacks)
    # find intersecting tiles
    stacksind = findall(Extents.intersects.(hstacks.GeoExtent, Ref(geotile.extent)))

    if isempty(stacksind)
        return nothing
    end

    gt = DataFrame()

    for sind in stacksind

        # extract data
        x = collect(hstacks[sind, :x]) * ones(1, length(hstacks[sind, :y]))
        y = ones(length(hstacks[sind, :x])) * collect(hstacks[sind, :y])'
        
        proj2geo = Proj.Transformation(hstacks[sind, :ProjString].val, "EPSG:4326"; ctx=Proj.proj_context_clone())
        lat_lon = proj2geo.(x, y)

        isin = GlobalGlacierAnalysis.within.(Ref(geotile.extent), getindex.(lat_lon, 2), getindex.(lat_lon, 1))

        if !any(isin)
            return nothing
        else
            rrange, crange = GlobalGlacierAnalysis.validrange(isin)
            isin = isin[rrange, crange]
            
            NCDataset(hstacks[sind, :path]) do ds
              
                ll = lat_lon[rrange, crange]
               
                t = ds["time"][:]::Vector{DateTime}
                if any(t .< Date(2000, 1, 1)) || any(t .> Date(2030, 1, 1))
                    error("data is being read incorrectly")
                end
 
                z = ds["z"][rrange, crange, :]::Array{Union{Missing,Float32},3}
                (m, n, p) = size(z)

                isin = repeat(isin, 1, 1, length(t))
                ll = repeat(ll, 1, 1, p)

                t = repeat(reshape(t, 1, 1, p), m, n, 1)

                ref_z = ds["ref_z"][rrange, crange]::Matrix{Union{Missing,Float32}}
                ref_z = repeat(ref_z, 1, 1, p)

                uncert = ds["uncert"][:]::Vector{Float32}
                uncert = repeat(reshape(uncert, 1, 1, p), m, n, 1)

                dem_names = ds["dem_names"][:]::Vector{String}
                dem_names = repeat(reshape(dem_names, 1, 1, p), m, n, 1)

                corr = ds["corr"][rrange, crange, :]::Array{Union{Missing,Int8},3}

                valid = .!ismissing.(z) .& isin .& .!isnan.(z) .& (corr .> 70)

                gt0 = DataFrame()
                gt0[!, :datetime] = t[valid]
                gt0[!, :height] = Float32.(z[valid])
                gt0[!, :latitude] = Float32.(getindex.(ll, 1)[valid])
                gt0[!, :longitude] = Float32.(getindex.(ll, 2)[valid])

                gt0[!, :granule_info] = dem_names[valid]
                gt0[!, :quality] = corr[valid] .> 70
                gt0[!, :height_error] .= uncert[valid]
                gt0[!, :height_reference] = Float32.(ref_z[valid])

                append!(gt, gt0)
            end
        end
    end
    return gt
end


"""
    geotile_build_hugonnet(geotile, geotile_dir, hstacks; force_remake=false)

Process and save Hugonnet elevation data for a specific geotile.

# Arguments
- `geotile`: Geotile information containing ID and boundaries
- `geotile_dir`: Directory where the processed geotile data will be saved
- `hstacks`: Hugonnet data stacks to process
- `force_remake`: Whether to recreate the file even if it already exists (default: false)

# Description
Converts Hugonnet elevation data stacks to geotile format and saves as an Arrow file.
The function tracks processing time and reports completion statistics.
"""
function geotile_build_hugonnet(geotile, geotile_dir, hstacks; force_remake=false)

    t1 = time()
    outfile = joinpath(geotile_dir, geotile[:id] * ".arrow")
    
    if !isfile(outfile) || force_remake
        gt = hstacks2geotile(geotile, hstacks)
        t2 = time()
        if !isnothing(gt) && !isempty(gt)
            tmp = tempname(dirname(outfile))
            Arrow.write(tmp, gt::DataFrame) # do not use compression... greatly slows read time
            mv(tmp, outfile; force=true)

            total_time = round((time()-t1)/60, digits = 2);
            save_time = round((time() - t2) / (time() - t1) * 100, digits=0)
            printstyled("    -> Hugonnet $(geotile.id) created: $(total_time) min [$(save_time)% for save] \n"; color=:light_black)
        else
            printstyled("    -> Hugonnet data does not overlap with $(geotile.id), skipping \n"; color=:light_red)
        end
    else
        printstyled("\n    -> Hugonnet $(geotile.id) already exists \n"; color=:light_green)
    end
end