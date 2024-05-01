using Altim
using GeoTiles
using Plots
using BinStatistics
using DataFrames
using Statistics
using LsqFit

project_id = :v01;
paths = project_paths(project_id = project_id);
products = project_products(; project_id=project_id);
dem = :cop30_v2; #:cop30_v2, :nasadem_v1, :rema_v2_10m, :arcticdem_v3_10m

products = (icesat=products.icesat,)

(extent, epsg) = Altim.region_extent(:RGI02);

extent = Extent(X=(-141.0, -141.5), Y=(58.0, 60.0))
extent = Extent(X=(-119.0, -118.5), Y=(39, 39.5))
#extent = Extent(X=(-131.8, -131.0), Y=(58.1, 58.2))
# define a dictionary mapping of user defined variable table name [used in geotile_binstats] 
# to file suffix, use Symbols

var_suffix = Dict(:dem => dem, :altim => :arrow, :masks => :masks, :canopyh => :canopyh);

# This function accepts a named tuple of geotile variables as defined by `var_suffix` and 
# retuns a BinStatistics DataFrame. 
function geotile_binstats_curv(gts)

    # convert data from rows of rows to vectors
    ice = reduce(vcat, gts.masks.landice)
    lat = reduce(vcat, gts.altim.latitude)
    height0 = reduce(vcat, gts.altim.height)

    dhddx = reduce(vcat, gts.dem.dhddx)
    dhddy = reduce(vcat, gts.dem.dhddy)
    height = reduce(vcat, gts.dem.height)

    dh = height0 .- height
    
    canopyh = reduce(vcat, gts.canopyh.canopyh)

    foo = Altim.meters2latlon_distance.(Ref(1), lat) # degree per meter
    lat_dist = getindex.(foo, 1)
    lon_dist = getindex.(foo, 2)

    curvature = -2(dhddx .* lon_dist .^ 2 + dhddy .* lat_dist .^ 2) * 100

    df = DataFrame(; dh, ice, curvature, height)
   
    curvature_bins = -1:0.05:1
    valid = (canopyh .<= 1) .& (abs.(dh) .< 100)
    idx = .!isnan.(dh) .& .!isnan.(curvature) .& .!ice .& valid
    if any(idx)
        bs1 = binstats(df[idx, :], :curvature, curvature_bins, :dh; col_function=[median, Altim.mad])
    else
        bs1 = DataFrame()
    end

    height_bins = 0:100:9000
    idx = .!isnan.(dh) .& ice
    if any(idx)
        bs2 = binstats(df[idx, :], :height, height_bins, :dh; col_function=[median, Altim.mad])
    else
        bs2 = DataFrame()
    end

    return bs1
end


function geotile_binstats_dhdt(gts)

    # convert data from rows of rows to vectors
    ice = reduce(vcat, gts.masks.landice)
    decyear = decimalyear.(reduce(vcat, gts.altim.datetime))
    height0 = reduce(vcat, gts.altim.height)
    height = reduce(vcat, gts.dem.height)

    dh = height0 .- height
    dt = decyear .- 2013

    dhdt = dh./dt

    df = DataFrame(; dhdt, ice, height)
   
    height_bins = 0:100:9000
    idx = .!isnan.(dhdt) .& ice .& (abs.(dhdt) .< 10)


    if any(idx)
        bs1 = binstats(df[idx, :], :height, height_bins, :dhdt; col_function=[median, Altim.mad])
    else
        bs1 = DataFrame()
    end

    return bs1
end


function binstats_combine(bsvec; fun=mean, weighted=false)
    # concatinate vector of data frames into a single dataframe
    ind = .!isempty.(bsvec);
    bsvec = reduce(vcat, bsvec[ind])

    # remove id column
    select!(bsvec, Not(:id))
    colname = names(bsvec)

    if weighted
        for cn in colname[3:end]
            bsvec[:, cn] = bsvec[:, cn] .* bsvec[:, :nrow]
        end
    end

    bsc = combine(groupby(bsvec, colname[1]),
        :nrow => sum,
        colname[3:end] .=> fun,
        nrow => :cnt;
        renamecols=false)

    if weighted
        for cn in colname[3:end]
            bsc[:, cn] = bsc[:, cn] ./ (bsc[:, :nrow] ./ bsc[:, :cnt])
        end
    end

    return bsc
end

function geotile_stats(path2dir, var_suffix, geotile_binstat; extent = nothing)
    # find requested geotiles within region, only returning geotiles where all suffixes exist
    gtfilelist = GeoTiles.listtiles_intersecting(path2dir, values(var_suffix); extent)
    
    bsvec = Vector{Vector{Any}}()
    for i in 1:Threads.nthreads()
        push!(bsvec, Any[])
    end

    # loop for each geotile
    Threads.@threads for gtinfo in eachrow(gtfilelist)
    #for gtinfo in eachrow(gtfilelist)    
        # read data
        gts = NamedTuple{Tuple(keys(var_suffix))}([GeoTiles.read(gtinfo[suffix]) for suffix in values(var_suffix)])
     
        # calculate binstatistics 
        foo = geotile_binstat(gts)
        
        # add tile ID
        foo[!, :id] .= gtinfo.id

        push!(bsvec[Threads.threadid()], foo)
    end;

    bsvec = reduce(vcat, bsvec);

    return bsvec
end

bsvec = Dict{String,Vector{Any}}();
#for mission in (:icesat, :icesat2, :gedi, :hugonnet)
mission = :hugonnet

    #bsvec = geotile_stats(path2dir, var_suffix, geotile_binstats_dhdt; extent)
    bsvec[String(mission)] = geotile_stats(paths[mission].geotile, var_suffix, geotile_binstats_curv; extent)
#end;


# create model to fit to data
@. model(x, p) = p[1] + x * p[2]
p0 = [0.,0.]


for key in keys(bsvec)
    if !isempty(bsvec[key])
        p1 = plot(legend=false, title="$key minus $dem")
        for bs in bsvec[key]
            x = BinStatistics.bincenter.(bs[:,1])
            y = bs[:,3]
            valid = bs[:,:nrow] .> 100
            plot!(x[valid],y[valid])
        end

        xlabel!("curvature")
        ylabel!("height anomaly")
        bs = binstats_combine(copy(bsvec[key]); fun=mean, weighted=true)
        x = BinStatistics.bincenter.(bs[:, 1])
        y = bs[:, 3]
        valid = bs[:, :nrow] .> 100

        plot!(x[valid], y[valid], linewidth=3, linecolor = :black)
        xlims!(-1, 1)

        p2 = bar(legend=false, x, bs.nrow / sum(bs.nrow))
        xlims!(-1, 1)
        xlabel!("curvature")
        ylabel!("density")
        p = plot(p1, p2, layout=(2, 1), size=(600, 800))
    end
end

p = plot();
for key in keys(bsvec)
    if !isempty(bsvec[key])
        bs = binstats_combine(bsvec[key]; fun=mean, weighted=true)
        bs[!, :bincenter] = BinStatistics.bincenter.(bs[:, 1])
        sort!(bs, :bincenter)
        x = bs[:, :bincenter]
        y = bs[:, 3]
        valid = bs[:, :nrow] .> 100
        plot!(x[valid], y[valid],label=key)
    end
end

p = scatter();
for key in keys(bsvec)
    p = plot(legend = false);
    for bs in bsvec[key]
        x = BinStatistics.bincenter.(bs[:,1])
        y = bs[:,3]
        valid = bs[:,:nrow] .> 100

        fit = curve_fit(model, Float32.(x[valid]), Float32.(y[valid]), p0)
        lat = mean(GeoTiles.extent(bs[1,:id]).Y)
        scatter!([lat],[fit.param[2]])
    end
end


p