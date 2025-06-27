using FileIO, DataFrames, Arrow
import GlobalGlacierAnalysis as GA

base_fn = "/mnt/bylot-r3/data/binned/2deg/geotiles_2deg.arrow"
df0 = DataFrame(Arrow.Table(base_fn))
base_names = names(df0)


fns = GA.allfiles("/mnt/bylot-r3/data/binned/2deg"; fn_startswith="geotile")

#for fn in fns
fn = fns[3]

#if fn == base_fn
#    continue
#else
    df2 = DataFrame(Arrow.Table(fn))
    @assert df0.id == df2.id

    for name0 in base_names
        df2[!, name0] = df0[:, name0]
    end


    geom = GGA.extent2rectangle.(getindex.(df.extent, :bounds))

    valid = sum.(df2.glacier_area_km2) .> 0

    plot(geom[valid], color=:red, alpha=0.5)

    #Arrow.write(fn, df2)

end

