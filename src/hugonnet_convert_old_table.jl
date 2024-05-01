using DataFrames;
using Altim;
using Arrow;

folder = "/mnt/bylot-r3/data/hugonnet/HSTACK/001/geotile/2deg/";
fns = allfiles(folder; fn_endswith=".arrow");

for outfile in fns
    reformat(outfile)
end


function reformat(outfile)
    t1 = time()
    df = DataFrame(Arrow.Table(outfile))

    T = SubArray{Float32, 1, Arrow.Primitive{Float32, Vector{Float32}}, Tuple{UnitRange{Int64}}, true}

    if any(typeof.((df[:,:latitude])) .== T)
        gt = DataFrame()
        for r in eachrow(df)
            n = length(r.latitude)
            df1 = DataFrame()

            df1[!, :latitude] = r.latitude
            df1[!, :longitude] = r.longitude
            df1[!, :datetime] = fill(r.datetime,n)
            df1[!, :height] = r.height
            df1[!, :quality] = r.quality
            df1[!, :height_error] = fill(r.height_error, n)
            df1[!, :id] = fill(r.granule_info[1], n)
            df1[!, :height_reference] = r.height_reference

            append!(gt,df1)
        end

        if !isempty(gt)
            tmp = tempname(dirname(outfile))
            Arrow.write(tmp, gt::DataFrame) # do not use compression... greatly slows read time
            mv(tmp, outfile; force=true)

            total_time = round((time() - t1))
            printstyled("    -> $(outfile) re-formatted: $(total_time) s] \n"; color=:light_black)
        end
    end
end
