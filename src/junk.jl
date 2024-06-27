files = allfiles("/mnt/bylot-r3/data/hugonnet/HSTACK/001/geotile/2deg/"; fn_endswith=".arrow")

using Arrow
using DataFrames
using Dates

for f in files
    println(f)

    df = DataFrame(Arrow.Table(f))

    df.datetime = df.datetime .- (Date(2000, 1, 1) - Date(1900, 1, 1))

    tmp = tempname(dirname(f))
    Arrow.write(tmp, df::DataFrame) # do not use compression... greatly slows read time
    mv(tmp, f; force=true)
end