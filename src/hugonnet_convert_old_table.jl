"""
    hugonnet_convert_old_table.jl

Convert Hugonnet elevation data from old table format to a standardized structure.

This script:
1. Processes Arrow files containing elevation data in the old format
2. Identifies files that need reformatting based on their data structure
3. Converts nested array structures into flat DataFrame records
4. Preserves all original data fields while standardizing the format
5. Writes the reformatted data back to the original files

The reformatting ensures consistent data access patterns for downstream analysis.
"""

using DataFrames;
using GlobalGlacierAnalysis;
using Arrow;

folder = "/mnt/bylot-r3/data/hugonnet/HSTACK/001/geotile/2deg/";
fns = allfiles(folder; fn_endswith=".arrow");

for outfile in fns
    reformat(outfile)
end

"""
    reformat(outfile)

Reformat a Hugonnet elevation data file from old nested format to flat structure.

This function:
1. Loads an Arrow file into a DataFrame
2. Detects if the file uses the old format with nested arrays
3. Converts each row with array fields into multiple flat records
4. Preserves all metadata and measurement values
5. Writes the reformatted data back to the original file

# Arguments
- `outfile`: Path to the Arrow file to be reformatted
"""
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
