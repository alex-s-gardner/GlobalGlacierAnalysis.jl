using Arrow
using DataFrames
using BinStatistics
using Statistics
using Dates

#path to file
fn = "/mnt/bylot-r3/data/icesat2/ATL06/006/geotile/2deg/lat[-74-72]lon[-064-062].arrow"
fn = "/mnt/bylot-r3/data/hugonnet/HSTACK/001/geotile/2deg/lat[+56+58]lon[-134-132].arrow"

# define some functions
"""
    decimalyear(datetime::DateTime)

converts DateTime to decimalyear [e.g. decimalyear(DateTime("2018-10-24")) = 2018.813698630137]
"""
function decimalyear(datetime)

    # check if date is a leapyear
    if isleapyear(datetime)
        number_of_days = 366
    else
        number_of_days = 365
    end

    # calculate decimal year
    decyear = year(datetime) + dayofyear(datetime) / number_of_days

    return decyear
end

"""
    madnorm_filter(x, threshold)

Filter data points based on their normalized median absolute deviation (MAD).

Takes a vector `x` and returns a boolean mask indicating which values fall within
`threshold` standard deviations of the median, using MAD as a robust estimate of spread.

The MAD is normalized by 1.4826 to make it consistent with standard deviation for normal distributions.

# Arguments
- `x`: Vector of values to filter
- `threshold`: Number of standard deviations to use as cutoff

# Returns
Boolean vector same length as `x`, with `true` for values within `threshold` standard deviations
"""
function madnorm(x)
    consistent_estimator = 1.4826 #mad to sigma
    x_abs = abs.(x .- median(x))
    return x_madnorm = x_abs ./ (median(x_abs) .* consistent_estimator)
end


# read the file
arrow_table = DataFrame(Arrow.Table(fn))

# convert to dataframe
df = DataFrame(arrow_table)

# calculate height difference (dh)
df[!, :dh] = df.height .- df.height_reference

#convert datetime to decimal year
df[!, :decyear] = decimalyear.(df.datetime)

# identify valid data
valid = abs.(df.height) .< 7000

# define bins
date_bins = 2000.:30/365:2025
height_bins = 0.:100.:10000.

# time = 42s
@time df0 = binstats(df[valid, :], [:decyear, :height_reference],
    [date_bins, height_bins], :dh; col_function=[madnorm], missing_bins=true)