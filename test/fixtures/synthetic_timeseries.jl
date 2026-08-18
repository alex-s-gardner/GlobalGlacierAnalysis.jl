"""
Synthetic time series generation for testing temporal analysis
"""

using Dates

"""
    generate_trend_seasonal_ts(; n_points=50, trend=0.5, amplitude=2.0, noise_sigma=0.1, start_year=2018)

Generate time series with linear trend + seasonal component + noise.

# Arguments
- `n_points`: Number of time points
- `trend`: Linear trend (units per year)
- `amplitude`: Seasonal amplitude
- `noise_sigma`: Standard deviation of Gaussian noise
- `start_year`: Starting year

# Returns
- Tuple of (dates, values)
"""
function generate_trend_seasonal_ts(;
    n_points=50,
    trend=0.5,
    amplitude=2.0,
    noise_sigma=0.1,
    start_year=2018
)
    # Generate quarterly dates
    dates = [DateTime(start_year, 1, 1) + Month(3*i) for i in 0:n_points-1]

    # Convert to decimal years
    t = [(Dates.value(date - DateTime(start_year,1,1))) / (365.25 * 24 * 60 * 60 * 1000) for date in dates]

    # Generate signal
    signal = trend .* t .+ amplitude .* sin.(2π .* t)

    # Add noise
    noise = noise_sigma .* randn(n_points)

    values = signal .+ noise

    return (dates, values)
end

"""
    generate_ts_with_gaps(; n_points=20, gap_fraction=0.2)

Generate time series with intentional gaps for interpolation testing.

# Returns
- Tuple of (dates, values) where some values are `missing`
"""
function generate_ts_with_gaps(;
    n_points=20,
    gap_fraction=0.2,
    trend=1.0,
    start_year=2018
)
    dates = [DateTime(start_year, 1, 1) + Month(i) for i in 0:n_points-1]
    t = [i / 12.0 for i in 0:n_points-1]  # Years

    # Generate underlying signal
    signal = trend .* t

    # Create gaps
    n_gaps = floor(Int, n_points * gap_fraction)
    gap_indices = rand(2:n_points-1, n_gaps)  # Don't gap first/last points

    values = Vector{Union{Float64, Missing}}(signal)
    values[gap_indices] .= missing

    return (dates, values)
end

"""
    generate_quadratic_seasonal_ts(; n_points=50, params=[5.0, 0.5, 0.1, 2.0])

Generate time series with quadratic trend + seasonal component.

# Arguments
- `params`: [offset, linear_coef, quadratic_coef, amplitude]

# Returns
- Tuple of (times, values) where times are in decimal years
"""
function generate_quadratic_seasonal_ts(;
    n_points=50,
    params=[5.0, 0.5, 0.1, 2.0],
    noise_sigma=0.05
)
    t = range(0, 2, length=n_points)  # 0 to 2 years

    # y = p[1] + p[2]*t + p[3]*t² + p[4]*sin(2πt)
    signal = params[1] .+ params[2] .* t .+ params[3] .* t.^2 .+ params[4] .* sin.(2π .* t)

    # Add noise
    values = signal .+ noise_sigma .* randn(n_points)

    return (collect(t), values)
end
