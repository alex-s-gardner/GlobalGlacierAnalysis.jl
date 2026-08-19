"""
Synthetic multi-mission altimetry data generation for testing synthesis algorithms.

This module generates realistic synthetic elevation change data for multiple satellite missions
(ICESat-2, ICESat, GEDI, Hugonnet) with known parameters, offsets, and trends to enable
validation of multi-mission synthesis and alignment algorithms.
"""

using Dates
using DataFrames
using DimensionalData
import DimensionalData as DD

"""
    generate_hypsometric_synthetic(;
        n_dates=24,
        n_heights=10,
        trend=-0.5,
        seasonal_amplitude=0.2,
        vertical_gradient=-0.001,
        start_date=DateTime(2018,1,1),
        height_range=(1000, 3000),
        noise_sigma=0.1
    )

Generate synthetic hypsometric (elevation-binned) time series with known parameters.

# Arguments
- `n_dates`: Number of time steps
- `n_heights`: Number of elevation bins
- `trend`: Linear trend in m/yr (default: -0.5)
- `seasonal_amplitude`: Amplitude of seasonal cycle in m (default: 0.2)
- `vertical_gradient`: Elevation-dependent trend in m/yr/m (default: -0.001)
- `start_date`: Starting date (default: 2018-01-01)
- `height_range`: Tuple of (min_height, max_height) in meters (default: 1000-3000m)
- `noise_sigma`: Standard deviation of Gaussian noise (default: 0.1 m)

# Returns
- Tuple of (dates, heights, dh_matrix) where dh_matrix[date, height] contains elevation changes
"""
function generate_hypsometric_synthetic(;
    n_dates=24,
    n_heights=10,
    trend=-0.5,
    seasonal_amplitude=0.2,
    vertical_gradient=-0.001,
    start_date=DateTime(2018,1,1),
    height_range=(1000, 3000),
    noise_sigma=0.1
)
    # Generate time vector
    dates = [start_date + Month(i) for i in 0:n_dates-1]
    t = [(Dates.value(date - start_date)) / (365.25 * 24 * 60 * 60 * 1000) for date in dates]

    # Generate elevation bins
    # `range(a, b, length=1)` errors with "endpoints differ", so a single bin is special-cased to
    # the midpoint of the requested range.
    heights = n_heights == 1 ? [(height_range[1] + height_range[2]) / 2] :
              collect(range(height_range[1], height_range[2], length=n_heights))

    # Create signal matrix
    dh = zeros(n_dates, n_heights)

    for (i, ti) in enumerate(t)
        for (j, h) in enumerate(heights)
            # Elevation-dependent trend: lower elevations thin faster
            local_trend = trend + vertical_gradient * (h - mean(heights))

            # Seasonal signal
            seasonal = seasonal_amplitude * sin(2π * ti)

            # Combined signal + noise
            dh[i, j] = local_trend * ti + seasonal + noise_sigma * randn()
        end
    end

    return (dates, collect(heights), dh)
end

"""
    generate_icesat2_synthetic(;
        dates,
        heights,
        baseline_dh,
        coverage=0.95,
        uncertainty_sigma=0.15
    )

Generate synthetic ICESat-2 elevation change data with realistic coverage and uncertainties.

# Arguments
- `dates`: Vector of observation dates
- `heights`: Vector of elevation bins
- `baseline_dh`: Matrix of true elevation changes [n_dates × n_heights]
- `coverage`: Fraction of space-time bins with observations (default: 0.95)
- `uncertainty_sigma`: Observation uncertainty (default: 0.15 m)

# Returns
- Named tuple (dh, nobs, sigma) with elevation changes, observation counts, and uncertainties
"""
function generate_icesat2_synthetic(;
    dates,
    heights,
    baseline_dh,
    coverage=0.95,
    uncertainty_sigma=0.15
)
    n_dates = length(dates)
    n_heights = length(heights)

    # Initialize arrays
    dh = fill(NaN, n_dates, n_heights)
    nobs = zeros(Int, n_dates, n_heights)
    sigma = fill(NaN, n_dates, n_heights)

    # Add data with specified coverage
    for i in 1:n_dates
        for j in 1:n_heights
            if rand() < coverage
                # Add observation with noise
                dh[i, j] = baseline_dh[i, j] + uncertainty_sigma * randn()
                nobs[i, j] = rand(10:100)  # ICESat-2 has many measurements
                sigma[i, j] = uncertainty_sigma
            end
        end
    end

    return (dh=dh, nobs=nobs, sigma=sigma)
end

"""
    generate_icesat_synthetic(;
        dates,
        heights,
        baseline_dh,
        coverage=0.70,
        uncertainty_sigma=0.25,
        time_range=(DateTime(2003,1,1), DateTime(2009,12,31))
    )

Generate synthetic ICESat elevation change data with temporal constraints and sparser coverage.

# Arguments
- `dates`: Vector of observation dates
- `heights`: Vector of elevation bins
- `baseline_dh`: Matrix of true elevation changes [n_dates × n_heights]
- `coverage`: Fraction of space-time bins with observations (default: 0.70)
- `uncertainty_sigma`: Observation uncertainty (default: 0.25 m)
- `time_range`: Tuple of (start_date, end_date) for ICESat operations

# Returns
- Named tuple (dh, nobs, sigma) with elevation changes, observation counts, and uncertainties
"""
function generate_icesat_synthetic(;
    dates,
    heights,
    baseline_dh,
    coverage=0.70,
    uncertainty_sigma=0.25,
    time_range=(DateTime(2003,1,1), DateTime(2009,12,31))
)
    n_dates = length(dates)
    n_heights = length(heights)

    # Initialize arrays
    dh = fill(NaN, n_dates, n_heights)
    nobs = zeros(Int, n_dates, n_heights)
    sigma = fill(NaN, n_dates, n_heights)

    # Add data only within ICESat operational period
    for i in 1:n_dates
        if dates[i] < time_range[1] || dates[i] > time_range[2]
            continue
        end

        for j in 1:n_heights
            if rand() < coverage
                # Add observation with noise
                dh[i, j] = baseline_dh[i, j] + uncertainty_sigma * randn()
                nobs[i, j] = rand(5:30)  # ICESat has fewer measurements than ICESat-2
                sigma[i, j] = uncertainty_sigma
            end
        end
    end

    return (dh=dh, nobs=nobs, sigma=sigma)
end

"""
    generate_gedi_synthetic(;
        dates,
        heights,
        baseline_dh,
        coverage=0.80,
        uncertainty_sigma=0.30,
        time_range=(DateTime(2019,4,1), DateTime(2025,12,31)),
        latitude_range=(-52, 52)
    )

Generate synthetic GEDI elevation change data with latitude constraints.

# Arguments
- `dates`: Vector of observation dates
- `heights`: Vector of elevation bins
- `baseline_dh`: Matrix of true elevation changes [n_dates × n_heights]
- `coverage`: Fraction of space-time bins with observations (default: 0.80)
- `uncertainty_sigma`: Observation uncertainty (default: 0.30 m)
- `time_range`: Tuple of (start_date, end_date) for GEDI operations
- `latitude_range`: Tuple of (min_lat, max_lat) for GEDI coverage (±52°)

# Returns
- Named tuple (dh, nobs, sigma) with elevation changes, observation counts, and uncertainties
"""
function generate_gedi_synthetic(;
    dates,
    heights,
    baseline_dh,
    coverage=0.80,
    uncertainty_sigma=0.30,
    time_range=(DateTime(2019,4,1), DateTime(2025,12,31)),
    latitude_range=(-52, 52),
    geotile_lat=45  # Default latitude for coverage check
)
    n_dates = length(dates)
    n_heights = length(heights)

    # Check if geotile is within GEDI coverage
    if geotile_lat < latitude_range[1] || geotile_lat > latitude_range[2]
        # Return empty data outside GEDI coverage
        return (
            dh=fill(NaN, n_dates, n_heights),
            nobs=zeros(Int, n_dates, n_heights),
            sigma=fill(NaN, n_dates, n_heights)
        )
    end

    # Initialize arrays
    dh = fill(NaN, n_dates, n_heights)
    nobs = zeros(Int, n_dates, n_heights)
    sigma = fill(NaN, n_dates, n_heights)

    # Add data only within GEDI operational period
    for i in 1:n_dates
        if dates[i] < time_range[1] || dates[i] > time_range[2]
            continue
        end

        for j in 1:n_heights
            if rand() < coverage
                # Add observation with noise
                dh[i, j] = baseline_dh[i, j] + uncertainty_sigma * randn()
                nobs[i, j] = rand(20:80)  # GEDI intermediate sample size
                sigma[i, j] = uncertainty_sigma
            end
        end
    end

    return (dh=dh, nobs=nobs, sigma=sigma)
end

"""
    generate_hugonnet_synthetic(;
        dates,
        heights,
        baseline_dh,
        coverage=0.85,
        uncertainty_sigma=0.20,
        time_range=(DateTime(2000,1,1), DateTime(2020,12,31)),
        bias_offset=0.0
    )

Generate synthetic Hugonnet elevation change data with potential systematic bias.

# Arguments
- `dates`: Vector of observation dates
- `heights`: Vector of elevation bins
- `baseline_dh`: Matrix of true elevation changes [n_dates × n_heights]
- `coverage`: Fraction of space-time bins with observations (default: 0.85)
- `uncertainty_sigma`: Observation uncertainty (default: 0.20 m)
- `time_range`: Tuple of (start_date, end_date) for Hugonnet data
- `bias_offset`: Systematic vertical bias in meters (default: 0.0)

# Returns
- Named tuple (dh, nobs, sigma) with elevation changes, observation counts, and uncertainties
"""
function generate_hugonnet_synthetic(;
    dates,
    heights,
    baseline_dh,
    coverage=0.85,
    uncertainty_sigma=0.20,
    time_range=(DateTime(2000,1,1), DateTime(2020,12,31)),
    bias_offset=0.0
)
    n_dates = length(dates)
    n_heights = length(heights)

    # Initialize arrays
    dh = fill(NaN, n_dates, n_heights)
    nobs = zeros(Int, n_dates, n_heights)
    sigma = fill(NaN, n_dates, n_heights)

    # Add data only within Hugonnet time period
    for i in 1:n_dates
        if dates[i] < time_range[1] || dates[i] > time_range[2]
            continue
        end

        for j in 1:n_heights
            if rand() < coverage
                # Add observation with noise + systematic bias
                dh[i, j] = baseline_dh[i, j] + bias_offset + uncertainty_sigma * randn()
                nobs[i, j] = rand(50:200)  # Hugonnet stacks have many DEM pairs
                sigma[i, j] = uncertainty_sigma
            end
        end
    end

    return (dh=dh, nobs=nobs, sigma=sigma)
end

"""
    generate_multimission_synthetic(;
        n_dates=60,
        n_heights=10,
        trend=-0.5,
        seasonal_amplitude=0.2,
        start_date=DateTime(2000,1,1),
        height_range=(1000, 3000),
        include_missions=[:icesat2, :icesat, :gedi, :hugonnet],
        mission_biases=Dict(:hugonnet => 0.5),  # Example: Hugonnet has 0.5m positive bias
        geotile_lat=45
    )

Generate complete multi-mission synthetic dataset with realistic characteristics.

# Arguments
- `n_dates`: Number of monthly time steps (default: 60 = 5 years)
- `n_heights`: Number of elevation bins (default: 10)
- `trend`: True linear trend in m/yr (default: -0.5)
- `seasonal_amplitude`: Seasonal amplitude in m (default: 0.2)
- `start_date`: Starting date (default: 2000-01-01 for full mission coverage)
- `height_range`: Elevation range in meters (default: 1000-3000m)
- `include_missions`: Missions to include (default: all four)
- `mission_biases`: Dictionary of systematic biases by mission (default: none)
- `geotile_lat`: Geotile latitude for GEDI coverage (default: 45°N)

# Returns
- Dictionary with keys for each mission, containing (dh, nobs, sigma) named tuples
- Also returns :dates, :heights, and :truth (true elevation changes) keys
"""
function generate_multimission_synthetic(;
    n_dates=60,
    n_heights=10,
    trend=-0.5,
    seasonal_amplitude=0.2,
    start_date=DateTime(2000,1,1),
    height_range=(1000, 3000),
    include_missions=[:icesat2, :icesat, :gedi, :hugonnet],
    mission_biases=Dict{Symbol, Float64}(),
    geotile_lat=45
)
    # Generate ground truth
    dates, heights, baseline_dh = generate_hypsometric_synthetic(
        n_dates=n_dates,
        n_heights=n_heights,
        trend=trend,
        seasonal_amplitude=seasonal_amplitude,
        start_date=start_date,
        height_range=height_range,
        noise_sigma=0.0  # No noise in baseline
    )

    # Initialize output dictionary
    data = Dict{Symbol, Any}(
        :dates => dates,
        :heights => heights,
        :truth => baseline_dh
    )

    # Generate mission-specific data
    if :icesat2 in include_missions
        bias = get(mission_biases, :icesat2, 0.0)
        data[:icesat2] = generate_icesat2_synthetic(
            dates=dates,
            heights=heights,
            baseline_dh=baseline_dh .+ bias
        )
    end

    if :icesat in include_missions
        bias = get(mission_biases, :icesat, 0.0)
        data[:icesat] = generate_icesat_synthetic(
            dates=dates,
            heights=heights,
            baseline_dh=baseline_dh .+ bias
        )
    end

    if :gedi in include_missions
        bias = get(mission_biases, :gedi, 0.0)
        data[:gedi] = generate_gedi_synthetic(
            dates=dates,
            heights=heights,
            baseline_dh=baseline_dh .+ bias,
            geotile_lat=geotile_lat
        )
    end

    if :hugonnet in include_missions
        bias = get(mission_biases, :hugonnet, 0.0)
        data[:hugonnet] = generate_hugonnet_synthetic(
            dates=dates,
            heights=heights,
            baseline_dh=baseline_dh .+ bias
        )
    end

    return data
end

"""
    synthetic_fill_params(dh; missions2align2=String[], n_model_params=9)

Build the per-mission `params_fill` table that the `hyps_*` gap-filling routines expect.

`dh` is the mission-keyed dictionary of `(geotile, date, height)` DimArrays. For each mission this
returns a DataFrame with one row per geotile, matching the structure built in
`geotile_binning`/`utilities_binning.jl`: bookkeeping counters, a per-geotile model-parameter
vector, reference offsets, and one offset column per alignment-reference mission.

# Why this exists

`hyps_model_fill!`, `hyps_align_dh!` and friends take this table as their **positional** `params`
argument -- it is where fitted coefficients are written back. It is *not* the bag of tuning knobs
(`bincount_min`, `smooth_n`, ...), which are keyword arguments. Passing a NamedTuple of knobs as
`params` fails with `MethodError: no method matching getindex(::@NamedTuple{...}, ::String)`.
"""
function synthetic_fill_params(dh; missions2align2=String[], n_model_params=9)
    params = Dict()
    for mission in keys(dh)
        geotiles = collect(dims(dh[mission], :geotile))
        n = length(geotiles)

        params[mission] = DataFrame(
            geotile=geotiles,
            nobs_raw=zeros(n), nbins_raw=zeros(n),
            nobs_final=zeros(n), nbins_filt1=zeros(n),
            param_m1=[fill(NaN, n_model_params) for _ in 1:n],
            h0=fill(NaN, n),
            t0=fill(NaN, n),
            dh0=fill(NaN, n),
            bin_std=fill(NaN, n),
            bin_anom_std=fill(NaN, n),
        )

        for mission_ref in missions2align2
            params[mission][!, "offset"] = zeros(n)
            params[mission][!, "offset_$mission_ref"] = fill(NaN, n)
            params[mission][!, "offset_nmad_$mission_ref"] = fill(NaN, n)
            params[mission][!, "offset_nobs_$mission_ref"] = zeros(Int64, n)
        end
    end
    return params
end
