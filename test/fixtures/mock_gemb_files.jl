"""
Mock GEMB (Glacier Energy and Mass Balance) file generation for testing.

GEMB model outputs are stored in MATLAB .mat files with specific variable names
and structure. This module generates minimal mock files for testing I/O and calibration.
"""

using Dates
using MAT

"""
    create_mock_gemb_mat(filepath; n_points=10, n_times=36, pscale=1.0)

Create a mock GEMB output file in MATLAB .mat format.

# Arguments
- `filepath`: Output path for the .mat file
- `n_points`: Number of spatial points (default: 10)
- `n_times`: Number of time steps (default: 36 = 3 years monthly)
- `pscale`: Precipitation scaling parameter (default: 1.0)

# Returns
- Nothing (writes file to disk)

# GEMB File Structure
GEMB .mat files contain these variables:
- lat, lon: Coordinates [n_points]
- time: Decimal years [n_times]
- SMB, EC, Accumulation, Runoff, Melt, Refreeze, Rain: **per-interval increments in mm**
  [n_points × n_times]
- FAC: Firn air content, in metres and not accumulated [n_points × n_times]
- FACtoDepth: FAC to depth conversion factor [scalar]
- H: Surface height [n_points × n_times]
- Ta: Air temperature [n_points × n_times]

!!! note "Flux variables are increments, not running totals"
    `gemb_read2` applies `cumsum(x ./ 1000, dims=2)` to SMB, EC, Accumulation, Runoff, Melt,
    Refreeze and Rain, i.e. it expects each entry to be that interval's increment in **mm** and
    produces cumulative **metres**. This fixture previously wrote already-accumulated metres, so
    reading it integrated twice and produced non-monotonic, sign-scrambled series.

    The flux increments here are also mutually consistent -- `Runoff = Melt - Refreeze` and
    `SMB = Accumulation - Runoff - EC` -- and only precipitation-driven terms scale with `pscale`,
    so raising `pscale` raises SMB (makes it less negative) as it should physically.
"""
function create_mock_gemb_mat(filepath; n_points=10, n_times=36, pscale=1.0)
    # Create spatial coordinates (scattered over a small region)
    lats = 60.0 .+ 2.0 .* rand(n_points)
    lons = -50.0 .+ 2.0 .* rand(n_points)

    # Create time vector (decimal years)
    start_year = 2015.0
    times = start_year .+ collect(0:n_times-1) ./ 12.0  # Monthly

    # Initialize data arrays [n_points × n_times]
    data = Dict{String, Any}()
    data["lat"] = lats
    data["lon"] = lons
    data["time"] = times

    # Flux variables are per-interval increments in mm; gemb_read2 turns them into cumulative
    # metres. One monthly step of a rate r m/yr is r/12 m = r/12*1000 mm.
    mm_per_step(rate_m_per_yr) = rate_m_per_yr / 12.0 * 1000.0

    # Only precipitation-driven terms scale with pscale; melt and refreeze do not.
    acc_step = mm_per_step(2.0 * pscale)        # snow input
    rain_step = mm_per_step(0.2 * pscale)
    melt_step = mm_per_step(2.8)                # surface melt
    refreeze_step = mm_per_step(0.3)            # refreezing of meltwater
    ec_step = mm_per_step(0.1)                  # sublimation / condensation

    # Small positive jitter, kept well below the mean so increments stay non-negative and the
    # cumulative series stay monotonic.
    jitter(scale) = scale .* abs.(randn(n_points, n_times))

    acc_inc = fill(acc_step, n_points, n_times) .+ jitter(0.02 * acc_step)
    rain_inc = fill(rain_step, n_points, n_times) .+ jitter(0.02 * rain_step)
    melt_inc = fill(melt_step, n_points, n_times) .+ jitter(0.02 * melt_step)

    # refreeze must never exceed melt
    refreeze_inc = min.(fill(refreeze_step, n_points, n_times) .+ jitter(0.02 * refreeze_step),
                        melt_inc)

    ec_inc = fill(ec_step, n_points, n_times) .+ jitter(0.02 * ec_step)

    # Keep the identities exact so tests can assert them:
    #   runoff = melt - refreeze          (rain excluded by design)
    #   smb    = accumulation - runoff - ec
    runoff_inc = melt_inc .- refreeze_inc
    smb_inc = acc_inc .- runoff_inc .- ec_inc

    data["Accumulation"] = acc_inc
    data["Rain"] = rain_inc
    data["Melt"] = melt_inc
    data["Refreeze"] = refreeze_inc
    data["Runoff"] = runoff_inc
    data["EC"] = ec_inc
    data["SMB"] = smb_inc

    # cumulative metres, used below for surface height
    ec_cumulative = cumsum(ec_inc ./ 1000, dims=2)

    # FAC: Firn air content (non-cumulative, slowly evolving)
    fac = zeros(n_points, n_times)
    for i in 1:n_points
        base_fac = 2.0 + 0.5 * randn()  # Base FAC ~2m
        for j in 1:n_times
            fac[i, j] = base_fac + 0.1 * sin(2π * times[j]) + 0.05 * randn()
        end
    end
    data["FAC"] = fac

    # FACtoDepth: Conversion factor (scalar)
    data["FACtoDepth"] = 0.85

    # H: Surface height (non-cumulative, changes with EC)
    h_initial = 1500.0  # meters
    height = zeros(n_points, n_times)
    for i in 1:n_points
        h0 = h_initial + 100 * randn()
        for j in 1:n_times
            height[i, j] = h0 + ec_cumulative[i, j]
        end
    end
    data["H"] = height

    # Ta: Air temperature (seasonal, elevation-dependent)
    temp = zeros(n_points, n_times)
    for i in 1:n_points
        base_temp = -10.0 + 5.0 * randn()  # Base temperature
        for j in 1:n_times
            seasonal = 15.0 * sin(2π * (times[j] - 0.25))  # Peak in summer
            temp[i, j] = base_temp + seasonal + 2.0 * randn()
        end
    end
    data["Ta"] = temp

    # Write to .mat file
    matwrite(filepath, data)

    return nothing
end

"""
    create_mock_gemb_ensemble(output_dir; n_members=5, pscale_range=(0.8, 1.2))

Create an ensemble of mock GEMB files with varying precipitation scaling.

# Arguments
- `output_dir`: Directory to write ensemble files
- `n_members`: Number of ensemble members (default: 5)
- `pscale_range`: Range of precipitation scaling parameters (default: 0.8 to 1.2)

# Returns
- Vector of file paths to created ensemble members
"""
function create_mock_gemb_ensemble(output_dir; n_members=5, pscale_range=(0.8, 1.2))
    mkpath(output_dir)

    pscales = range(pscale_range[1], pscale_range[2], length=n_members)
    filepaths = String[]

    for (i, pscale) in enumerate(pscales)
        filename = joinpath(output_dir, "gemb_pscale_$(round(pscale, digits=3)).mat")
        create_mock_gemb_mat(filename; pscale=pscale)
        push!(filepaths, filename)
    end

    return filepaths
end

"""
    create_mock_gemb_with_known_optimum(filepath;
        true_pscale=1.0,
        n_points=20,
        n_times=48
    )

Create a mock GEMB file where one precipitation scaling parameter is "correct".

This is useful for testing calibration algorithms - we know the true parameter value
and can verify that optimization recovers it.

# Arguments
- `filepath`: Output path for the .mat file
- `true_pscale`: The "true" precipitation scaling parameter (default: 1.0)
- `n_points`: Number of spatial points
- `n_times`: Number of time steps

# Returns
- Tuple of (filepath, true_pscale) for verification
"""
function create_mock_gemb_with_known_optimum(filepath;
    true_pscale=1.0,
    n_points=20,
    n_times=48
)
    create_mock_gemb_mat(filepath;
        n_points=n_points,
        n_times=n_times,
        pscale=true_pscale
    )

    return (filepath, true_pscale)
end

"""
    create_mock_gemb_with_extremes(output_dir)

Create mock GEMB files with extreme parameter values for testing physical constraints.

# Arguments
- `output_dir`: Directory to write test files

# Returns
- Dictionary with keys :low, :medium, :high pointing to file paths
"""
function create_mock_gemb_with_extremes(output_dir)
    mkpath(output_dir)

    files = Dict{Symbol, String}()

    # Very low precipitation scaling (unrealistic)
    files[:low] = joinpath(output_dir, "gemb_pscale_low.mat")
    create_mock_gemb_mat(files[:low]; pscale=0.3)

    # Medium precipitation scaling (realistic)
    files[:medium] = joinpath(output_dir, "gemb_pscale_medium.mat")
    create_mock_gemb_mat(files[:medium]; pscale=1.0)

    # Very high precipitation scaling (unrealistic)
    files[:high] = joinpath(output_dir, "gemb_pscale_high.mat")
    create_mock_gemb_mat(files[:high]; pscale=2.5)

    return files
end
