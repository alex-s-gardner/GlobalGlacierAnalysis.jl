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
- SMB: Surface mass balance (cumulative) [n_points × n_times]
- FAC: Firn air content [n_points × n_times]
- EC: Elevation change (cumulative) [n_points × n_times]
- Accumulation: Snow accumulation (cumulative) [n_points × n_times]
- Runoff: Meltwater runoff (cumulative) [n_points × n_times]
- Melt: Surface melt (cumulative) [n_points × n_times]
- Refreeze: Refreezing (cumulative) [n_points × n_times]
- Rain: Rainfall (cumulative) [n_points × n_times]
- FACtoDepth: FAC to depth conversion factor [scalar]
- H: Surface height [n_points × n_times]
- Ta: Air temperature [n_points × n_times]
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

    # Generate synthetic cumulative variables
    # SMB: -0.5 m/yr * pscale (cumulative, so linearly increasing in magnitude)
    smb_rate = -0.5 * pscale  # m/yr
    smb_cumulative = zeros(n_points, n_times)
    for i in 1:n_points
        for j in 1:n_times
            t_years = (j - 1) / 12.0
            smb_cumulative[i, j] = smb_rate * t_years + 0.1 * randn()
        end
    end
    data["SMB"] = smb_cumulative

    # Accumulation: ~2.0 m/yr * pscale (snow input)
    acc_rate = 2.0 * pscale
    acc_cumulative = zeros(n_points, n_times)
    for i in 1:n_points
        for j in 1:n_times
            t_years = (j - 1) / 12.0
            acc_cumulative[i, j] = acc_rate * t_years + 0.2 * randn()
        end
    end
    data["Accumulation"] = acc_cumulative

    # Runoff: ~2.5 m/yr * pscale (output)
    runoff_rate = 2.5 * pscale
    runoff_cumulative = zeros(n_points, n_times)
    for i in 1:n_points
        for j in 1:n_times
            t_years = (j - 1) / 12.0
            runoff_cumulative[i, j] = runoff_rate * t_years + 0.3 * randn()
        end
    end
    data["Runoff"] = runoff_cumulative

    # Melt: ~2.8 m/yr * pscale
    melt_rate = 2.8 * pscale
    melt_cumulative = zeros(n_points, n_times)
    for i in 1:n_points
        for j in 1:n_times
            t_years = (j - 1) / 12.0
            melt_cumulative[i, j] = melt_rate * t_years + 0.3 * randn()
        end
    end
    data["Melt"] = melt_cumulative

    # Refreeze: ~0.3 m/yr * pscale (refreezing of meltwater)
    refreeze_rate = 0.3 * pscale
    refreeze_cumulative = zeros(n_points, n_times)
    for i in 1:n_points
        for j in 1:n_times
            t_years = (j - 1) / 12.0
            refreeze_cumulative[i, j] = refreeze_rate * t_years + 0.05 * randn()
        end
    end
    data["Refreeze"] = refreeze_cumulative

    # Rain: ~0.2 m/yr * pscale
    rain_rate = 0.2 * pscale
    rain_cumulative = zeros(n_points, n_times)
    for i in 1:n_points
        for j in 1:n_times
            t_years = (j - 1) / 12.0
            rain_cumulative[i, j] = rain_rate * t_years + 0.05 * randn()
        end
    end
    data["Rain"] = rain_cumulative

    # EC: Elevation change (cumulative, related to SMB)
    # EC ≈ SMB / (ice_density / water_density) for solid ice
    ec_cumulative = smb_cumulative ./ 0.91  # Simplified conversion
    data["EC"] = ec_cumulative

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
