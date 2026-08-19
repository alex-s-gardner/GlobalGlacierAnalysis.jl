"""
Mock external data file generation for testing data readers.

This module generates minimal mock CSV/MAT files that match the format of external
datasets (glacier discharge, GRACE, GlaMBIE) for testing I/O functions without
requiring the full datasets.
"""

using CSV
using DataFrames
using MAT
using Dates

"""
    create_mock_discharge_csv(filepath; n_glaciers=5)

Create a mock glacier discharge CSV file matching Kochtitzky format.

# Arguments
- `filepath`: Output path for the CSV file
- `n_glaciers`: Number of mock glacier entries (default: 5)

# Returns
- Nothing (writes file to disk)

# Format
The CSV has these columns:
- RGIId: RGI glacier identifier
- Year: Observation year
- Discharge_Gt_yr: Ice discharge in Gt/yr
- Uncertainty_Gt_yr: Uncertainty in Gt/yr
"""
function create_mock_discharge_csv(filepath; n_glaciers=5)
    # Generate mock data
    rgi_ids = ["RGI60-05.$(10000+i)" for i in 1:n_glaciers]  # Greenland IDs
    years = 2000:2020
    n_years = length(years)

    # Discharge is kept comfortably away from zero, and the uncertainty is generated as a
    # *fraction* of it (~10%, bounded to 2-40%). Drawing the two independently -- as this used to --
    # let a near-zero discharge produce a relative uncertainty above 100%, so any test bounding
    # `uncertainty / discharge` failed intermittently depending on the RNG.
    n = n_glaciers * n_years
    discharge = max.(1.0, 5.0 .+ 2.0 .* randn(n))
    rel_uncertainty = clamp.(0.10 .+ 0.03 .* randn(n), 0.02, 0.40)

    # Create DataFrame
    df = DataFrame(
        RGIId = repeat(rgi_ids, inner=n_years),
        Year = repeat(years, outer=n_glaciers),
        Discharge_Gt_yr = discharge,
        Uncertainty_Gt_yr = discharge .* rel_uncertainty
    )

    # Write with header lines (Kochtitzky format has metadata).
    #
    # The CSV body is rendered to a String first and then written. Handing the open `io` straight to
    # `CSV.write` discarded the metadata lines already written to that stream, so the file ended up
    # with the column header on line 1 -- and every reader that skipped the documented 14 metadata
    # lines silently parsed data rows as column names.
    metadata = """
    # Glacier discharge data
    # Source: Mock data for testing
    # Citation: Test et al. (2024)
    #
    # Column descriptions:
    # RGIId: RGI glacier identifier
    # Year: Observation year
    # Discharge_Gt_yr: Ice discharge (Gt/yr)
    # Uncertainty_Gt_yr: Uncertainty (Gt/yr)
    #
    # Data start
    # ------------
    #
    #
    """

    open(filepath, "w") do io
        write(io, metadata)          # 14 metadata lines
        write(io, sprint(CSV.write, df))  # header on line 15, data from line 16
    end

    return nothing
end

"""
    create_mock_grace_mat(filepath; n_regions=5, n_times=250)

Create a mock GRACE mass change .mat file.

# Arguments
- `filepath`: Output path for the .mat file
- `n_regions`: Number of RGI regions (default: 5)
- `n_times`: Number of monthly time steps (default: 250 = ~20 years)

# Returns
- Nothing (writes file to disk)

# Format
GRACE .mat files contain:
- region_codes: RGI region abbreviations (e.g., "GRE", "PAT")
- time: Decimal years
- mass_change_Gt: Mass change time series [n_regions × n_times]
- uncertainty_Gt: Uncertainty time series [n_regions × n_times]
"""
function create_mock_grace_mat(filepath; n_regions=5, n_times=250)
    # Define mock regions (using actual GRACE region codes)
    region_codes = ["GRE", "PAT", "NAS", "ICE", "TRP"]

    # Create time vector (monthly from 2003 to ~2023)
    start_year = 2003.0
    times = start_year .+ collect(0:n_times-1) ./ 12.0

    # Generate synthetic mass change data
    # Each region has different trend + seasonal signal + noise
    mass_change = zeros(n_regions, n_times)
    uncertainty = zeros(n_regions, n_times)

    for i in 1:n_regions
        # Different trend per region (-300 to -50 Gt/yr)
        trend = -300.0 + 50.0 * i

        # Seasonal amplitude (20-40 Gt)
        seasonal_amp = 20.0 + 5.0 * i

        # Noise level (10-20 Gt)
        noise_sigma = 10.0 + 2.0 * i

        for j in 1:n_times
            t_years = (j - 1) / 12.0
            seasonal = seasonal_amp * sin(2π * (times[j] - 0.25))  # Peak in summer
            noise = noise_sigma * randn()

            # Cumulative mass change
            mass_change[i, j] = trend * t_years + seasonal + noise

            # Uncertainty (varies slightly with time)
            uncertainty[i, j] = noise_sigma + 3.0 * sin(2π * times[j])
        end
    end

    # Create MATLAB structure
    data = Dict{String, Any}(
        "region_codes" => region_codes,
        "time" => times,
        "mass_change_Gt" => mass_change,
        "uncertainty_Gt" => abs.(uncertainty)  # Ensure positive
    )

    # Write to .mat file
    matwrite(filepath, data)

    return nothing
end

"""
    create_mock_glambie_csv(filepath; n_years=25)

Create a mock GlaMBIE mass balance CSV file.

# Arguments
- `filepath`: Output path for the CSV file
- `n_years`: Number of years of data (default: 25 = 2000-2024)

# Returns
- Nothing (writes file to disk)

# Format
GlaMBIE CSV has these columns:
- Year: Observation year (2000-2024)
- RGI_Region: RGI region number (1-19, 98, 99)
- Mass_Balance_Gt: Cumulative mass balance (Gt)
- Uncertainty_Gt: Cumulative uncertainty (Gt)
"""
function create_mock_glambie_csv(filepath; n_years=25)
    years = 2000:(2000+n_years-1)
    rgi_regions = vcat(collect(1:19), [98, 99])  # All RGI regions + global
    n_regions = length(rgi_regions)

    # Generate data
    records = []
    for region in rgi_regions
        # Different trend per region
        if region == 98  # Global (sum of all)
            trend = -400.0  # Gt/yr
        elseif region == 99  # Global excluding Greenland/Antarctica
            trend = -150.0  # Gt/yr
        else
            trend = -20.0 * (region / 10)  # Varies by region
        end

        for (i, year) in enumerate(years)
            t_years = i - 1
            # Cumulative mass balance
            mass_balance = trend * t_years + 10.0 * randn()
            uncertainty = abs(5.0 * sqrt(t_years) + 2.0 * randn())

            push!(records, (
                Year = year,
                RGI_Region = region,
                Mass_Balance_Gt = mass_balance,
                Uncertainty_Gt = uncertainty
            ))
        end
    end

    df = DataFrame(records)

    # Write CSV
    CSV.write(filepath, df)

    return nothing
end

"""
    create_mock_data_directory(base_dir)

Create a directory structure with all mock data files for testing.

# Arguments
- `base_dir`: Base directory to create mock data files

# Returns
- Dictionary with paths to all created files
"""
function create_mock_data_directory(base_dir)
    mkpath(base_dir)

    paths = Dict{Symbol, String}()

    # Discharge data
    discharge_dir = joinpath(base_dir, "discharge")
    mkpath(discharge_dir)
    paths[:discharge_nh] = joinpath(discharge_dir, "kochtitzky_nh.csv")
    paths[:discharge_npi] = joinpath(discharge_dir, "kochtitzky_npi.csv")
    paths[:discharge_spi] = joinpath(discharge_dir, "kochtitzky_spi.csv")

    create_mock_discharge_csv(paths[:discharge_nh]; n_glaciers=10)
    create_mock_discharge_csv(paths[:discharge_npi]; n_glaciers=5)
    create_mock_discharge_csv(paths[:discharge_spi]; n_glaciers=5)

    # GRACE data
    grace_dir = joinpath(base_dir, "grace")
    mkpath(grace_dir)
    paths[:grace_rgi] = joinpath(grace_dir, "grace_rgi_mass_change.mat")

    create_mock_grace_mat(paths[:grace_rgi]; n_regions=5, n_times=250)

    # GlaMBIE data
    glambie_dir = joinpath(base_dir, "glambie")
    mkpath(glambie_dir)
    paths[:glambie_2024] = joinpath(glambie_dir, "glambie_2024_mass_balance.csv")

    create_mock_glambie_csv(paths[:glambie_2024]; n_years=25)

    return paths
end

"""
    create_mock_zemp2019_csv(filepath)

Create a mock Zemp et al. 2019 glacier mass balance CSV file.

# Arguments
- `filepath`: Output path for the CSV file

# Returns
- Nothing (writes file to disk)

# Format
Zemp 2019 format:
- YEAR: Observation year
- MEAN_BALANCE: Mean specific mass balance (m w.e./yr)
- LOWER_BOUND: Lower uncertainty bound
- UPPER_BOUND: Upper uncertainty bound
"""
function create_mock_zemp2019_csv(filepath)
    years = 1961:2019
    n_years = length(years)

    # Generate realistic mass balance trend
    # Becoming more negative over time
    mean_balance = -0.2 .- 0.015 .* (years .- 1961) .+ 0.1 .* randn(n_years)
    uncertainty = 0.05 .+ 0.002 .* (years .- 1961)

    df = DataFrame(
        YEAR = years,
        MEAN_BALANCE = mean_balance,
        LOWER_BOUND = mean_balance .- uncertainty,
        UPPER_BOUND = mean_balance .+ uncertainty
    )

    # Write with header
    open(filepath, "w") do io
        write(io, "# Global glacier mass balance data\n")
        write(io, "# Source: Zemp et al. (2019) - Mock version\n")
        write(io, "# Units: m w.e. yr-1 (meters water equivalent per year)\n")
        write(io, "#\n")
        CSV.write(io, df)
    end

    return nothing
end
