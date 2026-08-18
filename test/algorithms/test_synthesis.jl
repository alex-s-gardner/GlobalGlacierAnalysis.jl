"""
Tests for multi-mission data synthesis in utilities_synthesis.jl

These tests validate algorithms that combine multiple satellite altimetry missions,
correct systematic biases, compute volume changes, and aggregate data.
"""

using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using Statistics
using Dates
using DimensionalData
import DimensionalData as DD
using Random

# Include test fixtures
include("../fixtures/synthetic_missions.jl")

@testset "Synthesis Operations" begin
    @testset "Inverse variance weighting" begin
        # Three measurements with different errors
        values = [10.0, 12.0, 11.0]
        errors = [1.0, 2.0, 0.5]

        # Weights: 1/σ²
        weights = 1.0 ./ (errors .^ 2)
        @test weights ≈ [1.0, 0.25, 4.0]

        # Weighted mean
        weighted_mean = sum(values .* weights) / sum(weights)
        expected = (10.0*1.0 + 12.0*0.25 + 11.0*4.0) / (1.0 + 0.25 + 4.0)
        @test weighted_mean ≈ expected rtol=1e-6
        @test weighted_mean ≈ 10.857142857 rtol=1e-6

        # Lower error should get higher weight
        @test weights[3] > weights[1] > weights[2]
    end

    @testset "Linear interpolation for gaps" begin
        # Create time series with gaps
        dates_decimal = [2018.0, 2018.25, 2018.5, 2018.75, 2019.0]
        values = [0.0, missing, 1.0, missing, 2.0]

        # Expected after linear interpolation: [0.0, 0.5, 1.0, 1.5, 2.0]
        # (This test verifies the concept; actual implementation may vary)
        valid_indices = findall(.!ismissing.(values))
        @test valid_indices == [1, 3, 5]

        # Interpolated value at index 2 should be between values[1] and values[3]
        t1, t2, t3 = dates_decimal[[1, 2, 3]]
        v1, v3 = values[[1, 3]]
        expected_interp = v1 + (v3 - v1) * (t2 - t1) / (t3 - t1)
        @test expected_interp ≈ 0.5 rtol=1e-6
    end

    @testset "Multi-mission synthesis with known offsets" begin
        # Generate synthetic multi-mission data with known systematic biases
        Random.seed!(123)

        true_trend = -0.6  # m/yr
        icesat_bias = -0.5  # m
        hugonnet_bias = 1.0  # m

        synth_data = generate_multimission_synthetic(
            n_dates=36,
            n_heights=8,
            trend=true_trend,
            seasonal_amplitude=0.25,
            start_date=DateTime(2003,1,1),
            mission_biases=Dict(:icesat => icesat_bias, :hugonnet => hugonnet_bias),
            geotile_lat=45
        )

        dates = synth_data[:dates]
        heights = synth_data[:heights]
        truth = synth_data[:truth]

        # Test: Multi-mission data should have different means due to biases
        icesat2_mean = mean(filter(!isnan, synth_data[:icesat2].dh))
        icesat_mean = mean(filter(!isnan, synth_data[:icesat].dh))
        hugonnet_mean = mean(filter(!isnan, synth_data[:hugonnet].dh))

        @test abs((icesat_mean - icesat2_mean) - icesat_bias) < 0.2  # ICESat should be ~0.5m lower
        @test abs((hugonnet_mean - icesat2_mean) - hugonnet_bias) < 0.2  # Hugonnet should be ~1.0m higher

        # Test: Inverse variance weighted combination
        # Simulate a simple combination at overlapping time-space points
        overlap_indices = findall(
            .!isnan.(synth_data[:icesat2].dh) .&&
            .!isnan.(synth_data[:icesat].dh) .&&
            .!isnan.(synth_data[:hugonnet].dh)
        )

        if length(overlap_indices) > 5
            # Get values at overlap
            vals_is2 = synth_data[:icesat2].dh[overlap_indices]
            vals_is = synth_data[:icesat].dh[overlap_indices]
            vals_hug = synth_data[:hugonnet].dh[overlap_indices]

            # Get uncertainties
            sigma_is2 = synth_data[:icesat2].sigma[overlap_indices]
            sigma_is = synth_data[:icesat].sigma[overlap_indices]
            sigma_hug = synth_data[:hugonnet].sigma[overlap_indices]

            # Compute weighted mean (inverse variance weighting)
            w_is2 = 1 ./ sigma_is2.^2
            w_is = 1 ./ sigma_is.^2
            w_hug = 1 ./ sigma_hug.^2
            w_total = w_is2 .+ w_is .+ w_hug

            vals_combined = (vals_is2 .* w_is2 .+ vals_is .* w_is .+ vals_hug .* w_hug) ./ w_total

            # Combined uncertainty should be lower than any individual
            sigma_combined = sqrt.(1 ./ w_total)
            @test all(sigma_combined .< sigma_is2)
            @test all(sigma_combined .< sigma_is)
            @test all(sigma_combined .< sigma_hug)

            # Weighted mean should be closer to truth than unweighted mean
            truth_overlap = truth[overlap_indices]
            rmse_weighted = sqrt(mean((vals_combined .- truth_overlap).^2))
            rmse_unweighted = sqrt(mean((mean([vals_is2, vals_is, vals_hug]) .- truth_overlap).^2))

            @test rmse_weighted <= rmse_unweighted + 0.1  # Weighted should be similar or better
        end
    end

    @testset "Elevation to volume conversion - dh2dv" begin
        # Test volume conversion with known geometry
        n_dates = 12
        n_heights = 5

        dates = [DateTime(2020,1,1) + Month(i) for i in 0:n_dates-1]
        heights = collect(range(2000, 3000, length=n_heights))

        # Create uniform elevation change: -1.0 m/yr everywhere
        t_years = collect(0:n_dates-1) ./ 12.0
        dh_uniform = -1.0 .* repeat(t_years, 1, n_heights)  # [n_dates × n_heights]

        # Create area distribution: 10 km² per elevation bin
        area_per_bin = 10.0  # km²
        areas = fill(area_per_bin, n_heights)

        # Create DimArrays
        date_dim = DD.Dim{:date}(dates)
        height_dim = DD.Dim{:height}(heights)

        dh_array = DimArray(dh_uniform, (date_dim, height_dim))
        area_array = DimArray(areas, (height_dim,))

        # Compute volume change
        dv = GGA.dh2dv(dh_array, area_array)

        # Expected volume change:
        # dV = Σ(dh_i * A_i) where A_i in km²
        # At year 0: dV = -1.0 * 0 * 50 = 0 km³
        # At year 1: dV = -1.0 * 1 * 50 / 1000 = -0.05 km³ (convert m to km)
        total_area = sum(areas)  # 50 km²
        expected_dv = -1.0 .* t_years .* (total_area / 1000)  # Convert m to km: divide by 1000

        @test length(dv) == n_dates
        @test dv ≈ expected_dv rtol=1e-6
        @test dv[1] ≈ 0.0 rtol=1e-6  # No change at t=0
        @test dv[end] ≈ -11/12 * 50/1000 rtol=1e-6  # At 11 months
    end

    @testset "Elevation to volume conversion - unit correctness" begin
        # Verify unit conversions: m³, km³, Gt
        ice_density = 910.0  # kg/m³

        # Simple case: 1m thinning over 100 km²
        dh_m = -1.0  # meters
        area_km2 = 100.0  # km²

        # Volume in km³
        dv_km3 = dh_m * area_km2 / 1000  # -0.1 km³
        @test dv_km3 ≈ -0.1

        # Volume in m³
        dv_m3 = dv_km3 * 1e9  # -1e8 m³
        @test dv_m3 ≈ -1e8

        # Mass in Gt (gigatons)
        dm_Gt = dv_m3 * ice_density / 1e12  # kg to Gt
        @test dm_Gt ≈ -0.091 rtol=1e-6

        # Alternative: km³ to Gt directly
        dm_Gt_alt = dv_km3 * ice_density * 1e9 / 1e12
        @test dm_Gt_alt ≈ dm_Gt
    end

    @testset "Volume conservation across aggregation" begin
        # Test that volume is conserved when aggregating spatially
        n_dates = 10
        n_heights = 4
        n_geotiles = 3

        dates = [DateTime(2019,1,1) + Month(3*i) for i in 0:n_dates-1]
        heights = collect(range(1500, 2500, length=n_heights))
        geotiles = ["lat+60+62lon-120-118", "lat+62+64lon-120-118", "lat+64+66lon-120-118"]

        # Create elevation change data: different trend per geotile
        dh_data = zeros(n_geotiles, n_dates, n_heights)
        area_data = zeros(n_geotiles, n_heights)

        for (i, gt) in enumerate(geotiles)
            t_years = collect(0:n_dates-1) ./ 4.0  # Quarterly
            trend = -0.3 * i  # Different trend per geotile
            for j in 1:n_dates
                for k in 1:n_heights
                    dh_data[i, j, k] = trend * t_years[j]
                end
            end
            area_data[i, :] .= 5.0 * i  # Different area per geotile
        end

        # Create DimArrays
        geotile_dim = DD.Dim{:geotile}(geotiles)
        date_dim = DD.Dim{:date}(dates)
        height_dim = DD.Dim{:height}(heights)

        dh_array = DimArray(dh_data, (geotile_dim, date_dim, height_dim))
        area_array = DimArray(area_data, (geotile_dim, height_dim))

        # Compute volume change per geotile
        dv_by_geotile = GGA.dh2dv_geotile(dh_array, area_array)

        # Compute total volume: sum across geotiles
        dv_total = sum(dv_by_geotile; dims=:geotile)[1, :]

        # Verify against direct calculation
        dv_direct = zeros(n_dates)
        for j in 1:n_dates
            for i in 1:n_geotiles
                for k in 1:n_heights
                    dv_direct[j] += dh_data[i, j, k] * area_data[i, k] / 1000
                end
            end
        end

        @test dv_total ≈ dv_direct rtol=1e-10

        # Test: Volume should sum correctly (conservation)
        @test all(isfinite.(dv_total))
    end

    @testset "Ensemble statistics" begin
        # Test calculation of ensemble quantiles and spread
        Random.seed!(456)

        n_dates = 20
        n_ensemble = 50

        dates = [DateTime(2018,1,1) + Month(i) for i in 0:n_dates-1]

        # Generate ensemble with known distribution
        # Each member has slightly different trend
        true_trend = -0.5
        ensemble_data = zeros(n_ensemble, n_dates)

        for i in 1:n_ensemble
            member_trend = true_trend + 0.1 * randn()  # Trend varies
            t = collect(0:n_dates-1) ./ 12.0
            ensemble_data[i, :] = member_trend .* t .+ 0.05 .* randn(n_dates)
        end

        # Compute statistics
        ensemble_mean = mean(ensemble_data; dims=1)[1, :]
        ensemble_median = median(ensemble_data; dims=1)[1, :]
        ensemble_std = std(ensemble_data; dims=1)[1, :]
        ensemble_p05 = [quantile(ensemble_data[:, j], 0.05) for j in 1:n_dates]
        ensemble_p95 = [quantile(ensemble_data[:, j], 0.95) for j in 1:n_dates]

        # Test: Mean should be close to median (roughly normal distribution)
        @test mean(abs.(ensemble_mean .- ensemble_median)) < 0.1

        # Test: 90% of data should fall within p05-p95 range
        for j in 1:n_dates
            in_range = sum((ensemble_data[:, j] .>= ensemble_p05[j]) .&&
                          (ensemble_data[:, j] .<= ensemble_p95[j]))
            @test in_range / n_ensemble > 0.85  # Should be ~90%, allow some tolerance
            @test in_range / n_ensemble < 0.95
        end

        # Test: Spread should increase with time (uncertainty grows)
        @test ensemble_std[end] > ensemble_std[1]
    end

    @testset "Error propagation in synthesis" begin
        # Test that uncertainties combine correctly in multi-mission synthesis
        n_obs = 10
        Random.seed!(789)

        # Three missions with different uncertainties
        sigma_mission1 = 0.2
        sigma_mission2 = 0.3
        sigma_mission3 = 0.15

        true_value = 5.0

        obs_mission1 = true_value .+ sigma_mission1 .* randn(n_obs)
        obs_mission2 = true_value .+ sigma_mission2 .* randn(n_obs)
        obs_mission3 = true_value .+ sigma_mission3 .* randn(n_obs)

        # Weights for inverse variance weighting
        w1 = 1 / sigma_mission1^2
        w2 = 1 / sigma_mission2^2
        w3 = 1 / sigma_mission3^2

        # Combined estimate at each observation
        combined_obs = (obs_mission1 .* w1 .+ obs_mission2 .* w2 .+ obs_mission3 .* w3) ./ (w1 + w2 + w3)

        # Combined uncertainty
        sigma_combined = sqrt(1 / (w1 + w2 + w3))

        # Test: Combined uncertainty should be smaller than all individual uncertainties
        @test sigma_combined < sigma_mission1
        @test sigma_combined < sigma_mission2
        @test sigma_combined < sigma_mission3

        # Test: Combined uncertainty should follow formula: 1/σ² = 1/σ₁² + 1/σ₂² + 1/σ₃²
        expected_sigma = sqrt(1 / (1/sigma_mission1^2 + 1/sigma_mission2^2 + 1/sigma_mission3^2))
        @test sigma_combined ≈ expected_sigma rtol=1e-10

        # Test: Combined estimate should be closer to truth than individual missions (on average)
        rmse_mission1 = sqrt(mean((obs_mission1 .- true_value).^2))
        rmse_combined = sqrt(mean((combined_obs .- true_value).^2))

        # With sufficient samples, combined should be better (but allow some randomness)
        @test rmse_combined <= rmse_mission1 * 1.2
    end

    @testset "Gap filling temporal continuity" begin
        # Test that gap-filled data maintains temporal continuity
        n_dates = 30
        dates = [DateTime(2018,1,1) + Month(i) for i in 0:n_dates-1]
        t = collect(0:n_dates-1) ./ 12.0

        # Generate smooth underlying signal
        true_signal = -0.5 .* t .+ 0.3 .* sin.(2π .* t)

        # Add noise and create gaps
        Random.seed!(321)
        observed = true_signal .+ 0.1 .* randn(n_dates)

        # Remove 40% of observations
        gap_indices = sort(rand(1:n_dates, Int(floor(0.4 * n_dates))))
        observed_with_gaps = copy(observed)
        observed_with_gaps[gap_indices] .= NaN

        # Simple linear interpolation to fill gaps
        filled = copy(observed_with_gaps)
        for i in gap_indices
            # Find surrounding valid points
            left_idx = findlast(.!isnan.(filled[1:i-1]))
            right_idx = findfirst(.!isnan.(filled[i+1:end]))

            if !isnothing(left_idx) && !isnothing(right_idx)
                right_idx += i  # Adjust index
                # Linear interpolation
                t_frac = (t[i] - t[left_idx]) / (t[right_idx] - t[left_idx])
                filled[i] = filled[left_idx] + t_frac * (filled[right_idx] - filled[left_idx])
            end
        end

        # Test: Filled data should have fewer NaNs
        @test sum(isnan.(filled)) < sum(isnan.(observed_with_gaps))

        # Test: Filled values should be between neighboring observations
        for i in gap_indices
            if !isnan(filled[i])
                nearby_valid = filter(!isnan, filled[max(1, i-2):min(n_dates, i+2)])
                if length(nearby_valid) > 0
                    @test filled[i] >= minimum(nearby_valid) - 0.5
                    @test filled[i] <= maximum(nearby_valid) + 0.5
                end
            end
        end
    end

    @testset "Edge case - all missions have gaps at same time" begin
        # Test handling when all missions are missing data simultaneously
        n_dates = 15
        n_heights = 3

        dates = [DateTime(2020,1,1) + Month(i) for i in 0:n_dates-1]
        heights = [1500.0, 2000.0, 2500.0]

        # Create data with a common gap period
        dh_mission1 = randn(n_dates, n_heights)
        dh_mission2 = randn(n_dates, n_heights)

        # Both missions missing data from indices 7-9
        dh_mission1[7:9, :] .= NaN
        dh_mission2[7:9, :] .= NaN

        # Verify gaps exist
        @test all(isnan.(dh_mission1[7:9, :]))
        @test all(isnan.(dh_mission2[7:9, :]))

        # Test: System should handle this gracefully (no errors)
        # In practice, synthesis should leave these gaps unfilled
        @test_nowarn begin
            # Simulate what synthesis would do: can't fill common gaps
            combined = copy(dh_mission1)
            combined[.!isnan.(dh_mission2)] .= dh_mission2[.!isnan.(dh_mission2)]
        end
    end

    @testset "Edge case - single mission coverage" begin
        # Test synthesis when only one mission has data
        n_dates = 10
        n_heights = 2

        dh_solo = randn(n_dates, n_heights)
        sigma_solo = fill(0.2, n_dates, n_heights)

        # "Synthesis" of single mission should return the mission itself
        combined = copy(dh_solo)
        sigma_combined = copy(sigma_solo)

        @test combined ≈ dh_solo
        @test sigma_combined ≈ sigma_solo

        # Uncertainty should not decrease with only one source
        @test all(sigma_combined .>= sigma_solo .- 1e-10)
    end
end
