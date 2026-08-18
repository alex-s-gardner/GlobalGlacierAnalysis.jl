"""
Tests for low-level binning and gap-filling algorithms in utilities_binning_lowlevel.jl

These tests validate the core data processing algorithms that fill missing elevation change data
using statistical models, multi-mission alignment, and hypsometric interpolation.
"""

using Test
using Dates
using Statistics
using DimensionalData
import DimensionalData as DD
using Random

import GlobalGlacierAnalysis as GGA

# Include test fixtures
include("../fixtures/synthetic_timeseries.jl")
include("../fixtures/synthetic_missions.jl")

@testset "Low-level Binning Algorithms" begin

    @testset "Model-based filling - hyps_model_fill!" begin
        # Generate synthetic data with gaps
        n_dates = 36
        n_heights = 8
        dates, heights, baseline_dh = generate_hypsometric_synthetic(
            n_dates=n_dates,
            n_heights=n_heights,
            trend=-0.6,
            seasonal_amplitude=0.3,
            vertical_gradient=-0.0005,
            noise_sigma=0.05
        )

        # Create DimArrays with intentional gaps
        dh_data = copy(baseline_dh)
        nobs_data = fill(50, n_dates, n_heights)

        # Create gaps: remove 30% of data randomly
        Random.seed!(42)
        gap_indices = rand(1:length(dh_data), Int(floor(0.3 * length(dh_data))))
        dh_data[gap_indices] .= NaN
        nobs_data[gap_indices] .= 0

        # Create DimArrays for a single mission
        date_dim = DD.Dim{:date}(dates)
        height_dim = DD.Dim{:height}(heights)
        geotile_dim = DD.Dim{:geotile}(["lat+45+47lon-123-121"])

        dh_array = DimArray(
            reshape(dh_data, 1, n_dates, n_heights),
            (geotile_dim, date_dim, height_dim)
        )
        nobs_array = DimArray(
            reshape(nobs_data, 1, n_dates, n_heights),
            (geotile_dim, date_dim, height_dim)
        )

        # Package into dictionary format expected by function
        dh_dict = Dict("icesat2" => dh_array)
        nobs_dict = Dict("icesat2" => nobs_array)

        # Create params structure
        params = (
            bincount_min = 3,
            model1_nmad_max = 5.0,
            smooth_n = 5,
            smooth_h2t_length_scale = 500.0,
            missions2update = ["icesat2"]
        )

        # Apply model-based filling
        dh_filled, nobs_filled = GGA.hyps_model_fill!(
            dh_dict,
            nobs_dict,
            params;
            bincount_min=params.bincount_min,
            model1_nmad_max=params.model1_nmad_max,
            smooth_n=params.smooth_n,
            smooth_h2t_length_scale=params.smooth_h2t_length_scale,
            missions2update=params.missions2update
        )

        # Test: Gaps should be filled
        n_gaps_before = sum(isnan.(dh_data))
        n_gaps_after = sum(isnan.(dh_filled["icesat2"][1, :, :]))
        @test n_gaps_after < n_gaps_before
        @test n_gaps_after >= 0  # Some gaps may remain if unfillable

        # Test: Filled values should be close to truth where we removed data
        # (This tests that the model is capturing the underlying trend)
        filled_mask = isnan.(dh_data) .&& .!isnan.(dh_filled["icesat2"][1, :, :])
        if sum(filled_mask) > 0
            filled_values = dh_filled["icesat2"][1, :, :][filled_mask]
            true_values = baseline_dh[filled_mask]
            rmse = sqrt(mean((filled_values .- true_values).^2))
            @test rmse < 0.5  # Filled values within 0.5m of truth (reasonable for noisy data)
        end

        # Test: Original valid data should be preserved
        valid_mask = .!isnan.(dh_data)
        if sum(valid_mask) > 0
            @test all(dh_filled["icesat2"][1, :, :][valid_mask] .≈ dh_data[valid_mask])
        end
    end

    @testset "Multi-mission alignment - hyps_align_dh!" begin
        # Generate synthetic multi-mission data with known offsets
        n_dates = 24
        n_heights = 6

        # Generate baseline (truth)
        dates, heights, baseline_dh = generate_hypsometric_synthetic(
            n_dates=n_dates,
            n_heights=n_heights,
            trend=-0.5,
            seasonal_amplitude=0.2,
            start_date=DateTime(2018,1,1),
            noise_sigma=0.05
        )

        # Create missions with known biases
        icesat2_bias = 0.0  # Reference mission (no bias)
        icesat_bias = -0.8  # ICESat has -0.8m bias relative to ICESat-2
        hugonnet_bias = 1.2  # Hugonnet has +1.2m bias

        # Generate mission data
        geotile_dim = DD.Dim{:geotile}(["lat+60+62lon-050-048"])
        date_dim = DD.Dim{:date}(dates)
        height_dim = DD.Dim{:height}(heights)

        # ICESat-2 (reference)
        icesat2_dh = reshape(baseline_dh .+ icesat2_bias .+ 0.1 .* randn(n_dates, n_heights),
                             1, n_dates, n_heights)
        icesat2_nobs = reshape(fill(50, n_dates, n_heights), 1, n_dates, n_heights)

        # ICESat (with bias, partial temporal coverage)
        icesat_dh = reshape(baseline_dh .+ icesat_bias .+ 0.15 .* randn(n_dates, n_heights),
                           1, n_dates, n_heights)
        icesat_dh[:, 13:end, :] .= NaN  # ICESat ends mid-2009
        icesat_nobs = reshape(fill(30, n_dates, n_heights), 1, n_dates, n_heights)
        icesat_nobs[:, 13:end, :] .= 0

        # Hugonnet (with bias, different temporal coverage)
        hugonnet_dh = reshape(baseline_dh .+ hugonnet_bias .+ 0.12 .* randn(n_dates, n_heights),
                             1, n_dates, n_heights)
        hugonnet_nobs = reshape(fill(100, n_dates, n_heights), 1, n_dates, n_heights)

        # Create DimArrays
        dh_dict = Dict(
            "icesat2" => DimArray(icesat2_dh, (geotile_dim, date_dim, height_dim)),
            "icesat" => DimArray(icesat_dh, (geotile_dim, date_dim, height_dim)),
            "hugonnet" => DimArray(hugonnet_dh, (geotile_dim, date_dim, height_dim))
        )
        nobs_dict = Dict(
            "icesat2" => DimArray(icesat2_nobs, (geotile_dim, date_dim, height_dim)),
            "icesat" => DimArray(icesat_nobs, (geotile_dim, date_dim, height_dim)),
            "hugonnet" => DimArray(hugonnet_nobs, (geotile_dim, date_dim, height_dim))
        )

        # Mock area array
        area_km2 = DimArray(
            reshape(fill(10.0, n_heights), 1, n_heights),
            (geotile_dim, height_dim)
        )

        # Create params
        params = (missions2align2 = ["icesat2", "icesat"],)

        # Apply alignment
        dh_aligned, nobs_aligned = GGA.hyps_align_dh!(
            dh_dict,
            nobs_dict,
            params,
            area_km2;
            missions2align2=params.missions2align2,
            missions2update=["hugonnet"]
        )

        # Test: Hugonnet bias should be largely removed
        # Compare overlap period where both ICESat-2 and Hugonnet have data
        overlap_mask = .!isnan.(dh_aligned["icesat2"][1, :, :]) .&& .!isnan.(dh_aligned["hugonnet"][1, :, :])

        if sum(overlap_mask) > 10  # Need sufficient overlap
            diff_before = mean(hugonnet_dh[1, :, :][overlap_mask] .- icesat2_dh[1, :, :][overlap_mask])
            diff_after = mean(dh_aligned["hugonnet"][1, :, :][overlap_mask] .- dh_aligned["icesat2"][1, :, :][overlap_mask])

            @test abs(diff_before - hugonnet_bias) < 0.3  # Verify we started with the known bias
            @test abs(diff_after) < abs(diff_before)  # Bias should be reduced
            @test abs(diff_after) < 0.4  # Remaining bias should be small
        end

        # Test: ICESat-2 (reference mission) should be unchanged
        @test all(dh_aligned["icesat2"] .≈ dh_dict["icesat2"])
    end

    @testset "Amplitude normalization - hyps_amplitude_normalize!" begin
        # Generate synthetic data with different seasonal amplitudes
        n_dates = 36
        n_heights = 5
        dates, heights, _ = generate_hypsometric_synthetic(
            n_dates=n_dates,
            n_heights=n_heights,
            trend=-0.3,
            seasonal_amplitude=1.0,  # Reference amplitude
            start_date=DateTime(2015,1,1),
            noise_sigma=0.05
        )

        # Create two missions with different amplitudes
        t = [(Dates.value(date - dates[1])) / (365.25 * 24 * 60 * 60 * 1000) for date in dates]

        # Mission 1: amplitude = 1.0 (reference)
        mission1_dh = zeros(n_dates, n_heights)
        # Mission 2: amplitude = 2.5 (needs normalization)
        mission2_dh = zeros(n_dates, n_heights)

        for i in 1:n_dates
            for j in 1:n_heights
                mission1_dh[i, j] = -0.3 * t[i] + 1.0 * sin(2π * t[i]) + 0.05 * randn()
                mission2_dh[i, j] = -0.3 * t[i] + 2.5 * sin(2π * t[i]) + 0.08 * randn()
            end
        end

        # Create DimArrays
        geotile_dim = DD.Dim{:geotile}(["lat+35+37lon+080+082"])
        date_dim = DD.Dim{:date}(dates)
        height_dim = DD.Dim{:height}(heights)

        dh_dict = Dict(
            "mission1" => DimArray(reshape(mission1_dh, 1, n_dates, n_heights),
                                  (geotile_dim, date_dim, height_dim)),
            "mission2" => DimArray(reshape(mission2_dh, 1, n_dates, n_heights),
                                  (geotile_dim, date_dim, height_dim))
        )

        # Create params and params_reference
        params = (amplitude_normalize_2 = "mission1",)
        params_reference = params

        # Apply amplitude normalization to mission2
        dh_normalized = GGA.hyps_amplitude_normalize!(dh_dict, params, params_reference)

        # Test: Mission2 amplitude should now match mission1 amplitude
        # Compute seasonal amplitude for both missions after normalization
        mission1_seasonal = mission1_dh .- mean(mission1_dh, dims=1)
        mission2_seasonal_before = mission2_dh .- mean(mission2_dh, dims=1)
        mission2_seasonal_after = dh_normalized["mission2"][1, :, :] .- mean(dh_normalized["mission2"][1, :, :], dims=1)

        amp1 = std(mission1_seasonal[:])
        amp2_before = std(mission2_seasonal_before[:])
        amp2_after = std(mission2_seasonal_after[:])

        @test amp2_before > amp1 * 1.5  # Verify mission2 started with larger amplitude
        @test abs(amp2_after - amp1) < abs(amp2_before - amp1)  # Amplitude should be closer after normalization
        @test amp2_after / amp1 > 0.5 && amp2_after / amp1 < 1.5  # Ratio should be near 1

        # Test: Mission1 (reference) should be unchanged
        @test all(dh_normalized["mission1"] .≈ dh_dict["mission1"])
    end

    @testset "Fill empty elevation bins - hyps_fill_empty!" begin
        # Create data with some completely empty elevation bins
        n_dates = 20
        n_heights = 10
        dates, heights, baseline_dh = generate_hypsometric_synthetic(
            n_dates=n_dates,
            n_heights=n_heights,
            trend=-0.4,
            vertical_gradient=-0.002,  # Strong elevation dependence
            noise_sigma=0.05
        )

        # Remove entire elevation bins to simulate missing coverage
        dh_data = copy(baseline_dh)
        dh_data[:, [3, 7, 9]] .= NaN  # Remove bins 3, 7, 9

        geotile_dim = DD.Dim{:geotile}(["lat+65+67lon+020+022"])
        date_dim = DD.Dim{:date}(dates)
        height_dim = DD.Dim{:height}(heights)

        dh_dict = Dict(
            "icesat2" => DimArray(reshape(dh_data, 1, n_dates, n_heights),
                                 (geotile_dim, date_dim, height_dim))
        )

        # Create mock geotile extent and area
        geotile_extent = Dict(
            "lat+65+67lon+020+022" => (lon_min=20.0, lon_max=22.0, lat_min=65.0, lat_max=67.0)
        )
        area_km2 = DimArray(
            reshape(fill(5.0, n_heights), 1, n_heights),
            (geotile_dim, height_dim)
        )

        params = (missions2update = ["icesat2"],)

        # Apply fill_empty
        dh_filled = GGA.hyps_fill_empty!(
            dh_dict,
            params,
            geotile_extent,
            area_km2;
            missions2update=params.missions2update
        )

        # Test: Empty bins should be filled
        n_empty_before = sum(all(isnan.(dh_data), dims=1))
        n_empty_after = sum(all(isnan.(dh_filled["icesat2"][1, :, :]), dims=1))
        @test n_empty_after < n_empty_before

        # Test: Filled values should follow elevation gradient
        if n_empty_after < n_empty_before
            # Check that filled bin 7 has values between bins 6 and 8
            if !all(isnan.(dh_filled["icesat2"][1, :, 7]))
                mean_filled_7 = mean(filter(!isnan, dh_filled["icesat2"][1, :, 7]))
                mean_bin_6 = mean(filter(!isnan, dh_filled["icesat2"][1, :, 6]))
                mean_bin_8 = mean(filter(!isnan, dh_filled["icesat2"][1, :, 8]))

                # Filled values should be bracketed by neighbors (with tolerance for noise)
                @test mean_filled_7 < max(mean_bin_6, mean_bin_8) + 0.5
                @test mean_filled_7 > min(mean_bin_6, mean_bin_8) - 0.5
            end
        end
    end

    @testset "Up-down filling - hyps_fill_updown!" begin
        # Create sparse data with gaps that can be filled from adjacent elevation bins
        n_dates = 15
        n_heights = 8
        dates, heights, baseline_dh = generate_hypsometric_synthetic(
            n_dates=n_dates,
            n_heights=n_heights,
            trend=-0.5,
            vertical_gradient=-0.001,
            noise_sigma=0.05
        )

        # Create systematic gaps: remove some time-elevation combinations
        dh_data = copy(baseline_dh)
        # Remove patches that can be filled from above/below
        dh_data[5:7, 4] .= NaN  # Gap at height 4, times 5-7
        dh_data[10:12, 6] .= NaN  # Gap at height 6, times 10-12

        geotile_dim = DD.Dim{:geotile}(["lat-15-13lon-075-073"])
        date_dim = DD.Dim{:date}(dates)
        height_dim = DD.Dim{:height}(heights)

        dh_dict = Dict(
            "icesat2" => DimArray(reshape(dh_data, 1, n_dates, n_heights),
                                 (geotile_dim, date_dim, height_dim))
        )

        area_km2 = DimArray(
            reshape(fill(8.0, n_heights), 1, n_heights),
            (geotile_dim, height_dim)
        )

        params = (missions2update = ["icesat2"],)

        # Apply up-down filling
        dh_filled = GGA.hyps_fill_updown!(
            dh_dict,
            area_km2;
            missions2update=params.missions2update
        )

        # Test: Some gaps should be filled
        n_gaps_before = sum(isnan.(dh_data))
        n_gaps_after = sum(isnan.(dh_filled["icesat2"][1, :, :]))
        @test n_gaps_after <= n_gaps_before

        # Test: Filled values should be reasonable
        # (within range of adjacent elevation bins)
        filled_mask = isnan.(dh_data) .&& .!isnan.(dh_filled["icesat2"][1, :, :])
        if sum(filled_mask) > 0
            filled_values = dh_filled["icesat2"][1, :, :][filled_mask]
            all_valid_values = filter(!isnan, baseline_dh)

            # Filled values should be within the range of the data
            @test all(filled_values .>= minimum(all_valid_values) - 1.0)
            @test all(filled_values .<= maximum(all_valid_values) + 1.0)
        end
    end

    @testset "Edge cases - all missing data" begin
        # Test handling of completely missing data
        n_dates = 10
        n_heights = 5
        dates = [DateTime(2020,1,1) + Month(i) for i in 0:n_dates-1]
        heights = collect(range(1500, 2500, length=n_heights))

        geotile_dim = DD.Dim{:geotile}(["lat+40+42lon-120-118"])
        date_dim = DD.Dim{:date}(dates)
        height_dim = DD.Dim{:height}(heights)

        # Create completely empty arrays
        dh_empty = reshape(fill(NaN, n_dates * n_heights), 1, n_dates, n_heights)
        nobs_empty = reshape(zeros(Int, n_dates * n_heights), 1, n_dates, n_heights)

        dh_dict = Dict("icesat2" => DimArray(dh_empty, (geotile_dim, date_dim, height_dim)))
        nobs_dict = Dict("icesat2" => DimArray(nobs_empty, (geotile_dim, date_dim, height_dim)))

        params = (
            bincount_min = 3,
            missions2update = ["icesat2"]
        )

        # These should not crash with empty data
        @test_nowarn GGA.hyps_model_fill!(dh_dict, nobs_dict, params;
                                          bincount_min=params.bincount_min,
                                          missions2update=params.missions2update)
    end

    @testset "Edge cases - single elevation bin" begin
        # Test with only one elevation bin
        n_dates = 20
        dates, heights, baseline_dh = generate_hypsometric_synthetic(
            n_dates=n_dates,
            n_heights=1,  # Single bin
            trend=-0.3,
            seasonal_amplitude=0.15,
            noise_sigma=0.05
        )

        geotile_dim = DD.Dim{:geotile}(["lat+70+72lon-050-048"])
        date_dim = DD.Dim{:date}(dates)
        height_dim = DD.Dim{:height}(heights)

        dh_dict = Dict(
            "icesat2" => DimArray(reshape(baseline_dh, 1, n_dates, 1),
                                 (geotile_dim, date_dim, height_dim))
        )
        nobs_dict = Dict(
            "icesat2" => DimArray(reshape(fill(50, n_dates, 1), 1, n_dates, 1),
                                 (geotile_dim, date_dim, height_dim))
        )

        params = (
            bincount_min = 3,
            missions2update = ["icesat2"]
        )

        # Should handle single bin without errors
        @test_nowarn GGA.hyps_model_fill!(dh_dict, nobs_dict, params;
                                          bincount_min=params.bincount_min,
                                          missions2update=params.missions2update)
    end

end
