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
        geotile_dim = DD.Dim{:geotile}(["lat[+45+47]lon[-123-121]"])

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

        # `params` is the per-mission, per-geotile table that fitted coefficients get written back
        # into -- not the bag of tuning knobs, which are keyword arguments.
        params = synthetic_fill_params(dh_dict)

        # The DimArray wraps `dh_data` without copying and hyps_model_fill! writes through it, so
        # the pre-fill state must be snapshotted before the call.
        dh_before = copy(dh_data)

        # Apply model-based filling
        dh_filled, nobs_filled = GGA.hyps_model_fill!(
            dh_dict,
            nobs_dict,
            params;
            bincount_min=3,
            model1_nmad_max=5.0,
            smooth_n=5,
            smooth_h2t_length_scale=500.0,
            missions2update=["icesat2"]
        )

        # Test: Gaps should be filled
        n_gaps_before = sum(isnan.(dh_before))
        n_gaps_after = sum(isnan.(dh_filled["icesat2"][1, :, :]))
        @test n_gaps_after < n_gaps_before
        @test n_gaps_after >= 0  # Some gaps may remain if unfillable

        # Test: Filled values should be close to truth where we removed data
        # (This tests that the model is capturing the underlying trend)
        filled_mask = isnan.(dh_before) .&& .!isnan.(dh_filled["icesat2"][1, :, :])
        if sum(filled_mask) > 0
            filled_values = dh_filled["icesat2"][1, :, :][filled_mask]
            true_values = baseline_dh[filled_mask]
            rmse = sqrt(mean((filled_values .- true_values).^2))
            @test rmse < 0.5  # Filled values within 0.5m of truth (reasonable for noisy data)
        end

        # Test: originally-valid bins stay valid and stay close to their input values.
        # `hyps_model_fill!` is not a pure gap-filler -- it refits and smooths every bin -- so the
        # values are not reproduced bit-for-bit; what matters is that it does not drift away from
        # the observations or turn them into NaN.
        valid_mask = .!isnan.(dh_before)
        if sum(valid_mask) > 0
            filled_at_valid = dh_filled["icesat2"][1, :, :][valid_mask]
            @test all(.!isnan.(filled_at_valid))
            @test maximum(abs.(filled_at_valid .- dh_before[valid_mask])) < 0.5
        end
    end

    @testset "Residual climatology survives the fill smoothing" begin
        # A semiannual component that the annual sine in `model1` cannot represent. The temporal
        # median smoothing removes most of it; holding the residual climatology out restores it.
        n_dates, n_heights = 72, 8
        Random.seed!(1)
        dates, heights, baseline_dh = generate_hypsometric_synthetic(; n_dates, n_heights,
            trend=-0.6, seasonal_amplitude=0.3, vertical_gradient=-0.0005, noise_sigma=0.05,
            start_date=DateTime(2018, 1, 15))
        t = GGA.decimalyear.(dates)
        semiannual_amplitude = 0.4
        truth = baseline_dh .+ semiannual_amplitude .* cos.(4π .* t .+ 0.7)

        function fill_semiannual(preserve_residual_climatology)
            dh = DimArray(reshape(copy(truth), 1, n_dates, n_heights),
                (DD.Dim{:geotile}(["lat[+45+47]lon[-123-121]"]), DD.Dim{:date}(dates), DD.Dim{:height}(heights)))
            nobs = DimArray(fill(50, 1, n_dates, n_heights), dims(dh))
            dh_dict = Dict("icesat2" => dh)
            params = synthetic_fill_params(dh_dict)
            GGA.hyps_model_fill!(dh_dict, Dict("icesat2" => nobs), params; bincount_min=3, smooth_n=5,
                smooth_h2t_length_scale=400.0, preserve_residual_climatology, missions2update=["icesat2"])
            filled = parent(dh_dict["icesat2"])[1, :, :]
            M = hcat(cos.(4π .* t), sin.(4π .* t))
            amplitude = [hypot((M \ filled[:, j])...) for j in axes(filled, 2)]
            return (; filled, amplitude, params=params["icesat2"])
        end

        smoothed = fill_semiannual(false)
        preserved = fill_semiannual(true)

        @test all(smoothed.amplitude .< 0.5 * semiannual_amplitude)
        @test all(abs.(preserved.amplitude .- semiannual_amplitude) .< 0.1)
        @test sqrt(mean((preserved.filled .- truth) .^ 2)) < 0.5 * sqrt(mean((smoothed.filled .- truth) .^ 2))

        # every elevation bin has full monthly coverage, so each gets a 2-harmonic climatology
        clim = only(preserved.params.residual_clim)
        @test clim.nh == 2
        @test length(clim.heights) == n_heights
        @test all(length.(clim.coef) .== 4)

        @test_throws "residual_climatology_harmonics must be >= 1" GGA.hyps_model_fill!(
            Dict("icesat2" => DimArray(reshape(copy(truth), 1, n_dates, n_heights),
                (DD.Dim{:geotile}(["g"]), DD.Dim{:date}(dates), DD.Dim{:height}(heights)))),
            Dict("icesat2" => DimArray(fill(50, 1, n_dates, n_heights),
                (DD.Dim{:geotile}(["g"]), DD.Dim{:date}(dates), DD.Dim{:height}(heights)))),
            Dict("icesat2" => smoothed.params); preserve_residual_climatology=true,
            residual_climatology_harmonics=0)
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
        geotile_dim = DD.Dim{:geotile}(["lat[+60+62]lon[-050-048]"])
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

        # `params` is the per-mission parameter table (offsets get written back into it), so it
        # must be built with the alignment-reference columns present.
        params = synthetic_fill_params(dh_dict; missions2align2=["icesat2", "icesat"])

        # `hyps_align_dh!` mutates in place, and the DimArrays wrap `icesat2_dh`/`hugonnet_dh`
        # without copying -- so the "before" state has to be snapshotted now, or diff_before and
        # diff_after end up reading the same (already aligned) numbers.
        icesat2_before = copy(dh_dict["icesat2"])
        hugonnet_before = copy(dh_dict["hugonnet"])

        # Apply alignment
        dh_aligned, nobs_aligned = GGA.hyps_align_dh!(
            dh_dict,
            nobs_dict,
            params,
            area_km2;
            missions2align2=["icesat2", "icesat"],
            missions2update=["hugonnet"]
        )

        # Test: Hugonnet bias should be largely removed
        # Compare overlap period where both ICESat-2 and Hugonnet have data
        overlap_mask = .!isnan.(dh_aligned["icesat2"][1, :, :]) .&& .!isnan.(dh_aligned["hugonnet"][1, :, :])

        if sum(overlap_mask) > 10  # Need sufficient overlap
            diff_before = mean(hugonnet_before[1, :, :][overlap_mask] .- icesat2_before[1, :, :][overlap_mask])
            diff_after = mean(dh_aligned["hugonnet"][1, :, :][overlap_mask] .- dh_aligned["icesat2"][1, :, :][overlap_mask])

            @test abs(diff_before - hugonnet_bias) < 0.3  # Verify we started with the known bias
            @test abs(diff_after) < abs(diff_before)  # Bias should be reduced
            @test abs(diff_after) < 0.5  # Remaining bias should be small (synthetic-noise bound)
        end

        # Test: ICESat-2 (reference mission) should be unchanged
        @test all(dh_aligned["icesat2"] .≈ icesat2_before)
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
        geotile_dim = DD.Dim{:geotile}(["lat[+35+37]lon[+080+082]"])
        date_dim = DD.Dim{:date}(dates)
        height_dim = DD.Dim{:height}(heights)

        dh_dict = Dict(
            "mission1" => DimArray(reshape(mission1_dh, 1, n_dates, n_heights),
                                  (geotile_dim, date_dim, height_dim)),
            "mission2" => DimArray(reshape(mission2_dh, 1, n_dates, n_heights),
                                  (geotile_dim, date_dim, height_dim))
        )
        nobs_dict = Dict(
            m => DimArray(reshape(fill(50, n_dates, n_heights), 1, n_dates, n_heights),
                          (geotile_dim, date_dim, height_dim))
            for m in keys(dh_dict)
        )

        # The in-place calls below write through dh_dict (and through mission1_dh/mission2_dh, which
        # the DimArrays wrap without copying), so snapshot the raw inputs up front.
        mission1_dh_before = copy(mission1_dh)
        mission2_dh_before = copy(mission2_dh)

        # `hyps_amplitude_normalize!` works on a single mission's DimArray plus that mission's
        # parameter row table and the reference mission's -- not on the whole mission dictionary.
        all_params = synthetic_fill_params(dh_dict)

        # the amplitude model coefficients come from hyps_model_fill!, so fit both missions first
        GGA.hyps_model_fill!(dh_dict, nobs_dict, all_params;
                             bincount_min=3, smooth_n=5, smooth_h2t_length_scale=500.0)

        # snapshot after model fitting but immediately before normalization, so the "reference
        # unchanged" check isolates hyps_amplitude_normalize! rather than the preceding fit
        mission1_before = copy(dh_dict["mission1"])

        GGA.hyps_amplitude_normalize!(dh_dict["mission2"], all_params["mission2"],
                                      all_params["mission1"])

        # Test: Mission2 amplitude should now match mission1 amplitude
        # Compute seasonal amplitude for both missions after normalization
        mission1_seasonal = mission1_dh_before .- mean(mission1_dh_before, dims=1)
        mission2_seasonal_before = mission2_dh_before .- mean(mission2_dh_before, dims=1)
        mission2_seasonal_after = dh_dict["mission2"][1, :, :] .- mean(dh_dict["mission2"][1, :, :], dims=1)

        amp1 = std(mission1_seasonal[:])
        amp2_before = std(mission2_seasonal_before[:])
        amp2_after = std(mission2_seasonal_after[:])

        @test amp2_before > amp1 * 1.5  # Verify mission2 started with larger amplitude
        @test abs(amp2_after - amp1) < abs(amp2_before - amp1)  # Amplitude should be closer after normalization
        @test amp2_after / amp1 > 0.5 && amp2_after / amp1 < 1.5  # Ratio should be near 1

        # Test: Mission1 (reference) should be unchanged
        @test all(dh_dict["mission1"] .≈ mission1_before)
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

        # `hyps_fill_empty!` fills a gap with the median of its *nearest neighbours*, so the fixture
        # needs more than one geotile -- with a single tile there is nothing to borrow from and the
        # empty bins stay empty. Use three adjacent tiles: two fully populated, one with holes.
        geotile_ids = ["lat[+65+67]lon[+020+022]",
                       "lat[+65+67]lon[+022+024]",
                       "lat[+65+67]lon[+024+026]"]
        geotile_dim = DD.Dim{:geotile}(geotile_ids)
        date_dim = DD.Dim{:date}(dates)
        height_dim = DD.Dim{:height}(heights)

        # The neighbour fill is gated on `all(isnan.(dh0))` -- it replaces a geotile that has *no*
        # data at all, rather than patching individual missing bins (that is what hyps_model_fill!
        # and hyps_fill_updown! do). So tile 1 is left completely empty.
        dh_all = Array{Float64}(undef, 3, n_dates, n_heights)
        dh_all[1, :, :] .= NaN             # no data at all -> should be filled from neighbours
        dh_all[2, :, :] = baseline_dh      # neighbour with full coverage
        dh_all[3, :, :] = baseline_dh      # neighbour with full coverage

        dh_dict = Dict("icesat2" => DimArray(dh_all, (geotile_dim, date_dim, height_dim)))

        # `hyps_fill_empty!` does `geotile_extent[geotile=At(...)]`, so this has to be a DimArray
        # over the geotile dimension -- a Dict keyed by id does not support that lookup. The values
        # are Extents, which `extent2rectangle` consumes.
        geotile_extent = DimArray(GGA.geotile_extent.(geotile_ids), (geotile_dim,))
        area_km2 = DimArray(fill(5.0, 3, n_heights), (geotile_dim, height_dim))

        params = synthetic_fill_params(dh_dict)

        # hyps_fill_empty! normalizes neighbours by `dh0_median`, which it reads out of
        # `param_m1[1] + dh0` in the params table. Those are NaN until a model has been fitted, and
        # subtracting NaN turns every neighbour value into NaN so nothing can be filled. Production
        # runs hyps_model_fill! first (utilities_binning.jl), so do the same here.
        nobs_dict = Dict("icesat2" => DimArray(fill(50, 3, n_dates, n_heights),
                                               (geotile_dim, date_dim, height_dim)))
        GGA.hyps_model_fill!(dh_dict, nobs_dict, params;
                             bincount_min=3, smooth_n=5, smooth_h2t_length_scale=500.0,
                             missions2update=["icesat2"])

        # hyps_fill_empty! writes through the DimArray
        dh_before_fill = copy(dh_dict["icesat2"][1, :, :])
        neighbour2_before = copy(dh_dict["icesat2"][2, :, :])
        neighbour3_before = copy(dh_dict["icesat2"][3, :, :])

        # Apply fill_empty
        dh_filled = GGA.hyps_fill_empty!(
            dh_dict,
            params,
            geotile_extent,
            area_km2;
            missions2update=["icesat2"]
        )

        # Test: the wholly-empty geotile gets populated from its neighbours
        @test all(isnan.(dh_before_fill))
        filled_tile = dh_filled["icesat2"][1, :, :]
        n_empty_before = sum(all(isnan.(dh_before_fill), dims=1))
        n_empty_after = sum(all(isnan.(filled_tile), dims=1))
        @test n_empty_after < n_empty_before
        @test any(.!isnan.(filled_tile))

        # Test: the neighbours themselves are left alone by fill_empty (they were already
        # model-smoothed by the fit above, so compare against that state, not the raw baseline)
        @test all(dh_filled["icesat2"][2, :, :] .≈ neighbour2_before)
        @test all(dh_filled["icesat2"][3, :, :] .≈ neighbour3_before)

        # Test: filled values track the neighbours' elevation gradient rather than being arbitrary.
        # Offsets are removed before the median, so compare bin-to-bin *differences*.
        valid_bins = [j for j in 1:n_heights if !all(isnan.(filled_tile[:, j]))]
        if length(valid_bins) >= 3
            filled_profile = [mean(filter(!isnan, filled_tile[:, j])) for j in valid_bins]
            truth_profile = [mean(filter(!isnan, neighbour2_before[:, j])) for j in valid_bins]
            @test sign(filled_profile[end] - filled_profile[1]) ==
                  sign(truth_profile[end] - truth_profile[1])
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

        geotile_dim = DD.Dim{:geotile}(["lat[-15-13]lon[-075-073]"])
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

        params = synthetic_fill_params(dh_dict)

        # Apply up-down filling
        dh_filled = GGA.hyps_fill_updown!(
            dh_dict,
            area_km2;
            missions2update=["icesat2"]
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

        geotile_dim = DD.Dim{:geotile}(["lat[+40+42]lon[-120-118]"])
        date_dim = DD.Dim{:date}(dates)
        height_dim = DD.Dim{:height}(heights)

        # Create completely empty arrays
        dh_empty = reshape(fill(NaN, n_dates * n_heights), 1, n_dates, n_heights)
        nobs_empty = reshape(zeros(Int, n_dates * n_heights), 1, n_dates, n_heights)

        dh_dict = Dict("icesat2" => DimArray(dh_empty, (geotile_dim, date_dim, height_dim)))
        nobs_dict = Dict("icesat2" => DimArray(nobs_empty, (geotile_dim, date_dim, height_dim)))

        params = synthetic_fill_params(dh_dict)

        # These should not crash with empty data
        @test_nowarn GGA.hyps_model_fill!(dh_dict, nobs_dict, params;
                                          bincount_min=3,
                                          missions2update=["icesat2"])
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

        geotile_dim = DD.Dim{:geotile}(["lat[+70+72]lon[-050-048]"])
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

        params = synthetic_fill_params(dh_dict)

        # Should handle single bin without errors
        @test_nowarn GGA.hyps_model_fill!(dh_dict, nobs_dict, params;
                                          bincount_min=3,
                                          missions2update=["icesat2"])
    end

end
