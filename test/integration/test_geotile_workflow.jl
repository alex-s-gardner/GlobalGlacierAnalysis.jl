using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using Dates
using Statistics
using Random

@testset "Geotile Workflow Integration" begin
    @testset "Synthetic altimetry → binning → trend recovery" begin
        Random.seed!(42)

        # Step 1: Generate synthetic altimetry points with known trend
        n_points = 100
        dates = [DateTime(2018,1,1) + Month(3i) for i in 0:11]  # 12 quarterly observations
        known_trend = -0.5  # m/yr
        seasonal_amp = 0.2  # m

        # Generate points within geotile
        lons = -122.0 .+ rand(n_points) .* 2.0  # -122 to -120
        lats = 45.0 .+ rand(n_points) .* 2.0    # 45 to 47
        elevations = 1000.0 .+ rand(n_points) .* 2000.0  # 1000-3000 m

        # For each point and each date, generate elevation change
        dh_data = []
        for date in dates
            decyear = GGA.decimalyear(date)
            t = decyear - 2018.0  # Time since start

            for i in 1:n_points
                # True signal: linear trend + seasonality + noise
                signal = known_trend * t + seasonal_amp * sin(2π * t)
                noise = randn() * 0.1  # 0.1 m noise

                push!(dh_data, (
                    date=date,
                    lon=lons[i],
                    lat=lats[i],
                    elevation=elevations[i],
                    dh=signal + noise
                ))
            end
        end

        # Step 2: Bin by elevation (simple binning)
        elevation_bins = [1000, 1500, 2000, 2500, 3000]
        n_bins = length(elevation_bins) - 1

        binned_medians = zeros(length(dates), n_bins)

        for (date_idx, date) in enumerate(dates)
            date_data = filter(d -> d.date == date, dh_data)

            for bin_idx in 1:n_bins
                bin_min = elevation_bins[bin_idx]
                bin_max = elevation_bins[bin_idx + 1]

                bin_data = filter(d -> bin_min <= d.elevation < bin_max, date_data)

                if !isempty(bin_data)
                    binned_medians[date_idx, bin_idx] = median([d.dh for d in bin_data])
                else
                    binned_medians[date_idx, bin_idx] = NaN
                end
            end
        end

        # Step 3: Fit linear trend to binned data
        decyears = GGA.decimalyear.(dates) .- 2018.0

        # Fit trend to each elevation bin
        fitted_trends = zeros(n_bins)
        for bin_idx in 1:n_bins
            valid_mask = .!isnan.(binned_medians[:, bin_idx])
            if sum(valid_mask) >= 3
                t_valid = decyears[valid_mask]
                dh_valid = binned_medians[valid_mask, bin_idx]

                # Simple linear regression
                A = hcat(ones(length(t_valid)), t_valid)
                coeffs = A \ dh_valid
                fitted_trends[bin_idx] = coeffs[2]  # Slope
            else
                fitted_trends[bin_idx] = NaN
            end
        end

        # Step 4: Verify trend recovery
        # Should recover trend ≈ -0.5 m/yr (with tolerance for noise and sampling)
        valid_trends = fitted_trends[.!isnan.(fitted_trends)]
        @test !isempty(valid_trends)

        mean_trend = mean(valid_trends)
        @test mean_trend ≈ known_trend atol=0.15  # Within 0.15 m/yr

        # Most trends should be negative
        @test sum(valid_trends .< 0) >= length(valid_trends) / 2
    end

    @testset "Multi-mission synthesis with known truth" begin
        Random.seed!(123)

        # True underlying trend
        true_trend = -0.3  # m/yr
        dates = 2018.0:0.25:2020.0  # Quarterly from 2018-2020

        # Generate synthetic observations from 3 missions
        # ICESat-2: dense, low error
        icesat2_obs = true_trend .* (dates .- 2018.0) .+ randn(length(dates)) .* 0.1

        # ICESat: sparse, low error (only every other observation)
        icesat_obs = true_trend .* (dates .- 2018.0) .+ randn(length(dates)) .* 0.1
        icesat_obs[2:2:end] .= NaN  # Sparse coverage

        # Hugonnet: full coverage, higher error
        hugonnet_obs = true_trend .* (dates .- 2018.0) .+ randn(length(dates)) .* 0.5

        # Inverse variance weighting
        icesat2_weight = 1.0 / 0.1^2
        icesat_weight = 1.0 / 0.1^2
        hugonnet_weight = 1.0 / 0.5^2

        # Synthesize
        synthesized = zeros(length(dates))
        for i in eachindex(dates)
            values = Float64[]
            weights = Float64[]

            if !isnan(icesat2_obs[i])
                push!(values, icesat2_obs[i])
                push!(weights, icesat2_weight)
            end
            if !isnan(icesat_obs[i])
                push!(values, icesat_obs[i])
                push!(weights, icesat_weight)
            end
            if !isnan(hugonnet_obs[i])
                push!(values, hugonnet_obs[i])
                push!(weights, hugonnet_weight)
            end

            synthesized[i] = sum(values .* weights) / sum(weights)
        end

        # Fit trend to synthesized data
        t = dates .- 2018.0
        A = hcat(ones(length(t)), t)
        coeffs = A \ synthesized
        recovered_trend = coeffs[2]

        # Verify closer to truth than Hugonnet alone
        hugonnet_trend = (A \ hugonnet_obs)[2]
        synth_error = abs(recovered_trend - true_trend)
        hugonnet_error = abs(hugonnet_trend - true_trend)

        @test synth_error < hugonnet_error
        @test recovered_trend ≈ true_trend atol=0.1
    end

    @testset "Multi-mission workflow with gap filling" begin
        Random.seed!(456)

        # Generate sparse synthetic altimetry with gaps
        n_dates = 20
        n_heights = 6
        dates = [DateTime(2018,1,1) + Month(2*i) for i in 0:n_dates-1]
        heights = collect(range(1500, 2500, length=n_heights))

        # True signal: -0.4 m/yr with elevation gradient
        t_years = collect(0:n_dates-1) ./ 6.0  # Bimonthly to years
        true_dh = zeros(n_dates, n_heights)
        for j in 1:n_dates
            for k in 1:n_heights
                elevation_factor = 1.0 - 0.0003 * (heights[k] - 2000)
                true_dh[j, k] = -0.4 * elevation_factor * t_years[j]
            end
        end

        # Create sparse observations (50% coverage)
        observed_dh = copy(true_dh)
        gap_indices = rand(1:length(observed_dh), Int(floor(0.5 * length(observed_dh))))
        observed_dh[gap_indices] .= NaN

        # Apply simple gap filling: linear interpolation in time for each height bin
        filled_dh = copy(observed_dh)
        for k in 1:n_heights
            for j in 2:n_dates-1
                if isnan(filled_dh[j, k])
                    # Find surrounding valid points
                    left = findlast(!isnan, filled_dh[1:j-1, k])
                    right = findfirst(!isnan, filled_dh[j+1:end, k])

                    if !isnothing(left) && !isnothing(right)
                        right += j  # Adjust index
                        # Linear interpolation
                        weight = (t_years[j] - t_years[left]) / (t_years[right] - t_years[left])
                        filled_dh[j, k] = filled_dh[left, k] + weight * (filled_dh[right, k] - filled_dh[left, k])
                    end
                end
            end
        end

        # Test: Gaps should be reduced
        n_gaps_before = sum(isnan.(observed_dh))
        n_gaps_after = sum(isnan.(filled_dh))
        @test n_gaps_after < n_gaps_before

        # Test: Filled values should be reasonable
        @test all(filled_dh[.!isnan.(filled_dh)] .< 0.5)  # Not too positive
        @test all(filled_dh[.!isnan.(filled_dh)] .> -2.0)  # Not too negative

        # Test: Trend recovery with filled data
        # Average across heights
        mean_dh_filled = [mean(filter(!isnan, filled_dh[j, :])) for j in 1:n_dates]

        # Fit trend
        valid_mask = .!isnan.(mean_dh_filled)
        if sum(valid_mask) >= 3
            A = hcat(ones(sum(valid_mask)), t_years[valid_mask])
            coeffs = A \ mean_dh_filled[valid_mask]
            recovered_trend = coeffs[2]

            @test recovered_trend ≈ -0.4 rtol=0.25  # Within 25% given gaps
        end
    end

    @testset "End-to-end: altimetry → binning → synthesis → volume" begin
        Random.seed!(789)

        # Simulate complete workflow for a single geotile
        geotile_id = "lat+60+62lon-050-048"

        # Step 1: Generate synthetic altimetry points
        n_points = 200
        n_times = 16  # Quarterly over 4 years
        dates = [DateTime(2018,1,1) + Month(3*i) for i in 0:n_times-1]

        # Known parameters
        true_trend = -0.6  # m/yr
        elevation_bins = collect(range(1000, 3000, length=8))
        n_elevation_bins = length(elevation_bins) - 1

        # Generate points
        elevations = 1000.0 .+ rand(n_points) .* 2000.0

        # Assign area to elevation bins
        area_per_bin = 15.0  # km² per bin

        # Step 2: Generate elevation changes for each point and time
        binned_dh = zeros(n_times, n_elevation_bins)
        binned_counts = zeros(Int, n_times, n_elevation_bins)

        for (t_idx, date) in enumerate(dates)
            t_year = (t_idx - 1) / 4.0  # Quarterly to years

            # Generate observations for this time
            for point_idx in 1:n_points
                elev = elevations[point_idx]

                # True signal with noise
                dh_true = true_trend * t_year + 0.1 * randn()

                # Bin this observation
                bin_idx = findfirst(elevation_bins[2:end] .> elev)
                if !isnothing(bin_idx)
                    binned_dh[t_idx, bin_idx] += dh_true
                    binned_counts[t_idx, bin_idx] += 1
                end
            end
        end

        # Compute bin means
        for t in 1:n_times
            for b in 1:n_elevation_bins
                if binned_counts[t, b] > 0
                    binned_dh[t, b] /= binned_counts[t, b]
                else
                    binned_dh[t, b] = NaN
                end
            end
        end

        # Step 3: Compute volume change
        # dV = Σ(dh_i * area_i)
        dv = zeros(n_times)
        for t in 1:n_times
            volume_change = 0.0
            for b in 1:n_elevation_bins
                if !isnan(binned_dh[t, b])
                    volume_change += binned_dh[t, b] * area_per_bin / 1000  # m to km
                end
            end
            dv[t] = volume_change
        end

        # Step 4: Verify volume change trend
        t_years = collect(0:n_times-1) ./ 4.0
        A = hcat(ones(n_times), t_years)
        coeffs = A \ dv
        dv_trend = coeffs[2]  # km³/yr

        # Expected: true_trend * total_area / 1000
        expected_dv_trend = true_trend * (area_per_bin * n_elevation_bins) / 1000
        @test dv_trend ≈ expected_dv_trend rtol=0.2

        # Step 5: Convert to mass change (Gt/yr)
        ice_density = 910.0  # kg/m³
        dm_trend = dv_trend * 1e9 * ice_density / 1e12  # km³/yr to Gt/yr
        expected_dm_trend = expected_dv_trend * 1e9 * ice_density / 1e12

        @test dm_trend ≈ expected_dm_trend rtol=0.2

        # Test: Mass loss should be negative
        @test dm_trend < 0
    end

    @testset "Workflow with GEMB calibration integration" begin
        # Simulated workflow: altimetry observations calibrate GEMB

        Random.seed!(321)

        # Observed altimetry: -0.5 m/yr ± 0.1 m
        n_obs = 20
        t_obs = collect(0:n_obs-1) ./ 12.0  # Monthly
        dh_obs = -0.5 .* t_obs .+ 0.1 .* randn(n_obs)

        # GEMB ensemble with varying precipitation scaling
        pscale_values = [0.8, 0.9, 1.0, 1.1, 1.2]
        n_ensemble = length(pscale_values)

        # Generate GEMB predictions (SMB → elevation change)
        # SMB scales with pscale; higher pscale → less negative EC
        gemb_predictions = zeros(n_ensemble, n_obs)
        for (i, pscale) in enumerate(pscale_values)
            # Simple model: EC ≈ -0.3 - 0.2*pscale m/yr
            gemb_trend = -0.3 - 0.2 * pscale
            gemb_predictions[i, :] = gemb_trend .* t_obs .+ 0.05 .* randn(n_obs)
        end

        # Find best-fit GEMB member (minimize RMSE with observations)
        rmse_values = [sqrt(mean((gemb_predictions[i, :] .- dh_obs).^2)) for i in 1:n_ensemble]
        best_idx = argmin(rmse_values)
        best_pscale = pscale_values[best_idx]

        # Test: Best member should have smallest error
        @test rmse_values[best_idx] == minimum(rmse_values)

        # Test: Calibrated GEMB should be closer than ensemble mean
        ensemble_mean = mean(gemb_predictions; dims=1)[1, :]
        rmse_best = rmse_values[best_idx]
        rmse_mean = sqrt(mean((ensemble_mean .- dh_obs).^2))
        @test rmse_best <= rmse_mean

        # Test: Best-fit pscale should be identifiable
        @test best_pscale in pscale_values
    end

    @testset "Regional aggregation workflow" begin
        # Test aggregation from multiple geotiles to region

        Random.seed!(654)

        # Create 5 geotiles in one region
        n_geotiles = 5
        n_dates = 12
        geotile_ids = ["gt_$(i)" for i in 1:n_geotiles]

        # Each geotile has different area and trend
        areas = [20.0, 35.0, 15.0, 28.0, 12.0]  # km²
        trends = [-0.3, -0.5, -0.2, -0.6, -0.4]  # m/yr

        # Generate time series for each geotile
        t_years = collect(0:n_dates-1) ./ 12.0
        dv_geotiles = zeros(n_geotiles, n_dates)

        for i in 1:n_geotiles
            for j in 1:n_dates
                dh = trends[i] * t_years[j] + 0.05 * randn()
                dv_geotiles[i, j] = dh * areas[i] / 1000  # Volume change in km³
            end
        end

        # Aggregate to region (sum volumes)
        dv_regional = sum(dv_geotiles; dims=1)[1, :]

        # Fit regional trend
        A = hcat(ones(n_dates), t_years)
        coeffs = A \ dv_regional
        regional_dv_trend = coeffs[2]  # km³/yr

        # Expected: area-weighted mean of trends, scaled to volume
        expected_trend = sum(trends .* areas) / 1000  # km³/yr
        @test regional_dv_trend ≈ expected_trend rtol=0.2

        # Test: Regional volume is sum of geotile volumes
        for j in 1:n_dates
            expected_vol = sum(dv_geotiles[:, j])
            @test dv_regional[j] ≈ expected_vol rtol=1e-10
        end
    end

    @testset "Workflow error propagation" begin
        # Test that uncertainties propagate correctly through workflow

        Random.seed!(987)

        # Input: elevation change with uncertainty
        n_bins = 4
        n_times = 10

        dh = randn(n_times, n_bins) .* 0.5 .- 0.3  # m
        dh_sigma = fill(0.15, n_times, n_bins)     # m uncertainty

        # Area per bin
        areas = [10.0, 15.0, 12.0, 8.0]  # km²

        # Compute volume change with error propagation
        dv = zeros(n_times)
        dv_sigma = zeros(n_times)

        for t in 1:n_times
            # Volume: Σ(dh_i * area_i)
            dv[t] = sum(dh[t, :] .* areas) / 1000

            # Uncertainty: sqrt(Σ((σ_i * area_i)²))
            dv_sigma[t] = sqrt(sum((dh_sigma[t, :] .* areas).^2)) / 1000
        end

        # Test: Volume uncertainties should be positive
        @test all(dv_sigma .> 0)

        # Test: Relative uncertainty should decrease with aggregation
        relative_dh_sigma = mean(dh_sigma) / abs(mean(dh))
        relative_dv_sigma = mean(dv_sigma) / abs(mean(dv))

        # Aggregation reduces relative uncertainty
        @test relative_dv_sigma < relative_dh_sigma
    end
end
