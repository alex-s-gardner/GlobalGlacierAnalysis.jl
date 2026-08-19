"""
Tests for GEMB (Glacier Energy and Mass Balance) model integration.

These tests validate GEMB file reading, physical constraint enforcement,
model calibration, and volume change calculations.
"""

using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using Dates
using Random
using Statistics

# Include test fixtures
include("../fixtures/mock_gemb_files.jl")

@testset "GEMB Operations" begin
    # `gemb_rate_physical_constraints!` is keyed by *String* and requires "rain" alongside
    # acc/melt/refreeze/ec. It mutates the vectors in place and also writes "runoff" and "smb".
    @testset "gemb_rate_physical_constraints!" begin
        gemb = Dict(
            "acc" => [-1.0, 2.0],       # Negative accumulation (invalid)
            "rain" => [0.2, -0.3],      # Negative rain (invalid)
            "melt" => [3.0, 4.0],
            "refreeze" => [5.0, 2.0],   # Refreeze > melt at index 1 (invalid)
            "ec" => [0.5, -0.1],
        )

        GGA.gemb_rate_physical_constraints!(gemb)

        # Check corrections
        @test all(gemb["acc"] .>= 0)                        # No negative accumulation
        @test all(gemb["rain"] .>= 0)                       # No negative rain
        @test all(gemb["melt"] .>= 0)                       # No negative melt
        @test all(gemb["refreeze"] .>= 0)                   # No negative refreeze
        @test all(gemb["refreeze"] .<= gemb["melt"])        # Refreeze cannot exceed melt

        # Derived quantities are written back, and exclude rain by design
        @test gemb["runoff"] == gemb["melt"] .- gemb["refreeze"]
        @test gemb["smb"] == gemb["acc"] .- gemb["runoff"] .- gemb["ec"]

        # spot-check the clamped first element: acc 0, refreeze clamped to melt -> runoff 0
        @test gemb["acc"][1] == 0.0
        @test gemb["refreeze"][1] == 3.0
        @test gemb["runoff"][1] == 0.0
    end

    # `gemb_add_derived_vars!` is Symbol-keyed and built on `merge`, so it needs a NamedTuple
    # (not a Dict) and -- despite the `!` -- returns a new object rather than mutating.
    @testset "gemb_add_derived_vars!" begin
        dv_gemb = (; acc=10.0, melt=8.0, refreeze=3.0, ec=1.0, fac=50.0)

        result = GGA.gemb_add_derived_vars!(dv_gemb)

        # smb = acc - melt + refreeze - ec
        @test result[:smb] ≈ 10.0 - 8.0 + 3.0 - 1.0 rtol=1e-6

        # runoff = melt - refreeze
        @test result[:runoff] ≈ 8.0 - 3.0 rtol=1e-6

        # dv = smb + fac
        @test result[:dv] ≈ result[:smb] + 50.0 rtol=1e-6

        # the input is left untouched -- the function is not actually in-place
        @test !haskey(dv_gemb, :smb)

        # works elementwise on vectors too
        vec_gemb = (; acc=[10.0, 4.0], melt=[8.0, 1.0], refreeze=[3.0, 0.5],
                      ec=[1.0, 0.25], fac=[50.0, 20.0])
        vres = GGA.gemb_add_derived_vars!(vec_gemb)
        @test vres[:runoff] ≈ [5.0, 0.5]
        @test vres[:smb] ≈ [4.0, 3.25]
    end

    @testset "GEMB file reading - gemb_read2" begin
        mktempdir() do temp_dir
            # Create mock GEMB file
            gemb_file = joinpath(temp_dir, "test_gemb.mat")
            create_mock_gemb_mat(gemb_file; n_points=5, n_times=12, pscale=1.0)

            # Read the file
            gemb = GGA.gemb_read2(gemb_file)

            # Test: All expected variables should be present
            @test haskey(gemb, "latitude")
            @test haskey(gemb, "longitude")
            @test haskey(gemb, "date")
            @test haskey(gemb, "smb")
            @test haskey(gemb, "fac")
            @test haskey(gemb, "acc")
            @test haskey(gemb, "runoff")
            @test haskey(gemb, "melt")
            @test haskey(gemb, "refreeze")

            # Test: Dimensions should be correct
            n_points = length(gemb["latitude"])
            n_times = length(gemb["date"])
            @test n_points == 5
            @test n_times == 12

            # Test: Data arrays should have correct shape
            @test size(gemb["smb"]) == (n_points, n_times)
            @test size(gemb["acc"]) == (n_points, n_times)
            @test size(gemb["runoff"]) == (n_points, n_times)

            # Test: Coordinates should be reasonable
            @test all(60.0 .<= gemb["latitude"] .<= 62.0)
            @test all(-52.0 .<= gemb["longitude"] .<= -48.0)

            # Test: Dates should be in decimal years and increasing
            @test all(diff(gemb["date"]) .> 0)  # Monotonically increasing
            @test gemb["date"][1] >= 2015.0
            @test gemb["date"][end] <= 2020.0

            # Test: Variables should be rates (not cumulative after conversion)
            # The function converts cumulative to rates
            @test !all(gemb["smb"] .== 0)  # Should have non-zero values
        end
    end

    @testset "GEMB ensemble with varying pscale" begin
        mktempdir() do temp_dir
            # Create ensemble with different precipitation scaling
            n_members = 5
            pscale_values = range(0.8, 1.2, length=n_members)
            gemb_files = []

            for (i, pscale) in enumerate(pscale_values)
                filename = joinpath(temp_dir, "gemb_pscale_$(i).mat")
                create_mock_gemb_mat(filename; n_points=10, n_times=24, pscale=pscale)
                push!(gemb_files, filename)
            end

            # Read all ensemble members
            gemb_ensemble = [GGA.gemb_read2(f) for f in gemb_files]

            # Test: SMB should scale with pscale
            # Higher pscale should lead to less negative SMB (more accumulation)
            smb_means = [mean(g["smb"]) for g in gemb_ensemble]

            # With higher pscale, more snow accumulates, so SMB becomes less negative
            @test smb_means[end] > smb_means[1]  # pscale=1.2 > pscale=0.8
            @test issorted(smb_means)            # and it scales monotonically

            # Accumulation is what pscale actually scales
            acc_means = [mean(g["acc"]) for g in gemb_ensemble]
            @test issorted(acc_means)
            @test acc_means[end] > acc_means[1]

            # Runoff, by contrast, is melt-driven: melt is energy-limited rather than
            # precipitation-limited, so scaling precipitation leaves runoff essentially unchanged.
            runoff_means = [mean(g["runoff"]) for g in gemb_ensemble]
            @test runoff_means[end] ≈ runoff_means[1] rtol=0.05
        end
    end

    @testset "GEMB physical constraints enforcement" begin
        # String-keyed, and "rain" is required alongside the rest
        gemb_extreme = Dict(
            "acc" => [5.0, -10.0, 3.0],      # Large negative accumulation
            "rain" => [0.0, -1.0, 0.5],
            "melt" => [2.0, 1.0, 4.0],
            "refreeze" => [3.0, 5.0, 2.0],   # Refreeze > melt in multiple places
            "ec" => [1.0, -5.0, 0.5],        # Large negative elevation change
        )

        # Apply constraints
        GGA.gemb_rate_physical_constraints!(gemb_extreme)

        # Test: All values should be physically reasonable. Note `ec` is deliberately *not*
        # clamped -- sublimation/condensation is signed, and it enters smb as a subtraction.
        @test all(gemb_extreme["acc"] .>= 0)
        @test all(gemb_extreme["rain"] .>= 0)
        @test all(gemb_extreme["melt"] .>= 0)
        @test all(gemb_extreme["refreeze"] .>= 0)

        # Test: Mass balance constraints
        # Refreeze cannot exceed melt (can't refreeze more than what melted)
        for i in eachindex(gemb_extreme["melt"])
            @test gemb_extreme["refreeze"][i] <= gemb_extreme["melt"][i]
        end

        # Runoff is non-negative once refreeze is capped at melt
        @test all(gemb_extreme["runoff"] .>= 0)
    end

    @testset "GEMB derived variable calculations" begin
        # Symbol-keyed NamedTuple, and the result must be captured (merge, not mutation)
        dv = (;
            acc=2.0,       # 2 m/yr accumulation
            melt=1.5,      # 1.5 m/yr melt
            refreeze=0.3,  # 0.3 m/yr refreezes
            ec=0.1,        # 0.1 m/yr elevation change
            fac=10.0,      # 10 m firn air content
        )

        dv = GGA.gemb_add_derived_vars!(dv)

        # Test: SMB calculation
        # SMB = accumulation - melt + refreeze - elevation_change
        @test dv[:smb] ≈ 2.0 - 1.5 + 0.3 - 0.1 rtol=1e-10

        # Test: Runoff calculation
        # Runoff = melt - refreeze (water that leaves the system)
        @test dv[:runoff] ≈ 1.5 - 0.3 rtol=1e-10

        # Test: Runoff should always be >= 0 (after refreeze is capped)
        @test dv[:runoff] >= 0

        # dv = smb + fac
        @test dv[:dv] ≈ dv[:smb] + 10.0 rtol=1e-10
    end

    @testset "GEMB with extreme precipitation scaling" begin
        mktempdir() do temp_dir
            # Test very low and very high pscale values
            pscale_extremes = create_mock_gemb_with_extremes(temp_dir)

            gemb_low = GGA.gemb_read2(pscale_extremes[:low])
            gemb_med = GGA.gemb_read2(pscale_extremes[:medium])
            gemb_high = GGA.gemb_read2(pscale_extremes[:high])

            # Test: Accumulation should scale with pscale
            acc_low = mean(gemb_low["acc"])
            acc_med = mean(gemb_med["acc"])
            acc_high = mean(gemb_high["acc"])

            @test acc_low < acc_med < acc_high

            # Test: Higher pscale leads to more positive SMB
            smb_low = mean(gemb_low["smb"])
            smb_high = mean(gemb_high["smb"])
            @test smb_high > smb_low

            # Runoff is melt-driven and melt is energy-limited, so scaling precipitation leaves
            # runoff essentially unchanged rather than increasing it.
            runoff_low = mean(gemb_low["runoff"])
            runoff_high = mean(gemb_high["runoff"])
            @test runoff_high ≈ runoff_low rtol=0.05
        end
    end

    @testset "GEMB calibration - parameter recovery" begin
        # Simulate a simple calibration scenario
        # We know the "true" pscale and test if we can identify it

        mktempdir() do temp_dir
            # Create ensemble with known true parameter
            true_pscale = 1.1
            n_ensemble = 5
            pscale_range = range(0.8, 1.4, length=n_ensemble)

            # Generate ensemble
            gemb_files = []
            for pscale in pscale_range
                filename = joinpath(temp_dir, "gemb_pscale_$(pscale).mat")
                create_mock_gemb_mat(filename; n_points=15, n_times=36, pscale=pscale)
                push!(gemb_files, filename)
            end

            # Read ensemble
            gemb_ensemble = [GGA.gemb_read2(f) for f in gemb_files]

            # Generate synthetic "observed" SMB from the true parameter
            Random.seed!(42)
            true_file = joinpath(temp_dir, "gemb_true.mat")
            create_mock_gemb_mat(true_file; n_points=15, n_times=36, pscale=true_pscale)
            gemb_true = GGA.gemb_read2(true_file)

            # Add observation noise
            obs_smb = gemb_true["smb"] .+ 0.2 .* randn(size(gemb_true["smb"])...)

            # Compute RMSE for each ensemble member against observations
            rmse_values = Float64[]
            for gemb in gemb_ensemble
                rmse = sqrt(mean((gemb["smb"] .- obs_smb).^2))
                push!(rmse_values, rmse)
            end

            # Test: Best ensemble member should be closest to true pscale
            best_idx = argmin(rmse_values)
            best_pscale = pscale_range[best_idx]

            # Should identify the closest pscale value
            @test abs(best_pscale - true_pscale) <= 0.3  # Within one ensemble spacing

            # Test: RMSE should increase as we move away from true parameter
            # (at least for nearest neighbors)
            if best_idx > 1
                @test rmse_values[best_idx] < rmse_values[best_idx - 1]
            end
            if best_idx < n_ensemble
                @test rmse_values[best_idx] < rmse_values[best_idx + 1]
            end
        end
    end

    @testset "GEMB temporal rate conversion" begin
        # Test that cumulative variables are properly converted to rates
        mktempdir() do temp_dir
            gemb_file = joinpath(temp_dir, "test_rates.mat")
            create_mock_gemb_mat(gemb_file; n_points=3, n_times=24, pscale=1.0)

            gemb = GGA.gemb_read2(gemb_file)

            # Test: Rates should have reasonable magnitudes
            # SMB typically -2 to +2 m/yr for glaciers
            @test all(-5.0 .<= gemb["smb"] .<= 5.0)

            # Test: Accumulation should be positive (after constraint enforcement)
            @test all(gemb["acc"] .>= 0)

            # Test: Melt should be positive
            @test all(gemb["melt"] .>= 0)

            # Test: Runoff = melt - refreeze (should be consistent)
            for i in 1:size(gemb["runoff"], 1)
                for j in 1:size(gemb["runoff"], 2)
                    expected_runoff = gemb["melt"][i, j] - gemb["refreeze"][i, j]
                    @test gemb["runoff"][i, j] ≈ expected_runoff rtol=1e-6
                end
            end
        end
    end

    @testset "Edge case - zero melt (cold glacier)" begin
        # Test GEMB behavior for a very cold glacier with no melt
        gemb_cold = (;
            acc=[2.0, 2.5, 1.8],
            melt=[0.0, 0.0, 0.0],       # No melt
            refreeze=[0.0, 0.0, 0.0],   # No refreeze possible
            ec=[0.1, 0.05, 0.08],
            fac=[5.0, 5.2, 4.8],
        )

        gemb_cold = GGA.gemb_add_derived_vars!(gemb_cold)

        # Test: SMB should equal accumulation - ec (no melt/refreeze)
        @test gemb_cold[:smb] ≈ gemb_cold[:acc] .- gemb_cold[:ec] rtol=1e-10

        # Test: Runoff should be zero everywhere (no melt)
        @test all(gemb_cold[:runoff] .== 0.0)
    end

    @testset "Edge case - complete melt (warm glacier)" begin
        # Test GEMB for a glacier where all accumulation melts
        gemb_warm = (;
            acc=[1.0, 1.2, 0.9],
            melt=[1.5, 1.8, 1.4],        # More melt than accumulation
            refreeze=[0.1, 0.15, 0.08],
            ec=[-0.3, -0.4, -0.35],      # Surface lowering
            fac=[2.0, 1.8, 1.9],
        )

        gemb_warm = GGA.gemb_add_derived_vars!(gemb_warm)

        # Test: SMB should be negative everywhere (losing mass)
        @test all(gemb_warm[:smb] .< 0)

        # Test: Runoff should be substantial
        @test all(gemb_warm[:runoff] .> 1.0)  # Most meltwater runs off

        # Test: Mass balance: acc - melt + refreeze - ec
        expected_smb = gemb_warm[:acc] .- gemb_warm[:melt] .+ gemb_warm[:refreeze] .- gemb_warm[:ec]
        @test gemb_warm[:smb] ≈ expected_smb rtol=1e-10
    end
end
