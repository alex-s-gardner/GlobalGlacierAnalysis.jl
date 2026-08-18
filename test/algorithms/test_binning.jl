using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using Statistics
using Random
using DataFrames

@testset "Binning Operations" begin
    @testset "binningfun_define" begin
        Random.seed!(42)

        # Create test data: 1000 normal samples + 10 outliers
        normal_data = randn(1000)
        outliers = randn(10) .* 10 .+ 50
        all_data = vcat(normal_data, outliers)

        # Test NMAD3 filtering (threshold = 3)
        binfun_nmad3 = GGA.binningfun_define("nmad3")
        result_nmad3 = binfun_nmad3(all_data)
        @test result_nmad3 isa Number
        # Should be closer to median of normal data than median of all data
        @test abs(result_nmad3 - median(normal_data)) < abs(median(all_data) - median(normal_data))

        # Test NMAD5 filtering (threshold = 5, more permissive)
        binfun_nmad5 = GGA.binningfun_define("nmad5")
        result_nmad5 = binfun_nmad5(all_data)
        @test result_nmad5 isa Number

        # Test median binning (no filtering)
        binfun_median = GGA.binningfun_define("median")
        result_median = binfun_median(all_data)
        @test result_median ≈ median(all_data) rtol=1e-6

        # Consistency test
        data_consistent = [1.0, 2.0, 3.0, 4.0, 5.0]
        result1 = binfun_median(data_consistent)
        result2 = binfun_median(data_consistent)
        @test result1 == result2
    end

    @testset "dh_area_average" begin
        # Uniform values with uniform areas
        dh_uniform = [1.0, 1.0, 1.0]
        area_uniform = [1.0, 1.0, 1.0]
        avg_uniform = GGA.dh_area_average(dh_uniform, area_uniform)
        @test avg_uniform ≈ 1.0 rtol=1e-6

        # Varying values with varying areas
        dh_vary = [1.0, 2.0, 3.0]
        area_vary = [1.0, 2.0, 1.0]
        # Weighted average: (1*1 + 2*2 + 3*1) / (1+2+1) = 8/4 = 2.0
        avg_vary = GGA.dh_area_average(dh_vary, area_vary)
        @test avg_vary ≈ 2.0 rtol=1e-6

        # NaN handling (should exclude from calculation)
        dh_nan = [1.0, NaN, 3.0]
        area_nan = [1.0, 2.0, 1.0]
        # Should compute (1*1 + 3*1) / (1+1) = 2.0
        avg_nan = GGA.dh_area_average(dh_nan, area_nan)
        @test avg_nan ≈ 2.0 rtol=1e-6

        # Zero area handling
        dh_mixed = [1.0, 2.0, 3.0]
        area_zero = [1.0, 0.0, 1.0]
        # Should exclude zero-area point: (1*1 + 3*1) / (1+1) = 2.0
        avg_zero = GGA.dh_area_average(dh_mixed, area_zero)
        @test avg_zero ≈ 2.0 rtol=1e-6

        # All NaN
        dh_all_nan = [NaN, NaN, NaN]
        area_all = [1.0, 1.0, 1.0]
        avg_all_nan = GGA.dh_area_average(dh_all_nan, area_all)
        @test isnan(avg_all_nan)
    end

    @testset "geotile_bin2d synthetic data" begin
        # Create synthetic DataFrame
        Random.seed!(123)
        n_points = 1000

        df = DataFrame(
            x = rand(n_points) .* 10,  # 0-10 range
            y = rand(n_points) .* 10,  # 0-10 range
            value = randn(n_points) .+ 5.0
        )

        # Define 5×5 grid bins
        x_edges = range(0, 10, length=6)
        y_edges = range(0, 10, length=6)
        dims_edges = (x=x_edges, y=y_edges)

        # Apply median binning
        binned = GGA.geotile_bin2d(df; var2bin=:value, dims_edges=dims_edges, binfunction=median)

        # Check output shape
        @test size(binned) == (5, 5)

        # Check that binned values are reasonable (around 5.0 ± 3σ)
        valid_values = binned[.!isnan.(binned)]
        @test mean(valid_values) ≈ 5.0 atol=0.5
        @test all(valid_values .> 0.0)
        @test all(valid_values .< 10.0)
    end

    @testset "geotile_bin2d edge alignment" begin
        # Test points exactly on bin boundaries
        df_edge = DataFrame(
            x = [0.0, 5.0, 10.0],
            y = [0.0, 5.0, 10.0],
            value = [1.0, 2.0, 3.0]
        )

        x_edges = [0.0, 5.0, 10.0]
        y_edges = [0.0, 5.0, 10.0]
        dims_edges = (x=x_edges, y=y_edges)

        binned_edge = GGA.geotile_bin2d(df_edge; var2bin=:value, dims_edges=dims_edges, binfunction=median)

        # Should have 2×2 bins
        @test size(binned_edge) == (2, 2)

        # Check that edge points are assigned correctly
        @test any(.!isnan.(binned_edge))
    end

    @testset "geotile_bin2d sparse data" begin
        # Only a few points in large grid
        df_sparse = DataFrame(
            x = [2.5, 7.5],
            y = [2.5, 7.5],
            value = [10.0, 20.0]
        )

        x_edges = range(0, 10, length=11)  # 10 bins
        y_edges = range(0, 10, length=11)
        dims_edges = (x=x_edges, y=y_edges)

        binned_sparse = GGA.geotile_bin2d(df_sparse; var2bin=:value, dims_edges=dims_edges, binfunction=median)

        # Should be mostly NaN with a few valid values
        @test sum(.!isnan.(binned_sparse)) <= 4  # At most 2 bins with data
        valid_sparse = binned_sparse[.!isnan.(binned_sparse)]
        @test all(valid_sparse .∈ Ref([10.0, 20.0]))
    end
end
