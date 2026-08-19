using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using Statistics
using Random
using DataFrames
using DimensionalData
import DimensionalData as DD
using Dates

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
        # nmad3 takes the *mean* of the inlier subset, so the property worth asserting is
        # robustness: it lands closer to the clean mean than an unfiltered mean does. Comparing it
        # against `median(all_data)` is not a real test -- ten outliers barely move a median, so
        # that right-hand side is smaller than the filtered mean's own sampling error.
        @test abs(result_nmad3 - mean(normal_data)) < abs(mean(all_data) - mean(normal_data))
        # and the outliers near +50 are excluded outright
        @test result_nmad3 < 1.0

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

    # `dh_area_average(dh, area0)` expects `dh` to be a DimArray carrying a `:date` dimension and
    # averages over the remaining (spatial) dimension, weighting by `area0`. It returns a DimArray
    # indexed by date, not a scalar.
    @testset "dh_area_average" begin
        heights = 1000.0:100.0:1200.0          # three elevation bins
        dates = [DateTime(2019, 1, 1)]         # single date keeps the expected values obvious

        make_dh(vals) = DimArray(reshape(collect(vals), 1, 3),
                                 (Dim{:date}(dates), Dim{:height}(heights)))

        # Uniform values with uniform areas
        avg_uniform = GGA.dh_area_average(make_dh([1.0, 1.0, 1.0]), [1.0, 1.0, 1.0])
        @test avg_uniform[1] ≈ 1.0 rtol=1e-6

        # Varying values with varying areas: (1*1 + 2*2 + 3*1) / (1+2+1) = 8/4 = 2.0
        avg_vary = GGA.dh_area_average(make_dh([1.0, 2.0, 3.0]), [1.0, 2.0, 1.0])
        @test avg_vary[1] ≈ 2.0 rtol=1e-6

        # NaN values are excluded: (1*1 + 3*1) / (1+1) = 2.0
        avg_nan = GGA.dh_area_average(make_dh([1.0, NaN, 3.0]), [1.0, 2.0, 1.0])
        @test avg_nan[1] ≈ 2.0 rtol=1e-6

        # Zero-area points contribute nothing: (1*1 + 3*1) / (1+1) = 2.0
        avg_zero = GGA.dh_area_average(make_dh([1.0, 2.0, 3.0]), [1.0, 0.0, 1.0])
        @test avg_zero[1] ≈ 2.0 rtol=1e-6

        # All NaN -> NaN
        @test isnan(GGA.dh_area_average(make_dh([NaN, NaN, NaN]), [1.0, 1.0, 1.0])[1])

        # Multiple dates are averaged independently
        multi = DimArray([1.0 2.0 3.0; 4.0 4.0 4.0],
                         (Dim{:date}([DateTime(2019, 1, 1), DateTime(2019, 2, 1)]),
                          Dim{:height}(heights)))
        avg_multi = GGA.dh_area_average(multi, [1.0, 2.0, 1.0])
        @test length(avg_multi) == 2
        @test avg_multi[1] ≈ 2.0 rtol=1e-6
        @test avg_multi[2] ≈ 4.0 rtol=1e-6
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
        dims_edges = ("x" => x_edges, "y" => y_edges)

        # Apply median binning
        binned, _ = GGA.geotile_bin2d(df; var2bin="value", dims_edges=dims_edges, binfunction=median)

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
        dims_edges = ("x" => x_edges, "y" => y_edges)

        binned_edge, _ = GGA.geotile_bin2d(df_edge; var2bin="value", dims_edges=dims_edges, binfunction=median)

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
        dims_edges = ("x" => x_edges, "y" => y_edges)

        binned_sparse, _ = GGA.geotile_bin2d(df_sparse; var2bin="value", dims_edges=dims_edges, binfunction=median)

        # Should be mostly NaN with a few valid values
        @test sum(.!isnan.(binned_sparse)) <= 4  # At most 2 bins with data
        valid_sparse = binned_sparse[.!isnan.(binned_sparse)]
        @test all(valid_sparse .∈ Ref([10.0, 20.0]))
    end
end
