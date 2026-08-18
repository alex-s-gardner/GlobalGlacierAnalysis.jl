using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using Statistics
using Random

@testset "Statistical Functions" begin
    @testset "nmad - normalized median absolute deviation" begin
        # Uniform data (degenerate case)
        uniform = [0.0, 0.0, 0.0, 0.0]
        nmad_result = GGA.nmad(uniform)
        # All zeros since median is 0 and all deviations are 0
        @test all(isnan.(nmad_result) .| (nmad_result .== 0.0))

        # Symmetric distribution
        symmetric = [-2.0, -1.0, 0.0, 1.0, 2.0]
        nmad_sym = GGA.nmad(symmetric)
        # Median is 0, MAD is 1, so nmad should be deviations / (1.0 * 1.4826)
        expected = abs.(symmetric) ./ (1.0 * 1.4826)
        @test nmad_sym ≈ expected rtol=1e-6

        # Data with outliers
        with_outlier = [1.0, 2.0, 3.0, 4.0, 5.0, 100.0]
        nmad_outlier = GGA.nmad(with_outlier)
        # Outlier should have large NMAD value
        @test maximum(nmad_outlier) > 10.0

        # Normal distribution test (statistical property)
        Random.seed!(42)
        normal_data = randn(1000)
        nmad_normal = GGA.nmad(normal_data)
        # For N(0,1), approximately 68% should be within 1 NMAD
        within_1nmad = count(nmad_normal .< 1.0) / length(nmad_normal)
        @test within_1nmad ≈ 0.68 rtol=0.1
    end

    @testset "mad - median absolute deviation" begin
        # Test with simple sequence
        data = [1, 2, 3, 4, 5]
        # Median is 3, absolute deviations are [2, 1, 0, 1, 2], MAD is median of that = 1
        @test GGA.mad(data) == 1.0

        # Test empty array
        @test isnan(GGA.mad([]))

        # Test single value
        @test GGA.mad([42.0]) == 0.0

        # Test two values
        @test GGA.mad([1.0, 3.0]) == 1.0  # Median 2, deviations [1, 1], MAD = 1

        # Test with negative numbers
        neg_data = [-5, -3, -1, 1, 3]
        # Median is -1, deviations are [4, 2, 0, 2, 4], MAD is median = 2
        @test GGA.mad(neg_data) == 2.0

        # Test with duplicates
        dups = [1, 1, 1, 2, 2, 2, 3, 3, 3]
        # Median is 2, deviations are [1,1,1,0,0,0,1,1,1], MAD = 1
        @test GGA.mad(dups) == 1.0
    end

    @testset "Statistical consistency" begin
        # Test that NMAD = (x - median(x)) / (MAD * 1.4826)
        Random.seed!(123)
        data = randn(100) .* 10 .+ 50

        mad_val = GGA.mad(data)
        nmad_vals = GGA.nmad(data)

        # Manually calculate expected NMAD
        consistent_estimator = 1.4826
        deviations = abs.(data .- median(data))
        expected_nmad = deviations ./ (mad_val * consistent_estimator)

        @test nmad_vals ≈ expected_nmad rtol=1e-6
    end
end
