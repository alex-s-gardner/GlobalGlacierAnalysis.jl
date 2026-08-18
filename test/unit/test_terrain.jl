using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA

@testset "Terrain Derivatives" begin
    @testset "slope calculation" begin
        # Geographic coordinates (EPSG 4326)
        # Test at equator (lat=0)
        dhdx_eq = 0.001  # m/m gradient in X
        dhdy_eq = 0.001  # m/m gradient in Y
        slope_eq = GGA.slope(dhdx_eq, dhdy_eq, 4326; lat=0.0)
        # Slope magnitude should be sqrt(dhdx^2 + dhdy^2)
        expected_eq = sqrt(dhdx_eq^2 + dhdy_eq^2)
        @test slope_eq ≈ expected_eq rtol=1e-3

        # Test at 45° latitude (cosine correction)
        slope_45 = GGA.slope(0.001, 0.001, 4326; lat=45.0)
        # Should differ from equator due to latitude correction
        @test slope_45 isa Number
        @test slope_45 > 0

        # Projected coordinates (UTM, EPSG 32632)
        slope_proj = GGA.slope(0.001, 0.001, 32632; lat=45.0)
        # Should be similar to geographic at same gradient
        @test slope_proj isa Number
        @test slope_proj > 0

        # Zero gradient
        slope_flat = GGA.slope(0.0, 0.0, 4326; lat=0.0)
        @test slope_flat == 0.0

        # Large gradient
        slope_steep = GGA.slope(0.5, 0.5, 4326; lat=0.0)
        expected_steep = sqrt(0.5^2 + 0.5^2)
        @test slope_steep ≈ expected_steep rtol=1e-3
    end

    @testset "curvature calculation" begin
        # Geographic coordinates
        dhddx = 0.0001  # Second derivative in X
        dhddy = 0.0001  # Second derivative in Y
        curv_eq = GGA.curvature(dhddx, dhddy, 4326; lat=0.0)
        @test curv_eq isa Number

        # Test at different latitude
        curv_45 = GGA.curvature(0.0001, 0.0001, 4326; lat=45.0)
        @test curv_45 isa Number

        # Projected coordinates
        curv_proj = GGA.curvature(0.0001, 0.0001, 32632; lat=45.0)
        @test curv_proj isa Number

        # Zero curvature (flat surface)
        curv_flat = GGA.curvature(0.0, 0.0, 4326; lat=0.0)
        @test curv_flat == 0.0

        # Negative curvature (concave)
        curv_neg = GGA.curvature(-0.0001, -0.0001, 4326; lat=0.0)
        @test curv_neg < 0
    end

    @testset "Latitude correction consistency" begin
        # Verify that latitude correction affects results consistently
        dhdx, dhdy = 0.01, 0.01

        # At equator vs 60° latitude
        slope_0 = GGA.slope(dhdx, dhdy, 4326; lat=0.0)
        slope_60 = GGA.slope(dhdx, dhdy, 4326; lat=60.0)

        # Slopes should differ due to cosine correction
        # At 60°, cos(60°) = 0.5, so correction factor is larger
        @test slope_0 != slope_60
    end
end
