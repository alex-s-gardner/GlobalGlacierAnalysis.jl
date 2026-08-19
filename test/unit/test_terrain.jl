using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA

# `slope` returns the two scaled derivative *components* as a tuple -- it does not reduce them to a
# magnitude. `curvature` returns a scalar and carries a leading -2 factor, so a bowl (positive
# second derivatives) comes back negative, matching its "negative = concave" docstring.
@testset "Terrain Derivatives" begin
    @testset "slope calculation" begin
        # Geographic coordinates (EPSG 4326): inputs are per-degree and get scaled to per-metre,
        # so the components come back much smaller than the raw gradients.
        ydist_eq, xdist_eq = GGA.meters2lonlat_distance(1, 0.0)
        dhdx_scaled, dhdy_scaled = GGA.slope(0.001, 0.001, 4326; lat=0.0)
        @test dhdx_scaled ≈ 0.001 * xdist_eq rtol=1e-9
        @test dhdy_scaled ≈ 0.001 * ydist_eq rtol=1e-9
        @test dhdx_scaled > 0 && dhdy_scaled > 0

        # At 45° latitude the x component is inflated by 1/cos(lat); y is latitude-independent
        dhdx_45, dhdy_45 = GGA.slope(0.001, 0.001, 4326; lat=45.0)
        @test dhdx_45 ≈ dhdx_scaled / cosd(45.0) rtol=1e-9
        @test dhdy_45 ≈ dhdy_scaled rtol=1e-12
        @test dhdx_45 > dhdx_scaled

        # Projected coordinates (UTM) are already metres, so the inputs pass through untouched
        @test GGA.slope(0.001, 0.001, 32632; lat=45.0) == (0.001, 0.001)
        @test GGA.slope(0.5, -0.25, 32632; lat=0.0) == (0.5, -0.25)

        # Zero gradient stays zero in both components
        @test GGA.slope(0.0, 0.0, 4326; lat=0.0) == (0.0, 0.0)

        # Scaling is linear in the input gradient
        small = GGA.slope(0.001, 0.001, 4326; lat=0.0)
        large = GGA.slope(0.5, 0.5, 4326; lat=0.0)
        @test large[1] ≈ small[1] * 500 rtol=1e-9
        @test large[2] ≈ small[2] * 500 rtol=1e-9

        # Sign is preserved
        neg_x, neg_y = GGA.slope(-0.001, -0.001, 4326; lat=0.0)
        @test neg_x < 0 && neg_y < 0
    end

    @testset "curvature calculation" begin
        # Geographic coordinates
        @test GGA.curvature(0.0001, 0.0001, 4326; lat=0.0) isa Number
        @test GGA.curvature(0.0001, 0.0001, 4326; lat=45.0) isa Number

        # Projected coordinates: xdist = ydist = 1, so the result is just -2*(sum)*100
        @test GGA.curvature(0.0001, 0.0001, 32632; lat=45.0) ≈ -2 * 0.0002 * 100 rtol=1e-12

        # Zero curvature (flat surface)
        @test GGA.curvature(0.0, 0.0, 4326; lat=0.0) == 0.0

        # Sign convention: positive second derivatives (a bowl) are concave -> negative curvature
        @test GGA.curvature(0.0001, 0.0001, 4326; lat=0.0) < 0

        # negative second derivatives (a dome) are convex -> positive curvature
        @test GGA.curvature(-0.0001, -0.0001, 4326; lat=0.0) > 0

        # and the two are symmetric about zero
        @test GGA.curvature(0.0001, 0.0001, 4326; lat=0.0) ≈
              -GGA.curvature(-0.0001, -0.0001, 4326; lat=0.0) rtol=1e-12
    end

    @testset "Latitude correction consistency" begin
        # Verify that latitude correction affects results consistently
        dhdx, dhdy = 0.01, 0.01

        slope_0 = GGA.slope(dhdx, dhdy, 4326; lat=0.0)
        slope_60 = GGA.slope(dhdx, dhdy, 4326; lat=60.0)

        # Slopes differ because of the 1/cos(lat) longitude correction
        @test slope_0 != slope_60

        # cos(60°) = 0.5, so the x component roughly doubles while y is unchanged
        @test slope_60[1] ≈ slope_0[1] / cosd(60.0) rtol=1e-9
        @test slope_60[1] ≈ 2 * slope_0[1] rtol=1e-9
        @test slope_60[2] ≈ slope_0[2] rtol=1e-12

        # Projected coordinates carry no latitude dependence at all
        @test GGA.slope(dhdx, dhdy, 32632; lat=0.0) == GGA.slope(dhdx, dhdy, 32632; lat=60.0)
    end
end
