using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using Extents

@testset "Extent Conversions" begin
    @testset "nt2extent basic conversion" begin
        # Standard NamedTuple to Extent
        nt = (min_x=0.0, min_y=0.0, max_x=10.0, max_y=10.0)
        ext = GGA.nt2extent(nt)
        @test ext.X == (0.0, 10.0)
        @test ext.Y == (0.0, 10.0)

        # Negative coordinates
        nt_neg = (min_x=-120.0, min_y=-45.0, max_x=-110.0, max_y=-35.0)
        ext_neg = GGA.nt2extent(nt_neg)
        @test ext_neg.X == (-120.0, -110.0)
        @test ext_neg.Y == (-45.0, -35.0)
    end

    @testset "extent2nt basic conversion" begin
        # Standard Extent to NamedTuple
        ext = Extent(X=(0.0, 10.0), Y=(0.0, 10.0))
        nt = GGA.extent2nt(ext)
        @test nt.min_x == 0.0
        @test nt.max_x == 10.0
        @test nt.min_y == 0.0
        @test nt.max_y == 10.0

        # Geographic coordinates
        ext_geo = Extent(X=(-180.0, 180.0), Y=(-90.0, 90.0))
        nt_geo = GGA.extent2nt(ext_geo)
        @test nt_geo.min_x == -180.0
        @test nt_geo.max_x == 180.0
        @test nt_geo.min_y == -90.0
        @test nt_geo.max_y == 90.0
    end

    @testset "Round-trip conversion" begin
        # Test NamedTuple → Extent → NamedTuple
        nt_orig = (min_x=5.0, min_y=10.0, max_x=15.0, max_y=20.0)
        nt_roundtrip = GGA.extent2nt(GGA.nt2extent(nt_orig))
        @test nt_roundtrip.min_x == nt_orig.min_x
        @test nt_roundtrip.max_x == nt_orig.max_x
        @test nt_roundtrip.min_y == nt_orig.min_y
        @test nt_roundtrip.max_y == nt_orig.max_y

        # Test Extent → NamedTuple → Extent
        ext_orig = Extent(X=(100.0, 200.0), Y=(300.0, 400.0))
        ext_roundtrip = GGA.nt2extent(GGA.extent2nt(ext_orig))
        @test ext_roundtrip.X == ext_orig.X
        @test ext_roundtrip.Y == ext_orig.Y
    end

    @testset "Edge cases" begin
        # Zero-width extent
        nt_point = (min_x=5.0, min_y=5.0, max_x=5.0, max_y=5.0)
        ext_point = GGA.nt2extent(nt_point)
        @test ext_point.X == (5.0, 5.0)
        @test ext_point.Y == (5.0, 5.0)

        # Very large extent
        nt_large = (min_x=-1e6, min_y=-1e6, max_x=1e6, max_y=1e6)
        ext_large = GGA.nt2extent(nt_large)
        nt_back = GGA.extent2nt(ext_large)
        @test nt_back.min_x == nt_large.min_x
        @test nt_back.max_x == nt_large.max_x
    end
end
