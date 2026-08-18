using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using Extents

@testset "Coordinate Operations" begin
    @testset "geotile_extent from center" begin
        # Test standard geotile centered at (45, -120) with width 2
        ext = GGA.geotile_extent(45, -120, 2)
        @test ext.min_x == -121
        @test ext.max_x == -119
        @test ext.min_y == 44
        @test ext.max_y == 46

        # Test at equator
        ext_eq = GGA.geotile_extent(0, 0, 2)
        @test ext_eq.min_x == -1
        @test ext_eq.max_x == 1
        @test ext_eq.min_y == -1
        @test ext_eq.max_y == 1

        # Test at pole
        ext_pole = GGA.geotile_extent(89, 0, 2)
        @test ext_pole.min_y == 88
        @test ext_pole.max_y == 90
    end

    @testset "geotile_extent from ID string" begin
        # Test standard geotile ID parsing
        ext = GGA.geotile_extent("lat+45+47lon-123-121")
        @test ext.X == (-123, -121)
        @test ext.Y == (45, 47)

        # Test Southern hemisphere
        ext_south = GGA.geotile_extent("lat-45-43lon+010+012")
        @test ext_south.X == (10, 12)
        @test ext_south.Y == (-45, -43)

        # Test Western hemisphere negative longitude
        ext_west = GGA.geotile_extent("lat+01+03lon-180-178")
        @test ext_west.X == (-180, -178)
        @test ext_west.Y == (1, 3)

        # Test zero crossing
        ext_zero = GGA.geotile_extent("lat-01+01lon-001+001")
        @test ext_zero.X == (-1, 1)
        @test ext_zero.Y == (-1, 1)
    end

    @testset "within extent" begin
        # Test with NamedTuple extent
        ext_nt = (min_x=-121.0, min_y=44.0, max_x=-119.0, max_y=46.0)

        # Interior point
        @test GGA.within(ext_nt, -120.0, 45.0) == true

        # Boundary points (should be inclusive)
        @test GGA.within(ext_nt, -121.0, 44.0) == true  # min corner
        @test GGA.within(ext_nt, -119.0, 46.0) == true  # max corner

        # Exterior point
        @test GGA.within(ext_nt, -122.0, 45.0) == false
        @test GGA.within(ext_nt, -120.0, 47.0) == false

        # Test with Extent type
        ext_obj = Extent(X=(-121.0, -119.0), Y=(44.0, 46.0))
        @test GGA.within(ext_obj, -120.0, 45.0) == true
        @test GGA.within(ext_obj, -122.0, 45.0) == false
    end

    @testset "utm_epsg" begin
        # Standard UTM zones (Northern hemisphere)
        @test GGA.utm_epsg(0.0, 45.0) == 32631  # Zone 31N (prime meridian)
        @test GGA.utm_epsg(3.0, 45.0) == 32631  # Still zone 31N
        @test GGA.utm_epsg(9.0, 45.0) == 32632  # Zone 32N

        # Southern hemisphere
        @test GGA.utm_epsg(9.0, -45.0) == 32732  # Zone 32S

        # Norway special case (zone 32 extends west)
        @test GGA.utm_epsg(5.0, 60.0) == 32632  # Norway exception

        # Svalbard special case (zones 31 & 33 omitted)
        @test GGA.utm_epsg(15.0, 75.0) == 32633  # Svalbard zone 33X

        # North polar region (lat > 84.5)
        @test GGA.utm_epsg(0.0, 85.0) == 3413  # NSIDC North
        @test GGA.utm_epsg(0.0, 89.0) == 3413

        # South polar region (lat < -80)
        @test GGA.utm_epsg(0.0, -81.0) == 3031  # Antarctic Polar Stereographic
        @test GGA.utm_epsg(0.0, -85.0) == 3031
    end

    @testset "meters2lonlat_distance" begin
        # At equator, 1 degree ≈ 111,320 m
        lon_deg_eq = GGA.meters2lonlat_distance(111320.0, 0.0)
        @test lon_deg_eq ≈ 1.0 rtol=1e-3

        # At 60° latitude, cosine effect doubles the degrees needed
        lon_deg_60 = GGA.meters2lonlat_distance(111320.0, 60.0)
        @test lon_deg_60 ≈ 2.0 rtol=1e-2

        # At 45° latitude
        lon_deg_45 = GGA.meters2lonlat_distance(111320.0, 45.0)
        @test lon_deg_45 ≈ 1.0 / cos(deg2rad(45.0)) rtol=1e-2

        # Near pole (should handle gracefully)
        lon_deg_89 = GGA.meters2lonlat_distance(111320.0, 89.0)
        @test lon_deg_89 > 10.0  # Much larger factor near pole
    end

    @testset "vector_overlap" begin
        # Complete overlap
        @test GGA.vector_overlap([1, 5], [1, 5]) == true

        # Partial overlap
        @test GGA.vector_overlap([1, 4], [3, 6]) == true
        @test GGA.vector_overlap([3, 6], [1, 4]) == true

        # No overlap
        @test GGA.vector_overlap([1, 2], [3, 4]) == false
        @test GGA.vector_overlap([3, 4], [1, 2]) == false

        # Touching boundaries (check inclusive behavior)
        # Assuming boundaries are inclusive
        @test GGA.vector_overlap([1, 3], [3, 5]) == true

        # One range contains the other
        @test GGA.vector_overlap([1, 10], [3, 5]) == true
        @test GGA.vector_overlap([3, 5], [1, 10]) == true
    end
end
