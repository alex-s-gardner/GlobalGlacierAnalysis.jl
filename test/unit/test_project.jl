using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using DimensionalData

@testset "Project configuration" begin
    @testset "gemb_info" begin
        # Every registered run must be constructible: the pscale and Δheight grids carry one
        # label per value, and filename_gemb_combined resolves its interpolated id.
        for gemb_run_id in 1:6
            g = GGA.gemb_info(; gemb_run_id)

            @test length(g.precipitation_scale) == length(dims(g.precipitation_scale, :pscale))
            @test length(g.elevation_delta) == length(dims(g.elevation_delta, :Δheight))

            @test !occursin("\$", g.filename_gemb_combined)
            @test occursin(g.file_uniqueid, g.filename_gemb_combined)
            @test endswith(g.filename_gemb_combined, ".jld2")

            @test g.modify_melt_only isa Bool
        end

        @test_throws "unrecognized gemb_run_id" GGA.gemb_info(; gemb_run_id=7)
    end

    @testset "file_is_current" begin
        absent = joinpath(mktempdir(), "not_written.jld2")

        # A missing file is never current, and must not compare an mtime against `nothing`.
        @test GGA.file_is_current(absent, nothing) == false
        @test GGA.file_is_current(absent, DateTime(2020, 1, 1)) == false

        present = tempname()
        write(present, "x")
        try
            # No rebuild date means any existing file is reusable.
            @test GGA.file_is_current(present, nothing) == true

            # Written after the rebuild date -> reusable; written before -> rebuild.
            @test GGA.file_is_current(present, DateTime(2000, 1, 1)) == true
            @test GGA.file_is_current(present, DateTime(2999, 1, 1)) == false
        finally
            rm(present; force=true)
        end
    end

    @testset "project_paths covers every product" begin
        products = GGA.project_products(project_id=:v01)
        paths = GGA.project_paths(project_id=:v01)
        @test keys(paths) == keys(products)
        for mission in keys(products)
            @test occursin(string(products[mission].name), paths[mission].geotile)
        end
    end

    @testset "date bins" begin
        date_range, date_center = GGA.project_date_bins()
        decyear = GGA.project_decyear_bins()

        # One bin per centre, in both representations. Observations are grouped by the decimal-year
        # edges while the result is labelled with `date_center`, so a mismatch here would misalign
        # every binned array against its own date dimension. The two used to be independent literals
        # -- `Date(1990):Day(30):Date(2026,1,1)` and `1990:(30/365):2026` -- which is how the record
        # end silently stopped short of the archive.
        @test length(date_center) == length(date_range) - 1
        @test length(decyear) == length(date_range)

        # Edges keep the original spacing, so extending the record only appends bins. If this changes,
        # every observation moves to a different bin and published values shift.
        @test first(decyear) == 1990.0
        @test step(decyear) ≈ 30 / 365
        @test decyear[100] ≈ 1990.0 + 99 * (30 / 365)

        # The grid must reach the end of the longest archive: ICESat-2 ATL06 v7 runs to 2026-05-18.
        @test last(date_range) > Date(2026, 5, 18)

        # 30-day spacing, centres half a bin in
        @test all(diff(date_range) .== Day(30))
        @test date_center[1] == date_range[1] + Day(15)
    end
end
