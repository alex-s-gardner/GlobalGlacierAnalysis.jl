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
end
