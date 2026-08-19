"""
Tests for archive-building helpers in utilities_build_archive.jl

These cover the atomic-write, download-pooling, work-partitioning and retry logic used by
`geotile_build_archive`. All four are pure or filesystem-only, so they are exercised without
network access or real HDF5 granules.

Regression context: the ICESat-2 ATL06 archive was lost to an unbounded retry loop that wrote a
`tempname` stub into the output directory on every failed attempt, leaving 633,838 orphaned
8-byte files and no geotiles. The `atomic_write` and `download_with_retry!` tests below pin the two
behaviours that prevent a repeat.
"""

using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using Arrow
using DataFrames

# Minimal stand-in for a SpaceAltimetry/SpaceLiDAR granule. Only `id` and `url` are touched by
# the code under test, and `url` must be mutable because the download step rewrites it to a
# local path.
mutable struct MockGranule
    id::String
    url::String
end

@testset "Archive building helpers" begin

    @testset "atomic_write" begin
        mktempdir() do temp_dir
            outdir = joinpath(temp_dir, "geotile", "2deg")
            mkpath(outdir)
            target = joinpath(outdir, "lat[+60+62]lon[-140-138].arrow")

            # Test: a successful write lands at the target and leaves no scratch file behind
            GGA.atomic_write(target) do tmp
                @test dirname(tmp) == joinpath(outdir, "tmp")   # scratch is beside, not among
                Arrow.write(tmp, DataFrame(height=[1.0, 2.0]))
            end

            @test isfile(target)
            @test nrow(DataFrame(Arrow.Table(target))) == 2
            @test filter(!=("tmp"), readdir(outdir)) == ["lat[+60+62]lon[-140-138].arrow"]
            @test isempty(readdir(joinpath(outdir, "tmp")))

            # Test: a FAILED write leaves nothing behind at all -- not in the output directory and
            # not in tmp/. An 8-byte ARROW1 stub surviving here is what destroyed the ATL06
            # archive, so this is the invariant that matters most.
            @test_throws ErrorException GGA.atomic_write(target) do tmp
                write(tmp, "ARROW1\0\0")      # a partial write, as Arrow.write would leave
                error("simulated interruption")
            end

            @test isempty(readdir(joinpath(outdir, "tmp")))
            @test !any(startswith.(readdir(outdir), "jl_"))

            # Test: the pre-existing target is untouched by the failed write
            @test isfile(target)
            @test nrow(DataFrame(Arrow.Table(target))) == 2

            # Test: suffix is applied, for writers that dispatch on extension
            target2 = joinpath(outdir, "with_suffix.arrow")
            GGA.atomic_write(target2; suffix=".arrow") do tmp
                @test endswith(tmp, ".arrow")
                Arrow.write(tmp, DataFrame(a=[1]))
            end
            @test isfile(target2)
            @test isempty(readdir(joinpath(outdir, "tmp")))
        end
    end

    @testset "partition_rows" begin
        # Test: nothing selects everything
        @test collect(GGA.partition_rows(10, nothing)) == collect(1:10)
        @test collect(GGA.partition_rows(0, nothing)) == Int[]

        # Test: single partition is the identity
        @test collect(GGA.partition_rows(7, (1, 1))) == collect(1:7)

        # Test: partitions tile the row range exactly -- no gaps, no overlap. This is the
        # property that makes concurrent build processes safe.
        for nrows in (0, 1, 5, 8, 9, 1397), nparts in (1, 3, 8)
            covered = reduce(vcat, [collect(GGA.partition_rows(nrows, (i, nparts))) for i in 1:nparts])
            @test sort(covered) == collect(1:nrows)
            @test length(unique(covered)) == nrows
        end

        # Test: stride (not contiguous block) so each partition spans all longitudes. With
        # geotile_granules sorted by longitude, a contiguous block would hand one process every
        # slow high-latitude tile.
        @test collect(GGA.partition_rows(9, (1, 3))) == [1, 4, 7]
        @test collect(GGA.partition_rows(9, (2, 3))) == [2, 5, 8]
        @test collect(GGA.partition_rows(9, (3, 3))) == [3, 6, 9]

        # Test: partitions stay balanced to within one row
        sizes = [length(GGA.partition_rows(1397, (i, 8))) for i in 1:8]
        @test maximum(sizes) - minimum(sizes) <= 1
        @test sum(sizes) == 1397

        # Test: out-of-range and degenerate partitions are rejected rather than silently
        # skipping geotiles
        @test_throws ErrorException GGA.partition_rows(10, (0, 4))
        @test_throws ErrorException GGA.partition_rows(10, (5, 4))
        @test_throws ErrorException GGA.partition_rows(10, (1, 0))
    end

    @testset "pool_granules" begin
        g(name) = MockGranule(name, "https://example.org/$name")

        @testset "deduplicates granules shared between geotiles" begin
            # Adjacent geotiles share granules; the per-geotile loop this replaced requested the
            # same URL once per tile.
            lists = [[g("A.h5"), g("B.h5")], [g("B.h5"), g("C.h5")], [g("C.h5")]]
            granules, n_total = GGA.pool_granules(lists, Set{String}())

            @test n_total == 5                      # references seen
            @test length(granules) == 3             # unique files fetched
            @test [x.id for x in granules] == ["A.h5", "B.h5", "C.h5"]
        end

        @testset "skips files already on disk" begin
            lists = [[g("A.h5"), g("B.h5")], [g("C.h5")]]
            granules, n_total = GGA.pool_granules(lists, Set(["B.h5"]))

            @test n_total == 3
            @test [x.id for x in granules] == ["A.h5", "C.h5"]
        end

        @testset "returns nothing to do when everything is present" begin
            lists = [[g("A.h5")], [g("A.h5"), g("B.h5")]]
            granules, n_total = GGA.pool_granules(lists, Set(["A.h5", "B.h5"]))

            @test isempty(granules)
            @test n_total == 3
        end

        @testset "handles empty input" begin
            granules, n_total = GGA.pool_granules(Vector{Vector{MockGranule}}(), Set{String}())
            @test isempty(granules)
            @test n_total == 0

            granules, n_total = GGA.pool_granules([MockGranule[], MockGranule[]], Set{String}())
            @test isempty(granules)
            @test n_total == 0
        end

        @testset "download order is reproducible" begin
            # Chunk boundaries must be stable across runs so an interrupted download resumes
            # predictably. Geotiles are visited in order and granules within a geotile in order, so
            # input order is already deterministic -- no key sort needed.
            lists = [[g("Z.h5"), g("M.h5")], [g("A.h5"), g("M.h5")]]
            first_call = [x.id for x in GGA.pool_granules(lists, Set{String}())[1]]
            second_call = [x.id for x in GGA.pool_granules(lists, Set{String}())[1]]

            @test first_call == second_call
            @test first_call == ["Z.h5", "M.h5", "A.h5"]   # first-seen order, duplicates dropped
        end
    end

    @testset "download_with_retry!" begin
        granule = MockGranule("ATL06_20190101000000_00000000_006_01.h5", "https://example.org/x.h5")

        @testset "returns on first success" begin
            calls = Ref(0)
            ok(g, dir) = (calls[] += 1; nothing)

            @test GGA.download_with_retry!(granule, "/tmp"; downloader=ok, backoff=0) === granule
            @test calls[] == 1
        end

        @testset "retries a transient failure then succeeds" begin
            calls = Ref(0)
            flaky = (g, dir) -> begin
                calls[] += 1
                calls[] < 3 && error("simulated transient network failure")
                nothing
            end

            @test GGA.download_with_retry!(granule, "/tmp"; downloader=flaky, backoff=0) === granule
            @test calls[] == 3
        end

        @testset "gives up after max_attempts=$n and rethrows" for n in (1, 4)
            # The bug this replaces looped forever on a deterministic failure. The retry must be
            # bounded, must attempt exactly max_attempts times, and must surface the exception
            # rather than swallowing it.
            calls = Ref(0)
            always_fails = (g, dir) -> (calls[] += 1; error("granule withdrawn"))

            @test_throws ErrorException GGA.download_with_retry!(
                granule, "/tmp"; max_attempts=n, downloader=always_fails, backoff=0
            )
            @test calls[] == n
        end

        @testset "passes savedir through to the downloader" begin
            seen = Ref("")
            capture = (g, dir) -> (seen[] = dir; nothing)

            GGA.download_with_retry!(granule, "/data/raw"; downloader=capture, backoff=0)
            @test seen[] == "/data/raw"
        end
    end

    @testset "geotile_build_archive stage validation" begin
        # Stage names are validated before any path or network work, so a typo fails fast
        # instead of silently skipping a stage of a multi-day run.
        @test_throws ErrorException GGA.geotile_build_archive(; stages=(:bogus,))
        @test_throws ErrorException GGA.geotile_build_archive(; stages=(:search, :buidl))
    end
end
