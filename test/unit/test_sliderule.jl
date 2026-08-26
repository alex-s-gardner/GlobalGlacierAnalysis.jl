"""
Tests for the SlideRule ingest path in utilities_sliderule.jl

SlideRule returns its results as a native record stream carrying an Arrow IPC file, so the pieces
under test are the framing parser, the translation from SlideRule's column set to this project's
archive schema, and the incremental rule that decides which granules to ask for. HTTP is injected, so
nothing here touches the network.

The ATL06 schema mapping was validated against the HDF5 path on a real geotile
(`lat[-80-78]lon[+166+168]`, 3 granules, 45,466 coordinate-matched points): height, height_error,
height_reference, quality, track, strong_beam, detector_id all agreed exactly, and `datetime` differed
by exactly -18.000 s, which is the leap-second offset pinned below.

The GEDI mapping was validated the same way on `lat[+52+54]lon[-170-168]`: across 53
coordinate-matched points, `height`, `height_error`, `intensity`, `sensitivity`, `sun_angle`,
`height_reference`, `quality`, `surface`, `nmodes`, `track`, `strong_beam`, `classification` and the
granule `id` were all bit-identical, and `datetime` agreed to within 1 ms of the same 18 s offset.
"""

using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using Arrow
using DataFrames
using Dates
using Extents
using JSON
using Logging

# Build a SlideRule record: 8-byte big-endian header (version, type_size, data_size), the
# NUL-terminated type string, then the payload.
function mock_record(rectype::String, payload::Vector{UInt8})
    type_bytes = vcat(Vector{UInt8}(rectype), 0x00)
    out = UInt8[]
    append!(out, reinterpret(UInt8, [hton(Int16(1))]))
    append!(out, reinterpret(UInt8, [hton(Int16(length(type_bytes)))]))
    append!(out, reinterpret(UInt8, [hton(Int32(length(payload)))]))
    append!(out, type_bytes)
    append!(out, payload)
    return out
end

# arrowrec payloads open with a fixed 128-byte NUL-padded filename field
function mock_filename_field(name::String)
    field = zeros(UInt8, GGA.SLIDERULE_FILENAME_LEN)
    bytes = Vector{UInt8}(name)
    field[1:length(bytes)] = bytes
    return field
end

mock_meta(name::String, size::Integer) =
    mock_record("arrowrec.meta", vcat(mock_filename_field(name), reinterpret(UInt8, [Int64(size)])))

mock_data(name::String, payload::Vector{UInt8}) =
    mock_record("arrowrec.data", vcat(mock_filename_field(name), payload))

mock_except(text::String) =
    mock_record("exceptrec", vcat(reinterpret(UInt8, [Int32(0), UInt32(1)]), Vector{UInt8}(text), 0x00))

# A real Arrow IPC file, so the assembled bytes can be handed to Arrow.Table like the live path does.
function mock_feather(df::DataFrame)
    io = IOBuffer()
    Arrow.write(io, df)
    return take!(io)
end

"""
    gedi_raw_frame(; n=7, kwargs...) -> DataFrame

A `gedi02ax` response as SlideRule returns it: native columns plus the 21 ancillary columns, named
with the exact strings that were requested, slashes included.

With `n=7` (the default) row 1 passes the L3 quality filter and rows 2-7 each fail exactly one term,
in the order the filter applies them. With `n=1` only the passing row is produced, and any keyword
overrides that row's value -- which is how each term is exercised in isolation. `n=0` gives the
correctly typed empty frame, for a granule that returns no points at all.
"""
function gedi_raw_frame(; n=7,
    flags=GGA.GEDI_SURFACE_FLAG_MASK | GGA.GEDI_L2_QUALITY_FLAG_MASK,
    rx_quality=0x01, stale=0x00, maxamp=80.0f0, sd=1.0f0, algorithm=0x01,
    zcross=10.0f0, toploc=20.0f0, dem=NaN32)

    df = DataFrame(
        longitude=[86.5], latitude=[28.5],
        elevation_lm=Float32[100.0], elevation_hr=Float32[110.0],
        solar_elevation=Float32[12.0], sensitivity=Float32[0.9],
        time_ns=[DateTime(2019, 1, 1)],
        shot_number=UInt64[1753_00_00_3_00518838],
        orbit=UInt32[1753], track=UInt16[1683], beam=UInt8[5],
        flags=UInt8[flags], srcid=Int32[0],
        elevation_bin0_error=Float32[0.5], energy_total=Float32[400.0],
        num_detectedmodes=UInt8[3], digital_elevation_model=Float32[dem],
    )
    df[!, Symbol("rx_assess/quality_flag")] = UInt8[rx_quality]
    df[!, Symbol("geolocation/stale_return_flag")] = UInt8[stale]
    df[!, Symbol("rx_assess/rx_maxamp")] = Float32[maxamp]
    df[!, Symbol("rx_assess/sd_corrected")] = Float32[sd]
    df[!, :selected_algorithm] = UInt8[algorithm]
    for a in 1:6
        df[!, Symbol("rx_processing_a$(a)/zcross")] = Float32[zcross]
        df[!, Symbol("rx_processing_a$(a)/toploc")] = Float32[toploc]
    end

    n == 0 && return df[1:0, :]
    n == 1 && return df

    # one row per filter term, each a copy of the passing row with that term's input spoiled
    rows = [df]
    for spoil in (
        d -> d[1, Symbol("rx_assess/quality_flag")] = 0x00,
        d -> d[1, Symbol("geolocation/stale_return_flag")] = 0x01,
        d -> d[1, Symbol("rx_assess/rx_maxamp")] = 7.0f0,
        d -> d[1, Symbol("rx_processing_a1/zcross")] = 0.0f0,
        d -> d[1, Symbol("rx_processing_a1/toploc")] = 0.0f0,
        d -> d[1, :flags] = GGA.GEDI_SURFACE_FLAG_MASK | GGA.GEDI_DEGRADE_FLAG_MASK,
    )
        row = copy(df)
        spoil(row)
        push!(rows, row)
    end
    return reduce(vcat, rows)[1:n, :]
end

"""
    anonymous(f)

Run `f` with no SlideRule credentials reachable: environment variable unset, token file pointed at a
path that does not exist, and the cached bearer token cleared.

Credentials are process-global and cached, so without this the suite's behaviour depends on whether
the machine running it has a real `~/.sliderule_pat` -- and if it does, tests that expect anonymous
access reach out to the live login endpoint instead. The cache is restored afterwards so an
authenticated session in the same process is not disturbed.
"""
function anonymous(f)
    saved = GGA.SLIDERULE_TOKEN[]
    GGA.SLIDERULE_TOKEN[] = nothing
    try
        withenv("SLIDERULE_GITHUB_TOKEN" => nothing,
            "SLIDERULE_PAT_FILE" => joinpath(tempdir(), "no-such-sliderule-pat-$(getpid())")) do
            f()
        end
    finally
        GGA.SLIDERULE_TOKEN[] = saved
    end
end

@testset "SlideRule ingest" begin

    @testset "host selection" begin
        withenv("SLIDERULE_ORGANIZATION" => nothing) do
            @test GGA.sliderule_host() == GGA.SLIDERULE_PUBLIC_HOST
        end
        withenv("SLIDERULE_ORGANIZATION" => "") do
            @test GGA.sliderule_host() == GGA.SLIDERULE_PUBLIC_HOST
        end
        withenv("SLIDERULE_ORGANIZATION" => "my-cluster") do
            @test GGA.sliderule_host() == "my-cluster.slideruleearth.io"
        end
    end

    @testset "anonymous access needs no token" begin
        # The public cluster works unauthenticated; absence of a PAT must not raise.
        #
        # `SLIDERULE_PAT_FILE` is pointed at nothing and the token cache cleared, because otherwise
        # this depends on whether the machine happens to have `~/.sliderule_pat`. It did not when this
        # was written, and once a real token was installed the test exchanged it for a JWT against the
        # live login endpoint -- a network call from a suite that promises not to make any, with the
        # bearer token printed into the failure output.
        anonymous() do
            @test GGA.sliderule_pat() == ""
            @test GGA.sliderule_token() === nothing
        end
    end

    @testset "record framing" begin
        payload = Vector{UInt8}("hello")
        stream = vcat(mock_record("arrowrec.data", payload), mock_record("exceptrec", Vector{UInt8}("x")))
        records = GGA.sliderule_records(stream)
        @test length(records) == 2
        @test records[1][1] == "arrowrec.data"
        @test records[1][2] == payload
        @test records[2][1] == "exceptrec"

        # Test: a truncated transfer is an error, not a short read. Returning the records that
        # happened to parse would present partial data as a complete result.
        truncated = stream[1:end-3]
        @test_throws GGA.NonRetryable GGA.sliderule_records(truncated)

        # Test: trailing bytes too short to be a header are also rejected
        @test_throws GGA.NonRetryable GGA.sliderule_records(vcat(stream, UInt8[0x01, 0x02]))

        # Test: a nonsense header is caught rather than used to index out of bounds
        bogus = vcat(reinterpret(UInt8, [hton(Int16(1)), hton(Int16(-4))]), reinterpret(UInt8, [hton(Int32(0))]))
        @test_throws GGA.NonRetryable GGA.sliderule_records(bogus)
    end

    @testset "arrow reassembly" begin
        name = "out.feather"
        df = DataFrame(a=[1, 2, 3], b=["x", "y", "z"])
        file = mock_feather(df)

        # Test: single-chunk transfer round-trips through Arrow
        stream = vcat(mock_meta(name, length(file)), mock_data(name, file))
        assembled = GGA.sliderule_arrow(stream; filename=name)
        @test assembled == file
        @test DataFrame(Arrow.Table(IOBuffer(assembled))) == df

        # Test: a chunked transfer concatenates in order. The live service splits large files, so
        # this is the normal case, not an edge case.
        half = length(file) ÷ 2
        chunked = vcat(mock_meta(name, length(file)),
            mock_data(name, file[1:half]), mock_data(name, file[half+1:end]))
        @test GGA.sliderule_arrow(chunked; filename=name) == file

        # Test: no output file at all -- every requested granule missed the region -- is `nothing`,
        # not an error. This is common: CMR adds geolocation margin to granule polygons.
        empty_stream = vcat(mock_except("no data for granule"), mock_except("still nothing"))
        @test GGA.sliderule_arrow(empty_stream) === nothing

        # Test: a size that disagrees with what arrived is an error. Silently accepting it would
        # write a half-granule geotile that looks complete.
        @test_throws GGA.NonRetryable GGA.sliderule_arrow(
            vcat(mock_meta(name, length(file) + 500), mock_data(name, file)); filename=name)

        # Test: a filename we did not ask for means the payload layout is not what we assume
        @test_throws GGA.NonRetryable GGA.sliderule_arrow(
            vcat(mock_meta("other.feather", length(file)), mock_data("other.feather", file)); filename=name)

        # Test: interleaved exceptrec records do not corrupt the assembled file
        interleaved = vcat(mock_except("note"), mock_meta(name, length(file)),
            mock_except("another"), mock_data(name, file))
        @test GGA.sliderule_arrow(interleaved; filename=name) == file
    end

    @testset "alert text" begin
        @test GGA.sliderule_alert(vcat(reinterpret(UInt8, [Int32(0), UInt32(1)]),
            Vector{UInt8}("granule empty"), 0x00)) == "granule empty"
        @test GGA.sliderule_alert(UInt8[]) == ""
    end

    @testset "request parameters" begin
        extent = Extent(X=(8.0, 10.0), Y=(46.0, 48.0))
        parms = GGA.sliderule_atl06_parms(extent; granules=["ATL06_a.h5", "ATL06_b.h5"],
            t0=DateTime(2019, 1, 1), t1=DateTime(2019, 3, 1), filename="f.feather")

        @test parms["asset"] == "icesat2-atl06"
        @test parms["output"]["format"] == "feather"   # Arrow IPC, readable without a parquet dep
        @test parms["output"]["path"] == "f.feather"
        @test parms["resources"] == ["ATL06_a.h5", "ATL06_b.h5"]
        @test parms["t0"] == "2019-01-01T00:00:00Z"
        @test parms["t1"] == "2019-03-01T00:00:00Z"
        # height_reference comes from an ancillary field, not the default x-series column set
        @test parms["atl06_fields"] == ["dem/dem_h"]

        # Test: polygon is a closed ring covering the extent
        poly = parms["poly"]
        @test length(poly) == 5
        @test poly[1] == poly[end]
        @test extrema(p["lon"] for p in poly) == (8.0, 10.0)
        @test extrema(p["lat"] for p in poly) == (46.0, 48.0)

        # Test: without granules the server does its own CMR query, so no resources key is sent
        @test !haskey(GGA.sliderule_atl06_parms(extent), "resources")
    end

    @testset "granule id parsing" begin
        @test GGA._atl06_key("ATL06_20181017143645_02900102_007_01.h5") == (290, 1, 2)
        @test GGA._atl06_key("/some/dir/ATL06_20250101000000_12345614_007_02.h5") == (1234, 56, 14)
        @test GGA._atl06_key("not-a-granule.h5") === nothing
        @test GGA._atl06_key("ATL06_20181017143645_0290010_007_01.h5") === nothing
    end

    @testset "beam mapping" begin
        # SlideRule reports gt as 10, 20 ... 60; the archive stores the beam name
        @test GGA._beam_name.(10:10:60) == ["gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"]
        @test GGA._beam_name(0) == ""
        @test GGA._beam_name(70) == ""
    end

    @testset "ancillary DEM column name" begin
        # SlideRule labels an ancillary column with the exact string requested, so this arrives as
        # "dem/dem_h". Reading `:dem_h` instead yields an all-NaN height_reference that looks like
        # missing data rather than a bug -- which is exactly how it went unnoticed once.
        @test GGA.SLIDERULE_ATL06_FIELDS == ["dem/dem_h"]
        @test first(GGA.SLIDERULE_DEM_COLUMNS) == Symbol("dem/dem_h")

        n = 2
        as_requested = DataFrame(Symbol("dem/dem_h") => Float32[10.0, 20.0])
        @test GGA._dem_column(as_requested) == Float32[10.0, 20.0]

        # a bare `dem_h` is still accepted, in case the server ever normalises the name
        @test GGA._dem_column(DataFrame(dem_h=Float32[30.0, 40.0])) == Float32[30.0, 40.0]

        # absent entirely: NaN, but noisily
        missing_col = DataFrame(longitude=[1.0, 2.0])
        @test all(isnan, (@test_logs (:warn,) match_mode = :any GGA._dem_column(missing_col)))
        @test length(GGA._dem_column(missing_col)) == n
    end

    @testset "schema mapping" begin
        granules = ["ATL06_20181017143645_02900102_007_01.h5", "ATL06_20181017152951_02900110_007_01.h5"]
        raw = DataFrame(
            longitude=[8.5, 8.6, 8.7],
            latitude=[46.5, 46.6, 46.7],
            h_li=Float32[1000.0, 2000.0, GGA.ATL06_FILL_VALUE],
            h_li_sigma=Float32[3.0, 0.0, 1.0],
            sigma_geo_h=Float32[4.0, 0.0, 1.0],
            time_ns=[DateTime(2019, 1, 1, 0, 0, 0), DateTime(2019, 1, 1, 0, 0, 1), DateTime(2019, 1, 1, 0, 0, 2)],
            atl06_quality_summary=Int8[0, 1, 0],
            spot=UInt8[1, 2, 3],
            gt=UInt8[10, 20, 30],
            rgt=UInt16[290, 290, 999],
            cycle=UInt16[1, 1, 1],
            region=UInt8[2, 10, 2],
        )
        # named as SlideRule returns it, not as a bare `dem_h`
        raw[!, Symbol("dem/dem_h")] = Float32[900.0, 1900.0, GGA.ATL06_FILL_VALUE]

        out = GGA.sliderule2archive(raw, granules)

        @test names(out) == names(GGA.sliderule_empty_table())
        @test eltype.(eachcol(out)) == eltype.(eachcol(GGA.sliderule_empty_table()))

        # Test: height_error combines both terms the way the HDF5 reader does
        @test out.height_error[1] ≈ sqrt(3.0f0^2 + 4.0f0^2)

        # Test: fill values become NaN, since the pipeline tests for NaN not for 3.4e38
        @test isnan(out.height[3])
        @test isnan(out.height_reference[3])
        @test out.height_reference[1] == 900.0f0

        # Test: quality is inverted -- ATL06 uses 0 for "no problems", the archive uses true for good
        @test out.quality == [true, false, true]

        # Test: strong beams are the odd spots, independent of spacecraft orientation
        @test out.strong_beam == [true, false, true]
        @test out.detector_id == Int8[1, 2, 3]
        @test out.track == ["gt1l", "gt1r", "gt2l"]

        # Test: granule id recovered from (rgt, cycle, region); an unmatched triple yields ""
        @test out.id == [granules[1], granules[2], ""]

        # Test: timestamps land on the archive's GPS timescale, 18 s ahead of SlideRule's UTC
        @test out.datetime[1] == DateTime(2019, 1, 1, 0, 0, 18)
        @test GGA.sliderule2archive(raw, granules; gps_time=false).datetime[1] == DateTime(2019, 1, 1)

        # Test: an empty result still has the archive's columns and types, so it vcats and serialises
        empty_out = GGA.sliderule2archive(DataFrame(), granules)
        @test names(empty_out) == names(GGA.sliderule_empty_table())
        @test isempty(empty_out)
    end

    @testset "query with injected transport" begin
        extent = Extent(X=(8.0, 10.0), Y=(46.0, 48.0))
        granules = ["ATL06_20181017143645_02900102_007_01.h5"]
        raw = DataFrame(
            longitude=[8.5], latitude=[46.5], h_li=Float32[100.0], h_li_sigma=Float32[1.0],
            sigma_geo_h=Float32[1.0], time_ns=[DateTime(2019, 1, 1)],
            atl06_quality_summary=Int8[0], spot=UInt8[1], gt=UInt8[10],
            rgt=UInt16[290], cycle=UInt16[1], region=UInt8[2])
        raw[!, Symbol("dem/dem_h")] = Float32[90.0]

        requests = String[]
        function poster(url, headers, body)
            push!(requests, body)
            file = mock_feather(raw)
            stream = vcat(mock_meta("gga_atl06.feather", length(file)),
                mock_data("gga_atl06.feather", file))
            return (; status=200, body=stream)
        end

        df, failed, partial, unreadable, retried = GGA.sliderule_atl06(extent; granules, poster)
        @test nrow(df) == 1
        @test df.id == granules
        @test isempty(failed)
        @test length(requests) == 1
        @test occursin("icesat2-atl06", requests[1])

        # Test: an HTTP error is not retried into oblivion -- it is a request problem
        bad_poster(url, headers, body) = (; status=400, body=Vector{UInt8}("bad request"))
        @test_throws GGA.NonRetryable GGA.sliderule_atl06(extent; granules, poster=bad_poster)

        # Test: resources are split into batches, so one geotile's several hundred granules do not
        # go out as a single request
        many = ["ATL06_2018101714364$(i)_0290010$(i % 10)_007_01.h5" for i in 1:(GGA.SLIDERULE_RESOURCE_CHUNK+5)]
        empty!(requests)
        GGA.sliderule_atl06(extent; granules=many, poster)
        @test length(requests) == 2
    end

    @testset "incremental build" begin
        mktempdir() do dir
            extent = Extent(X=(8.0, 10.0), Y=(46.0, 48.0))
            granules = ["ATL06_20181017143645_02900102_007_01.h5",
                "ATL06_20181017152951_02900110_007_01.h5"]
            geotiles = DataFrame(
                id=["lat[+46+48]lon[+008+010]"],
                extent=[extent],
                granules=[[(id=g, url=g) for g in granules]],
            )

            # Serve points only for the first granule; the second is a granule that intersects the
            # CMR polygon but has no data inside the region.
            asked = Vector{String}[]
            function poster(url, headers, body)
                push!(asked, String.(JSON.parse(body)["parms"]["resources"]))
                rows = DataFrame(
                    longitude=[8.5], latitude=[46.5], h_li=Float32[100.0], h_li_sigma=Float32[1.0],
                    sigma_geo_h=Float32[1.0], time_ns=[DateTime(2019, 1, 1)],
                    atl06_quality_summary=Int8[0], spot=UInt8[1], gt=UInt8[10],
                    rgt=UInt16[290], cycle=UInt16[1], region=UInt8[2])
                rows[!, Symbol("dem/dem_h")] = Float32[90.0]
                file = mock_feather(rows)
                return (; status=200,
                    body=vcat(mock_meta("gga_atl06.feather", length(file)),
                        mock_data("gga_atl06.feather", file)))
            end

            GGA.geotile_build_sliderule(geotiles, dir; poster, ntasks=1)
            outfile = joinpath(dir, "lat[+46+48]lon[+008+010].arrow")
            @test isfile(outfile)
            built = DataFrame(Arrow.Table(outfile))

            # one real point plus a placeholder for the granule that returned nothing
            @test nrow(built) == 2
            @test Set(built.id) == Set(granules)
            @test length(asked) == 1 && Set(asked[1]) == Set(granules)

            # Test: rerunning asks for nothing and changes nothing. The placeholder row is what stops
            # the empty granule from being requested forever. `isequal`, not `==`, because the
            # placeholder row is NaN and NaN != NaN.
            empty!(asked)
            GGA.geotile_build_sliderule(geotiles, dir; poster, ntasks=1)
            @test isempty(asked)
            @test isequal(DataFrame(Arrow.Table(outfile)), built)

            # Test: a newly acquired granule is the only one requested on the next pass
            new_granule = "ATL06_20260518000000_02900103_007_01.h5"
            geotiles.granules = [[(id=g, url=g) for g in vcat(granules, new_granule)]]
            empty!(asked)
            GGA.geotile_build_sliderule(geotiles, dir; poster, ntasks=1)
            @test length(asked) == 1 && asked[1] == [new_granule]
            @test nrow(DataFrame(Arrow.Table(outfile))) > nrow(built)

            # Test: no scratch files left among the outputs
            @test isempty(filter(f -> !endswith(f, ".arrow"), readdir(joinpath(dir, "tmp"))))
        end
    end

    @testset "build rejects an empty granule list" begin
        mktempdir() do dir
            geotiles = DataFrame(id=["lat[+46+48]lon[+008+010]"],
                extent=[Extent(X=(8.0, 10.0), Y=(46.0, 48.0))],
                granules=[NamedTuple[]])
            @test_throws GGA.NonRetryable GGA.geotile_build_sliderule(geotiles, dir)
        end
    end

    @testset "build rejects an unsupported mission" begin
        mktempdir() do dir
            geotiles = DataFrame(id=["lat[+46+48]lon[+008+010]"],
                extent=[Extent(X=(8.0, 10.0), Y=(46.0, 48.0))],
                granules=[[(id="g.h5", url="g.h5")]])
            # :icesat and :hugonnet have no SlideRule query; failing here beats issuing a request
            # that cannot produce the archive schema.
            @test_throws GGA.NonRetryable GGA.geotile_build_sliderule(geotiles, dir; mission=:icesat)
        end
    end

    @testset "token from a credentials file" begin
        # A file keeps the PAT out of shell history and the process table for a build that runs days.
        mktempdir() do dir
            path = joinpath(dir, "pat")
            write(path, "\n  ghp_fromfile  \nignored second line\n")
            chmod(path, 0o600)

            withenv("SLIDERULE_GITHUB_TOKEN" => nothing, "SLIDERULE_PAT_FILE" => path) do
                # leading blank line skipped, surrounding whitespace stripped
                @test GGA.sliderule_pat() == "ghp_fromfile"
            end

            # Test: the environment wins, so a one-off run can override the stored token
            withenv("SLIDERULE_GITHUB_TOKEN" => "ghp_fromenv", "SLIDERULE_PAT_FILE" => path) do
                @test GGA.sliderule_pat() == "ghp_fromenv"
            end

            # Test: a loose-permission file is still used, but noisily. Ignoring it would silently
            # downgrade the build to anonymous, which is the more expensive failure.
            #
            # `geotile_build_sliderule(warnings=false)` calls `Logging.disable_logging(Warn)`, which
            # is global and outlives the testset that triggered it, so re-enable before asserting on
            # a warning -- otherwise this passes or fails depending on testset order.
            Logging.disable_logging(Logging.BelowMinLevel)
            chmod(path, 0o644)
            withenv("SLIDERULE_GITHUB_TOKEN" => nothing, "SLIDERULE_PAT_FILE" => path) do
                @test (@test_logs (:warn,) match_mode = :any GGA.sliderule_pat()) == "ghp_fromfile"
            end

            # Test: absent or empty file is anonymous, not an error. The cache has to be cleared as
            # well, or a token minted by an earlier assertion answers in its place.
            saved = GGA.SLIDERULE_TOKEN[]
            GGA.SLIDERULE_TOKEN[] = nothing
            withenv("SLIDERULE_GITHUB_TOKEN" => nothing,
                "SLIDERULE_PAT_FILE" => joinpath(dir, "nope")) do
                @test GGA.sliderule_pat() == ""
                @test GGA.sliderule_token() === nothing
            end
            GGA.SLIDERULE_TOKEN[] = saved
            empty_path = joinpath(dir, "empty")
            write(empty_path, "\n\n")
            withenv("SLIDERULE_GITHUB_TOKEN" => nothing, "SLIDERULE_PAT_FILE" => empty_path) do
                @test GGA.sliderule_pat() == ""
            end
        end
    end

    @testset "dedicated capacity" begin
        # `user_service` provisioning without SlideRule's Python client: three JSON posts, mirroring
        # Session.scaleout. Requests still go to the ordinary public hostname -- the gateway routes to
        # our nodes on the bearer token -- so only the provisioning calls are new.
        @test GGA.SLIDERULE_MAX_TTL == 720

        # Test: anonymous is not an error. The request path is identical either way, so a run without
        # credentials still works; it just shares the public cluster.
        withenv("SLIDERULE_GITHUB_TOKEN" => nothing, "SLIDERULE_PAT_FILE" => tempname()) do
            GGA.SLIDERULE_TOKEN[] = nothing
            @test GGA.sliderule_service() === nothing
            @test GGA.sliderule_capacity() == 0
            @test (@test_logs (:warn,) match_mode = :any GGA.sliderule_scaleout()) == 0
        end

        # Pretend we hold a token for identity "me", so the provisioning calls can be inspected.
        GGA.SLIDERULE_TOKEN[] = (token="jwt", renew_at=time() + 3600, service="me")
        try
            @test GGA.sliderule_service() == "me"

            calls = Tuple{String,Any}[]
            nodes = Ref(0)
            function poster(url, headers, body)
                parsed = JSON.parse(body)
                push!(calls, (url, parsed))
                payload = if occursin("discovery/status", url)
                    Dict("nodes" => nodes[])
                else
                    nodes[] = 4       # a deploy brings nodes up
                    Dict("status" => "ok")
                end
                # Test: the bearer token is attached, or the provisioner cannot tell who is asking
                @test any(h -> first(h) == "Authorization" && last(h) == "Bearer jwt", headers)
                return (; status=200, body=Vector{UInt8}(JSON.json(payload)))
            end

            # Test: no capacity yet means deploy, against the public cluster under our identity
            available = GGA.sliderule_scaleout(; node_capacity=4, ttl=60, block=false, poster)
            deploy = only(filter(c -> occursin("deploy", first(c)), calls))
            @test first(deploy) == "https://$(GGA.SLIDERULE_PROVISIONER_HOST)/deploy/me"
            @test last(deploy)["node_capacity"] == 4
            @test last(deploy)["ttl"] == 60
            @test last(deploy)["is_public"] === false
            @test last(deploy)["cluster"] == GGA.SLIDERULE_CLUSTER
            @test available == 0     # block=false reports what was up before the request

            # Test: already at capacity extends the TTL instead of deploying again
            empty!(calls)
            GGA.sliderule_scaleout(; node_capacity=4, ttl=60, block=false, poster)
            extend = only(filter(c -> occursin("extend", first(c)), calls))
            @test first(extend) == "https://$(GGA.SLIDERULE_PROVISIONER_HOST)/extend/me"
            @test last(extend)["ttl"] == 60
            @test isempty(filter(c -> occursin("deploy", first(c)), calls))

            # Test: the TTL is clamped to what the provisioner will grant, rather than being rejected
            empty!(calls)
            nodes[] = 0
            GGA.sliderule_scaleout(; node_capacity=4, ttl=10_000, block=false, poster)
            @test only(filter(c -> occursin("deploy", first(c)), calls))[2]["ttl"] == GGA.SLIDERULE_MAX_TTL

            # Test: capacity is read from the discovery API
            @test GGA.sliderule_capacity(; poster) == 4
            status = last(filter(c -> occursin("discovery/status", first(c)), calls))
            @test last(status) == Dict("service" => "me")

            # Test: a provisioner that will not cooperate downgrades the run to the shared cluster
            # rather than aborting it. Dedicated nodes are throughput, not correctness -- and a stuck
            # stack on the service side really did make every deploy fail, which must not cost hours
            # of ingest that would otherwise have completed.
            nodes[] = 0
            err_poster(url, headers, body) = (; status=200,
                body=Vector{UInt8}(JSON.json(Dict("error" => "nope", "error_description" => "quota"))))
            errored = @test_logs (:warn,) match_mode = :any GGA.sliderule_scaleout(; node_capacity=4, block=false, poster=err_poster)
            @test errored == 0

            # Test: an HTTP-level failure from the provisioner is handled the same way. This is the
            # real case: a stack stuck server-side made every deploy return 500.
            http_err(url, headers, body) = (; status=500,
                body=Vector{UInt8}("""{"error": "internal error", "error_description": "AlreadyExistsException"}"""))
            http_failed = @test_logs (:warn,) match_mode = :any GGA.sliderule_scaleout(; node_capacity=4, block=false, poster=http_err)
            @test http_failed == 0

            # Test: and the build still runs when provisioning failed
            ran = @test_logs (:warn,) match_mode = :any GGA.sliderule_keepalive(() -> :ran; node_capacity=4, poster=http_err)
            @test ran == :ran

            # Test: keepalive runs the body and returns its value
            nodes[] = 4
            @test GGA.sliderule_keepalive(() -> :done; node_capacity=4, every=3600, poster) == :done

            # Test: and it still runs the body when there is no dedicated capacity to keep alive
            GGA.SLIDERULE_TOKEN[] = nothing
            withenv("SLIDERULE_GITHUB_TOKEN" => nothing, "SLIDERULE_PAT_FILE" => tempname()) do
                @test (@test_logs (:warn,) match_mode = :any GGA.sliderule_keepalive(() -> :ran)) == :ran
            end
        finally
            GGA.SLIDERULE_TOKEN[] = nothing
        end
    end

    @testset "failed resource detection" begin
        # A per-beam read can fail while the request still succeeds, returning that beam with zero
        # rows. Indistinguishable in the data from a beam with no points in the region, so it has to
        # be recovered from the message stream or the granule is silently frozen into the archive.
        stream = vcat(
            mock_except("(dataframe.lua:108) request <X> on A.h5 generated dataframe [beam3] with 61 rows"),
            mock_except("Failure on resource A.h5 beam beam11: H5Coro::Future read failure on BEAM1011/shot_number"),
            mock_except("Failure on resource B.h5 beam beam5: H5Coro::Future read failure on BEAM0101/quality_flag"),
            mock_except("Successfully completed processing resource [1 out of 2]: A.h5"))
        @test GGA.sliderule_failed_resources(stream) == Set(["A.h5", "B.h5"])

        # Test: an informative message is not a failure. Most exceptrecs are "no data for this
        # granule" notices, which are normal -- CMR pads granule polygons.
        @test isempty(GGA.sliderule_failed_resources(
            vcat(mock_except("generated dataframe [beam0] with 0 rows"),
                mock_except("Successfully completed processing resource [2 out of 2]: B.h5"))))
        @test isempty(GGA.sliderule_failed_resources(UInt8[]))

        @test GGA._sliderule_failed_resource("Failure on resource g.h5 beam beam1: boom") == "g.h5"
        @test GGA._sliderule_failed_resource("Failure on resource g.h5") === nothing

        # Test: both consumers accept an already-parsed record stream, so a response is walked once
        # instead of three times. Re-walking a 20 MB response for each consumer measured 5x the
        # response size in allocations, over tens of thousands of requests.
        name = "out.feather"
        file = mock_feather(DataFrame(a=[1, 2, 3]))
        raw = vcat(mock_meta(name, length(file)), mock_data(name, file),
            mock_except("Failure on resource A.h5 beam beam3: boom"))
        records = GGA.sliderule_records(raw)
        @test GGA.sliderule_failed_resources(records) == GGA.sliderule_failed_resources(raw)
        @test GGA.sliderule_arrow(records; filename=name) == GGA.sliderule_arrow(raw; filename=name)
    end

    @testset "GEDI request parameters" begin
        extent = Extent(X=(86.0, 88.0), Y=(28.0, 30.0))
        parms = GGA.sliderule_gedi_parms(extent; granules=["a.h5", "b.h5"],
            t0=DateTime(2024, 5, 1), t1=DateTime(2025, 7, 1), filename="g.feather")

        @test parms["asset"] == "gedil2a"
        @test parms["output"]["format"] == "feather"
        @test parms["resources"] == ["a.h5", "b.h5"]
        @test parms["t0"] == "2024-05-01T00:00:00Z"

        # Test: degrade and surface are pushed server-side -- the archive's own filter drops those
        # shots anyway, so filtering there saves transferring them.
        @test parms["degrade_filter"] === true
        @test parms["surface_filter"] === true

        # Test: l2_quality_filter is NOT set. The archive keeps `quality` as a column (13% of its
        # rows are false), so filtering on it here would discard data the HDF5 path retains.
        @test !haskey(parms, "l2_quality_filter")

        # Test: every field the schema and the quality filter need is requested, and rx_algrunflag is
        # not -- see `_gedi_l3_filter` for why that term is deliberately absent.
        for field in ("elevation_bin0_error", "energy_total", "num_detectedmodes",
            "digital_elevation_model", "rx_assess/quality_flag", "geolocation/stale_return_flag",
            "rx_assess/rx_maxamp", "rx_assess/sd_corrected", "selected_algorithm")
            @test field in parms["anc_fields"]
        end
        for a in 1:6
            @test "rx_processing_a$(a)/zcross" in parms["anc_fields"]
            @test "rx_processing_a$(a)/toploc" in parms["anc_fields"]
        end
        @test !any(f -> occursin("algrunflag", f), parms["anc_fields"])
        @test length(parms["anc_fields"]) == 21

        poly = parms["poly"]
        @test length(poly) == 5 && poly[1] == poly[end]
        @test extrema(p["lon"] for p in poly) == (86.0, 88.0)
        @test !haskey(GGA.sliderule_gedi_parms(extent), "resources")
    end

    @testset "GEDI beam mapping" begin
        # SlideRule's `beam` value is the binary number the beam is named after. Strong beams were
        # confirmed against 9.8M rows of the existing archive: exactly these four.
        @test GGA._gedi_beam_name.([0, 1, 2, 3, 5, 6, 8, 11]) ==
              ["BEAM0000", "BEAM0001", "BEAM0010", "BEAM0011",
            "BEAM0101", "BEAM0110", "BEAM1000", "BEAM1011"]
        # values that are not GEDI beams yield "", not a plausible-looking name
        @test GGA._gedi_beam_name.([4, 7, 9, 10, 12, 16, -1]) == fill("", 7)
        @test GGA.GEDI_STRONG_BEAMS == (0x05, 0x06, 0x08, 0x0b)
    end

    @testset "GEDI granule id from shot_number" begin
        # GEDI02_A_<datetime>_O<orbit>_<granule>_T<track>_...
        @test GGA._gedi_key("GEDI02_A_2019094180555_O01753_03_T01683_02_003_01_V002.h5") == (1753, 3, 1683)
        # v003 names have the same field layout, so the key is version-independent
        @test GGA._gedi_key("GEDI02_A_2025190223019_O37238_01_T09057_02_004_02_V003.h5") == (37238, 1, 9057)
        @test GGA._gedi_key("/dir/GEDI02_A_2022175044201_O19995_04_T07648_02_003_03_V002.h5") == (19995, 4, 7648)
        @test GGA._gedi_key("not-a-granule.h5") === nothing
        @test GGA._gedi_key("GEDI02_A_2019094180555_X01753_03_T01683_02_003_01_V002.h5") === nothing

        granules = ["GEDI02_A_2019094180555_O01753_03_T01683_02_003_01_V002.h5",
            "GEDI02_A_2022175044201_O19995_04_T07648_02_003_03_V002.h5"]

        # shot_number = orbit*10^13 + beam*10^11 + reserved*10^9 + granule*10^8 + shot_index.
        # The layout is right-anchored, so a 5-digit orbit gives an 18-digit shot number: both widths
        # must decode, which is the whole reason for dividing rather than slicing the digits.
        df = DataFrame(
            shot_number=UInt64[1753_00_00_3_00518838, 19995_00_00_4_00000017, 1753_00_00_9_00000001],
            track=UInt16[1683, 7648, 1683])
        @test ndigits(df.shot_number[1]) == 17
        @test ndigits(df.shot_number[2]) == 18

        ids = GGA._gedi_granule_ids(df, granules)
        @test ids[1] == granules[1]
        @test ids[2] == granules[2]
        # Test: the granule number matters. Row 3 shares orbit and track with row 1 and differs only
        # in granule number, which is exactly the case `(orbit, track)` alone cannot separate -- it
        # collides for 11,785 of the archive's granules.
        @test ids[3] == ""
    end

    @testset "GEDI schema mapping" begin
        granules = ["GEDI02_A_2019094180555_O01753_03_T01683_02_003_01_V002.h5"]

        # One row per outcome: passes everything, then one killed by each filter term.
        raw = gedi_raw_frame()
        out = GGA.sliderule2archive_gedi(raw, granules)

        @test names(out) == names(GGA.sliderule_gedi_empty_table())
        @test eltype.(eachcol(out)) == eltype.(eachcol(GGA.sliderule_gedi_empty_table()))

        # only the first row survives the L3 filter
        @test nrow(out) == 1
        @test out.height[1] == 100.0f0
        @test out.height_error[1] == 0.5f0
        @test out.intensity[1] == 400.0f0
        @test out.nmodes[1] == 0x03
        @test out.sun_angle[1] == 12.0f0
        @test out.track[1] == "BEAM0101"
        @test out.strong_beam[1] == true
        @test out.classification[1] == "ground"
        @test out.id[1] == granules[1]

        # Test: quality comes from the flags bitmask and is NOT inverted, unlike ATL06 -- GEDI's
        # quality_flag is already 1 for good.
        @test out.quality[1] == true
        @test out.surface[1] == true

        # Test: timestamps land on the archive's GPS timescale, 18 s ahead of SlideRule's UTC
        @test out.datetime[1] == DateTime(2019, 1, 1, 0, 0, 18)
        @test GGA.sliderule2archive_gedi(raw, granules; gps_time=false).datetime[1] == DateTime(2019, 1, 1)

        # Test: the DEM's -999999 fill becomes NaN. 11.5% of the archive is NaN here, so an all-NaN
        # column would not look wrong -- which is why the field's absence throws instead.
        @test isnan(GGA.sliderule2archive_gedi(gedi_raw_frame(dem=-999999.0f0), granules).height_reference[1])
        @test GGA.sliderule2archive_gedi(gedi_raw_frame(dem=42.0f0), granules).height_reference[1] == 42.0f0

        # Test: an empty result still carries the archive's columns and types
        empty_out = GGA.sliderule2archive_gedi(DataFrame(), granules)
        @test names(empty_out) == names(GGA.sliderule_gedi_empty_table())
        @test isempty(empty_out)
    end

    @testset "GEDI quality filter" begin
        # The existing archive was built through SpaceLiDAR.points at its default filtered=true, so
        # every point in it passed this filter. Appending unfiltered data would make the new points
        # systematically noisier than the old with nothing downstream to reveal it.
        keep = GGA._gedi_l3_filter(gedi_raw_frame())
        @test keep == [true, false, false, false, false, false, false]

        # Test: each term in isolation, so a term that stops working is not masked by another
        @test GGA._gedi_l3_filter(gedi_raw_frame(n=1))[1]
        @test !GGA._gedi_l3_filter(gedi_raw_frame(n=1, rx_quality=0x00))[1]
        @test !GGA._gedi_l3_filter(gedi_raw_frame(n=1, stale=0x01))[1]
        @test !GGA._gedi_l3_filter(gedi_raw_frame(n=1, maxamp=7.0f0, sd=1.0f0))[1]
        @test GGA._gedi_l3_filter(gedi_raw_frame(n=1, maxamp=8.0f0, sd=1.0f0))[1]
        @test !GGA._gedi_l3_filter(gedi_raw_frame(n=1, zcross=0.0f0))[1]
        @test !GGA._gedi_l3_filter(gedi_raw_frame(n=1, toploc=-1.0f0))[1]
        # surface and degrade come from the flags byte; also enforced server-side
        @test !GGA._gedi_l3_filter(gedi_raw_frame(n=1, flags=GGA.GEDI_L2_QUALITY_FLAG_MASK))[1]
        @test !GGA._gedi_l3_filter(gedi_raw_frame(n=1,
            flags=GGA.GEDI_SURFACE_FLAG_MASK | GGA.GEDI_L2_QUALITY_FLAG_MASK | GGA.GEDI_DEGRADE_FLAG_MASK))[1]

        # Test: zcross/toploc are read from the group named by selected_algorithm, not a fixed one.
        # Reading the wrong group is silent -- the values are plausible floats either way.
        for a in 1:6
            good = gedi_raw_frame(n=1, algorithm=UInt8(a), zcross=0.0f0)
            good[!, Symbol("rx_processing_a$(a)/zcross")] = Float32[5.0]
            @test GGA._gedi_l3_filter(good)[1]
            bad = gedi_raw_frame(n=1, algorithm=UInt8(a))
            bad[!, Symbol("rx_processing_a$(a)/zcross")] = Float32[0.0]
            @test !GGA._gedi_l3_filter(bad)[1]
        end
        # an algorithm outside 1:6 has no group to read, so the row is dropped
        @test !GGA._gedi_l3_filter(gedi_raw_frame(n=1, algorithm=0x07))[1]

        # Test: a missing filter input throws rather than passing rows through. Silently keeping
        # everything would produce an archive that is inhomogeneous in a way no plot would show.
        for field in ("rx_assess/quality_flag", "geolocation/stale_return_flag",
            "rx_assess/rx_maxamp", "selected_algorithm", "rx_processing_a1/zcross")
            missing_field = select(gedi_raw_frame(n=1), Not(Symbol(field)))
            @test_throws GGA.NonRetryable GGA._gedi_l3_filter(missing_field)
        end
        @test_throws GGA.NonRetryable GGA.sliderule2archive_gedi(
            select(gedi_raw_frame(n=1), Not(:digital_elevation_model)),
            ["GEDI02_A_2019094180555_O01753_03_T01683_02_003_01_V002.h5"])
    end

    @testset "GEDI query with injected transport" begin
        extent = Extent(X=(86.0, 88.0), Y=(28.0, 30.0))
        granules = ["GEDI02_A_2019094180555_O01753_03_T01683_02_003_01_V002.h5"]

        requests = String[]
        function poster(url, headers, body)
            push!(requests, body)
            file = mock_feather(gedi_raw_frame(n=1))
            return (; status=200, body=vcat(mock_meta("gga_gedi.feather", length(file)),
                mock_data("gga_gedi.feather", file)))
        end

        df, failed, partial, unreadable, retried = GGA.sliderule_gedi(extent; granules, poster)
        @test nrow(df) == 1
        @test df.id == granules
        @test isempty(failed)
        @test occursin("gedil2a", requests[1])

        # Test: resources are chunked. GEDI needs a smaller chunk than ATL06 -- eight beams and 21
        # ancillary datasets per granule is a lot of server-side reading to lose to one failure.
        @test GGA.SLIDERULE_GEDI_RESOURCE_CHUNK < GGA.SLIDERULE_RESOURCE_CHUNK
        many = ["GEDI02_A_20190941805$(lpad(i, 2, '0'))_O0$(lpad(i, 4, '0'))_03_T01683_02_003_01_V002.h5"
                for i in 1:(GGA.SLIDERULE_GEDI_RESOURCE_CHUNK+3)]
        empty!(requests)
        GGA.sliderule_gedi(extent; granules=many, poster)
        @test length(requests) == 2
    end

    @testset "a granule that yields nothing is not recorded" begin
        # A granule whose read broke and produced no points must be left out of the file entirely, so
        # the incremental rule asks for it again. Writing it -- with or without a placeholder -- would
        # lose its data permanently, because the id being present is what stops a re-request.
        mktempdir() do dir
            extent = Extent(X=(86.0, 88.0), Y=(28.0, 30.0))
            good = "GEDI02_A_2019094180555_O01753_03_T01683_02_003_01_V002.h5"
            bad = "GEDI02_A_2022175044201_O19995_04_T07648_02_003_03_V002.h5"
            geotiles = DataFrame(id=["lat[+28+30]lon[+086+088]"], extent=[extent],
                granules=[[(id=g, url=g) for g in [good, bad]]])

            asked = Vector{String}[]
            function poster(url, headers, body)
                resources = String.(JSON.parse(body)["parms"]["resources"])
                push!(asked, resources)
                # `bad` never returns a row of its own, only a failure
                rows = good in resources ? gedi_raw_frame(n=1) : gedi_raw_frame(n=0)
                body = if isempty(rows)
                    UInt8[]
                else
                    file = mock_feather(rows)
                    vcat(mock_meta("gga_gedi.feather", length(file)),
                        mock_data("gga_gedi.feather", file))
                end
                return (; status=200, body=vcat(body,
                    mock_except("Failure on resource $(bad) beam beam5: H5Coro::Future read failure")))
            end

            GGA.geotile_build_sliderule(geotiles, dir; mission=:gedi, poster, ntasks=1)
            built = DataFrame(Arrow.Table(joinpath(dir, "lat[+28+30]lon[+086+088].arrow")))

            # the good granule is recorded; the failed one appears nowhere, not even as a placeholder
            @test Set(built.id) == Set([good])
            @test !(bad in built.id)

            # Test: the failure is retried in place, and only the broken granule is asked for again
            @test length(asked) == GGA.SLIDERULE_BEAM_ATTEMPTS
            @test asked[1] == [good, bad]
            @test all(==([bad]), asked[2:end])

            # Test: and it is still outstanding on the next pass
            empty!(asked)
            GGA.geotile_build_sliderule(geotiles, dir; mission=:gedi, poster, ntasks=1)
            @test first(asked) == [bad]
        end
    end

    @testset "a retried granule is recorded once" begin
        # The retry must not double-count: the first attempt returns the granule's readable beams and
        # reports the failure, the second returns all of them. Keeping both would duplicate points.
        extent = Extent(X=(86.0, 88.0), Y=(28.0, 30.0))
        granule = "GEDI02_A_2019094180555_O01753_03_T01683_02_003_01_V002.h5"

        attempt = Ref(0)
        function poster(url, headers, body)
            attempt[] += 1
            file = mock_feather(gedi_raw_frame(n=1))
            body = vcat(mock_meta("gga_gedi.feather", length(file)),
                mock_data("gga_gedi.feather", file))
            # fail on the first attempt only
            attempt[] == 1 && (body = vcat(body,
                mock_except("Failure on resource $(granule) beam beam5: H5Coro::Future read failure")))
            return (; status=200, body)
        end

        df, failed, partial, unreadable, retried = GGA.sliderule_gedi(extent; granules=[granule], poster)
        @test attempt[] == 2
        @test isempty(failed)
        @test nrow(df) == 1          # not 2 -- the failed attempt's rows were discarded
        @test df.id == [granule]

        # Test: a granule that yields nothing at all is a failure after the attempts are used up, so
        # it is left unrecorded and asked for again next pass
        always_fails(url, headers, body) = (; status=200,
            body=mock_except("Failure on resource $(granule) beam beam5: read failure"))
        df2, failed2, partial2, retried2 = GGA.sliderule_gedi(extent; granules=[granule], poster=always_fails)
        @test failed2 == Set([granule])
        @test isempty(partial2)
        @test isempty(df2)
        @test names(df2) == names(GGA.sliderule_gedi_empty_table())
    end

    @testset "unreadable versus empty-in-box" begin
        # The distinction that decides whether a build can converge. A granule outside the requested
        # box reports every readable beam with zero rows; if two of its beams are also unreadable, the
        # result looks identical to "could not read this granule" if you only count rows. Calling that
        # a failure means it is never recorded and every later pass asks again.
        g = "GEDI02_A_2024125115322_O30546_03_T00000_02_004_02_V002.h5"

        # six beams finished with zero rows, two failed -> the granule answered; not unreadable
        answered = vcat(
            [mock_except("Failure on resource $g beam beam$(b): H5Coro::Future read failure on BEAM/lat_lowestmode")
             for b in (8, 11)]...,
            [mock_except("request <X> on $g generated dataframe [beam$(b)] with 0 rows and 30 columns")
             for b in (0, 1, 2, 3, 5, 6, 8, 11)]...)
        @test GGA.sliderule_beam_verdicts(GGA.sliderule_records(answered))[g] === :reported

        # every beam that reported also failed -> genuinely unreadable
        dead = vcat(
            [mock_except("Failure on resource $g beam beam$(b): H5Coro::Future read failure on BEAM/quality_flag")
             for b in (0, 1, 2, 3, 5, 6, 8, 11)]...,
            [mock_except("request <X> on $g generated dataframe [beam$(b)] with 0 rows and 9 columns")
             for b in (0, 1, 2, 3, 5, 6, 8, 11)]...)
        @test GGA.sliderule_beam_verdicts(GGA.sliderule_records(dead))[g] === :unreadable

        # Test: failures with no dataframe lines at all establish nothing about the granule. This is
        # what a service-level error looks like, and calling it unreadable would blank readable granules
        # out of the archive on the strength of an outage.
        @test GGA.sliderule_beam_verdicts(GGA.sliderule_records(
            mock_except("Failure on resource $g beam beam0: boom")))[g] === :inconclusive

        # Test: a granule with no failures at all does not appear
        @test !haskey(GGA.sliderule_beam_verdicts(GGA.sliderule_records(
            mock_except("request <X> on $g generated dataframe [beam0] with 12 rows and 30 columns"))), g)

        # Test: end to end -- a granule whose readable beams are empty gets a placeholder and is not
        # asked for again, which is the whole point.
        mktempdir() do dir
            extent = Extent(X=(-74.0, -72.0), Y=(6.0, 8.0))
            geotiles = DataFrame(id=["lat[+06+08]lon[-074-072]"], extent=[extent],
                granules=[[(id=g, url=g)]])
            asked = Vector{String}[]
            function poster(url, headers, body)
                push!(asked, String.(JSON.parse(body)["parms"]["resources"]))
                return (; status=200, body=answered)
            end

            GGA.geotile_build_sliderule(geotiles, dir; mission=:gedi, poster, ntasks=1)
            built = DataFrame(Arrow.Table(joinpath(dir, "lat[+06+08]lon[-074-072].arrow")))
            @test built.id == [g]          # recorded, as a placeholder
            @test all(isnan, built.height)
            # empty in this box, not unreadable -- the two placeholder kinds stay distinguishable
            @test built.track == [GGA.PLACEHOLDER_TRACK_EMPTY]

            empty!(asked)
            GGA.geotile_build_sliderule(geotiles, dir; mission=:gedi, poster, ntasks=1)
            @test isempty(asked)           # and never requested again
        end
    end

    @testset "an unreadable granule is recorded and stops being requested" begin
        # The convergence property. A granule nothing can be read from fails identically on every pass,
        # so leaving it unrecorded means requesting it forever and the build never reports itself done.
        # Recording it as a placeholder marked unreadable retires it while keeping it findable.
        g = "GEDI02_A_2022132212910_O19339_03_T01650_02_003_02_V002.h5"
        dead = vcat(
            [mock_except("Failure on resource $g beam beam$(b): H5Coro::Future read failure on BEAM/quality_flag")
             for b in (0, 1, 2, 3, 5, 6, 8, 11)]...,
            [mock_except("request <X> on $g generated dataframe [beam$(b)] with 0 rows and 9 columns")
             for b in (0, 1, 2, 3, 5, 6, 8, 11)]...)

        mktempdir() do dir
            extent = Extent(X=(86.0, 88.0), Y=(28.0, 30.0))
            gt = "lat[+28+30]lon[+086+088]"
            geotiles = DataFrame(id=[gt], extent=[extent], granules=[[(id=g, url=g)]])
            asked = Vector{String}[]
            function poster(url, headers, body)
                push!(asked, String.(JSON.parse(body)["parms"]["resources"]))
                return (; status=200, body=dead)
            end

            GGA.geotile_build_sliderule(geotiles, dir; mission=:gedi, poster, ntasks=1)
            built = DataFrame(Arrow.Table(joinpath(dir, gt * ".arrow")))
            @test built.id == [g]
            @test built.track == [GGA.PLACEHOLDER_TRACK_UNREADABLE]
            @test all(isnan, built.height)

            empty!(asked)
            GGA.geotile_build_sliderule(geotiles, dir; mission=:gedi, poster, ntasks=1)
            @test isempty(asked)

            # Test: logged under its own kind, so the frozen granules can be listed
            log = readlines(joinpath(dir, GGA.SLIDERULE_ANOMALY_FILE))
            @test any(l -> occursin("\tunreadable\t" * g, l), log)
        end
    end

    @testset "an inconclusive failure is left for the next pass" begin
        # Beam failures with no per-beam accounting are what a service outage looks like. Recording such
        # a granule would freeze real data out of the archive, so it must stay unrecorded and be asked
        # for again -- the opposite of the unreadable case above.
        g = "GEDI02_A_2024125115322_O30546_03_T00000_02_004_02_V002.h5"
        blind = mock_except("Failure on resource $g beam beam0: H5Coro::Future read failure")

        mktempdir() do dir
            extent = Extent(X=(86.0, 88.0), Y=(28.0, 30.0))
            gt = "lat[+28+30]lon[+086+088]"
            geotiles = DataFrame(id=[gt], extent=[extent], granules=[[(id=g, url=g)]])
            asked = Vector{String}[]
            function poster(url, headers, body)
                push!(asked, String.(JSON.parse(body)["parms"]["resources"]))
                return (; status=200, body=blind)
            end

            GGA.geotile_build_sliderule(geotiles, dir; mission=:gedi, poster, ntasks=1)
            @test !isfile(joinpath(dir, gt * ".arrow"))   # nothing recorded at all

            empty!(asked)
            GGA.geotile_build_sliderule(geotiles, dir; mission=:gedi, poster, ntasks=1)
            @test g in reduce(vcat, asked)                # and asked for again
        end
    end

    @testset "anomalies are recorded by granule id" begin
        # Counts are not enough: a granule recorded without all its beams is indistinguishable in the
        # archive from a complete one, so nothing downstream can find it and the incremental rule will
        # never re-request it. The ids are the only way to audit or repair afterwards.
        mktempdir() do dir
            GGA.record_sliderule_anomalies(dir, "lat[+28+30]lon[+086+088]";
                failed=Set(["bad2.h5", "bad1.h5"]), partial=Set(["thin.h5"]),
                unreadable=Set(["dead.h5"]))
            GGA.record_sliderule_anomalies(dir, "lat[+30+32]lon[+086+088]"; partial=Set(["x.h5"]))

            path = joinpath(dir, GGA.SLIDERULE_ANOMALY_FILE)
            lines = readlines(path)
            @test lines[1] == "timestamp\tgeotile\tkind\tgranule"
            body = lines[2:end]
            @test length(body) == 5
            @test count(l -> occursin("\tfailed\t", l), body) == 2
            @test count(l -> occursin("\tpartial\t", l), body) == 2
            # the three kinds are distinguishable, which is what makes the log auditable
            @test count(l -> occursin("\tunreadable\t", l), body) == 1
            # ids present and attributed to the right geotile
            @test any(l -> occursin("lat[+28+30]lon[+086+088]\tpartial\tthin.h5", l), body)
            @test any(l -> occursin("lat[+28+30]lon[+086+088]\tunreadable\tdead.h5", l), body)
            @test any(l -> occursin("lat[+30+32]lon[+086+088]\tpartial\tx.h5", l), body)
            # failed ids sorted, so a diff between passes is readable
            @test findfirst(l -> occursin("bad1.h5", l), body) < findfirst(l -> occursin("bad2.h5", l), body)

            # Test: nothing to report writes nothing at all
            other = mktempdir()
            GGA.record_sliderule_anomalies(other, "g")
            @test !isfile(joinpath(other, GGA.SLIDERULE_ANOMALY_FILE))
        end
    end

    @testset "a permanently unreadable beam is not a failed granule" begin
        # Some granules have beams that never read, however many attempts are made -- release-004
        # granules whose BEAM1000 and BEAM1011 always fail are the case that surfaced this, at ~1.8% of
        # granules over a calibration run. Treating those as transient would reject them on every pass
        # forever and silently discard the beams that did read.
        extent = Extent(X=(86.0, 88.0), Y=(28.0, 30.0))
        granule = "GEDI02_A_2024123163515_O30518_02_T00000_02_004_02_V002.h5"

        attempts = Ref(0)
        function poster(url, headers, body)
            attempts[] += 1
            # six beams return data, two fail, every single time
            raw = gedi_raw_frame(n=1)
            raw.shot_number = UInt64[30518_00_00_2_00000042]
            raw.track = UInt16[0]
            file = mock_feather(raw)
            return (; status=200, body=vcat(
                mock_meta("gga_gedi.feather", length(file)),
                mock_data("gga_gedi.feather", file),
                mock_except("Failure on resource $(granule) beam beam8: H5Coro::Future read failure on BEAM1000/lat_lowestmode"),
                mock_except("Failure on resource $(granule) beam beam11: H5Coro::Future read failure on BEAM1011/lat_lowestmode"),
                # the six beams that did read report too -- a real response always carries these, and
                # they are what distinguishes "some beams failed" from "nothing could be read"
                [mock_except("request <X> on $(granule) generated dataframe [beam$(b)] with 1200 rows and 30 columns")
                 for b in (0, 1, 2, 3, 5, 6)]...,
                [mock_except("request <X> on $(granule) generated dataframe [beam$(b)] with 0 rows and 9 columns")
                 for b in (8, 11)]...))
        end

        df, failed, partial, unreadable, retried = GGA.sliderule_gedi(extent; granules=[granule], poster)

        # Test: retried to exhaustion, then accepted rather than discarded
        @test attempts[] == GGA.SLIDERULE_BEAM_ATTEMPTS
        @test isempty(failed)
        @test partial == Set([granule])
        @test nrow(df) == 1
        @test df.id == [granule]

        # Test: a track number of zero, which release 004 uses, still resolves the granule id
        @test GGA._gedi_key(granule) == (30518, 2, 0)

        # Test: and the geotile records it, so it is not requested again forever
        mktempdir() do dir
            geotiles = DataFrame(id=["lat[+28+30]lon[+086+088]"], extent=[extent],
                granules=[[(id=granule, url=granule)]])
            GGA.geotile_build_sliderule(geotiles, dir; mission=:gedi, poster, ntasks=1)
            built = DataFrame(Arrow.Table(joinpath(dir, "lat[+28+30]lon[+086+088].arrow")))
            @test built.id == [granule]

            asked = Ref(0)
            counting(url, headers, body) = (asked[] += 1; poster(url, headers, body))
            GGA.geotile_build_sliderule(geotiles, dir; mission=:gedi, poster=counting, ntasks=1)
            @test asked[] == 0
        end
    end
end
