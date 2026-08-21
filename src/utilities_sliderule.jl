"""
SlideRule ingest for ICESat-2 ATL06.

SlideRule (https://slideruleearth.io) subsets granules server-side in the same AWS region the NASA
archive lives in and returns only the points inside a region of interest. For a geotile archive that
turns an update into work proportional to the *new* data, instead of downloading whole 40-100 MB
granules to keep the few thousand points that fall inside a 2 degree tile.

This file is the alternative source for the same per-geotile `.arrow` archive that
[`geotile_build`](@ref) writes from local HDF5. It produces identical columns, keys incremental
updates on granule id exactly as the HDF5 path does, and writes through [`atomic_write`](@ref), so the
two sources can populate one archive interchangeably.

Verified against the HDF5 path on `lat[-80-78]lon[+166+168]`: across 47,930 coordinate-matched points,
`height`, `height_error`, `height_reference`, `quality`, `track`, `strong_beam`, `detector_id` and the
granule `id` were all identical, and `datetime` agreed to the millisecond once
[`SLIDERULE_GPS_UTC_OFFSET`](@ref) was applied. The one difference is coverage: SlideRule returned
99.5% of the in-box points and never a point the HDF5 path lacked -- see
[`SLIDERULE_POLY_BUFFER`](@ref).

# Which endpoint

`atl06x` -- the "x-series" ATL06 Dataframe endpoint. It subsets the **standard** ATL06 product, which
is what this project's calibration is built on. Two near-misses worth naming:

  - `atl06p` runs SlideRule's own surface-fit algorithm over ATL03 photons to make ATL06-SR. That is a
    different product with different corrections, not a substitute.
  - `atl06sp` also subsets standard ATL06, but does not support `atl06_fields`, which is how
    `height_reference` (`land_ice_segments/dem/dem_h`) is obtained.

# Wire format

Responses are SlideRule's native record stream, regardless of `Accept` header. Each record is

    | Int16 version | Int16 type_size | Int32 data_size |   <- 8 byte header, big endian
    | type_size bytes of NUL-terminated type string   |
    | data_size bytes of payload                      |   <- fields are native (little) endian

Requesting `output = {"format": "feather", ...}` makes the server emit the result as an Arrow IPC
file, delivered as `arrowrec.meta` (filename + total size), one or more `arrowrec.data` chunks, and
`arrowrec.eof` (checksum). Concatenating the `data` payloads -- each prefixed by a
`SLIDERULE_FILENAME_LEN` byte filename field -- reconstructs a file `Arrow.jl` reads directly, so no
parquet dependency is needed.

`exceptrec` records are interleaved. Most are the expected "no data for this granule" notices: CMR
adds geolocation margin to ICESat-2 granule polygons, so a search returns granules whose data does not
actually reach the region. Those are informational, not failures.

# Cluster and credentials

Defaults to the public cluster, which needs no credentials but applies the strictest rate limits. A
GitHub personal access token authenticates the request -- which relaxes rate limits on the public
cluster too, not only on private ones -- and is read from `SLIDERULE_GITHUB_TOKEN`, or failing that
from the file named by `SLIDERULE_PAT_FILE` (default `~/.sliderule_pat`, see
[`sliderule_token`](@ref)). Set `SLIDERULE_ORGANIZATION` to target a private cluster, which
additionally requires membership of the SlideRuleEarth GitHub organization and a cluster deployed
through `provisioner.slideruleearth.io` with a time-to-live.
"""

const SLIDERULE_PUBLIC_HOST = "sliderule.slideruleearth.io"
const SLIDERULE_LOGIN_HOST = "login.slideruleearth.io"
const SLIDERULE_ATL06_ENDPOINT = "atl06x"

# Length of the fixed `filename` field that prefixes every arrowrec payload; `FILE_NAME_MAX_LEN` in
# SlideRule's OutputLib.h. Read out of the arrowrec.meta record and cross-checked against the
# filename we asked for, so a server-side change in this layout fails loudly rather than silently
# shifting every byte of the assembled file.
const SLIDERULE_FILENAME_LEN = 128

# exceptrec payload: Int32 code, UInt32 level, then NUL-terminated text (EventLib.h alert_t).
const SLIDERULE_ALERT_TEXT_OFFSET = 8

# ATL06 fill value for the height fields. SpaceAltimetry's HDF5 reader maps these to NaN before use
# and the binning code assumes NaN, so do the same here rather than letting 3.4e38 into the archive.
const ATL06_FILL_VALUE = 3.4028235f38

# GPS minus UTC, in seconds.
#
# ATL06 stores `delta_time` against the ATLAS SDP GPS epoch. SlideRule converts to true UTC;
# SpaceAltimetry adds the epoch plus `gps_offset` without subtracting leap seconds, so every
# `datetime` already in the archive is GPS time -- 18 s ahead of UTC. Measured directly by
# coordinate-matching a geotile built both ways: the offset is -18.000 s on every one of 45,466
# points.
#
# The archive's convention wins, so that a geotile fed from both sources carries one timescale. 18 s
# is nothing against the monthly bins downstream, but two conventions inside one file would be a trap.
# The count has been 18 since 2017 and so covers all of ICESat-2; a future leap second would apply only
# to data after it, and this constant would need to become time-dependent.
const SLIDERULE_GPS_UTC_OFFSET = Second(18)

# Granules per request. Each geotile can intersect several hundred granules; splitting keeps a single
# request from monopolising cluster nodes and bounds what has to be retried after a failure.
const SLIDERULE_RESOURCE_CHUNK = 100

# Ancillary ATL06 field backing the archive's `height_reference` column, and the name it comes back
# under. SlideRule labels ancillary columns with the exact string that was requested, so this is
# "dem/dem_h" and not "dem_h" -- checking for the latter silently yields an all-NaN column.
const SLIDERULE_ATL06_FIELDS = ["dem/dem_h"]
const SLIDERULE_DEM_COLUMNS = (Symbol("dem/dem_h"), :dem_h)

# Degrees by which the requested polygon is grown beyond the geotile. Zero by default -- see below.
#
# SlideRule returns slightly fewer points than reading the granule locally: on `lat[-80-78]lon[+166+168]`
# with 6 granules, 47,930 of the HDF5 path's 48,186 in-box points came back, the 256 missing ones lying
# in a 0.07 degree band against the tile's northern edge. The shortfall is one-directional -- zero
# points were returned that the HDF5 path did not also have -- so this loses a little edge data rather
# than inventing any.
#
# Growing the polygon and clipping locally makes it *worse*, not better (buffer=0.1 returned a 59,085
# point superset yet only 47,056 in-box, losing 1,130), so the mechanism is not a simple inset boundary
# and a buffer is not the fix. Left at zero as the measured best, and exposed as a keyword so the
# behaviour can be re-tested if SlideRule's subsetting changes.
const SLIDERULE_POLY_BUFFER = 0.0

# Cached bearer token: the JWT, the epoch seconds at which to renew it, and the `sub` claim, which is
# the identity the provisioner books dedicated capacity against.
const SLIDERULE_TOKEN = Ref{Union{Nothing,NamedTuple{(:token, :renew_at, :service),Tuple{String,Float64,String}}}}(nothing)

# Default location of the personal access token, used when `SLIDERULE_GITHUB_TOKEN` is unset. A file
# keeps the secret out of shell history, `env` dumps and the process table, which matters for a job
# that runs for days; `SLIDERULE_PAT_FILE` overrides it.
const SLIDERULE_PAT_FILE = "~/.sliderule_pat"

"""
    sliderule_pat() -> String

The GitHub personal access token, or `""` when none is configured.

`SLIDERULE_GITHUB_TOKEN` wins when set. Otherwise the first non-empty line of the file named by
`SLIDERULE_PAT_FILE` (default [`SLIDERULE_PAT_FILE`](@ref)) is used, so a long-running or scheduled
build does not need the secret in its environment. A token file readable by anyone but its owner
warns rather than being ignored: refusing it would silently downgrade a build to anonymous, which is
the more expensive failure.
"""
function sliderule_pat()
    token = get(ENV, "SLIDERULE_GITHUB_TOKEN", "")
    isempty(token) || return strip(token)

    path = expanduser(get(ENV, "SLIDERULE_PAT_FILE", SLIDERULE_PAT_FILE))
    isfile(path) || return ""

    # 0o077 covers group and other; the token is a credential, not a config value.
    mode = filemode(path)
    if mode & 0o077 != 0
        @warn "SlideRule token file is readable beyond its owner; chmod 600 it" path mode = string(mode & 0o777, base=8) maxlog = 1
    end

    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) || return stripped
    end
    return ""
end

"""
    sliderule_host() -> String

Host to send requests to: `<SLIDERULE_ORGANIZATION>.slideruleearth.io` when that variable is set,
otherwise the public cluster.
"""
function sliderule_host()
    org = get(ENV, "SLIDERULE_ORGANIZATION", "")
    return isempty(org) ? SLIDERULE_PUBLIC_HOST : "$(org).slideruleearth.io"
end

"""
    sliderule_token(; force=false) -> Union{Nothing,String}

Bearer token for authenticated requests, or `nothing` when running anonymously.

Exchanges the GitHub personal access token from [`sliderule_pat`](@ref) for a JWT at
`login.slideruleearth.io/auth/github/pat`, and caches it until half its lifetime has elapsed. Returns
`nothing` when no token is configured, which is a valid (rate-limited) mode on the public cluster.

A failed exchange warns and falls back to anonymous rather than aborting: an expired PAT should
degrade a request, not kill a multi-hour build.
"""
function sliderule_token(; force=false)
    # Cache first, then the PAT. A live token is a live token: re-deriving whether one *could* be
    # minted before using one we already hold would drop a valid session if the token file moved
    # mid-run, and would make the cached identity unreachable to `sliderule_service`.
    cached = SLIDERULE_TOKEN[]
    if !force && !isnothing(cached) && time() < cached.renew_at
        return cached.token
    end

    pat = sliderule_pat()
    isempty(pat) && return nothing

    try
        response = HTTP.post("https://$(SLIDERULE_LOGIN_HOST)/auth/github/pat",
            ["Content-Type" => "application/json"], JSON.json(Dict("pat" => pat));
            status_exception=false)
        result = JSON.parse(String(response.body))
        if get(result, "status", "") != "success"
            @warn "SlideRule PAT login failed; continuing anonymously" status = get(result, "status", "unknown")
            return nothing
        end
        token = result["token"]
        metadata = get(result, "metadata", Dict())
        # Renew at half the remaining lifetime, matching the Python client, so a long build never
        # runs a request with a token that expires mid-flight.
        exp = get(metadata, "exp", time() + 3600)
        SLIDERULE_TOKEN[] = (token=token, renew_at=(exp - time()) / 2 + time(),
            service=String(get(metadata, "sub", "")))
        return token
    catch e
        @warn "SlideRule PAT login errored; continuing anonymously" exception = e
        return nothing
    end
end

"""
    sliderule_service() -> Union{Nothing,String}

The identity dedicated capacity is booked against -- the `sub` claim of the bearer token.

`nothing` when running anonymously. This is what the Python client calls the "service" in
`user_service` mode, and it is neither the cluster name nor the hostname: dedicated nodes sit behind
the ordinary public hostname and the gateway routes to them based on who is asking.
"""
function sliderule_service()
    isnothing(sliderule_token()) && return nothing
    cached = SLIDERULE_TOKEN[]
    (isnothing(cached) || isempty(cached.service)) && return nothing
    return cached.service
end

# ---------------------------------------------------------------------------------------------------
# Dedicated capacity
#
# A run of this size -- tens of thousands of granule reads -- does not belong on the shared public
# cluster: it is slower for us because of the anonymous rate limits, and it crowds out everyone else.
# With dedicated capacity provisioned against our identity the same requests go to nodes nobody else
# is queueing behind.
#
# This does not need SlideRule's Python client. `user_service` mode keeps talking to the ordinary
# public hostname with the bearer token we already mint from the PAT; the only additions are three
# JSON posts, mirroring `Session.scaleout` in `clients/python/sliderule/session.py`:
#
#   deploy   POST provisioner.slideruleearth.io/deploy/<service>
#   extend   POST provisioner.slideruleearth.io/extend/<service>
#   capacity POST <host>/discovery/status  {"service": <service>}  -> {"nodes": n}
#
# Clusters carry a time-to-live and shut themselves down when it lapses, so nothing has to be torn
# down by hand -- but a build outliving its TTL loses the nodes underneath it, hence
# [`sliderule_keepalive`](@ref).
# ---------------------------------------------------------------------------------------------------

const SLIDERULE_PROVISIONER_HOST = "provisioner.slideruleearth.io"

# Cluster name to provision against. Dedicated capacity is booked under the public cluster; the
# `service` path segment is what makes it ours.
const SLIDERULE_CLUSTER = "sliderule"

# Longest TTL the provisioner grants, in minutes.
const SLIDERULE_MAX_TTL = 720

# How long to wait for nodes to come up before giving up and using whatever is available.
const SLIDERULE_SCALEOUT_TIMEOUT = 600

"""
    sliderule_gateway(api, data; host=SLIDERULE_PROVISIONER_HOST, poster=HTTP.post) -> Dict

POST `data` to a SlideRule gateway service and parse the JSON reply.

Used for the provisioner and discovery APIs, which are plain JSON rather than the record stream
[`sliderule_post`](@ref) handles. A cold AWS API Gateway commonly answers the first request with 401,
so that one status is retried once before being treated as a real authorization failure.
"""
function sliderule_gateway(api, data; host=SLIDERULE_PROVISIONER_HOST, poster=HTTP.post)
    url = "https://$(host)/$(api)"
    headers = ["Content-Type" => "application/json"]
    token = sliderule_token()
    isnothing(token) || push!(headers, "Authorization" => "Bearer $(token)")
    body = JSON.json(data)

    response = poster(url, headers, body)
    if response.status == 401
        sleep(1)
        response = poster(url, headers, body)
    end
    if response.status != 200
        throw(NonRetryable("SlideRule $(api) returned HTTP $(response.status): " *
                           first(String(response.body), 300)))
    end
    return JSON.parse(String(response.body))
end

"""
    sliderule_capacity(; poster=HTTP.post) -> Int

Nodes currently serving our dedicated capacity, or 0 when there are none (or we are anonymous).
"""
function sliderule_capacity(; poster=HTTP.post)
    service = sliderule_service()
    isnothing(service) && return 0
    try
        result = sliderule_gateway("discovery/status", Dict("service" => service);
            host=sliderule_host(), poster)
        return Int(get(result, "nodes", 0))
    catch e
        # No cluster yet is the normal answer before the first deploy, and reads as an error here.
        @debug "SlideRule capacity check failed" exception = e
        return 0
    end
end

"""
    sliderule_stack_status(; poster=HTTP.post) -> Union{Nothing,String}

CloudFormation status of our dedicated cluster's stack, or `nothing` if there is none.

Capacity alone cannot tell a cluster that is still starting from one that failed and rolled back --
both report zero nodes -- so a build would otherwise wait out the whole scaleout timeout before
falling back, every time, with no indication of why.
"""
function sliderule_stack_status(; poster=HTTP.post)
    service = sliderule_service()
    isnothing(service) && return nothing
    try
        result = sliderule_gateway("status/$(service)", Dict("cluster" => SLIDERULE_CLUSTER);
            host=SLIDERULE_PROVISIONER_HOST, poster)
        response = get(result, "response", Dict())
        return response isa Dict ? get(response, "StackStatus", nothing) : nothing
    catch e
        @debug "SlideRule stack status check failed" exception = e
        return nothing
    end
end

"""
    sliderule_scaleout(; node_capacity=10, ttl=SLIDERULE_MAX_TTL, block=true, timeout=SLIDERULE_SCALEOUT_TIMEOUT, poster=HTTP.post) -> Int

Provision `node_capacity` dedicated nodes for `ttl` minutes and return the nodes available.

Deploys when short of `node_capacity`, extends the TTL when already at it, then (with `block`) waits
for the nodes to come up. Returns 0 and warns when running anonymously, which is not fatal: the
request path is identical either way, so the build still runs against the public cluster.
"""
function sliderule_scaleout(; node_capacity=10, ttl=SLIDERULE_MAX_TTL, block=true,
    timeout=SLIDERULE_SCALEOUT_TIMEOUT, poster=HTTP.post)

    service = sliderule_service()
    if isnothing(service)
        @warn "no SlideRule credentials, so no dedicated capacity; falling back to the shared " *
              "public cluster and its anonymous rate limits"
        return 0
    end

    ttl = min(ttl, SLIDERULE_MAX_TTL)
    available = sliderule_capacity(; poster)

    # A stack left in ROLLBACK_COMPLETE blocks every future deploy with AlreadyExistsException, and
    # clearing it needs the `destroy` endpoint, which an `affiliate` org role is not permitted to call.
    # So this is not something a retry or a longer wait can fix: say so and fall back, rather than
    # burning the scaleout timeout on a cluster that cannot come up.
    if available == 0
        status = sliderule_stack_status(; poster)
        if !isnothing(status) && occursin("ROLLBACK", status)
            @warn "SlideRule dedicated cluster is stuck in $(status) and cannot be redeployed; " *
                  "clearing it needs the provisioner's destroy endpoint, which requires a higher org " *
                  "role. Falling back to the shared cluster -- ask the SlideRule team to delete the " *
                  "stack." service
            return 0
        end
    end

    # Provisioning is an optimisation, never a precondition: the request path is identical with or
    # without dedicated nodes. So a provisioner that will not cooperate downgrades the run to the
    # shared cluster instead of aborting it -- otherwise a service-side fault costs hours of ingest
    # that would have completed fine, which is exactly what a stuck stack did here.
    try
        if available < node_capacity
            result = sliderule_gateway("deploy/$(service)",
                Dict("cluster" => SLIDERULE_CLUSTER, "is_public" => false,
                    "node_capacity" => node_capacity, "ttl" => ttl, "version" => "latest"); poster)
            haskey(result, "error") && error("deploy: $(get(result, "error_description", result))")
            printstyled("    -> requested $(node_capacity) dedicated node(s) for $(ttl) min\n"; color=:light_black)
        else
            result = sliderule_gateway("extend/$(service)",
                Dict("cluster" => SLIDERULE_CLUSTER, "ttl" => ttl); poster)
            haskey(result, "error") && error("extend: $(get(result, "error_description", result))")
            printstyled("    -> extended existing cluster to $(ttl) min\n"; color=:light_black)
        end
    catch e
        e isa InterruptException && rethrow(e)
        @warn "could not provision dedicated SlideRule capacity; continuing on the shared cluster" *
              " (this costs throughput, not correctness)" exception = (e, catch_backtrace())
        return 0
    end

    block || return available

    start = time()
    while available < node_capacity
        if time() - start > timeout
            @warn "SlideRule cluster did not reach the requested capacity in time; " *
                  "continuing with what is up" available node_capacity
            break
        end
        sleep(10)
        available = sliderule_capacity(; poster)

        # Watch the stack, not just the node count: a failed deploy also reports zero nodes, and
        # without this the wait runs to timeout before anyone learns the stack rolled back.
        if available == 0
            status = sliderule_stack_status(; poster)
            if !isnothing(status) && occursin("ROLLBACK", status)
                @warn "SlideRule cluster deployment rolled back ($(status)); falling back to the " *
                      "shared cluster" service
                return 0
            end
        end

        printstyled("    -> $(available)/$(node_capacity) node(s) up after $(round(Int, time() - start)) s\n";
            color=:light_black)
    end
    return available
end

"""
    sliderule_keepalive(f; node_capacity=10, ttl=SLIDERULE_MAX_TTL, every=1800, poster=HTTP.post)

Run `f()` with the dedicated cluster's TTL refreshed in the background.

A cluster shuts down when its TTL lapses, which mid-build means every request starts failing. Rather
than betting the run finishes inside one TTL, extend periodically for as long as `f` is running. The
extension is best-effort: a failed extend warns and is retried at the next tick, because losing the
refresh is not a reason to abandon work already done.
"""
function sliderule_keepalive(f; node_capacity=10, ttl=SLIDERULE_MAX_TTL, every=1800, poster=HTTP.post)
    nodes = sliderule_scaleout(; node_capacity, ttl, poster)
    nodes == 0 && return f()

    running = Ref(true)
    refresher = @async while running[]
        # short sleeps so the task notices `f` finishing promptly instead of holding the process open
        for _ in 1:(every÷5)
            running[] || break
            sleep(5)
        end
        running[] || break
        try
            sliderule_scaleout(; node_capacity, ttl, block=false, poster)
        catch e
            @warn "SlideRule TTL extension failed; will retry" exception = e
        end
    end

    try
        return f()
    finally
        running[] = false
        # let the task observe the flag rather than leaving it dangling
        try
            wait(refresher)
        catch
        end
    end
end

"""
    sliderule_post(endpoint, parms; host=sliderule_host(), poster=HTTP.post, max_attempts=5, backoff=1) -> Vector{UInt8}

POST `parms` to a SlideRule `endpoint` and return the raw response body.

`poster` is the function actually called as `poster(url, headers, body)`; it exists so the request
path can be exercised without network access. Transient failures retry through
[`with_retry`](@ref); an HTTP error status is a [`NonRetryable`](@ref) request problem and fails
immediately.
"""
function sliderule_post(endpoint, parms; host=sliderule_host(), poster=HTTP.post, max_attempts=5, backoff=1)
    url = "https://$(host)/source/$(endpoint)"
    headers = ["Content-Type" => "application/json"]
    token = sliderule_token()
    isnothing(token) || push!(headers, "Authorization" => "Bearer $(token)")
    body = JSON.json(Dict("parms" => parms))

    return with_retry("sliderule $(endpoint)"; max_attempts, backoff) do
        response = poster(url, headers, body)
        if response.status != 200
            throw(NonRetryable("SlideRule $(endpoint) returned HTTP $(response.status): " *
                               first(String(response.body), 300)))
        end
        return Vector{UInt8}(response.body)
    end
end

"""
    sliderule_records(bytes) -> Vector{Tuple{String,Vector{UInt8}}}

Split a SlideRule response into `(record type, payload)` pairs.

Walks the 8-byte big-endian header of each record. Throws if a header describes a record that runs
past the end of the buffer, which is what a truncated transfer looks like -- silently returning the
records that did parse would present partial data as complete.
"""
function sliderule_records(bytes::Vector{UInt8})
    records = Tuple{String,Vector{UInt8}}[]
    i = firstindex(bytes)
    stop = lastindex(bytes)

    while i + 7 <= stop
        type_size = Int(ntoh(reinterpret(Int16, bytes[i+2:i+3])[1]))
        data_size = Int(ntoh(reinterpret(Int32, bytes[i+4:i+7])[1]))
        if type_size <= 0 || data_size < 0
            throw(NonRetryable("SlideRule response has a malformed record header at byte $(i - firstindex(bytes)) " *
                               "(type_size=$type_size, data_size=$data_size)"))
        end
        record_end = i + 7 + type_size + data_size
        if record_end > stop
            throw(NonRetryable("SlideRule response is truncated: record at byte $(i - firstindex(bytes)) " *
                               "needs $(type_size + data_size) bytes, only $(stop - i - 7) remain"))
        end
        # type string is NUL-terminated, so drop the final byte
        rectype = String(bytes[i+8:i+6+type_size])
        payload = bytes[i+8+type_size:record_end]
        push!(records, (rectype, payload))
        i = record_end + 1
    end

    if i <= stop
        throw(NonRetryable("SlideRule response has $(stop - i + 1) trailing bytes that are too few for a record header"))
    end
    return records
end

"""
    sliderule_alert(payload) -> String

Text of an `exceptrec` record.
"""
function sliderule_alert(payload::Vector{UInt8})
    length(payload) <= SLIDERULE_ALERT_TEXT_OFFSET && return ""
    text = payload[(begin+SLIDERULE_ALERT_TEXT_OFFSET):end]
    stop = findfirst(iszero, text)
    return String(isnothing(stop) ? text : text[1:stop-1])
end

"""
    sliderule_arrow(bytes; filename=nothing, warnings=false) -> Union{Nothing,Vector{UInt8}}

Reassemble the Arrow IPC file carried by a SlideRule response.

Returns `nothing` when the response contains no output file, which is the normal answer when none of
the requested granules had data inside the region. Verifies the assembled length against the size
declared in `arrowrec.meta`, and -- when `filename` is given -- that the server is talking about the
file that was asked for. `warnings=true` echoes the `exceptrec` messages, which are usually just
granules that did not intersect after all.
"""
sliderule_arrow(bytes::Vector{UInt8}; filename=nothing, warnings=false) =
    sliderule_arrow(sliderule_records(bytes); filename, warnings)

function sliderule_arrow(records::Vector{Tuple{String,Vector{UInt8}}}; filename=nothing, warnings=false)
    chunks = UInt8[]
    declared = nothing

    for (rectype, payload) in records
        if rectype == "arrowrec.meta"
            length(payload) >= SLIDERULE_FILENAME_LEN + 8 ||
                throw(NonRetryable("SlideRule arrowrec.meta is $(length(payload)) bytes, too short for " *
                                   "a $(SLIDERULE_FILENAME_LEN) byte filename plus size"))
            name = _sliderule_filename(payload)
            if !isnothing(filename) && name != filename
                throw(NonRetryable("SlideRule returned output file \"$name\" but \"$filename\" was requested"))
            end
            declared = reinterpret(Int64, payload[(begin+SLIDERULE_FILENAME_LEN):(begin+SLIDERULE_FILENAME_LEN+7)])[1]
        elseif rectype == "arrowrec.data"
            length(payload) > SLIDERULE_FILENAME_LEN ||
                throw(NonRetryable("SlideRule arrowrec.data carries no bytes past its filename field"))
            append!(chunks, payload[(begin+SLIDERULE_FILENAME_LEN):end])
        elseif rectype == "exceptrec" && warnings
            message = sliderule_alert(payload)
            isempty(message) || printstyled("    -> sliderule: $message\n"; color=:light_black)
        end
    end

    isempty(chunks) && return nothing

    if !isnothing(declared) && declared > 0 && length(chunks) != declared
        throw(NonRetryable("SlideRule output file is incomplete: assembled $(length(chunks)) of $declared bytes"))
    end
    return chunks
end

"""
    sliderule_failed_resources(bytes) -> Set{String}

Granules whose server-side read failed part-way, recovered from the `exceptrec` stream.

A per-beam read can fail (`H5Coro::Future read failure on BEAM0001/quality_flag`) while the request
as a whole still succeeds: the beam's dataframe is returned with zero rows and the output file is
assembled from whatever did come back. That is indistinguishable, in the data, from a beam that
genuinely had no points in the region -- so without this the granule is written to the archive as
complete, the incremental rule never asks for it again, and the missing beams are lost permanently.

Returning the affected resources lets the build refuse to record them, so the next pass retries. Seen
in practice: 2 of 16 beam-dataframes on the first ATL06 request made against the public cluster, and
8 of 8 on every GEDI02_A v003 granule, which is why v003 cannot be ingested at all.
"""
sliderule_failed_resources(bytes::Vector{UInt8}) = sliderule_failed_resources(sliderule_records(bytes))

function sliderule_failed_resources(records::Vector{Tuple{String,Vector{UInt8}}})
    failed = Set{String}()
    for (rectype, payload) in records
        rectype == "exceptrec" || continue
        resource = _sliderule_failed_resource(sliderule_alert(payload))
        isnothing(resource) || push!(failed, resource)
    end
    return failed
end

# "Failure on resource <resource> beam <beam>: <reason>" -- the alert text SlideRule's readers emit
# when a dataset read throws. Matched on the fixed prefix rather than the reason, which varies.
function _sliderule_failed_resource(message::AbstractString)
    m = match(r"Failure on resource (\S+) beam ", message)
    return isnothing(m) ? nothing : String(m.captures[1])
end

# filename field is a fixed-width, NUL-padded string
function _sliderule_filename(payload::Vector{UInt8})
    field = payload[begin:(begin+SLIDERULE_FILENAME_LEN-1)]
    stop = findfirst(iszero, field)
    return String(isnothing(stop) ? field : field[1:stop-1])
end

"""
    sliderule_atl06_parms(extent; granules=nothing, t0=nothing, t1=nothing, buffer=SLIDERULE_POLY_BUFFER, filename="gga.feather") -> Dict

Request parameters for an `atl06x` query over `extent`.

Passing `granules` (a vector of ATL06 filenames) pins the request to exactly those resources, which is
what makes an incremental update exact: only granules a geotile does not already hold are asked for.
With `granules=nothing` the server queries CMR itself from the polygon and the optional `t0`/`t1`
bounds.

`buffer` grows the requested polygon in degrees; see [`SLIDERULE_POLY_BUFFER`](@ref). Callers are
expected to clip the result to the exact extent, which [`geotile_build_sliderule`](@ref) does.
"""
function sliderule_atl06_parms(extent::Extent; granules=nothing, t0=nothing, t1=nothing,
    buffer=SLIDERULE_POLY_BUFFER, filename="gga.feather")

    x0, x1 = extent.X .+ (-buffer, buffer)
    y0, y1 = extent.Y .+ (-buffer, buffer)
    # keep the request inside valid coordinates at the poles
    y0 = max(y0, -90.0)
    y1 = min(y1, 90.0)
    # closed ring, as SlideRule expects
    poly = [Dict("lat" => y0, "lon" => x0), Dict("lat" => y0, "lon" => x1),
        Dict("lat" => y1, "lon" => x1), Dict("lat" => y1, "lon" => x0),
        Dict("lat" => y0, "lon" => x0)]

    parms = Dict{String,Any}(
        "asset" => "icesat2-atl06",
        "poly" => poly,
        # `height_reference` in the archive is ATL06's own DEM height, which is not part of the
        # default x-series field set.
        "atl06_fields" => SLIDERULE_ATL06_FIELDS,
        # feather is Arrow IPC, which Arrow.jl reads without another dependency
        "output" => Dict("format" => "feather", "path" => filename, "open_on_complete" => false),
    )

    isnothing(granules) || (parms["resources"] = collect(granules))
    isnothing(t0) || (parms["t0"] = _sliderule_time(t0))
    isnothing(t1) || (parms["t1"] = _sliderule_time(t1))
    return parms
end

_sliderule_time(t::DateTime) = Dates.format(t, "yyyy-mm-ddTHH:MM:SS") * "Z"
_sliderule_time(t::AbstractString) = String(t)

# Attempts at a batch before its still-broken granules are given up on and left for the next pass.
#
# Per-beam read failures are common and transient: 5 of 25 granules failed on one live GEDI request,
# and the same 5 succeeded when asked for again. Retrying in place rather than only across passes is
# what keeps a build converging -- at a 20% per-attempt failure rate, three attempts leaves under 1%.
const SLIDERULE_BEAM_ATTEMPTS = 3

"""
    _sliderule_batches(endpoint, filename, make_parms, to_archive, granules, chunk; kwargs...) -> (; frames, failed, partial)

Request `granules` in batches, retrying granules whose server-side read failed, and return the
translated frames along with the granules that never yielded anything (`failed`) and those recorded
without every beam (`partial`).

Shared by the ATL06 and GEDI queries, which differ only in endpoint, parameters and translation.

Rows the translation could not attribute to a requested granule (`id == ""`) are dropped, loudly. An
unattributed row cannot be recognised by the incremental rule, so writing it would put points in the
archive that no later pass can account for -- the same failure mode as the `id == "0"` placeholder rows
this archive had to be repaired for.

With `granules=nothing` the server runs its own CMR query, there is no resource list to narrow a retry
to, so failures are reported without one.

# Transient failures versus beams that are simply not there

A beam read can fail two ways that look identical in the message stream, and they need opposite
handling.

Transient: the read glitches, the beam comes back with zero rows, the request still reports success.
Keeping those rows and moving on would freeze an incomplete granule into the archive, because the id
would be present and no pass would ask again. So while attempts remain, rows for a broken granule are
dropped and only that granule is re-requested.

Permanent: the beam is not readable at all, and no number of retries will change it. Some granules are
simply like this -- release-004 granules where BEAM1000 and BEAM1011 fail every attempt are the case
that surfaced it, though it is per-granule rather than per-release, and ran at 1.8% of granules over a
4,511-granule calibration. Treating those as transient is a trap: the granule would be retried,
rejected and re-requested on every future pass forever, and the beams that *did* read (~8,000 points
each) discarded every time, silently.

The two are separated by whether the failure survives `beam_attempts`. A glitch is unlikely to recur
on the same beam three times while its siblings succeed; an absent beam always will. So on the last
attempt the rows that did arrive are kept, and the granule is reported as `partial` rather than
`failed` -- which is what the HDF5 path does anyway, since `SpaceLiDAR.points` filters the beam list
to the groups the file actually contains. A granule that yields nothing at all across every attempt
stays `failed` and is left unrecorded.
"""
function _sliderule_batches(endpoint, filename, make_parms, to_archive, granules, chunk;
    warnings=false, gps_time=true, poster=HTTP.post, beam_attempts=SLIDERULE_BEAM_ATTEMPTS)

    batches = isnothing(granules) ? [nothing] :
              collect(Iterators.partition(collect(granules), chunk))

    frames = DataFrame[]
    failed = Set{String}()
    partial = Set{String}()
    retried = 0

    for batch in batches
        pending = batch
        attempts = isnothing(batch) ? 1 : beam_attempts

        for attempt in 1:attempts
            last_attempt = attempt == attempts || isnothing(batch)

            parms = make_parms(pending, filename)
            bytes = sliderule_post(endpoint, parms; poster)

            # Parse the record stream once. Both the failure scan and the file reassembly used to take
            # the raw bytes and re-walk them independently, each copying every record out again --
            # measured at 5x the response size in allocations for a 20 MB response, on a run that
            # issued tens of thousands of requests.
            records = sliderule_records(bytes)
            broken = sliderule_failed_resources(records)

            frame = nothing
            payload = sliderule_arrow(records; filename, warnings)
            if !isnothing(payload)
                frame = _sliderule_drop_unattributed(
                    to_archive(DataFrame(Arrow.Table(IOBuffer(payload))), pending; gps_time))
            end

            if !isempty(broken)
                returned = isnothing(frame) || isempty(frame) ? Set{String}() : Set(unique(frame.id))
                if last_attempt
                    # Out of attempts: keep what arrived for a granule whose beams are evidently not
                    # coming, and record it. Only a granule that produced nothing is a failure.
                    union!(partial, intersect(broken, returned))
                    union!(failed, setdiff(broken, returned))
                elseif !isnothing(frame)
                    frame = frame[.!in.(frame.id, Ref(broken)), :]
                end
            end

            isnothing(frame) || isempty(frame) || push!(frames, frame)

            if isempty(broken) || last_attempt
                pending = String[]
                break
            end

            pending = sort(collect(broken))
            retried += length(pending)
        end
    end

    return (; frames, failed, partial, retried)
end

# Rows the translation could not tie back to a requested granule. The archive keys its incremental
# rule on `id`, so a row without one is unreachable: no pass can tell whether it is present, stale or
# duplicated. Dropping is the only safe answer, and it is worth a warning because in production it
# should not happen at all -- it means the returned rows did not come from the resources requested.
function _sliderule_drop_unattributed(frame::DataFrame)
    isempty(frame) && return frame
    unattributed = frame.id .== ""
    n = count(unattributed)
    n == 0 && return frame
    @warn "SlideRule returned $n row(s) that match no requested granule; dropping them" maxlog = 5
    return frame[.!unattributed, :]
end

"""
    sliderule_atl06(extent; granules=nothing, t0=nothing, t1=nothing, buffer=SLIDERULE_POLY_BUFFER, warnings=false, gps_time=true, poster=HTTP.post) -> (DataFrame, Set{String}, Set{String})

Query `atl06x` and return the points inside `extent` in this project's archive schema, together with
the granules whose read failed part-way (see [`sliderule_failed_resources`](@ref)).

Returns an empty, correctly typed DataFrame when the request produced no points. Granules are
requested in batches of `SLIDERULE_RESOURCE_CHUNK`. The result covers `extent` grown by `buffer`, so
callers that need the exact extent must clip -- [`geotile_build_sliderule`](@ref) does.

# Columns
Matches what `SpaceLiDAR.points` plus `SpaceLiDAR.add_id` produce for a local ATL06 granule, so
geotiles from either source are interchangeable: `longitude`, `latitude`, `height`, `height_error`,
`datetime`, `quality`, `track`, `strong_beam`, `detector_id`, `height_reference`, `id`.
"""
function sliderule_atl06(extent::Extent; granules=nothing, t0=nothing, t1=nothing,
    buffer=SLIDERULE_POLY_BUFFER, warnings=false, gps_time=true, poster=HTTP.post,
    beam_attempts=SLIDERULE_BEAM_ATTEMPTS)

    result = _sliderule_batches(SLIDERULE_ATL06_ENDPOINT, "gga_atl06.feather",
        (batch, filename) -> sliderule_atl06_parms(extent; granules=batch, t0, t1, buffer, filename),
        sliderule2archive, granules, SLIDERULE_RESOURCE_CHUNK;
        warnings, gps_time, poster, beam_attempts)

    df = isempty(result.frames) ? sliderule_empty_table() : reduce(vcat, result.frames)
    return (df, result.failed, result.partial, result.retried)
end

"""
    sliderule_empty_table() -> DataFrame

Empty DataFrame with the archive's column names and types.

Needed so a geotile whose granules all returned nothing still produces a table that `vcat`s and
serialises like a populated one.
"""
sliderule_empty_table() = DataFrame(
    longitude=Float64[], latitude=Float64[], height=Float32[], height_error=Float32[],
    datetime=DateTime[], quality=Bool[], track=String[], strong_beam=Bool[],
    detector_id=Int8[], height_reference=Float32[], id=String[])

"""
    sliderule2archive(df, granules; gps_time=true) -> DataFrame

Translate an `atl06x` result into the archive schema.

`granules` is the resource list the points came from, used to recover each row's granule filename.
`gps_time=true` shifts timestamps onto the archive's GPS timescale (see
[`SLIDERULE_GPS_UTC_OFFSET`](@ref)); pass `false` to keep SlideRule's UTC.

# Mapping notes
  - `datetime` is shifted by `SLIDERULE_GPS_UTC_OFFSET` so it agrees with geotiles built from local
    HDF5, which are 18 s ahead of UTC.
  - `quality` is inverted: ATL06 stores `atl06_quality_summary` as 0 for "no problems", while the
    archive's `quality` is true when a point is good.
  - `height_error` combines the two error terms the same way the HDF5 reader does,
    `sqrt(sigma_geo_h^2 + h_li_sigma^2)`.
  - `track` is the beam name (`"gt1l"` ...). SlideRule reports `gt` as 10, 20 ... 60, which indexes
    the beam tuple directly.
  - `strong_beam` derives from `spot`, not from `gt`: which beams are strong swaps with spacecraft
    orientation, but spots 1, 3 and 5 are always the strong ones.
  - `id` is matched on `(rgt, cycle, region)`, which identifies an ATL06 granule uniquely. The
    alternative, SlideRule's `srcid` column, indexes a source table that the feather output does not
    carry.
"""
function sliderule2archive(df::DataFrame, granules; gps_time=true)
    isempty(df) && return sliderule_empty_table()

    height = _nanfill(df.h_li)
    h_li_sigma = _nanfill(df.h_li_sigma)
    sigma_geo_h = _nanfill(df.sigma_geo_h)
    spot = Int8.(coalesce.(df.spot, Int8(0)))
    gt = Int.(coalesce.(df.gt, 0))

    out = DataFrame(
        longitude=Float64.(coalesce.(df.longitude, NaN)),
        latitude=Float64.(coalesce.(df.latitude, NaN)),
        height=height,
        height_error=sqrt.(sigma_geo_h .^ 2 .+ h_li_sigma .^ 2),
        datetime=gps_time ? DateTime.(df.time_ns) .+ SLIDERULE_GPS_UTC_OFFSET : DateTime.(df.time_ns),
        quality=.!Bool.(coalesce.(df.atl06_quality_summary, Int8(1))),
        track=_beam_name.(gt),
        strong_beam=isodd.(spot),
        detector_id=spot,
        height_reference=_dem_column(df),
        id=_granule_ids(df, granules),
    )
    return out
end

"""
    _dem_column(df) -> Vector{Float32}

ATL06 DEM height for each row, for the archive's `height_reference` column.

Warns rather than quietly returning NaNs when the ancillary column is absent: an all-NaN
`height_reference` looks like legitimately missing data, so a rename on the server side would
otherwise pass unnoticed.
"""
function _dem_column(df::DataFrame)
    for name in SLIDERULE_DEM_COLUMNS
        hasproperty(df, name) && return _nanfill(getproperty(df, name))
    end
    @warn "SlideRule response has no ATL06 DEM column $(SLIDERULE_DEM_COLUMNS); " *
          "height_reference will be NaN" maxlog = 1
    return fill(NaN32, nrow(df))
end

# ATL06 marks missing heights with a large sentinel; the rest of the pipeline tests for NaN.
function _nanfill(column)
    values = Float32.(coalesce.(column, NaN32))
    values[values.>=ATL06_FILL_VALUE] .= NaN32
    return values
end

# SlideRule's `gt` is 10, 20 ... 60 for gt1l, gt1r ... gt3r. `icesat2_tracks` is SpaceAltimetry's beam
# tuple, reused rather than restated so the names cannot drift from the HDF5 reader's.
function _beam_name(gt::Integer)
    tracks = SpaceLiDAR.icesat2_tracks
    index = gt ÷ 10
    return (1 <= index <= length(tracks)) ? String(tracks[index]) : ""
end

"""
    _granule_ids(df, granules) -> Vector{String}

Granule filename for each row, matched on `(rgt, cycle, region)`.

An ATL06 filename encodes those three numbers as `ATL06_<datetime>_<rgt><cycle><region>_<ver>_<rev>`,
so they identify the source granule without needing the server's source table. Rows that match no
supplied resource get `""`; with `granules=nothing` (server-side CMR query) the triple is formatted
directly, since no filename is available to match against.
"""
function _granule_ids(df::DataFrame, granules)
    rgt = Int.(coalesce.(df.rgt, 0))
    cycle = Int.(coalesce.(df.cycle, 0))
    region = hasproperty(df, :region) ? Int.(coalesce.(df.region, 0)) : zeros(Int, nrow(df))

    if isnothing(granules) || isempty(granules)
        return [@sprintf("ATL06_%04d%02d%02d", r, c, g) for (r, c, g) in zip(rgt, cycle, region)]
    end

    lookup = Dict{Tuple{Int,Int,Int},String}()
    for granule in granules
        key = _atl06_key(granule)
        isnothing(key) || (lookup[key] = String(granule))
    end
    return [get(lookup, (r, c, g), "") for (r, c, g) in zip(rgt, cycle, region)]
end

"""
    _atl06_key(filename) -> Union{Nothing,Tuple{Int,Int,Int}}

`(rgt, cycle, region)` parsed from an ATL06 granule filename, or `nothing` if it does not parse.
"""
function _atl06_key(filename::AbstractString)
    parts = split(basename(String(filename)), "_")
    length(parts) >= 3 || return nothing
    field = parts[3]
    length(field) == 8 || return nothing
    rgt = tryparse(Int, field[1:4])
    cycle = tryparse(Int, field[5:6])
    region = tryparse(Int, field[7:8])
    (isnothing(rgt) || isnothing(cycle) || isnothing(region)) && return nothing
    return (rgt, cycle, region)
end

# ---------------------------------------------------------------------------------------------------
# GEDI L2A
#
# `gedi02ax` -- the x-series GEDI 2A Dataframe endpoint, which subsets the standard GEDI02_A product.
# Same role for GEDI that `atl06x` plays for ICESat-2, and the same reason for preferring it: the
# alternative `gedi02ap` returns a fixed point record without ancillary field support, and this
# archive needs six ancillary datasets for its columns alone.
#
# Only v002 can be ingested. GEDI02_A v003 removed the per-shot `quality_flag` dataset, which
# SlideRule's L2A reader reads unconditionally, so every beam of every v003 granule throws and the
# request returns zero rows. Nothing here is version-specific -- the filename parsing and the shot
# number layout are identical in v003 -- so this path will pick v003 up unchanged if SlideRule stops
# requiring that dataset.
#
# # Granules with unreadable beams
#
# A minority of granules have one or more beams that cannot be read at all: the request reports a
# failure for that beam on every attempt while the others return normally. First seen on release 004
# granules (acquired after the 2024-04 restart) where BEAM1000 and BEAM1011 fail, but it is a
# per-granule property, not a per-release one -- measured across 4,511 granules of a calibration run,
# 1.8% were affected, and 5.6M post-2024 points from unaffected granules carry all eight beams with the
# same ~0.52 `strong_beam` fraction as the pre-gap archive. See [`_sliderule_batches`](@ref) for why
# such a granule is recorded with the beams that did read rather than rejected.
# ---------------------------------------------------------------------------------------------------

const SLIDERULE_GEDI_ENDPOINT = "gedi02ax"
const SLIDERULE_GEDI_ASSET = "gedil2a"

# Granules per request. Far smaller than the ATL06 chunk: a GEDI granule means eight beams and, with
# the L3 filter's inputs, twenty-one ancillary datasets each, so a hundred of them in one request is
# a great deal of server-side reading to lose to a single failure.
const SLIDERULE_GEDI_RESOURCE_CHUNK = 25

# `flags` bitmask, from `GediParameters::flags_t`. SlideRule packs the three per-shot flags into one
# byte rather than returning them as columns.
const GEDI_DEGRADE_FLAG_MASK = 0x01
const GEDI_L2_QUALITY_FLAG_MASK = 0x02
const GEDI_SURFACE_FLAG_MASK = 0x80

# `digital_elevation_model` fill. Unlike the height fields, which use the shared 3.4e38 sentinel that
# [`_nanfill`](@ref) handles, the ATL06-style DEM field marks missing data with -999999.
const GEDI_DEM_FILL_VALUE = -999999.0f0

# Beams carrying full laser power, as SlideRule's `beam` values. Confirmed against 9.8M rows of the
# existing archive: `strong_beam` is true for exactly BEAM0101, BEAM0110, BEAM1000 and BEAM1011.
# Unlike ICESat-2 this does not depend on spacecraft orientation.
const GEDI_STRONG_BEAMS = (0x05, 0x06, 0x08, 0x0b)

# `shot_number` field layout, right-anchored so it is independent of the orbit's digit count:
#
#     orbit * 10^13 + beam * 10^11 + reserved * 10^9 + granule * 10^8 + shot_index
#
# Verified against a live response: the decoded orbit and beam equal the `orbit` and `beam` columns on
# every row, and the decoded granule number equals the granule field of the filename for every source.
# This is the only per-row route to the source granule -- `srcid` keys a table the feather output does
# not carry, and the `orbit`/`track` columns alone are ambiguous, since the four sub-orbit granules of
# one orbit share a track number (11,785 collisions across the archive's 35,002 granules).
const GEDI_SHOT_ORBIT_DIVISOR = 10^13
const GEDI_SHOT_BEAM_DIVISOR = 10^11
const GEDI_SHOT_GRANULE_DIVISOR = 10^8

# Waveform processing algorithms; `selected_algorithm` names which one's results apply to a shot.
const GEDI_ALGORITHMS = 1:6

# `rx_maxamp / sd_corrected` threshold of the L3 quality filter.
const GEDI_RX_MAXAMP_THRESHOLD = 8

# Ancillary datasets, named relative to the beam group. Six back archive columns; the rest are inputs
# to the L3 quality filter. SlideRule labels an ancillary column with the exact string requested, so
# these are the column names in the response, slashes included.
const SLIDERULE_GEDI_FIELDS = [
    # archive columns
    "elevation_bin0_error", "energy_total", "num_detectedmodes", "digital_elevation_model",
    # L3 filter inputs
    "rx_assess/quality_flag", "geolocation/stale_return_flag",
    "rx_assess/rx_maxamp", "rx_assess/sd_corrected", "selected_algorithm",
    ("rx_processing_a$(a)/zcross" for a in GEDI_ALGORITHMS)...,
    ("rx_processing_a$(a)/toploc" for a in GEDI_ALGORITHMS)...,
]

"""
    sliderule_gedi_parms(extent; granules=nothing, t0=nothing, t1=nothing, buffer=SLIDERULE_POLY_BUFFER, filename="gga.feather") -> Dict

Request parameters for a `gedi02ax` query over `extent`.

Mirrors [`sliderule_atl06_parms`](@ref): `granules` pins the request to exactly those resources, which
is what makes an incremental update exact, and `buffer` grows the polygon in degrees.

`degrade_filter` and `surface_filter` are set because the archive's own filter drops those shots
anyway (see [`_gedi_l3_filter`](@ref)), and applying them server-side means the points are never
transferred. `l2_quality_filter` is deliberately **not** set: the archive keeps `quality` as a column
rather than filtering on it -- 13% of its rows are `quality == false` -- so filtering here would
silently discard data the HDF5 path retains.
"""
function sliderule_gedi_parms(extent::Extent; granules=nothing, t0=nothing, t1=nothing,
    buffer=SLIDERULE_POLY_BUFFER, filename="gga.feather")

    x0, x1 = extent.X .+ (-buffer, buffer)
    y0, y1 = extent.Y .+ (-buffer, buffer)
    y0 = max(y0, -90.0)
    y1 = min(y1, 90.0)
    poly = [Dict("lat" => y0, "lon" => x0), Dict("lat" => y0, "lon" => x1),
        Dict("lat" => y1, "lon" => x1), Dict("lat" => y1, "lon" => x0),
        Dict("lat" => y0, "lon" => x0)]

    parms = Dict{String,Any}(
        "asset" => SLIDERULE_GEDI_ASSET,
        "poly" => poly,
        "anc_fields" => SLIDERULE_GEDI_FIELDS,
        # equivalent to two terms of the archive's own filter, but applied before transfer
        "degrade_filter" => true,
        "surface_filter" => true,
        "output" => Dict("format" => "feather", "path" => filename, "open_on_complete" => false),
    )

    isnothing(granules) || (parms["resources"] = collect(granules))
    isnothing(t0) || (parms["t0"] = _sliderule_time(t0))
    isnothing(t1) || (parms["t1"] = _sliderule_time(t1))
    return parms
end

"""
    sliderule_gedi(extent; granules=nothing, t0=nothing, t1=nothing, buffer=SLIDERULE_POLY_BUFFER, warnings=false, gps_time=true, poster=HTTP.post) -> (DataFrame, Set{String}, Set{String})

Query `gedi02ax` and return the points inside `extent` in this project's archive schema, together with
the granules whose read failed part-way (see [`sliderule_failed_resources`](@ref)).

Granules are requested in batches of `SLIDERULE_GEDI_RESOURCE_CHUNK`. The result covers `extent` grown
by `buffer`, so callers needing the exact extent must clip -- [`geotile_build_sliderule`](@ref) does.

# Columns
Matches what `SpaceLiDAR.points` plus `SpaceLiDAR.add_id` produce for a local GEDI02_A granule, so
geotiles from either source are interchangeable: `longitude`, `latitude`, `height`, `height_error`,
`datetime`, `intensity`, `sensitivity`, `surface`, `quality`, `nmodes`, `track`, `strong_beam`,
`classification`, `sun_angle`, `height_reference`, `id`.
"""
function sliderule_gedi(extent::Extent; granules=nothing, t0=nothing, t1=nothing,
    buffer=SLIDERULE_POLY_BUFFER, warnings=false, gps_time=true, poster=HTTP.post,
    beam_attempts=SLIDERULE_BEAM_ATTEMPTS)

    result = _sliderule_batches(SLIDERULE_GEDI_ENDPOINT, "gga_gedi.feather",
        (batch, filename) -> sliderule_gedi_parms(extent; granules=batch, t0, t1, buffer, filename),
        sliderule2archive_gedi, granules, SLIDERULE_GEDI_RESOURCE_CHUNK;
        warnings, gps_time, poster, beam_attempts)

    df = isempty(result.frames) ? sliderule_gedi_empty_table() : reduce(vcat, result.frames)
    return (df, result.failed, result.partial, result.retried)
end

"""
    sliderule_gedi_empty_table() -> DataFrame

Empty DataFrame with the GEDI archive's column names and types.
"""
sliderule_gedi_empty_table() = DataFrame(
    longitude=Float64[], latitude=Float64[], height=Float32[], height_error=Float32[],
    datetime=DateTime[], intensity=Float32[], sensitivity=Float32[], surface=Bool[],
    quality=Bool[], nmodes=UInt8[], track=String[], strong_beam=Bool[],
    classification=String[], sun_angle=Float32[], height_reference=Float32[], id=String[])

"""
    sliderule2archive_gedi(df, granules; gps_time=true) -> DataFrame

Translate a `gedi02ax` result into the archive schema, applying the same quality filter the HDF5 path
applies (see [`_gedi_l3_filter`](@ref)).

`granules` is the resource list the points came from, used to recover each row's granule filename.
`gps_time=true` shifts timestamps onto the archive's GPS timescale; pass `false` to keep SlideRule's
UTC.

# Mapping notes
  - `datetime` is shifted by [`SLIDERULE_GPS_UTC_OFFSET`](@ref) for the same reason as ATL06:
    SlideRule converts GEDI's `delta_time` to true UTC, while SpaceAltimetry adds the epoch without
    subtracting leap seconds, so every `datetime` in the archive is 18 s ahead of UTC.
  - `quality` is **not** inverted, unlike ATL06: GEDI's `quality_flag` is already 1 for good.
  - `height` is `elev_lowestmode` and the coordinates are `lat/lon_lowestmode`, matching the archive's
    `classification == "ground"`, which is the only classification it contains.
  - `track` is the beam name (`"BEAM0000"` ...), recovered from SlideRule's 4-bit `beam` value, whose
    bits *are* the name's digits.
  - `surface` and `quality` are unpacked from the `flags` bitmask.
"""
function sliderule2archive_gedi(df::DataFrame, granules; gps_time=true)
    isempty(df) && return sliderule_gedi_empty_table()

    keep = _gedi_l3_filter(df)
    any(keep) || return sliderule_gedi_empty_table()
    df = df[keep, :]

    flags = _gedi_flag_column(df, "flags")
    beam = _gedi_flag_column(df, "beam")

    return DataFrame(
        longitude=Float64.(coalesce.(df.longitude, NaN)),
        latitude=Float64.(coalesce.(df.latitude, NaN)),
        height=_nanfill(df.elevation_lm),
        height_error=_nanfill(_gedi_field(df, "elevation_bin0_error")),
        datetime=gps_time ? DateTime.(df.time_ns) .+ SLIDERULE_GPS_UTC_OFFSET : DateTime.(df.time_ns),
        intensity=_nanfill(_gedi_field(df, "energy_total")),
        sensitivity=_nanfill(df.sensitivity),
        surface=(flags .& GEDI_SURFACE_FLAG_MASK) .!= 0,
        quality=(flags .& GEDI_L2_QUALITY_FLAG_MASK) .!= 0,
        nmodes=UInt8.(coalesce.(_gedi_field(df, "num_detectedmodes"), 0x00)),
        track=_gedi_beam_name.(beam),
        strong_beam=in.(beam, Ref(GEDI_STRONG_BEAMS)),
        classification=fill("ground", nrow(df)),
        sun_angle=_nanfill(df.solar_elevation),
        height_reference=_gedi_dem_column(df),
        id=_gedi_granule_ids(df, granules),
    )
end

"""
    _gedi_l3_filter(df) -> Vector{Bool}

Rows to keep, reproducing the quality filter `SpaceLiDAR.points` applies to GEDI02_A by default.

The existing archive was built through `getpoints`, which calls `SpaceLiDAR.points` without
`filtered=false`, so every point in it has passed this filter. Not reapplying it here would make the
appended data systematically noisier than what it is appended to, which no downstream binning would
reveal as a bug -- so a missing input dataset is an error, not a reason to pass rows through.

The terms, after `degrade_flag` and `surface_flag` (both also enforced server-side):

  - `rx_assess/quality_flag != 0`
  - `geolocation/stale_return_flag == 0`
  - `rx_assess/rx_maxamp / rx_assess/sd_corrected >= GEDI_RX_MAXAMP_THRESHOLD`
  - `zcross > 0` and `toploc > 0`, read from the `rx_processing_a<n>` group named by
    `selected_algorithm`

One term is deliberately omitted. `SpaceLiDAR.points` also means to filter on `rx_algrunflag`, but
the line reads `m .& algrun .!= 0` -- no in-place assignment, so the result is discarded and the
filter has never had any effect. The archive therefore contains rows that term would have dropped,
and reproducing the *intent* here would make the two sources disagree. This mirrors the archive, and
`rx_algrunflag` is left out of [`SLIDERULE_GEDI_FIELDS`](@ref) accordingly.
"""
function _gedi_l3_filter(df::DataFrame)
    n = nrow(df)
    keep = trues(n)

    keep .&= _gedi_flag_column(df, "rx_assess/quality_flag") .!= 0
    keep .&= _gedi_flag_column(df, "geolocation/stale_return_flag") .== 0

    maxamp = _nanfill(_gedi_field(df, "rx_assess/rx_maxamp"))
    sd_corrected = _nanfill(_gedi_field(df, "rx_assess/sd_corrected"))
    keep .&= (maxamp ./ sd_corrected) .>= GEDI_RX_MAXAMP_THRESHOLD

    flags = _gedi_flag_column(df, "flags")
    keep .&= (flags .& GEDI_SURFACE_FLAG_MASK) .!= 0
    keep .&= (flags .& GEDI_DEGRADE_FLAG_MASK) .== 0

    # zcross/toploc live in a per-algorithm group, so gather the selected algorithm's values. Rows
    # naming an algorithm outside 1:6 keep NaN and are dropped, which is the safe reading of a value
    # the product should never produce.
    algorithm = _gedi_flag_column(df, "selected_algorithm")
    zcross = fill(NaN32, n)
    toploc = fill(NaN32, n)
    for a in GEDI_ALGORITHMS
        rows = algorithm .== a
        any(rows) || continue
        zcross[rows] = _nanfill(_gedi_field(df, "rx_processing_a$(a)/zcross"))[rows]
        toploc[rows] = _nanfill(_gedi_field(df, "rx_processing_a$(a)/toploc"))[rows]
    end
    keep .&= zcross .> 0
    keep .&= toploc .> 0

    return keep
end

"""
    _gedi_field(df, name) -> AbstractVector

An ancillary or native column of a `gedi02ax` response, by its exact requested name.

Absence is a [`NonRetryable`](@ref) error rather than a warning: every caller either builds an archive
column from the result or filters on it, so a renamed field would otherwise turn into an all-NaN
column or a filter term that silently passes everything.
"""
function _gedi_field(df::DataFrame, name::AbstractString)
    column = Symbol(name)
    hasproperty(df, column) || throw(NonRetryable(
        "SlideRule gedi02ax response has no column \"$name\"; the field set or its naming has " *
        "changed and the archive schema or quality filter cannot be reproduced"))
    return getproperty(df, column)
end

# Small unsigned columns -- flags, beam id, algorithm id, mode count -- widened to Int so masking and
# comparison cannot overflow or trip on Union{Missing,UInt8}.
_gedi_flag_column(df::DataFrame, name::AbstractString) = Int.(coalesce.(_gedi_field(df, name), 0))

"""
    _gedi_dem_column(df) -> Vector{Float32}

`height_reference` from GEDI's `digital_elevation_model`, with its -999999 fill mapped to NaN.

11.5% of the existing archive's rows are NaN here, so an all-NaN column would not look wrong -- hence
[`_gedi_field`](@ref) throwing on absence rather than warning.
"""
function _gedi_dem_column(df::DataFrame)
    values = _nanfill(_gedi_field(df, "digital_elevation_model"))
    values[values.<=GEDI_DEM_FILL_VALUE] .= NaN32
    return values
end

# SlideRule's `beam` value is the binary number the beam is named after: 5 -> "BEAM0101", 11 ->
# "BEAM1011". Checked against `gedi_tracks`, SpaceAltimetry's beam tuple, so the names cannot drift
# from the HDF5 reader's.
function _gedi_beam_name(beam::Integer)
    (0 <= beam <= 15) || return ""
    name = "BEAM" * string(beam; base=2, pad=4)
    return name in SpaceLiDAR.gedi_tracks ? name : ""
end

"""
    _gedi_granule_ids(df, granules) -> Vector{String}

Granule filename for each row, matched on `(orbit, granule number, track)`.

Orbit and granule number are decoded from `shot_number` (see [`GEDI_SHOT_ORBIT_DIVISOR`](@ref)) and
the track comes from the `track` column. A GEDI02_A filename encodes the three as
`GEDI02_A_<datetime>_O<orbit>_<granule>_T<track>_...`, and the triple is unique across the archive's
35,002 granules, where `(orbit, track)` alone is not.

Rows matching no supplied resource get `""`; with `granules=nothing` (server-side CMR query) the
triple is formatted directly, since no filename is available to match against.
"""
function _gedi_granule_ids(df::DataFrame, granules)
    shot = UInt64.(coalesce.(df.shot_number, UInt64(0)))
    orbit = Int.(shot .÷ GEDI_SHOT_ORBIT_DIVISOR)
    granule_number = Int.((shot .÷ GEDI_SHOT_GRANULE_DIVISOR) .% 10)
    track = _gedi_flag_column(df, "track")

    keys = zip(orbit, granule_number, track)
    if isnothing(granules) || isempty(granules)
        return [@sprintf("GEDI02_A_O%05d_%02d_T%05d", o, g, t) for (o, g, t) in keys]
    end

    lookup = Dict{Tuple{Int,Int,Int},String}()
    for granule in granules
        key = _gedi_key(granule)
        isnothing(key) || (lookup[key] = String(granule))
    end
    return [get(lookup, key, "") for key in keys]
end

"""
    _gedi_key(filename) -> Union{Nothing,Tuple{Int,Int,Int}}

`(orbit, granule number, track)` parsed from a GEDI02_A granule filename, or `nothing` if it does not
parse. Version-independent: v002 and v003 names have identical field layouts.
"""
function _gedi_key(filename::AbstractString)
    parts = split(basename(String(filename)), "_")
    length(parts) >= 6 || return nothing
    (startswith(parts[4], "O") && startswith(parts[6], "T")) || return nothing
    orbit = tryparse(Int, parts[4][2:end])
    granule_number = tryparse(Int, parts[5])
    track = tryparse(Int, parts[6][2:end])
    (isnothing(orbit) || isnothing(granule_number) || isnothing(track)) && return nothing
    return (orbit, granule_number, track)
end

# Per-mission query function. Each returns `(DataFrame in archive schema, failed resources)` for one
# extent and granule list, which is the only thing the build differs on between missions.
const SLIDERULE_QUERIES = Dict{Symbol,Function}(
    :icesat2 => sliderule_atl06,
    :gedi => sliderule_gedi,
)

"""
    geotile_build_sliderule(geotile_granules, geotile_dir; mission=:icesat2, warnings=false, fmt=:arrow, ntasks=8, poster=HTTP.post)

Build or update geotiles from SlideRule instead of local HDF5 granules.

Same contract as [`geotile_build`](@ref): each geotile is independent, granules already present in a
geotile are not requested again, requested granules that return no data leave a placeholder row so
they are not requested on the next pass, and every file is written through [`atomic_write`](@ref).
Only the source of the points differs, so an archive can be populated by either path.

`mission` selects the query through [`SLIDERULE_QUERIES`](@ref) -- `:icesat2` for ATL06, `:gedi` for
GEDI02_A. Everything after the query is mission-independent, including the incremental rule and the
handling of granules whose server-side read failed part-way, which are left unrecorded so the next
pass retries them (see [`sliderule_failed_resources`](@ref)).

# Arguments
- `geotile_granules`: DataFrame of geotiles with `id`, `extent` and `granules` columns, as
  [`granules_load`](@ref) returns
- `geotile_dir`: directory holding the per-geotile files

# Keywords
- `mission=:icesat2`: which SlideRule query to use
- `warnings=false`: echo SlideRule's per-granule notices
- `fmt=:arrow`: output format
- `ntasks=8`: geotiles in flight at once. Unlike the HDF5 path this is task-parallel rather than
  process-parallel -- the reading happens on SlideRule's nodes, so there is no libhdf5 lock to route
  around -- and the ceiling is the cluster's tolerance, not local CPU.
- `poster=HTTP.post`: injection point for tests
"""
function geotile_build_sliderule(geotile_granules, geotile_dir; mission=:icesat2, warnings=false,
    fmt=:arrow, ntasks=8, poster=HTTP.post)

    query = get(SLIDERULE_QUERIES, mission) do
        throw(NonRetryable("source=:sliderule is not implemented for mission :$mission " *
                           "(have $(join(sort(collect(keys(SLIDERULE_QUERIES))), ", ")))"))
    end

    printstyled("building $(mission) geotiles from sliderule\n"; color=:blue, bold=true)

    geotile_granules = geotile_granules[.!isempty.(geotile_granules.granules), :]
    if isempty(geotile_granules)
        throw(NonRetryable("no granules for any requested geotile -- run stages=(:search,) to build " *
                           "the remote granule list first"))
    end

    if !warnings
        Logging.disable_logging(Logging.Warn)
    end

    progress = Progress(nrow(geotile_granules); dt=1, desc="Building geotiles from sliderule...")

    asyncmap(eachrow(geotile_granules); ntasks) do row
        try
            _build_geotile_sliderule(row, geotile_dir; query, fmt, warnings, poster)
        catch e
            e isa InterruptException && rethrow(e)
            printstyled("\n    -> $(row.id): failed [$(sprint(showerror, e))]\n"; color=:light_red)
        finally
            next!(progress)
        end
    end
    finish!(progress)
    return nothing
end

function _build_geotile_sliderule(row, geotile_dir; query=sliderule_atl06, fmt=:arrow,
    warnings=false, poster=HTTP.post)
    outfile = joinpath(geotile_dir, row.id * ".$fmt")
    wanted = [g.id for g in row.granules]

    df0 = nothing
    if isfile(outfile)
        df0 = fmt == :arrow ? DataFrame(Arrow.Table(outfile)) : FileIO.load(outfile, "df")
        if !isempty(df0)
            have = Set(unique(df0.id))
            wanted = filter(!in(have), wanted)
        end
    end

    if isempty(wanted)
        printstyled("\n    -> $(row.id): no new granules to add to exisitng GeoTile\n"; color=:light_green)
        return nothing
    end

    t1 = time()
    requested = length(wanted)
    df, failed, partial, retried = query(row.extent; granules=wanted, warnings, poster)

    # Subsetting happens server-side against the same polygon, but keep the client-side clip so both
    # sources apply identical bounds to the archive.
    if !isempty(df)
        keep = within.(Ref(row.extent), df.longitude, df.latitude) .| isnan.(df.longitude)
        deleteat!(df, .!keep)
    end

    # A granule that produced nothing across every attempt is left out of the file entirely, so the
    # next pass asks for it again. Recording it -- even as a placeholder -- would freeze it into the
    # archive, since the incremental rule keys on the id being present.
    if !isempty(failed)
        isempty(df) || deleteat!(df, in.(df.id, Ref(failed)))
        wanted = filter(!in(failed), wanted)
        printstyled("\n    -> $(row.id): $(length(failed)) granule(s) returned nothing after " *
                    "$(SLIDERULE_BEAM_ATTEMPTS) attempts, left for the next pass\n"; color=:light_yellow)
    end

    # Granules recorded without every beam -- see `_sliderule_batches`. Ran at 1.8% of granules on a
    # calibration run. Reported rather than silent so a rate that climbs, which would mean the archive
    # is being quietly thinned, is visible.
    if !isempty(partial)
        printstyled("\n    -> $(row.id): $(length(partial)) granule(s) recorded without all beams\n";
            color=:light_black)
    end

    # Record granules that came back empty, so the next pass does not ask for them again. This is the
    # same bookkeeping `geotile_build` does with `emptyrow`.
    returned = isempty(df) ? Set{String}() : Set(unique(df.id))
    missing_ids = filter(!in(returned), wanted)
    if !isempty(missing_ids)
        placeholder = emptyrow(df)
        # By name, not `[end]` -- a placeholder carrying `emptyrow`'s "0" rather than its granule id
        # never matches the already-present check, so that granule would be re-requested every pass.
        id_column = columnindex(df, :id)
        id_column == 0 && error("point table has no :id column; cannot record placeholders")
        for id in missing_ids
            placeholder[id_column] = id
            push!(df, placeholder)
        end
    end

    if !isnothing(df0) && !isempty(df0)
        df = vcat(df0::DataFrame, df::DataFrame)
    end
    read_time = round((time() - t1) / 60, digits=1)

    # Nothing new to write: every requested granule failed, so leave the file as it was rather than
    # rewriting it identically.
    if isempty(wanted)
        printstyled("\n    -> $(row.id): nothing recorded, all requested granules failed\n"; color=:light_red)
        return nothing
    end

    t1 = time()
    atomic_write(outfile; suffix=".$fmt") do tmp
        if fmt == :arrow
            Arrow.write(tmp, df::DataFrame)
        else
            save(tmp, Dict("df" => df::DataFrame))
        end
    end
    write_time = round((time() - t1) / 60, digits=1)
    # Seconds per granule is the number that predicts a full run, so report it rather than leaving it
    # to be divided out of the log afterwards.
    per_granule = round(read_time * 60 / max(length(wanted), 1), digits=1)
    # Retry rate is the number that predicts throughput. A run at 12 concurrent geotiles spent 28% of
    # its granule-requests on retries, each one a fresh request with its own backoff, which is where
    # 12.7 s/granule came from -- so surface it per geotile rather than leaving it to be counted out
    # of the log afterwards.
    retry_pct = round(100 * retried / max(requested, 1), digits=1)
    printstyled("\n    -> $(row.id): generation complete [$(length(wanted)) granules, " *
                "read: $read_time min ($(per_granule) s/granule), retries: $(retry_pct)%, " *
                "write: $write_time min]\n"; color=:light_black)
    return nothing
end
