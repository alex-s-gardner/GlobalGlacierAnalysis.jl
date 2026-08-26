# update_gedi_sliderule.jl
#
# Refresh the GEDI02_A geotile archive with recently acquired data, using SlideRule.
#
# Usage:
#   julia --project -t 16 src/update_gedi_sliderule.jl                # search, then update
#   julia --project      src/update_gedi_sliderule.jl --no-search     # reuse the granule list on disk
#   julia --project      src/update_gedi_sliderule.jl --geotiles=5    # pilot on 5 geotiles
#   julia --project      src/update_gedi_sliderule.jl --nodes=0       # shared public cluster instead
#
# What it does
#
#   1. Re-runs the CMR search, so the per-geotile granule lists include everything published since the
#      last pass. This is the only stage that needs threads.
#   2. Asks SlideRule for the granules each geotile does not already hold, and appends the points.
#
# Nothing is downloaded: SlideRule reads the granules in its own AWS region and returns only the points
# inside each geotile. That matters more for GEDI than for any other mission here -- the local raw
# GEDI02_A archive is already 79 TB, and this update adds to it not at all. Interrupting is safe;
# geotiles already holding a granule are skipped next run and every file is written atomically.
#
# Dedicated capacity
#
# By default this provisions `--nodes` dedicated nodes for the duration of the run and keeps their
# time-to-live refreshed. Tens of thousands of granule reads do not belong on the shared public
# cluster: the anonymous rate limits make it slower for us and it crowds out everyone else. Pass
# `--nodes=0` to use the public cluster anyway. Requests go to the same hostname either way -- the
# gateway routes to our nodes based on the bearer token -- so nothing else about the run changes.
#
# Before the first run
#
#   * Run `src/repair_gedi_placeholders.jl` FIRST, and before any search. The archive's "no data here"
#     placeholder rows are missing their granule ids, so without the repair this update re-requests
#     43,846 granules it already knows are empty -- roughly doubling the work.
#   * Put a GitHub personal access token in `~/.sliderule_pat` (`chmod 600`). Authentication lifts the
#     public cluster's anonymous rate limits, which is the binding constraint on a job this size. Set
#     `SLIDERULE_ORGANIZATION` as well to use a private cluster.
#
# Expect several passes. A per-beam read can fail server-side; those granules are retried in place
# (`SLIDERULE_BEAM_ATTEMPTS`) and, if still broken, left for the next run rather than recorded as
# complete. Re-run until it reports nothing outstanding.
#
# A granule every one of whose reporting beams fails is the exception: nothing in it can be read, so it
# is recorded as an unreadable placeholder and stops being requested, and a run is not waiting on it. To
# put such granules back in front of the incremental rule -- if SlideRule's L2A reader gains the ability
# to read them -- delete the archive rows whose `track` is `PLACEHOLDER_TRACK_UNREADABLE` and re-run.
#
# Version: this updates GEDI02_A v002, which the archive is built on. v003 exists but SlideRule cannot
# read it -- its L2A reader requires the per-shot `quality_flag` dataset that v003 removed, so every
# beam fails. v003 is a reprocessing rather than new coverage (both end 2025-07), so a v002 update
# reaches the same observations.

using Arrow
using Dates
using Printf
using DataFrames
import GlobalGlacierAnalysis as GGA

# Redirected to a file, Julia block-buffers stdout and ProgressMeter goes quiet on a non-TTY, so a
# run logged with `> log 2>&1` shows nothing at all for its first half hour and looks hung. Since this
# job is measured in hours and will normally be logged, flush at every checkpoint.
say(args...) = (println(args...); flush(stdout))

"""
    parse_args(args) -> (do_search, geotile_limit, node_capacity, ntasks)

Read the command line.

In a function, not a top-level loop: assigning to a global from inside a top-level `for` creates a new
local instead, so every flag parses into nothing and is silently discarded. That is not a hypothetical
-- it sent a run that was asked for 5 geotiles on the public cluster off across all 256 at 30-way
concurrency instead. An unrecognised flag is an error for the same reason: a mistyped `--geotiles`
must not quietly become a full run.
"""
function parse_args(args)
    do_search = true
    geotile_limit = nothing
    node_capacity = 10          # nodes the account is configured for
    ntasks = 30                 # concurrent geotiles; ~3 requests per node
    ntasks_given = false

    for arg in args
        if arg == "--no-search"
            do_search = false
        elseif startswith(arg, "--geotiles=")
            geotile_limit = parse(Int, split(arg, "=")[2])
        elseif startswith(arg, "--nodes=")
            node_capacity = parse(Int, split(arg, "=")[2])
        elseif startswith(arg, "--ntasks=")
            ntasks = parse(Int, split(arg, "=")[2])
            ntasks_given = true
        else
            error("unrecognized argument \"$arg\" -- expected --no-search, --geotiles=N, " *
                  "--nodes=N or --ntasks=N")
        end
    end

    # Without dedicated nodes the default 30 is too many to point at the shared cluster, so back it
    # off -- but only when it was left at the default. An explicit --ntasks is a decision, and quietly
    # overriding it would be the same class of bug as the flags that used to be silently discarded.
    if node_capacity == 0 && !ntasks_given
        ntasks = min(ntasks, 8)
    end
    return (do_search, geotile_limit, node_capacity, ntasks)
end

do_search, geotile_limit, node_capacity, ntasks = parse_args(ARGS)

project_id = :v01
geotile_width = 2
domain = :glacier

products = GGA.project_products(; project_id)
paths = GGA.project_paths(; project_id)

printstyled("GEDI $(products.gedi.name) v$(products.gedi.version) update via SlideRule\n";
    color=:blue, bold=true)
say("  cluster:       ", GGA.sliderule_host())
say("  authenticated: ", !isnothing(GGA.sliderule_token()))
say("  identity:      ", something(GGA.sliderule_service(), "anonymous"))
say("  nodes:         ", node_capacity == 0 ? "shared public cluster" : "$(node_capacity) dedicated")
say("  concurrency:   ", ntasks, " geotiles")
say("  geotiles:      ", paths.gedi.geotile)

if !isfile(paths.gedi.granules_remote * ".pre-update")
    printstyled("\n  !! src/repair_gedi_placeholders.jl has not been run.\n"; color=:light_yellow)
    printstyled("     Without it this update re-requests ~43,846 granules already known to be " *
                "empty here.\n"; color=:light_yellow)
end

geotiles = GGA.geotiles_w_mask(geotile_width; remake=false)
geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :]
say("  $(domain) geotiles: ", nrow(geotiles))

before = count(endswith(".arrow"), readdir(paths.gedi.geotile))
t_start = time()

if do_search
    printstyled("\nsearching CMR for new granules\n"; color=:blue, bold=true)
    # No `after=`: the list is rebuilt wholesale, so bounding the search would drop every granule
    # acquired before that date and make the archive look complete when it is not.
    GGA.geotile_search_granules(geotiles, :gedi, products.gedi.name, products.gedi.version,
        paths.gedi.granules_remote)
else
    println("\nskipping search; using the granule list already on disk")
end

geotile_granules = GGA.granules_load(paths.gedi.granules_remote, :gedi; geotiles)

# Report the work outstanding, so a run that is about to do nothing says so before spending an hour
# discovering it, and so the pilot can be pointed at geotiles that actually have something to fetch.
function outstanding(row)
    file = joinpath(paths.gedi.geotile, row.id * ".arrow")
    wanted = Set(g.id for g in row.granules)
    isfile(file) || return length(wanted)
    return length(setdiff(wanted, Set(String.(unique(Arrow.Table(file).id)))))
end

printstyled("\nassessing outstanding granules\n"; color=:blue, bold=true); flush(stdout)
geotile_granules[!, :outstanding] = outstanding.(eachrow(geotile_granules))
@printf("  %d granule-geotile fetches across %d geotiles\n",
    sum(geotile_granules.outstanding), count(>(0), geotile_granules.outstanding))

geotile_granules = geotile_granules[geotile_granules.outstanding.>0, :]
if isempty(geotile_granules)
    printstyled("\nnothing to do -- the archive already holds every granule in the list\n"; color=:green)
    exit()
end

if !isnothing(geotile_limit)
    # Pilot: sample evenly across the outstanding-granule distribution rather than taking the head of
    # it. Taking the smallest picks the tiles above GEDI's +/-51.6 degree coverage limit, which return
    # nothing at all -- a pilot that "passed" against those measured no throughput and would have hidden
    # a bug that silently discarded every post-2024 granule. A spread exercises real work.
    sort!(geotile_granules, :outstanding)
    n = nrow(geotile_granules)
    take = min(geotile_limit, n)
    picks = unique(round.(Int, range(1, n; length=take)))
    geotile_granules = geotile_granules[picks, :]
    @warn "!!!!!!!!!!!!!! PILOT RUN: $(nrow(geotile_granules)) geotiles sampled across the " *
          "distribution, $(sum(geotile_granules.outstanding)) granules " *
          "($(minimum(geotile_granules.outstanding))-$(maximum(geotile_granules.outstanding)) each) " *
          "!!!!!!!!!!!!!!"
end

build() = GGA.geotile_build_sliderule(geotile_granules, paths.gedi.geotile;
    mission=:gedi, warnings=false, ntasks)

if node_capacity > 0
    printstyled("\nprovisioning dedicated capacity\n"; color=:blue, bold=true)
    GGA.sliderule_keepalive(build; node_capacity)
else
    build()
end

after = count(endswith(".arrow"), readdir(paths.gedi.geotile))
@printf("\nupdate finished in %.2f h; geotiles %d -> %d\n", (time() - t_start) / 3600, before, after)
println("Re-run until no geotile reports granules left for the next pass.")
