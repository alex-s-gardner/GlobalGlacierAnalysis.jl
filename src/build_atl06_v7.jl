# build_atl06_v7.jl
#
# Build the ICESat-2 ATL06 geotile archive from granules already sitting in
# `icesat2/ATL06/<version>/raw`, spread over several processes.
#
# Usage:
#   julia --project src/build_atl06_v7.jl          # 8 partitions
#   julia --project src/build_atl06_v7.jl 16       # 16 partitions
#
# Run the search and download stages first -- this script only builds:
#   julia --project -t 16 -e 'import GlobalGlacierAnalysis as GGA; GGA.geotile_build_archive(; missions=(:icesat2,), stages=(:search,))'
#   julia --project     -e 'import GlobalGlacierAnalysis as GGA; GGA.geotile_build_archive(; missions=(:icesat2,), stages=(:download,))'
#
# Why processes rather than threads: HDF5.jl routes every libhdf5 call through one global
# `ReentrantLock`, so threaded reads serialize no matter how many threads are available. Each
# partition is a separate process with its own libhdf5 state, and `partition=(i, N)` hands each one a
# disjoint stride of geotiles. Every output goes through `atomic_write`, so partitions cannot corrupt
# each other and an interrupted run leaves no stubs -- rerunning picks up where it stopped and only
# appends granules a geotile does not already hold.

using Dates
using Printf
import GlobalGlacierAnalysis as GGA

nparts = isempty(ARGS) ? 8 : parse(Int, ARGS[1])
nparts >= 1 || error("nparts must be >= 1, got $nparts")

project_id = :v01
products = GGA.project_products(; project_id)
paths = GGA.project_paths(; project_id)
version = products.icesat2.version
geotile_dir = paths.icesat2.geotile
raw_dir = paths.icesat2.raw_data

printstyled("ICESat-2 ATL06 v$(version) archive build\n"; color=:blue, bold=true)
println("  raw:      ", raw_dir)
println("  geotiles: ", geotile_dir)
println("  partitions: ", nparts)

# ---- preflight ------------------------------------------------------------------------------
# `geotile_build` reads the *local* granule list, which the download stage writes only after every
# transfer finishes. Its absence means the download has not completed, and building now would send
# `getpoints` looking for files that are not there yet.
if !isfile(paths.icesat2.granules_local)
    printstyled("\nno local granule list at $(paths.icesat2.granules_local)\n"; color=:red, bold=true)
    println("The download stage writes it when it finishes. Check progress with:")
    @printf("  ls %s | wc -l\n", raw_dir)
    println("and run the download stage if it is not already running:")
    println("  julia --project -e 'import GlobalGlacierAnalysis as GGA; GGA.geotile_build_archive(; missions=(:icesat2,), stages=(:download,))'")
    exit(1)
end

n_raw = length(readdir(raw_dir))
n_built_before = count(endswith(".arrow"), readdir(geotile_dir))
println("  granules on disk: ", n_raw)
println("  geotiles already built: ", n_built_before)

# ---- launch ---------------------------------------------------------------------------------
logdir = mkpath(joinpath(dirname(geotile_dir), "build_logs"))
julia = Base.julia_cmd()[1]
project = dirname(Base.active_project())

println("  logs: ", logdir, "\n")
t_start = time()

procs = map(1:nparts) do i
    code = """
    import GlobalGlacierAnalysis as GGA
    GGA.geotile_build_archive(; missions=(:icesat2,), stages=(:build,), partition=($i, $nparts))
    """
    log = joinpath(logdir, @sprintf("partition_%02d_of_%02d.log", i, nparts))
    cmd = `$julia --project=$project -t 1 -e $code`
    printstyled("    -> launched partition $i of $nparts\n"; color=:light_black)
    run(pipeline(cmd; stdout=log, stderr=log); wait=false)
end

# ---- monitor --------------------------------------------------------------------------------
# Progress is the count of finished geotiles, which is what the partitions are actually producing.
monitor = @async while any(process_running, procs)
    sleep(300)
    built = count(endswith(".arrow"), readdir(geotile_dir))
    elapsed = (time() - t_start) / 3600
    @printf("    [%5.2f h] %d geotiles built (+%d), %d/%d partitions running\n",
        elapsed, built, built - n_built_before, count(process_running, procs), nparts)
    flush(stdout)
end

for p in procs
    wait(p)
end

# ---- report ---------------------------------------------------------------------------------
elapsed = (time() - t_start) / 3600
built = count(endswith(".arrow"), readdir(geotile_dir))
failed = [i for (i, p) in enumerate(procs) if !success(p)]

println()
printstyled(@sprintf("build finished in %.2f h\n", elapsed); color=:blue, bold=true)
println("  geotiles built: ", built, " (was ", n_built_before, ")")

if isempty(failed)
    printstyled("  all $nparts partitions exited cleanly\n"; color=:green)
else
    printstyled("  partitions that failed: " * join(failed, ", ") * "\n"; color=:red, bold=true)
    println("  inspect: ", logdir)
    println("  rerunning is safe -- finished geotiles are skipped, so only the gaps are rebuilt.")
end
