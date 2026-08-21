# update_atl06_sliderule.jl
#
# Refresh the ICESat-2 ATL06 geotile archive with recently acquired data, using SlideRule.
#
# Usage:
#   julia --project -t 16 src/update_atl06_sliderule.jl            # search, then update
#   julia --project     src/update_atl06_sliderule.jl --no-search   # reuse the existing granule list
#
# What it does
#
#   1. Re-runs the CMR search, so the per-geotile granule lists include anything published since the
#      last pass. This is the only stage that needs threads.
#   2. Asks SlideRule for the granules each geotile does not already hold, and appends the points.
#
# Nothing is downloaded: SlideRule reads the granules in its own region and returns only the points
# inside each geotile, so an update costs roughly what the new data is worth rather than a granule at a
# time. Interrupting is safe -- geotiles already holding a granule are skipped on the next run, and
# every file is written atomically.
#
# Credentials are optional but worth setting: with `SLIDERULE_GITHUB_TOKEN` in the environment the
# request is authenticated, which lifts the public cluster's anonymous rate limits. Set
# `SLIDERULE_ORGANIZATION` as well to use a private cluster.

using Dates
using Printf
import GlobalGlacierAnalysis as GGA

do_search = !("--no-search" in ARGS)

project_id = :v01
geotile_width = 2
domain = :glacier

products = GGA.project_products(; project_id)
paths = GGA.project_paths(; project_id)

printstyled("ICESat-2 ATL06 v$(products.icesat2.version) update via SlideRule\n"; color=:blue, bold=true)
println("  cluster:  ", GGA.sliderule_host())
println("  authenticated: ", !isnothing(GGA.sliderule_token()))
println("  geotiles: ", paths.icesat2.geotile)

geotiles = GGA.geotiles_w_mask(geotile_width; remake=false)
geotiles = geotiles[geotiles[!, "$(domain)_frac"].>0, :]
println("  $(domain) geotiles: ", nrow(geotiles))

before = count(endswith(".arrow"), readdir(paths.icesat2.geotile))
t_start = time()

if do_search
    printstyled("\nsearching CMR for new granules\n"; color=:blue, bold=true)
    GGA.geotile_search_granules(geotiles, :icesat2, products.icesat2.name, products.icesat2.version,
        paths.icesat2.granules_remote)
else
    println("\nskipping search; using the granule list already on disk")
end

GGA.geotile_build_archive(; project_id, geotile_width, domain, missions=(:icesat2,),
    stages=(:build,), source=:sliderule)

after = count(endswith(".arrow"), readdir(paths.icesat2.geotile))
@printf("\nupdate finished in %.2f h; geotiles %d -> %d\n", (time() - t_start) / 3600, before, after)
