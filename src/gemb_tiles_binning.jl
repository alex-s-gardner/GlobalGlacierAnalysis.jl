# =============================================================================
# Build the GEMB forcing ensemble from GEMB_GlacierSims tile NetCDF.
# =============================================================================
#
# Replaces `gemb_classes_binning.jl`. That script's job was mostly to repair sparse point output:
# buffer each geotile's search extent until the elevation profile filled, interpolate and extrapolate
# across height bins, and fake elevation classes by scaling melt. A tile file is complete over its own
# geotile's height range and carries a real (ΔT × precipitation scaling) forcing matrix, so all of that
# is gone and this script only reads, area-weights and rebins. Height interpolation survives in one
# place: a geotile the sweep does not cover borrows the nearest tile, whose bands need not line up with
# the target's ice.
#
# The forcing matrix's second axis is a temperature offset in K (`:ΔT`), where the `.mat` path used a
# melt multiplier (`:mscale`). Downstream code detects which it is from the sign of the axis values.
#
# Output: one DimStack per surface mask over (geotile, date, pscale, ΔT) in km³ of ice equivalent,
# area-weighted by that mask's hypsometry, saved as "gemb_dv" with its "surface_mask" to
# `gemb_geotile_dv_file(gembinfo, surface_mask)` and read back by
# `gemb_ensemble_dv(; gemb_run_id, surface_mask)`. Altimetry members are calibrated against the ensemble
# of their own mask.
# GEMB output is in metres of ice equivalent at the density each tile records as `model_density_ice`.

begin
    import GlobalGlacierAnalysis as GGA
    using DataFrames
    using Dates
    using FileIO
    using DimensionalData
    using Statistics

    # Restrict to a single RGI region, or `nothing` for every geotile holding ice.
    rgi_subset = nothing

    single_geotile_test = nothing # e.g. "lat[+62+64]lon[-152-150]"

    project_id = :v01
    geotile_width = 2
    gemb_surface_masks = [:glacier, :glacier_rgi7]
    gemb_run_id = 8

    gembinfo = GGA.gemb_info(; gemb_run_id)

    # Date bins: the project's 30-day spacing, extended back to 1970 while keeping the existing edges
    # so no bin already in use moves.
    date_range, _ = GGA.project_date_bins()
    Δd = 30
    date_range = reverse(last(date_range):-Day(Δd):Date(1970, 1, 1))
    date_center = date_range[1:end-1] .+ Day(Δd / 2)

    length(date_range) == length(date_center) + 1 ||
        error("date_range and date_center lengths are inconsistent")

    # Every geotile holding ice is processed. Those the sweep covers are read from their own tile; the
    # rest borrow the nearest one and are rebinned onto their own hypsometry by `process_gemb_tiles`,
    # which reports each substitution. Only tiles spanning the whole date axis are eligible, as donor or
    # otherwise: a partial one would arrive partly NaN and fail the checks below with no indication of
    # why, and the directory also accumulates shorter development sweeps.
    coverage = GGA.gemb_tile_coverage(gembinfo.tile_dir)
    available = keys(GGA.gemb_tile_paths(gembinfo.tile_dir;
                                        covering=(first(date_center), last(date_center))))
    println("tile files: $(length(coverage)), of which $(length(available)) span the date axis")
    println("dates: $(length(date_center)) bins, $(first(date_center)) .. $(last(date_center))")
end;

for surface_mask in gemb_surface_masks
    geotiles = GGA._geotile_load_align(; surface_mask, geotile_order=nothing,
                                      only_geotiles_w_area_gt_0=true)
    area_km2 = GGA._geotile_area_km2(; surface_mask, geotile_width)

    println("$(surface_mask): geotiles with a tile of their own: $(sum(in(id, available) for id in geotiles.id)) of $(nrow(geotiles))")

    if !isnothing(rgi_subset)
        GGA.add_single_rgi_column!(geotiles)
        geotiles = geotiles[geotiles.rgi.==rgi_subset, :]
    end

    if !isnothing(single_geotile_test)
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)], OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!"
        geotiles = geotiles[geotiles.id.==single_geotile_test, :]
    end

    nrow(geotiles) == 0 && error("$(surface_mask): no geotiles left to process")
    area_km2 = area_km2[geotile=At(geotiles.id)]

    println("$(surface_mask): geotiles to process: $(nrow(geotiles))  ($(round(Int, sum(area_km2))) km² of ice)")

    gemb_dv0 = GGA.process_gemb_tiles(geotiles, area_km2; gembinfo.tile_dir, date_center)

    gemb_dv0 = GGA.gemb_add_derived_vars!(gemb_dv0)

    # The same gates the `.mat` path asserts. They are cheap and they have caught real problems.
    for k in keys(gemb_dv0)
        if any(isnan.(gemb_dv0[k][date=Near(Date(2000, 1, 1))]))
            error("$(surface_mask): NaN found in $(k) for nearest Date(2000,1,1)")
        end
    end

    for k in [:refreeze, :rain, :acc, :melt, :runoff]
        if any(diff(gemb_dv0[k], dims=:date) .< -1E-9)
            error("$(surface_mask): negative values found in $(k)")
        end
    end

    if any(diff(gemb_dv0[:refreeze], dims=:date) .- diff(gemb_dv0[:melt], dims=:date) .> 1E-9)
        error("$(surface_mask): refreeze is greater than melt")
    end

    println("$(surface_mask): all checks passed")

    if isnothing(single_geotile_test)
        gemb_geotile_filename_dv = GGA.gemb_geotile_dv_file(gembinfo, surface_mask)
        save(gemb_geotile_filename_dv, Dict("gemb_dv" => gemb_dv0, "surface_mask" => String(surface_mask)))
        println("wrote $(gemb_geotile_filename_dv)")
    end
end
