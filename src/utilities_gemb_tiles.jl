# Read GEMB_GlacierSims tile NetCDF output into the geotile volume-change ensemble the synthesis
# consumes.
#
# The tile files carry band-resolved output on the same 100 m grid as `project_height_bins`, over a
# (ΔT × precipitation scaling) forcing matrix, for most glacierized 2° tiles. That makes the
# hypsometric gap-filling in `process_gemb_geotiles` unnecessary: there is no simulated elevation
# class, and the only interpolation across the forcing grid is the one `gemb_dv_sample` already
# performs downstream.
#
# The sweep does not reach every geotile this package finds ice in, so a geotile with no tile file of
# its own borrows the nearest one ([`gemb_tile_donors`](@ref)) and is rebinned onto its own hypsometry.
# The borrowed series is a worse estimate the further the donor is and the less its band range overlaps
# the target's ice, which is why `process_gemb_tiles` reports both per substitution.
#
# Bands are weighted by this package's own `_geotile_area_km2`, not the tile file's `band_area`. The two
# are different inventories on the same grid and disagree by a few percent per tile; the altimetry these
# series are calibrated against is integrated over `_geotile_area_km2`, so using it for both is what
# keeps the comparison free of an area bias.

# The ice density a tile was run with, kg m-3, read from the file rather than assumed. Mass fluxes are
# stored in kg m-2 and the height layers in metres, so a mass layer crosses this once on its way to a
# volume -- and it must be the same value GEMB divided by to build `dh_mass`, or the budget will not
# close. Read rather than taken from `δice` so a tile run at a non-default density still converts
# correctly; the two agree at GEMB's default of 917.
function _gemb_density_ice(ds)
    haskey(ds.attrib, "model_density_ice") ||
        error("tile file has no model_density_ice attribute; cannot convert mass fluxes to volume " *
              "consistently with its height layers")
    return Float64(ds.attrib["model_density_ice"])
end

# km³ per (km² × m).
const _GEMB_KM3_PER_KM2_M = 1e-3

# Epoch the tile time axes are reduced to before resampling, so a file written against a different
# `days since` origin still lands on the same axis.
const _GEMB_TIME_EPOCH = DateTime(1950, 1, 1)

# Layers `process_gemb_tiles` builds directly. `smb`, `runoff` and `dv` are derived from these by
# `gemb_add_derived_vars!`, which the driver applies before saving so the file matches the `.mat`
# path's contents.
const GEMB_TILE_LAYERS = (:fac, :acc, :refreeze, :melt, :rain, :ec)

# Band layers read from each tile file. The `dh` layers are `-cumsum(ice_flux)` and already cumulative;
# the mass fluxes are per-interval sums in kg m-2.
const _GEMB_BAND_HEIGHT_LAYERS = ("dh", "dh_mass")
const _GEMB_BAND_MASS_LAYERS = ("melt", "runoff", "refreeze", "rain", "precipitation",
                                "evaporation_condensation")

"""
    gemb_tile_coverage(tile_dir) -> Dict{String,@NamedTuple{path::String, first::DateTime, last::DateTime, n::Int}}

Geotile id, path and time span of every `.nc` in `tile_dir`.

Ids come from each file's `geotile_id` attribute, which `GEMB_GlacierSims` already writes in this
package's own form (`"lat[+62+64]lon[-152-150]"`), so no filename parsing is involved. Only the time
axis and the attributes are read.

The directory accumulates runs over different windows and output intervals -- a tile from a short
development sweep sits beside a full-record one -- so a caller must select on the span rather than
assume every file is usable. [`gemb_tile_paths`](@ref) does that selection.

Throws if a file carries no `geotile_id`, or if two files claim the same geotile.
"""
function gemb_tile_coverage(tile_dir)
    isdir(tile_dir) || error("no GEMB tile directory at $tile_dir")

    coverage = Dict{String,@NamedTuple{path::String, first::DateTime, last::DateTime, n::Int}}()
    for filename in readdir(tile_dir)
        endswith(filename, ".nc") || continue
        path = joinpath(tile_dir, filename)

        entry = NCDataset(path, "r") do ds
            geotile_id = get(ds.attrib, "geotile_id", nothing)
            isnothing(geotile_id) && error("no geotile_id attribute in $path")
            time = DateTime.(ds["time"][:])
            (geotile_id, (; path, first=first(time), last=last(time), n=length(time)))
        end

        geotile_id, span = entry
        if haskey(coverage, geotile_id)
            error("two tile files claim geotile $geotile_id: $(coverage[geotile_id].path) and $path")
        end
        coverage[geotile_id] = span
    end

    isempty(coverage) && error("no .nc tile files found in $tile_dir")
    return coverage
end

"""
    gemb_tile_paths(tile_dir; covering=nothing) -> Dict{String,String}

Map geotile id to tile file path, keeping only tiles whose time axis spans `covering`.

`covering` is a pair or tuple of the first and last date the ensemble needs. Tiles that do not reach
both ends are dropped: the ensemble is built on one date axis, so a tile covering part of it could only
contribute NaN over the rest. Pass `nothing` to keep every tile.

See [`gemb_tile_coverage`](@ref) for why this filter is not optional in practice.
"""
function gemb_tile_paths(tile_dir; covering=nothing)
    coverage = gemb_tile_coverage(tile_dir)
    isnothing(covering) && return Dict(id => span.path for (id, span) in coverage)

    lo, hi = DateTime(first(covering)), DateTime(last(covering))
    return Dict(id => span.path for (id, span) in coverage
                if span.first <= lo && span.last >= hi)
end

# Geotile center as (longitude, latitude), the order `Distances.Haversine` takes.
function _geotile_center(geotile_id)
    extent = geotile_extent(geotile_id)
    return ((extent.X[1] + extent.X[2]) / 2, (extent.Y[1] + extent.Y[2]) / 2)
end

"""
    gemb_tile_donors(geotile_ids, paths; max_distance_km=2000) -> Dict{String,@NamedTuple{donor::String, distance_km::Float64}}

For each of `geotile_ids`, the geotile whose tile file supplies its GEMB series.

A geotile present in `paths` is its own donor at zero distance. One that is not borrows the tile whose
geotile center is nearest by great-circle distance; the borrowed series is then rebinned onto the
target's own hypsometry by [`_gemb_band_weights`](@ref), which interpolates across the donor's bands in
height and holds the outermost band constant beyond them.

A donor is only defensible where the two geotiles share a climate and an elevation range, and neither
condition is checked here -- distance is the whole criterion. `max_distance_km` is the backstop: it is
loose enough for the isolated sub-polar islands the sweep omits, whose nearest neighbor is over a
thousand kilometres away, and tight enough that a sweep which dropped a whole region throws rather than
smearing one surviving tile across it.

Throws if `paths` is empty or if any geotile's nearest donor is further than `max_distance_km`.
"""
function gemb_tile_donors(geotile_ids, paths; max_distance_km=2000)
    isempty(paths) && error("no GEMB tile files to draw donors from")

    donor_ids = collect(keys(paths))
    donor_centers = _geotile_center.(donor_ids)
    # Kilometres, so `max_distance_km` and the reported distances share units.
    haversine = Haversine(6371.0)

    donors = Dict{String,@NamedTuple{donor::String, distance_km::Float64}}()
    too_far = Tuple{String,String,Float64}[]

    for geotile_id in geotile_ids
        if haskey(paths, geotile_id)
            donors[geotile_id] = (; donor=geotile_id, distance_km=0.0)
            continue
        end

        center = _geotile_center(geotile_id)
        distance_km, j = findmin(donor_center -> haversine(center, donor_center), donor_centers)

        if distance_km > max_distance_km
            push!(too_far, (geotile_id, donor_ids[j], distance_km))
        else
            donors[geotile_id] = (; donor=donor_ids[j], distance_km)
        end
    end

    if !isempty(too_far)
        report = join(("$(id) (nearest $(donor) at $(round(Int, d)) km)" for (id, donor, d) in too_far),
                      ", ")
        error("$(length(too_far)) geotiles have no GEMB tile within $(max_distance_km) km: $report")
    end

    return donors
end

"""
    gemb_tile_forcing_grid(path) -> (; pscale, ΔT)

The precipitation-scaling and temperature-offset axes of a tile file, as ascending vectors.

`gemb_dv_sample` interpolates over both with `Gridded(Linear())`, which requires ascending knots, so
the ensemble is built on sorted axes and the tile data is permuted into that order on read.
"""
function gemb_tile_forcing_grid(path)
    return NCDataset(path, "r") do ds
        (; pscale = sort(Float64.(ds["precipitation_scaling"][:])),
           ΔT = sort(Float64.(ds["delta_temperature"][:])))
    end
end

"""
    _gemb_band_weights(band_centers, area_km2) -> (w, extrapolated_fraction)

Collapse "fill the ice-bearing height bins GEMB has no band for, weight by glacier area, sum over
height" into one weight per band, so an area-weighted total becomes a single `series * w`.

`area_km2` is one geotile's hypsometry over `project_height_bins` centers. A bin holding ice that no
band covers is filled by linear interpolation in height across the bands that do, held constant beyond
the outermost band -- the rule `process_gemb_geotiles` applies. That fill is linear in the band values
and the area weighting is linear in the filled values, so the composition is a fixed linear functional
of the band series and `w` is that functional.

`w` is in km³ per metre of ice equivalent. `extrapolated_fraction` is the share of the geotile's ice
area lying outside the band range, which a caller should report rather than silently accept.
"""
function _gemb_band_weights(band_centers, area_km2)
    issorted(band_centers) || error("band centers are not ascending: $(band_centers)")

    height_centers = collect(val(dims(area_km2, :height)))
    areas = Float64.(collect(area_km2))
    ice = findall(>(0), areas)
    isempty(ice) && error("geotile hypsometry holds no ice")

    n_band = length(band_centers)
    w = zeros(Float64, n_band)
    extrapolated_area = 0.0

    for i in ice
        height = height_centers[i]
        weight = areas[i] * _GEMB_KM3_PER_KM2_M

        j = searchsortedfirst(band_centers, height)
        if j <= n_band && band_centers[j] == height
            # A band covers this bin exactly, which holds for all but a handful of bins.
            w[j] += weight
        elseif j == 1 || j > n_band
            # Outside the band range: hold the nearest band constant.
            w[clamp(j, 1, n_band)] += weight
            extrapolated_area += areas[i]
        else
            lo, hi = band_centers[j-1], band_centers[j]
            f = (height - lo) / (hi - lo)
            w[j-1] += weight * (1 - f)
            w[j] += weight * f
        end
    end

    return w, extrapolated_area / sum(areas[ice])
end

# Cumulative levels resampled onto `target_days`. Linear interpolation of the level is the right
# reduction: the increment over a target interval is then the integral of the flux across it.
function _gemb_resample(series, source_days, target_days)
    return DataInterpolations.LinearInterpolation(series, source_days).(target_days)
end

"""
    read_gemb_tile(path, area_km2, target_days; n_pscale, n_ΔT)

One tile's volume-change ensemble as `(; layers, extrapolated_fraction, closure_km3, melt_offset_km3)`.

`layers` is a `Dict{Symbol,Array{Float64,3}}` over `(date, pscale, ΔT)` in km³ of ice equivalent,
cumulative and referenced to the first target date. `area_km2` is the geotile's row of
`_geotile_area_km2`; `target_days` are the output date centers as days since `_GEMB_TIME_EPOCH`.

# Layer construction

`dv` and `dv_mass` are the area-weighted band `dh` and `dh_mass`. Every mass flux is split by where it
ended up before being accumulated, over each output interval:

    rain_retained = max(rain - runoff, 0)     -- rain that never left, so it refroze in the ice
    rain_shed     = rain - rain_retained      -- rain that ran off
    melt_runoff   = max(runoff - rain, 0)     -- the part of runoff that is meltwater

giving the stored layers

    acc      = cumulative (precipitation - rain + rain_retained)   -- snowfall plus refrozen rain
    rain     = cumulative rain_shed
    refreeze = cumulative refreeze                                 -- GEMB's, unaltered
    melt     = cumulative (melt_runoff + refreeze)
    fac      = dv - dv_mass
    ec       = -cumulative evaporation_condensation

This matches the convention the altimetry products are on: **accumulation is snowfall plus the rain that
refroze**, and **runoff excludes rain**, so `gemb_add_derived_vars!`'s `runoff = melt - refreeze` comes out
as `melt_runoff` alone. Neither can be read straight off the file. GEMB's `runoff` carries rain as well
as meltwater despite its `long_name`, which the water balance settles:
`melt + rain - runoff - refreeze - Δstored` closes while the rain-free version misses by exactly the
rain. GEMB's `melt` is separately net of rain (`melt_total = max(0, melt_sum - rain)` in
`calculate_melt.jl`); `melt_offset_km3` reports how far the layer above sits from it, which is a percent
or so.

The split loses nothing. `max(0, a-b) - max(0, b-a) == a - b` for any interval, so the retained and shed
halves cancel in `smb`, which reduces to `precipitation - runoff + evaporation_condensation` -- GEMB's
own surface mass balance. Hence

  * `smb == dv_mass`, and `ec` is the real evaporation and condensation rather than a closure residual,
  * `acc`, `rain`, `refreeze`, `melt` and the derived `runoff` are all monotonic by construction, with
    nothing clamped away,
  * `dv == Σ dh·area`, and `dm = (dv - fac)·δice` is exactly GEMB's mass term.

`closure_km3` is the largest `|smb - dv_mass|` over the forcing grid and should sit at the rounding
floor; a nonzero value means one of the assumptions above no longer holds for that tile.

!!! note "Ice density"
    Mass fluxes are divided by the tile's own `model_density_ice`, which is what GEMB divided by to build
    `dh_mass`. Using any other value leaves a closure error of order a percent of the accumulated fluxes
    -- tens of km³ on a large tile. `δice` is set to GEMB's default of 917 so that the volumes here and
    the `volume2mass = δice/1000` applied downstream describe the same ice.
"""
function read_gemb_tile(path, area_km2, target_days; n_pscale, n_ΔT)
    return NCDataset(path, "r"; maskingvalue = NaN) do ds
        band_centers = Float64.(ds["band_center"][:])
        w, extrapolated_fraction = _gemb_band_weights(band_centers, area_km2)
        w_mass = w ./ _gemb_density_ice(ds)

        source_days = Dates.value.(DateTime.(ds["time"][:]) .- _GEMB_TIME_EPOCH) ./ 86_400_000
        if first(source_days) > first(target_days) || last(source_days) < last(target_days)
            error("$(basename(path)) spans $(first(source_days))..$(last(source_days)) days since " *
                  "$(_GEMB_TIME_EPOCH) but $(first(target_days))..$(last(target_days)) was requested")
        end

        # The ensemble axes are ascending; the file's are not necessarily, so read through a permutation.
        iΔT = sortperm(Float64.(ds["delta_temperature"][:]))
        ipscale = sortperm(Float64.(ds["precipitation_scaling"][:]))
        length(iΔT) == n_ΔT && length(ipscale) == n_pscale ||
            error("$(basename(path)) forcing grid is $(length(iΔT))x$(length(ipscale)), expected " *
                  "$(n_ΔT)x$(n_pscale)")

        # Read each band layer once, whole. Slicing a deflated variable per forcing node instead costs
        # one chunk decompression per slice and dominates the runtime -- 49 nodes over 7 layers is two
        # orders of magnitude more decompression than the file needs. Kept at the stored `Float32`:
        # about 45 MB per layer, and the arithmetic below promotes against the `Float64` weights.
        band = Dict(name => ds[name][:, :, :, :] for name in
                    (_GEMB_BAND_HEIGHT_LAYERS..., _GEMB_BAND_MASS_LAYERS...))
        n_time = length(target_days)
        layers = Dict(k => fill(NaN, n_time, n_pscale, n_ΔT)
                      for k in (GEMB_TILE_LAYERS..., :dv, :dv_mass))
        closure_km3 = 0.0
        melt_offset_km3 = 0.0

        for (jp, p) in enumerate(ipscale), (jt, t) in enumerate(iΔT)
            # `w` folds the height fill and the area weighting into one vector, so each layer's
            # area-weighted total is a single matrix-vector product over bands.
            total(name, weights) = view(band[name], :, :, t, p) * weights

            precipitation = total("precipitation", w_mass)
            rain = total("rain", w_mass)
            runoff = total("runoff", w_mass)
            refreeze = total("refreeze", w_mass)

            # Split the rain by where it ended up. GEMB's `runoff` carries rain as well as meltwater
            # despite its name, so rain in excess of what left the column over an interval must have
            # stayed in the ice: that part is accumulation, and the remainder ran off. Both halves are
            # per-interval quantities, hence the split before any cumulative sum.
            rain_retained = max.(rain .- runoff, 0.0)
            rain_shed = rain .- rain_retained
            melt_runoff = max.(runoff .- rain, 0.0)

            native = Dict{Symbol,Vector{Float64}}(
                :dv => total("dh", w),
                :dv_mass => total("dh_mass", w),
                # Snowfall plus the rain that refroze in the ice.
                :acc => cumsum(precipitation .- rain .+ rain_retained),
                # Only the rain that left; the retained part is in `acc`, not counted twice.
                :rain => cumsum(rain_shed),
                :refreeze => cumsum(refreeze),
                # So `gemb_add_derived_vars!`'s `melt - refreeze` is runoff net of rain.
                :melt => cumsum(melt_runoff .+ refreeze),
                # Positive `ec` is mass loss here, the reverse of GEMB's flux convention.
                :ec => -cumsum(total("evaporation_condensation", w_mass)),
                # GEMB's own melt layer, kept only to report how far the above sits from it.
                :melt_gemb => cumsum(total("melt", w_mass)),
            )

            series = Dict(k => _gemb_resample(v, source_days, target_days) for (k, v) in native)
            for v in values(series)
                v .-= first(v)
            end

            series[:fac] = series[:dv] .- series[:dv_mass]
            smb = series[:acc] .- series[:melt] .+ series[:refreeze] .- series[:ec]
            closure_km3 = max(closure_km3, maximum(abs, smb .- series[:dv_mass]))
            melt_offset_km3 = max(melt_offset_km3, abs(last(series[:melt]) - last(series[:melt_gemb])))

            for k in (GEMB_TILE_LAYERS..., :dv, :dv_mass)
                layers[k][:, jp, jt] = series[k]
            end
        end

        return (; layers, extrapolated_fraction, closure_km3, melt_offset_km3)
    end
end

"""
    process_gemb_tiles(geotiles, area_km2; tile_dir, date_center, show_stats=true, max_donor_distance_km=2000)

Build the GEMB volume-change ensemble from `GEMB_GlacierSims` tile NetCDF, as a drop-in for
`process_gemb_geotiles`.

Returns a `DimStack` of [`GEMB_TILE_LAYERS`](@ref) over `(geotile, date, pscale, ΔT)` in km³ of ice
equivalent. `smb`, `runoff` and `dv` are added later by `gemb_add_derived_vars!`, as with the `.mat`
path.

Every band is weighted by `area_km2`, so the ensemble and the altimetry it is calibrated against share
one hypsometry. A geotile with no tile file of its own borrows the nearest one within
`max_donor_distance_km` -- see [`gemb_tile_donors`](@ref) -- and is rebinned onto its own hypsometry, so
every geotile passed in gets a series rather than an all-NaN row that would reach the ensemble's NaN
check with no indication of why.

`show_stats` prints the substitutions worst-area first, the geotiles whose ice extends beyond the band
range they were read against, and the largest `ec` closure residual. All three are silent data-quality
problems otherwise, and the first two compound: a borrowed tile whose bands miss the target's ice
entirely is held constant from its nearest band, which is the weakest estimate this produces.
"""
function process_gemb_tiles(geotiles, area_km2; tile_dir, date_center, show_stats=true,
                            max_donor_distance_km=2000)
    paths = gemb_tile_paths(tile_dir; covering=extrema(date_center))
    donors = gemb_tile_donors(geotiles.id, paths; max_distance_km=max_donor_distance_km)

    grid = gemb_tile_forcing_grid(paths[donors[first(geotiles.id)].donor])
    ddate = Dim{:date}(collect(date_center))
    dgeotile = Dim{:geotile}(collect(geotiles.id))
    dpscale = Dim{:pscale}(grid.pscale)
    dΔT = Dim{:ΔT}(grid.ΔT)

    gemb_dv = DimStack([DimArray(fill(NaN, (dgeotile, ddate, dpscale, dΔT)); name=k)
                        for k in GEMB_TILE_LAYERS]...)

    target_days = Dates.value.(DateTime.(date_center) .- _GEMB_TIME_EPOCH) ./ 86_400_000
    extrapolated = fill(NaN, dgeotile)
    melt_offset = fill(NaN, dgeotile)
    closure = fill(NaN, dgeotile)

    @showprogress desc = "Reading GEMB tiles" Threads.@threads for row in collect(eachrow(geotiles))
        tile = read_gemb_tile(paths[donors[row.id].donor], area_km2[geotile=At(row.id)], target_days;
                              n_pscale=length(grid.pscale), n_ΔT=length(grid.ΔT))

        for k in GEMB_TILE_LAYERS
            gemb_dv[k][geotile=At(row.id)] = tile.layers[k]
        end

        extrapolated[At(row.id)] = tile.extrapolated_fraction
        melt_offset[At(row.id)] = tile.melt_offset_km3
        # `dv` reconstructed through this package's identities against the one the tile file implies.
        smb = tile.layers[:acc] .- tile.layers[:melt] .+ tile.layers[:refreeze] .- tile.layers[:ec]
        closure[At(row.id)] = maximum(abs, (smb .+ tile.layers[:fac]) .- tile.layers[:dv])
    end

    if show_stats
        printstyled("GEMB tiles read: $(length(dgeotile))\n"; color=:light_green)
        @printf("  max |dv| closure error : %.3g km3\n", maximum(closure))
        @printf("  melt vs GEMB's melt    : up to %.1f km3 (%s) -- the rain GEMB nets off its melt\n",
                maximum(melt_offset), val(dgeotile)[argmax(melt_offset)])

        borrowed = [id for id in val(dgeotile) if donors[id].donor != id]
        if isempty(borrowed)
            println("  every geotile has a GEMB tile of its own")
        else
            ice_km2 = Dict(id => sum(area_km2[geotile=At(id)]) for id in val(dgeotile))
            sort!(borrowed; by=id -> ice_km2[id], rev=true)
            total_km2 = sum(values(ice_km2))
            @printf("  %d geotiles have no tile of their own and borrow the nearest (%.0f km2, %.3f%% of the ice):\n",
                    length(borrowed), sum(ice_km2[id] for id in borrowed),
                    100 * sum(ice_km2[id] for id in borrowed) / total_km2)
            for id in first(borrowed, 10)
                @printf("      %-28s %8.2f km2 <- %-28s %5.0f km, %5.1f%% of ice beyond its bands\n",
                        id, ice_km2[id], donors[id].donor, donors[id].distance_km,
                        100 * extrapolated[At(id)])
            end
            length(borrowed) > 10 && println("      ... $(length(borrowed) - 10) more, all smaller")
        end

        beyond = findall(>(0), collect(extrapolated))
        if isempty(beyond)
            println("  every ice-bearing height bin is covered by a GEMB band")
        else
            order = beyond[sortperm(collect(extrapolated)[beyond]; rev=true)]
            @printf("  %d geotiles have ice beyond the band range (median %.2f%%):\n",
                    length(beyond), 100 * median(collect(extrapolated)[beyond]))
            for i in first(order, 10)
                @printf("      %-28s %6.2f%% of ice area\n",
                        val(dgeotile)[i], 100 * extrapolated[i])
            end
        end
    end

    return gemb_dv
end
