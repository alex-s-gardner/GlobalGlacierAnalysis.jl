#density of glacier ice
# Ice density [kg m-3]. Matches the `density_ice` the GEMB runs use, which the tile files record as
# their `model_density_ice` attribute and divide by to build the height-change decomposition. The two
# have to agree: a GEMB volume converted to mass with a different density carries that ratio as a bias.
const δice = 917;
const local2utc = Hour(7) # LA timezone to UTC
const seasonality_weight = 85/100
# Strength of the forcing-prior penalty (`origin_penalty_mode = :prior`): the cost is multiplied by
# 1 + wd × d, d the Mahalanobis distance from `gemb_forcing_prior`. 0.35 is the largest value within 1% of the
# best held-out skill in temporal cross-validation; see notes/methods_2026-10_seasonal_cycle_and_calibration.md.
const distance_from_origin_penalty = 35 / 100
# Only used by the `:legacy` penalty mode.
const ΔT_to_pscale_weight = 50/100

# Empirical prior on the GEMB forcing corrections for `origin_penalty_mode = :prior`: the area-weighted
# centre, spread and correlation of the unpenalized fits of the well-constrained geotile groups (sd of
# log pscale < 0.15 and of ΔT < 0.5 K within 5% of the minimum cost; 88% of ice area) for the reference
# ensemble member, fill set 6, GEMB run 8. See notes/methods_2026-10_seasonal_cycle_and_calibration.md.
const gemb_forcing_prior = (pscale=1.60, log_pscale_sd=0.45, ΔT=2.20, ΔT_sd=1.65, corr=0.50)
const ocean_area_km2 = 362.5 * 1E6
const reference_ensemble_file = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_rgi7_dh_cop30_v2_cc_nmad5_v01_filled_ac_p2_aligned.jld2"; 

"""
    project_products(; project_id = :v01)

Get the elevation products configuration for a specific project.

# Arguments
- `project_id::Symbol`: Project identifier (default: :v01)

# Returns
- Named tuple containing configured ElevationProduct instances for different missions
- NOTE: latitude_limits are the latitude limits of the mission's data product and need to be specified as Integers

# Examples
```julia
julia> product = project_products(; project_id=:v01)
julia> icesat2_cfg = product.icesat2
```
"""
function project_products(; project_id = :v01)
    if project_id == :v01
        product = (
            # ATL06 v7. NSIDC retired v006 from CMR when v007 was released, so a v6 archive can no
            # longer be searched or downloaded -- `search(:ICESat2, :ATL06; version=6)` returns
            # nothing from either NSIDC_CPRD or the (now decommissioned) NSIDC_ECS provider. The
            # version is part of the data path, so v7 builds into `icesat2/ATL06/007/` and leaves
            # the existing v6 raw granules in `006/raw` untouched. That directory is now the only
            # copy of the v6 data used for the published analysis: do not delete it.
            icesat2=(mission=:icesat2, name=:ATL06, version=7, id="I206", error_sigma=0.1, halfwidth=11 / 4, kernel=:gaussian, apply_quality_filter=false, coregister=true, latitude_limits=[-88, 88], longitude_limits=[-180, 180]),

            icesat=(mission=:icesat, name=:GLAH06, version=34, id="I106", error_sigma=0.1, halfwidth=35 / 4, kernel=:gaussian, apply_quality_filter=false, coregister=true, latitude_limits=[-86, 86], longitude_limits=[-180, 180]),

            gedi=(mission=:gedi, name=:GEDI02_A, version=2, id="G02A", error_sigma=0.1, halfwidth=22 / 4, kernel=:gaussian, apply_quality_filter=false, coregister=true, latitude_limits=[-52, 52], longitude_limits=[-180, 180]),

            hugonnet=(mission=:hugonnet, name=:HSTACK, version=1, id="HS01", error_sigma=5, halfwidth=100 / 2, kernel=:gaussian, apply_quality_filter=false, coregister=true, latitude_limits=[-90, 90], longitude_limits=[-180, 180])
        )
    end
    return product
end


"""
    project_paths(; project_id = :v01)

Get the file paths configuration for a specific project.

# Arguments
- `project_id::Symbol`: Project identifier (default: :v01)

# Returns
- Named tuple containing configured paths for different data products

# Examples
```julia
julia> paths = project_paths(; project_id=:v01)
julia> hugonnet_path = paths.hugonnet
```
"""
function project_paths(; project_id = :v01)
    if project_id == :v01
        geotile_width = 2 #geotile width [degrees]

        p = project_products(project_id = project_id)

        # One entry per registered product, so a new mission needs adding only to project_products.
        paths = map(v -> setpaths(geotile_width, v.mission, "$(v.name)", lpad("$(v.version)", 3, '0')), p)
    end
    return paths
end

"""
    project_geotiles(; geotile_width = 2, domain = :all, extent=nothing)

Define geotiles for the project, optionally filtered by domain.

# Arguments
- `geotile_width::Int`: Width of geotiles in degrees (default: 2)
- `domain::Symbol`: Domain filter (:all or :landice) (default: :all)
- `extent`: Optional bounding box to limit geotile creation

# Returns
- DataFrame of geotiles with their properties

# Examples
```julia
julia> geotiles = project_geotiles(; geotile_width=2, domain=:all)
julia> geotiles_ice = project_geotiles(; geotile_width=2, domain=:landice)
```
"""
function project_geotiles(; geotile_width = 2, domain = :all, extent=nothing)
    paths = setpaths()
    geotiles = GeoTiles.define(geotile_width; extent)

    if domain == :all

    elseif domain == :landice
        icemaskfn = paths.icemask
        icemask = GeoArrays.read(icemaskfn);
        geotilemask = GeoArrays.crop.(Ref(icemask), extent2nt.(geotiles.extent))
        #geotilemask = GeoArrays.crop.(Ref(icemask), geotiles.extent)
        hasice = [any(mask.A.==1) for mask in geotilemask];
        geotiles = geotiles[hasice,:];
    else
        error("unrecognized domain = $domain")
    end
    return geotiles
end

"""
    analysis_paths(; geotile_width = 2)

Create and return paths for analysis outputs.

# Arguments
- `geotile_width::Int`: Width of geotiles in degrees (default: 2)

# Returns
- Named tuple containing paths for analysis outputs

# Examples
```julia
julia> paths = analysis_paths(; geotile_width=2)
julia> binned_path = paths.binned
```
"""
function analysis_paths(; geotile_width = 2)
    paths = (
        binned = joinpath(setpaths().data_dir, "binned", "$(geotile_width)deg"),
    )
    
    for p in paths
        if !isdir(p)
            mkpath(p)
        end
    end
    return paths
end


"""
    project_date_bins()

Define temporal bins for the project.

The end date is the one place the length of the record is set; [`project_decyear_bins`](@ref) derives
its own extent from this, so the two cannot drift apart. Extend it when a mission's archive grows past
the last bin, and re-run the binning: everything downstream inherits the grid from here.

The bins run to `2026-07-01`, which covers ICESat-2 ATL06 v7 (to 2026-05-18) with one bin to spare.
Other inputs stop earlier -- GEDI in 2025-07, Hugonnet in 2019-10, ICESat in 2009-10, and the GEMB
runs in 2024 -- and simply hold no data in later bins, which the synthesis already handles.

Binned products carry `date_center` as their `:date` dimension. `geotile_binning` builds that dimension
only when the output file does not yet exist, so a product updated in place keeps the dimension it was
first written with; rebuild every mission (`missions2update=nothing`) whenever this grid changes.

# Returns
- Tuple containing (date_range, date_center) where:
  - date_range: `Date` range of bin edges with 30-day intervals
  - date_center: `Date` values at the center of each bin, one per bin

# Examples
```julia
julia> date_range, date_center = project_date_bins()
julia> ddate = Dim{:date}(date_center)
```
"""
function project_date_bins()
        Δd = 30
        date_range = Date(1990):Day(Δd):Date(2026, 7, 1)
        date_center = date_range[1:end-1] .+ Day(Δd / 2)

    return date_range, date_center
end

"""
    project_decyear_bins()

Bin edges used to group observations by date, in decimal years.

Same count as [`project_date_bins`](@ref)'s `date_range`, so a binned array always has one date bin per
`date_center`. Taking only the count from there, rather than converting the dates, is deliberate: the
edges keep their original `1990 + k * 30/365` spacing, so extending the record appends bins without
moving any existing one.

Note that `30/365` is not thirty calendar days, so the window an observation is grouped into runs late
relative to the date it is labelled with. The offset accumulates: zero in 1990, 5 days by 2009, and
9 days by 2025 -- largest over exactly the years with the most data. Deriving these edges from the
dates instead would move about 30% of observations into an adjacent bin (41% for 2025 onward), which
changes published values, so the spacing is left as it is.
"""
function project_decyear_bins()
    date_range, _ = project_date_bins()
    return range(1990.0; step=30 / 365, length=length(date_range))
end

"""
    project_height_bins()

Define elevation bins for the project.

# Returns
- Tuple containing (height_range, height_center) where:
  - height_range: Range of elevation bin edges from 0 to 10000m at 100m intervals
  - height_center: Values at the center of each elevation bin

# Examples
```julia
julia> height_range, height_center = project_height_bins()
julia> dheight = Dim{:height}(height_center)
```
"""
function project_height_bins()
    Δh = 100;
    height_range = 0:100:10000;
    height_center = height_range[1:end-1] .+ Δh / 2;

    return height_range, height_center
end

# An alternative binning must be a distinctly named function, not a second `project_mscale_bins()`
# method: a duplicate empty signature is silently shadowed, and method overwriting makes the module
# unprecompilable, forcing a from-source rebuild on every load.

"""
    project_mscale_bins()

Define melt-scaling bins (log-style spacing).

Used only by `process_gemb_geotiles`, whose forcing axis is a melt multiplier centred on 1. The tile
path's axis is a temperature offset in K centred on 0 and takes its values from the data, so it does not
come through here.

# Returns
- Tuple containing (mscale_range, mscale_center) with values [1/6, 1/4, 1/2, 2, 4, 6] and centers.

# Examples
```julia
julia> mscale_range, mscale_center = project_mscale_bins()
```
"""
function project_mscale_bins()

    mscale_range = [1/6, 1/4, 1/2, 2, 4, 6]
    mscale_center = [1/5, 1/3, 1, 3, 5]

    return mscale_range, mscale_center
end



"""
    mission_land_trend()

Get the land elevation trend correction for each mission.

# Returns
- DimensionalArray with trend values (m/yr) for each mission

# Examples
```julia
julia> mission_trend = mission_land_trend()
julia> gedi_trend = mission_trend[At("gedi")]
```
"""
function mission_land_trend()
    missions0 = ["icesat", "icesat2", "gedi", "hugonnet"]
    dmission = Dim{:mission}(missions0)
    mission_trend_myr = fill(0.0, dmission)
    mission_trend_myr[At(["gedi"])] .= -0.144

    return mission_trend_myr
end

const geotiles_golden_test = [
    "lat[+30+32]lon[+078+080]",
    "lat[+60+62]lon[-142-140]", 
    "lat[+62+64]lon[-052-050]", 
    "lat[-68-66]lon[-070-068]", 
    "lat[-44-42]lon[-074-072]", 
    "lat[-34-32]lon[-070-068]",
    "lat[-74-72]lon[-080-078]",
    "lat[+34+36]lon[+076+078]",
    "lat[+40+42]lon[+078+080]",
    "lat[+46+48]lon[+008+010]",
    "lat[+76+78]lon[+016+018]",
    "lat[+60+62]lon[+006+008]",
    "lat[+64+66]lon[-018-016]",
    "lat[+66+68]lon[-052-050]",
    "lat[+78+80]lon[-076-074]",
    "lat[+68+70]lon[-070-068]",
    "lat[+56+58]lon[-134-132]",
    "lat[-48-46]lon[-074-072]",
    "lat[-34-32]lon[-070-068]",
    "lat[+64+66]lon[+058+060]",
    "lat[+50+52]lon[-126-124]"
    ]

"""
    gemb_info(; gemb_run_id=4)

Get GEMB (Glacier Energy and Mass Balance) model configuration for a specific run.

# Arguments
- `gemb_run_id::Int`: GEMB run identifier (1-4, default: 4)

# Returns
- Named tuple containing GEMB configuration parameters:
  - `gemb_folder`: Path(s) to GEMB data folders
  - `file_uniqueid`: Unique file identifier string
  - `elevation_delta`: Array of elevation adjustments [m]
  - `precipitation_scale`: Array of precipitation scaling factors
  - `filename_gemb_combined`: Output file path for combined data

# Throws
- `ErrorException`: If gemb_run_id is not recognized

# Examples
```julia
julia> gembinfo = gemb_info(; gemb_run_id=4)
julia> dpscale = gembinfo.precipitation_scale
```
"""
function gemb_info(; gemb_run_id = 4)

    if gemb_run_id == 1
        dpscale = Dim{:pscale}(["p1"]) # do not change order as these are lookup values
        dΔheight = Dim{:Δheight}(["t1"])
        elevation_delta = DimArray([0], dΔheight) # do not change order as these are lookup values
        precipitation_scale = DimArray([1], dpscale) # do not change order as these are lookup values
        file_uniqueid = "rv1_0_19500101_20231231"
        gemb_info = (;
            gemb_folder = ["/home/schlegel/Share/GEMBv1/"],
            file_uniqueid,
            elevation_delta,
            precipitation_scale,
            filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/$file_uniqueid.jld2",
            modify_melt_only = false
        )
    elseif gemb_run_id == 2
        dpscale = Dim{:pscale}(["p1", "p2", "p3", "p4", "p5", "p6"]) # do not change order as these are lookup values
        dΔheight = Dim{:Δheight}(["t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9"])
        elevation_delta = DimArray([-1000, -750, -500, -250, 0, 250, 500, 750, 1000], dΔheight) # do not change order as these are lookup values
        precipitation_scale = DimArray([0.5, 1, 1.5, 2, 5, 10], dpscale) # do not change order as these are lookup values
        file_uniqueid = "1979to2023_820_40_racmo_grid_lwt"
        gemb_info = (;
            gemb_folder = "/home/schlegel/Share/GEMBv1/Alaska_sample/v1/",
            file_uniqueid,
            elevation_delta,
            precipitation_scale,
            filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/$file_uniqueid.jld2",
            modify_melt_only = false
        )
    elseif gemb_run_id == 3
        dpscale = Dim{:pscale}(["p1", "p2", "p3", "p4", "p5"]) # do not change order as these are lookup values
        dΔheight = Dim{:Δheight}(["t1", "t2", "t3", "t4", "t5"])
        elevation_delta = DimArray([-200, 0, 200, 500, 1000], dΔheight) # do not change order as these are lookup values
        precipitation_scale = DimArray([0.75, 1, 1.25, 1.5, 2], dpscale) # do not change order as these are lookup values

        gemb_info = (;
            gemb_folder = ["/mnt/bylot-r3/data/gemb/mat/no_lw_correction/"],
            file_uniqueid="1979to2023_820_40_",
            elevation_delta,
            precipitation_scale,
            filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0.jld2",
            modify_melt_only = false
        )
    elseif gemb_run_id == 4
        dpscale = Dim{:pscale}(["p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8"]) # do not change order as these are lookup values
        dΔheight = Dim{:Δheight}(["t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8"])
        elevation_delta = DimArray([-200, 0, 200, 500, 1000, -2000, -500, 2000], dΔheight) # do not change order as these are lookup values
        precipitation_scale = DimArray([0.75, 1, 1.25, 1.5, 2, 0.25, 3, 4], dpscale) # do not change order as these are lookup values

        gemb_info = (;
            gemb_folder=["/mnt/bylot-r3/data/gemb/mat/no_lw_correction/"],
            file_uniqueid="1979to2024_820_40_",
            elevation_delta,
            precipitation_scale,
            filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2024_820_40_racmo_grid_lwt_e97_0.jld2",
            modify_melt_only = false
        )
    elseif gemb_run_id == 5
        dpscale = Dim{:pscale}(["p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8"]) # do not change order as these are lookup values
        dΔheight = Dim{:Δheight}(["t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8"])
        elevation_delta = DimArray([-200, 0, 200, 500, 1000, -2000, -500, 2000], dΔheight) # do not change order as these are lookup values
        precipitation_scale = DimArray([0.75, 1, 1.25, 1.5, 2, 0.25, 3, 4], dpscale) # do not change order as these are lookup values

        gemb_info = (;
            gemb_folder=["/mnt/bylot-r3/data/gemb/mat/lw_correction/"],
            file_uniqueid="1979to2024_820_40_",
            elevation_delta,
            precipitation_scale,
            filename_gemb_combined="/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2024_820_40_lwt_e97_0_corrected_dmelt.jld2",
            modify_melt_only = true
        )
    elseif gemb_run_id == 6
        dpscale = Dim{:pscale}(["p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8"]) # do not change order as these are lookup values
        dΔheight = Dim{:Δheight}(["t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8"])
        elevation_delta = DimArray([-200, 0, 200, 500, 1000, -2000, -500, 2000], dΔheight) # do not change order as these are lookup values
        precipitation_scale = DimArray([0.75, 1, 1.25, 1.5, 2, 0.25, 3, 4], dpscale) # do not change order as these are lookup values

        gemb_info = (;
            gemb_folder=["/mnt/bylot-r3/data/gemb/mat/lw_correction/"],
            file_uniqueid="1979to2024_820_40_",
            elevation_delta,
            precipitation_scale,
            filename_gemb_combined="/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2024_820_40_lwt_e97_0_corrected.jld2",
            modify_melt_only=false
        )
    elseif gemb_run_id == 7
        # `GEMB_GlacierSims` tile NetCDF rather than point `.mat` files. The forcing axes live in the
        # data and are read from it by `gemb_tile_forcing_grid`, so there is nothing to declare here:
        # the tiles are complete over the height range, so `elevation_delta` -- which existed to fake
        # elevation classes from sparse point runs -- has no counterpart. `tile_dir` replaces
        # `gemb_folder`, and `read_gemb_files` is not used on this path.
        gemb_info = (;
            tile_dir=joinpath(get(ENV, "CLIMATE_CACHE", "/mnt/bylot-r3/data/era5land"), "tile_runs"),
            precipitation_scale=nothing,
            filename_gemb_combined="/mnt/bylot-r3/data/gemb/raw/gemb_glacier_sims_tiles_1950to2026.jld2",
            modify_melt_only=false
        )
    elseif gemb_run_id == 8
        # Same tile-based path as run 7, over a wider ΔT axis: `[-3, -1, -0.5, 0, 0.5, 1, 3, 4, 5, 6]`
        # against run 7's `[-3, -1, -0.5, 0, 0.5, 1, 3]`. `merge_tile_perturbations` built each tile here
        # by joining a +4/+5/+6 supplementary sweep onto the corresponding run-7 tile, so the two trees
        # share every point run 7 has and this one only adds coverage above +3 K, where `gemb_calibration`
        # pins a fitted ΔT at the grid ceiling for a sizeable share of global glacier area under run 7.
        gemb_info = (;
            tile_dir=joinpath(get(ENV, "CLIMATE_CACHE", "/mnt/bylot-r3/data/era5land"), "tile_runs_merged"),
            precipitation_scale=nothing,
            filename_gemb_combined="/mnt/bylot-r3/data/gemb/raw/gemb_glacier_sims_tiles_1950to2026_merged.jld2",
            modify_melt_only=false
        )
    else
        error("unrecognized gemb_run_id: $gemb_run_id")
    end

    return gemb_info
end


"""
    geotile_groups_forced()

Return manual overrides for geotile groupings where large glaciers span multiple tiles.

Used when glaciers that cross tile boundaries should be grouped or treated separately.
Groupings are only applied when geotile grouping is recomputed (e.g. force_remake_before is set).

# Returns
- Vector of vectors of geotile ID strings; each inner vector is one forced group.

# Examples
```julia
julia> forced = geotile_groups_forced()
julia> geotile_grouping!(geotiles0, glaciers, 100; geotile_groups_manual=forced)
```
"""
function geotile_groups_forced() 
    out = [
        ["lat[+62+64]lon[-148-146]"],
        ["lat[+62+64]lon[-146-144]"],
        ["lat[+60+62]lon[-138-136]"],
        ["lat[+58+60]lon[-138-136]"],
        ["lat[+58+60]lon[-136-134]"],
        ["lat[+58+60]lon[-134-132]"],
        ["lat[+56+58]lon[-134-132]"],
        ["lat[+78+80]lon[-084-082]"],
        ["lat[-52-50]lon[-074-072]"],
        ["lat[+56+58]lon[-130-128]"],
        ["lat[-74-72]lon[-080-078]", "lat[-74-72]lon[-078-076]"],
        ["lat[+78+80]lon[+010+012]", "lat[+78+80]lon[+012+014]", "lat[+78+80]lon[+014+016]"],
        ["lat[+76+78]lon[+062+064]", "lat[+76+78]lon[+064+066]", "lat[+76+78]lon[+066+068]", "lat[+76+78]lon[+068+070]", "lat[+74+76]lon[+062+064]", "lat[+74+76]lon[+064+066]", "lat[+74+76]lon[+066+068]", "lat[+74+76]lon[+066+068]"]
    ]
    return out
end

const plot_order = Dict("missions" => ["hugonnet", "icesat", "gedi", "icesat2"], "synthesis" => ["hugonnet", "ICESat & ICESat 2", "gedi", "Synthesis"])

"""
    gemb_altim_cost(x, dv_altim, dv_gemb, kwargs)

Compute cost between altimetry-derived volume change and GEMB-sampled volume change.

Samples GEMB at (pscale, ΔT) from x, subtracts from dv_altim, then evaluates
model_fit_cost_function on the residuals with the given kwargs.

# Arguments
- `x`: Parameter vector [pscale, ΔT]
- `dv_altim`: DimArray of altimetry-derived volume change (geotile, date)
- `dv_gemb`: GEMB volume change array used for sampling
- `kwargs`: Keyword arguments passed to model_fit_cost_function (e.g. seasonality_weight, distance_from_origin_penalty)

# Returns
- Scalar cost value.

# Examples
```julia
julia> cost = gemb_altim_cost([1.0, 0.0], dv_altim, dv_gemb, (; seasonality_weight=0.85, distance_from_origin_penalty=0.2, ΔT_to_pscale_weight=0.5))
```
"""
function gemb_altim_cost(x, dv_altim, dv_gemb, kwargs)
    pscale = x[1]
    ΔT = x[2]

    res = dv_altim .- gemb_dv_sample(pscale, ΔT, dv_gemb)
    res .-= mean(res)

    cost = model_fit_cost_function(res, pscale, ΔT; kwargs...)

    return cost
end

"""
    model_fit_cost_function(res, pscale, ΔT; seasonality_weight, distance_from_origin_penalty, ΔT_to_pscale_weight, calibrate_to=:all, is_scaling_factor, origin_penalty_mode=:prior, forcing_prior=gemb_forcing_prior)

Compute a composite cost function for fitting a model to altimetry data, incorporating trend, seasonality, and parameter penalties.

# Arguments
- `res`: Residuals between observed and modeled values. Should be a DimArray or array-like object with a :date dimension.
- `pscale`: Precipitation scaling factor (numeric).
- `ΔT`: Air temperature offset applied to the GEMB forcing (numeric, in K). An additive parameter whose
  no-op is 0, unlike `pscale`; `is_scaling_factor` selects which distance the `:legacy` penalty uses.
- `seasonality_weight`: Weight (0–1) for the seasonal amplitude in the cost function. Higher values emphasize seasonality.
- `distance_from_origin_penalty`: Strength `wd` of the penalty on the forcing corrections.
- `ΔT_to_pscale_weight`: Weight of the temperature-offset penalty against the precipitation-scaling one
  (`:legacy` only).
- `calibrate_to`` = [:all, :annual, :five_year, :trend]
- `origin_penalty_mode`: `:prior` (default) multiplies the fit cost by `1 + wd * d`, `d` the Mahalanobis
  distance of (log pscale, ΔT) from `forcing_prior`, i.e. the distance in prior standard deviations
  allowing for the correlation between the two. `:legacy` multiplies it by `1 + wd * d`, `d` the weighted
  distance from the no-op point (1, 0).
- `forcing_prior`: `(; pscale, log_pscale_sd, ΔT, ΔT_sd, corr)`, the centre, spreads and correlation used
  by `origin_penalty_mode = :prior`; defaults to `gemb_forcing_prior`.

# Returns
A tuple `cost` where:
- `cost`: Composite cost metric, combining RMSE, seasonal amplitude, and penalties for parameter deviation.

# Details
- The function fits a seasonal model to the residuals using `ts_seasonal_model`.
- If `calibrate_to_annual_change_only` is `true`, the cost is computed using only the annual change (residuals at the seasonal minimum).
- If `calibrate_to_trend_only` is `true`, the cost is based only on the absolute value of the linear trend.
- Otherwise, the cost is a weighted sum of RMSE and the amplitude of the seasonal cycle.
- The penalty measures deviation of (`pscale`, `ΔT`) from `forcing_prior`, or from the no-op point (1, 0)
  when `origin_penalty_mode = :legacy`.

# Examples
```julia
julia> cost = model_fit_cost_function(res, 1.0, 0.0; seasonality_weight=0.85, distance_from_origin_penalty=0.35, ΔT_to_pscale_weight=0.5)
```
"""
function model_fit_cost_function(res, pscale, ΔT; seasonality_weight, distance_from_origin_penalty, ΔT_to_pscale_weight, calibrate_to=:all, is_scaling_factor=Dict("pscale" => true, "ΔT" => false), origin_penalty_mode=:prior, forcing_prior=gemb_forcing_prior)

    origin_penalty_mode in (:legacy, :prior) || error("origin_penalty_mode must be :legacy or :prior, got $origin_penalty_mode")

    # remove linear trend to emphasize seasonality
    if calibrate_to != :five_year
        fit = ts_seasonal_model(res; interval=nothing);
    else
        res0 = groupby(res, :date => Bins(year, 5))
        res0 = mean.(res0)
        rmse_cost = sqrt(mean(res0 .^ 2))
    end

    # calibrate to annual change only
    if calibrate_to == :annual
        seasonal_min =  mod1(fit.phase_peak_month + 6, 12)
        res = res[date=Near(DateTime(minimum(year.(dims(res, :date))), seasonal_min, 15):Year(1):DateTime(maximum(year.(dims(res, :date))), seasonal_min, 15))]
        rmse_cost = sqrt(mean(res .^ 2))
    end

    if calibrate_to == :all
        rmse = sqrt(mean(res .^ 2))
        cost = ((1 - seasonality_weight) * rmse) + (seasonality_weight * fit.amplitude)
    elseif calibrate_to == :trend
        cost = (1 - seasonality_weight) * abs(fit.trend)
    else
        cost = (1 - seasonality_weight) * rmse_cost
    end

    if origin_penalty_mode == :prior
        (is_scaling_factor["pscale"] && !is_scaling_factor["ΔT"]) ||
            error("origin_penalty_mode = :prior assumes a multiplicative pscale and an additive ΔT")
        return cost * (1 + distance_from_origin_penalty * _forcing_prior_distance(pscale, ΔT, forcing_prior))
    else
        # Distance of a scaling factor from its no-op value of 1, measured symmetrically so that
        # halving and doubling are penalized equally.
        dp = _distance_from_origin(pscale, is_scaling_factor["pscale"])
        dT = _distance_from_origin(ΔT, is_scaling_factor["ΔT"])
        return cost * (1 + sqrt((dT * ΔT_to_pscale_weight)^2 + (dp * (1 - ΔT_to_pscale_weight))^2) * distance_from_origin_penalty)
    end
end

# Mahalanobis distance of (log pscale, ΔT) from the prior centre, in prior standard deviations.
function _forcing_prior_distance(pscale, ΔT, prior)
    z1 = (log(pscale) - log(prior.pscale)) / prior.log_pscale_sd
    z2 = (ΔT - prior.ΔT) / prior.ΔT_sd
    ρ = prior.corr
    return sqrt(max((z1^2 - 2ρ * z1 * z2 + z2^2) / (1 - ρ^2), 0.0))
end

function _distance_from_origin(scale, is_scaling_factor)
    is_scaling_factor || return scale
    return scale < 1 ? 1 / scale - 1 : scale - 1
end
