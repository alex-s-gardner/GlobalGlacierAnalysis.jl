#density of glacier ice
const δice = 910; #kg m-3
const local2utc = Hour(7) # LA timezone to UTC
const seasonality_weight = 85/100
const distance_from_origin_penalty = 5/100
const ΔT_to_pscale_weight = 50/100




"""
    project_products(; project_id = :v01)

Get the elevation products configuration for a specific project.

# Arguments
- `project_id::Symbol`: Project identifier (default: :v01)

# Returns
- Named tuple containing configured ElevationProduct instances for different missions
- NOTE: latitude_limits are the latitude limits of the mission's data product and need to be specified as Integers
"""
function project_products(; project_id = :v01)
    if project_id == :v01
        product = (
            icesat2=(mission=:icesat2, name=:ATL06, version=6, id="I206", error_sigma=0.1, halfwidth=11 / 4, kernel=:gaussian, apply_quality_filter=false, coregister=true, latitude_limits=[-88, 88], longitude_limits=[-180, 180]),

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
"""
function project_paths(; project_id = :v01)
    if project_id == :v01
        geotile_width = 2 #geotile width [degrees]

        p = project_products(project_id = project_id)

        paths = (
            icesat2 = setpaths(geotile_width, :icesat2, "$(p.icesat2.name)", lpad("$(p.icesat2.version)", 3, '0')),
            icesat = setpaths(geotile_width, :icesat, "$(p.icesat.name)", lpad("$(p.icesat.version)", 3, '0')),
            gedi = setpaths(geotile_width, :gedi, "$(p.gedi.name)", lpad("$(p.gedi.version)", 3, '0')),
            hugonnet = setpaths(geotile_width, :hugonnet, "$(p.hugonnet.name)", lpad("$(p.hugonnet.version)", 3, '0'))
        );
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

# Returns
- Tuple containing (date_range, date_center) where:
  - date_range: DateTime range with 30-day intervals
  - date_center: DateTime values at the center of each bin
"""
function project_date_bins()
        Δd = 30
        date_range = Date(1990):Day(Δd):Date(2026, 1, 1)
        date_center = date_range[1:end-1] .+ Day(Δd / 2)

    return date_range, date_center
end

"""
    project_height_bins()

Define elevation bins for the project.

# Returns
- Tuple containing (height_range, height_center) where:
  - height_range: Range of elevation bin edges from 0 to 10000m at 100m intervals
  - height_center: Values at the center of each elevation bin
"""
function project_height_bins()
    Δh = 100;
    height_range = 0:100:10000;
    height_center = height_range[1:end-1] .+ Δh / 2;

    return height_range, height_center
end


function project_ΔT_bins()
    ΔT = 1
    ΔT_range = -10.5:ΔT:10.5
    ΔT_center = ΔT_range[1:end-1] .+ ΔT / 2

    return ΔT_range, ΔT_center
end

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
"""
function mission_land_trend()
    missions0 = ["icesat", "icesat2", "gedi", "hugonnet"]
    dmission = Dim{:mission}(missions0)
    mission_trend_myr = fill(0.0, dmission)
    mission_trend_myr[At(["gedi"])] .= -0.144

    return mission_trend_myr
end

geotiles_golden_test = [
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
"""
function gemb_info(; gemb_run_id = 4)

    if gemb_run_id == 1
        dpscale = Dim{:pscale}(["p1"]) # do not change order as these are lookup values
        dΔheight = Dim{:Δheight}(["t1"])
        elevation_delta = DimArray([0], dΔheight) # do not change order as these are lookup values
        precipitation_scale = DimArray([1], dpscale) # do not change order as these are lookup values
        gemb_info = (;
            gemb_folder = ["/home/schlegel/Share/GEMBv1/"],
            file_uniqueid = "rv1_0_19500101_20231231",
            elevation_delta,
            precipitation_scale,
            filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/$file_uniqueid.jld2",
            modify_melt_only = false
        )
    elseif gemb_run_id == 2
        dpscale = Dim{:pscale}(["p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9"]) # do not change order as these are lookup values
        dΔheight = Dim{:Δheight}(["t1", "t2", "t3", "t4", "t5", "t6"])
        elevation_delta = DimArray([-1000, -750, -500, -250, 0, 250, 500, 750, 1000], dΔheight) # do not change order as these are lookup values
        precipitation_scale = DimArray([0.5, 1, 1.5, 2, 5, 10], dpscale) # do not change order as these are lookup values
        gemb_info = (;
            gemb_folder = "/home/schlegel/Share/GEMBv1/Alaska_sample/v1/",
            file_uniqueid = "1979to2023_820_40_racmo_grid_lwt",
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
    else
        error("unrecognized gemb_run_id: $gemb_run_id")
    end

    return gemb_info
end


# Manual override for specific geotile groups
# This handles cases where large glaciers cross multiple tiles but should be treated separately
# NOTE: groupings are only updated if force_remake == true
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

plot_order = Dict("missions" => ["hugonnet", "icesat", "gedi", "icesat2"], "synthesis" => ["hugonnet", "ICESat & ICESat 2", "gedi", "Synthesis"])



function gemb_altim_cost(x, dv_altim, dv_gemb, kwargs)
    pscale = x[1]
    ΔT = x[2]
    res = dv_altim .- gemb_dv_sample(pscale, ΔT, dv_gemb)
    
    cost = model_fit_cost_function(res, pscale, ΔT; kwargs...)

    return cost
end

function gemb_altim_cost!(res::DimArray, gemb_dv_out, x, dv_altim, dv_gemb, kwargs)
    res[:] = dv_altim .- gemb_dv_sample!(gemb_dv_out, x[1], x[2], dv_gemb)
    res .-= mean(res)


    cost = model_fit_cost_function(res, x[1], x[2]; kwargs...)

    return cost
end


"""
    model_fit_cost_function(res, pscale, ΔT; seasonality_weight, distance_from_origin_penalty, calibrate_to_trend_only=false, calibrate_to_annual_change_only=true)

Compute a composite cost function for fitting a model to altimetry data, incorporating trend, seasonality, and parameter penalties.

# Arguments
- `res`: Residuals between observed and modeled values. Should be a DimArray or array-like object with a :date dimension.
- `pscale`: Precipitation scaling factor (numeric).
- `ΔT`: Elevation offset (numeric, in meters).
- `seasonality_weight`: Weight (0–1) for the seasonal amplitude in the cost function. Higher values emphasize seasonality.
- `distance_from_origin_penalty`: Penalty factor for deviation of parameters from their reference values.
- `ΔT_to_pscale_weight`: Weight for the difference between melt scaling and precipitation scaling in the cost function.
- `calibrate_to`` = [:all, :annual, :five_year, :trend]

# Returns
A tuple `cost` where:
- `cost`: Composite cost metric, combining RMSE, seasonal amplitude, and penalties for parameter deviation.

# Details
- The function fits a seasonal model to the residuals using `ts_seasonal_model`.
- If `calibrate_to_annual_change_only` is `true`, the cost is computed using only the annual change (residuals at the seasonal minimum).
- If `calibrate_to_trend_only` is `true`, the cost is based only on the absolute value of the linear trend.
- Otherwise, the cost is a weighted sum of RMSE and the amplitude of the seasonal cycle.
- Penalty terms are applied for deviation of `pscale` from 1 and `ΔT` from 0 (scaled to kilometers).
"""
function model_fit_cost_function(res, pscale, ΔT; seasonality_weight, distance_from_origin_penalty, ΔT_to_pscale_weight, calibrate_to = :all, is_scaling_factor=Dict("pscale" => true, "ΔT" => true))
    
    # remove linear trend to emphasize seasonality
    if calibrate_to != :five_year
        fit = ts_seasonal_model(res; interval=nothing);
    else
        res0 = groupby(res, :date => Bins(year, 5))
        res0 = ts_seasonal_model.(res0)
        interval_trends = map(p -> p.trend, res0)
        rmse_cost =sqrt(mean(interval_trends .^ 2))
    end
    
    # calibrate to annual change only
    if calibrate_to == :annual
        seasonal_min =  mod1(fit.phase_peak_month + 6, 12)
        res = res[date=Near(DateTime(minimum(year.(dims(res, :date))), seasonal_min, 15):Year(1):DateTime(maximum(year.(dims(res, :date))), seasonal_min, 15))]
        rmse_cost = sqrt(mean(res .^ 2))
    end

    if is_scaling_factor["pscale"]

        if pscale < 1
            dp = 1/pscale - 1
        else
            dp = pscale - 1
        end
    else
        dp = pscale;
    end

    if is_scaling_factor["ΔT"]
         if ΔT < 1
            dT = 1/ΔT - 1
        else
            dT = ΔT - 1
        end
    else
        dT = ΔT;
    end

    if calibrate_to == :all
        rmse = sqrt(mean(res .^ 2)) 
        cost = ((1 - seasonality_weight) * rmse + seasonality_weight * fit.amplitude) * (1 + (sqrt((ΔT * (ΔT_to_pscale_weight))^2 + (dp * (1 - ΔT_to_pscale_weight))^2) * distance_from_origin_penalty))
    elseif calibrate_to == :trend
        cost = ((1 - seasonality_weight) * abs(fit.trend)) * (1 + (sqrt((ΔT * (ΔT_to_pscale_weight))^2 + (dp * (1 - ΔT_to_pscale_weight))^2) * distance_from_origin_penalty)
    else
        cost = ((1 - seasonality_weight) * rmse_cost) * (1 + (sqrt((ΔT * (ΔT_to_pscale_weight))^2 + (dp * (1 - ΔT_to_pscale_weight))^2) * distance_from_origin_penalty)
    end

    return cost
end


function model_fit_cost_function2(res, pscale_idx, ΔT_idx; seasonality_weight, distance_from_origin_penalty, ΔT_to_pscale_weight, pscale_idx2dist=nothing, ΔT_idx2dist=nothing)

    # remove linear trend to emphasize seasonality

    fit = ts_seasonal_model(res; interval=nothing)
    rmse = sqrt(mean(res .^ 2))
    cost = ((1 - seasonality_weight) * rmse + seasonality_weight * fit.amplitude) * (1 + (sqrt((ΔT_idx2dist(ΔT_idx) * (ΔT_to_pscale_weight))^2 + (pscale_idx2dist(pscale_idx) * (1 - ΔT_to_pscale_weight))^2) * distance_from_origin_penalty))

    return cost
end
