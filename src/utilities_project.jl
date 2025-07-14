#density of glacier ice
const δice = 910; #kg m-3
const local2utc = Hour(7) # LA timezone to UTC

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
        date_range = DateTime(1990):Day(Δd):DateTime(2026, 1, 1)
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
    "lat[+64+66]lon[+058+060]"
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
            filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/$file_uniqueid.jld2"
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
            filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/$file_uniqueid.jld2"
        )
    elseif gemb_run_id == 3
        dpscale = Dim{:pscale}(["p1", "p2", "p3", "p4", "p5"]) # do not change order as these are lookup values
        dΔheight = Dim{:Δheight}(["t1", "t2", "t3", "t4", "t5"])
        elevation_delta = DimArray([-200, 0, 200, 500, 1000], dΔheight) # do not change order as these are lookup values
        precipitation_scale = DimArray([0.75, 1, 1.25, 1.5, 2], dpscale) # do not change order as these are lookup values

        gemb_info = (;
            gemb_folder = ["/home/schlegel/Share/GEMBv1/NH_sample/", "/home/schlegel/Share/GEMBv1/SH_sample/"],
            file_uniqueid="1979to2023_820_40_",
            elevation_delta,
            precipitation_scale,
            filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0.jld2"
        )
    elseif gemb_run_id == 4
        dpscale = Dim{:pscale}(["p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8"]) # do not change order as these are lookup values
        dΔheight = Dim{:Δheight}(["t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8"])
        elevation_delta = DimArray([-200, 0, 200, 500, 1000, -2000, -500, 2000], dΔheight) # do not change order as these are lookup values
        precipitation_scale = DimArray([0.75, 1, 1.25, 1.5, 2, 0.25, 3, 4], dpscale) # do not change order as these are lookup values

        gemb_info = (;
            gemb_folder = ["/home/schlegel/Share/GEMBv1/NH_sample", "/home/schlegel/Share/GEMBv1/SH_sample"],
            file_uniqueid="1979to2024_820_40_",
            elevation_delta,
            precipitation_scale,
            filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2024_820_40_racmo_grid_lwt_e97_0.jld2"
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



const seasonality_weight = 85/100
const distance_from_origin_penalty = 5/100

"""
    model_fit_cost_function(res, pscale, Δheight)

Compute the cost function for fitting the model to altimetry data.

# Arguments
- `res`: Residuals between observed and modeled values (DimArray or array-like with :date dimension)
- `pscale`: Precipitation scaling factor (numeric)
- `Δheight`: Elevation offset (numeric)

# Returns
A tuple `(rmse, cost)` where:
- `rmse`: Root mean square error of the residuals
- `cost`: Composite cost metric including RMSE, seasonality, and parameter penalties

# Details
- Removes the linear trend from the residuals to emphasize seasonality.
- The cost metric is the sum of the RMSE, the RMSE of the detrended residuals, and penalty terms for deviation of `pscale` from 1 and `Δheight` from 0 (scaled by 6.5/1000).
"""
function model_fit_cost_function(res, pscale, Δheight; seasonality_weight, distance_from_origin_penalty) 

    decyear = decimalyear.(val(dims(res, :date)))
    decyear .-= mean(decyear)

    # remove linear trend to emphasize seasonality
    fit = ts_seasonal_model(res; interval=nothing);

    rmse = sqrt(mean(res .^ 2))

    if pscale < 1
        dp = 1/pscale - 1
    else
        dp = pscale - 1
    end

    dh = abs(Δheight / 1000)

    cost = ((1 - seasonality_weight) * rmse + seasonality_weight * fit.amplitude) * (1 + sqrt(dh^2 + dp^2) * distance_from_origin_penalty)

    return (rmse, cost)
end

