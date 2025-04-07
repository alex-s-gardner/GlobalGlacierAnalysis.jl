#density of glacier ice
const δice = 910; #kg m-3

""" 
    ElevationProduct
Defines ElevationProduct parameter type
"""
@kwdef struct ElevationProduct
    mission::Symbol
    name::Symbol
    version::Int64
    id::String       
    error_sigma::Float32
    halfwidth::Float32
    kernel::Symbol
    apply_quality_filter::Bool
    coregister::Bool
end


# project specific helper funcitons
function project_products(; project_id = :v01)
    if project_id == :v01
        product = (
            icesat2=ElevationProduct(mission=:icesat2, name=:ATL06, version=6, id="I206", error_sigma=0.1, halfwidth=11 / 4, kernel=:gaussian, apply_quality_filter=false, coregister=true),

            icesat=ElevationProduct(mission=:icesat, name=:GLAH06, version=34, id="I106", error_sigma=0.1, halfwidth=35 / 4, kernel=:gaussian, apply_quality_filter=false, coregister=true),

            gedi=ElevationProduct(mission=:gedi, name=:GEDI02_A, version=2, id="G02A", error_sigma=0.1, halfwidth=22 / 4, kernel=:gaussian, apply_quality_filter=false, coregister=true),

            hugonnet=ElevationProduct(mission=:hugonnet, name=:HSTACK, version=1, id="HS01", error_sigma=5, halfwidth=100 / 2, kernel=:gaussian, apply_quality_filter=false, coregister=true)
        )
    end
    return product
end


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

function project_geotiles(; geotile_width = 2, domain = :all, extent=nothing)
    paths = setpaths()
    geotiles = GeoTiles.define(geotile_width; extent)

    if domain == :all

    elseif domain == :landice
        icemaskfn = paths.icemask
        icemask = GeoArrays.read(icemaskfn);
        geotilemask = GeoArrays.crop.(Ref(icemask), Altim.extent2nt.(geotiles.extent))
        #geotilemask = GeoArrays.crop.(Ref(icemask), geotiles.extent)
        hasice = [any(mask.A.==1) for mask in geotilemask];
        geotiles = geotiles[hasice,:];
    else
        error("unrecognized domain = $domain")
    end
    return geotiles
end

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


function project_date_bins()
        Δd = 30
        date_range = DateTime(1990):Day(Δd):DateTime(2026, 1, 1)
        date_center = date_range[1:end-1] .+ Day(Δd / 2)

    return date_range, date_center
end

function project_height_bins()
    Δh = 100;
    height_range = 0:100:10000;
    height_center = height_range[1:end-1] .+ Δh / 2;

    return height_range, height_center
end


function mission_land_trend()
    missions = ["icesat", "icesat2", "gedi", "hugonnet"]
    dmission = Dim{:mission}(missions)
    mission_trend_myr = fill(0.0, dmission)
    mission_trend_myr[At(["gedi"])] .= -0.144

    return mission_trend_myr
end
