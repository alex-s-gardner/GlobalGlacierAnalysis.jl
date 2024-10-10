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

#=
# regional_dischage (Gt/yr)
discharge_gtyr = Dict(
    "rgi1" => 17.1, # Burgess, 2013
    "rgi2" => 0.75, # Menounos personal communication
    "rgi3" => 2.2, # van Wychen 2016 
    "rgi4" => 0.06, # van Wychen 2015
    "rgi5" => 5.1, # Bollen et al, 2021 AGU poster
    "rgi6" => 2.1, # Tómas Jóhannesson, 2020 (geothermal and lake calving)
    "rgi7" => 10, # moholdt personal comunication [6.75; % Baszczyk, 2009]
    "rgi8" => 0,
    "rgi9" => 10, # moholdt personal comunication  (0.5+2.26+0.65) * dh2dm; % 0.5 Novaya Zemlya, Melkonian 2016; 
    "rgi10" => 0,
    "rgi11" => 0,
    "rgi12" => -4, # this is done as the model generates WAY too much melt that a match can't be found... # -4 Gt/yr is the mean rate for the icesat period with no smb scaling
    "rgi13" => 0,
    "rgi14" => 0,
    "rgi15" => 0,
    "rgi16" => 0,
    "rgi17" => 61.3, # Schaefer, 2015
    "rgi18" => 0,
    "rgi19" => 70, #75 Gt/yr is the mean rate for the icesat period -6 gt with no smb scaling.  If set to NaN then deteremined inside of gem_Δvolume!(df)
)
=#



function mission_land_trend()
    missions = ["icesat", "icesat2", "gedi", "hugonnet"]
    dmission = Dim{:mission}(missions)
    mission_trend_myr = fill(0.0, dmission)
    mission_trend_myr[At(["gedi"])] .= -0.144

    return mission_trend_myr
end
