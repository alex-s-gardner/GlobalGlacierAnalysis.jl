# This script is to explore the seasonality of the glaciers for Chad's paper
begin
    using DataFrames, CSV, Plots, Dates, CairoMakie
    using NearestNeighbors
    import GeometryOps as GO
    using FlexiJoins
    using Altim
    using FileIO
    using DimensionalData

    # seasonality data from Chad
    master_file = "/mnt/devon-r2/data/seasonality/data/flowline_seasonality_compiler_2025-02-06.csv"

    # glacier smb data
    paths = Altim.pathlocal
    reference_run = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_dh_best_cc_meanmadnorm3_v01_filled_ac_p1_synthesized.jld2";
    perglacier_reference_run = replace(reference_run, ".jld2" => "_perglacier.jld2")
    globaldischarge_fn = joinpath(paths.data_dir, "GlacierOutlines/GlacierDischarge/global_glacier_discharge.jld2")
end

# find nearest point in glacier_flux to seasonality data
# It takes 10 min to load and manipulate the data
@time begin
    seasonality = CSV.read(master_file, DataFrame)
    glacier_flux = load(perglacier_reference_run, "glaciers") 
    discharge = load(globaldischarge_fn, "discharge")

    # only consider main flowlines
    seasonality = seasonality[seasonality.is_main.==1, :]


    # add geometry column to seasonality and discharge
    seasonality[!, :geometry] = tuple.(seasonality.lon, seasonality.lat)
    discharge[!, :geometry] = tuple.(discharge.longitude, discharge.latitude)
    rename!(glacier_flux, :geom => :geometry)
    # fit trends to glacier flux data
    # find the index of valid data
    varnames = [:runoff, :smb, :dh]

    for varname in varnames
        date = collect(dims(glacier_flux[1, varname], :date))
        index = .!isnan.(glacier_flux[1, varname])
        datelimits = extrema(date[index])
        DataFrames.colmetadata!(glacier_flux, varname, "date", date; style=:note)
        @time Altim.df_tsfit!(glacier_flux, [varname]; progress=true, datelimits)
    end

    # remove flux time series
    glacier_flux = glacier_flux[:, Not(:runoff, :smb, :dh, :fac)]
    # this takes 5 min
    glacier_flux[!, :geometry_buffer] = GO.buffer.(glacier_flux.geometry, 0.00001)
end;

begin
    # join glacier flux and discharge
    glacier_flux = FlexiJoins.rightjoin(
        (discharge, glacier_flux),
        by_pred(:geometry, GO.within, :geometry_buffer)
    )

    # join seasonality and glacier flux
    seasonality = FlexiJoins.leftjoin(
        (seasonality, glacier_flux),
        by_pred(:geometry, GO.within, :geometry_buffer)
    )


    vars2remove = [
        "geometry"
        "latitude"
        "longitude"
        "geometry_2"
        "geometry_1"
        "RGIId"
        "GLIMSId"
        "BgnDate"
        "EndDate"
        "CenLon"
        "CenLat"
        "O1Region"
        "O2Region"
        "Area"
        "Name"
        "layer"
        "path"
        "geotile"
        "area_km2_1"
        "runoff_offset"
        "smb_offset"
        "dh_offset"
        "geometry_buffer"
        ]

    seasonality = seasonality[:, Not(vars2remove)]

    # convert from m to Gtyr
    vas2convert = ["runoff_trend"
        "runoff_amplitude"
        "smb_trend"
        "smb_amplitude"
        "dh_trend"
        "dh_amplitude"
    ]

    for varname in vas2convert
        seasonality[:, varname] = seasonality[:, varname] .* seasonality.area_km2 / 1000
        rename!(seasonality, varname => Symbol(varname, "_gtyr"))
    end

    CSV.write(replace(master_file, ".csv" => "_plus.csv"), seasonality)
end

var0 = "dh_phase"
valid = .!(isnan.(seasonality[:,var0]) .| ismissing.(seasonality[:,var0]))
CairoMakie.hist(seasonality[valid,var0])
