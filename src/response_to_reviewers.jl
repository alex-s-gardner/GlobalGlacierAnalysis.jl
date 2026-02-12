# Load packages and data paths
begin
    using CSV
    using DataFrames
    import GlobalGlacierAnalysis as GGA
    using DimensionalData
    using Dates
    using CairoMakie
    using Statistics
    using ProgressMeter
    using NCDatasets
    using FileIO
    using CFTime
    using GeoDataFrames
    
    # Import functions from utilities
    include("utilities_response2reviewers.jl")

    # data paths
    paths = (
        wgms_mb = "/mnt/bylot-r3/data/glacier_mb/WGMS/DOI-WGMS-FoG-2025-02b/data/mass_balance.csv",
        wgms_glacier = "/mnt/bylot-r3/data/glacier_mb/WGMS/DOI-WGMS-FoG-2025-02b/data/glacier.csv",

        # GEMB dataset paths
        gemb_corrected = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_NH_1979to2024_820_40_racmo_grid_lwt_e97_0_corrected_geotile_dv.jld2",
        gemb_original = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2024_820_40_racmo_grid_lwt_e97_0_geotile_dv.jld2",

        path2glacier = GGA.pathlocal[Symbol("glacier_individual")],
    )

    # Analysis parameters
    study_start_year = 2000
    min_neighbors_variogram = 3
    geotile_thresholds = [1, 20, 65]
    histogram_bounds = (-5, 3)
    historic_analysis_mincount_range = 5:15

    # RGI glacier counts for context
    rgi6_count = 215547
    rgi7_count = 274531
end




# load hugonnet data
# load in all regions
begin
    hugonnet = DataFrame()
    for rgi = 1:19
        rgi_str = lpad(string(rgi), 2, '0')
        path2hugonnet = "/mnt/bylot-r3/data/glacier_mb/Hugonnet2021/time_series/time_series_$(rgi_str)/dh_$(rgi_str)_rgi60_pergla_rates.csv"
        hugonnet = vcat(hugonnet, CSV.read(path2hugonnet, DataFrame))
    end

    # find annual balance and format similar to wgms_mb for spatial variogram
    hugonnet[!, "start_date"] = DateTime.(getindex.(split.(hugonnet.period, "_"), 1))
    hugonnet[!, "end_date"] = DateTime.(getindex.(split.(hugonnet.period, "_"), 2))

    hugonnet[!, "end_year"] = year.(hugonnet[!, "end_date"])
    hugonnet[!, "start_year"] = year.(hugonnet[!, "start_date"])
    hugonnet[!, "dt_years"] = hugonnet[!, "end_year"] - hugonnet[!, "start_year"]
    hugonnet[!, "year"] = hugonnet[!, "end_year"]
    index = hugonnet[!, "dt_years"] .== 1
    hugonnet = hugonnet[index, :]
    hugonnet[!, :annual_balance] = hugonnet.dhdt * 0.85


    rgi6 = GeoDataFrames.read(GGA.pathlocal.glacier_individual)
    index = indexin(hugonnet.rgiid, rgi6.RGIId)

    # it seems some ids do not match... I think this is because hugonnet updated glaciers that are not in rgi6
    index_not_match = isnothing.(index)
    hugonnet = hugonnet[.!index_not_match, :]

    hugonnet[!, :longitude] = rgi6.CenLon[index[.!index_not_match]]
    hugonnet[!, :latitude] = rgi6.CenLat[index[.!index_not_match]]

subset_random = rand(1:nrow(hugonnet), round(Int, nrow(hugonnet)/100))
variogram, ddistance = calculate_spatial_variogram(deepcopy(hugonnet[subset_random, [:annual_balance, :longitude, :latitude, :year]]); mincount=min_neighbors_variogram);

f_var_hugonnet = plot_variogram(ddistance, variogram, 0, 0)
display(f_var_hugonnet)
save(joinpath(GGA.pathlocal.figures, "wgms_mb_variogram_hugonnet.png"), f_var_hugonnet)

end


# link hugonnet glaciers to wgms_mb
index_into_wgms = indexin(hugonnet.rgiid, wgms_glacier.rgi60_ids)
index_valid = .!isnothing.(index_into_wgms)
hugonnet[!, :wgms_id] .= 0
hugonnet[index_valid, :wgms_id] = wgms_glacier.id[index_into_wgms[index_valid]]


wgms_mb[!, :annual_balance_hugonnet] .= NaN
for r in eachrow(wgms_mb)
    index = findfirst((hugonnet.dt_years .== 1) .& (hugonnet.wgms_id .== r.glacier_id) .& (hugonnet.year .== r.year))
    if !isnothing(index)
        r.annual_balance_hugonnet = hugonnet.dhdt[index] .* 0.85
    end
end

delta = wgms_mb.annual_balance_hugonnet .- wgms_mb.annual_balance
index_valid = (.!ismissing.(delta)) .& (.!isnan.(delta))
f = Figure();
ax = Axis(f[1, 1]; xlabel="wgms minus hugonnet [m w.e. yr⁻¹]", ylabel="count");
hist!(ax, delta[index_valid])
f
save(joinpath(GGA.pathlocal.figures, "wgms_mb_histogram_hugonnet.png"), f)


# Create precipitation raster of maximum monthly precipitation
begin
    path2precip = "/mnt/bylot-r3/data/glacier_mb/GPCP/precip.mon.ltm.1991-2020.nc"

    precip = GGA.Rasters.Raster(path2precip, name=:precip);
    precip_max = dropdims(getindex.(argmax(precip, dims=3), 3), dims=:Ti)

    output_path = joinpath(splitpath(path2precip)[1:end-1]..., "precip.mon_of_max.ltm.1991-2020.tif")
    precip_max = reverse(precip_max, dims=:Y)

    precip_max = rebuild(precip_max, data = UInt8.(precip_max.data))
    precip_max = set(precip_max, X => dims(precip_max, X).val .-180)

    GGA.Rasters.write(output_path, precip_max, force=true)
end;



