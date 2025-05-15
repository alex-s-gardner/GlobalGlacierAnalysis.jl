"""
    gmax_global.jl

Calculate the climatological maximum monthly fractional contribution of glacier runoff to total river flux (gmax).

This script:
1. Loads glacier runoff and river flux data from NetCDF files
2. Calculates the fraction of total river flux contributed by glacier runoff
3. Determines the maximum monthly contribution (gmax) for each river segment
4. Processes river network data to identify terminal points (sinks)
5. Adds glacier discharge data from external sources
6. Writes results to GeoPackage files for further analysis

The script handles time series data from 2000-2024 and focuses on peak runoff periods.
"""

begin
    using NCDatasets
    using GeoDataFrames
    using DimensionalData
    using GlobalGlacierAnalysis
    using Dates
    using Statistics
    import GeometryOps as GO
    import GeoInterface as GI
    using FileIO
    using DataFrames
    using CairoMakie

    dates4trend = [DateTime(2000, 3, 1), DateTime(2024, 12, 15)]

    paths = GlobalGlacierAnalysis.pathlocal

    glacier_summary_file = joinpath(paths[:project_dir], "gardner2025_glacier_summary.nc")
    glacier_summary_riverflux_file = replace(glacier_summary_file, ".nc" => "_riverflux.nc")
    glacier_summary_gmax_file = replace(glacier_summary_file, ".nc" => "_gmax.pkg")

    glacier_rivers_land_flux_path = joinpath(paths[:river], "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier_Qs_acc_Qsb_acc.nc")
    glacier_rivers_snow_flux_path = joinpath(paths[:river], "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier_Qsm_acc.nc")

    path2discharge = joinpath(paths[:data_dir], "GlacierOutlines/GlacierDischarge/global_glacier_discharge.jld2")
    rivers_paths = GlobalGlacierAnalysis.allfiles(paths[:river]; fn_endswith="MERIT_Hydro_v07_Basins_v01.shp", fn_startswith="riv_pfaf")
    glacier_rivers_path = joinpath(paths[:river], "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.arrow")

    dates4trend = (Date(2000,1,1), Date(2024,10,1))
    SECONDS_PER_DAY = 86400
end

"""
Calculate the maximum monthly fractional contribution of glacier runoff to total river flux.

This function:
1. Loads glacier runoff and river flux data from NetCDF files
2. Aligns time series data and ensures consistent dimensions
3. Focuses on peak runoff periods (±2 months around maximum) for each river segment
4. Calculates the fraction of total river flux contributed by glacier runoff
5. Determines monthly averages and extremes (min/max) of glacier contribution
6. Adds results to river network data for spatial analysis

Returns a GeoDataFrame with river segments and their glacier runoff contribution metrics.
"""
@time begin
    # load data from netcdf files
    glacier_flux_nc = NCDataset(glacier_summary_riverflux_file)
    river_flux_nc = NCDataset(glacier_rivers_land_flux_path)

    # convert to DimensionalData
    glacier_runoff = GlobalGlacierAnalysis.nc2dd(glacier_flux_nc["runoff"]) / (1000 * u"kg/m^3")
    river_flux = GlobalGlacierAnalysis.nc2dd(river_flux_nc["flux"])
   
    # having issues with sortslices
    #@time glacier_flux = sortslices(glacier_flux, dims=:COMID)
    #@time river_flux = sortslices(river_flux, dims=:COMID)

    @assert issorted(dims(glacier_runoff,:COMID))

    p = sortperm(dims(river_flux,:COMID).val)
    river_flux = river_flux[:,p]

    # do this little hack to get COMID as ForwardOrdered
    river_flux = DimArray(parent(river_flux), (dims(river_flux, :Ti), Dim{:COMID}(val(dims(river_flux, :COMID).val))))

    # subset the data to the date interval
    glacier_runoff = glacier_runoff[dates4trend[1]..dates4trend[2], :]
    river_flux = river_flux[dates4trend[1]..dates4trend[2], :]

    glacier_Ti = dims(glacier_runoff, :Ti)
    river_Ti = dims(river_flux, :Ti)

    glacier_index = (glacier_Ti .>= minimum(dates4trend)) .& (glacier_Ti .<= maximum(dates4trend))
    river_index = (river_Ti .>= minimum(dates4trend)) .& (river_Ti .<= maximum(dates4trend))


    #TODO: WHY DOES RUNOFF NOT GO TO ZERO? This needs further investigation... 
    # COMID = 81020092 

    # require that gmax must be calculated with +/- 2 months peak runoff
    glacier_monthly = groupby(glacier_runoff, Ti => Bins(month, 1:12))
    glacier_monthly = mean.(glacier_monthly; dims=:Ti)
    dmonth = Dim{:Month}((collect(val(val((dims.(glacier_monthly, :Ti)))))[1]))

    glacier_monthly = vcat(parent.(glacier_monthly)...)    
    glacier_monthly = DimArray(glacier_monthly, (dmonth, dims(glacier_runoff, :COMID)))

    glacier_max_month = mapslices(argmax, glacier_monthly; dims=:Month)
    glacier_max_month = dropdims(glacier_max_month; dims=:Month)

    glacier_max = mapslices(maximum, glacier_monthly; dims=:Month)
    glacier_max = dropdims(glacier_max; dims=:Month)

    ti_month = month.(dims(glacier_runoff, :Ti))
    dt = 2;
    for comid in dims(glacier_monthly, :COMID)
    #comid = dims(glacier_monthly, :COMID)[2]
        month_max = glacier_max_month[At(comid)] 
        month_max = collect((month_max - dt):(month_max + dt))
        month_max[month_max .>12] .= month_max[month_max .>12] .- 12
        month_max[month_max .<1] .= month_max[month_max .<1] .+ 12
        index = .!in.(ti_month, Ref(month_max))
        glacier_runoff[index, At(comid)] .= 0 * u"m^3/s"
    end

    # calcualte the fraction of flux that is from glacier runoff
    glacier_fraction = @d glacier_runoff[findall(glacier_index), :] ./ (glacier_runoff[findall(glacier_index), :] .+ river_flux[findall(river_index), :])
    # calcualte the fraction of flux that is from glacier runoff
    glacier_fraction = @d glacier_runoff[findall(glacier_index), :] ./ (glacier_runoff[findall(glacier_index), :] .+ river_flux[findall(river_index), :])
    # calcualte the fraction of flux that is from glacier runoff
    glacier_fraction = @d glacier_runoff[findall(glacier_index), :] ./ (glacier_runoff[findall(glacier_index), :] .+ river_flux[findall(river_index), :])
    
    glacier_fraction[isnan.(glacier_fraction)] .= 0
    glacier_fraction = round.(Int8, glacier_fraction.*100)
    glacier_fraction_monthly = groupby(glacier_fraction, Ti => Bins(month, 1:12))
    glacier_fraction_monthly = mean.(glacier_fraction_monthly; dims=:Ti)
    glacier_fraction_monthly = cat(glacier_fraction_monthly...; dims=dims(glacier_fraction_monthly, :Ti))
    glacier_fraction_monthly = round.(Int8, glacier_fraction_monthly)
    dmonth = Dim{:month}(month.(dims(glacier_fraction_monthly, :Ti)))
    glacier_fraction_monthly = DimArray(glacier_fraction_monthly, (dmonth, dims(glacier_fraction_monthly, :COMID)))

    gmax_average = maximum(glacier_fraction_monthly, dims=:1)

    gmax_month = deepcopy(gmax_average)
    for i in eachindex(gmax_average)
        gmax_month[i] = findfirst(glacier_fraction_monthly[:,i] .== gmax_average[i])
    end

    glacier_fraction_monthly = groupby(glacier_fraction, Ti => month)
    gmax_maximum = deepcopy(gmax_average)
    gmax_minimum = deepcopy(gmax_average)

    for i in eachindex(gmax_month)
        foo = glacier_fraction_monthly[At(gmax_month[i])][:,i]
        gmax_maximum[i] = maximum(foo)
        gmax_minimum[i] = minimum(foo)
    end

    rivers = GeoDataFrames.read(glacier_rivers_path)
    rivers = rivers[:, Not("NextDownID", "up1", "glacier_table_row", "RGIId")]
    sort!(rivers, :COMID)
    @assert collect(dims(gmax_maximum , :COMID)) == rivers.COMID

    # add the glacier fraction maximum to the rivers dataframe
    rivers[!, :runoff_max_avg] = ustrip.(vec(parent(glacier_max)))
    rivers[!, :runoff_max_month] = ustrip.(vec(parent(glacier_max_month)))
    rivers[!, :gmax_max] = ustrip.(vec(parent(gmax_maximum)))
    rivers[!, :gmax_min] = ustrip.(vec(parent(gmax_minimum)))
    rivers[!, :gmax_range] = rivers[!, :gmax_max] .- rivers[!, :gmax_min]
    rivers[!, :gmax_avg] = ustrip.(vec(parent(gmax_average)))
    rivers[!, :gmax_month] = ustrip.(vec(parent(gmax_month)))

    # write to disk
    GeoDataFrames.write(glacier_summary_gmax_file , rivers)
end

"""
Create a point dataset of river sinks with discharge information.

This function:
1. Loads glacier river network data
2. Identifies terminal rivers (sinks) where rivers end
3. Converts river geometries to point centroids
4. Adds runoff information from the rivers dataframe
5. Incorporates discharge data from external sources
6. Categorizes sinks by type (ocean terminating, inland, etc.)
7. Calculates marker sizes based on runoff volume
8. Saves the resulting point dataset to disk

The output is a GeoPackage file containing point geometries for all river sinks
with associated attributes for visualization and analysis.
"""
@time begin #[90s]
    # load glacier river paths
    glacier_rivers = GeoDataFrames.read(glacier_rivers_path)

    # subset to terminal rivers
    sinks = glacier_rivers[glacier_rivers[!, :NextDownID] .== 0, :]
    sinks.geometry = GO.centroid.(sinks.geometry)

    # create a point dataset of river sinks
    @assert issorted(rivers.COMID)

    river_idex = searchsortedfirst.(Ref(rivers.COMID), sinks.COMID)
    sinks[!, :runoff_max_avg] = rivers[river_idex, :runoff_max_avg]

    # remove unneeded columns
    sinks = sinks[!, Not(["NextDownID", "up1", "glacier_table_row", "RGIId", "continent", "country", "lengthkm"])]
    sinks[:, :discharge] .= false;

    discharge = load(path2discharge, "discharge")
    discharge0 = DataFrame()
    discharge0[!, :geometry] = tuple.(discharge.longitude, discharge.latitude)
    discharge0[!, :COMID] .= 0
    discharge0[!, :discharge] .= true
    discharge0[!, :ocean_terminating] .= true
    discharge0[!, :runoff_max_avg] .= GlobalGlacierAnalysis.gt_per_yr_to_m3_per_s.(discharge.discharge_gtyr)

    # join the discharge to the river sinks
    sinks = vcat(sinks, discharge0)
    sinks[isnan.(sinks.runoff_max_avg), :runoff_max_avg] .= 0

    # sixe of marker
    sinks[!, :size] = log.(sinks.runoff_max_avg .+ 1)

    # type of marker
    sinks[!, :type] .= Int8(1)
    sinks[sinks.discharge, :type] .= Int8(3)
    sinks[.!sinks.ocean_terminating, :type] .= Int8(2)
    # write to disk
    sinks_path = replace(glacier_rivers_path, ".arrow" => "_sinks.pkg")
    GeoDataFrames.write(sinks_path, sinks)
end