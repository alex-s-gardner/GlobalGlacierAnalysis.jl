# Calculate the climatological maximum monthly fractional contribution of glacier runoff to total river flux [gmax]
begin
    using NCDatasets
    using GeoDataFrames
    using DimensionalData
    using Altim
    using Dates
    using Statistics
    import GeometryOps as GO

    paths = Altim.pathlocal

    glacier_rivers_runoff_qout_path = joinpath(paths[:river], "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier_runoff_qout.nc")
    glacier_rivers_land_flux_path = joinpath(paths[:river], "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier_Qs_acc_Qsb_acc.nc")
    
    glacier_rivers_path = joinpath(paths[:river], "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.arrow")

    sdate = Date(2000,1,1)
    edate = Date(2024,1,1)
    date_interval =  sdate..edate;

    SECONDS_PER_DAY = 86400
end


@time begin
    # load glacier river paths
    glacier_rivers = GeoDataFrames.read(glacier_rivers_path)
    glacier_latitude = getindex.(GO.centroid.(glacier_rivers.geometry), 2)
    dcomid = Dim{:COMID}(glacier_rivers.COMID)
    glacier_latitude = DimArray(glacier_latitude, (dcomid,))

    # load data from netcdf files
    glacier_flux_nc = NCDataset(glacier_rivers_runoff_qout_path)
    river_flux_nc = NCDataset(glacier_rivers_land_flux_path)

    # convert to DimensionalData
    glacier_runoff = Altim.nc2dd(glacier_flux_nc, "runoff")
    river_flux = Altim.nc2dd(river_flux_nc, "flux")

    # having issues with sortslices
    #@time glacier_flux = sortslices(glacier_flux, dims=:COMID)
    #@time river_flux = sortslices(river_flux, dims=:COMID)

    @assert issorted(dims(glacier_runoff,:COMID))

    p = sortperm(dims(river_flux,:COMID).val)
    river_flux = river_flux[:,p]

    # do this little hack to get COMID as ForwardOrdered
    river_flux = DimArray(parent(river_flux), (dims(river_flux, :Ti), Dim{:COMID}(val(dims(river_flux, :COMID).val))))

    # subset the data to the date interval
    glacier_runoff = glacier_runoff[date_interval, :]
    river_flux = river_flux[date_interval, :]

    glacier_Ti = dims(glacier_runoff, :Ti)
    river_Ti = dims(river_flux, :Ti)

    glacier_index = (glacier_Ti .>= minimum(river_Ti)) .& (glacier_Ti .<= maximum(river_Ti))
    river_index = (river_Ti .>= minimum(glacier_Ti)) .& (river_Ti .<= maximum(glacier_Ti))


    #TODO: WHY DOES RUNOFF NOT GO TO ZERO? This needs further investigation... 
    # COMID = 81020092 

    # require that gmax must be calculated with +/- 2 months peak runoff
    glacier_monthly = groupby(glacier_runoff, Ti => month)
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
        glacier_runoff[index, At(comid)] .= 0
    end

    # calcualte the fraction of flux that is from glacier runoff
    glacier_fraction = @d glacier_runoff[findall(glacier_index), :] ./ (glacier_runoff[findall(glacier_index), :] .+ river_flux[findall(river_index), :])
    # calcualte the fraction of flux that is from glacier runoff
    glacier_fraction = @d glacier_runoff[findall(glacier_index), :] ./ (glacier_runoff[findall(glacier_index), :] .+ river_flux[findall(river_index), :])
    # calcualte the fraction of flux that is from glacier runoff
    glacier_fraction = @d glacier_runoff[findall(glacier_index), :] ./ (glacier_runoff[findall(glacier_index), :] .+ river_flux[findall(river_index), :])
    glacier_fraction[isnan.(glacier_fraction)] .= 0
    glacier_fraction = round.(Int8, glacier_fraction.*100)
    glacier_fraction_monthly = groupby(glacier_fraction, Ti => month)
    glacier_fraction_monthly = mean.(glacier_fraction_monthly; dims=:Ti)
    glacier_fraction_monthly = vcat(glacier_fraction_monthly...)
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
    rivers[!, :runoff_max_avg] = vec(parent(glacier_max ))
    rivers[!, :runoff_max_month] = vec(parent(glacier_max_month))
    rivers[!, :gmax_max] =vec(parent(gmax_maximum))
    rivers[!, :gmax_min] =vec(parent(gmax_minimum))
    rivers[!, :gmax_range] = rivers[!, :gmax_max] .- rivers[!, :gmax_min]
    rivers[!, :gmax_avg] = vec(parent(gmax_average))
    rivers[!, :gmax_month] = vec(parent(gmax_month))
   

    alpha = zeros(size(rivers[!, :runoff_max_avg]))
    alpha[rivers[!, :runoff_max_avg] .> 1000] .= 1
    alpha[rivers[!, :runoff_max_avg] .< 1000] .= 0.8
    alpha[rivers[!, :runoff_max_avg] .< 100] .= 0.6
    alpha[rivers[!, :runoff_max_avg] .< 10] .= 0.4
    alpha[rivers[!, :runoff_max_avg] .< 1] .= 0.2
    rivers[!, :alpha] = Int8.(alpha * 100)


    # write to disk
    rivers_gmax_path = replace(glacier_rivers_path, ".arrow" => "_gmax.gpkg")
    GeoDataFrames.write(rivers_gmax_path, rivers)
end

