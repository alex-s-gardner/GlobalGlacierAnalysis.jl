# Calculate gmaxclimatological maximum monthly fractional contribution of glacier runoff to total river flux
begin
    using NCDatasets
    using GeoDataFrames
    using DimensionalData
    using Altim
    using Dates
    using Statistics
    using CairoMakie

    paths = Altim.pathlocal

    glacier_rivers_runoff_qout_path = "/mnt/bylot-r3/data/rivers/MERIT_Hydro_v07_Basins_v01/riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier_runoff_qout.nc"
    
    glacier_rivers_land_flux_path = "/mnt/bylot-r3/data/rivers/MERIT_Hydro_v07_Basins_v01/riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier_Qs_acc_Qsb_acc.nc"
    glacier_rivers_snow_flux_path = "/mnt/bylot-r3/data/rivers/MERIT_Hydro_v07_Basins_v01/riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier_Qsm_acc.nc"
    land_path = joinpath(paths.data_dir, "rivers/MERIT_Hydro_v07_Basins_v01")
    glacier_rivers_path = joinpath(rivers_path, "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.arrow")
end

seperate_out_snow = true;
#@time begin
    # load data from netcdf files
    glacier_flux = NCDataset(glacier_rivers_runoff_qout_path)
    land_flux = NCDataset(glacier_rivers_land_flux_path)
    snow_flux = NCDataset(glacier_rivers_snow_flux_path)

    # convert to DimensionalData
    glacier_flux = Altim.nc2dd(glacier_flux, "flux")
    land_flux = Altim.nc2dd(land_flux, "flux")
    snow_flux = Altim.nc2dd(snow_flux, "flux")
    # having issues with sortslices
    #@time glacier_flux = sortslices(glacier_flux, dims=:COMID)
    #@time land_flux = sortslices(land_flux, dims=:COMID)

    p = sortperm(dims(glacier_flux,:COMID).val)
    glacier_flux = glacier_flux[:,p]

    p = sortperm(dims(land_flux,:COMID).val)
    land_flux = land_flux[:,p]

    p = sortperm(dims(snow_flux,:COMID).val)
    snow_flux = snow_flux[:,p]

    # do this little hack to get COMID as ForwardOrdered
    land_flux = DimArray(parent(land_flux), (dims(land_flux, :Ti), Dim{:COMID}(val(dims(land_flux, :COMID).val))))
    snow_flux = DimArray(parent(snow_flux), (dims(snow_flux, :Ti), Dim{:COMID}(val(dims(snow_flux, :COMID).val))))

    #TODO: WHY DOES RUNOFF NOT GO TO ZERO? This needs further investigation... 
    # COMID = 81020092 

    locations = [
        (45047742, "Kerachi, Pakistan"),
        (71017956, "Edmonton, Canada"),
        (71031599, "Calgary, Canada"),
        (71039426, "Lethbridge, Canada"),
        (78013077, "Vancouver, Canada"),
        (78016783, "Tacoma, USA"),
        (81014312,"Fairbanks, USA"),
        (67007399, "Lima, Peru"),
        (66010882,"Copiapó, Chile"),
        (66015561, "Temuco, Chile"),
        (65003342,"Mendoza, Argentina"),
        (66013257, "Santiago, Chile"),
        (65001686, "San Juan, Argentina"),
        (67007407,"Lima, Peru"), 
        (45011579, "Lahore, Pakistan"),
        (45008062, "Jammu, India"),
        (45022890, "Delhi, India"),
        (45043858, "Kanpur, India "),
        (45004712,"Peshawar, Pakistan"),
        (46044113, "Dushanbe, Tajikistan"),
        (46047000, "Tashkent, Uzbekistan"),
        (46041153, "Samarkand, Uzbekistan"),
        (46037217, "Tashkent, Uzbekistan"), 
        (47017377, "Ürümqi, China"),
        (47032572, "Jiayuguan, China"),
        (21001926, "Aosta, Italy"),    
        (21001122, "Geneva, Switzerland"),
        (21001840  , "Lyon, France"),
        (21002710, "Grenoble, France"),
        (22027761, "Innsbruck, Austria"),
        (24003651  , "Luleå, Sweden"),
        (25005999, "Mo i Rana, Norway")
    ]

    sdate = Date(2000,1,1)
    edate = Date(2023,1,1)
    date_interval =  sdate..edate;

for location in locations
#location = last(locations)
    COMID = location[1]
    name = location[2]

    fig = Figure()

    y_snow = collect(snow_flux[date_interval, At(COMID)]) ./ 1000

    ax1 = Axis(fig[1:2, 1]; xticklabelsvisible=false, title = name)
    
    # plot land flux
    y = collect(land_flux[date_interval,At(COMID)])./1000

    if seperate_out_snow
        y .-= y_snow
    end

    x = Altim.decimalyear.(collect(dims(land_flux, :Ti)[date_interval]))
    x = vcat(x[1],x, x[end])
    y = vcat(0.,y, 0.)
    max_y = maximum(y)
    p = poly!(ax1, Point.(x, y); color = (:orange, 1), label="land")

    # plot snow flux
    if seperate_out_snow
        y = y_snow;
        x = Altim.decimalyear.(collect(dims(land_flux, :Ti)[date_interval]))
        x = vcat(x[1], x, x[end])
        y = vcat(0.0, y, 0.0)
        max_y = max(maximum(y),max_y)
        p = poly!(ax1, Point.(x, y); color=(:gray, 0.5), label="snow")
    end

    y = collect(glacier_flux[date_interval,At(COMID)])./1000
    x = Altim.decimalyear.(collect(dims(glacier_flux, :Ti)[date_interval]))
    x = vcat(x[1],x, x[end])
    y = vcat(0.,y, 0.)
    poly!(ax1, Point.(x, y); color = (:blue, 0.5), label = "glacier")

    axislegend(ax1; merge = true, backgroundcolor=(:white,0), framevisible = false)
    #xlabel!(p.axis[1], "Year")
    ax1.ylabel = "river flux (m³/s) x 10³"
    xlims!(ax1, Altim.decimalyear(sdate), Altim.decimalyear(edate))
    ylims!(ax1, 0, max_y)

    ax2 = Axis(fig[3, 1])
    glacier_frac = glacier_flux[date_interval,At(COMID)] ./ (land_flux[date_interval,At(COMID)] .+ glacier_flux[date_interval,At(COMID)])*100
    x = Altim.decimalyear.(collect(dims(glacier_flux, :Ti)[date_interval]))
    gmax = groupby(glacier_frac, Ti => month)
    gmax = mean.(gmax; dims=:Ti)
    gmax_month = month.(collect(getindex.(val.(val.(dims.(gmax, :Ti))),1)))
    gmax = vcat(gmax...)
    gmax_month = gmax_month[argmax(gmax)]

    gmax = round(Int8, maximum(gmax))

    foo = glacier_frac[month.(dims(glacier_frac, :Ti)) .== gmax_month]
    gmax_min = round(Int8, minimum(foo))
    gmax_max = round(Int8, maximum(foo))

    x = Altim.decimalyear.(collect(dims(glacier_flux, :Ti)[date_interval]))
    lines!(ax2, x, collect(glacier_frac); color = (:blue, 0.5))
    lines!(ax2, [Altim.decimalyear(sdate), Altim.decimalyear(edate)], [gmax_min, gmax_min]; color = (:gray, 0.5), linestyle = :dash)
    lines!(ax2, [Altim.decimalyear(sdate), Altim.decimalyear(edate)], [gmax_max, gmax_max]; color=(:gray, 0.5), linestyle=:dash)
    lines!(ax2, [Altim.decimalyear(sdate), Altim.decimalyear(edate)], [gmax, gmax]; color=(:black, 0.8))
    text!(Altim.decimalyear(edate), mean((gmax_max, gmax)), text="gmax = $gmax ", align=[:right, :center], color=(:black, 0.8))
    ax2.ylabel = "glacier fraction [%]"

    xlims!(ax2, Altim.decimalyear(sdate), Altim.decimalyear(edate))
    display(fig)
end
