# ============================================================================
# gmax_point_figure.jl
#
# Calculate and visualize the climatological maximum monthly fractional contribution of glacier runoff to total river flux (gmax).
#
# This script:
# 1. Loads glacier runoff and river flux data from NetCDF files
# 2. Processes time series data for specific geographic locations
# 3. Creates visualizations showing:
#    - River flux components (land, snow, glacier)
#    - Glacier runoff as a fraction of total river flux
#    - The maximum monthly contribution (gmax) with seasonal patterns
# 4. Saves the resulting figures for each location
#
# The script handles data from 2000-2024 and focuses on key river basins near glaciated regions.
# ============================================================================

begin
    using NCDatasets
    using GeoDataFrames
    using DimensionalData
    import GlobalGlacierAnalysis as GGA
    using Dates
    using Statistics
    using CairoMakie
    using Unitful

    dates4plot = (Date(2000,1,1), Date(2024,10,1))
    seperate_out_snow = false;

    paths = GGA.pathlocal

    glacier_flux_path = joinpath(paths[:project_dir], "gardner2025_glacier_summary_riverflux.nc")

    glacier_rivers_land_flux_path = joinpath(paths[:river],"riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier_Qs_acc_Qsb_acc.nc")
    glacier_rivers_snow_flux_path = joinpath(paths[:river],"riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier_Qsm_acc.nc")
    land_path = paths[:river]
    glacier_rivers_path = joinpath(paths[:river], "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.arrow")
end


begin #[50s]
    # load data from netcdf files
    glacier_flux = NCDataset(glacier_flux_path)
    land_flux = NCDataset(glacier_rivers_land_flux_path)
    seperate_out_snow && (snow_flux = NCDataset(glacier_rivers_snow_flux_path))

    # convert to DimensionalData
    glacier_flux = GGA.nc2dd(glacier_flux["runoff"]) / (1000 * u"kg/m^3")
    land_flux = GGA.nc2dd(land_flux["flux"])
    seperate_out_snow && (snow_flux = GGA.nc2dd(snow_flux["flux"]))
    # having issues with sortslices
    #@time glacier_flux = sortslices(glacier_flux, dims=:COMID)
    #@time land_flux = sortslices(land_flux, dims=:COMID)

    p = sortperm(dims(glacier_flux,:COMID).val)
    glacier_flux = glacier_flux[:,p]

    p = sortperm(dims(land_flux,:COMID).val)
    land_flux = land_flux[:,p]

    seperate_out_snow && (p = sortperm(dims(snow_flux,:COMID).val))
    seperate_out_snow && (snow_flux = snow_flux[:,p])

    # do this little hack to get COMID as ForwardOrdered
    land_flux = DimArray(parent(land_flux), (dims(land_flux, :Ti), Dim{:COMID}(val(dims(land_flux, :COMID).val))))
    seperate_out_snow && (snow_flux = DimArray(parent(snow_flux), (dims(snow_flux, :Ti), Dim{:COMID}(val(dims(snow_flux, :COMID).val)))))

    #TODO: WHY DOES RUNOFF NOT GO TO ZERO? This needs further investigation... 
    # COMID = 81020092 

    locations = [
        (45047742, "Karachi, Pakistan"),
        (71017956, "Edmonton, Canada"),
        (71031599, "Calgary, Canada"),
        (71039426, "Lethbridge, Canada"),
        (78013077, "Vancouver, Canada"),
        (78016783, "Tacoma, USA"),
        (81014312, "Fairbanks, USA"),
        (67007399, "Lima, Peru"),
        (66010882, "Copiapó, Chile"),
        (66015561, "Temuco, Chile"),
        (65003342, "Mendoza, Argentina"),
        (66013257, "Santiago, Chile"),
        (65001686, "San Juan, Argentina"),
        (67007407, "Lima, Peru"), 
        (45011579, "Lahore, Pakistan"),
        (45008062, "Jammu, India"),
        (45022890, "Delhi, India"),
        (45043858, "Kanpur, India "),
        (45004712, "Peshawar, Pakistan"),
        (46044113, "Dushanbe, Tajikistan"),
        (46047000, "Tashkent, Uzbekistan"),
        (46041153, "Samarkand, Uzbekistan"),
        (46037217, "Tashkent, Uzbekistan"), 
        (47017377, "Ürümqi, China"),
        (47032572, "Jiayuguan, China"),
        (21001926, "Aosta, Italy"),    
        (21001122, "Geneva, Switzerland"),
        (21001840, "Lyon, France"),
        (21002710, "Grenoble, France"),
        (22027761, "Innsbruck, Austria"),
        (24003651, "Luleå, Sweden"),
        (25005999, "Mo i Rana, Norway")
    ]


    date_interval = dates4plot[1]..dates4plot[2]
end


#for location in locations
begin
location = locations[1]

    COMID = location[1]
    name = location[2]

    fig = Figure()

    xticks = collect(year(dates4plot[1]):5:(year(dates4plot[2]) + year(1)));

    
    seperate_out_snow && (y_snow = collect(snow_flux[date_interval, At(COMID)]) ./ 1000)

    ax1 = Axis(fig[1:2, 1]; xticklabelsvisible=false, xticksvisible=false, title = name, xticks = xticks)
    
    # plot land flux
    y = collect(land_flux[date_interval, At(COMID)])./1000

    seperate_out_snow && (y .-= y_snow)

    x = GGA.decimalyear.(collect(dims(land_flux, :Ti)[date_interval]))
    x = vcat(x[1],x, x[end])
    yunits = Unitful.unit(land_flux[1])
    y = vcat(0*yunits, y, 0*yunits)
    max_y = maximum(y)
    p = poly!(ax1, GGA.GI.Point.(ustrip(x), ustrip(y)); color=(:peru, 1), label="land")

    # plot snow flux
    if seperate_out_snow
        y = y_snow;
        x = GGA.decimalyear.(collect(dims(land_flux, :Ti)[date_interval]))
        x = vcat(x[1], x, x[end])
        yunits = Unitful.unit(y[1])
        y = vcat(0.0*yunits , y, 0.0*yunits)
        max_y = max(maximum(y),max_y)
        p = poly!(ax1, GGA.GI.Point.(x, ustrip(y)); color=(:gray, 0.5), label="snow")
    end

    y = collect(glacier_flux[date_interval,At(COMID)])./1000
    yunits = Unitful.unit(y[1])
    x = GGA.decimalyear.(collect(dims(glacier_flux, :Ti)[date_interval]))
    x = vcat(x[1],x, x[end])
    y = vcat(0*yunits, y, 0*yunits)
    poly!(ax1, GGA.GI.Point.(ustrip(x), ustrip(y)); color=(:blue, 0.5), label="glacier")

    axislegend(ax1; merge=true, backgroundcolor=(:white, 0), framevisible=false, position=:lt)
    #xlabel!(p.axis[1], "Year")
    ax1.ylabel = "river flux (m³/s) x 10⁻³"
    xlims!(ax1, minimum(xticks), maximum(xticks))
    ylims!(ax1, 0, ustrip(max_y))

    ax2 = Axis(fig[3, 1],  xticks = xticks)
    glacier_frac = glacier_flux[date_interval,At(COMID)] ./ (land_flux[date_interval,At(COMID)] .+ glacier_flux[date_interval,At(COMID)])*100
    x = GGA.decimalyear.(collect(dims(glacier_flux, :Ti)[date_interval]))
    gmax = groupby(glacier_frac, Ti => Bins(month, 1:12))
    gmax = mean.(gmax; dims=:Ti)
    gmax = cat(gmax...; dims=dims(gmax, :Ti))
    dTi = dims(gmax, :Ti)
    gmax_month = dTi[argmax(gmax)]
    gmax = round(Int8, maximum(gmax))

    foo = glacier_frac[month.(dims(glacier_frac, :Ti)) .== gmax_month]
    gmax_min = round(Int8, minimum(foo))
    gmax_max = round(Int8, maximum(foo))

    x = GGA.decimalyear.(collect(dims(glacier_flux, :Ti)[date_interval]))
    lines!(ax2, collect(x), collect(glacier_frac); color = (:blue, 0.5))
    lines!(ax2, [GGA.decimalyear(dates4plot[1]), GGA.decimalyear(dates4plot[end])], [gmax_min, gmax_min]; color = (:gray, 0.5), linestyle = :dash)
    lines!(ax2, [GGA.decimalyear(dates4plot[1]), GGA.decimalyear(dates4plot[end])], [gmax_max, gmax_max]; color=(:gray, 0.5), linestyle=:dash)
    lines!(ax2, [GGA.decimalyear(dates4plot[1]), GGA.decimalyear(dates4plot[end])], [gmax, gmax]; color=(:black, 1))
    ax2.ylabel = "glacier fraction [%]"

    xlims!(ax2, minimum(xticks), maximum(xticks))
    text!(GGA.decimalyear(dates4plot[end]), mean((gmax_max, gmax)), text="gmax = $gmax [$(Dates.monthname(gmax_month))] ", align=[:right, :center], color=(:black, 0.8), overdraw=true)
    display(fig)

    out_name = replace(split(name, ",")[1], " " => "_")
    
    output_path = joinpath(paths[:figures], "riverflux_$out_name.png")
    save(output_path, fig)
end
