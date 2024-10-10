#begin
    using Altim
    using FileIO
    using CairoMakie
    using DataFrames
    using Statistics
    geotile_width = 2;
    project_id = :v01
    daterange = (2000, 2024)

    final_data_dir = joinpath(Altim.pathlocal.data_dir, "altim_final", "$(geotile_width)_$(project_id)")
    final_filename = "dvdm.jld2"
    final_figure_dir = joinpath(Altim.pathlocal.data_dir, "altim_final", "$(geotile_width)_$(project_id)", "figures")

    df = load(joinpath(final_data_dir, final_filename), "df")

    df = Altim.dvdm_crop2dates!(df, daterange)

    # recompute trends after crop2dates!(df, daterange)
    df = Altim.geotile_dvdm_add_trend!(df; iterations=1000)

    dates = DataFrames.metadata(df, "date")
    decyear = Altim.decimalyear.(dates)

    rgi = "global"
    mission = "synthesis"
    var = "dm"
    index = findfirst((df.rgi .== rgi) .& (df.mission .== mission) .& (df.var .== var))
    nonnan = .!isnan.(df[index, :].mid)
    dmass_mid = round(Int, df[index, :].mid[findlast(nonnan)] - df[index, :].mid[findfirst(nonnan)])
    err = round(Int, sqrt((df[index, :].low[findfirst(nonnan)] - df[index, :].mid[findfirst(nonnan)]).^2 + (df[index, :].low[findlast(nonnan)] - df[index, :].mid[findlast(nonnan)]).^2))

# total glacier loss
println("mission = $mission, rgi = $rgi, variable = $var: $dmass_mid ± $err $(df[index,:unit]) between $(dates[findfirst(nonnan)]) & $(dates[findlast(nonnan)])")


# acceleration
println("mission = $mission, rgi = $rgi, variable = $var: $(round(df[index, :acceleration], digits = 1)) ± $(round(df[index, :acceleration_err], digits = 1)) $(df[index,:unit]) yr⁻² between $(dates[findfirst(nonnan)]) & $(dates[findlast(nonnan)])")


#
rgi = "global"
mission = "gemb"
var = "fac"
index = findfirst((df.rgi .== rgi) .& (df.mission .== mission) .& (df.var .== var))
nonnan = .!isnan.(df[index, :].mid)
dmass_mid = round(Int, df[index, :].mid[findlast(nonnan)] - df[index, :].mid[findfirst(nonnan)])
err = round(Int, sqrt((df[index, :].low[findfirst(nonnan)] - df[index, :].mid[findfirst(nonnan)]) .^ 2 + (df[index, :].low[findlast(nonnan)] - df[index, :].mid[findlast(nonnan)]) .^ 2))

# total loss of fac
println("mission = $mission, rgi = $rgi, variable = $var: $dmass_mid ± $err $(df[index,:unit]) between $(dates[findfirst(nonnan)]) & $(dates[findlast(nonnan)])")

# offset by fac gain in antarctic
rgi = "rgi19"
index = findfirst((df.rgi .== rgi) .& (df.mission .== mission) .& (df.var .== var))
nonnan = .!isnan.(df[index, :].mid)
dmass_mid = round(Int, df[index, :].mid[findlast(nonnan)] - df[index, :].mid[findfirst(nonnan)])
err = round(Int, sqrt((df[index, :].low[findfirst(nonnan)] - df[index, :].mid[findfirst(nonnan)]) .^ 2 + (df[index, :].low[findlast(nonnan)] - df[index, :].mid[findlast(nonnan)]) .^ 2))



# 


# compare to ice sheets
begin
    grace = Altim.read_grace_rgi()
    df = load(joinpath(final_data_dir, final_filename), "df")
    dates_grace = Altim.datenum2date.(grace["rgi5"]["dM_gt_mdl_fill_date"]);
    dates_synth = DataFrames.metadata(df, "date")
    notnan = .!isnan.(df[findfirst(df.mission .== "synthesis"), :].mid)
    synth_range = extrema(dates_synth[notnan])

    notnan = .!isnan.(grace["rgi5"]["dM_gt_mdl_fill"])
    grace_range = extrema(dates_grace[notnan])

    index_synth = (dates_synth .> max(synth_range[1], grace_range[1])) .& (dates_synth .< min(synth_range[2], grace_range[2]))

    index_grace = (dates_grace .> max(synth_range[1], grace_range[1])) .& (dates_grace .< min(synth_range[2], grace_range[2]))


    daterange_synth = Altim.decimalyear.(extrema(dates_synth[index_synth]))

    y = grace["rgi5"]["dM_gt_mdl_fill"][index_grace] .+ grace["rgi19"]["dM_gt_mdl_fill"][index_grace]
    x = Altim.decimalyear.(Altim.datenum2date.(grace["rgi5"]["dM_gt_mdl_fill_date"][index_grace]))
    x = x .- mean(x)
    
    icesheet_fit = Altim.curve_fit(Altim.model3, x, y, Altim.p3).param
    icesheet_trend = icesheet_fit[2]
    icesheet_acceleration = icesheet_fit[3]

    df = Altim.dvdm_crop2dates!(df, daterange_synth)
    df = df[(df.mission.=="synthesis").& (df.var.=="dm"), :]
    df = Altim.geotile_dvdm_add_trend!(df; iterations=1000)


    icesheet_trend = round(Int, icesheet_trend - (df[findfirst(df.rgi .== "rgi5"), :trend] + df[findfirst(df.rgi .== "rgi19"), :trend]))
    icesheet_acceleration = round(icesheet_acceleration - (df[findfirst(df.rgi .== "rgi5"), :acceleration] + df[findfirst(df.rgi .== "rgi19"), :acceleration]), digits = 2)


    glacier_trend = round(Int, df[findfirst(df.rgi .== "global"), :trend])
    glacier_acceleration = round(df[findfirst(df.rgi .== "global"), :acceleration], digits=2)

    println("ice sheet $(icesheet_trend) vs. glacier $(glacier_trend) Gt yr⁻¹")
    println("ice sheet $(icesheet_acceleration) vs. glacier $(glacier_acceleration) Gt yr⁻²")
end