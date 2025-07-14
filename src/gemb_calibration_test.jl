import GlobalGlacierAnalysis as GGA
using DimensionalData
using CairoMakie
using DataFrames
using Dates
using Statistics
using ProgressMeter

path2runs_synthesized = ["/mnt/bylot-r3/data/binned/2deg/glacier_dh_best_nmad5_v01_filled_ac_p2_synthesized.jld2"]

discharge = GGA.global_discharge_filled(;
    surface_mask="glacier",
    discharge_global_fn=GGA.pathlocal[:discharge_global],
    gemb_run_id=4,
    discharge2smb_max_latitude=-60,
    discharge2smb_equilibrium_period=(Date(1979), Date(2000)),
    pscale=1,
    Δheight=0,
    geotile_width=2,
    force_remake_before=DateTime("2025-07-08T14") + GGA.local2utc
);

vars2load = ["runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze"]
gemb = GGA.gemb_ensemble_dv(; gemb_run_id=4, vars2load)

# load geotiles
geotile_grouping_min_feature_area_km2 = 100
force_remake_before = nothing
geotile_width = 2

seasonality_weights = 0.8:0.02:0.98
distance_from_origin_penalties = 0.01:0.01:0.2

date_interval = DateTime(2000, 4, 1).. DateTime(2024, 12, 31)

dws = Dim{:ws}(seasonality_weights)
dwd = Dim{:wd}(distance_from_origin_penalties)
amp = DimArray(zeros(dws, dwd))
trend = DimArray(zeros(dws, dwd))


#seasonality_weight = seasonality_weights[1]
seasonality_weight = GGA.seasonality_weight
 #  for distance_from_origin_penalty in distance_from_origin_penalties
    distance_from_origin_penalty = GGA.distance_from_origin_penalty

        begin
            GGA.gemb_calibration(
                path2runs_synthesized, 
                discharge,
                gemb;
                geotile_width,
                geotile_grouping_min_feature_area_km2,
                single_geotile_test=nothing,
                plots_show=false,
                plots_save=false,
                seasonality_weight,
                distance_from_origin_penalty,
                force_remake_before = DateTime("2026-07-09T14") + GGA.local2utc
            )

            # [14] Calibrate GEMB model to altimetry data for each synthesized geotile dataset [3 hrs]
            GGA.geotile_synthesis_gembfit_dv(
                path2runs_synthesized, 
                discharge,
                gemb;
                geotile_grouping_min_feature_area_km2=100, 
                geotile_width=2, 
                force_remake_before=DateTime("2026-07-09T14") + GGA.local2utc
            )

            path2ref = replace(path2runs_synthesized[1], ".jld2" => "_gembfit_dv.jld2")

            # let's inspect that data
            df = GGA.FileIO.load(path2ref, "geotiles")
            dgeotile = Dim{:geotile}(df.id)
        
        
            if true
                f = Figure();
                ax = Axis(f[1, 1]);
                date_interval = DateTime(2000, 4, 1).. DateTime(2024, 12, 31)
                ylims = (0,0)

                for varname in ["dm", "dm_altim"]
                    ddate = Dim{:date}(DataFrames.colmetadata(df, varname, "date"))
                    da = DimArray(reduce(hcat, df[!, varname])', (dgeotile, ddate))
                    ts = dropdims(sum(da, dims=:geotile), dims=:geotile)
                    ts .-= mean(ts[date=date_interval])
                    ylims = extrema(vcat(collect(extrema(ts[date=date_interval])), collect(ylims)))
                    
                    fit = GGA.ts_seasonal_model(ts, interval=(date_interval.left, date_interval.right))

                    label = "$(varname) - (trend: $(round(Int, fit.trend)) Gt/yr, amp: $(round(Int, fit.amplitude)) Gt, phase peak month: $(fit.phase_peak_month))"
                    lines!(ts; label=label)
                end

                axislegend(ax, position=:rt)
                xlims!(ax, DateTime(2000, 1, 1), DateTime(2024, 12, 31))
                ylims!(ax, ylims)

                display(f)
            else
                varname = "dm_altim"
                ddate = Dim{:date}(DataFrames.colmetadata(df, varname, "date"))
                da = DimArray(reduce(hcat, df[!, varname])', (dgeotile, ddate))
                ts = dropdims(sum(da, dims=:geotile), dims=:geotile)
                ts .-= mean(ts[date=date_interval])
                fit_altim = GGA.ts_seasonal_model(ts, interval=(date_interval.left, date_interval.right))

                varname = "dm"
                ddate = Dim{:date}(DataFrames.colmetadata(df, varname, "date"))
                da = DimArray(reduce(hcat, df[!, varname])', (dgeotile, ddate))
                ts = dropdims(sum(da, dims=:geotile), dims=:geotile)
                ts .-= mean(ts[date=date_interval])
                fit = GGA.ts_seasonal_model(ts, interval=(date_interval.left, date_interval.right))
                amp[ws=At(seasonality_weight), wd=At(distance_from_origin_penalty)] = fit.amplitude - fit_altim.amplitude
                trend[ws=At(seasonality_weight), wd=At(distance_from_origin_penalty)] = fit.trend - fit_altim.trend
            end
        end
    end
end

save(joinpath(GGA.pathlocal[:project_dir], "gemb_calibration_test.jld2"), Dict("amp" => amp, "trend" => trend))

amp0 = abs.(amp)
trend0 = abs.(trend)

cost0 = (amp0 /2) .+ trend0;
idx = argmin(cost0)
DimPoints(cost0)[idx]
heatmap(cost0, colorrange = (0, 30))
heatmap(trend0)





