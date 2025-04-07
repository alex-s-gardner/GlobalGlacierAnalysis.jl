"""
    regional_results.jl

This script processes and analyzes glacier mass change data at regional scales.

# Main components:
1. Data aggregation and uncertainty estimation:
   - Loads multiple model runs and aggregates geotile data into regional time series
   - Calculates uncertainties from ensemble spread (95% confidence intervals)
   - Incorporates GRACE satellite data for comparison
   - Reformats data for consistent plotting and analysis

2. Glacier discharge and mass change analysis:
   - Processes glacier discharge data by region
   - Maps discharge points to RGI regions using spatial joins
   - Identifies endorheic (inland-draining) glaciers
   - Computes mass change rates for each glacier and region

3. Trend and acceleration analysis:
   - Fits trends and seasonal cycles to time series
   - Identifies regions with significant acceleration/deceleration
   - Calculates mass change under different scenarios (with/without firn air content correction)
   - Compares altimetry-derived trends with GRACE satellite measurements

4. Visualization and reporting:
   - Generates multi-region plots of mass change
   - Calculates statistics on glacier contribution to sea level rise - Analyzes endorheic basin contributions
   - Reports regional mass change trends with uncertainties

# Key outputs:
- Regional time series of glacier mass change
- Trend and acceleration estimates with uncertainties
- Comparison between different measurement techniques
- Analysis of endorheic vs. ocean-draining glacier contributions
"""

# NOTE: errors are at the 95% confidence interval from the reference run for this study and 95% confidence interval from GRACE
begin
    using Altim
    using DataFrames
    using FileIO
    using CairoMakie
    using DimensionalData
    using Statistics
    using Dates
    using DataInterpolations
    using LsqFit
    using ColorSchemes
    using GeoDataFrames
    using StatsBase
    using FlexiJoins
    import GeometryOps as GO

    # to include in uncertainty
    project_id = ["v01"]
    surface_mask = ["glacier", "glacier_rgi7"]
    dem_id = ["best", "cop30_v2"]
    curvature_correct = [false, true]
    amplitude_correct = [true]
    binning_method = ["median", "meanmadnorm5", "meanmadnorm3"]
    paramater_set = [1, 2, 3, 4]
    binned_folder = ["/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"]

    dates2plot =DateTime(1980, 1, 1):Month(1):DateTime(2025, 1, 15) # also dates used for trend fitting
    dates4trend = DateTime(2000, 1, 1):Month(1):DateTime(2025, 1, 15) # use these date for comparing individual model runs to GRACE (only need for rmse anlysis)
    #dates2plot = DateTime(2002, 4, 1):Month(1):DateTime(2025, 1, 15) # use these date for comparing individual model runs to GRACE (only need for rmse anlysis)
    
    paths = Altim.pathlocal
    reference_run = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_dh_best_cc_meanmadnorm3_v01_filled_ac_p1_synthesized.jld2";
    perglacier_reference_run = replace(reference_run, ".jld2" => "_perglacier.jld2")
    path2discharge = joinpath(paths.data_dir, "GlacierOutlines/GlacierDischarge/global_glacier_discharge.jld2")
    path2rgi_regions = paths.rgi6_regions_shp

    glacier_rivers_path = joinpath(paths.river, "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.arrow")

    param_nt = (;project_id, surface_mask, dem_id, curvature_correct, amplitude_correct, binning_method, paramater_set, binned_folder)
    params = Altim.ntpermutations(param_nt)

    # only include files that exist
    path2runs = String[]
    for param in params
        binned_aligned_file = Altim.binned_aligned_filepath(; param...)
        if isfile(binned_aligned_file)
            push!(path2runs, binned_aligned_file)
        end
    end
    path2runs = replace.(path2runs, "aligned.jld2" => "synthesized.jld2")

end


# load discharge for each RGI [<1s]
@time discharge_rgi = Altim.discharge_rgi(path2discharge, path2rgi_regions)

# load results for all runs for each RGI [18s]
@time runs_rgi = Altim.runs2rgi(path2runs)

# fit trends and acceleration for each model and RGI [33s]
@time runs_rgi_fits = Altim.rgi_trends(runs_rgi, discharge_rgi, dates4trend)


# create df for reference run

varnames = ["fac"  "refreeze"  "rain"  "acc"  "melt"  "runoff"  "dm_altim"  "dm"  "dv_altim"  "smb"  "dv"  "ec"]

    binned_synthesized_file = reference_run
    binned_synthesized_dv_file = replace(binned_synthesized_file, ".jld2" => "_gembfit_dv.jld2")

    geotiles0 = FileIO.load(binned_synthesized_dv_file, "geotiles")

    # group by rgi and sum
    geotiles_reg = groupby(geotiles0, :rgiid)

    # sum across timeseries
    regions = DataFrames.combine(geotiles_reg, varnames .=> Ref ∘ sum; renamecols=false)

    # add an HMA region 
    index = (regions[:, :rgiid] .<= 15) .& (regions[:, :rgiid] .>= 13)
    foo = DataFrames.combine(regions[index, :], varnames .=> Ref ∘ sum; renamecols=false)
    foo.rgiid .= 98
    regions = append!(regions, foo)

    # add global region 
    index = (regions[:, :rgiid] .<= 19) .& (regions[:, :rgiid] .>= 1)
    foo = DataFrames.combine(regions[index, :], vcat("rgiid", varnames) .=> Ref ∘ sum; renamecols=false)
    foo.rgiid .= 99
    regions = append!(regions, foo)

    # copy metadata
    for varname in names(regions)[2:end]
        cmetadata = colmetadata(geotiles0, varname, "date")
        regions = colmetadata!(regions, varname, "date", cmetadata; style=:note)
    end

    # add error columns to regions
    for varname in varnames
    #varname = "acc"
        varname_err = "$(varname)_err"
        var0 = regions0[varname]

        # loop for each rgi and remove mean to center data
        for rgi in drgi
            var1 = var0[:, At(rgi), :]
            rowrange, colrange = Altim.validrange(.!isnan.(var1))

            delta = collect(mean(var1[:, colrange], dims=2))
            var0[:, At(rgi), :] .-= delta * ones(1, size(var1, 2))
        end

        # remove reference run from ensemble for error calculation
        var1 = deepcopy(var0)
        for run in drun
            var1[At(run), :, :] .-= var0[At(reference_run), :, :]
        end

        ddate = dims(var1, :date)
        err0 = var1[1,:,:]
        err0 .= NaN
        for rgi in drgi
            for date in ddate
                foo = abs.(var1[:, At(rgi), At(date)]);
                valid = .!isnan.(foo)
                if any(valid)

                    ## THIS IS WHERE THE ERROR IS DEFINED [symetric 95% confidence interval from reference run]
                    err0[At(rgi), At(date)] = quantile(abs.(var1[valid, At(rgi), At(date)]), 0.95)
                    ##
                end
            end
        end
  
        regions[!, varname_err] = eachrow(collect(err0))
        regions = colmetadata!(regions, varname_err, "date", colmetadata(regions, varname, "date"); style=:note)
    end

    # read in GRACE data
    grace = Altim.read_grace_rgi(; datadir=setpaths()[:grace_rgi])

    regions[!, :dm_grace] = [fill(NaN, length(grace["rgi1"]["dM_gt_mdl_fill"])) for _ in 1:nrow(regions)]
    gdate = vec(Altim.datenum2date.(grace["rgi1"]["dM_gt_mdl_fill_date"]))
    regions = colmetadata!(regions, :dm_grace, "date", gdate; style=:note)
    regions[!, :dm_grace_err] = [fill(NaN, length(grace["rgi1"]["dM_gt_mdl_fill"])) for _ in 1:nrow(regions)]
    regions = colmetadata!(regions, :dm_grace_err, "date", gdate; style=:note)

    for r in eachrow(regions)
        #r = first(eachrow(regions))
        if r.rgiid == 98
            rgi = "HMA"
        elseif r.rgiid == 99
            rgi = "Global"
        else
            rgi = "rgi$(r.rgiid)"
        end

        haskey(grace, rgi) || continue

        r.dm_grace = vec(grace[rgi]["dM_gt_mdl_fill"])

        ## multiply by 2 to get 95% confidence interval
        r.dm_grace_err = vec(grace[rgi]["dM_sigma_gt_mdl_fill"]) * 2
    end

    derror = Dim{:error}([false, true])
    dvarname = Dim{:varname}(collect(keys(regions0_fit)))
    regions_fit = fill(NaN, dvarname, drgi, dparameter, derror)
    for varname in dvarname
        for rgi in drgi
            for param in dparameter
                v = regions0_fit[varname][At(binned_synthesized_file), At(rgi), At(param)]
                regions_fit[At(varname), At(rgi), At(param), At(false)] = v

                err0 = quantile(abs.(regions0_fit[varname][:, At(rgi), At(param)] .- v), 0.95)
                regions_fit[At(varname), At(rgi), At(param), At(true)] = err0
            end
        end
    end

    # reformat for plotting
    varnames = setdiff(names(regions), ["rgiid"])
    dates_decyear = Altim.decimalyear.(dates2plot)
    df = DataFrame()

    for varname in varnames
        #varname = "discharge"

        date0 = Altim.decimalyear.(colmetadata(regions, varname, "date"))
        v = fill(NaN, length(dates_decyear))

        if Base.contains(varname, "altim")
            mission = "synthesis"
        elseif Base.contains(varname, "grace")
            mission = "grace"
        else
            mission = "gemb"
        end

        if Base.contains(varname, "dm")
            if !Base.contains(varname, "err")
                varname_out = "dm"
            else
                varname_out = "dm_err"
            end
        elseif Base.contains(varname, "dv")
            if !Base.contains(varname, "err")
                varname_out = "dv"
            else
                varname_out = "dv_err"
            end
        else
            varname_out = varname
        end

        for r in eachrow(regions)
            #r = first(eachrow(regions))
            model = DataInterpolations.LinearInterpolation(r[varname], date0)
            index = (dates_decyear .> minimum(date0)) .& (dates_decyear .< maximum(date0))
            v[index] = model(dates_decyear[index])
            varname_mid = "$(varname)"
            df = vcat(df, DataFrame("rgi" => r[:rgiid], "mission" => mission, "var" => varname_out, "mid" => Ref(copy(v))))
        end
    end

    # set GRACE to NaN for Greenland and Antarctica
    index = findfirst((df.rgi .== 5) .& (df.mission .== "grace"))
    df[index, :mid] .= NaN
    index = findfirst((df.rgi .== 19) .& (df.mission .== "grace"))
    df[index, :mid] .= NaN

    df[!, :low] = deepcopy(df[!, :mid])
    df[!, :high] = deepcopy(df[!, :mid])
    df[!, :unit] .= "gt"
    metadata!(df, "date", dates2plot; style=:note)

    gdf = groupby(df, "rgi")

    for g in gdf
        #g = gdf[3]
        for mission in unique(g.mission)
            #mission = unique(g.mission)[1]
            mission_index = g.mission .== mission

            vars = unique(g[mission_index, :var])
            vars = vars[.!Base.contains.(vars, Ref("err"))]

            for var0 in vars
                #var0 = vars[2]
                var_index = g.var .== var0
                ind_mid = findfirst(var_index .& mission_index)
                ind_err = findfirst((g.var .== "$(var0)_err") .& mission_index)
                if !isnothing(ind_err)
                    g[ind_mid, :low] .-= g[ind_err, :mid]
                    g[ind_mid, :high] .+= g[ind_err, :mid]
                end
            end
        end
    end

    # remove "dm_grace_err"
    df = df[.!Base.contains.(df.var, Ref("err")), :]

    # align to altim
    gdf = groupby(df, "rgi")
    for g in gdf
        #g = gdf[1]
        index_var = g.var .== "dm"

        grace0 = g[findfirst((g.mission .== "grace") .& index_var), :mid]
        altim0 = g[findfirst((g.mission .== "synthesis") .& index_var), :mid]
        gemb0 = g[findfirst((g.mission .== "gemb") .& index_var), :mid]

        index = .!isnan.(altim0)
        for ts in (grace0, gemb0)
            index2 = .!isnan.(ts)
            if any(index2)
                index = index .& index2
            end
        end

        delta = mean(altim0[index])
        g[findfirst((g.mission .== "synthesis") .& index_var), :mid] .-= delta
        g[findfirst((g.mission .== "synthesis") .& index_var), :low] .-= delta
        g[findfirst((g.mission .== "synthesis") .& index_var), :high] .-= delta

        delta = mean(grace0[index])
        g[findfirst((g.mission .== "grace") .& index_var), :mid] .-= delta
        g[findfirst((g.mission .== "grace") .& index_var), :low] .-= delta
        g[findfirst((g.mission .== "grace") .& index_var), :high] .-= delta

        delta = mean(gemb0[index])
        g[findfirst((g.mission .== "gemb") .& index_var), :mid] .-= delta
        g[findfirst((g.mission .== "gemb") .& index_var), :low] .-= delta
        g[findfirst((g.mission .== "gemb") .& index_var), :high] .-= delta
    end
end


# 2. Processes glacier outlines and river data:
#    - Loads glacier outlines and river network data
#    - Identifies endorheic (inland-draining) glaciers
#    - Computes mass change rates for each glacier
#    - Aggregates results by RGI region
#
# 3. Calculates mass change trends:
#    - Computes trends during GRACE/altimetry overlap period
#    - Calculates mass change for different scenarios:
#      - With/without static firn air content correction
#      - With/without endorheic basin correction
#    - Stores results in multi-dimensional array
#
# Key outputs:
# - discharge_rgi: DataFrame with discharge by RGI region
# - glacier_dm_gt: DataFrame with mass change by RGI region
# - dm_gt: Array of mass change trends under different scenarios

begin

    # load in glacier outlines
    glaciers = load(perglacier_reference_run, "glaciers")

    # load in river flux data
    river_flux = GeoDataFrames.read(glacier_rivers_path)

    # select only terminating rivers
    river_flux = river_flux[river_flux[:, :NextDownID] .== 0, :]

    # add endorheic column to glaciers
    glaciers[!, :endorheic] .= false

    ## THIS SECTION IS VERY VERY SLOW
    Threads.@threads for g in eachrow(glaciers)
        for r in eachrow(river_flux)
            if g.RGIId in r.RGIId
                g.endorheic = !r.ocean_terminating
            end
        end
    end

    # compute mass change for each region for endorheic and glaciers

    # now compute dh rate for all glaciers
    varname = "dh"
    date = collect(dims(glaciers[1, varname], :date))

    # find the index of valid data
    index = .!isnan.(glaciers[1, varname])
    datelimits = extrema(date[index]);
    DataFrames.colmetadata!(glaciers, varname, "date", date; style=:note)

    # fit dh trend to all glaciers
    Altim.df_tsfit!(glaciers, [varname]; progress=true, datelimits)

    # fit dh trend to all glaciers
    varname = "fac"
    date = collect(dims(glaciers[1, varname], :date))
    DataFrames.colmetadata!(glaciers, varname, "date", date; style=:note)
    Altim.df_tsfit!(glaciers, [varname]; progress=true, datelimits)

    # add rgi column
    glaciers[!, :rgi] = [parse(Int16, g[7:8]) for g in glaciers.RGIId]

    # (dh - fac)/1000 * area = dm_gt
    glaciers[!, :dm_gt] = (glaciers.dh_trend .- glaciers.fac_trend) / 1000 .* sum.(glaciers.area_km2)

    # derive statistics for each region
    dendorheic_only = Dim{:endorheic_correction}([false,true])
    glacier_dm_gt = zeros(drgi, dendorheic_only)

    for rgi in setdiff(collect(drgi), [98,99])
        index = glaciers.rgi .== rgi
        endorheic = glaciers[:, :endorheic] .& index;
        glacier_dm_gt[At(rgi), At(false)] = sum(glaciers[index, "dm_gt"])
        glacier_dm_gt[At(rgi), At(true)] = sum(glaciers[endorheic, "dm_gt"]) 
    end

    # add HMA and global
    glacier_dm_gt[At(98),:] = sum(glacier_dm_gt[13..15, :], dims = :rgi)

    # add global
    glacier_dm_gt[At(99),:] = sum(glacier_dm_gt[1..19, :], dims = :rgi)
end

begin
    # compute trends on df values for overlaping grace and synthesis period
    valid = .!isnan.(df[findfirst(df.mission .== "synthesis"), :mid]) .& .!isnan.(df[findfirst(df.mission .== "grace"), :mid])
    dates = DataFrames.metadata(df, "date")
    colmetadata!(df, "mid", "date", dates; style=:note)
    datelimits = extrema(DataFrames.metadata(df, "date")[valid])
    Altim.df_tsfit!(df, [:mid]; progress=true, datelimits)

    # compute dm_gt for each rgi, static_fac, and endorheic_correction
    drgi = Dim{:rgi}(unique(df.rgi))
    dstatic_fac = Dim{:fac_static}([false, true])
    dendorheic_correction = Dim{:endorheic_correction}([false,true])

    dm_gt = fill(NaN,(dendorheic_correction, dstatic_fac, drgi))
    grace_dm_gt = fill(NaN,(drgi))

    for rgi in drgi
    #rgi = first(drgi)
        index2 = findfirst((df.mission .== "grace") .& (df.rgi .== rgi) .& (df.var .== "dm"))
        if !isnothing(index2)
            grace_dm_gt[At(rgi)] = df[index2, :mid_trend]
        end

        for static_fac in dstatic_fac
            for endorheic_correction in dendorheic_correction
                
                scale = 1.0
                if endorheic_correction
                    scale = (glacier_dm_gt[At(rgi), At(false)] - glacier_dm_gt[At(rgi), At(true)] ) / glacier_dm_gt[At(rgi), At(false)]
                end

                if static_fac
                    index1 = findfirst((df.mission .== "synthesis") .& (df.rgi .== rgi) .& (df.var .== "dv"))
                    scale = scale * 0.85
                else
                    index1 = findfirst((df.mission .== "synthesis") .& (df.rgi .== rgi) .& (df.var .== "dm"))
                end

                println(rgi)
                dm_gt[At(endorheic_correction), At(static_fac), At(rgi)] = df[index1, :mid_trend] * scale
            end
        end
    end
end

# find regions with significant acceleration in mass change
regions_with_acceleration = collect(dims(regions_fit, :rgi)[regions_fit[At("dm_altim"), :, At("acceleration"), 1] .<  -(regions_fit[At("dm_altim"), :, At("acceleration"), 2])])
println("Regions with significant acceleration:")
for rgi in regions_with_acceleration
    println("$(Altim.rgi2label[Altim.rginum2txt[rgi]])")
end

regions_with_decceleration = collect(dims(regions_fit, :rgi)[regions_fit[At("dm_altim"), :, At("acceleration"), 1] .> (regions_fit[At("dm_altim"), :, At("acceleration"), 2])])
println("Regions with decceleration:")
for rgi in regions_with_decceleration
    println("$(Altim.rgi2label[Altim.rginum2txt[rgi]])")
end


# plot altimetry results for RGI regions
outfile = joinpath(paths[:figures], "regional_dm_altim.png")
exclude_mission = (df.mission .== "gemb") .| (df.mission .== "grace")
exclude_regions = [13, 14, 15, 99]

f, plot_order_rgi, offset, ylims = Altim.plot_multiregion_dvdm(df[.!exclude_mission, :];
    variable="dm",
    featured_mission="synthesis",
    regions=setdiff(unique(df.rgi), exclude_regions),
    showlines=false,
    showmissions=true,
    fontsize=15,
    cmap=:Dark2_4,
    regions_ordered=false,
    region_offsets=nothing,
    ylims=nothing,
    title=nothing,
    palette=nothing,
    delta_offset=nothing,
)
display(f)
CairoMakie.save(outfile, f)

# plot altimetry results with GRACE data for RGI regions
outfile = joinpath(paths[:figures], "regional_dm_altim_grace.png")
exclude_mission = (df.mission .== "gemb")
exclude_regions = [13, 14, 15, 99]

f, plot_order_rgi, offset, ylims = Altim.plot_multiregion_dvdm(df[.!exclude_mission, :];
    variable="dm",
    featured_mission="synthesis",
    regions=setdiff(unique(df.rgi), exclude_regions),
    showlines=false,
    showmissions=true,
    fontsize=15,
    cmap=:Dark2_4,
    regions_ordered=false,
    region_offsets=nothing,
    ylims=nothing,
    title=nothing,
    palette=nothing,
    delta_offset=nothing,
)
display(f)
CairoMakie.save(outfile, f)

# plot altimetry results with gemb data for RGI regions
outfile = joinpath(paths[:figures], "regional_dm_altim_gemb.png")
exclude_mission = (df.mission .== "grace")
exclude_regions = [13, 14, 15, 99]

f, plot_order_rgi, offset, ylims = Altim.plot_multiregion_dvdm(df[.!exclude_mission, :];
    variable="dm",
    featured_mission="synthesis",
    regions=setdiff(unique(df.rgi), exclude_regions),
    showlines=false,
    showmissions=true,
    fontsize=15,
    cmap=:Dark2_4,
    regions_ordered=false,
    region_offsets=nothing,
    ylims=nothing,
    title=nothing,
    palette=nothing,
    delta_offset=nothing,
)
display(f)
CairoMakie.save(outfile, f)


n = 6;
large_region_glaciers = reverse(plot_order_rgi[(end-(n-1)):end]);

large_region_dm_gt = sum(glacier_dm_gt[At(large_region_glaciers), At(false)]) / glacier_dm_gt[At(99), At(false)]
println("$(n) largest regions comprise $(round(large_region_dm_gt, digits=2)) of total glacier loss. Regions: $(large_region_glaciers)")

# plot all iterations for a single region
var0 = "acc"
rgi = 1

p = lines(regions0[var0][1, At(rgi), :]);
for i = 2:length(regions0[var0][:, At(rgi), 1])
    lines!(regions0[var0][i, At(rgi), :])
end
display(p)


# create table with rgi results
begin
    exclude_regions = [13, 14, 15]
    varnames = ["acc", "runoff", "dm_altim", "dm", "ec", "dv_altim", "fac", "runoff_eq", "runoff_uns_perc"]
    varnames = ["dm_altim", "dm", "acc", "runoff", "ec", "fac", "runoff_eq", "runoff_uns_perc"]
    params = ["trend", "acceleration", "amplitude", "phase"]

    # constuct a DataFrame with summary results
    rgis = setdiff(drgi, exclude_regions)
    rgi_labels = Altim.rginum2label.(rgis)

    regional_results = DataFrames.DataFrame(rgi = rgis, region_name = rgi_labels)

    for varname in varnames
        println("    $(varname)")
        for param in params
            if param in ["acceleration"]
                unit = "Gt_per_decade"
            elseif param in ["trend"]
                unit = "Gt_yr"
            elseif param in ["amplitude"]
                unit = "Gt"
            elseif param in ["phase"]
                unit = "day_of_max"
            end

            column_label = Symbol("$(varname)_$(param)_[$unit]")
            regional_results[!, column_label] .= ""
        end
    end

    # add non standard columns
    regional_results[!, :discharge_Gt_yr] .= ""
    regional_results[!, :equilibrium_runoff_Gt_yr] .= ""
    regional_results[!, :non_sustainable_runoff_perc] .= ""

    for (i, rgi) in enumerate(rgis)
        println("$(Altim.rgi2label[Altim.rginum2txt[rgi]])")
        for varname in varnames
            println("    $(varname)")
            for param in params
                v = regions_fit[At(varname), At(rgi), At(param), At(false)]
                err = regions_fit[At(varname), At(rgi), At(param), At(true)]
                if param in ["phase"]
                    v = round(Int, v)
                    err = round(Int, err)
                else
                    v = round(v, digits=1)
                    err = round(err, digits=1)
                end

                if param in ["acceleration"]
                    v = v *10
                    err = err * 10
                    unit = "Gt_per_decade"
                elseif param in ["trend"]
                    unit = "Gt_yr"
                elseif param in ["amplitude"]
                    unit = "Gt"
                elseif param in ["phase"]
                    unit = "day_of_max"
                end
               
                column_label = Symbol("$(varname)_$(param)_[$unit]")
                regional_results[i, column_label] = "$(v) ± $(err)"
            end
        end
    end

    regional_results[!, :discharge_Gt_yr] .= string.(round.(discharge_rgi[At(rgis)], digits=1))
end


# dump table to console
colnames = names(regional_results)
last_set = ceil(Int, (length(colnames) - 2) / 4 - 1)
for i = 0:last_set
    printstyled("\n------------------------------------------------------- $(colnames[3 .+ (4 * i)])-------------------------------------------------------\n", color=:yellow)
    if (i == last_set) && (length(colnames) > (last_set+4))
        show(regional_results[:, [1:2; (3 .+ (4*i)):end]],allrows=true, allcols=true)
    else
        show(regional_results[:,[1:2; (3:6).+(4*i)]],allrows=true, allcols=true)
    end
end


# fraction of glacier loss that is endorheic
endorheic_frac = glacier_dm_gt[At(99), At(true)] /  glacier_dm_gt[At(99), At(false)]
println("glacier loss that flows to ocean = $(round(Int,(1-endorheic_frac)*100))%")

hma_frac_of_total_endorheic = glacier_dm_gt[At(98), At(true)] / glacier_dm_gt[At(99), At(true)];
println("HMA comprises $(round(Int,hma_frac_of_total_endorheic*100))% of total endorheic loss")

amp_static_fac = round(Int, regions_fit[At("dv"), At(99), At("amplitude"), At(false)]) * 0.85;
amp_dynamic_fac = round(Int, regions_fit[At("dm"), At(99), At("amplitude"), At(false)]);
println("Global amplitude would be $(round(Int,((amp_static_fac - amp_dynamic_fac) ./ amp_dynamic_fac *100)))% larger if static fac was applied")

# GRACE rgi regions
rgis = collect([1:4; 6:12; 16:18; 98])

function resample(colorscheme::Symbol, n)
    newscheme = ColorSchemes.ColorScheme(
        get(ColorSchemes.colorschemes[colorscheme], (0:n) / n))
    return newscheme
end

colormap = :Dark2_4;
cmap = resample(colormap, length(plot_order_rgi));
# set_theme!(Theme(palette=(color=resample(colormap, length(rgis)),),))

begin
    f = Figure()
    output_dir = joinpath(paths[:figures], "regional_dm_altim_grace.png")
    title = "glacier mass change [Gt/yr]"
    ax = Axis(f[1, 1]; 
        ylabel = "this study [Gt/yr]",
        xlabel = "GRACE /FO [Gt/yr]",
        title=title
    )

    for (i, rgi) in enumerate(plot_order_rgi)
        if (rgi == 5) || (rgi == 19)
            continue
        end
        scatter!(grace_dm_gt[At(rgi)], dm_gt[At(true), At(false), At(rgi)]; label=Altim.rgi2label[Altim.rginum2txt[rgi]], markersize=15, color=cmap[i])
    end

    for rgi = 98
        scatter!(grace_dm_gt[At(rgi)], dm_gt[At(false), At(false), At(rgi)]; label=nothing, markersize=15, color=cmap[findfirst(plot_order_rgi .== rgi)], strokecolor=:black, strokewidth=3)

        scatter!(-20, -75; label=nothing, markersize=15, color=(cmap[findfirst(plot_order_rgi .== rgi)], 0.0), strokecolor=:black, strokewidth=3)

        txt = "including endorheic"
        text!(-17, -75, text=txt, align=(:left, :center))
        
        #text!(grace_dm_gt[At(rgi)], dm_gt[At(false), At(true), At(rgi)], text=txt, align=(:left, :top))
    end

    xlims!(ax, -80, 20)
    ylims!(ax, -80, 20)
    lines!(ax, [-80, 20], [-80, 20], color=:black, linestyle=:dash)

    rmse = sqrt(mean((grace_dm_gt[At(setdiff(plot_order_rgi, [5, 19]))] .- dm_gt[At(true), At(false), At(setdiff(plot_order_rgi, [5, 19]))]) .^ 2))

    text!(-75, 15, text="RMSE = $(round(rmse, digits=1)) Gt/yr", align=(:left, :top))

    leg = Legend(f[1, 2], ax)
    display(f)

    save(output_dir, f)
    
    println("RMSE = $(round(rmse, digits=1)) Gt/yr")

    rmse_incl_endorheic = sqrt(mean((grace_dm_gt[At(setdiff(plot_order_rgi, [5, 19]))] .- dm_gt[At(false), At(false), At(setdiff(plot_order_rgi, [5, 19]))]) .^ 2))
    println("RMSE including endorheic = $(round(rmse_incl_endorheic, digits=1)) Gt/yr")

    rmse_static_fac = sqrt(mean((grace_dm_gt[At(setdiff(plot_order_rgi, [5, 19]))] .- dm_gt[At(true), At(true), At(setdiff(plot_order_rgi, [5, 19]))]) .^ 2))
    println("RMSE including static fac = $(round(rmse_static_fac, digits=1)) Gt/yr")
end





#### THIS COMPUTES THE RMSE BETWEEN GRACE AND ALL MODELLED REGIONAL TRENDS
# ---> NOTE: THIS SHOULD BE USED THE dates4trend(at top of script) ARE THE SAME AS THE GRACE DATES <---
f = Figure()
ax = Axis(f[1, 1])

nonendorheic_scale = collect(1 .- glacier_dm_gt[At(setdiff(plot_order_rgi, [5, 19])), At(true)] ./ glacier_dm_gt[At(setdiff(plot_order_rgi, [5, 19])), At(false)])

grace_atlim_rmse = zeros(drun)
for run0 in drun
    x = collect(grace_dm_gt[At(setdiff(plot_order_rgi, [5, 19]))]);
    y = collect(regions0_fit["dm_altim"][At(run0), At(setdiff(plot_order_rgi, [5, 19])), At("trend")]) .* nonendorheic_scale
    scatter!(x, y, label=run0)

    delta = x .- y;
    grace_atlim_rmse[At(run0)] = sqrt(mean(delta .^ 2))
    println("RMSE = $(round(grace_atlim_rmse[At(run0)], digits=1)) Gt/yr: $(run0)")
    
end
f

min_rmse = minimum(grace_atlim_rmse)
println("minimum RMSE with GRACE = $(round(min_rmse, digits=1)) Gt/yr for run $(drun[findfirst(minimum(grace_atlim_rmse) .== grace_atlim_rmse)])")


