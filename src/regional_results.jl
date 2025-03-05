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
    surface_mask=["glacier", "glacier_rgi7"]
    dem_id=["best", "cop30_v2"]
    curvature_correct=[false, true]
    amplitude_correct=[true]
    binning_method = ["median" "meanmadnorm10" "meanmadnorm5" "meanmadnorm3"]
    paramater_set=[1, 2, 3, 4]
    binned_folder=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg")

    paths = Altim.pathlocal
    reference_run = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_dh_best_cc_meanmadnorm3_v01_filled_ac_p1_synthesized.jld2";
    perglacier_reference_run = replace(reference_run, ".jld2" => "_perglacier.jld2")
    globaldischarge_fn = joinpath(paths.data_dir, "GlacierOutlines/GlacierDischarge/global_glacier_discharge.jld2")

    rivers_path = joinpath(paths.data_dir, "rivers/MERIT_Hydro_v07_Basins_v01")
    glacier_rivers_path = joinpath(rivers_path, "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.arrow")

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


# aggrigate individual geotile data to regions, for each model run and all variables
# This section takes 45s:
# 1. Loads and aggregates geotile data into regional time series for each model run
# 2. Calculates regional uncertainties from ensemble spread
# 3. Incorporates GRACE satellite data for comparison
# 4. Reformats data for plotting, including:
#    - Interpolating to common time grid
#    - Adding High Mountain Asia (HMA) and global regions
#    - Calculating error bounds
#    - Aligning different data sources
# Key outputs:
# - regions0: Dictionary containing raw regional time series for each variable/run
# - regions: DataFrame with processed regional time series and metadata
# - df: Final DataFrame formatted for plotting, containing:
#   - Regional time series from altimetry, GRACE and GEMB
#   - Error bounds
#   - Common time grid and units

begin
    # read in example file
    regions = FileIO.load(replace(reference_run, ".jld2" => "_gembfit_dv.jld2"), "geotiles")
    varnames = setdiff(names(regions), ["id", "extent", "glacier_frac", "area_km2", "discharge", "landice_frac", "floating_frac", "geometry", "group", "pscale", "Δheight", "mie2cubickm", "rgiid"])

    # drun needs full path as some runs are distiguished by folder not just filename
    drun = Dim{:run}(path2runs)
    drgi = Dim{:rgi}([1:19; 98; 99])
    dparameter = Dim{:parameter}(["trend", "acceleration", "amplitude", "phase"])
    dvarname = Dim{:varname}(varnames)
    regions0 = Dict() # this needs to be a Dict since variables have different date dimensions
    regions0_fit = Dict()
    
    for varname in varnames
        ddate = Dim{:date}(colmetadata(regions, varname, "date"))
        regions0[varname] = fill(NaN, drun, drgi, ddate)
        regions0_fit[varname] = fill(NaN, drun, drgi,dparameter)
    end

    for binned_synthesized_file in path2runs
    #binned_synthesized_file = path2runs[50]

        binned_synthesized_dv_file = replace(binned_synthesized_file, ".jld2" => "_gembfit_dv.jld2")
        geotiles0 = FileIO.load(binned_synthesized_dv_file, "geotiles")

        # group by rgi and sum
        geotiles_reg = groupby(geotiles0, :rgiid)

        # sum across timeseries
        regions = DataFrames.combine(geotiles_reg, varnames .=> Ref ∘ sum; renamecols=false)

        # add a HMA region 
        index = (regions[:, :rgiid] .<= 15) .& (regions[:, :rgiid] .>= 13)
        foo = DataFrames.combine(regions[index, :], vcat("rgiid", varnames) .=> Ref ∘ sum; renamecols=false)
        foo.rgiid .= 98
        regions = append!(regions, foo)

        # add a Global region 
        index = (regions[:, :rgiid] .<= 19) .& (regions[:, :rgiid] .>= 1)
        foo = DataFrames.combine(regions[index, :], vcat("rgiid", varnames) .=> Ref ∘ sum; renamecols=false)
        foo.rgiid .= 99
        regions = append!(regions, foo)

        for varname in varnames
        #varname = varnames[1]
            
            regions0[varname][At(binned_synthesized_file), At(regions.rgiid), :] = reduce(hcat, regions[!, varname])'

            d = Altim.decimalyear.(dims(regions0[varname], :date))
            d = d .- mean(d)
            for rgi in drgi
                valid = .!isnan.(regions0[varname][At(binned_synthesized_file), At(rgi), :])

                
                # first fit linear trend and seasonal [this creates less spread in linear fit as acceleration fixed to zero]
                dm_fit = curve_fit(Altim.offset_trend_seasonal2, d[valid], regions0[varname][At(binned_synthesized_file), At(rgi), :][valid], Altim.p3)
                regions0_fit[varname][At(binned_synthesized_file), At(rgi), At("trend")] = dm_fit.param[2]

                # then fit linear trend and acceleration and seasonal
                dm_fit = curve_fit(Altim.offset_trend_acceleration_seasonal2, d[valid], regions0[varname][At(binned_synthesized_file), At(rgi), :][valid], Altim.p3)
                regions0_fit[varname][At(binned_synthesized_file), At(rgi), At("acceleration")] = dm_fit.param[3]
                regions0_fit[varname][At(binned_synthesized_file), At(rgi), At("amplitude")] = hypot(dm_fit.param[4], dm_fit.param[5])
                regions0_fit[varname][At(binned_synthesized_file), At(rgi), At("phase")] = 365.25 * (mod(0.25 - atan(dm_fit.param[5], dm_fit.param[4]) / (2π), 1))
            end
        end
    end

    # create df for reference run
    binned_synthesized_file = reference_run
    binned_synthesized_dv_file = replace(binned_synthesized_file, ".jld2" => "_gembfit_dv.jld2")

    geotiles0 = FileIO.load(binned_synthesized_dv_file, "geotiles")

    # group by rgi and sum
    geotiles_reg = groupby(geotiles0, :rgiid)

    # sum across timeseries
    regions = DataFrames.combine(geotiles_reg, varnames .=> Ref ∘ sum; renamecols=false)

    # add an HMA region 
    index = (regions[:, :rgiid] .<= 15) .& (regions[:, :rgiid] .>= 13)
    foo = DataFrames.combine(regions[index, :], vcat("rgiid", varnames) .=> Ref ∘ sum; renamecols=false)
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
        #varname = "dv_altim"
        varname_err = "$(varname)_err"
        var0 = regions0[varname]

        # loop for each rgi and remove mean to center data
        for rgi in drgi
            var1 = var0[:, At(rgi), :]
            rowrange, colrange = Altim.validrange(.!isnan.(var1))

            delta = collect(mean(var1[:, colrange], dims=2))
            var0[:, At(rgi), :] .-= delta * ones(1, size(var1, 2))
        end

        # sanity check
        #=
        region = 1
        p = lines(var0[1,At(region),:])
        for i = 2:size(var0, 1)
            lines!(var0[i,At(region),:])
        end 
        p
        =#

        # santiy check
        #CairoMakie.heatmap(regions0[varname][:,9,:])

        # remove reference run from ensemble for error calculation
        var1 = deepcopy(var0)
        for run in drun
            var1[At(run), :, :] .-= var0[At(reference_run), :, :]
        end

        # sanity check
        #=
        p = lines(var1[1,At(region),:])
        for i = 2:size(var0, 1)
            lines!(var1[i,At(region),:])
        end 
        p
        =#

        ddate = dims(var1, :date)
        err = var1[1,:,:]
        err .= NaN
        for rgi in drgi
            for date in ddate
                foo = abs.(var1[:, At(rgi), At(date)]);
                valid = .!isnan.(foo)
                if any(valid)

                    ## THIS IS WHERE THE ERROR IS DEFINED [symetric 95% confidence interval from reference run]
                    err[At(rgi), At(date)] = quantile(abs.(var1[valid, At(rgi), At(date)]), 0.95)
                    ##
                end
            end
        end
  
        regions[!, varname_err] = eachrow(collect(err))
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
    regions_fit = fill(NaN, dvarname, drgi, dparameter, derror)
    for varname in varnames
        for rgi in drgi
            for param in dparameter
                v = regions0_fit[varname][At(binned_synthesized_file), At(rgi), At(param)]
                regions_fit[At(varname), At(rgi), At(param), At(false)] = v

                err = quantile(abs.(regions0_fit[varname][:, At(rgi), At(param)] .- v), 0.95)
                regions_fit[At(varname), At(rgi), At(param), At(true)] = err
            end
        end
    end

    # reformat for plotting
    varnames = setdiff(names(regions), ["rgiid"])
    dates = DateTime(2000, 1, 15):Month(1):DateTime(2023, 1, 15)
    dates_decyear = Altim.decimalyear.(dates)
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
            model = LinearInterpolation(r[varname], date0)
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
    metadata!(df, "date", dates; style=:note)

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


# This section processes glacier discharge and mass change data by region:
#
# 1. Loads and processes glacier discharge data:
#    - Loads discharge data from file
#    - Maps discharge points to RGI regions using spatial join
#    - Aggregates discharge by RGI region
#    - Adds HMA and global totals
#
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
# - glacier_rgi: DataFrame with mass change by RGI region
# - dm_gt: Array of mass change trends under different scenarios

begin
    # load in discharge data     
    discharge = load(globaldischarge_fn, "discharge")

    # determine wich region discharge is within
    rgi_regions = GeoDataFrames.read(paths.rgi6_regions_shp)

    discharge[!, :geometry] = tuple.(discharge.longitude, discharge.latitude)

    discharge = FlexiJoins.innerjoin(
        (discharge, rgi_regions),
        by_pred(:geometry, GO.within, :geometry)
    )
    discharge = rename(discharge, "RGI_CODE" => "rgi" )


    # derive statistics for each region
    discharge_rgi = DataFrame()
    discharge_rgi[!, :rgi] = 1:19
    discharge_rgi[!, :discharge_gtyr] = fill(NaN, 19)

    for g in eachrow(discharge_rgi)
        index = discharge.rgi .== g.rgi
        g.discharge_gtyr = sum(discharge[index, "discharge_gtyr"])
    end

    # add HMA and global
    push!(discharge_rgi, deepcopy(discharge_rgi[1,:]))
    discharge_rgi[end,:rgi] = 98
    index = (discharge_rgi.rgi .>= 13) .& (discharge_rgi.rgi .<= 15)
    discharge_rgi[end,Not(:rgi)] = sum.(eachcol(discharge_rgi[index, Not(:rgi)]))

    push!(discharge_rgi, deepcopy(discharge_rgi[1,:]))
    discharge_rgi[end,:rgi] = 99
    index = discharge_rgi.rgi .<= 19
    discharge_rgi[end,Not(:rgi)] = sum.(eachcol(discharge_rgi[index, Not(:rgi)]))


    # load in glacier outlines
    glaciers = load(perglacier_reference_run, "glaciers")

    # load in river flux data
    river_flux = GeoDataFrames.read(glacier_rivers_path)

    # select only terminating rivers
    river_flux = river_flux[river_flux[:, :NextDownID] .== 0, :]

    # add endorheic column to glaciers
    glaciers[!, :endorheic] .= false

    ## THIS SECTION IS VERY VERY SLOW
    for g in eachrow(glaciers)
        g.RGIId
        for r in eachrow(river_flux)
            if any(r.RGIId .== g.RGIId)
                g.endorheic = !r.ocean_terminating
                break
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
    glacier_rgi = DataFrame()
    glacier_rgi[!, :rgi] = 1:19
    glacier_rgi[!, :dm_gt] = fill(NaN, 19)
    glacier_rgi[!, :endorheic_dm_gt] = fill(NaN, 19)

    for g in eachrow(glacier_rgi)
        index = glaciers.rgi .== g.rgi
        endorheic = glaciers[:, :endorheic] .& index;
        g.dm_gt = sum(glaciers[index, "dm_gt"])
        g.endorheic_dm_gt = sum(glaciers[endorheic, "dm_gt"]) 
    end

    # add HMA and global
    push!(glacier_rgi, deepcopy(glacier_rgi[1,:]))
    glacier_rgi[end,:rgi] = 98
    index = (glacier_rgi.rgi .>= 13) .& (glacier_rgi.rgi .<= 15)
    glacier_rgi[end,Not(:rgi)] = sum.(eachcol(glacier_rgi[index, Not(:rgi)]))

    push!(glacier_rgi, deepcopy(glacier_rgi[1,:]))
    glacier_rgi[end,:rgi] = 99
    index = glacier_rgi.rgi .<= 19
    glacier_rgi[end,Not(:rgi)] = sum.(eachcol(glacier_rgi[index, Not(:rgi)]))

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
        index2 = findfirst((df.mission .== "grace") .& (df.rgi .== rgi) .& (df.var .== varname))
        if !isnothing(index2)
            grace_dm_gt[At(rgi)] = df[index2, :mid_trend]
        end

        for static_fac in dstatic_fac
            for endorheic_correction in dendorheic_correction
                
                scale = 1.0
                if endorheic_correction
                    idx = findfirst(glacier_rgi.rgi .== rgi)
                    scale = (glacier_rgi[idx, :dm_gt] - glacier_rgi[idx, :endorheic_dm_gt]) / glacier_rgi[idx, :dm_gt]
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


# fit trend and acceleration to regional data
# this approach assumes errors are normally distributed and uncorrelated which is not 
# correct... it's more accurate to calculate trends and accelerations directly on the 
# ensemble members ... then take the mean and std of the trends and accelerations

function tsfit_regions(regions; daterange = (nothing, nothing), acceleration = false)
    if isnothing(daterange[1])
        sdate = DateTime(0)
    else
        sdate = daterange[1]
    end

    if isnothing(daterange[2])
        edate = DateTime(2500)
    else
        edate = daterange[2]
    end

    acceleration && (m = (x, p) -> p[1] .+ p[2] * x .+ p[3] * x.^2)
    acceleration || (m = (x, p) -> p[1] .+ p[2] * x)

    df = DataFrame(rgiid = regions[:, :rgiid])
    for varname in names(regions)[.!occursin.("err", names(regions)).&.!occursin.("rgiid", names(regions))]
        #varname = first(varnames)
        df[!, "$(varname)_trend"] = fill(NaN, nrow(df))
        df[!, "$(varname)_trend_err"] = fill(NaN, nrow(df))
        acceleration && (df[!, "$(varname)_acc"] = fill(NaN, nrow(df)))
        acceleration && (df[!, "$(varname)_acc_err"] = fill(NaN, nrow(df)))

        for i in 1:nrow(regions)
            #i = 1
            x = colmetadata(regions, varname, "date")
            y = regions[i,varname]
            index = (x .>= sdate) .& (x .<= edate) .& .!(isnan.(y))

            if !any(index)
                continue
            end

            y_err = regions[i, "$(varname)_err"][index]
            x = x[index]
            y = y[index]

            # normalize x and Y
            x = Altim.decimalyear.(x)

            # normalize
            x = x .- mean(x)
            y = y .- mean(y)

            wt = 1 ./ (y_err .^ 2)
            acceleration && (p0 = [0.5, 0.5, 0.5])
            acceleration || (p0 = [0.5, 0.5])
            # p: model parameters
            fit = curve_fit(m, x, y, wt, p0)

            cf = coef(fit)
            ci = confint(fit; level = 0.95)    # 5% significance level

            df[i, "$(varname)_trend"] = cf[2]
            df[i, "$(varname)_trend_err"] = abs(ci[2][1] - cf[2])
            acceleration && (df[i, "$(varname)_acc"] = cf[3])
            acceleration && (df[i, "$(varname)_acc_err"] = abs(ci[3][1] - cf[3]))
        end
    end
    return df
end

# NOTE: This is not used as error are not properly propogated 
regions_fit_lsq = tsfit_regions(regions; daterange=(Date(2003), Date(2023)), acceleration=true)



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

# make table with regional results
regional = DataFrame()
rgi2include = [1:4; 6:12; 16:18; 98]

Altim.rgi2label[Altim.rginum2txt[rgi]]







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
f


n = 6;
index = falses(nrow(glacier_rgi))
large_region_glaciers = reverse(plot_order_rgi[(end-(n-1)):end]);
for rgi in large_region_glaciers
    index = index .| (glacier_rgi.rgi .== rgi)
end

large_region_dm_gt = sum(glacier_rgi.dm_gt[index]) / glacier_rgi.dm_gt[findfirst(glacier_rgi.rgi .== 99)]
println("$(n) largest glaciers comprise $(round(large_region_dm_gt, digits=2)) of total glacier loss. Regions: $(large_region_glaciers)")



# plot all iterations for a single region
var0 = "dv_altim"
rgi = 98

p = lines(regions0[var0][1, At(rgi), :])
for i = 2:length(regions0[var0][:, At(rgi), 1])
    lines!(regions0[var0][i, At(rgi), :])
end
p


exclude_regions = [13, 14, 15]
varnames = ["acc", "runoff", "dm_altim", "dm", "ec", "dv_altim", "fac"]
varnames = [ "dm_altim", "dm"]
params = ["trend", "acceleration", "amplitude", "phase"]

begin
    for rgi in drgi
        if rgi in exclude_regions
            continue
        end
        println("$(Altim.rgi2label[Altim.rginum2txt[rgi]])")
        for varname in varnames
            println("    $(varname)")
            for param in params
                v = round(Int,regions_fit[At(varname), At(rgi), At(param), At(false)])
                err = round(Int,regions_fit[At(varname), At(rgi), At(param), At(true)])

                println("        $(param): $(v) ± $(err)")
            end
        end

        # estimate fractional increase in runoff
        discharge0 = discharge_rgi[findfirst(discharge_rgi.rgi .== rgi), :discharge_gtyr]
        runoff = round(Int, regions_fit[At("runoff"), At(rgi), At("trend"), At(false)])
        acc = round(Int, regions_fit[At("acc"), At(rgi), At("trend"), At(false)])
        ec = round(Int, regions_fit[At("ec"), At(rgi), At("trend"), At(false)])
        increase_in_runoff = runoff - (acc - discharge0 - ec)

        println()
        println("    total discharge = $(round(Int,discharge0)) Gt/yr")
        println("    increase in runoff = $(round(Int,increase_in_runoff)) Gt/yr")
        println("    equilibrium runoff = $(round(Int,(acc - discharge0 - ec))) Gt/yr")
        frac_increase_in_runoff = increase_in_runoff / (acc - discharge0 - ec)
        println("    fractional increase in runoff = $(round(frac_increase_in_runoff, digits=2))")
    end
end


# fraction of glacier loss that is endorheic
endorheic_frac = glacier_rgi[findfirst(glacier_rgi.rgi .== 99), :endorheic_dm_gt] / glacier_rgi[findfirst(glacier_rgi.rgi .== 99), :dm_gt];
println("fraction of glacier loss that flows to ocean = $(1-round(endorheic_frac, digits=3))")

hma_frac_of_total_endorheic = glacier_rgi[findfirst(glacier_rgi.rgi .== 98), :endorheic_dm_gt] / sum(glacier_rgi[glacier_rgi.rgi.<=19, :endorheic_dm_gt]);
println("HMA comprises $(round(hma_frac_of_total_endorheic, digits=3)) of total endorheic loss")

amp_static_fac = round(Int, regions_fit[At("dv"), At(99), At("amplitude"), At(false)]) * 0.85;
amp_dynamic_fac = round(Int, regions_fit[At("dm"), At(99), At("amplitude"), At(false)]);
println("Global amplitude would be $(round((amp_static_fac - amp_dynamic_fac) ./ amp_dynamic_fac, digits=2)) larger if static fac was applied")

# GRACE rgi regions
rgis = collect([1:4; 6:12; 16:18; 98])

function resample(colorscheme, n)
    newscheme = ColorSchemes.ColorScheme(
        get(ColorSchemes.colorschemes[colorscheme], (0:n) / n))
    return newscheme
end

colormap = :Dark2_4
cmap = resample(colormap, length(plot_order_rgi))
# set_theme!(Theme(palette=(color=resample(colormap, length(rgis)),),))


begin
    f = Figure()
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
        scatter!(grace_dm_gt[At(rgi)], dm_gt[At(true), At(true), At(rgi)]; label=Altim.rgi2label[Altim.rginum2txt[rgi]], markersize=15, color=cmap[i])
    end

    for rgi = 98
        scatter!(grace_dm_gt[At(rgi)], dm_gt[At(false), At(true), At(rgi)]; label=nothing, markersize=15, color=cmap[findfirst(plot_order_rgi .== rgi)], strokecolor=:black, strokewidth=3)

        scatter!(-17, -75; label=nothing, markersize=15, color=(cmap[findfirst(drgi .== rgi)], 0.), strokecolor=:black, strokewidth=3)

        txt = "including endorheic"
        text!(-14, -75, text=txt, align=(:left, :center))
        
        #text!(grace_dm_gt[At(rgi)], dm_gt[At(false), At(true), At(rgi)], text=txt, align=(:left, :top))
    end

    xlims!(ax, -80, 20)
    ylims!(ax, -80, 20)
    lines!(ax, [-80, 20], [-80, 20], color=:black, linestyle=:dash)

    rmse = sqrt(mean((grace_dm_gt .- dm_gt[At(true), At(true), :]) .^ 2))

    text!(-75, 15, text="RMSE = $(round(rmse, digits=1)) Gt/yr", align=(:left, :top))

    leg = Legend(f[1, 2], ax)
    display(f)
end


