
"""
    discharge_rgi(path2discharge, path2rgi_regions)

Calculate glacier discharge for each RGI (Randolph Glacier Inventory) region.

# Arguments
- `path2discharge`: Path to the glacier discharge data file (.jld2)
- `path2rgi_regions`: Path to the RGI regions shapefile

# Returns
- A DimensionalArray containing discharge values (in Gt/yr) for each RGI region,
  including aggregated values for High Mountain Asia (RGI ID 98) and global (RGI ID 99)

This function loads glacier discharge data, assigns each discharge point to its 
corresponding RGI region using spatial operations, and sums the discharge values
for each region. It also calculates aggregated values for High Mountain Asia 
(regions 13-15) and global (regions 1-19).
"""
function discharge_rgi(path2discharge, path2rgi_regions; fractional_error = 0.15)

    drgi = Dim{:rgi}([1:19; 98; 99])
    derror = Dim{:error}([false, true])
    
    discharge = load(path2discharge, "discharge")

    # determine wich region discharge is within
    rgi_regions = GeoDataFrames.read(path2rgi_regions)

    discharge[!, :geometry] = tuple.(discharge.longitude, discharge.latitude)

    discharge = FlexiJoins.innerjoin(
        (discharge, rgi_regions),
        by_pred(:geometry, GO.within, :geometry)
    )

    discharge = rename(discharge, "RGI_CODE" => "rgi" )

    # sum discharge for each region
    discharge_rgi = zeros(drgi, derror; name = "frontal ablation [Gt yr⁻¹]")
    for rgi in drgi
        index = discharge.rgi .== rgi
        discharge_rgi[At(rgi), At(false)] = sum(discharge[index, "discharge_gtyr"])

    end

    # add HMA and global
    discharge_rgi[At(98), At(false)] = sum(discharge_rgi[13..15, At(false)])
    
    # add global
    discharge_rgi[At(99), At(false)] = sum(discharge_rgi[1..19, At(false)])

    # add error (simple fraction of total)
    discharge_rgi[:, At(true)] = discharge_rgi[:, At(false)] * fractional_error
    return discharge_rgi
end



"""
    runs2rgi(path2runs, vars2sum=["dv_altim", "runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze", "dv", "dm", "dm_altim"])

Read in individual model runs, sum to regional level, and add HMA and global aggregations.

# Arguments
- `path2runs`: Array of paths to model run files (.jld2)
- `vars2sum`: Array of variable names to sum across regions (default includes common glaciological variables)

# Returns
- Dictionary where keys are variable names and values are arrays with dimensions [run, region, date]

This function processes multiple model run files, aggregates data by RGI (Randolph Glacier Inventory) 
regions, and adds special aggregations for High Mountain Asia (RGI ID 98) and global (RGI ID 99).
"""
function runs2rgi(path2runs, vars2sum=["dv_altim", "runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze", "dv", "dm", "dm_altim"])

    # read in example file
    drgi = Dim{:rgi}([1:19; 98; 99])
   
    # drun needs full path as some runs are distiguished by folder not just filename
    drun = Dim{:run}(path2runs)

    regional_sum = Dict() # this needs to be a Dict since variables have different date dimensions

    example_data = FileIO.load(replace(path2runs[1], ".jld2" => "_gembfit_dv.jld2"), "geotiles")
    for varname in vars2sum
        ddate = Dim{:date}(colmetadata(example_data, varname, "date"))
        regional_sum[varname] = fill(NaN, drun, drgi, ddate)
    end

    Threads.@threads for binned_synthesized_file in path2runs
    #for binned_synthesized_file in path2runs
        #binned_synthesized_file = path2reference

        binned_synthesized_dv_file = replace(binned_synthesized_file, ".jld2" => "_gembfit_dv.jld2")
        geotiles0 = FileIO.load(binned_synthesized_dv_file, "geotiles")

        # group by rgi and sum
        geotiles_reg = groupby(geotiles0, :rgiid)

        # sum across timeseries
        regions0 = DataFrames.combine(geotiles_reg, vars2sum .=> Ref ∘ sum; renamecols=false)

        # add a HMA region 
        index = (regions0[:, :rgiid] .<= 15) .& (regions0[:, :rgiid] .>= 13)
        foo = DataFrames.combine(regions0[index, :], vcat("rgiid", vars2sum) .=> Ref ∘ sum; renamecols=false)
        foo.rgiid .= 98
        regions0 = append!(regions0, foo)

        # add a Global region 
        index = (regions0[:, :rgiid] .<= 19) .& (regions0[:, :rgiid] .>= 1)
        foo = DataFrames.combine(regions0[index, :], vcat("rgiid", vars2sum) .=> Ref ∘ sum; renamecols=false)
        foo.rgiid .= 99
        regions0 = append!(regions0, foo)

        for varname in vars2sum
            # reduce all region time series vectors to a matrix
            regional_sum[varname][At(binned_synthesized_file), At(regions0.rgiid), :] = reduce(hcat, regions0[!, varname])'
        end
    end

    return regional_sum
end



"""
    rgi_trends(regional_sum, discharge_rgi, daterange)

Calculate trends, accelerations, amplitudes, and phases for regional glacier mass change variables.

This function fits time series models to regional glacier data and computes various parameters
including linear trends, accelerations, seasonal amplitudes, and phases. It also calculates
derived variables like equivalent runoff and runoff excess percentage.

# Arguments
- `regional_sum`: Dictionary containing regional time series data for different variables
- `discharge_rgi`: Discharge data by RGI (Randolph Glacier Inventory) region
- `daterange`: Date range to use for trend fitting

# Returns
- A DimensionalData array with dimensions for variable names, model runs, RGI regions, and fit parameters
  (trend, acceleration, amplitude, phase)

# Notes
- First fits a linear trend with seasonal cycle, then fits a model with acceleration included
- Calculates equivalent runoff as accumulation minus evaporation/condensation minus discharge
- Runoff excess percentage is calculated as ((runoff - runoff_eq) / runoff_eq) * 100
"""
function rgi_trends(regional_sum::AbstractDict, discharge_rgi, daterange)
    varnames = vcat(collect(keys(regional_sum)), ["net_acc", "gsi"])
    dparameter = Dim{:parameter}(["trend", "acceleration", "amplitude", "phase"])
    dvarname = Dim{:varname}(varnames)
    drun = dims(regional_sum[varnames[1]], :run)
    drgi = dims(regional_sum[varnames[1]], :rgi)
    region_fit = zeros(dvarname, drun, drgi, dparameter)

    Threads.@threads for binned_synthesized_file in drun

        for varname in setdiff(varnames, ["gsi"])
           
            if varname == "net_acc"
                ddate = dims(regional_sum["runoff"][:,:,minimum(daterange)..maximum(daterange)], :date)
            else
                ddate = dims(regional_sum[varname][:,:,minimum(daterange)..maximum(daterange)], :date)
            end

            x = Altim.decimalyear.(ddate)
            x = x .- mean(x)

            
            for rgi in drgi
                if varname == "net_acc"
                    y = regional_sum["acc"][At(binned_synthesized_file), At(rgi), minimum(daterange)..maximum(daterange)] .- regional_sum["ec"][At(binned_synthesized_file), At(rgi), minimum(daterange)..maximum(daterange)]
                else
                    y = regional_sum[varname][At(binned_synthesized_file), At(rgi), minimum(daterange)..maximum(daterange)]
                end

                # first fit linear trend and seasonal [this creates less spread in linear fit as acceleration fixed to zero]
                dm_fit = curve_fit(Altim.offset_trend_seasonal2, x, y, Altim.p3)
                region_fit[At(varname), At(binned_synthesized_file), At(rgi), At("trend")] = dm_fit.param[2]

                # then fit linear trend and acceleration and seasonal
                dm_fit = curve_fit(Altim.offset_trend_acceleration_seasonal2, x, y, Altim.p3)
                region_fit[At(varname), At(binned_synthesized_file), At(rgi), At("acceleration")] = dm_fit.param[3]
                region_fit[At(varname), At(binned_synthesized_file), At(rgi), At("amplitude")] = hypot(dm_fit.param[4], dm_fit.param[5])
                region_fit[At(varname), At(binned_synthesized_file), At(rgi), At("phase")] = 365.25 * (mod(0.25 - atan(dm_fit.param[5], dm_fit.param[4]) / (2π), 1))
            end
        end
    end

    varname = "gsi"
    for rgi in drgi
        for param in ["trend", "acceleration"]

            if param == "trend"
                D = discharge_rgi[At(rgi), At(false)];
            else
                D = 0
            end

            region_fit[At(varname), :, At(rgi), At(param)] = ((region_fit[At("runoff"), :, At(rgi), At(param)] .+ D) ./ region_fit[At("net_acc"), :, At(rgi), At(param)])
        end
    end

    return region_fit
end


function rgi_trends(da::AbstractDimArray, daterange)

    dparameter = Dim{:parameter}(["trend", "acceleration", "amplitude", "phase"])
    drgi = dims(da, :rgi)

    ddate = dims(da[:,minimum(daterange)..maximum(daterange)], :date)

    d = Altim.decimalyear.(ddate)
    d = d .- mean(d)

    region_fit = zeros(drgi, dparameter)

    for rgi in drgi

        # first fit linear trend and seasonal [this creates less spread in linear fit as acceleration fixed to zero]
        dm_fit = curve_fit(Altim.offset_trend_seasonal2, d, da[At(rgi), minimum(daterange)..maximum(daterange)], Altim.p3)
        region_fit[At(rgi), At("trend")] = dm_fit.param[2]

        # then fit linear trend and acceleration and seasonal
        dm_fit = curve_fit(Altim.offset_trend_acceleration_seasonal2, d, da[At(rgi), minimum(daterange)..maximum(daterange)], Altim.p3)
        region_fit[At(rgi), At("acceleration")] = dm_fit.param[3]
        region_fit[At(rgi), At("amplitude")] = hypot(dm_fit.param[4], dm_fit.param[5])
        region_fit[At(rgi), At("phase")] = 365.25 * (mod(0.25 - atan(dm_fit.param[5], dm_fit.param[4]) / (2π), 1))
    end

    return region_fit
end



"""
    region_fit_ref_and_err(region_fit, path2reference; p = 0.95)

Extract reference run values and calculate error estimates for regional fit parameters.

This function takes a multi-dimensional array of regional fit parameters from multiple model runs
and extracts the values from a reference run while calculating error estimates based on the
spread across all runs.

# Arguments
- `region_fit`: A DimensionalData array with dimensions for variable name, model run, RGI region, and parameter.
- `path2reference`: The identifier for the reference model run to extract values from.
- `p`: Probability level for quantile-based error estimation (default: 0.95 for 95% confidence).

# Returns
- A DimensionalData array with dimensions for variable name, RGI region, parameter, and error,
  where the error dimension contains the reference value (false) and the error estimate (true).

# Note
Error estimates are calculated as the p-th quantile of the absolute deviations from the reference value.
"""
function region_fit_ref_and_err(region_fit, path2reference; p = 0.95, discharge = nothing)
    derror = Dim{:error}([false, true])
    (dvarname, _, drgi, dparameter) = dims(region_fit)
    regions_fit = fill(NaN, dvarname, drgi, dparameter, derror)
    for varname in dvarname
        for rgi in drgi
            for param in dparameter
                v = region_fit[At(varname), At(path2reference), At(rgi), At(param)]
                regions_fit[At(varname), At(rgi), At(param), At(false)] = v

                err0 = quantile(abs.(region_fit[At(varname), :, At(rgi), At(param)] .- v), p)
                regions_fit[At(varname), At(rgi), At(param), At(true)] = err0
            end
        end
    end

    if "net_acc" in dvarname
        if isnothing(discharge)
            error("gsi is an exisiting variable, so discharge_error must be provided as a kwarg")
        end

        A = regions_fit[At("runoff"), :, At("trend"), At(false)] .+ discharge[:, At(false)]
        σA = sqrt.(regions_fit[At("runoff"), :, At("trend"), At(true)].^ 2 .+ discharge[:, At(true)] .^ 2)

        B = regions_fit[At("net_acc"), :, At("trend"), At(false)]
        σB = regions_fit[At("net_acc"), :, At("trend"), At(true)]

        err = abs.(regions_fit[At("gsi"), :, At("trend"), At(true)]) .* sqrt.((σA ./ A) .^ 2 .+ (σB ./ B) .^ 2)
        regions_fit[At("gsi"), :, At("trend"), At(true)] = err
    end

    return regions_fit
end


function runs_center!(runs_rgi::Dict, daterange)
 
    for k in keys(runs_rgi)
    #k = "dm_altim"
        mean_rgi = dropdims(mean(runs_rgi[k][:,:,minimum(daterange)..maximum(daterange)], dims = :date), dims = :date)
        runs_rgi[k] .-= mean_rgi
    end

    return runs_rgi
end

function regions_center!(runs_rgi, daterange)

    #k = "dm_altim"
    mean_rgi = dropdims(mean(runs_rgi[:, minimum(daterange)..maximum(daterange), At(false)], dims=:date), dims=:date)
    runs_rgi .-= mean_rgi

    return runs_rgi
end


function runs_delta!(runs_rgi, path2reference)
    for k in keys(runs_rgi)
        runs_rgi[k] = @d runs_rgi[k] .- runs_rgi[k][At(path2reference), :, :]
    end
    return runs_rgi
end


function runs_quantile(runs_rgi, p; on_abs=true)
    runs_quantile = Dict()
    on_abs ? (fun=abs) : (fun = identity)

    for k in keys(runs_rgi)
        (drun, drgi, ddate) = dims(runs_rgi[k])
        runs_quantile[k] = fill(NaN, drgi, ddate)
        for rgi in drgi
            for date in ddate
                v0 = @view runs_rgi[k][:, At(rgi), At(date)]
                if any(isnan.(v0))
                    continue
                else
                    runs_quantile[k][At(rgi), At(date)] = quantile(fun.(v0), p)
                end
            end
        end
    end
    return runs_quantile
end


function runs_select(runs_rgi, runs2select)
    runs_selected = Dict()
    for k in keys(runs_rgi)
        runs_selected[k] = runs_rgi[k][At(runs2select), :, :]
    end
    return runs_selected
end


function runs_ref_and_err(runs_rgi, path2reference; p = 0.95)
    # remove reference run
    runs_rgi_delta = Altim.runs_delta!(deepcopy(runs_rgi), path2reference)

    # get 95% confidence interval
    regions_err = Altim.runs_quantile(runs_rgi_delta, p; on_abs=true)

    # selected run as reference
    regions_ref = Altim.runs_select(runs_rgi, path2reference)

    regions = Dict()
    derror = Dim{:error}([false, true])
    for k in keys(regions_ref)
        regions[k] = zeros(dims(regions_ref[k])..., derror)
        regions[k][:,:,At(false)] = regions_ref[k]
        regions[k][:,:,At(true)] = regions_err[k]
    end

    return regions
end


function rgi_endorheic(path2river_flux, path2perglacier; dates4trend = nothing)

    # We need to map between river sink and glacier source
    volume2mass = Altim.δice / 1000
    
    drgi = Dim{:rgi}([collect(1:19)..., 98, 99])

    # load river flux data
    river_flux = GeoDataFrames.read(path2river_flux)

    # load in glacier outlines
    glaciers = load(path2perglacier, "glaciers")

    # select only terminal rivers
    river_flux = river_flux[river_flux[:, :NextDownID].==0, :]

    # match glaciers to river sink type (endorheic or ocean terminating)
    gid1 = reduce(vcat, river_flux.RGIId)
    endorheic = falses(length(gid1))

    c = 1
    for r in eachrow(river_flux)
        n = length(r.RGIId)
        endorheic[c:c+n-1] .= !(r.ocean_terminating)
        c += n
    end

    p1 = sortperm(gid1)
    gid1 = gid1[p1]
    endorheic = endorheic[p1]

    indices = [searchsortedfirst(gid1, x) for x in glaciers[:, "RGIId"]]

    # some glaciers are not in the river_flux file (no rivers in the Antarctic)
    valid_indices = indices .<= length(gid1)

    # add endorheic column to glaciers
    glaciers[!, :endorheic] .= false
    glaciers[valid_indices, :endorheic] .= endorheic[indices[valid_indices]]

    # compute mass change for each region for endorheic and glaciers

    # now compute dh rate for all glaciers
    varname = "dh"  # unit of [m i.e.]
    date = collect(dims(glaciers[1, varname], :date))

    # find the index of valid data
    index = .!isnan.(glaciers[1, varname])
    DataFrames.colmetadata!(glaciers, varname, "date", date; style=:note)

    # fit dh trend to all glaciers
    Altim.df_tsfit!(glaciers, [varname]; progress=true, datelimits=dates4trend)

    # fit fac trend to all glaciers
    varname = "fac"
    date = collect(dims(glaciers[1, varname], :date))
    DataFrames.colmetadata!(glaciers, varname, "date", date; style=:note)
    Altim.df_tsfit!(glaciers, [varname]; progress=true, datelimits=dates4trend)

    # fit runoff trend to all glaciers
    varname = "runoff"
    date = collect(dims(glaciers[1, varname], :date))
    DataFrames.colmetadata!(glaciers, varname, "date", date; style=:note)
    Altim.df_tsfit!(glaciers, [varname]; progress=true, datelimits=dates4trend)

    # add rgi column
    glaciers[!, :rgi] = [parse(Int16, g[7:8]) for g in glaciers.RGIId]

    # (dh - fac)/1000 * area * volume2mass = dm_gt
    glaciers[!, :dm] = (glaciers.dh_trend .- glaciers.fac_trend) / 1000 .* sum.(glaciers.area_km2) * volume2mass
    glaciers[!, :runoff] = (glaciers.runoff_trend) / 1000 .* sum.(glaciers.area_km2)

    # derive statistics for each region
    dvarnames = Dim{:varname}(["dm", "runoff"])
    dendorheic_only = Dim{:endorheic_correction}([false, true])
    glacier_dm_gt = zeros(dvarnames, drgi, dendorheic_only; name="all vs. endorheic only", metadata=Dict("units" => "Gt/yr", "daterange" => dates4trend))

    for rgi in setdiff(collect(drgi), [98, 99])
        index = glaciers.rgi .== rgi
        endorheic = glaciers[:, :endorheic] .& index
        for varname in dvarnames
            glacier_dm_gt[At(varname), At(rgi), At(false)] = sum(glaciers[index, varname])
            glacier_dm_gt[At(varname), At(rgi), At(true)] = sum(glaciers[endorheic, varname])
        end
    end

    # add HMA and global
    glacier_dm_gt[:, At(98), :] = sum(glacier_dm_gt[:, 13..15, :], dims=:rgi)

    # add global
    glacier_dm_gt[:, At(99), :] = sum(glacier_dm_gt[:, 1..19, :], dims=:rgi)

    return glacier_dm_gt
end


function error_bar_table(region_fits, varnames, params, rgis; digits=2)

    rgi_labels = Altim.rginum2label.(rgis)

    regional_results = DataFrames.DataFrame(rgi = rgis, region_name = rgi_labels)

    for varname in varnames
        #println("    $(varname)")
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

    for (i, rgi) in enumerate(rgis)
        #println("$(Altim.rginum2label(rgi))")
        for varname in varnames
            #println("    $(varname)")
            for param in params
                v = region_fits[At(varname), At(rgi), At(param), At(false)]
                err = region_fits[At(varname), At(rgi), At(param), At(true)]
                

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
               
                if param in ["phase"]
                    v = round(Int, v)
                    err = round(Int, err)
                else
                    v = round(v; digits)
                    err = round(err; digits)
                end
                
                column_label = Symbol("$(varname)_$(param)_[$unit]")
                regional_results[i, column_label] = "$(v) ± $(err)"
            end
        end
    end
    return regional_results
end


function geotile2dimarray_kgm2(
    binned_synthesized_dv_file; 
    vars2extract=["runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze", "dm"], 
    trim=true, 
    reference_period=nothing
    )

    geotiles0 = load(binned_synthesized_dv_file, "geotiles")

    geotiles0[!, :area_km2] = sum.(geotiles0.area_km2)

    # trim data to valid range and update date metadata
    if trim
        valid, = Altim.validrange(.!isnan.(geotiles0[1, vars2extract[1]]))
    else
        valid = 1:length(geotiles0[1, vars2extract[1]])
    end

    ddate = colmetadata(geotiles0, vars2extract[1], "date")[valid]
    dTi = Ti(ddate)

    varunits = Dict()
    for varname in vars2extract
        if varname == "dv"
            varunits[varname] = u"m"
        else
            varunits[varname] = u"kg/m^2"
        end
    end
    geotile_area = Dict()
    for (id, area) in zip(geotiles0.id, geotiles0.area_km2)
        geotile_area[id] = area * u"km^2"
    end

    dvarname = Dim{:varname}(vars2extract)
    dgeotile = Dim{:geotile}(geotiles0.id)

    # TODO: I've hardcoded the units here to be the units of the first variable in vars2extract ... this will error for mixed units... fixing this would require splitting geotile_out into multiple arrays. 
    geotile_out = fill(NaN * varunits[first(vars2extract)], dvarname, dgeotile, dTi; name=joinpath(splitpath(binned_synthesized_dv_file)[end-2:end]), metadata=Dict("area" => geotile_area))

    # remove unneeded columns
    geotiles0 = geotiles0[:, vcat("geometry", "area_km2", vars2extract)]

    for varname in vars2extract

        if varname == "dv"
            # convert from km3 to m
            geotile_out[At(varname), :, :] = reduce(hcat, getindex.(geotiles0[:, varname], (valid,)) ./ geotiles0[:, :area_km2] * 1000 * 1000 * varunits[varname])'

        elseif varname == "dm"
            # convert from Gt to kg/m²
            geotile_out[At(varname), :, :] = reduce(hcat, getindex.(geotiles0[:, varname], (valid,)) ./ geotiles0[:, :area_km2] * 1000 * 1000 * varunits[varname])'
        else
            # convert from km3(assumed ice density of 910 kg/m³) to kg/m²
            geotile_out[At(varname), :, :] = reduce(hcat, getindex.(geotiles0[:, varname], (valid,)) ./ geotiles0[:, :area_km2] * 910 * 1000 * varunits[varname])'
        end
    end

    if !isnothing(reference_period)
        geotile_out .-= mean(geotile_out[Ti=reference_period[1]..reference_period[2]], dims=Ti)
    end

    return geotile_out
end


function geotiles_mean_error(path2ref, path2files; p=0.95, reference_period=nothing)

    geotiles_ref = Altim.geotile2dimarray_kgm2(path2ref; reference_period)
    geotiles0 = Altim.geotile2dimarray_kgm2.(path2files; reference_period)

    # calculate the deviation relative to a reference run
    for gt in geotiles0
        gt .-= geotiles_ref
    end

    geotiles0 = cat(geotiles0..., dims=Dim{:run}(DimensionalData.name.(geotiles0)))

    # calculate the 95% quantile of the absolute deviation relative to the reference run
    geotiles_ref_err = copy(geotiles_ref)
    (dvarname, dgeotile, dTi) = dims(geotiles_ref_err)
    Threads.@threads for varname in dvarname
        for geotile in dgeotile
            for Ti in dTi
                geotiles_ref_err[At(varname), At(geotile), At(Ti)] = quantile(abs.(geotiles0[At(varname), At(geotile), At(Ti), :]), p)
            end
        end
    end

    metadata0 = DimensionalData.metadata(geotiles_ref)
    geotiles_out = cat(geotiles_ref, geotiles_ref_err; dims=Dim{:error}([false, true]))

    return geotiles_out
end


function geotiles_mean_error_glaciers(glaciers, geotiles; show_progress=true, vars2downscale=nothing)

    drgiid = Dim{:rgiid}(glaciers.rgiid)
    derror = Dim{:error}([false, true])
    (dvarname, dgeotile, dTi, derror) = dims(geotiles)
    
    if !isnothing(vars2downscale)
        dvarname = Dim{:varname}(vars2downscale)
    end

    glacier_out = DimArray(zeros(length(dvarname), length(drgiid), length(dTi), length(derror)) * 1u"Gt", (dvarname, drgiid, dTi, derror))

    prog = Progress(nrow(glaciers)*length(dvarname); desc="Downscaling from geotiles to glaciers", enabled=show_progress)

    for varname in dvarname
        Threads.@threads for r in eachrow(glaciers)
            glacier_out[At(varname), At(r.rgiid), :, At(false)] = geotiles[At(varname), At(r.geotile), :, At(false)] .* r.area_km2

            glacier_out[At(varname), At(r.rgiid), :, At(true)] = geotiles[At(varname), At(r.geotile), :, At(true)] .* r.area_km2

            next!(prog)
        end
      
    end
    finish!(prog)

    return glacier_out
end