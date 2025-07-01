"""
    discharge_rgi(path2discharge, path2rgi_regions; fractional_error = 0.15)

Aggregate glacier discharge data by RGI regions and calculate uncertainty.

# Arguments
- `path2discharge`: Path to file containing discharge data
- `path2rgi_regions`: Path to file containing RGI region geometries
- `fractional_error`: Fractional error to apply to discharge values (default: 0.15)

# Returns
- Array with dimensions [region, error] containing discharge values and uncertainties
  for each RGI region, including special aggregations for HMA (98) and global (99)
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

Aggregate glacier model run data by RGI (Randolph Glacier Inventory) regions.

This function loads data from multiple model runs, aggregates variables by RGI regions,
and creates special aggregations for HMA (regions 13-15) and global (regions 1-19).

# Arguments
- `path2runs`: Vector of paths to model run files
- `vars2sum`: List of variable names to aggregate (default includes common glacier variables)

# Returns
- Dictionary containing aggregated data for each variable, with dimensions [run, rgi, date]
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
    rgi_trends(regional_sum::AbstractDict, discharge_rgi, daterange)

Calculate trends, accelerations, amplitudes, and phases for regional glacier variables.

This function fits time series models to various glacier variables across different regions
and model runs, computing key parameters that characterize their temporal evolution.

# Arguments
- `regional_sum`: Dictionary containing time series data for different glacier variables
- `discharge_rgi`: Discharge values by RGI region, used for GSI calculations
- `daterange`: Time period over which to calculate trends

# Returns
- A DimensionalArray with dimensions for variable name, model run, RGI region, and parameter
  (trend, acceleration, amplitude, phase)

The function handles special calculations for derived variables:
- `net_acc`: Calculated as accumulation minus evaporation/condensation
- `gsi`: Glacier sustainability index, calculated as (runoff + discharge) / net_acc
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

            x = decimalyear.(ddate)
            x = x .- mean(x)

            
            for rgi in drgi
                if varname == "net_acc"
                    y = regional_sum["acc"][At(binned_synthesized_file), At(rgi), minimum(daterange)..maximum(daterange)] .- regional_sum["ec"][At(binned_synthesized_file), At(rgi), minimum(daterange)..maximum(daterange)]
                else
                    y = regional_sum[varname][At(binned_synthesized_file), At(rgi), minimum(daterange)..maximum(daterange)]
                end

                # first fit linear trend and seasonal [this creates less spread in linear fit as acceleration fixed to zero]
                dm_fit = curve_fit(offset_trend_seasonal2, x, y, p3)
                region_fit[At(varname), At(binned_synthesized_file), At(rgi), At("trend")] = dm_fit.param[2]

                # then fit linear trend and acceleration and seasonal
                dm_fit = curve_fit(offset_trend_acceleration_seasonal2, x, y, p3)
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
"""
    rgi_trends(da::AbstractDimArray, daterange)

Calculate trend, acceleration, amplitude, and phase parameters for each RGI region in a time series.

This function fits temporal models to glacier data for each RGI region within the specified date range.
It performs two curve fits: first a linear trend with seasonal component, then a model with
linear trend, acceleration, and seasonal components.

# Arguments
- `da`: A DimensionalData array with dimensions for RGI region and date
- `daterange`: Date range to analyze (can be specified as a range of DateTime objects)

# Returns
- A DimensionalData array with dimensions for RGI region and parameter, where parameters include
  trend, acceleration, amplitude, and phase of the fitted curves.

# Note
The function first centers the decimal years around their mean to improve numerical stability.
"""
function rgi_trends(da::AbstractDimArray, daterange)

    dparameter = Dim{:parameter}(["trend", "acceleration", "amplitude", "phase"])
    drgi = dims(da, :rgi)

    ddate = dims(da[:,minimum(daterange)..maximum(daterange)], :date)

    d = decimalyear.(ddate)
    d = d .- mean(d)

    region_fit = zeros(drgi, dparameter)

    for rgi in drgi

        # first fit linear trend and seasonal [this creates less spread in linear fit as acceleration fixed to zero]
        dm_fit = curve_fit(offset_trend_seasonal2, d, da[At(rgi), minimum(daterange)..maximum(daterange)], p3)
        region_fit[At(rgi), At("trend")] = dm_fit.param[2]

        # then fit linear trend and acceleration and seasonal
        dm_fit = curve_fit(offset_trend_acceleration_seasonal2, d, da[At(rgi), minimum(daterange)..maximum(daterange)], p3)
        region_fit[At(rgi), At("acceleration")] = dm_fit.param[3]
        region_fit[At(rgi), At("amplitude")] = hypot(dm_fit.param[4], dm_fit.param[5])
        region_fit[At(rgi), At("phase")] = 365.25 * (mod(0.25 - atan(dm_fit.param[5], dm_fit.param[4]) / (2π), 1))
    end

    return region_fit
end

"""
    region_fit_ref_and_err(region_fit, path2reference; p = 0.95, discharge = nothing)

Calculate reference values and error estimates for regional fit parameters.

This function extracts reference values from a specified path and computes error estimates
based on the quantile of absolute differences between all values and the reference value.

# Arguments
- `region_fit`: A DimensionalData array with dimensions for variable name, path, RGI region, and parameter
- `path2reference`: The reference path to use for extracting reference values
- `p`: Quantile level for error estimation (default: 0.95)
- `discharge`: Discharge values with error estimates, required when "net_acc" is present in variables

# Returns
- A DimensionalData array with dimensions for variable name, RGI region, parameter, and error,
  where the error dimension contains the reference values (false) and error estimates (true)

# Note
When "net_acc" is present, the function also calculates error propagation for the glacier state index (gsi)
using standard error propagation formulas.
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


"""
    runs_center!(runs_rgi::Dict, daterange)

Center glacier model run data by subtracting the temporal mean for each variable.

# Arguments
- `runs_rgi`: Dictionary containing time series data for different glacier variables
- `daterange`: Time period over which to calculate the mean

# Returns
- The modified `runs_rgi` dictionary with centered data

This function modifies the input dictionary in-place, subtracting the temporal mean
of each variable over the specified date range from all values of that variable.
"""
function runs_center!(runs_rgi::Dict, daterange)
 
    for k in keys(runs_rgi)
        mean_rgi = dropdims(mean(runs_rgi[k][:,:,minimum(daterange)..maximum(daterange)], dims = :date), dims = :date)
        runs_rgi[k] .-= mean_rgi
    end

    return runs_rgi
end

"""
    regions_center!(runs_rgi, daterange)

Center regional glacier data by subtracting the temporal mean.

# Arguments
- `runs_rgi`: DimensionalArray containing regional glacier data
- `daterange`: Time period over which to calculate the mean

# Returns
- The modified `runs_rgi` array with centered data

This function modifies the input array in-place, subtracting the temporal mean
over the specified date range from all values in the array.
"""
function regions_center!(runs_rgi, daterange)
    mean_rgi = dropdims(mean(runs_rgi[:, minimum(daterange)..maximum(daterange), At(false)], dims=:date), dims=:date)
    runs_rgi .-= mean_rgi

    return runs_rgi
end


"""
    runs_delta!(runs_rgi::Dict, path2reference) -> Dict

Calculate differences between each run and a reference run for all variables.

# Arguments
- `runs_rgi`: Dictionary containing time series data for different glacier variables
- `path2reference`: Path to the reference run to subtract from all other runs

# Returns
- The modified `runs_rgi` dictionary with delta values

This function modifies the input dictionary in-place, subtracting the reference run
from all runs for each variable, creating anomalies relative to the reference.
"""
function runs_delta!(runs_rgi, path2reference)
    for k in keys(runs_rgi)
        runs_rgi[k] = @d runs_rgi[k] .- runs_rgi[k][At(path2reference), :, :]
    end
    return runs_rgi
end
"""
    runs_quantile(runs_rgi, p; on_abs=true) -> Dict

Calculate quantiles of glacier data across multiple model runs.

# Arguments
- `runs_rgi`: Dictionary containing DimensionalArrays with model run data
- `p`: Quantile to calculate (e.g., 0.95 for 95th percentile)
- `on_abs`: Whether to calculate quantiles on absolute values (default: true)

# Returns
- Dictionary with same keys as input, containing quantile values for each region and date

This function computes quantiles across the run dimension for each variable, region, and date.
When `on_abs=true`, quantiles are calculated on absolute values, useful for error estimation.
"""
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

"""
    runs_select(runs_rgi, runs2select) -> Dict

Extract a subset of model runs from a dictionary of glacier data.

# Arguments
- `runs_rgi`: Dictionary containing DimensionalArrays with model run data
- `runs2select`: Vector of run identifiers to select

# Returns
- Dictionary with same keys as input, containing only the selected runs

This function creates a new dictionary with the same structure as the input,
but including only the specified model runs.
"""
function runs_select(runs_rgi, runs2select)
    runs_selected = Dict()
    for k in keys(runs_rgi)
        runs_selected[k] = runs_rgi[k][At(runs2select), :, :]
    end
    return runs_selected
end


"""
    runs_ref_and_err(runs_rgi, path2reference; p = 0.95) -> Dict

Create a dictionary containing reference values and error estimates for glacier model runs.

# Arguments
- `runs_rgi`: Dictionary containing DimensionalArrays with model run data
- `path2reference`: Path identifier for the reference run
- `p`: Quantile level for error estimation (default: 0.95)

# Returns
- Dictionary with same keys as input, containing DimensionalArrays with an additional
  error dimension where false=reference values and true=error estimates

This function calculates deviations from a reference run, computes quantile-based
error estimates, and combines them into a single data structure with an error dimension.
"""
function runs_ref_and_err(runs_rgi, path2reference; p = 0.95)
    # remove reference run
    runs_rgi_delta = runs_delta!(deepcopy(runs_rgi), path2reference)

    # get 95% confidence interval
    regions_err = runs_quantile(runs_rgi_delta, p; on_abs=true)

    # selected run as reference
    regions_ref = runs_select(runs_rgi, path2reference)

    regions = Dict()
    derror = Dim{:error}([false, true])
    for k in keys(regions_ref)
        regions[k] = zeros(dims(regions_ref[k])..., derror)
        regions[k][:,:,At(false)] = regions_ref[k]
        regions[k][:,:,At(true)] = regions_err[k]
    end

    return regions
end


"""
    rgi_endorheic(path2river_flux, glacier_summary_file; dates4trend=nothing)

Calculate glacier mass change and runoff trends for endorheic and non-endorheic basins by RGI region.

# Arguments
- `path2river_flux`: Path to file containing river flux data with terminal status information
- `glacier_summary_file`: Path to NetCDF file containing glacier runoff and mass change data
- `dates4trend`: Optional date range over which to calculate trends

# Returns
- A DimensionalArray with dimensions for variable name (dm, runoff), RGI region, and endorheic correction
  where false=all glaciers and true=endorheic glaciers only

This function identifies endorheic (landlocked) glaciers by matching glacier IDs to terminal river basins,
calculates mass change and runoff trends, and aggregates results by RGI region including special
aggregations for HMA (region 98) and global (region 99).
"""
function rgi_endorheic(path2river_flux, glacier_summary_file; dates4trend=nothing)
    ds = NCDataset(glacier_summary_file)
    glaciers = DataFrame(rgiid = ds["rgiid"][:], runoff = eachrow(ds["runoff"][:,:]), dm = eachrow(ds["dm"][:,:]))
    DataFrames.colmetadata!(glaciers,"runoff", "date", ds["Ti"][:]; style=:note)
    DataFrames.colmetadata!(glaciers,"dm", "date", ds["Ti"][:]; style=:note)

    
    drgi = Dim{:rgi}([collect(1:19)..., 98, 99])

    # load river flux data
    river_flux = GeoDataFrames.read(path2river_flux)

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

    indices = [searchsortedfirst(gid1, x) for x in glaciers[:, "rgiid"]]

    # some glaciers are not in the river_flux file (no rivers in the Antarctic)
    valid_indices = indices .<= length(gid1)

    # add endorheic column to glaciers
    glaciers[!, :endorheic] .= false
    glaciers[valid_indices, :endorheic] .= endorheic[indices[valid_indices]]

    # compute mass change for each region for endorheic and glaciers

    # now compute rates for all glaciers
    varnames = ["dm", "runoff"]  # unit of Gt

    # fit trend to all glaciers
    df_tsfit!(glaciers, varnames; progress=true, datelimits=dates4trend)

    
    # add rgi column
    glaciers[!, :rgi] = [parse(Int16, g[7:8]) for g in glaciers.rgiid]

    # (dh - fac)/1000 * area * volume2mass = dm_gt
    glaciers[!, :dm] = glaciers.dm_trend
    glaciers[!, :runoff] = glaciers.runoff_trend

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
"""
    error_bar_table(region_fits, varnames, params, rgis; digits=2)

Create a formatted table of regional glacier parameter values with error estimates.

# Arguments
- `region_fits`: DimensionalArray containing fitted parameters with dimensions for variable name, 
                 RGI region, parameter type, and error
- `varnames`: List of variable names to include in the table
- `params`: List of parameters to include (e.g., "trend", "acceleration", "amplitude", "phase")
- `rgis`: List of RGI region IDs to include in the table
- `digits`: Number of decimal places for rounding numeric values (default: 2)

# Returns
- DataFrame with columns for RGI ID, region name, and formatted parameter values with error estimates
  in the format "value ± error" for each variable and parameter combination

The function automatically determines appropriate units for each parameter type and includes
them in the column headers.
"""
function error_bar_table(region_fits, varnames, params, rgis; digits=2)

    rgi_labels = rginum2label.(rgis)

    regional_results = DataFrames.DataFrame(rgi = rgis, region_name = rgi_labels)

    for varname in varnames
        #println("    $(varname)")
        for param in params
            if param in ["acceleration"]
                unit = "Gt/yr^2"
            elseif param in ["trend"]
                unit = "Gt/yr"
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
        #println("$(rginum2label(rgi))")
        for varname in varnames
            #println("    $(varname)")
            for param in params
                v = region_fits[At(varname), At(rgi), At(param), At(false)]
                err = region_fits[At(varname), At(rgi), At(param), At(true)]
                

                if param in ["acceleration"]
                    v = v
                    err = err
                    unit = "Gt/yr^2"
                elseif param in ["trend"]
                    unit = "Gt/yr"
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


"""
    geotile2dimarray_kgm2(binned_synthesized_dv_file; vars2extract=["runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze", "dm"], trim=true, reference_period=nothing)

Convert geotile data to a dimensional array with standardized units (kg/m²).

# Arguments
- `binned_synthesized_dv_file`: Path to the file containing geotile data
- `vars2extract`: List of variable names to extract from the file (default includes common glacier variables)
- `trim`: Whether to trim data to valid range (default: true)
- `reference_period`: Optional time period to use as reference for anomaly calculation

# Returns
- A DimensionalArray with dimensions for variable name, geotile ID, and time, with units in kg/m²
  (except for 'dv' which is in meters)

This function loads geotile data, converts variables to consistent units (kg/m²), and optionally
calculates anomalies relative to a reference period. It handles special unit conversions for
volume ('dv') and mass ('dm') variables.
"""
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
        valid, = validrange(.!isnan.(geotiles0[1, vars2extract[1]]))
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


"""
    geotiles_mean_error(path2ref, path2files; p=0.95, reference_period=nothing)

Calculate reference values and error estimates for geotile data across multiple model runs.

This function computes the deviation of each model run from a reference run, then calculates
error estimates based on the quantile of absolute differences.

# Arguments
- `path2ref`: Path to the reference model run file
- `path2files`: Vector of paths to model run files to compare against the reference
- `p`: Quantile level for error estimation (default: 0.95)
- `reference_period`: Optional time period for centering data (subtracting temporal mean)

# Returns
- A DimensionalArray with dimensions for variable name, geotile, time, and error,
  where the error dimension contains reference values (false) and error estimates (true)
"""
function geotiles_mean_error(path2ref, path2files; p=0.95, reference_period=nothing)

    geotiles_ref = geotile2dimarray_kgm2(path2ref; reference_period)
    geotiles0 = geotile2dimarray_kgm2.(path2files; reference_period)

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


"""
    geotiles_mean_error_glaciers(glaciers, geotiles; show_progress=true, vars2downscale=nothing)

Downscale geotile-level data to individual glaciers, preserving reference values and error estimates.

# Arguments
- `glaciers`: DataFrame containing glacier data with columns for rgiid, geotile, and area_km2
- `geotiles`: DimensionalArray with dimensions for variable name, geotile, time, and error
- `show_progress`: Whether to display a progress bar during computation (default: true)
- `vars2downscale`: Optional subset of variables to downscale (default: all variables in geotiles)

# Returns
- A DimensionalArray with dimensions for variable name, glacier RGI ID, time, and error,
  where values are scaled by glacier area and the error dimension contains reference values (false)
  and error estimates (true)
"""
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


"""
    _simrun_init(nsamples, single_geotile_test)

Initialize a dictionary of DimensionalArrays for simulation runs.

# Arguments
- `nsamples`: Number of simulation runs
- `single_geotile_test`: Geotile identifier for testing

# Returns
- Dictionary with mission keys containing DimensionalArrays filled with NaN values
  with dimensions for run, geotile, date, and height
"""
function _simrun_init(;nsamples, missions2include, single_geotile_test)
    height_range, height_center = project_height_bins()
    date_range, date_center = project_date_bins()
    dheight = Dim{:height}(height_center)
    ddate = Dim{:date}(date_center)
    drun = Dim{:run}(1:nsamples)
    dgeotile = Dim{:geotile}([single_geotile_test])
    dmission = Dim{:mission}(keys(project_products()))

    dh = Dict{String,DimArray}()
    for mission in vcat(missions2include, "Synthesis")
        dh[mission] = fill(NaN, drun, dgeotile, ddate, dheight)
    end

    return dh
end

"""
    _simrun2areaaverage(dh; surface_mask, geotile_width=2, geotile2extract=nothing, sigma_error=2)

Calculate area-averaged height anomalies from simulation runs.

# Arguments
- `dh`: Dictionary of DimensionalArrays containing height data for each mission
- `surface_mask`: Surface mask type for area calculations
- `geotile_width`: Width of geotiles in degrees (default: 2)
- `geotile2extract`: Specific geotile to extract (default: nothing)
- `sigma_error`: Multiplier for error calculation (default: 2)

# Returns
- Tuple of (median_values, error_estimates) dictionaries containing area-averaged statistics
"""
function _simrun2areaaverage(dh; surface_mask, geotile_width=2, geotile2extract=nothing, sigma_error=2)
    # Calculate area-averaged height anomalies
    dh0_area_average = Dict()
    area_km2 = _geotile_area_km2(surface_mask, geotile_width)

    # Extract dimensions and convert dates to decimal years
    (drun, dgeotile, ddate, dheight) = dims(dh[first(keys(dh))])
    decyear = decimalyear.(collect(ddate))
    ddate = Dim{:date}(decyear)

    # Extract specific geotile data for each mission
    for mission in keys(dh)
        dh[mission] = DimArray(dh[mission][geotile=At(geotile2extract)][:, :, :], (drun, ddate, dheight))
    end

    # Initialize area-averaged arrays
    for mission in keys(dh)
        dh0_area_average[mission] = fill(NaN, drun, ddate)
    end

    # Calculate area-averaged values for each mission
    for mission in keys(dh)
        valid_run, valid_date, valid_height = validrange(.!isnan.(dh[mission]))
        for i in drun
            dh0_area_average[mission][At(i), valid_date] = dh_area_average(dh[mission][At(i), valid_date, valid_height], area_km2[At(geotile2extract), valid_height])
        end
    end

    # Calculate error estimates
    dh0_area_average_error = Dict()
    for mission in keys(dh)
        valid_run, valid_date = validrange(.!isnan.(dh0_area_average[mission]))
        dh0_area_average_error[mission] = dropdims(std(dh0_area_average[mission][:, valid_date], dims=:run), dims=:run) * sigma_error
    end

    # Calculate median values
    dh0_area_average_median = Dict()
    for mission in keys(dh)
        valid_run, valid_date = validrange(.!isnan.(dh0_area_average[mission]))
        dh0_area_average_median[mission] = dropdims(median(dh0_area_average[mission][:, valid_date], dims=:run), dims=:run)
    end

    return dh0_area_average_median, dh0_area_average_error
end
