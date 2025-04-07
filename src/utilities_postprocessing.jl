
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
function discharge_rgi(path2discharge, path2rgi_regions)

    drgi = Dim{:rgi}([1:19; 98; 99])
    
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
    discharge_rgi = zeros(drgi)
    for rgi in drgi
        index = discharge.rgi .== rgi
        discharge_rgi[At(rgi)] = sum(discharge[index, "discharge_gtyr"])
    end

    # add HMA and global
    discharge_rgi[At(98),Not(:rgi)] = sum(discharge_rgi[13..15])
    
    # add global
    discharge_rgi[At(99),Not(:rgi)] = sum(discharge_rgi[1..19])

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
        #binned_synthesized_file = reference_run

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



function rgi_trends(regional_sum, discharge_rgi, daterange)
    varnames = vcat(collect(keys(regional_sum)), ["runoff_eq", "runoff_excess_perc"])
    dparameter = Dim{:parameter}(["trend", "acceleration", "amplitude", "phase"])
    dvarname = Dim{:varname}(varnames)
    drun = dims(regional_sum[varnames[1]], :run)
    drgi = dims(regional_sum[varnames[1]], :rgi)
    region_fit = zeros(dvarname, drun, drgi, dparameter)

    Threads.@threads for binned_synthesized_file in drun

        for varname in setdiff(varnames, ["runoff_eq", "runoff_excess_perc"])
            #varname = varnames[1]
            ddate = dims(regional_sum[varname], :date)
            date_index = (ddate .>= minimum(daterange)) .&& (ddate .<= maximum(daterange))

            d = Altim.decimalyear.(dims(regional_sum[varname], :date))
            d = d .- mean(d)
            for rgi in drgi

                valid = .!isnan.(regional_sum[varname][At(binned_synthesized_file), At(rgi), :]) .& date_index

                # first fit linear trend and seasonal [this creates less spread in linear fit as acceleration fixed to zero]
                dm_fit = curve_fit(Altim.offset_trend_seasonal2, d[valid], regional_sum[varname][At(binned_synthesized_file), At(rgi), valid], Altim.p3)
                region_fit[At(varname), At(binned_synthesized_file), At(rgi), At("trend")] = dm_fit.param[2]

                # then fit linear trend and acceleration and seasonal
                dm_fit = curve_fit(Altim.offset_trend_acceleration_seasonal2, d[valid], regional_sum[varname][At(binned_synthesized_file), At(rgi), :][valid], Altim.p3)
                region_fit[At(varname), At(binned_synthesized_file), At(rgi), At("acceleration")] = dm_fit.param[3]
                region_fit[At(varname), At(binned_synthesized_file), At(rgi), At("amplitude")] = hypot(dm_fit.param[4], dm_fit.param[5])
                region_fit[At(varname), At(binned_synthesized_file), At(rgi), At("phase")] = 365.25 * (mod(0.25 - atan(dm_fit.param[5], dm_fit.param[4]) / (2π), 1))
            end
        end
    end

    ## THERE MIGHT BE A BUG HERE with runoff_eq acceleration
    for rgi in drgi
        for param in ["trend", "acceleration"]
            region_fit[At("runoff_eq"), :, At(rgi), At(param)] = region_fit[At("acc"), :, At(rgi), At(param)] .- region_fit[At("ec"), :, At(rgi), At(param)] .- discharge_rgi[At(rgi)]

            # same as (A-B)/B
            region_fit[At("runoff_excess_perc"), :, At(rgi), At(param)] = (region_fit[At("runoff"), :, At(rgi), At(param)] .- region_fit[At("runoff_eq"), :, At(rgi), At(param)]) ./ region_fit[At("runoff_eq"), :, At(rgi), At(param)] * 100
        end
    end

    return region_fit
end