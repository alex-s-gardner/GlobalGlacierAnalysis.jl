"""
    gemb_read(gemb_files; vars=[...])

Read GEMB (Glacier Energy and Mass Balance) data from multiple files.

# Arguments
- `gemb_files`: Array of file paths to GEMB data files
- `vars`: Array of variable names to read from the files. Default includes all common variables.

# Returns
A DataFrame containing the requested variables from all files concatenated together.
"""
function gemb_read(
    gemb_files;
    vars=[
        "latitude",
        "longitude",
        "date",
        "smb", "fac",
        "ec",
        "acc",
        "runoff",
        "melt",
        "fac_to_depth",
        "height_ref",
        "refreeze",
        "t_air",
        "rain"
    ]
)

    any(vars .== "latitude") && (lat = Vector{Vector{Float64}}())
    any(vars .== "longitude") && (lon = Vector{Vector{Float64}}())
    any(vars .== "date") && (t = Vector{Vector{Float64}}())
    any(vars .== "smb") && (smb = Vector{Vector{Float64}}())
    any(vars .== "fac") && (fac = Vector{Vector{Float64}}())
    any(vars .== "ec") && (ec = Vector{Vector{Float64}}())
    any(vars .== "acc") && (acc = Vector{Vector{Float64}}())
    any(vars .== "runoff") && (runoff = Vector{Vector{Float64}}())
    any(vars .== "melt") && (melt = Vector{Vector{Float64}}())
    any(vars .== "fac_to_depth") && (fac_to_depth = Vector{Vector{Float64}}())
    any(vars .== "height_ref") && (h = Vector{Vector{Float64}}())
    any(vars .== "refreeze") && (refreeze = Vector{Vector{Float64}}())
    any(vars .== "t_air") && (ta = Vector{Vector{Float64}}())
    any(vars .== "rain") && (rain = Vector{Vector{Float64}}())


    for fn in gemb_files

        matopen(fn) do file
            any(vars .== "latitude") && (lat0 = read(file, "lat"))
            any(vars .== "longitude") && (lon0 = read(file, "lon"))
            any(vars .== "date") && (t0 = read(file, "time"))
            any(vars .== "fac") && (fac0 = read(file, "FAC"))
            any(vars .== "smb") && (smb0 = read(file, "SMB"))
            any(vars .== "height_ref") && (h0 = read(file, "H"))
            any(vars .== "ec") && (ec0 = read(file, "EC"))
            any(vars .== "melt") && (melt0 = read(file, "Melt"))
            any(vars .== "fac_to_depth") && (fac_to_depth0 = read(file, "FACtoDepth"))
            any(vars .== "refreeze") && (refreeze0 = read(file, "Refreeze"))
            any(vars .== "t_air") && (ta0 = read(file, "Ta"))
            any(vars .== "rain") && (rain0 = read(file, "Rain"))
            any(vars .== "runoff") && (runoff0 = read(file, "Runoff"))
            any(vars .== "acc") && (acc0 = read(file, "Accumulation"))

            # this is going to cause issues at some point 
            m, n = size(smb0)

            any(vars .== "latitude") && (push!(lat, vec(repeat(lat0, 1, n)[:])))
            any(vars .== "longitude") && (push!(lon, vec(repeat(lon0, 1, n)[:])))
            any(vars .== "height_ref") && (push!(h, vec(repeat(h0, 1, n)[:])))
            any(vars .== "date") && (push!(t, vec(repeat(t0, m, 1)[:])))
            any(vars .== "fac") && (push!(fac, fac0[:]))
            any(vars .== "smb") && (push!(smb, smb0[:]))
            any(vars .== "ec") && (push!(ec, ec0[:]))
            any(vars .== "melt") && (push!(melt, melt0[:]))
            any(vars .== "fac_to_depth") && (push!(fac_to_depth, fac_to_depth0[:]))
            any(vars .== "refreeze") && (push!(refreeze, refreeze0[:]))
            any(vars .== "t_air") && (push!(ta, ta0[:]))
            any(vars .== "rain") && (push!(rain, rain0[:]))
            any(vars .== "runoff") && (push!(runoff, runoff0[:]))
            any(vars .== "acc") && (push!(acc, acc0[:]))
        end
    end
    gemb = DataFrame()

    any(vars .== "latitude") && (gemb[!, :latitude] = vcat(lat...))
    any(vars .== "longitude") && (gemb[!, :longitude] = vcat(lon...))
    any(vars .== "date") && (gemb[!, :date] = decimalyear2datetime.(vcat(t...)))
    any(vars .== "height_ref") && (gemb[!, :height_ref] = vcat(h...))
    any(vars .== "fac") && (gemb[!, :fac] = vcat(fac...))
    any(vars .== "smb") && (gemb[!, :smb] = vcat(smb...))
    any(vars .== "ec") && (gemb[!, :ec] = vcat(ec...))
    any(vars .== "melt") && (gemb[!, :melt] = vcat(melt...))
    any(vars .== "fac_to_depth") && (gemb[!, :fac_to_depth] = vcat(fac_to_depth...))
    any(vars .== "refreeze") && (gemb[!, :refreeze] = vcat(refreeze...))
    any(vars .== "t_air") && (gemb[!, :t_air] = vcat(ta...))
    any(vars .== "rain") && (gemb[!, :rain] = vcat(rain...))
    any(vars .== "runoff") && (gemb[!, :runoff] = vcat(runoff...))
    any(vars .== "acc") && (gemb[!, :acc] = vcat(acc...))

    return gemb
end


"""
    gemb_read2(gemb_file; vars="all", datebin_edges=nothing)

Read GEMB data from a single file with optional date binning.

# Arguments
- `gemb_file`: Path to a GEMB data file
- `vars`: Either "all" to read all variables, or an array of specific variable names
- `datebin_edges`: Optional array of date edges for binning the data

# Returns
A dictionary containing the requested variables with optional date binning applied.
"""
function gemb_read2(gemb_file; vars = "all", datebin_edges = nothing)
       
    vars0=["latitude","longitude", "date", "smb", "fac", "ec", "acc", "runoff", "melt", "fac_to_depth", "height", "refreeze", "t_air", "rain"]

    # variables that need to be converted from cumulitive outputs to rates

    var_cumulitive= ["smb", "ec", "acc", "runoff", "melt", "refreeze", "rain"]

    varmap= Dict(
        "latitude" => "lat", 
        "longitude" => "lon", 
        "date" => "time", 
        "smb" => "SMB", 
        "fac" => "FAC", 
        "ec" => "EC", 
        "acc" => "Accumulation", 
        "runoff" => "Runoff", 
        "melt" => "Melt", 
        "fac_to_depth" => "FACtoDepth", 
        "height" => "H", 
        "refreeze" => "Refreeze", 
        "t_air" => "Ta", 
        "rain" => "Rain"
    )

    gemb = Dict()
    if vars == "all"
        foo = matread(gemb_file)
        for v in varmap
            f = foo[v[2]];
            f[isnan.(f)] .= 0; #set any random nan to zero
            gemb[v[1]] = f
        end
    else
        if length(intersect(vars0, vars)) != length(vars)
            error("one or more vars not recognized")
        end

        gemb = matopen(gemb_file) do file
            foo = MAT.read.(Ref(file), [varmap[v] for v in vars])
            for (i, v) in enumerate(vars)
                f = foo[i];
                f[isnan.(f)] .= 0 #set any random nan to zero
                gemb[v] = f
            end
            return gemb
        end
    end

    # inputs are in cumulative amounts ... e.g. mm per 5 days
    # convert to cumulitive meters
    for k in keys(gemb)
        if any(var_cumulitive .== k)
            gemb[k] = cumsum(gemb[k] ./ 1000, dims=2)
        end
    end

    if .!isnothing(datebin_edges)
        ind = Vector{Union{UnitRange{Int64}, Nothing}}()
        t = vec(gemb["date"])
        for i = eachindex(datebin_edges)[1:end-1]
            a = findfirst((t .>= datebin_edges[i]) .& (t .< datebin_edges[i+1]))
            b = findlast((t .>= datebin_edges[i]) .& (t .< datebin_edges[i+1]))

            if isnothing(a)
                push!(ind, nothing)
            else
                push!(ind,a:b)
            end
        end

        for k in keys(gemb)
            if k == "date"
                continue
            end
            
            if size(gemb[k],2) == size(gemb["date"],2)
                foo = fill(NaN, (size(gemb[k],1), length(ind)))
                for i in eachindex(ind)
                    if !isnothing(ind[i])
                        foo[:, i] = mean(gemb[k][:, ind[i]], dims=2)
                    end
                end
                gemb[k] = foo;
            end
        end
        gemb["date"] = reshape((datebin_edges[1:end-1] .+ datebin_edges[2:end]) ./ 2, (1, length(datebin_edges) - 1))
    end
    
    return gemb
end


"""
    gemb2dataframe(; path2file)

Convert GEMB data from a JLD2 file to a DataFrame format.

# Arguments
- `path2file`: Path to a JLD2 file containing GEMB data

# Returns
A DataFrame with GEMB variables organized by RGI region, height adjustment, and precipitation scale.
"""
function gemb2dataframe(;
    path2file="/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0_geotile_filled_d_reg.jld2"
)

    gemb = FileIO.load(path2file)

    vars = setdiff(keys(gemb), ["area_km2", "nobs"])

    drgi = dims(gemb[first(vars)], :rgi)
    dΔheight = dims(gemb[first(vars)], :Δheight)
    dpscale = dims(gemb[first(vars)], :pscale)
    ddate = dims(gemb[first(vars)], :date)

    df = DataFrame()

    area = gemb["area_km2"]
    for var0 in vars
        #var0 = first(vars)

        v0 = gemb[var0]
        nobs0 = gemb["nobs"]
        for rgi in drgi
            for Δheight in dΔheight
                for pscale in dpscale
        
                    vout = vec(v0[At(rgi), :, At(pscale), At(Δheight)])
                    nobsout = vec(nobs0[At(rgi), :, At(pscale), At(Δheight)])

                    append!(df,
                        DataFrame(; rgi, var="$(var0)_km3", mission="gemb", Δheight, pscale, surface_mask="glacier", area_km2=area[At(rgi)], val=[vout], nobs=[nobsout], filename=path2file)
                    )
                end
            end
        end
    end

    DataFrames.metadata!(df, "date", collect(ddate); style=:note)
    return df
end

"""
    gemb_classes_densify!(df; n_densify=4)

Densify GEMB parameter classes by interpolating between existing parameter values.

# Arguments
- `df`: DataFrame containing GEMB data with pscale and Δheight parameters
- `n_densify`: Number of interpolated points to add between existing parameter values (default: 4)

# Returns
The input DataFrame with additional interpolated parameter combinations added.
"""
function gemb_classes_densify!(df; n_densify = 4)

    # interpolate between gemb classes
    classes_Δheight = sort(unique(df.Δheight))
    classes_pscale = sort(unique(df.pscale))
    rgis = unique(df.rgi)
    vars = unique(df.var)
    cols_interp = ["val", "nobs"]
    #n_densify = 4 # how many "between" runs to add

    df_interp = DataFrame()
    for var0 in vars
        index_var = df.var .== var0


        for rgi in rgis
            #rgi = first(rgis)
            index_rgi = df.rgi .== rgi
            for h = eachindex(classes_Δheight)
                h0 = df.Δheight .== classes_Δheight[h]

                for p = eachindex(classes_pscale)[1:end-1]
                    p1 = df.pscale .== classes_pscale[p]
                    p2 = df.pscale .== classes_pscale[p+1]

                    index = h0 .&& index_rgi .& index_var
                    v1 = df[index.&p1, :]
                    v2 = df[index.&p2, :]

                    if nrow(v2) != 1
                        error("nrow(v2) != 1")
                    end

                    foo = copy(v1)
                    for w in (1:n_densify) / (n_densify + 1)
                        for col in cols_interp
                            foo[:, col][1] = (w .* v1[:, col][1]) .+ ((1 - w) .* v2[:, col][1])
                            foo.pscale[1] = (w .* classes_pscale[p]) .+ ((1 - w) .* classes_pscale[p+1])
                            df_interp = append!(df_interp, foo)
                        end
                    end
                end
            end

            for p = eachindex(classes_pscale)
                p0 = df.pscale .== classes_pscale[p]

                for h = eachindex(classes_Δheight)[1:end-1]
                    h1 = df.Δheight .== classes_Δheight[h]
                    h2 = df.Δheight .== classes_Δheight[h+1]

                    index = p0 .&& index_rgi .& index_var
                    v1 = df[index.&h1, :]
                    v2 = df[index.&h2, :]

                    if nrow(v2) != 1
                        error("nrow(v2) != 1")
                    end

                    foo = copy(v1)
                    for w in (1:n_densify) / (n_densify + 1)
                        for col in cols_interp
                            foo[:, col][1] = (w .* v1[:, col][1]) .+ ((1 - w) .* v2[:, col][1])
                            foo.Δheight[1] = (w .* classes_Δheight[h]) .+ ((1 - w) .* classes_Δheight[h+1])
                            df_interp = append!(df_interp, foo)
                        end
                    end
                end
            end
        end
    end
    df = append!(df, df_interp)

    return df
end


"""
    _matrix_shift_ud!(matrix, shift; exclude_zeros_in_extrapolation=false, npts_linear_extrapolation=Inf)

Shifts a matrix up or down by a specified number of rows, with linear extrapolation at boundaries.

# Arguments
- `matrix`: Input matrix to be shifted in-place
- `shift`: Integer number of rows to shift. Positive shifts down, negative shifts up.
- `exclude_zeros_in_extrapolation`: If true, zero values are excluded when calculating extrapolation
- `npts_linear_extrapolation`: Maximum number of points to use for linear extrapolation

# Returns
The modified input matrix
"""
function _matrix_shift_ud!(matrix, shift; exclude_zeros_in_extrapolation=false, npts_linear_extrapolation=Inf)

    n = size(matrix, 2)
    h0 = (1.:n)./n .- 0.5
    M = hcat(ones(n),h0)

    # shift with linear extrapolation
    if shift > 0 # shift values down 
        # shift values down
        matrix[:, 1:(end-shift)] = matrix[:, (shift+1):end]

        # extrapolate upper boundary
        M0 = @view M[1:(end-shift), :]

        for i in axes(matrix, 1)
            y0 = @view matrix[i, 1:(end-shift)]
            if exclude_zeros_in_extrapolation
                valid = y0 .!= 0
                if any(valid)
                    if sum(valid) > npts_linear_extrapolation
                        valid = findall(valid)[(end-npts_linear_extrapolation+1):end]
                    end
                    param1 = M0[valid, :] \ y0[valid] 
                else
                    matrix[i, (end-shift+1):end] .= 0
                    continue
                end
            else
                param1 = M0 \ y0
            end
            
            matrix[i, (end-shift+1):end] = param1[1] .+ param1[2] .* h0[(end-shift+1):end]
        end

    elseif shift < 0 # shift values up
       
        # shift values up
        matrix[:, (1-shift):end] = matrix[:, 1:(end+shift)]

        # extrapolate lower boundary
        M0 = @view M[(1-shift):end, :]

        for i in axes(matrix, 1)
            y0 = @view matrix[i, (1-shift):end]
            if exclude_zeros_in_extrapolation
                valid = y0 .!= 0
                if any(valid)
                    if sum(valid) > npts_linear_extrapolation
                        valid = findall(valid)[1:npts_linear_extrapolation]
                    end
                    param1 = M0[valid, :] \ y0[valid]
                else
                    matrix[i, 1:(-shift)] .= 0
                    continue
                end
            else
                param1 = M0 \ y0
            end

            matrix[i, 1:(-shift)] = param1[1] .+ param1[2] .* h0[1:(-shift)]
        end
    end
    return matrix
end


"""
    gemb_rate_physical_constraints!(gemb)

Apply physical constraints to GEMB rate variables to ensure consistency.

# Arguments
- `gemb`: Dictionary containing GEMB variables

# Returns
The input dictionary with physically consistent values:
- Negative values for accumulation, rain, melt, and refreeze are set to zero
- Refreeze cannot exceed melt
- Runoff equals melt minus refreeze
- SMB equals accumulation minus runoff minus evaporation/condensation
"""
function gemb_rate_physical_constraints!(gemb)
    acc = gemb["acc"]
    acc[acc.<0] .= 0

    rain = gemb["rain"]
    rain[rain.<0] .= 0

    melt = gemb["melt"]
    melt[melt.<0] .= 0

    refreeze = gemb["refreeze"]
    refreeze[refreeze.<0] .= 0

    # refreeze can not exceed melt
    refreeze[refreeze.>melt] = melt[refreeze.>melt]

    # runoff must equal melt - refreeze [rain is not included in runoff calculation]
    gemb["runoff"] = melt .- refreeze;

    # SMB must equal acc - runoff - ec [rain is not included in SMB calculation]
    gemb["smb"] = acc .- gemb["runoff"] .- gemb["ec"];

    return gemb
end



"""
    gemb_bestfit_grouped(dv_altim, smb, fac, discharge, geotiles; single_geotile_test = nothing)

Calibrate the SMB model to grouped geotiles to handle glaciers that cross multiple geotiles.

# Arguments
- `dv_altim`: Volume change from altimetry data
- `smb`: Surface mass balance data
- `fac`: Firn air content data
- `discharge`: Discharge data
- `geotiles`: DataFrame containing geotile information
- `single_geotile_test`: Optional geotile ID to examine model fits for debugging 

# Returns
DataFrame with calibrated parameters for each geotile: id, extent, pscale, Δheight, and rmse
"""

function gemb_bestfit_grouped(dv_altim, smb, fac, discharge0, geotiles0; 
    distance_from_origin_penalty=2/100,
    seasonality_weight = 95/100, 
    single_geotile_test = nothing
    )

    # there are issues with calibrating the smb model to individual geotiles as glaciers
    # can cross multiple geotiles. To tackle this issue we to disconected geotile groups

    ddate = dims(dv_altim, :date)
    dΔheight = dims(smb, :Δheight)
    dpscale = dims(smb, :pscale)
    ddate_gemb = dims(smb, :date)
    dgeotile = dims(smb, :geotile)

    volume2mass = δice / 1000

    if !isnothing(single_geotile_test)
        geotile_ref = single_geotile_test

        # extract all geotiles with overlapping groups
        geotile_index = findfirst(==(single_geotile_test), geotiles0.id)
        group_index = findall(==(geotiles0[geotile_index, :group]), geotiles0.group)
        geotiles0 = geotiles0[group_index, :]

        # extract all geotiles with overlapping groups
    else
        geotile_ref = dgeotile[1]
    end

    # find common overlap 
    index = .!isnan.(smb[At(geotile_ref), :, 1, 1])
    ex_gemb = extrema(ddate_gemb[index])
    index = .!isnan.(dv_altim[At(geotile_ref), :])
    ex_dv = extrema(ddate[vec(index)])
    ex = (max(ex_gemb[1], ex_dv[1]), min(ex_gemb[2], ex_dv[2]))
    index_dv = (ddate .>= ex[1]) .& (ddate .<= ex[2])
    index_gemb = (ddate_gemb .>= ex[1]) .& (ddate_gemb .<= ex[2])

    decyear = decimalyear.(ddate_gemb[index_gemb])
    Δdecyear = decyear .- mean(decyear)

    # loop for each geotile
    geotiles = DataFrame()
    geotiles[!, :id] = geotiles0[:, :id]
    geotiles[!, :group] = geotiles0[:, :group]
    geotiles[!, :extent] = geotiles0[:, :extent]
    geotiles[!, :area_km2] = geotiles0[:, :area_km2]

    geotiles[!, :pscale] .= 1.0
    geotiles[!, :Δheight] .= 0.0
    geotiles[!, :rmse] .= 0.0

    if !isnothing(single_geotile_test)
        cost_metric_minimum = Inf
        best_rmse = Inf
        best_fit = zeros(ddate)
        best_smb = zeros(ddate)
        best_fac = zeros(ddate)
        best_discharge = zeros(ddate)
        best_offset = 0.

        # Can not yet modify DateTime auto ticks in Makie... convert to decimal year
        ddate = dims(smb, :date)
        decyear = decimalyear.(collect(val(ddate)))[index_gemb]
        xtickspacing = 5
        xlims = (floor(Int, minimum(decyear ./ xtickspacing)) * xtickspacing, ceil(Int, maximum(decyear ./ xtickspacing)) * xtickspacing)

        fig = []
    end

    # loop for each group
    # NOTE: do not do a groupby on geotiles as you need to keep the original geotiles order

    # Threads only buys you a 2x speedup with 128 cores... This could likely be made more efficient by better grouping the data for Threads
    for grp in unique(geotiles.group)

        gindex = grp .== geotiles.group

        area_km2 = sum(sum(geotiles[gindex, :area_km2]))
        discharge_km3yr = sum(discharge0[varname=At("discharge_gtyr"), geotile=At(geotiles.id[gindex])]) ./ volume2mass

        pscale0 = 1
        Δheight0 = 0

        dv0 = dv_altim[At(geotiles[gindex, :id]), index_dv]
        dv0 = dropdims(sum(dv0, dims=:geotile), dims=:geotile)

        cost_metric = fill(Inf, (dpscale, dΔheight))
        rmse = fill(Inf, (dpscale, dΔheight))

        Threads.@threads for pscale in dpscale
            #pscale = first(dpscale)

            for Δheight in dΔheight
                #Δheight = 0.

                smb0 = smb[At(geotiles[gindex, :id]), index_gemb, At(pscale), At(Δheight)] # smb is in unit of km3 assuming an ice density of 910 kg/m3
                smb0 = dropdims(sum(smb0, dims=:geotile), dims=:geotile)

                fac0 = fac[At(geotiles[gindex, :id]), index_gemb, At(pscale), At(Δheight)] # fac is in unit of km3
                fac0 = dropdims(sum(fac0, dims=:geotile), dims=:geotile)

                res = dv0 .- (smb0 .+ fac0 .- (discharge_km3yr .* Δdecyear))

                if any(isnan.(res))
                    println("dv0 has NaNs: $(any(isnan.(dv0)))")
                    println("smb0 has NaNs: $(any(isnan.(smb0)))")
                    println("fac0 has NaNs: $(any(isnan.(fac0)))")
                    println("discharge_km3yr has NaNs: $(any(isnan.(discharge_km3yr)))")
                    #display(discharge_km3yr)
                    #println(cname)

                    error("NANs in residual")
                end

                # align the center of the residual distribution to zero
                res0 = median(res)
                res = res .- res0

                if any(isnan.(res))
                    #display(dv0)
                    #display(discharge_km3yr)
                    #println(cname)

                    error("NANs in residual")
                end

                # try using the root mean square error instead of MAD...MAD did a poor job of caturing seasonality

                # ----------------- Objective function definition -[this should be passed as a kwarg]- ----------------
                # fit a model to the residuals... use the amplitude of the residuals as a weighting of the objective function
                # fit1 = curve_fit(model3, Δdecyear, res, p3)

                # remove trend for improved fit to seasonality

                rmse[pscale=At(pscale), Δheight=At(Δheight)], cost_metric[pscale=At(pscale), Δheight=At(Δheight)] = model_fit_cost_function(res, pscale, Δheight; seasonality_weight, distance_from_origin_penalty)
                # --------------------------------------------------------------

                if .!isnothing(single_geotile_test)
                    if cost_metric_minimum > cost_metric[pscale=At(pscale), Δheight=At(Δheight)]
                        cost_metric_minimum = cost_metric[pscale=At(pscale), Δheight=At(Δheight)]

                        best_rmse = rmse[pscale=At(pscale), Δheight=At(Δheight)]
                        best_fit = (smb0 .+ fac0 .- (discharge_km3yr .* Δdecyear)) .+ res0
                        best_smb = smb0
                        best_fac = fac0
                        # thi is +/- is done to get dates right for ploting
                        best_discharge = discharge_km3yr .* Δdecyear .+ best_fac .- best_fac
                        best_offset = res0

                        pscale0 = pscale
                        Δheight0 = Δheight
                    end
                end
            end
        end

        # optimal fit (minimizes objective function)
        index_minimum = argmin(cost_metric)
        geotiles[gindex, :pscale] .= DimPoints(cost_metric)[index_minimum][1]
        geotiles[gindex, :Δheight] .= DimPoints(cost_metric)[index_minimum][2]
        geotiles[gindex, :rmse] .= rmse[index_minimum]

        if .!isnothing(single_geotile_test)
            v2h = 1000 / area_km2;
            fig = plot_model_fit(best_smb .* v2h, best_fac .* v2h, best_discharge .* v2h, best_fit .* v2h, dv0 .* v2h, cost_metric; xlims, xtickspacing, colormap=:thermal)
            printstyled("best fit for $single_geotile_test [pscale = $pscale0, Δheight = $Δheight0, RMSE = $(round(best_rmse; digits=2))m]"; color=:light_blue)
        end
    end

    if .!isnothing(single_geotile_test)
        return geotiles[:, [:id, :extent, :pscale, :Δheight, :rmse]], fig
    else
        return geotiles[:, [:id, :extent, :pscale, :Δheight, :rmse]]
    end
end

"""
    add_pscale_classes(var0, dpscale_new; allow_negative=false)

Interpolate and extrapolate values across a new set of precipitation scale factors.

# Arguments
- `var0`: Input dimensional array with dimensions (:date, :height, :pscale, :Δheight)
- `dpscale_new`: New precipitation scale factor dimension to interpolate/extrapolate to
- `allow_negative`: Whether to allow negative values in the output (default: false)

# Returns
A new dimensional array with the same structure as `var0` but with the new pscale dimension.

This function performs linear interpolation for pscale values within the range of the original
data, and linear extrapolation for values outside that range. For extrapolation, it uses either
the first/last three points or all available points depending on data availability.
"""
function add_pscale_classes(var0, dpscale_new; allow_negative=false)

    ddate = dims(var0, :date)
    dheight = dims(var0, :height)
    dpscale = dims(var0, :pscale)
    dΔheight = dims(var0, :Δheight)

    interpolate_index = (dpscale_new .>= minimum(dpscale)) .& (dpscale_new .<= maximum(dpscale))
    extrapolate_lower_index = .!interpolate_index
    extrapolate_lower_index[findfirst(interpolate_index):end] .= false
    extrapolate_upper_index = .!interpolate_index
    extrapolate_upper_index[1:findlast(interpolate_index)] .= false
    extrapolate_index = extrapolate_lower_index .| extrapolate_upper_index
    x = collect(dpscale)
    x_new = collect(dpscale_new)

    # this can be removed after code rerun
    x_sort_perm = sortperm(x)
    x = x[x_sort_perm]

    M = hcat(ones(length(x)), x)

    gemb_new = fill(NaN, (ddate, dheight, dpscale_new, dΔheight))

    valid_range, = validrange(.!isnan.(var0[:, 1, 1, 1]))

    Threads.@threads for date in ddate[valid_range]
        #date = ddate[820]

        y_new = fill(NaN, length(x_new))

        for height in dheight
            #height = 1050

            for Δheight in dΔheight
                #Δheight = 0

                y = var0[At(date), At(height), :, At(Δheight)]

                # this can be removed after code rerun
                y[:] = y[x_sort_perm]

                if !allow_negative
                    valid_interp = y .>= 0
                end

                if sum(valid_interp) >= 2

                    # Loess interpolation on so few points == linear interpolation but is much slower
                    #model = Loess.loess(x, y, span=0.3)
                    #y_new[interpolate_index] = Loess.predict(model, x_new[interpolate_index])

                    model = DataInterpolations.LinearInterpolation(y, x)
                    y_new[interpolate_index] = model(x_new[interpolate_index])

                    x0 = x[valid_interp]
                    y0 = y[valid_interp]
                    M0 = M[valid_interp, :]

                    if length(x0) > 3
                        lowest_pts = 1:3
                        highest_pts = length(x0)-2:length(x0)

                        param1 = M0[lowest_pts, :] \ y0[lowest_pts]
                        y_new[extrapolate_lower_index] = param1[1] .+ x_new[extrapolate_lower_index] .* param1[2]

                        param2 = M0[highest_pts, :] \ y0[highest_pts]
                        y_new[extrapolate_upper_index] = param2[1] .+ x_new[extrapolate_upper_index] .* param2[2]
                    else
                        param1 = M0 \ y0
                        y_new[extrapolate_index] = param1[1] .+ x_new[extrapolate_index] .* param1[2]
                    end
                
                    if !allow_negative
                        y_new[y_new.<0] .= 0
                    end

                    # sanity check
                    # p = lines(x_new, y_new);
                    # CM.plot!(x,collect(y))
                else
                    y_new[:] .= 0
                end

                gemb_new[At(date), At(height), :, At(Δheight)] = y_new
            end
        end
    end

    return gemb_new
end


"""
    gemb_calibration(path2runs; surface_mask_default="glacier", gemb_run_id, geotile_grouping_min_feature_area_km2=100, single_geotile_test=nothing, force_remake_before)

Calibrate GEMB (Glacier Energy and Mass Balance) model to altimetry data by finding optimal fits for grouped geotiles.

# Arguments
- `path2runs`: Array of file paths to binned synthesized data files
- `surface_mask_default`: Default surface mask type (default: "glacier")
- `gemb_run_id`: GEMB run identifier
- `geotile_grouping_min_feature_area_km2`: Minimum feature area in km² for geotile grouping (default: 100)
- `single_geotile_test`: Single geotile ID for testing (default: nothing)
- `force_remake_before`: Date to force file regeneration (default: nothing)

# Description
Loads GEMB data and altimetry volume change data, then finds optimal model fits for grouped geotiles.
Saves calibration results as Arrow files with geometry information.
"""
function gemb_calibration(
    path2runs, 
    discharge,
    gemb;
    geotile_width,
    geotile_grouping_min_feature_area_km2=100, 
    single_geotile_test=nothing, 
    plots_show=false,
    plots_save=false,
    seasonality_weight=95/100,
    distance_from_origin_penalty=2/100,
    force_remake_before=nothing
    )

    # Load GEMB data
    # Note: GEMB volume change is for a single surface mask - this is a hack since making GEMB dv for multiple surface masks is onerous
    # However, since GEMB data is calibrated to altimetry that uses different surface masks, this is mostly a non-issue

    if !isnothing(single_geotile_test)
        # NOTE: don't prematurely trim data to `single_geotile_test` as fit optimization is done by geotile groups
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)], OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!"
    end
   
    dgeotile = dims(gemb["smb"], :geotile)
    params = binned_filled_fileparts.(path2runs)
    surface_masks = unique(getindex.(params, :surface_mask))

    # Process each surface mask to load aligned geotiles and calculate hypsometry
    geotiles = Dict()
    area_km2 = Dict()
    discharge0 = Dict()
    for surface_mask in surface_masks # takes 7 seconds

        # Load and align geotiles with the synthesized data dimensions
        geotiles[surface_mask] = _geotile_load_align(;
            surface_mask,
            geotile_width,
            only_geotiles_w_area_gt_0=true,
            geotile_order=collect(dgeotile)
        )

        area_km2[surface_mask] = _geotile_area_km2(; surface_mask, geotile_width)[geotile=At(val(dgeotile))]

        # Add in groups
        geotiles_groups = geotile_grouping(; surface_mask, min_area_km2=geotile_grouping_min_feature_area_km2, geotile_width, force_remake_before)

        # Add groups to geotiles
        _, grp_index = intersectindices(geotiles[surface_mask].id, geotiles_groups.id)
        geotiles[surface_mask][!, :group] = geotiles_groups[grp_index, :group]

        discharge0[surface_mask] = discharge2geotile(discharge, geotiles[surface_mask])
    end
   
    @showprogress desc = "Finding optimal GEMB fits to altimetry data..."  Threads.@threads for binned_synthesized_file in path2runs

        synthesized_gemb_fit = replace(binned_synthesized_file, ".jld2" => "_gemb_fit.arrow")

        if isfile(synthesized_gemb_fit) && (isnothing(force_remake_before) || Dates.unix2datetime(mtime(synthesized_gemb_fit)) > force_remake_before) && isnothing(single_geotile_test)
            printstyled("    -> Skipping $(synthesized_gemb_fit) because it was created after force_remake_before:$force_remake_before\n"; color=:light_green)
            continue
        else
            
            t1 = time()

            run_parameters = binned_filled_fileparts(synthesized_gemb_fit)
            surface_mask = run_parameters.surface_mask

            dh = FileIO.load(binned_synthesized_file, "dh_hyps")

            # Convert elevation change to volume change
            dv_altim = dh2dv_geotile(dh, area_km2[surface_mask])
            all_nans = dropdims(all(isnan.(dv_altim), dims=:date), dims=:date)
            
            if any(all_nans)
                error("geotile volume change contains all NaNs")
            end

            # Find optimal fit to GEMB data
            # There are issues with calibrating the SMB model to individual geotiles since glaciers 
            # can cross multiple geotiles, therefore we calibrate the model for groups of 
            # distinct geotiles.

            if !isnothing(single_geotile_test)
                @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)], OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!"
                df, fig = gemb_bestfit_grouped(dv_altim, gemb["smb"], gemb["fac"], discharge0[surface_mask], geotiles[surface_mask]; single_geotile_test, seasonality_weight, distance_from_origin_penalty)

                param = binned_filled_fileparts(binned_synthesized_file)
                fig_folder = pathlocal[:figures]

                if occursin("binned_unfiltered", param.binned_folder)
                    fig_folder = joinpath(fig_folder, "binned_unfiltered")
                else
                    fig_folder = joinpath(fig_folder, "binned")
                end

                binned_synthesized_file_parts = splitpath(binned_synthesized_file)
                plot_save_path_prefix = joinpath(fig_folder, replace(binned_synthesized_file_parts[end], ".jld2" => ""))

                plots_show && display(fig)
                plots_save && save(plot_save_path_prefix * "_$(single_geotile_test)_model_fit.png", fig)

                return
    
            else
                df = gemb_bestfit_grouped(dv_altim, gemb["smb"], gemb["fac"], discharge0[surface_mask], geotiles[surface_mask]; single_geotile_test, seasonality_weight, distance_from_origin_penalty)
            end

            df[!, :area_km2] = sum.(geotiles[surface_mask].area_km2)
            df.rmse = df.rmse ./ df.area_km2 * 1000
            rename!(df, "rmse" => "rmse_m")
            df.geometry = extent2rectangle.(df.extent)
            df = df[:, Not(:extent)]

            GeoDataFrames.write(synthesized_gemb_fit, df)
            println("\n $binned_synthesized_file optimal GEMB fit found: $(round(Int,time() -t1))s")
        end
    end
end

"""
    gemb2geotile(gemb, geotile_row; geotile_buffer=1000, min_gemb_coverage=0.75, ddate::Dim{:date})

Aggregate GEMB model output to a geotile by averaging over grid points within the geotile extent and height bins, for each precipitation scaling factor.

# Arguments
- `gemb`: Dictionary of GEMB output arrays.
- `geotile_row`: Row from geotile table, containing extent and area information.
- `geotile_buffer`: Buffer distance (meters) to expand geotile extent if coverage is insufficient (default: 1000).
- `min_gemb_coverage`: Minimum fraction of geotile area that must be covered by GEMB data (default: 0.75).
- `ddate`: Date dimension (Dim{:date}).

# Returns
- `Dict` with variables averaged within the geotile, binned by height and precipitation scale.

# Notes
- Buffers geotile extent iteratively until minimum coverage is achieved.
- Returns zeros if geotile area is zero.
"""
function gemb2geotile(gemb, geotile_row; geotile_buffer=1000, min_gemb_coverage=0.75, ddate::Dim{:date})

    if min_gemb_coverage > 1
        error("min_gemb_coverage must be less than 1 but is $min_gemb_coverage")
    end

    unique_precipitation_scale = unique(gemb["precipitation_scale"])
    vars = setdiff(collect(keys(gemb)), ("precipitation_scale", "elevation_delta", "latitude", "longitude", "height", "date", "height_effective"))

    height_range, height_center = project_height_bins()
    dheight = Dim{:height}(height_center)
    dpscale = Dim{:pscale}(sort(unique_precipitation_scale))

    gemb0 = Dict()
    for k in vcat(vars, "nobs")
        gemb0[k] = fill(NaN, (ddate, dheight, dpscale))
    end

    geotile_extent = geotile_row.extent
    latitude_distance, longitude_distance = meters2lonlat_distance(geotile_buffer, mean(geotile_extent.Y))
    has_data = falses(size(height_center))
    in_geotile = falses(size(gemb["longitude"]))

    total_area = sum(geotile_row.area_km2)
    if total_area == 0
        for k in vars
            gemb0[k][:] .= 0.
        end
        return gemb0
    end

    cnt = 0
    while (sum(geotile_row.area_km2[has_data]) / total_area) < min_gemb_coverage
        in_geotile = vec([within(geotile_extent, x, y) for (x, y) in zip(gemb["longitude"], gemb["latitude"])])
        h = gemb["height_effective"][in_geotile]

        for i in 1:(length(height_range)-1)
            has_data[i] = any((h .>= height_range[i]) .& (h .< height_range[i+1]))
        end

        geotile_extent = Extents.buffer(geotile_extent, (X=longitude_distance, Y=latitude_distance))
        cnt += 1
    end

    if cnt > 10
        @warn "Buffering $(geotile_row.id) extent by $(round(Int,geotile_buffer*cnt/1000)) km, to achieve minimum GEMB coverage of $min_gemb_coverage"
    end

    Threads.@threads for i in 1:(length(height_range)-1)
        index_height = (gemb["height_effective"] .>= height_range[i]) .& (gemb["height_effective"] .< height_range[i+1])
        height0 = height_center[i]

        if !any(index_height)
            continue
        end

        for pscale in dpscale
            index_precipitation = gemb["precipitation_scale"] .== pscale
            index = index_precipitation .& in_geotile .& index_height
            sum_index = sum(index)

            if sum_index == 0
                continue
            elseif sum_index == 1
                for k in vars
                    gemb0[k][:, At(height0), At(pscale)] = gemb[k][index, :]
                end
            else
                for k in vars
                    gemb0[k][:, At(height0), At(pscale)] = mean(gemb[k][index, :], dims=1)
                end
            end

            valid_date = .!isnan.(gemb0[vars[1]][:, At(height0), At(pscale)])
            gemb0["nobs"][valid_date, At(height0), At(pscale)] .= sum(index)
        end
    end

    return gemb0
end

"""
    gemb_fill_gaps_and_add_Δheight_classes(gemb, dΔheight)

Fill gaps in GEMB geotile data by interpolating and extrapolating missing values, then generate new classes for elevation change (Δheight) by shifting the filled data accordingly.

# Arguments
- `gemb`: Dictionary of GEMB output arrays, indexed by variable name.
- `dΔheight`: Range or dimension object specifying Δheight classes to generate.

# Returns
- `Dict` containing filled and Δheight-shifted GEMB variables, with physical constraints enforced and cumulative variables reconstructed.

# Notes
- "smb" and "runoff" are excluded from filling and are reconstructed after processing to maintain mass conservation.
- Filling is performed on rates (not cumulative values) for physical consistency.
- Elevation classes are created by shifting the filled data vertically according to Δheight.
"""
function gemb_fill_gaps_and_add_Δheight_classes(gemb, dΔheight)
    # Exclude "smb" and "runoff" from filling; reconstruct later for mass conservation
    vars = setdiff(collect(keys(gemb)), ["nobs", "smb", "runoff"])

    k = first(keys(gemb))
    ddate = dims(gemb[k], :date)
    dheight = dims(gemb[k], :height)
    dpscale = dims(gemb[k], :pscale)
    x = collect(dheight)
    npts_linear_extrapolation = 7
    height_bin_interval = dΔheight.val.step.hi

    Threads.@threads for k in vars
        for pscale in dpscale
            # Interpolate/extrapolate on rates for physical consistency
            var0 = gemb[k][:, :, At(pscale)]
            dheight = dims(var0, :height)
            x0 = (1:length(dheight)) ./ length(dheight) .- 0.5
            M0 = hcat(ones(length(x0)), x0)
            npts = npts_linear_extrapolation

            # Skip if all values are zero
            if all(var0 .== 0)
                continue
            end

            # Calculate rates (except for "fac")
            if k == "fac"
                var0_rate = var0[2:end, :]
            else
                ind = findfirst(.!vec(all(isnan.(var0), dims=:height))) - 1
                var0[ind, :] .= 0
                var0_rate = var0[2:end, :] .- collect(var0[1:end-1, :])
            end

            for date in dims(var0_rate, :date)
                y = var0_rate[At(date), :]
                valid = .!isnan.(collect(y))

                if all(valid) || !any(valid)
                    continue
                end

                if k in ("acc", "refreeze", "melt", "rain", "fac")
                    valid = valid .& (y .>= 0)
                end

                if all(valid) || !any(valid)
                    y[:] .= 0
                else
                    validgap = validgaps(valid)
                    if sum(valid) > 4
                        model = loess(x[valid], y[valid], span=0.3)
                        y[valid .| validgap] = predict(model, x[valid .| validgap])
                    end

                    fill_index = isnan.(collect(y))
                    if k == "ec"
                        valid_interp_index = .!fill_index
                    else
                        valid_interp_index = .!fill_index .& (y .!= 0)
                    end

                    if any(fill_index) && (sum(valid_interp_index) > 3)
                        if sum(valid_interp_index) < npts + 2
                            M = @view M0[valid_interp_index, :]
                            param1 = M \ y[valid_interp_index]
                            y[fill_index] = param1[1] .+ x0[fill_index] .* param1[2]
                        else
                            fill_index_lower = copy(fill_index)
                            fill_index_lower[findfirst(.!fill_index):end] .= false

                            fill_index_upper = copy(fill_index)
                            fill_index_upper[1:findfirst(.!fill_index_upper)] .= false

                            valid_interp_index_lower = findall(valid_interp_index)[1:npts]
                            valid_interp_index_upper = findall(valid_interp_index)[end-npts+1:end]

                            if any(fill_index_lower)
                                M = @view M0[valid_interp_index_lower, :]
                                param1 = M \ y[valid_interp_index_lower]
                                y[fill_index_lower] = param1[1] .+ x0[fill_index_lower] .* param1[2]
                            end

                            if any(fill_index_upper)
                                M = @view M0[valid_interp_index_upper, :]
                                param1 = M \ y[valid_interp_index_upper]
                                y[fill_index_upper] = param1[1] .+ x0[fill_index_upper] .* param1[2]
                            end
                        end

                        if k in ("acc", "refreeze", "melt", "rain", "fac")
                            y[y .< 0] .= 0
                        end
                    else
                        y[fill_index] .= 0
                    end
                end

                gemb[k][At(date), :, At(pscale)] = y
            end
        end
    end

    # Enforce physical constraints
    gemb = gemb_rate_physical_constraints!(gemb)
    gemb["fac"][gemb["fac"] .< 0] .= 0

    gemb0 = Dict()
    for k in vars
        gemb0[k] = deepcopy(gemb[k][:, :, :])
    end

    gemb_mie = Dict()
    for k in vars
        gemb_mie[k] = fill(NaN, (ddate, dheight, dpscale, dΔheight))
    end

    Threads.@threads for k in vars
        exclude_zeros_in_extrapolation = (k != "ec")
        for Δheight in dΔheight
            for pscale in dpscale
                if k == "nobs"
                    gemb_mie[k][:, :, At(pscale), At(Δheight)] = gemb0[k][:, :, At(pscale)]
                else
                    height_shift = round(Int16(Δheight / height_bin_interval))
                    gemb_mie[k][:, :, At(pscale), At(Δheight)] =
                        _matrix_shift_ud!(deepcopy(gemb0[k][:, :, At(pscale)]), height_shift;
                            exclude_zeros_in_extrapolation, npts_linear_extrapolation)
                end
            end
        end
    end

    gemb_mie = gemb_rate_physical_constraints!(gemb_mie)
    gemb_mie["fac"][gemb_mie["fac"] .< 0] .= 0

    # Reconstruct cumulative variables from rates
    for k in setdiff(collect(keys(gemb_mie)), ["nobs", "fac"])
        sindex = findfirst(.!isnan.(gemb_mie[k][:, 1, 1, 1]))
        gemb_mie[k][sindex:end, :, :, :] = cumsum(gemb_mie[k][sindex:end, :, :, :], dims=:date)
    end

    return gemb_mie
end


"""
    read_gemb_files(gemb_files, gembinfo; vars2extract=["acc", "refreeze", "melt", "rain", "ec", "fac", "smb", "dv"], date_range)

Read and merge multiple GEMB output files into a single dictionary of arrays, extracting specified variables and metadata.

# Arguments
- `gemb_files`: Array of GEMB output file paths.
- `gembinfo`: Struct or dictionary containing precipitation scale and elevation delta information.
- `vars2extract`: Variables to extract from each file (default: common GEMB variables).
- `date_range`: Tuple or array specifying the date range to extract.

# Returns
- `Dict` containing merged arrays for each variable, with additional metadata fields and effective height.
"""
function read_gemb_files(gemb_files, gembinfo; vars2extract=["acc", "refreeze", "melt", "rain", "ec", "fac", "smb", "dv"], date_range)
    datebin_edges = decimalyear.(date_range)

    gemb = []
    @showprogress desc = "Reading raw GEMB output into memory, this will take ~2 min [peak memory usage: 150GB]" Threads.@threads for gemb_file in gemb_files
        gemb0 = gemb_read2(gemb_file; vars=vars2extract, datebin_edges)

        # Determine precipitation scale and elevation delta indices
        if (length(gembinfo.elevation_delta) == 1) && (length(gembinfo.precipitation_scale) == 1)
            precip_scale_ind = "p1"
            elev_delta_ind = "t1"
        else
            ind1 = findlast("_p", gemb_file)
            ind2 = findlast("_t", gemb_file)
            ind3 = findlast(".mat", gemb_file)
            precip_scale_ind = gemb_file[ind1[end]:ind2[1]-1]
            elev_delta_ind = gemb_file[ind2[end]:ind3[1]-1]
        end

        gemb0["precipitation_scale"] = gembinfo.precipitation_scale[At(precip_scale_ind)]
        gemb0["elevation_delta"] = gembinfo.elevation_delta[At(elev_delta_ind)]

        push!(gemb, gemb0)
    end

    # If only one class, merge all into a single dictionary
    if (length(gembinfo.elevation_delta) == 1) && (length(gembinfo.precipitation_scale) == 1)
        gemb0 = Dict()
        for k in keys(gemb[1])
            if k !== "date"
                gemb0[k] = reduce(vcat, [g[k] for g in gemb])
            else
                gemb0[k] = gemb[1][k]
            end
        end
        gemb = [gemb0]
    end

    # Merge all simulations into a dictionary of arrays
    dates = vec(gemb[1]["date"])
    latitude = vec(reduce(vcat, getindex.(gemb, Ref("latitude"))))
    longitude = vec(reduce(vcat, getindex.(gemb, Ref("longitude"))))
    height0 = vec(reduce(vcat, getindex.(gemb, Ref("height"))))

    # Ensure longitude is in -180 to 180 range
    index = longitude .< -180
    longitude[index] = longitude[index] .+ 360
    index = longitude .> 180
    longitude[index] = longitude[index] .- 360

    precipitation_scale0 = Float64[]
    elevation_delta0 = Float64[]
    for g in gemb
        precipitation_scale0 = vcat(precipitation_scale0, fill(g["precipitation_scale"], length(g["latitude"])))
        elevation_delta0 = vcat(elevation_delta0, fill(g["elevation_delta"], length(g["latitude"])))
    end

    gemb0 = Dict(
        "date" => decimalyear2datetime.(dates),
        "latitude" => latitude,
        "longitude" => longitude,
        "height" => height0,
        "precipitation_scale" => precipitation_scale0,
        "elevation_delta" => elevation_delta0
    )

    foo = Dict()
    for k in setdiff(vars2extract, ["date", "latitude", "longitude", "height"])
        foo[k] = reduce(vcat, getindex.(gemb, Ref(k)))
    end
    gemb = merge(foo, gemb0)

    # Add effective height (height + elevation_delta)
    gemb["height_effective"] = gemb["height"] .+ gemb["elevation_delta"]

    return gemb
end


"""
    gemb_dv(gemb, area_km2, dpscale_new)

Compute volume change (ΔV) for GEMB output, interpolated to new precipitation scaling factors.

# Arguments
- `gemb`: Dictionary of GEMB output arrays, each keyed by variable name.
- `area_km2`: Array of grid cell areas in km².
- `dpscale_new`: Array of new precipitation scaling factors to interpolate to.

# Returns
- `Dict` mapping each variable to a 3D array of volume change (date × pscale × Δheight).
"""
function gemb_dv(gemb, area_km2, dpscale_new)
    k = first(keys(gemb))
    dpscale = dims(gemb[k], :pscale)
    dΔheight = dims(gemb[k], :Δheight)
    ddate = dims(gemb[k], :date)

    # Ensure all old pscales are present in new pscales
    if !isempty(setdiff(dpscale, dpscale_new))
        error("old pscales [$(dpscale)] must exist in new pscales [$(dpscale_new)]")
    end

    # Initialize output dictionary
    gemb_dv = Dict()
    for k in keys(gemb)
        gemb_dv[k] = fill(NaN, (ddate, dpscale_new, dΔheight))
    end

    Threads.@threads for k in collect(keys(gemb))
        valid = .!isnan.(gemb[k])
        (date_range, height_range, pscale_range, Δheight_range) = validrange(valid)
        for date in ddate[date_range]
            for Δheight in dΔheight
                dv = @d dropdims(sum(gemb[k][date=At(date), Δheight=At(Δheight)][height_range, :] .* area_km2[height_range] ./ 1000, dims=:height), dims=:height)
                dv_interp = DataInterpolations.LinearInterpolation(dv, val(dpscale))
                gemb_dv[k][date=At(date), Δheight=At(Δheight)] = dv_interp(val(dpscale_new))
            end
        end
    end
    return gemb_dv
end


"""
    gemb_ensemble_dv(; gemb_run_id, vars2load=["runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze"])

Load GEMB geotile ensemble data for a given run ID. Returns a dictionary of variables specified in `vars2load` from the corresponding GEMB geotile volume change file. If `vars2load` is `nothing`, loads all variables from the file.

# Arguments
- `gemb_run_id`: Identifier for the GEMB run.
- `vars2load`: Array of variable names to load (default: common mass balance variables). If `nothing`, loads all variables.

# Returns
- `Dict` containing the requested variables.
"""
function gemb_ensemble_dv(; gemb_run_id, vars2load = ["runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze"])
    gembinfo = gemb_info(; gemb_run_id)
    gemb_geotile_filename_dv = replace(gembinfo.filename_gemb_combined, ".jld2" => "_geotile_dv.jld2")
    if isnothing(vars2load)
        gemb = load(gemb_geotile_filename_dv)
    else
        gemb = Dict()
        for varname in vars2load
            gemb[varname] = load(gemb_geotile_filename_dv, varname)
        end
    end
    return gemb
end
