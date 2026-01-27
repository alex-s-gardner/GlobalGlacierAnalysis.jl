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
function gemb_read2(gemb_file; vars="all", datebin_edges=nothing)

    vars0 = ["latitude", "longitude", "date", "smb", "fac", "ec", "acc", "runoff", "melt", "fac_to_depth", "height", "refreeze", "t_air", "rain"]

    # variables that need to be converted from cumulitive outputs to rates

    var_cumulitive = ["smb", "ec", "acc", "runoff", "melt", "refreeze", "rain"]

    varmap = Dict(
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
            f = foo[v[2]]
            f[isnan.(f)] .= 0 #set any random nan to zero
            gemb[v[1]] = f
        end
    else
        if length(intersect(vars0, vars)) != length(vars)
            error("one or more vars not recognized")
        end

        gemb = matopen(gemb_file) do file
            foo = MAT.read.(Ref(file), [varmap[v] for v in vars])
            for (i, v) in enumerate(vars)
                f = foo[i]
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
        ind = Vector{Union{UnitRange{Int64},Nothing}}()
        t = vec(gemb["date"])
        for i = eachindex(datebin_edges)[1:end-1]
            a = findfirst((t .>= datebin_edges[i]) .& (t .< datebin_edges[i+1]))
            b = findlast((t .>= datebin_edges[i]) .& (t .< datebin_edges[i+1]))

            if isnothing(a)
                push!(ind, nothing)
            else
                push!(ind, a:b)
            end
        end

        for k in keys(gemb)
            if k == "date"
                continue
            end

            if size(gemb[k], 2) == size(gemb["date"], 2)
                foo = fill(NaN, (size(gemb[k], 1), length(ind)))
                for i in eachindex(ind)
                    if !isnothing(ind[i])
                        foo[:, i] = mean(gemb[k][:, ind[i]], dims=2)
                    end
                end
                gemb[k] = foo
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
    dΔT = dims(gemb[first(vars)], :ΔT)
    dpscale = dims(gemb[first(vars)], :pscale)
    ddate = dims(gemb[first(vars)], :date)

    df = DataFrame()

    area = gemb["area_km2"]
    for var0 in vars
        #var0 = first(vars)

        v0 = gemb[var0]
        nobs0 = gemb["nobs"]
        for rgi in drgi
            for ΔT in dΔT
                for pscale in dpscale

                    vout = vec(v0[At(rgi), :, At(pscale), At(ΔT)])
                    nobsout = vec(nobs0[At(rgi), :, At(pscale), At(ΔT)])

                    append!(df,
                        DataFrame(; rgi, var="$(var0)_km3", mission="gemb", ΔT, pscale, surface_mask="glacier", area_km2=area[At(rgi)], val=[vout], nobs=[nobsout], filename=path2file)
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
- `df`: DataFrame containing GEMB data with pscale and ΔT parameters
- `n_densify`: Number of interpolated points to add between existing parameter values (default: 4)

# Returns
The input DataFrame with additional interpolated parameter combinations added.
"""
function gemb_classes_densify!(df; n_densify=4)

    # interpolate between gemb classes
    classes_ΔT = sort(unique(df.ΔT))
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
            for h = eachindex(classes_ΔT)
                h0 = df.ΔT .== classes_ΔT[h]

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

                for h = eachindex(classes_ΔT)[1:end-1]
                    h1 = df.ΔT .== classes_ΔT[h]
                    h2 = df.ΔT .== classes_ΔT[h+1]

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
                            foo.ΔT = (w .* classes_ΔT[h]) .+ ((1 - w) .* classes_ΔT[h+1])
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
    h0 = (1.:n) ./ n .- 0.5
    M = hcat(ones(n), h0)

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
    gemb["runoff"] = melt .- refreeze

    # SMB must equal acc - runoff - ec [rain is not included in SMB calculation]
    gemb["smb"] = acc .- gemb["runoff"] .- gemb["ec"]

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
DataFrame with calibrated parameters for each geotile: id, extent, pscale, ΔT, and rmse
"""

function gemb_bestfit_grouped(dv_altim, dv_gemb, geotiles0;
    distance_from_origin_penalty=20 / 100,
    ΔT_to_pscale_weight=1,
    seasonality_weight=85 / 100,
    calibrate_to=:all,
    is_scaling_factor=Dict("pscale" => true, "ΔT" => true)
)

    if ΔT_to_pscale_weight > 1 || ΔT_to_pscale_weight < 0
        error("ΔT_to_pscale_weight must be between 0 and 1")
    end

    if seasonality_weight > 1 || seasonality_weight < 0
        error("seasonality_weight must be between 0 and 1")
    end
    
    # there are issues with calibrating the dv_gemb model to individual geotiles as glaciers
    # can cross multiple geotiles. To tackle this issue one needs to process mutually 
    # exclusive geotile groups

    # loop for each geotile
    geotiles = DataFrame()
    geotiles[!, :id] = geotiles0[:, :id]
    geotiles[!, :group] = geotiles0[:, :group]
    geotiles[!, :extent] = geotiles0[:, :extent]
    geotiles[!, :area_km2] = geotiles0[:, :area_km2]

    geotiles[!, :pscale] .= 1.0
    geotiles[!, :ΔT] .= 1.0
    geotiles[!, :cost] .= Inf

    dpscale = dims(dv_gemb, :pscale)
    dΔT = dims(dv_gemb, :ΔT)
 
    bounds = boxconstraints(lb=[minimum(dpscale.val), minimum(dΔT.val)].+.01, ub=[maximum(dpscale.val), maximum(dΔT.val)].-.01)


    # loop for each group
    # NOTE: do not do a groupby on geotiles as you need to keep the original geotiles order

    # Threads only buys you a 2x speedup with 128 cores... This could likely be made more efficient by better grouping the data for Threads
    if nrow(geotiles0) != 1
        Threads.@threads for grp in unique(geotiles.group)
            #@warn "threads turned off"
            #for grp in unique(geotiles.group)

            gindex = grp .== geotiles.group

            dv_altim0 = dropdims(sum(dv_altim[geotile=At(geotiles[gindex, :id])], dims=:geotile), dims=:geotile)
            dv_gemb0 = dropdims(sum(dv_gemb[geotile=At(geotiles[gindex, :id])], dims=:geotile), dims=:geotile)

            kwargs2 = (seasonality_weight=seasonality_weight, distance_from_origin_penalty=distance_from_origin_penalty, ΔT_to_pscale_weight=ΔT_to_pscale_weight, calibrate_to=calibrate_to, is_scaling_factor=is_scaling_factor)
            f2 = Base.Fix{1}(Base.Fix{2}(Base.Fix{4}(Base.Fix{5}(Base.Fix{6}(gemb_altim_cost!, kwargs2), dv_gemb0), dv_altim0), copy(dv_altim0)), copy(dv_altim0))

            # black box optimization IS NOT THREAD SAFE !!!!!!!
            result = optimize(f2, bounds, ECA(N=0, K=4, p_bin=0.05))

            (pscale_best, ΔT_best) = minimizer(result)

            if (!(minimum(dpscale.val) <= pscale_best <= maximum(dΔT.val))) || isnan(pscale_best)
                error("pscale_best is out of bounds -or is NaN")
            end

            if (!(minimum(dΔT.val) <= ΔT_best <= maximum(dΔT.val))) || isnan(ΔT_best)
                error("ΔT_best is out of bounds -or is NaN")
            end

            geotiles[gindex, :pscale] .= pscale_best
            geotiles[gindex, :ΔT] .= ΔT_best

            #results = bboptimize(f2, [1., 1.]; SearchRange=[(extrema(dpscale)), (extrema(dΔT))], MaxFuncEvals=1500, TraceMode=:silent, Method=:adaptive_de_rand_1_bin_radiuslimited)

            # (pscale_best, ΔT_best) = best_candidate(results)
            # geotiles[gindex, :pscale] .= pscale_best
            # geotiles[gindex, :ΔT] .= ΔT_best
            #println("pscale_fit: $(round(pscale_best, digits=1)), ΔT_fit: $(round(ΔT_best, digits=1))")
        end

        return geotiles[:, [:id, :extent, :pscale, :ΔT]]


    else
        gindex = geotiles.group[1] .== geotiles.group

        dv_altim0 = dropdims(sum(dv_altim[geotile=At(geotiles[gindex, :id])], dims=:geotile), dims=:geotile)
        dv_gemb0 = dropdims(sum(dv_gemb[geotile=At(geotiles[gindex, :id])], dims=:geotile), dims=:geotile)

        kwargs2 = (seasonality_weight=seasonality_weight, distance_from_origin_penalty=distance_from_origin_penalty, ΔT_to_pscale_weight=ΔT_to_pscale_weight, calibrate_to=calibrate_to, is_scaling_factor=is_scaling_factor)
        f2 = Base.Fix{1}(Base.Fix{2}(Base.Fix{4}(Base.Fix{5}(Base.Fix{6}(gemb_altim_cost!, kwargs2), dv_gemb0), dv_altim0), copy(dv_altim0)), copy(dv_altim0))

        step_size = 0.1
        s0 = -1 / minimum(dpscale.val)
        e0 = maximum(dpscale.val)
        pscale0 = s0:step_size:e0
        dpscale_search = vcat(-1 ./ pscale0[pscale0.<-1], pscale0[pscale0.>=1])

        if all(dΔT .> 0)
            step_size = 0.1
            s0 = -1 / minimum(dΔT.val)
            e0 = maximum(dΔT.val)
            ΔT0 = s0:step_size:e0
            dΔT_search = vcat(-1 ./ ΔT0[ΔT0.<-1], ΔT0[ΔT0.>=1])
        else
            dΔT_search = minimum(dΔT.val):((maximum(dΔT.val)-minimum(dΔT.val))/20):maximum(dΔT.val)
        end

        dpscale_search = Dim{:pscale}(dpscale_search)
        dΔT_search = Dim{:ΔT}(dΔT_search)

        cost = Inf;
        pscale_best = NaN;
        ΔT_best = NaN;

        cost = fill(NaN, dpscale_search, dΔT_search)
        for pscale in dpscale_search
            for  ΔT in dΔT_search
                cost[pscale=At(pscale), ΔT=At(ΔT)] = f2([pscale, ΔT])
            end
        end

        #println("pscale_best: $(round(pscale_best, digits=1)), ΔT_best: $(round(ΔT_best, digits=1))")
        return cost
    end
end

"""
    add_pscale_classes(var0, dpscale_new; allow_negative=false)

Interpolate and extrapolate values across a new set of precipitation scale factors.

# Arguments
- `var0`: Input dimensional array with dimensions (:date, :height, :pscale, :ΔT)
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
    dΔT = dims(var0, :ΔT)

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

    gemb_new = fill(NaN, (ddate, dheight, dpscale_new, dΔT))

    valid_range, = validrange(.!isnan.(var0[:, 1, 1, 1]))

    Threads.@threads for date in ddate[valid_range]
        #date = ddate[820]

        y_new = fill(NaN, length(x_new))

        for height in dheight
            #height = 1050

            for ΔT in dΔT
                #ΔT

                y = var0[At(date), At(height), :, At(ΔT)]

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

                gemb_new[At(date), At(height), :, At(ΔT)] = y_new
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
    gemb;
    geotile_width=2,
    geotile_grouping_min_feature_area_km2=100,
    single_geotile_test=nothing,
    seasonality_weight=95/100,
    distance_from_origin_penalty=2 / 100,
    ΔT_to_pscale_weight=1,
    force_remake_before=nothing,
    calibrate_to = :all
)

    # Load GEMB data
    # Note: GEMB volume change is for a single surface mask - this is a hack since making GEMB dv for multiple surface masks is onerous
    # However, since GEMB data is calibrated to altimetry that uses different surface masks, this is mostly a non-issue

    if !isnothing(single_geotile_test)
        # NOTE: don't prematurely trim data to `single_geotile_test` as fit optimization is done by geotile groups
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)], OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!"
    end

    dv_gemb = gemb[:dv];

    dgeotile = dims(dv_gemb, :geotile)
    dΔT = dims(dv_gemb, :ΔT)
    dpscale = dims(dv_gemb, :pscale)
    params = binned_filled_fileparts.(path2runs)
    surface_masks = unique(getindex.(params, :surface_mask))

    # Process each surface mask to load aligned geotiles and calculate hypsometry
    geotiles = Dict()
    area_km2 = Dict()


    if all(dΔT .> 0)
        ΔT_is_scaling = true;
    else
        ΔT_is_scaling = false;
    end

    if all(dpscale .> 0)
        pscale_is_scaling = true;
    else
        pscale_is_scaling = false;
    end

    is_scaling_factor = Dict("pscale" => pscale_is_scaling, "ΔT" => ΔT_is_scaling)


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
    
        if !isnothing(single_geotile_test)
            index_geotile = findfirst(geotiles[surface_mask].id .== single_geotile_test)
            geotiles[surface_mask] = geotiles[surface_mask][(geotiles[surface_mask].group .== geotiles[surface_mask][index_geotile, :group]), :]
        end
    end


    # align dv_altim and dv_gemb
    begin
        synthesized_gemb_fit = replace(path2runs[1], ".jld2" => "_gemb_fit.arrow")
        dh = FileIO.load(first(path2runs), "dh_hyps")
        ddate = dims(dh, :date)
        ddate_gemb = dims(dv_gemb, :date)
        index = .!isnan.(dv_gemb[geotile=1, pscale=1, ΔT=1])
        ex_gemb = extrema(ddate_gemb[index])
        index = .!isnan.(dh[geotile=1, height=1])
        ex_dv = extrema(ddate[vec(index)])
        ex = max(ex_gemb[1], ex_dv[1]), min(ex_gemb[2], ex_dv[2])
        index_dv = (ddate .>= ex[1]) .& (ddate .<= ex[2])
        index_gemb = (ddate_gemb .>= ex[1]) .& (ddate_gemb .<= ex[2])
        sum(index_dv) == sum(index_gemb) ? nothing : error("index_dv and index_gemb do not have the same length: $(ex)")

        dv_gemb = dv_gemb[date=At(ddate_gemb[index_gemb].val)]
    end

    # I'm getting an "geotile volume change contains all NaNs" when run using Threads.. no clue why
    @showprogress desc = "Calibrating GEMB model to altimetry data" for binned_synthesized_file in path2runs

        synthesized_gemb_fit = replace(binned_synthesized_file, ".jld2" => "_gemb_fit.arrow")

        if isfile(synthesized_gemb_fit) && (isnothing(force_remake_before) || Dates.unix2datetime(mtime(synthesized_gemb_fit)) > force_remake_before) && isnothing(single_geotile_test)
            printstyled("    -> Skipping $(synthesized_gemb_fit) because it was created after force_remake_before:$force_remake_before\n"; color=:light_green)
            continue
        else

            run_parameters = binned_filled_fileparts(synthesized_gemb_fit)
            surface_mask = run_parameters.surface_mask

            dh = FileIO.load(binned_synthesized_file, "dh_hyps")

            # Convert elevation change to volume change
            dv_altim = dh2dv_geotile(dh, area_km2[surface_mask])
            dv_altim = dv_altim[date=At(ddate[index_dv].val)]

            all_nans = dropdims(all(isnan.(dv_altim), dims=:date), dims=:date)

            if any(all_nans)
                error("geotile volume change contains all NaNs: $(binned_synthesized_file)")
            end

            # Find optimal fit to GEMB data
            # There are issues with calibrating the SMB model to individual geotiles since glaciers 
            # can cross multiple geotiles, therefore we calibrate the model for groups of 
            # distinct geotiles.

            if isnothing(single_geotile_test)
                df = gemb_bestfit_grouped(dv_altim, dv_gemb, geotiles[surface_mask]; seasonality_weight, distance_from_origin_penalty, ΔT_to_pscale_weight, calibrate_to, is_scaling_factor)

                df[!, :area_km2] = sum.(geotiles[surface_mask].area_km2)
                df.geometry = extent2rectangle.(df.extent)
                df = df[:, Not(:extent)]
                GeoDataFrames.write(synthesized_gemb_fit, df)
            else
                cost = gemb_bestfit_grouped(dv_altim, dv_gemb, geotiles[surface_mask]; seasonality_weight, distance_from_origin_penalty, ΔT_to_pscale_weight, calibrate_to, is_scaling_factor)

                return (cost, dv_altim)
            end
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

    geotile_total_area = sum(geotile_row.area_km2)
    if geotile_total_area == 0
        for k in vars
            gemb0[k][:] .= 0.
        end
        return gemb0
    end

    cnt = 0
    while (sum(geotile_row.area_km2[has_data]) / geotile_total_area) < min_gemb_coverage
        in_geotile = vec([within(geotile_extent, x, y) for (x, y) in zip(gemb["longitude"], gemb["latitude"])])
        h = gemb["height_effective"][in_geotile]

        for i in 1:(length(height_range)-1)
            has_data[i] = any((h .>= height_range[i]) .& (h .< height_range[i+1]))
        end

        geotile_extent = Extents.buffer(geotile_extent, (X=longitude_distance, Y=latitude_distance))
        cnt += 1

        if cnt > 100
            error("Buffering $(geotile_row.id) extent by $(round(Int,geotile_buffer*cnt/1000)) km, to achieve minimum GEMB coverage of $min_gemb_coverage, has taken too long. Check the code (maybe check coverage of gemb data) for potential issues.")
        end
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
    gemb_fill_gaps_and_add_ΔTsses(gemb, dΔT

Fill gaps in GEMB geotile data by interpolating and extrapolating missing values, then generate new classes for elevation change (ΔT shifting the filled data accordingly.

# Arguments
- `gemb`: Dictionary of GEMB output arrays, indexed by variable name.
- `dΔT`: Range or dimension object specifying ΔT classes to generate.

# Returns
- `Dict` containing filled and ΔT-shifted GEMB variables, with physical constraints enforced and cumulative variables reconstructed.

# Notes
- "smb" and "runoff" are excluded from filling and are reconstructed after processing to maintain mass conservation.
- Filling is performed on rates (not cumulative values) for physical consistency.
- Elevation classes are created by shifting the filled data vertically according to ΔT.
"""
function gemb_fill_gaps_and_add_ΔT_classes(gemb, dΔT; modify_melt_only=true)
    # --------------------------------------------------------------------
    # Fill missing data in GEMB outputs, then create ΔT-shifted classes
    # --------------------------------------------------------------------
    # Exclude SMB and runoff from gap-filling to preserve mass conservation

    vars = setdiff(collect(keys(gemb)), ["nobs", "smb", "runoff"])

    # Grab dimension info from one variable (assume all share dims)
    k = first(keys(gemb))
    ddate = dims(gemb[k], :date)
    dheight = dims(gemb[k], :height)
    dpscale = dims(gemb[k], :pscale)
    x = collect(dheight)  # height coordinates (for fitting)
    npts_linear_extrapolation = 7  # number of points for linear edge extrapolation

    if !modify_melt_only
        height_bin_interval = dΔT.val.step.hi  # spacing between ΔTsses
    end

    # Fill missing values for each variable and precipitation scaling scenario
    Threads.@threads for k in vars
        for pscale in dpscale
            # Interpolate/extrapolate on rates (not cumulative values)
            var0 = gemb[k][:, :, At(pscale)]
            dheight = dims(var0, :height)
            # Normalize x coordinate for fitting (center at 0, range -0.5:0.5)
            x0 = (1:length(dheight)) ./ length(dheight) .- 0.5
            M0 = hcat(ones(length(x0)), x0)
            npts = npts_linear_extrapolation

            # Skip if variable is zero everywhere
            if all(var0 .== 0)
                continue
            end

            # Compute rates (except for fac, which does not have meaningful rates)
            if k == "fac"
                var0_rate = var0[2:end, :]
            else
                # Find the first valid (not all-NaN) "time" slice for each height
                ind = findfirst(.!vec(all(isnan.(var0), dims=:height))) - 1
                var0[ind, :] .= 0  # ensure initial value is set
                var0_rate = var0[2:end, :] .- collect(var0[1:end-1, :])
            end

            # Fill missing data for each time slice ("date")
            for date in dims(var0_rate, :date)
                y = var0_rate[date=At(date)]  # 1D vector over heights for this date
                valid = .!isnan.(collect(y))  # Boolean array: where data is present

                # Skip if all values are valid or all are NaN
                if all(valid) || !any(valid)
                    continue
                end

                # For certain variables, treat only positive values as valid
                if k in ("acc", "refreeze", "melt", "rain", "fac")
                    valid = valid .& (y .>= 0)
                end

                # If all all missing, set to zero
                if !any(valid)
                    y[:] .= 0
                else
                    # Identify gaps surrounded by valid values (for loess interpolation)
                    validgap = validgaps(valid)

                    # If enough valid points, use LOESS smoothing/interpolation
                    if sum(valid) > 4
                        model = loess(x[valid], y[valid], span=0.3)
                        y[valid.|validgap] = predict(model, x[valid.|validgap])
                    end

                    # Fill any remaining NaNs (may be at edges) using linear fits
                    fill_index = isnan.(collect(y))  # locations still missing

                    # Choose which points are usable for edge regression
                    if k == "ec"
                        valid_interp_index = .!fill_index
                    else
                        valid_interp_index = .!fill_index .& (y .!= 0)
                    end

                    # Fill edges with regression if enough "good" points, else just set to zero
                    if any(fill_index) && (sum(valid_interp_index) > 3)
                        if sum(valid_interp_index) < npts + 2
                            # Not enough for two-sided regression, use all
                            M = @view M0[valid_interp_index, :]
                            param1 = M \ y[valid_interp_index]  # linear regression
                            y[fill_index] = param1[1] .+ x0[fill_index] .* param1[2]
                        else
                            # Do separate regression for lower/upper edge NaNs
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

                        # For accumulation-type variables: ensure non-negative
                        if k in ("acc", "refreeze", "melt", "rain", "fac")
                            y[y.<0] .= 0
                        end
                    else
                        # For hard-to-fill cases, just fill with 0 (conservative)
                        y[fill_index] .= 0
                    end
                end

                if any(isnan.(y))
                    error("NANs in y")
                end

                # Insert back into main array (after filling/imputation)
                gemb[k][At(date), :, At(pscale)] = y
            end
        end
    end

    # Apply physical/mass balance constraints directly to the filled arrays
    gemb = gemb_rate_physical_constraints!(gemb)
    gemb["fac"][gemb["fac"].<0] .= 0

    # Deepcopy filled rates to temporary store before ΔT-class construction
    gemb0 = Dict()
    for k in vars
        gemb0[k] = deepcopy(gemb[k][:, :, :])
    end

    # Allocate output array: will hold all ΔT-shifted classes
    gemb_mie = Dict()
    for k in vars
        gemb_mie[k] = fill(NaN, (ddate, dheight, dpscale, dΔT))
    end

    # this is a hack for testing
    if modify_melt_only
        for k in vars
            for ΔT in dΔT
                for pscale in dpscale
                    if k == "melt"
                        gemb_mie[k][:, :, At(pscale), At(ΔT)] = gemb0[k][:, :, At(pscale)] * ΔT
                    else
                        # Just copy across classes all other variables
                        gemb_mie[k][:, :, At(pscale), At(ΔT)] = gemb0[k][:, :, At(pscale)]
                    end
                end
            end
        end
    else
        # For each variable, for each ΔT class and pscale, create shifted arrays
        Threads.@threads for k in vars
            # Exclude zeros in edge extrapolation for everything except "ec"
            exclude_zeros_in_extrapolation = (k != "ec")
            for ΔT in dΔT
                for pscale in dpscale
                    if k == "nobs"
                        # For "nobs", just copy across classes (no height shifting)
                        gemb_mie[k][:, :, At(pscale), At(ΔT)] = gemb0[k][:, :, At(pscale)]
                    else
                        height_shift = round(Int16(ΔTeight_bin_interval))
                        gemb_mie[k][:, :, At(pscale), At(ΔT)] =
                            _matrix_shift_ud!(deepcopy(gemb0[k][:, :, At(pscale)]), height_shift;
                                exclude_zeros_in_extrapolation, npts_linear_extrapolation)
                    end
                end
            end
        end
    end

    # Enforce physical/mass conservation constraints again (this time to ΔTemble)
    gemb_mie = gemb_rate_physical_constraints!(gemb_mie)
    gemb_mie["fac"][gemb_mie["fac"].<0] .= 0

    # Reconstruct cumulative variables ("acc", "melt", etc) by time-integrating the rates
    for k in setdiff(collect(keys(gemb_mie)), ["nobs", "fac"])
        sindex = findfirst(.!isnan.(gemb_mie[k][:, 1, 1, 1]))
        gemb_mie[k][sindex:end, :, :, :] = cumsum(gemb_mie[k][sindex:end, :, :, :], dims=:date)
    end

    # Output dictionary: each variable [date, height, pscale, ΔTss]
    return gemb_mie
end

"""
    read_gemb_files(
        gemb_files,
        gembinfo;
        vars2extract = ["acc", "refreeze", "melt", "rain", "ec", "fac", "smb", "dv"],
        date_range,
        date_center,
        path2land_sea_mask = nothing,
        minimum_land_coverage_fraction = 0.70,
    )

Read and combine data from multiple GEMB output files (.mat), extracting requested variables and metadata, and standardizing into a single dictionary suitable for geotile binning.

All mass variables are returned in units of meters ice equivalent, except "fac" which is in meters of air.

# Arguments
- `gemb_files::Vector{String}`: Paths to GEMB output files to read and merge.
- `gembinfo`: Structure or dictionary containing precipitation scale (`precipitation_scale`) and elevation delta (`elevation_delta`) information to assign to outputs.
- `vars2extract::Vector{String}`: List of variables to extract from each file (default includes main GEMB outputs).
- `date_range`: Either two-element tuple or array, used to convert times to decimal years and select matching slices.
- `date_center`: Central value for time dimension for all files (must be provided). !!! Needed to ensure exact date matching with other data sources !!!
- `path2land_sea_mask::Union{String,Nothing}`: Path to land-sea mask raster file to screen data (optional).
- `minimum_land_coverage_fraction::Float64`: Minimum fraction of land required for a point to be included (default 0.70).

# Returns
- `Dict{String, Any}`: Dictionary containing merged variables. Includes:
    - All selected variables as arrays concatenated across all files
    - Coordinates: "latitude", "longitude", "height", "date"
    - "precipitation_scale", "elevation_delta": Array of assigned class values per data point
    - "height_effective": height plus elevation delta (for use in further gridding).
    - If land mask is used: "land_fraction"
"""
function read_gemb_files(gemb_files, gembinfo; vars2extract=["acc", "refreeze", "melt", "rain", "ec", "fac", "smb", "dv"], date_range, date_center, path2land_sea_mask=nothing, minimum_land_coverage_fraction=0.70)
    datebin_edges = decimalyear.(date_range)

    gemb = Vector{Any}(undef, length(gemb_files))
    pscale = gembinfo.precipitation_scale
    Δelevation = gembinfo.elevation_delta
    @showprogress desc = "Reading raw GEMB output into memory, this will take ~2 min [peak memory usage: 150GB]" Threads.@threads for i in eachindex(gemb_files)
        gemb_file = gemb_files[i]
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

        gemb0["precipitation_scale"] = pscale[At(precip_scale_ind)]
        gemb0["elevation_delta"] = Δelevation[At(elev_delta_ind)]

        gemb[i] = gemb0
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
    latitude = vec(reduce(vcat, getindex.(gemb, Ref("latitude"))))
    longitude = vec(reduce(vcat, getindex.(gemb, Ref("longitude"))))
    height0 = vec(reduce(vcat, getindex.(gemb, Ref("height"))))

    # Ensure longitude is in -180 to 180 range
    index = (longitude .< -180)
    longitude[index] .+= 360
    index = (longitude .> 180)
    longitude[index] .-= 360

    precipitation_scale0 = Float64[]
    elevation_delta0 = Float64[]
    for g in gemb
        precipitation_scale0 = vcat(precipitation_scale0, fill(g["precipitation_scale"], length(g["latitude"])))
        elevation_delta0 = vcat(elevation_delta0, fill(g["elevation_delta"], length(g["latitude"])))
    end

    gemb1 = Dict(
        "date" => date_center,
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
    gemb = merge(foo, gemb1)

    # Add effective height (height + elevation_delta)
    gemb["height_effective"] = gemb["height"] .+ gemb["elevation_delta"]

    # add land fraction to gemb
    if !isnothing(path2land_sea_mask)
        land_sea_mask = Raster(path2land_sea_mask)[:, :, 1]

        xpt = copy(gemb["longitude"])
        index = xpt .< 0
        xpt[index] .+= 360
        ypt = gemb["latitude"]

        pt = tuple.(X.(Near.(xpt)), Y.(Near.(ypt)))
        gemb["land_fraction"] = land_sea_mask[pt].data

        index_valid = gemb["land_fraction"] .>= minimum_land_coverage_fraction
        n = size(gemb["land_fraction"], 1)
        for k in keys(gemb)
            if size(gemb[k], 1) == n
                if ndims(gemb[k]) == 1
                    gemb[k] = gemb[k][index_valid]
                else
                    gemb[k] = gemb[k][index_valid, :]
                end
            end
        end

    end
    return gemb
end


"""
    dv_gemb(gemb, area_km2, dpscale_new)

Compute volume change (ΔV) for GEMB output, interpolated to new precipitation scaling factors.

# Arguments
- `gemb`: Dictionary of GEMB output arrays, each keyed by variable name.
- `area_km2`: Array of grid cell areas in km².
- `dpscale_new`: Array of new precipitation scaling factors to interpolate to.

# Returns
- `Dict` mapping each variable to a 3D array of volume change (date × pscale × ΔT
"""
function dv_gemb(gemb, area_km2, dpscale_new)
    k = first(keys(gemb))
    dpscale = dims(gemb[k], :pscale)
    dΔT = dims(gemb[k], :ΔT)
    ddate = dims(gemb[k], :date)

    # Ensure all old pscales are present in new pscales
    if !isempty(setdiff(dpscale, dpscale_new))
        error("old pscales [$(dpscale)] must exist in new pscales [$(dpscale_new)]")
    end

    # Initialize output dictionary
    dv_gemb = Dict()
    for k in keys(gemb)
        dv_gemb[k] = fill(NaN, (ddate, dpscale_new, dΔT))
    end

    Threads.@threads for k in collect(keys(gemb))
        valid = .!isnan.(gemb[k])
        (date_range, height_range, pscale_range, ΔT_range) = validrange(valid)
        for date in ddate[date_range]
            for ΔT in dΔT[ΔT_range]
                dv = @d dropdims(sum(gemb[k][date=At(date), ΔT=At(ΔT)][height_range, :] .* area_km2[height_range] ./ 1000, dims=:height), dims=:height)
                dv_interp = DataInterpolations.LinearInterpolation(dv, val(dpscale))
                dv_gemb[k][date=At(date), ΔT=At(ΔT)] = dv_interp(val(dpscale_new))
            end
        end
    end
    return dv_gemb
end


"""
    gemb_ensemble_dv(; gemb_run_id, vars2load=["runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze"])

Load GEMB geotile ensemble data for a given run ID. Returns a dictionary of variables specified in `vars2load` from the corresponding GEMB geotile volume change file. If `vars2load` is `nothing`, loads all variables from the file.

# Arguments
- `gemb_run_id`: Identifier for the GEMB run.
- `dpscale_expanded_increment`: Step size for increasing the resolution or range of `pscale`.
- `dΔT_expanded_increment`: Step size for increasing the resolution or range of `ΔT`.
- `vars2load`: Array of variable names to load (default: common mass balance variables). If `nothing`, loads all variables.

# Returns
- `Dict` containing the requested variables.
"""
function gemb_ensemble_dv(; gemb_run_id=4)
    gembinfo = gemb_info(; gemb_run_id)
    gemb_geotile_filename_dv = replace(gembinfo.filename_gemb_combined, ".jld2" => "_geotile_dv.jld2")
    dv_gemb = FileIO.load(gemb_geotile_filename_dv, "gemb_dv")
    
    gemb_add_derived_vars!(dv_gemb)

    return dv_gemb

end

function gemb_dv_interpolater(geotile, pscale, ΔT, gemb_dv)
    @warn "this is really slow due to large number of allocations"
    ddate = dims(gemb_dv, :date)
    dgeotile = dims(gemb_dv, :geotile)
    dpscale = dims(gemb_dv, :pscale)
    dΔT = dims(gemb_dv, :ΔT)

    itp = Dict()
    for k in keys(gemb_dv)
        itp[k] = Interpolations.interpolate((eachindex(dgeotile), eachindex(ddate), dpscale.val, dΔT.val), gemb_dv[k].data, (NoInterp(), NoInterp(), Gridded(Linear()), Gridded(Linear())))
    end

    gemb_dv_out = copy(gemb_dv[geotile = 1, pscale = 1, ΔT = 1])
    geotile_index = findfirst(dgeotile.val .== geotile)

    for k in keys(gemb_dv)
        gemb_dv_out[k][:] = itp[k].(Ref(geotile_index), eachindex(ddate), Ref(pscale), Ref(ΔT))
    end

    return gemb_dv_out
end

function gemb_add_derived_vars!(dv_gemb)
    dv_gemb = merge(dv_gemb, (; smb=dv_gemb[:acc] .- dv_gemb[:melt] .+ dv_gemb[:refreeze] .- dv_gemb[:ec], runoff=dv_gemb[:melt] .- dv_gemb[:refreeze]))
    dv_gemb = merge(dv_gemb, (; dv=dv_gemb[:smb] .+ dv_gemb[:fac]))

    return dv_gemb
end

function gemb_dv_sample(pscale, ΔT, dv_gemb)

    ddate = dims(dv_gemb, :date)
    dpscale = dims(dv_gemb, :pscale)
    dΔT = dims(dv_gemb, :ΔT)

    itp = Interpolations.interpolate((eachindex(ddate), dpscale.val, dΔT.val), dv_gemb.data, (NoInterp(), Gridded(Linear()), Gridded(Linear())))
    gemb_dv_out = DimArray(itp.(eachindex(ddate), Ref(pscale), Ref(ΔT)),ddate)

    return gemb_dv_out  
end

function gemb_dv_sample!(gemb_dv_out::DimArray, pscale, ΔT, dv_gemb::DimArray)

    ddate = dims(dv_gemb, :date)
    dpscale = dims(dv_gemb, :pscale)
    dΔT = dims(dv_gemb, :ΔT)

    itp = Interpolations.interpolate((eachindex(ddate), dpscale.val, dΔT.val), dv_gemb.data, (NoInterp(), Gridded(Linear()), Gridded(Linear())))
    gemb_dv_out[:] .= itp.(eachindex(ddate), Ref(pscale), Ref(ΔT))

    return gemb_dv_out
end

function gemb_dv_sample!(gemb_dv_out, pscale_index, ΔT_index, dv_gemb)

    itp = Interpolations.interpolate((axes(dv_gemb, 1), axes(dv_gemb, 2), axes(dv_gemb, 3)), dv_gemb, (NoInterp(), Gridded(Linear()), Gridded(Linear())))
    gemb_dv_out[:] .= itp.(axes(dv_gemb, 1), pscale_index, ΔT_index)

    return gemb_dv_out
end





"""
    gemb_dv_interpolate(dv_gemb, dpscale_expanded_increment, dΔT_expanded_increment)

Interpolate a GEMB DimStack of volume change (`dv_gemb`) across extended precipitation scaling (`pscale`) and temperature perturbation (`ΔT`) classes.

# Arguments
- `dv_gemb`: DimStack or dictionary containing GEMB volume change data, with dimensions including :geotile, :date, :pscale, and :ΔT.
- `dpscale_expanded_increment`: Step size for increasing the resolution or range of `pscale`.
- `dΔT_expanded_increment`: Step size for increasing the resolution or range of `ΔT`.

# Returns
- `gemb_dv_new`: Interpolated `DimStack` structure with expanded/adjusted :pscale and :ΔT dimensions.

# Description
- This function first determines the expanded/desired ranges for `pscale` and `ΔT` using the provided increments.
- It then initializes a new DimStack with NaN-filled arrays matching the desired output dimensions.
- For each variable in `dv_gemb`, a linear gridded interpolation is constructed across its original (`pscale`, `ΔT`) grid for each geotile and date, and evaluated at the new grid, storing the results in `gemb_dv_new`.

# Example
```julia
gemb_dv_new = gemb_dv_interpolate(dv_gemb; dpscale_expanded_increment=0.25, dΔT_expanded_increment=0.5)
```
"""
function gemb_dv_interpolate(dv_gemb; dpscale_expanded_increment=0.25, dΔT_expanded_increment=0.5)

    dgeotile = dims(dv_gemb, :geotile)
    ddate = dims(dv_gemb, :date)
    dpscale = dims(dv_gemb, :pscale)
    dΔT = dims(dv_gemb, :ΔT)

    pscale_start = minimum(dpscale)
    pscale_start < 1 ? pscale_start = -1 / pscale_start : pscale_start
    pscale_start = ceil(pscale_start, digits=0)
    pscale_end = maximum(dpscale)
    pscale_end < 1 ? pscale_end = -1 / pscale_end : pscale_end
    pscale_end = floor(pscale_end, digits=0)
    pscale_new = pscale_start:dpscale_expanded_increment:pscale_end
    pscale_new = pscale_new[(1 .<= pscale_new).|(pscale_new.<-1)]
    pscale_new[pscale_new.<0] = -1 ./ pscale_new[pscale_new.<1]
    dpscale_new = Dim{:pscale}(pscale_new)

    ΔT_start = minimum(dΔT)
    ΔT_start = ceil(ΔT_start, digits=0)
    ΔT_end = maximum(dΔT)
    ΔT_end = floor(ΔT_end, digits=0)
    ΔT_new = ΔT_start:dΔT_expanded_increment:ΔT_end
    dΔT_new = Dim{:ΔT}(ΔT_new)

    gemb_dv_new = DimStack([DimArray(fill(NaN, (dgeotile, ddate, dpscale_new, dΔT_new); name=Symbol(k))) for k in keys(dv_gemb)]...)

    index_date_valid = dropdims(any(.!isnan.(dv_gemb[first(keys(dv_gemb))][pscale=At(1), ΔT=At(0)]), dims=:geotile), dims=:geotile)
    index_date_range, = validrange(index_date_valid)
    daterange = ddate[first(index_date_range)] .. ddate[last(index_date_range)]

    if false
        k = :melt
        begin
            geotile_id = geotiles_golden_test[1]
            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], title="$geotile_id: $(k)")
            cmap_bar = Makie.resample_cmap(:thermal, (length(dpscale)))
            cnt = 1
            for pscale in dpscale
                lines!(dv_gemb[k][date=daterange, geotile=At(geotile_id), pscale=At(pscale), ΔT=At(0)]; color=cmap_bar[cnt], label="pscale: $(pscale)")
                cnt += 1
            end
            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)

            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], title="$geotile_id: $(k)")
            cmap_bar = Makie.resample_cmap(:thermal, (length(dΔT)))
            cnt = 1


            for ΔT in dΔT
                lines!(dv_gemb[k][date=daterange, geotile=At(geotile_id), pscale=At(1), ΔT=At(ΔT)]; color=cmap_bar[cnt], label="ΔT: $(ΔT)")
                cnt += 1
            end

            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)

            dates2plot = [DateTime(2000, 5, 1), DateTime(2010, 12, 30), DateTime(2015, 9, 25)]
            geotile_id = geotiles_golden_test[1]
            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], title="$geotile_id: $(k) as a function of pscale")
            cmap_bar = Makie.resample_cmap(:thermal, (length(dpscale)))
            cnt = 1
            for date in dates2plot
                lines!(dv_gemb[k][date=Near(date), geotile=At(geotile_id), ΔT=At(0)]; color=cmap_bar[cnt], label="date: $(date)")
                cnt += 1
            end
            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)

            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], title="$geotile_id: $(k) as a function of ΔT")
            cmap_bar = Makie.resample_cmap(:thermal, (length(dpscale)))
            cnt = 1
            for date in dates2plot
                lines!(dv_gemb[k][date=Near(date), geotile=At(geotile_id), pscale=At(1)]; color=cmap_bar[cnt], label="date: $(date)")
                cnt += 1
            end
            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)
        end
    end

    # linearly interpolate between pscale and ΔT... maybe this should be done in fitting function?
    Threads.@threads for k in keys(dv_gemb)
        for date in index_date_range
            for gtidx = eachindex(dgeotile)
                itp = Interpolations.interpolate((dpscale.val, dΔT.val), dv_gemb[k][date=date, geotile=gtidx].data, Gridded(Linear()))
                gemb_dv_new[k][geotile=gtidx, date=date] = itp(dpscale_new.val, dΔT_new.val)
            end
        end
    end

    if false
        k = :melt
        begin
            geotile_id = geotiles_golden_test[1]
            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], title="$geotile_id: $(k)")
            cmap_bar = Makie.resample_cmap(:thermal, (length(dpscale_new)))
            cnt = 1
            for pscale in dpscale_new
                lines!(gemb_dv_new[k][date=daterange, geotile=At(geotile_id), pscale=At(pscale), ΔT=At(0)]; color=cmap_bar[cnt], label="pscale: $(pscale)")
                cnt += 1
            end
            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)

            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], title="$geotile_id: $(k)")
            cmap_bar = Makie.resample_cmap(:thermal, (length(dΔT_new)))
            cnt = 1


            for ΔT in dΔT_new
                lines!(gemb_dv_new[k][date=daterange, geotile=At(geotile_id), pscale=At(1), ΔT=At(ΔT)]; color=cmap_bar[cnt], label="ΔT: $(ΔT)")
                cnt += 1
            end
            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)

            dates2plot = [DateTime(2000, 5, 1), DateTime(2010, 12, 30), DateTime(2015, 9, 25)]
            geotile_id = geotiles_golden_test[1]
            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], title="$geotile_id: $(k) as a function of pscale")
            cmap_bar = Makie.resample_cmap(:thermal, (length(dpscale)))
            cnt = 1
            for date in dates2plot
                lines!(gemb_dv_new[k][date=Near(date), geotile=At(geotile_id), ΔT=At(0)]; color=cmap_bar[cnt], label="date: $(date)")
                cnt += 1
            end
            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)

            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], title="$geotile_id: $(k) as a function of ΔT")
            cmap_bar = Makie.resample_cmap(:thermal, (length(dpscale)))
            cnt = 1
            for date in dates2plot
                lines!(gemb_dv_new[k][date=Near(date), geotile=At(geotile_id), pscale=At(1)]; color=cmap_bar[cnt], label="date: $(date)")
                cnt += 1
            end
            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)
        end
    end

    return gemb_dv_new
end


"""
    process_gemb_geotiles(gemb, geotiles; maximum_search_iterations=20, height_bins_extrapolate=1, maximum_extrap_fraction=0.1, search_buffer_increment=100_000, show_interp_extrap_plots=false, show_interp_extrap_stats=false)

Process GEMB data and organize it into geotiles with elevation and precipitation classes.

This function processes GEMB model output by:
1. Organizing data into geotiles with consistent spatial and temporal dimensions
2. Applying spatial buffering to ensure minimum coverage requirements
3. Filling data gaps through interpolation and extrapolation across elevation bands
4. Calculating volume changes for each geotile

# Arguments
- `gemb::Dict`: Dictionary containing GEMB model output with keys including variables and metadata
- `geotiles::DataFrame`: DataFrame with geotile definitions, hypsometry, and extent information
- `maximum_search_iterations::Int`: Maximum number of search iterations to populate elevation profile (default: 20)
- `height_bins_extrapolate::Int`: Number of bins above and below valid height range to extrapolate (default: 1)
- `maximum_extrap_fraction::Float64`: Maximum fraction of area below/above valid height range for extrapolation (default: 0.1)
- `search_buffer_increment::Int`: Increment of search buffer in meters with each search iteration (default: 100_000)
- `show_interp_extrap_plots::Bool`: Whether to show interpolation/extrapolation plots (default: false)
- `show_interp_extrap_stats::Bool`: Whether to show interpolation/extrapolation statistics (default: false)


# Returns
- `dv_gemb::DimStack`: DimStack containing processed GEMB data with dimensions (:geotile, :date, :pscale, :ΔT)
"""
function process_gemb_geotiles(
    gemb,
    geotiles;
    maximum_search_iterations=20,
    height_bins_extrapolate=1,
    maximum_extrap_fraction=0.1,
    search_buffer_increment=100_000,
    show_interp_extrap_plots=false,
    show_interp_extrap_stats=false,
    single_geotile_test=nothing,
    elevation_classes_method=:none,
)

    unique_precipitation_scale = unique(gemb["precipitation_scale"])
    unique_elevation_delta = unique(gemb["elevation_delta"])
    vars = setdiff(collect(keys(gemb)), ("precipitation_scale", "elevation_delta", "latitude", "longitude", "height", "date", "height_effective", "land_fraction"))

    height_range, height_center = project_height_bins()
    ddate = Dim{:date}(gemb["date"])
    dheight = Dim{:height}(height_center)
    dpscale = Dim{:pscale}(sort(unique_precipitation_scale))
    dΔheight = Dim{:Δheight}(sort(unique_elevation_delta))
    dgeotile = Dim{:geotile}(geotiles.id)

    index_date_valid = .!isnan.(gemb[vars[1]][1, :])
    index_date_range, = validrange(index_date_valid)
    daterange = ddate[first(index_date_range)] .. ddate[last(index_date_range)]

    cmap_bar = Makie.resample_cmap(:dense, (maximum_search_iterations * 3 + 1))

    if !isnothing(single_geotile_test)
        geotiles = geotiles[geotiles.id.==single_geotile_test, :]
        show_interp_extrap_plots = true
    end

    if elevation_classes_method == :none
        gemb_dv0 = DimStack([DimArray(fill(NaN, (dgeotile, ddate, dpscale, dΔheight)); name=k) for k in vars]...)
        gemb_dv0 = DimStack([DimArray(fill(NaN, (dgeotile, ddate, dpscale, dΔheight)); name=k) for k in vars]...)
    else
        ΔT_range, ΔT_center = project_ΔT_bins()
        dΔT = Dim{:ΔT}(ΔT_center)
        dh1 = dheight.val.step.hi
        max_height_offset = ceil(abs((1000 / -6.5) * maximum(abs.(dΔT.val))) / dh1) * abs(dh1)

        height_center1 = (minimum(dheight.val)+minimum(dΔheight.val)):dh1:((maximum(dheight.val)+maximum(dΔheight.val)))
        height_range1 = (minimum(height_center1)-dh1/2):dh1:(maximum(height_center1)+dh1/2)
        dheight1 = Dim{:height}(height_center1)

        if elevation_classes_method == :Δelevation
            gemb_dv0 = DimStack([DimArray(fill(NaN, (dgeotile, ddate, dpscale, dΔT)); name=k) for k in vars]...)
        elseif elevation_classes_method == :mscale
            mscale_range, mscale_center = project_mscale_bins()
            dmscale = Dim{:ΔT}(mscale_center)
            gemb_dv0 = DimStack([DimArray(fill(NaN, (dgeotile, ddate, dpscale, dmscale)); name=k) for k in vars]...)
        end
    end

    @showprogress desc = "Populating geotiles with GEMB data, this will take 2 to 10 min on 128 threads, depending on method" Threads.@threads for geotile_row in eachrow(geotiles)
    #@warn "threads turned off"
    #for geotile_row in eachrow(geotiles)

        geotile_hyps_area_km2 = DimArray(geotile_row.area_km2, dheight; name="area_km2")
        geotile_total_area = sum(geotile_hyps_area_km2)

        if geotile_total_area == 0
            for k in Symbol.(vars)
                gemb_dv0[k][geotile=At(geotile_row.id), date=daterange] .= 0.
            end
            continue
        end

        index_has_ice = geotile_hyps_area_km2 .> 0
        index_has_ice_range, = validrange(index_has_ice)

        search_extent = geotile_row.extent
        cnt = 0
        bar = 1
        search_buffer = 0
        coverage_height_fraction = 0.
        points_in_geotile_count = 0

        if show_interp_extrap_plots
            f = _publication_figure(; columns=1, rows=2)
            ax = CairoMakie.Axis(f[1, 1], yaxisposition=:right, title="$(geotile_row.id)")
            hidedecorations!(ax)  # hides ticks, grid and labels
            hidespines!(ax)  # hide the frame
        end

        # Populate classes through interpolation and extrapolation... elevation classes are not used to pupulate elevation range
        if elevation_classes_method == :none

            gemb0 = DimStack([DimArray(fill(NaN, (ddate, dheight, dpscale, dΔheight)); name=k) for k in vars]...)

            while isnan(coverage_height_fraction) || (round(coverage_height_fraction, digits=3) < 1.0)
                cnt += 1
                if cnt > 1
                    search_buffer += search_buffer_increment
                    latitude_distance, longitude_distance = meters2lonlat_distance(search_buffer, mean(search_extent.Y))
                    search_extent = Extents.buffer(search_extent, (X=longitude_distance, Y=latitude_distance))
                end

                # points in geotiles 
                index_geotile = within.(Ref(search_extent), gemb["longitude"], gemb["latitude"])

                if !any(index_geotile)
                    bar += 3
                    continue
                end

                index_geotile = findall(index_geotile)

                height0 = gemb["height"][index_geotile]
                elevation_delta0 = gemb["elevation_delta"][index_geotile]
                precipitation_scale0 = gemb["precipitation_scale"][index_geotile]

                index_existing_heights = .!isnan.(gemb0[:fac][pscale=At(1), Δheight=At(0)])
                index_existing_heights = any(index_existing_heights, dims=:date)[:]

                # loop is not thread safe
                for i in eachindex(height_center)[.!index_existing_heights]
                    index_height = height_range[i] .<= height0 .< height_range[i+1]

                    if !any(index_height)
                        continue
                    end

                    for Δheight in dΔheight
                        index_Δheight = Δheight .== elevation_delta0

                        for pscale in dpscale
                            index_precipitation = precipitation_scale0 .== pscale
                            index_valid = index_geotile[index_height.&index_Δheight.&index_precipitation]

                            for k in vars
                                gemb0[Symbol(k)][height=At(height_center[i]), pscale=At(pscale), Δheight=At(Δheight)] = mean(gemb[k][index_valid, :], dims=1)
                            end
                        end
                    end
                end
            
                valid_ensemble = .!isnan.(gemb0[Symbol(vars[1])][pscale=At(1), Δheight=At(0)])
                index_height_valid = any(valid_ensemble, dims=:date)[:]

                if any(index_height_valid)
                    # Ask questions about the ensemble coverage

                    # plot the ensemble coverage
                    if show_interp_extrap_plots
                        bar_index = (index_height_valid .& .!index_existing_heights)
                        bar_labels = repeat(["raw-$(cnt)"], sum(index_has_ice))
                        bar_labels[.!bar_index[index_has_ice]] .= ""
                        barplot!(ax, dheight[index_has_ice].val, bar_index[index_has_ice]; direction=:x, color=(cmap_bar[bar]), bar_labels=bar_labels, label_position=:center)
                        bar += 1
                    end

                    if false
                        for k in [:fac, :acc, :melt]
                            f = _publication_figure(; columns=1, rows=2)
                            ax = CairoMakie.Axis(f[1, 1]; palette=(:color => Makie.resample_cmap(:thermal, length(findall(index_height_valid)))), title="$(geotile_row.id): $(k)")

                            for height_index in findall(index_height_valid)
                                lines!(ax, gemb0[k][dates=daterange, height=At(dheight[height_index]), pscale=At(1), Δheight=At(0)]; label="$(dheight[height_index])")
                            end
                            
                            f[1, 2] = Legend(f, ax, "height [m]", framevisible=false)
                            display(f)
                        end
                    end

                    # Ask a bunch of questions to determine where extrapolation is appropriate
                    index_height_valid_extrapolate = dilate(index_height_valid, height_bins_extrapolate)
                    index_height_range, = validrange(index_height_valid_extrapolate)

                    if 1 < first(index_height_range) > first(index_has_ice_range)
                        area_below_height_range = sum(geotile_hyps_area_km2[1:(first(index_height_range)-1)])
                        coverage_below_height_range = area_below_height_range ./ geotile_total_area
                    else
                        area_below_height_range = 0
                        coverage_below_height_range = 0
                    end

                    if last(index_height_range) < last(index_has_ice_range)
                        area_above_height_range = sum(geotile_hyps_area_km2[(last(index_height_range)+1):end])
                        coverage_above_height_range = area_above_height_range ./ geotile_total_area
                    else
                        coverage_above_height_range = 0
                    end

                    if (0 < coverage_below_height_range < maximum_extrap_fraction)
                        index_height_range = first(index_has_ice_range):last(index_height_range)
                    end

                    if (0 < coverage_above_height_range < maximum_extrap_fraction)
                        index_height_range = first(index_height_range):last(index_has_ice_range)
                    end

                    index_extrapolate = falses(length(dheight))
                    index_extrapolate[index_height_range] .= true
                    index_interp_range, = validrange(index_height_valid)
                    index_extrapolate[index_interp_range] .= false

                    index_fill = falses(length(dheight))
                    index_fill[index_height_range] .= true
                    index_fill = (index_fill .& .!index_height_valid)

                    # limit size of interpolation gap
                    index_height_range_interpolate, = validrange(index_height_valid)
                    index_height_interpolate = falses(length(dheight))
                    index_height_interpolate[index_height_range_interpolate] .= true
                    index_height_interpolate = (index_height_interpolate .& .!index_height_valid)

                    if any(index_height_interpolate)
                        index_valid_interpolate = .!(true_block_size(index_fill) .> (height_bins_extrapolate * 2 + 1))
                        index_valid_interpolate = dilate(index_valid_interpolate, height_bins_extrapolate)
                        index_fill = index_fill .& index_valid_interpolate
                    end

                    if sum(index_height_valid) > 1
                        for k in Symbol.(vars)
                            itp = extrapolate(Interpolations.interpolate((1:length(ddate), dheight.val[index_height_valid], 1:length(dpscale), 1:length(dΔheight)), gemb0[k][height=findall(index_height_valid)], (NoInterp(), Gridded(Linear()), NoInterp(), NoInterp())), Flat())
                            for date in index_date_range
                                for Δheight in eachindex(dΔheight)
                                    for pscale in eachindex(dpscale)
                                        gemb0[k][date=date, height=index_fill, pscale=pscale, Δheight=Δheight] = itp.(Ref(date), dheight.val[index_fill], Ref(pscale), Ref(Δheight))
                                    end
                                end
                            end
                        end
                    else
                        for k in Symbol.(vars)
                            fill_values = gemb0[k][height=findall(index_height_valid)]

                            for height0 in dheight[index_fill]
                                gemb0[Symbol(k)][height=At(height0)] = fill_values
                            end
                        end
                    end

                    if show_interp_extrap_plots
                        bar_index = index_fill .& .!index_extrapolate
                        bar_labels = repeat(["interp-$(cnt)"], sum(index_has_ice))
                        bar_labels[.!bar_index[index_has_ice]] .= ""
                        barplot!(ax, dheight[index_has_ice].val, bar_index[index_has_ice]; direction=:x, color=(cmap_bar[bar]), bar_labels=bar_labels, label_position=:center)
                        bar += 1

                        bar_index = index_extrapolate
                        bar_labels = repeat(["extrap-$(cnt)"], sum(index_has_ice))
                        bar_labels[.!bar_index[index_has_ice]] .= ""
                        barplot!(ax, dheight[index_has_ice].val, bar_index[index_has_ice]; direction=:x, color=(cmap_bar[bar]), bar_labels=bar_labels, label_position=:center)
                        bar += 1
                    end

                    valid_ensemble = .!isnan.(gemb0[:fac][pscale=At(1), Δheight=At(0)])
                    index_height_valid = any(valid_ensemble, dims=:date)[:]
                    coverage_height_fraction = sum(geotile_hyps_area_km2[index_height_valid]) / geotile_total_area
                end

                if cnt > maximum_search_iterations
                    error("MAXIMUM SEARCH ITERATIONS REACHED: Buffering $(geotile_row.id) search extent by $(round(Int,search_buffer/1000)) km to populate elevation profile")
                end
            end

            # calculate volume change in unit of km3 [GEMB assumes an ice density of 910 kg/m3]
            for k in Symbol.(vars)
                dv = @d gemb0[k][date=daterange] .* geotile_hyps_area_km2 ./ 1000
                dv[isnan.(dv)] .= 0
                dv = sum(dv; dims=:height)
                gemb_dv0[k][geotile=At(geotile_row.id), date=daterange] = dropdims(dv, dims=:height)
            end
        else
            # ensure that height range covers desired ΔT classes
            gemb1 = DimStack([DimArray(fill(NaN, (ddate, dheight1, dpscale)); name=k) for k in vars]...)

            index_has_ice1 = (dheight1.val .>= (dheight[first(index_has_ice_range)] - max_height_offset)) .& (dheight1.val .<= (dheight[last(index_has_ice_range)] + max_height_offset))
            index_has_ice_range1, = validrange(index_has_ice1)

            # use simulated elevation classes to populate elevation range
            while isnan(coverage_height_fraction) || (round(coverage_height_fraction, digits=3) < 1.0)

                cnt += 1
                if cnt > 1
                    search_buffer += search_buffer_increment
                    latitude_distance, longitude_distance = meters2lonlat_distance(search_buffer, mean(search_extent.Y))
                    search_extent = Extents.buffer(search_extent, (X=longitude_distance, Y=latitude_distance))
                end

                # points in geotiles 
                index_geotile = within.(Ref(search_extent), gemb["longitude"], gemb["latitude"])
                points_in_geotile_count = sum(index_geotile)
                if points_in_geotile_count == 0
                    bar += 3
                    continue
                end

                index_geotile = findall(index_geotile)

                height0 = gemb["height_effective"][index_geotile]
                elevation_delta0 = gemb["elevation_delta"][index_geotile]
                precipitation_scale0 = gemb["precipitation_scale"][index_geotile]

                index_existing_heights = .!isnan.(gemb1[:fac][pscale=At(1)])
                index_existing_heights = any(index_existing_heights, dims=:date)[:]

                # loop is not thread safe
                for i in eachindex(height_center1)[.!index_existing_heights]
                    index_height = height_range1[i] .<= height0 .< height_range1[i+1]

                    if !any(index_height)
                        continue
                    end

                    for pscale in dpscale
                        index_precipitation = (precipitation_scale0 .== pscale)
                        index_valid = index_geotile[index_height .& index_precipitation]

                        for k in vars
                            gemb1[Symbol(k)][height=At(height_center1[i]), pscale=At(pscale)] = mean(gemb[k][index_valid, :], dims=1)
                        end
                    end
                end

                valid_ensemble = .!isnan.(gemb1[Symbol(vars[1])][pscale=At(1)])
                index_height_valid = any(valid_ensemble, dims=:date)[:]

                if sum(index_height_valid) != 0
                    # Ask questions about the ensemble coverage

                    # plot the ensemble coverage
                    if show_interp_extrap_plots
                        bar_index = (index_height_valid .& .!index_existing_heights)
                        bar_labels = repeat(["raw-$(cnt)"], sum(index_has_ice1))
                        bar_labels[.!bar_index[index_has_ice1]] .= ""
                        barplot!(ax, dheight1[index_has_ice1].val, bar_index[index_has_ice1]; direction=:x, color=(cmap_bar[bar]), bar_labels=bar_labels, label_position=:center)
                        bar += 1
                    end

                    if false
                        for k in [:fac, :acc, :melt]
                            f = _publication_figure(; columns=1, rows=2)
                            ax = CairoMakie.Axis(f[1, 1]; palette=(:color => Makie.resample_cmap(:thermal, length(findall(index_height_valid)))), title="$(geotile_row.id): $(k)")

                            for height_index in findall(index_height_valid)
                                lines!(ax, gemb1[k][dates=daterange, height=At(dheight1[height_index]), pscale=At(1)]; label="$(dheight1[height_index])")
                            end
                            f[1, 2] = Legend(f, ax, "height [m]", framevisible=false)
                            display(f)
                        end
                    end


                    # Ask a bunch of questions to determine where extrapolation is appropriate
                    index_height_valid_extrapolate = dilate(index_height_valid, height_bins_extrapolate)
                    index_height_range, = validrange(index_height_valid_extrapolate)

                    if 1 < first(index_height_range) > first(index_has_ice_range1)
                        area_below_height_range = sum(geotile_hyps_area_km2[height=-Inf .. (dheight1[first(index_height_range)-1] + max_height_offset)])

                        coverage_below_height_range = area_below_height_range ./ geotile_total_area
                    else
                        area_below_height_range = 0
                        coverage_below_height_range = 0
                    end

                    if last(index_height_range) < last(index_has_ice_range1)
                        area_above_height_range = sum(geotile_hyps_area_km2[height=(dheight1[last(index_height_range)+1] - max_height_offset) .. Inf])
                        coverage_above_height_range = area_above_height_range ./ geotile_total_area
                    else
                        coverage_above_height_range = 0
                    end

                    if (0 < coverage_below_height_range < maximum_extrap_fraction)
                        index_height_range = first(index_has_ice_range1):last(index_height_range)
                    end

                    if (0 < coverage_above_height_range < maximum_extrap_fraction)
                        index_height_range = first(index_height_range):last(index_has_ice_range1)
                    end

                    index_extrapolate = falses(length(dheight1))
                    index_extrapolate[index_height_range] .= true

                    index_interp_range, = validrange(index_height_valid)
                    index_extrapolate[index_interp_range] .= false

                    index_fill = falses(length(dheight1))
                    index_fill[index_height_range] .= true

                    # limit size of interpolation gap
                    index_height_range_interpolate, = validrange(index_height_valid)
                    index_height_interpolate = falses(length(dheight1))
                    index_height_interpolate[index_height_range_interpolate] .= true

                    if any(index_height_interpolate)
                        index_fill = index_fill .| index_height_interpolate
                    end

                    if sum(index_height_valid) > 2
                        for k in Symbol.(vars)
                            for date in ddate.val[index_date_range]
                                for pscale in dpscale
                                    itp = DataInterpolations.QuadraticInterpolation(gemb1[k][height=findall(index_height_valid), pscale=At(pscale), date=At(date)], dheight1.val[index_height_valid]; extrapolation=ExtrapolationType.Constant)
                                    gemb1[k][date=At(date), height=At(dheight1.val[index_fill]), pscale=At(pscale)] = itp.(dheight1.val[index_fill])
                                end
                            end
                        end
                    elseif sum(index_height_valid) > 1
                        for k in Symbol.(vars)
                            for date in ddate.val[index_date_range]
                                for pscale in dpscale
                                    itp = DataInterpolations.LinearInterpolation(gemb1[k][height=findall(index_height_valid), pscale=At(pscale), date=At(date)], dheight1.val[index_height_valid]; extrapolation=ExtrapolationType.Constant)
                                    gemb1[k][date=At(date), height=At(dheight1.val[index_fill]), pscale=At(pscale)] = itp.(dheight1.val[index_fill])
                                end
                            end
                        end
                    else
                        for k in Symbol.(vars)
                            fill_values = gemb1[k][height=findall(index_height_valid)]

                            for height0 in dheight1[index_fill]
                                gemb1[Symbol(k)][height=At(height0)] = fill_values
                            end
                        end
                    end

                    if show_interp_extrap_plots
                        bar_index = index_fill .& .!index_extrapolate .& .!index_height_valid
                        bar_labels = repeat(["interp-$(cnt)"], sum(index_has_ice1))
                        bar_labels[.!bar_index[index_has_ice1]] .= ""
                        barplot!(ax, dheight1[index_has_ice1].val, bar_index[index_has_ice1]; direction=:x, color=(cmap_bar[bar]), bar_labels=bar_labels, label_position=:center)
                        bar += 1

                        bar_index = index_extrapolate
                        bar_labels = repeat(["extrap-$(cnt)"], sum(index_has_ice1))
                        bar_labels[.!bar_index[index_has_ice1]] .= ""
                        barplot!(ax, dheight1[index_has_ice1].val, bar_index[index_has_ice1]; direction=:x, color=(cmap_bar[bar]), bar_labels=bar_labels, label_position=:center)
                        bar += 1
                    end

                    valid_ensemble = .!isnan.(gemb1[:fac][pscale=At(1)])
                    index_height_valid = any(valid_ensemble, dims=:date)[:]

                    coverage_height_fraction = 0
                    for ΔT in dΔT
                        Δheight = round((ΔT * (1000 / -6.5)) / dh1, digits=0) * dh1

                        height_effective = dheight1.val .- Δheight
                        index_height_effective = (height_effective .>= dheight.val[1]) .& (height_effective .<= dheight.val[end])

                        coverage_height_fraction += sum(geotile_hyps_area_km2[height=At(height_effective[index_height_effective])]) / geotile_total_area
                    end
                    coverage_height_fraction = coverage_height_fraction / length(dΔT)
                end

                if cnt > maximum_search_iterations
                    error("MAXIMUM SEARCH ITERATIONS REACHED: Buffering $(geotile_row.id) search extent by $(round(Int,search_buffer/1000)) km to populate elevation profile")
                end
            end

            # calculate volume change in unit of km3 [GEMB assumes an ice density of 910 kg/m3]
            index_height = (dheight1.val .>= dheight.val[1]) .& (dheight1.val .<= dheight.val[end])
            dΔT = dims(gemb_dv0, :ΔT)
          
            for k in Symbol.(vars)

                if (k in [:refreeze, :rain, :acc, :melt])
                    # these variables can not be negative change
                    v = diff(gemb1[k][date=index_date_range], dims=:date)
                    v[v.<0] .= 0

                    v1 = gemb1[k][date=index_date_range[1]:index_date_range[1]].data
                    v1[v1.<0] .= 0

                    v = cat(v1, v.data; dims=1)
                    v = cumsum(v; dims=1)
                     
                    gemb1[k][date=index_date_range] = v
                end

                for pscale in dpscale
                    v = gemb1[k][date=daterange, pscale=At(pscale)]
                    v0 = zeros(ddate, dheight)
                    v0 = v0[date=daterange]

                    for ΔT in dΔT
                        if elevation_classes_method == :Δelevation
                            Δheight = round((ΔT * (1000 / -6.5)) / dh1, digits=0) * dh1
                            height_effective = dheight1.val .- Δheight
                            index_height_effective = (height_effective .>= dheight.val[1]) .& (height_effective .<= dheight.val[end])
                            v0[:, :] = v[:, index_height_effective]

                        elseif elevation_classes_method == :mscale   
            
                            if k == :melt
                                v0[:, :] = v[:, index_height] .* ΔT
                            else
                                v0[:, :] = v[:, index_height]
                            end             
                        else
                            error("unknown elevation_classes_method: $(elevation_classes_method)")
                        end

                        dv = @d v0 .* geotile_hyps_area_km2 ./ 1000
                        dv[isnan.(dv)] .= 0
                        dv = sum(dv; dims=:height)
                        gemb_dv0[k][geotile=At(geotile_row.id), date=daterange, ΔT=At(ΔT), pscale=At(pscale)] = dropdims(dv, dims=:height)
                    end
                end
            end
        end

        show_interp_extrap_stats && printstyled("$(geotile_row.id): search buffer = $(round(Int,search_buffer/1000)) km, gemb nodes in search extent = $(round(Int,points_in_geotile_count/length(dΔheight)/length(dpscale)))\n"; color=:white)

        if show_interp_extrap_plots
            ax2 = CairoMakie.Axis(f[1, 1], xlabel="area [km²]", ylabel="height [m]")
            lines!(ax2, geotile_hyps_area_km2[index_has_ice].data, dheight[index_has_ice].val)
            display(f)
        end
    end

    # change from Δheight to ΔT
    if elevation_classes_method == :none
        # must be sorted in descending order for interpolation to work
        dΔT = Dim{:ΔT}(reverse(dΔheight.val .* (-6.5 / 1000)))
        gemb_dv0 = DimStack([DimArray(reverse(gemb_dv0[k].data, dims=4), (dgeotile, ddate, dpscale, dΔT); name=k) for k in Symbol.(vars)]...)
    end

    gemb_dv0 = gemb_add_derived_vars!(gemb_dv0)

    return gemb_dv0
end
