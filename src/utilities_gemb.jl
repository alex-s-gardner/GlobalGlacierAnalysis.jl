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
            f[isnan.(f)] .= 0;
            gemb[v[1]] = f
        end
    else
        if length(intersect(vars0, vars)) != length(vars)
            error("one or more vars not recognized")
        end

        gemb = matopen(gemb_file) do file
            foo = read.(Ref(file), [varmap[v] for v in vars])
            for (i, v) in enumerate(vars)
                f = foo[i];
                f[isnan.(f)] .= 0
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
    gemb_bestfit_grouped(dv_altim, smb, fac, discharge, geotiles; examine_model_fits = nothing)

Calibrate the SMB model to grouped geotiles to handle glaciers that cross multiple geotiles.

# Arguments
- `dv_altim`: Volume change from altimetry data
- `smb`: Surface mass balance data
- `fac`: Firn air content data
- `discharge`: Discharge data
- `geotiles`: DataFrame containing geotile information
- `examine_model_fits`: Optional geotile ID to examine model fits for debugging 

# Returns
DataFrame with calibrated parameters for each geotile: id, extent, pscale, Δheight, and mad
"""
function gemb_bestfit_grouped(dv_altim, smb, fac, discharge, geotiles; examine_model_fits = nothing)

    # there are issues with calibrating the smb model to individual geotiles as glaciers
    # can cross multiple geotiles. To tackle this issue we to disconected geotile groups

    ddate = dims(dv_altim, :date)
    dΔheight = dims(smb, :Δheight)
    dpscale = dims(smb, :pscale)
    ddate_gemb = dims(smb, :date)
    dgeotile = dims(smb, :geotile)

    volume2mass = δice / 1000

    if .!isnothing(examine_model_fits)
        geotile_ref = examine_model_fits
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
    geotiles = copy(geotiles)
    geotiles[!, :pscale] .= 1.0;
    geotiles[!, :Δheight] .= 0.0;
    geotiles[!, :mad] .= 0.0;

    if .!isnothing(examine_model_fits)
        verbose = true
        ind = geotiles.id .== examine_model_fits;
        geotiles = geotiles[geotiles.group .== geotiles[ind, :group], :]
    else
        verbose = false
    end

    # loop for each group
    # NOTE: do not do a groupby on geotiles as you need to keep the original geotiles order

    # !!! POSSABLY NOT SAVE TO MULTI THREAD ON THIS LOOP !!!
    #Threads.@threads for grp in unique(geotiles.group)
    for grp in unique(geotiles.group)

        discharge_km3yr = 0.0
        gindex = grp .== geotiles.group

        for geotile = eachrow(geotiles[gindex, :])

            # total discharge D in Gt/yr converted to km3/yr
            index = within.(Ref(geotile.extent), discharge.longitude, discharge.latitude)

            if any(index)
                discharge_km3yr += sum(discharge[index, :discharge_gtyr]) ./ volume2mass
            end
        end

        pscale0 = 1
        Δheight0 = 0
        mad0 = Inf

        dv0 = dv_altim[At(geotiles[gindex, :id]), index_dv]
        dv0 = dropdims(sum(dv0, dims=:geotile), dims=:geotile)
        
        if verbose
            best_fit = zeros(size(dv0))
            best_smb = zeros(size(dv0))
            best_fac = zeros(size(dv0))
            best_discharge = zeros(size(dv0))
            best_offset = 0.

            f = Figure(size=(1000, 1000))
            f[1, 1] = CairoMakie.Axis(f)
        end

        for pscale in dpscale
            #pscale = first(dpscale)

            for Δheight in dΔheight
                #Δheight = first(dΔheight)
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
                #res0 = round(Int,mean(res))
                res0 = res[ceil(Int, length(res) / 2)]
                res0 = median(res)
                res = res .- res0;

                verbose && CairoMakie.lines!((smb0 .+ fac0 .- (discharge_km3yr .* Δdecyear)) .+ res0)

                if any(isnan.(res))
                    #display(dv0)
                    #display(discharge_km3yr)
                    #println(cname)

                    error("NANs in residual")
                end

                # try using the root mean square error instead of MAD...MAD did a poor job of caturing seasonality

                # ----------------- Objective function definition -[this should be passed as a kwarg]- ----------------
                # fit a model to the residuals... use the amplitude of the residuals as a weighting of the objective function
                fit1 = curve_fit(model3, Δdecyear, res, p3)
                mad1 = sqrt(mean(res .^ 2)) + (2*abs(fit1.param[4])) 
                # --------------------------------------------------------------

                if mad1 < mad0
                    mad0 = mad1
                    pscale0 = pscale
                    Δheight0 = Δheight

                    if verbose
                        best_fit = (smb0 .+ fac0 .- (discharge_km3yr .* Δdecyear)) .+ res0
                        best_smb = smb0
                        best_fac = fac0
                        # thi is +/- is done to get dates right for ploting
                        best_discharge = discharge_km3yr .* Δdecyear .+ best_fac .- best_fac 
                        best_offset = res0
                    end
                end
            end
        end


        if verbose
            CairoMakie.lines!(dv0; color = :black)
            CairoMakie.lines!(best_fit; color=:red)
            CairoMakie.lines!(best_fit; color=:red)
            display(f)

            f2 = Figure(size=(1000, 1000))
            p2 = f2[1, 1] = CairoMakie.Axis(f2, title="best fit for $examine_model_fits [pscale = $pscale0, Δheight = $Δheight0, mad = $(round(mad0; digits=2))]", xlabel="year", ylabel="km3")

            CairoMakie.lines!(best_smb; label="smb")
            CairoMakie.lines!(best_fac; label="fac")
            CairoMakie.lines!(best_discharge; label="discharge")
            CairoMakie.lines!(dv0 .- best_offset; label="dh - obs", color = :black, linewidth = 2)
            CairoMakie.lines!(best_fit .- best_offset; label="dh - model", color=:red, linewidth=2)

            f2[1, 2] = Legend(f2[1, 1], p2, framevisible=false)
            
            display(f2)
        end

        geotiles[gindex, :pscale] .= pscale0
        geotiles[gindex, :Δheight] .= Δheight0
        geotiles[gindex, :mad] .= mad0
    end
    return geotiles[:, [:id, :extent, :pscale, :Δheight, :mad]]
    
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