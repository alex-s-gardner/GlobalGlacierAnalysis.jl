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
    any(vars .== "date") && (gemb[!, :date] = Altim.decimalyear2datetime.(vcat(t...)))
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

function gem_Δvolume!_old(df)

    volume2mass = Altim.δice / 1000
    dates = DataFrames.metadata(df, "date")
    decyear = Altim.decimalyear.(dates)
    Δdecyear = decyear .- decyear[1]

    index_smb = df.var .== "smb_km3"
    index_fac = df.var .== "fac_km3"

    df_smb = df[index_smb, :]
    df_fac = df[index_fac, :]

    df_dv = copy(df_smb)
    df_dv[!, :var] .= "dv_km3"
    df_dv[!, :val] .= [fill(NaN, length(decyear))]
    df_dv[!, :nobs] .= [fill(0, length(decyear))]

    df_dm = copy(df_dv)
    df_dm[!, :var] .= "dm_gt"

    Threads.@threads for i in 1:nrow(df_dv)

        smb = df_smb[i, :val]
        fac = df_fac[i, :val]
        rgi = df_smb[i, :rgi]
        nobs0 = df_smb[i, :nobs]

        discharge = Altim.discharge_gtyr[rgi]

        if isnan(discharge) && (rgi == "rgi19")
            # discharge set equal to mass balance from icesat period
            index = (decyear .> 2003) .& (decyear .<= 2009)
            x = Δdecyear[index] .- mean(Δdecyear[index])
            y = smb[index] .- mean(smb[index])
            discharge = x \ y * volume2mass

            println(discharge)

            # println("Δheight = $(r.Δheight), pscale = $(r.pscale): discharge = $discharge Gt")
        elseif isnan(discharge) && (rgi == "rgi12")
            # discharge set equal to mass balance from icesat period
            index = (decyear .> 2000) .& (decyear .<= 2005)
            x = Δdecyear[index] .- mean(Δdecyear[index])
            y = smb[index] .- mean(smb[index])
            discharge = x \ y * volume2mass

            #println("Δheight = $(r.Δheight), pscale = $(r.pscale): discharge = $discharge Gt")
        end

        discharge = discharge * Δdecyear

        dv_km3 = (smb / volume2mass) .+ fac .- (discharge / volume2mass)
        df_dv[i, :val] = dv_km3
        df_dv[i, :nobs] = nobs0

        dm_gt = smb .- discharge
        df_dm[i, :val] = dm_gt
        df_dm[i, :nobs] = nobs0
    end

    df = append!(df, df_dv)
    df = append!(df, df_dm)

    return df
end


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


function gembscale(var0, gemb_fit)
    # initialize fields and group
    height_step = step(dims(var0, :height))
    var2 = similar(var0[:, :, :, 1])

    # doing as DimensionalData adds NO increase in processing time
    for geotile in dims(var0, :geotile)
        ind = findfirst(gemb_fit.id .== geotile)
        Δheight = gemb_fit[ind, :Δheight]
        pscale = gemb_fit[ind, :pscale]

        # select appropriate precipitation scaling 
        var2[At(geotile), :, :] = var0[At(geotile), :, :, At(pscale)]

        # shift elevation to simulate lowering and raising the elevation of the model
        height_shift = round(Int16(Δheight / height_step))

        if height_shift > 0
            var2[At(geotile), :, 1:(end-height_shift)] = var2[At(geotile), :, (height_shift+1):end]
        elseif height_shift < 0
            var2[At(geotile), :, (-height_shift+1):end] = var2[At(geotile), :, 1:(end+height_shift)]
        end
    end
    return var2
end