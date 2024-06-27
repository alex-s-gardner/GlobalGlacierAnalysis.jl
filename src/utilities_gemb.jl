function gemb_read(
    gemb_files;
    vars=[
        "latitude",
        "longitude",
        "datetime",
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
    any(vars .== "datetime") && (t = Vector{Vector{Float64}}())
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
            any(vars .== "datetime") && (t0 = read(file, "time"))
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
            any(vars .== "datetime") && (push!(t, vec(repeat(t0, m, 1)[:])))
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
    any(vars .== "datetime") && (gemb[!, :datetime] = Altim.decimalyear2datetime.(vcat(t...)))
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
       
    vars0=["latitude","longitude", "datetime", "smb", "fac", "ec", "acc", "runoff", "melt", "fac_to_depth", "height_ref", "refreeze", "t_air", "rain"]

    # variables that need to be converted from cumulitive outputs to rates

    var_cumulitive= ["smb", "ec", "acc", "runoff", "melt", "refreeze", "rain"]

    varmap= Dict(
        "latitude" => "lat", 
        "longitude" => "lon", 
        "datetime" => "time", 
        "smb" => "SMB", 
        "fac" => "FAC", 
        "ec" => "EC", 
        "acc" => "Accumulation", 
        "runoff" => "Runoff", 
        "melt" => "Melt", 
        "fac_to_depth" => "FACtoDepth", 
        "height_ref" => "H", 
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
        t = vec(gemb["datetime"])
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
            if k == "datetime"
                continue
            end
            
            if size(gemb[k],2) == size(gemb["datetime"],2)
                foo = fill(NaN, (size(gemb[k],1), length(ind)))
                for i in eachindex(ind)
                    if !isnothing(ind[i])
                        foo[:, i] = mean(gemb[k][:, ind[i]], dims=2)
                    end
                end
                gemb[k] = foo;
            end
        end
        gemb["datetime"] = reshape((datebin_edges[1:end-1] .+ datebin_edges[2:end]) ./ 2, (1, length(datebin_edges) - 1))
    end
    
    return gemb
end