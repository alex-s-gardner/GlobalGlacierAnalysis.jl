
function gemb_read(gemb_files)

    lat = [];
    t = [];
    smb = [];
    fac = [];
    ec = [];
    acc = [];
    lon = [];
    runoff = [];
    melt = [];
    fac_to_depth = [];
    h = [];
    refreeze = [];
    ta = [];
    rain = [];


    for fn in gemb_files

        matopen(fn) do file
            lat0 = read(file, "lat");
            lon0 = read(file, "lon");
            t0 = read(file, "time");
            fac0 = read(file, "FAC");
            smb0 = read(file, "SMB")
            h0 = read(file, "H")
            ec0 = read(file, "EC")
            melt0 = read(file, "Melt")
            fac_to_depth0 = read(file, "FACtoDepth")
            refreeze0 = read(file, "Refreeze")
            ta0  = read(file, "Ta")
            rain0 = read(file, "Rain")

            runoff0 = read(file, "Runoff")
            acc0 = read(file, "Accumulation")
            m, n = size(fac0)

            push!(lat, vec(repeat(lat0,1, n)[:]))
            push!(lon, vec(repeat(lon0,1, n)[:]))
            push!(h, vec(repeat(h0, 1, n)[:]))
            push!(t, vec(repeat(t0,m, 1)[:]))
            push!(fac, fac0[:])
            push!(smb, smb0[:])
            push!(ec, ec0[:])
            push!(melt, melt0[:])
            push!(fac_to_depth, fac_to_depth0[:])
            push!(refreeze, refreeze0[:])
            push!(ta, ta0[:])
            push!(rain, rain0[:])
            push!(runoff, runoff0[:])
            push!(acc, acc0[:])
        end
    end
        gemb = DataFrame(
            datetime=Altim.decimalyear2datetime.(vcat(t...)), 
            latitude=vcat(lat...), 
            longitude=vcat(lon...), 
            height_ref=vcat(h...), 
            fac=vcat(fac...), 
            smb=vcat(smb...),
            ec=vcat(ec...),
            acc=vcat(acc...),
            runoff=vcat(runoff...),
            melt=vcat(melt...),
            fac_to_depth=vcat(fac_to_depth...),
            refreeze=vcat(refreeze...),
            t_air=vcat(ta...),
            rain=vcat(rain...),
        )
    return gemb
end