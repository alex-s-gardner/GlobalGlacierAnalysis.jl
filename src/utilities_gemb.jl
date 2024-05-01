
function gemb_fac(fac_files)
    lat = [];
    lon = [];
    fac1 = [];
    t = [];
    h = []

    for fn in fac_files
        file = matopen(fn)
        lat0 = read(file, "lat");
        lon0 = read(file, "lon");
        t0 = read(file, "time");
        fac0 = read(file, "FAC");
        h0 = read(file, "H")
        m, n = size(fac0)

        push!(lat, vec(repeat(lat0,1, n)[:]))
        push!(lon, vec(repeat(lon0,1, n)[:]))
        push!(h, vec(repeat(h0, 1, n)[:]))
        push!(t, vec(repeat(t0,m, 1)[:]))
        push!(fac1, fac0[:])
    end
    fac = DataFrame(datetime=Altim.decimalyear2datetime.(vcat(t...)), latitude=vcat(lat...), longitude=vcat(lon...), height_ref=vcat(h...), fac=vcat(fac1...))
    return fac
end
