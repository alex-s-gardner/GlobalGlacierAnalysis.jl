using LsqFit
using Statistics
using CSV
using DataFrames
using NCDatasets
using Dates
using GeoStats
using Arrow

# density of ice
const di = 910; #kg m^3
# Define model that will be fit to all data binned by hypsometry
#model::Function = model(t, h; t_intercept=2010) = hcat(ones(size(t)), (t .- t_intercept), h, h.^2, cos.(2 * pi * t), sin.(2 * pi * t))
model1::Function = model1(x, p) =
    p[1] .+
    p[2] .* x[:, 1] .+
    p[3] .* x[:, 1] .^ 2 .+
    p[4] .* x[:, 2] .+
    p[5] .* x[:, 2] .^ 2 .+
    sin.(2 .* pi .* (x[:, 1] .+ p[6])) .* (p[7] .+ p[8] .* x[:, 2] .* p[9] .* x[:, 2] .^ 2)

const p1 = zeros(9);
const lb1 = [-10.0, -3.0, -2.0, -0.05, -0.0001, -1.0, -7.0, -0.05, -0.001];
const ub1 = [+10.0, +3.0, +2.0, +0.05, +0.0001, +1.0, +7.0, +0.05, +0.001];

# seasonal only 
model1_seasonal::Function = model1_seasonal(x, p) = sin.(2 .* pi .* (x[:, 1] .+ p[6])) .* (p[7] .+ p[8] .* x[:, 2] .* p[9] .* x[:, 2] .^ 2)

# including quadratic for seasonal does not improve std(anom)
#(p[6] .* cos.(2 .* pi .* x[:, 1]) .+  p[7].* sin.(2 .* pi .* x[:, 1])) .* (1 .+ p[8] .* x[:, 2] .+ p[9] .* x[:, 2] .^ 2)
#p1 = zeros(9);
#lb1 = [-10., -3., -2., -.05, -0.0001, -10., -10., -0.01, -0.0001];
#ub1 = [+10., +3., +2., +.05, +0.0001, +10., +10., +0.01, +0.0001];

# model fit across all geotiles for a region for a given year
model2::Function = model2(h, p) = p[1] .+ p[2] .* h .+ p[3] .* h .^ 2;
const p2 = zeros(3);
const lb2 = [-30.0, -0.1, -0.01];
const ub2 = [+30.0, +0.1, 0.01];

model3::Function = model3(t, p) = p[1] .+ p[2] .* t .+ p[3] .* t .^ 2 .+ p[4] .* sin.(2 .* pi .* (t .+ p[5]))
const p3 = zeros(5)

function geotiles_mutually_exclusive_rgi!(geotiles) 
    reg = names(geotiles);
    reg = reg[contains.(reg, "rgi")];
    rgi_id = replace.(reg, Ref("rgi" => ""));
    rgi_id = parse.(Int64, rgi_id);

    # make regions mutually exclusive by assigning geotiles to region of largest overlap
    for geotile in eachrow(geotiles[!, reg])
        maxind = findfirst(maximum(geotile) .== collect(geotile))
        for r in eachindex(reg)
            if r == maxind
                continue
            else
                geotile[r] = 0
            end
        end
    end
    return geotiles, reg
end


function hyps_model_fill!(dh1, nobs1, params; bincount_min=5, model1_madnorm_max=5, smooth_n=9, smooth_h2t_length_scale=800, variogram_range_ratio = false)

    t = Altim.decimalyear.(dims(dh1[first(keys(dh1))], :date))
    t = repeat(t, 1, length(dims(dh1[first(keys(dh1))], :height)))

    h = val(dims(dh1[first(keys(dh1))], :height))'
    h = repeat(h, length(dims(dh1[first(keys(dh1))], :date)), 1)
    
    dgeotiels = dims(dh1[first(keys(dh1))], :geotile)
    foo = DimArray(fill(NaN, length(dgeotiels)), dgeotiels)

    if variogram_range_ratio
        range_ratio = Dict();
        for mission in keys(dh1)
            push!(range_ratio, String(mission) => copy(foo))
        end
    end

    for mission in keys(dh1)

        # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><>
        #mission = "hugonnet"
        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

        # find valid range for entire mission so that all filled geotiles cover the same date range
        rrange, = Altim.validrange(vec(any(.!isnan.(dh1[mission]), dims=(1, 3))))

        Threads.@threads for geotile in dims(dh1[mission], :geotile)
            #for geotile in dims(dh1[mission], :geotile)
            # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><>
            #geotile =       "lat[+80+82]lon[+058+060]"
            # geotile = first(dims(dh1[mission], :geotile))

            # geotile = geotiles[findfirst((geotiles.rgi1 .> 0.) .& (geotiles.glacier_frac .> 0.3)),:]

            # geotile = "lat[-72-70]lon[-014-012]"
            # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

            # println("$mission: $geotile")

            k = findfirst(params[mission].geotile .== geotile)
            df = @view params[mission][k, :]

            dh0 = dh1[mission][At(geotile), :, :]
            nobs0 = nobs1[mission][At(geotile), :, :]
            df.nobs_raw = sum(nobs0)
            df.nbins_raw = sum(nobs0 .> 0)

            ###################################### FILTER 1 ################################
            valid1 = .!isnan.(dh0) .& (nobs0 .> bincount_min[mission]) .& (abs.(dh0) .< 200)
            ################################################################################
            dh0[.!valid1] .= NaN
            dh1[mission][At(geotile), :, :] = dh0
            nobs0[.!valid1] .= 0
            nobs1[mission][At(geotile), :, :] = nobs0

            # if there are not enough points to fit a model the set all to NaNs
            va = sum(valid1)
            if va <= (length(p1) + 2)
                dh1[mission][At(geotile), :, :] .= NaN
                nobs1[mission][At(geotile), :, :] .= 0
                continue
            end

            # determine valid range of data
            (_, crange) = Altim.validrange(valid1)
            valid0 = valid1[rrange, crange]

            # if there are not enough points to fit a model the set all to NaNs
            va = sum(valid0)
            if va <= (length(p1) + 2)
                dh1[mission][At(geotile), :, :] .= NaN
                nobs1[mission][At(geotile), :, :] .= 0
                continue
            end

            dh0 = dh0[rrange, crange];
            nobs0 = nobs0[rrange, crange];
            t0 = t[rrange, crange]
            h0 = h[rrange, crange]

            # center predictors and observations
            t0_mean = round(mean(t0)) # remove an exact integer to keep phase from shifting
            df.t0 = t0_mean
            t0 = t0 .- t0_mean

            h0_mean = mean(h0[:])
            df.h0 = h0_mean
            h0 = h0 .- h0_mean

            dh0_median = median(dh0[valid0]) + 0.0000001 # add a small offset to prevent numerical instability
            df.dh0 = dh0_median
            dh0 = dh0 .- dh0_median

            df.bin_std = std(dh0[valid0])

            # fit global model 
            fit1 = []
            try
                fit1 = curve_fit(model1, hcat(t0[valid0], h0[valid0]), dh0[valid0], nobs0[valid0], p1; lower=lb1, upper=ub1)
            catch
                fit1 = curve_fit(model1, hcat(t0[valid0], h0[valid0]), dh0[valid0], nobs0[valid0], p1)
            end

            dh0_mdl = model1(hcat(t0[valid0], h0[valid0]), fit1.param)
            dh0_anom = dh0[valid0] .- dh0_mdl

            ###################################### FILTER 2 ####################################
            # filter model1_madnorm_max sigma outliers
            valid0[valid0] = Altim.madnorm(dh0_anom) .<= model1_madnorm_max
            vb = sum(valid0)
            df.nbins_filt1 = vb

            # if there are not enough points to fit a model then set all to NaNs
            if vb <= (length(p1) + 2)
                dh1[mission][At(geotile), :, :] .= NaN
                nobs1[mission][At(geotile), :, :] .= 0
                continue
                continue
            end

            if vb < va
                dh0[.!valid0] .= NaN
                nobs0[.!valid0] .= 0
                fit1 = curve_fit(model1, hcat(t0[valid0], h0[valid0]), dh0[valid0], nobs0[valid0], p1; lower=lb1, upper=ub1)

                dh0_mdl = model1(hcat(t0[valid0], h0[valid0]), fit1.param)
                dh0_anom = dh0[valid0] .- dh0_mdl
            end

            ####################################################################################
            if variogram_range_ratio
                coord = [(a,b) for (a,b) in zip(t0[valid0], h0[valid0])]
                table = (; Z=vec(dh0_anom))
                geotable = georef(table, coord)

                try
                    vg1 = DirectionalVariogram((1.0, 0.0), geotable, :Z)
                    f1 = fit(SphericalVariogram, vg1)
                    vg2 = DirectionalVariogram((0.0, 1.0), geotable, :Z)
                    f2 = fit(SphericalVariogram, vg2)
                    range_ratio[mission][At(geotile)] = round(Int, (f2.ball.radii ./ f1.ball.radii)[1])
                    println("elevation / height range: $(range_ratio[mission][At(geotile)])")

                catch
                end


                continue
            end

            df.bin_anom_std = std(dh0_anom)
            nobs1[mission][At(geotile), rrange, crange] = nobs0

            # add final parameters to DataFrame
            df.nobs_final = sum(nobs0)
            df.param_m1 = fit1.param

            # take the median of the x closest neighbors 
            if sum(valid0) < smooth_n[mission]
                anom_smooth = zeros(length(dh0))
            else
                # scale height distance relative to time (i.e. length-scale)
                pts = hcat(t0[valid0], h0[valid0] / smooth_h2t_length_scale)'
                kdtree = KDTree(pts)
                (idxs, _) = knn(kdtree, pts, smooth_n[mission])
                anom0 = map(ind -> median(dh0_anom[ind]), idxs)

                # extrema(anom0)
                # interpolate anomalies using weighted distance (Shepard(2))
                itp = ScatteredInterpolation.interpolate(Shepard(2), pts, anom0)
                pts = hcat(t0[:], h0[:] / smooth_h2t_length_scale)'
                anom_smooth = vec(ScatteredInterpolation.evaluate(itp, pts))
            end

            # fill out valid range (no extraploation beyond (rrange,crange) of geotile with the model only
            dh1[mission][At(geotile), rrange, crange] = model1(hcat(t0[:], h0[:]), fit1.param) .+ dh0_median .+ anom_smooth
            # println("granule interp: $(mission) - $(geotile)")
        end
    end

    if variogram_range_ratio
        return range_ratio
    end
end


function hyps_amplitude_normalize!(dh1, params; mission_reference = "icesat2")

    t = Altim.decimalyear.(dims(dh1[first(keys(dh1))], :date))
    t = repeat(t, 1, length(dims(dh1[first(keys(dh1))], :height)))

    h = val(dims(dh1[first(keys(dh1))], :height))'
    h = repeat(h, length(dims(dh1[first(keys(dh1))], :date)), 1)

    for mission in keys(dh1)
        #mission = "hugonnet"

        if mission == mission_reference
            continue
        end

        for geotile in dims(dh1[mission], :geotile)
            #geotile = "lat[-72-70]lon[-014-012]"
            #geotile = first(dims(dh1[mission], :geotile))
            k = findfirst(params[mission].geotile .== geotile)
            df0 = params[mission][k, :]
            dfr = params[mission_reference][k, :]

            if any(isnan.(df0.param_m1)) || any(isnan.(dfr.param_m1))
                continue
            end

            dh0 = dh1[mission][At(geotile), :, :]
            valid = .!isnan.(dh0)
        
            (rrange, crange) = Altim.validrange(valid)

            dh0 = dh0[rrange, crange]
            valid = valid[rrange, crange]
            t0 = t[rrange, crange] .- df0.t0
            h0x = h[rrange, crange] .- df0.h0

            p0 = df0.param_m1
            pr = dfr.param_m1
            Δp = pr .- p0
            Δp[1:6] = p0[1:6]

            delta = model1_seasonal(hcat(t0[:], h0x[:]), Δp)
            if !any(isnan.(delta)) # there seem to be rare cases where model1_seasonal returns nans.
                dh1[mission][At(geotile), rrange, crange] = dh0[:] .+ delta
            end
        end
    end
end

function hyps_fill_empty!(dh1, params, geotiles; mask = :glacier)

    lon = mean.(getindex.(geotiles.extent, :X))
    lat = mean.(getindex.(geotiles.extent, :Y))
    ll = [(a, b) for (a, b) in zip(lon, lat)]

    for mission in keys(dh1)
        # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><
        #mission = "hugonnet"
        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

         # find valid range for entire mission so that all filled geotiles cover the same date range
        rrange, = Altim.validrange(vec(any(.!isnan.(dh1[mission]), dims=(1, 3))))
        dh0_median = getindex.(params[mission][:, :param_m1], 1) .+ params[mission][:, :dh0]
        
        # copy dh as it gets morphed inside of parallel loop
        dh2 = copy(dh1[mission])

        Threads.@threads for geotile in (dims(dh1[mission], :geotile))
            #for geotile in dims(dh1[mission], :geotile)
            # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><>
            #geotile = first(dims(dh1[mission], :geotile)[empty_geotiles])
            #geotile = first(dims(dh1[mission], :geotile))

            # geotile = geotiles[findfirst((geotiles.rgi1 .> 0.) .& (geotiles.glacier_frac .> 0.3)),:]
            #geotile = "lat[-80-78]lon[+166+168]"
            # <><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

            # fill with median of nearest neighbors if no data
            k = findfirst(geotiles.id .== geotile)
            has_ice = geotiles[k, "$(mask)_area_km2"] .> 0

            if !any(has_ice)
                continue
            else
                crange, = Altim.validrange(has_ice)
            end

            dh0 = @view dh1[mission][At(geotile), rrange, crange]

            if all(isnan.(dh0))

                # find distance between goetiles
                dist2geotile = haversine.(Ref(ll[k]), ll, 6371000)

                # find closest X valid geotiles
                nnearest = 5;
                
                #sorted index
                idx = sortperm(dist2geotile)

                # has data in desired elevation range
                has_data = vec(any(.!isnan.(dh2[:, rrange, crange]), dims=(2, 3)));
                has_data = has_data[idx]

                # check that there are enough geotiles with overlapping elevation ranges
                nnearest0 = sum(has_data)
                if nnearest0 < nnearest
                    if nnearest0 == 0
                        printstyled("$mission $geotile: no geotiles with overlapping elevation ranges to fill geotile, setting dh .= 0 \n"; color=:red, bold=true)
                        dh0 .= 0;
                        continue
                    else
                        printstyled("$mission $geotile: less than $nnearest geotiles with overlapping elevation ranges to fill geotile, using closest $nnearest0 geotiles instead\n"; color=:light_red, bold=false)
                        nnearest = nnearest0
                    end
                end

                idx = idx[has_data][1:nnearest]
                f1 = dh2[idx, rrange, crange]

                # remove mean offset to normalize between regions
                dh0_median0 = dh0_median[idx]
                for g in eachindex(dh0_median0)
                    f1[g, :, :] .-= dh0_median0[g]
                end

                for i in 1:length(dims(dh0, :date))
                    for j in 1:length(dims(dh0, :height))
                        f2 = vec(f1[:,i, j])
                        valid00 = .!isnan.(f2)
                        if any(valid00)
                            dh0[i, j] = median(f2[valid00])
                        end
                    end
                end

                # if there are any gaps along an elevation profile then they need to be filled
                (rrange0, crange0) = Altim.validrange(.!isnan.(dh0))
                if any(isnan.(dh0[rrange0, crange0]))
                    x = 1:size(dh0,2)
                    for i in rrange0
                        y = vec(dh0[i,:]);
                        valid = .!isnan.(y)
                        if any(.!valid)
                            itp = DataInterpolations.LinearInterpolation(y[valid], x[valid])
                            y[.!valid] = itp(x[.!valid])
                            dh0[i,:] = y;
                        end
                    end
                end

                # add back median model offset 
                dh0 .+= median(dh0_median0)
            end
        end
    end
end

function hyps_fill_updown!(dh1, geotiles; mask = :glacier)

    for mission in keys(dh1)
        # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><
        # mission = "hugonnet"
        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

        rrange, = Altim.validrange(vec(any(.!isnan.(dh1[mission]), dims=(1, 3))))

        Threads.@threads for geotile in (dims(dh1[mission], :geotile))

            # fill with median of nearest neighbors if not data
            k = findfirst(geotiles.id .== geotile)
            valid = geotiles[k, "$(mask)_area_km2"] .> 0

            if !any(valid)
                continue
            else
                crange, = Altim.validrange(valid)
            end
            dh0 = @view dh1[mission][At(geotile), rrange, crange]

            # fill lowest hand highest elevations with first and last values
            valid = .!isnan.(dh0)
            if .!any(valid) || all(valid)
                continue
            end

            for i in 1:length(dims(dh0, 1))
                if any(valid[i, :])
                    f = findfirst(vec(valid[i, :]))
                    l = findlast(vec(valid[i, :]))
                    dh0[i, 1:f] .= dh0[i, f]
                    dh0[i, l:end] .= dh0[i, l]
                end
            end
        end
    end
end

function plot_dvdm(
    dv, 
    dm, 
    nobs;
    title=nothing,
    fontsize=18,
    colors = palette(:Set1_4, length(keys(dh1))),
    date_intercept = 2012,
    area = NaN,
    dmdm_flag=false,
    δ_effective=""
    )
    
    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98), size=(1000, 700), fontsize=fontsize)

    ga = f[1:5, 1] = GridLayout()
    gb = f[6, 1] = GridLayout()

    if !isnothing(title)
        Label(ga[1, 1, Top()],
            title, valign=:bottom,
            font=:bold,
            padding=(0, 0, 5, 0)
        )
    end

    decdate = Altim.decimalyear.(dims(dv, :date).val)
    kwargs = (; xminorgridvisible=true, xminorticks=IntervalsBetween(5))
    if dmdm_flag
        ylabel="mass/volume anomaly [Gt]"
    else
        ylabel="mass/volume anomaly [Gt/km³]"
    end

    axmain = Axis(ga[1, 1]; ylabel, kwargs...)

    hidexdecorations!(axmain; grid=false, minorgrid=false)

    ylabel = dims(nobs, :mission).val;
    ylabel = replace.(ylabel, "hugonnet" => "ASTER")
    ylabel = replace.(ylabel, "icesat" => "ICESat")
    ylabel = replace.(ylabel, "icesat2" => "ICESat-2")
    ylabel = replace.(ylabel, "gedi" => "GEDI")

    axbottom = Axis(gb[1, 1]; yticks=(1:4, lowercase.(ylabel)), yticklabelcolor = colors, 
        yminorgridvisible=true, yminorticks=IntervalsBetween(2), 
        yminorgridwidth=3, kwargs...)

    hideydecorations!(axbottom; grid=true, minorgrid=false, ticklabels=false)

    linkxaxes!(axmain, axbottom)

    # fit model to average across all missions
    CairoMakie.xlims!(axbottom, 2000, 2024)

    dm0, dm_fit = region_fit(dm, date_intercept)
    dv0, dv_fit = region_fit(dv, date_intercept)

    # random error
    rand_err = std(dm_fit.resid) *2

    # correlated error [take as a fraction of the fac correction]
    correlated_err = 0.3 .* vec(abs.(dv0 .- dm0))

    # combined error 
    err = correlated_err.+ rand_err

    valid = .!isnan.(dm0.data)
    if !dmdm_flag
        CairoMakie.band!(axmain, decdate[valid], dm0.data[valid] .- err[valid] , dm0.data[valid] .+ err[valid], color=(:black, 0.05))
    end

    if hasdim(dm, :mission)
        
        for (i, mission) in enumerate(dims(dm,:mission))
            dv0 = dv[At(mission),:]
            dm0 = dm[At(mission),:]
            if dmdm_flag
                CairoMakie.lines!(axmain, decdate, vec(dv0), color=(colors[i], 0.5), linewidth=2, label="850 kg m³")
            else
                CairoMakie.lines!(axmain, decdate, vec(dv0), color=(colors[i], 0.2), linewidth=2)
            end
            CairoMakie.lines!(axmain, decdate, vec(dm0), label="$(ylabel[i])", color=(colors[i], 1), linewidth=2)
        end
    else
        if dmdm_flag
            CairoMakie.lines!(axmain, decdate, vec(dv), color=(:black, 0.5), linewidth=2)
        else
            CairoMakie.lines!(axmain, decdate, vec(dv), color=(:black, 0.2), linewidth=2)
        end
        CairoMakie.lines!(axmain, decdate, vec(dm), color=(:black, 1), linewidth=2)
    end

    if dmdm_flag
        text = "δ:                   $(δ_effective) kg m³
trend:            $(round(dm_fit.param[2], digits = 1)) Gt yr⁻¹
acceleration: $(round(dm_fit.param[3], digits = 1)) Gt yr⁻²
amplitude:     $(abs(round(dm_fit.param[4], digits = 1))) Gt"
    else
        text = "trend:            $(round(dm_fit.param[2], digits = 1)) Gt yr⁻¹
acceleration: $(round(dm_fit.param[3], digits = 1)) Gt yr⁻²
amplitude:     $(abs(round(dm_fit.param[4], digits = 1))) Gt
area:             $(round(Int64,area)) km²"
    end

    text!(
        axmain, 0.70, 0.95,
        text = text,
        #font=:bold,
        color=(:black, 1),
        align=(:left, :top),
        space=:relative,
        fontsize=fontsize
    )

    if dmdm_flag
     text!(
        axmain, 0.70, 0.75,
            text="δ:                   850 kg m³
trend:            $(round(dv_fit.param[2], digits = 1)) Gt yr⁻¹
acceleration: $(round(dv_fit.param[3], digits = 1)) Gt yr⁻²
amplitude:     $(abs(round(dv_fit.param[4], digits = 1))) Gt",
        #font=:bold,
        color=(:black, 0.5),
        align=(:left, :top),
        space=:relative,
        fontsize=fontsize
    )
    end

    var0 = log.(nobs.data');
    var0[isinf.(var0)] .= NaN

    hm = CairoMakie.heatmap!(axbottom, decdate, 1:4, var0, colormap=Reverse(:magma),)

    rowgap!(ga, 1)
    rowgap!(gb, 1)

    CairoMakie.Colorbar(f[6, 2], hm; label="log(count)"); #, vertical=false, bbox=axbottom.scene.px_area)

    return f
end

function plot_dm_permutations(t, dm, title; colors=colorschemes[:tab20c], date_intercept=2012)
    f = Figure(size=(1500, 500); fontsize);
    #title = "Randolph Glacier Inventory: Region $(rgi[4:end])"
    Label(f[0, :], text=title)
    ax = Axis(f[1, 1:3], ylabel="mass anomaly Gt")
    l1 = []
    dm_fit = []
    
    for (i, r) in enumerate(eachrow(dm))
        dm0 = DimArray(r.dm_gt, Dim{:date}(Altim.decimalyear2datetime.(t)))
        _, dm_fit0 = region_fit(dm0, date_intercept)
        push!(dm_fit,dm_fit0)
        label = "bin = $(r.binning_method), dem = $(r.dem_id), curve = $(r.curvature_correct), amp = $(r.amplitude_correct)"
        l1 = lines!(ax, t, r.dm_gt; label, color = colors[i], linewidth=2)
    end

    tr = round.(extrema([f.param[2] for f in dm_fit]), digits = 1)
    acc = round.(extrema([f.param[3] for f in dm_fit]), digits = 2)
    amp = round.(extrema([abs(f.param[4]) for f in dm_fit]), digits = 1)

     text!(
        0.60, 0.95,
        text=
"trend:            $(tr[1]) to $(tr[2]), Gt yr⁻¹
acceleration: $(acc[1]) to $(acc[2]) Gt yr⁻²
amplitude:     $(amp[1]) to $(amp[2]) Gt",
        #font=:bold,
        color=(:black, 0.5),
        align=(:left, :top),
        space=:relative,
        fontsize=fontsize
    )

    Legend(f[1, 4], ax, framevisible=false)
    return f
end


function region_fit(dm, date_intercept)
    decdate = Altim.decimalyear.(dims(dm, :date).val)

    if hasdim(dm, :mission)
        valid = .!isnan.(dm)
        dmm =  dropdims(sum(dm .* valid, dims=1) ./ sum(valid, dims=1), dims = :mission)
        valid = vec(.!isnan.(dmm))
    else
        valid = .!isnan.(dm)
        dmm = dm
    end
    dm_fit = curve_fit(model3, decdate[valid] .- date_intercept, dmm[valid], p3)
    return dmm, dm_fit
end


function hyps_volume_change(dh1, nobs1, geotiles; mask = :glacier)
    (ngeotile, ndate, nheight) = size(dh1[first(keys(dh1))])
    missions = collect(keys(dh1))
    nmissions = length(missions)
    dmission = Dim{:mission}(missions)
    
    dv = DimArray(fill(NaN, nmissions, ngeotile, ndate), (dmission, dims(dh1[first(missions)], :geotile), dims(dh1[first(missions)], :date)))
    nobs = DimArray(fill(NaN, nmissions, ngeotile, ndate), (dmission, dims(dh1[first(missions)], :geotile), dims(dh1[first(missions)], :date)))
                
    var_area = reduce(hcat, geotiles[:, "$(mask)_area_km2"])'
    var_area = permutedims(repeat(var_area, 1, 1, ndate), (1, 3, 2))

    # if there is no glacier ice then set dhdt to zero
    for mission in keys(dh1)
        #mission = "icesat2"
        rrange, = Altim.validrange(vec(any(.!isnan.(dh1[mission]), dims=(1, 3))))
        area_mask = var_area .== 0

        if rrange.start > 1
            area_mask[:, 1:(rrange.start-1), :] .= false
        end

        if rrange.stop < ndate
            area_mask[:, (rrange.stop+1):end, :] .= false
        end

        dh1[mission][area_mask] .= 0
    end

    for mission in keys(dh1)
        dv[At(mission), :, :] = dropdims(sum(dh1[mission] .* var_area, dims=3), dims=3) / 1000
    end

    for mission in keys(nobs1)
        nobs[At(mission), :, :] = dropdims(sum(nobs1[mission], dims=3), dims=3)
    end

    return dv, nobs
end

function hyps_geotile_aggrigate(dv, geotiles, reg; fun=sum::Function)
    dmission, dgeotile, ddate = dims(dv)
    drgi = Dim{:rgi}(reg)

    dv_reg = DimArray(fill(NaN, length(dmission), length(reg), length(ddate)), (dmission, drgi, ddate))
    
    for mission in dmission
        for rgi in drgi
            rgi_ind = geotiles[:, rgi] .> 0
            geotile_ids = geotiles[rgi_ind, :].id
            dv_reg[At(mission), At(rgi), :] = vec(fun(dv[At(mission), At(geotile_ids), :], dims=1))
        end
    end

    return dv_reg
end

function hyps_offset2reference!(dv_reg, mission_reference)
    offset = copy(dv_reg[:, :, 1])
    offset[:] .= 0;

    for mission in dims(offset, :mission)
        if mission ==  mission_reference 
            continue
        end
        
        for rgi in dims(offset, :rgi)

            dv1 = dv_reg[At(mission_reference),At(rgi),:]
            valid1 = .!isnan.(dv1)

            dv2 = dv_reg[At(mission),At(rgi),:]
            valid2 = .!isnan.(dv2)

            overlap = valid1 .& valid2

            if any(overlap)
                offset[At(mission),At(rgi)] = median(dv1[overlap]) - median(dv2[overlap])
            end
        end
    end

    for mission in dims(offset, :mission)
        for rgi in dims(offset, :rgi)
            dv_reg[At(mission), At(rgi), :] .+= offset[At(mission), At(rgi)]
        end
    end

    return offset, dv_reg
end






function hyps_fac_correction_old(fac, smb, nobs_gemb, dh1, geotiles, reg; mask = :glacier)
    
    _, _, dheight = dims(smb)
    ddate = dims(dh1[first(keys(dh1))], :date)
    drgi = Dim{:rgi}(reg)
    h0 = collect(dheight)
    
    facv_reg = DimArray(fill(0.0, length(drgi), length(ddate), length(dheight)), (drgi, ddate, dheight))
    smbv_reg = copy(facv_reg)

    # map gemb date into dh date
    d1 = Altim.decimalyear.(collect(dims(smb,:date)))
    kdtree = KDTree(d1')
    d2 = Altim.decimalyear.(collect(ddate))
    idxs, dists = nn(kdtree, d2')

    # find first day of valid data 
    valid_date = falses(length(ddate))
    for mission in keys(dh1)
        valid_date .|= vec(any(.!isnan.(dh1[mission]), dims=(1, 3)))
    end
    first_valid = findfirst(valid_date)
    

    for rgi in reg
        #  <><><><><><><><><><><><><><><><><><> FOR TESTING  <><><><><><><><><><><><><><><><><><>
        # rgi = "rgi10"
        #  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

        rgi_ind = geotiles[:, rgi] .> 0
        geotile = geotiles[rgi_ind, :]

        # remove climatolotgy (assume some rate of vertical flow (divergence))

        fac1 = fac[At(geotile.id),idxs,:]
        nobs_gemb1 = nobs_gemb[At(geotile.id),idxs,:]
        smb1 = smb[At(geotile.id), idxs, :]
        valid1 = .!isnan.(fac1)

        if !any(valid1)
            continue
        end
        area_rgi = vec(sum(reduce(hcat, geotile[:, "$(mask)_area_km2"]), dims=2))

        # average accross all geotiles
        fac_foo = fac1 .* nobs_gemb1;
        smb_foo = smb1 .* nobs_gemb1
        fac_foo[.!valid1] .= 0;
        smb_foo[.!valid1] .= 0;

        nobs_gemb_foo = valid1 .* nobs_gemb1
        fac_rgi = dropdims(sum(fac_foo, dims=1) ./ sum(nobs_gemb_foo, dims=1), dims=1)
        smb_rgi = dropdims(sum(smb_foo, dims=1) ./ sum(nobs_gemb_foo, dims=1), dims=1)

        valid = .!isnan.(fac_rgi)

        # remove first values for each elevation band to get fac anomaly
        for i = 1:size(fac_rgi,2)
            if any(valid[:, i])
                f = findfirst(vec(valid[:, i]))
                fac_rgi[:, i] .-= fac_rgi[f, i]
            end
        end

        # take cumulitive sum of smb
        for i = 1:size(smb_rgi, 2)
            if any(valid[:, i])
                smb_rgi[:, i] = cumsum(smb_rgi[:, i])
            end
        end

        # set unmodeled elevations == to first and last modeled levelations
        valid = .!isnan.(fac_rgi);
        (crange,_) = Altim.validrange(valid)

        for i in crange
            if any(valid[i,:])
                f = findfirst(vec(valid[i,:]))
                l = findlast(vec(valid[i, :]))

                fac_rgi[i, 1:f] .= fac_rgi[i,f]
                fac_rgi[i, l:end] .= fac_rgi[i,f]

                smb_rgi[i, 1:f] .= smb_rgi[i,f]
                smb_rgi[i, l:end] .= smb_rgi[i,f]
            end
        end

        # convert from m to km.^3
        for r in eachrow(fac_rgi)
            r[:] = r .* area_rgi / 1000
        end

        # convert from mm to km.^3
        for r in eachrow(smb_rgi)
            r[:] = r .* area_rgi ./ (1000^2)
        end

        # set unmodeled times == to first and last modeled elevations
        valid = .!isnan.(fac_rgi)
        (_, rrange) = Altim.validrange(valid)
        for i in rrange
            if any(valid[:,i])
                f = findfirst(vec(valid[:, i]))
                l = findlast(vec(valid[:, i]))

                fac_rgi[1:f,i] .= fac_rgi[f,i]
                fac_rgi[l:end,i] .= fac_rgi[l,i]

                smb_rgi[1:f,i] .= smb_rgi[f,i]
                smb_rgi[l:end,i] .= smb_rgi[l,i]
            end
        end

        # fill gap if it exists
        for r in eachrow(fac_rgi)
            if any(isnan.(r))
                valid = .!isnan.(r)
                itp = DataInterpolations.LinearInterpolation(r[valid].data, h0[valid])
                r[.!valid] = itp(h0[.!valid])
            end
        end

        for r in eachrow(smb_rgi)
            if any(isnan.(r))
                valid = .!isnan.(r)
                itp = DataInterpolations.LinearInterpolation(r[valid].data, h0[valid])
                r[.!valid] = itp(h0[.!valid])
            end
        end

        facv_reg[At(rgi), :, :] = fac_rgi 
        smbv_reg[At(rgi), :, :] = smb_rgi 

    end

    return facv_reg, smbv_reg
end

function hyps_scale_fac!(facv_reg, smbv_reg, dv_reg)
   
    fac_scale = fill(0., dims(facv_reg, :rgi))

    for rgi in dims(facv_reg, :rgi)
        fac1 = facv_reg[At(rgi), :]
        smb1 = smbv_reg[At(rgi), :]
        geb_dv1 = (smb1 ./ (di/1000)) .+ fac1
        dv1 = dv_reg[At(mission_reference), At(rgi), :]

        # last valid date of gemb that was not interpolated
        if any(isnan.(facv_reg))
            lastvalid = dims(facv_reg, :date)[findlast(any(.!isnan.(facv_reg), dims=(1, 2)))[2]]
        else
            lastvalid = last(dims(facv_reg, :date))
        end


        # findoverlap
        valid1 = .!isnan.(fac1)
        valid2 = .!isnan.(dv1)
        overlap = vec((valid1 .& valid2) .& (dims(dv1, :date) .<= lastvalid))
        t = decimalyear.(dims(dv1, :date))
        
        # find amplitude
        dv_fit = curve_fit(model3, t[overlap], dv1[overlap], p3)

        dgemb_fit = curve_fit(model3, t[overlap], geb_dv1[overlap], p3)

        sf = dv_fit.param[4] / dgemb_fit.param[4]
        if isinf(sf)
            sf = 0
        end
        fac_scale[At(rgi)] = sf
    end

    # apply scale factor
    for rgi in dims(facv_reg, :rgi)
        facv_reg[At(rgi), :, :] .*= fac_scale[At(rgi)]
        smbv_reg[At(rgi), :, :] .*= fac_scale[At(rgi)] ./ (di / 1000) # convert to volume
    end

    return fac_scale, facv_reg, smbv_reg
end


function hyps_mass_change(dv_reg, facv_reg)

    dm_reg = copy(dv_reg)
    for rgi in dims(dm_reg, :rgi)

        dfac = facv_reg[At(rgi), :]
            
        for mission in dims(dm_reg, :mission)

            dv0 = dv_reg[At(mission), At(rgi), :]
            dm_reg[At(mission), At(rgi), :] = (vec(dv0) .- vec(dfac)) * (di/1000)
        end
    end
    return dm_reg
end

function plot_altim_grace(
    altim,
    grace;
    zemp = nothing,
    title=nothing,
    fontsize=18,
    date_intercept=2012,
    area=NaN
)

    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98), size=(1000, 700), fontsize=fontsize)

    ga = f[1, 1] = GridLayout()

    if !isnothing(title)
        Label(ga[1, 1, Top()],
            title, valign=:bottom,
            font=:bold,
            padding=(0, 0, 5, 0)
        )
    end

    kwargs = (; xminorgridvisible=true, xminorticks=IntervalsBetween(5))
    axmain = Axis(ga[1, 1]; ylabel="mass anomaly [Gt]", kwargs...)

    #hidexdecorations!(axmain; grid=false, minorgrid=false)

    # fit model to average across all missions
    xlim = (2000, 2024)
    CairoMakie.xlims!(axmain, xlim)

    # make sure overlaping dates are used for model fit
    valid = .!isnan.(altim["dm_gt"])
    aex = extrema(altim["date"][valid])

    valid = .!isnan.(grace["dm_gt"])
    gex = extrema(grace["date"][valid])

    if isnothing(zemp)
        ex = (max(aex[1], gex[1]), min(aex[2], gex[2]))
    else
        valid = .!isnan.(zemp["dm_gt"])
        zex = extrema(zemp["date"][valid])
        ex = (max(aex[1], gex[1], zex[1]), min(aex[2], gex[2], zex[2]))
    end

    valid = .!isnan.(altim["dm_gt"]) .& (altim["date"] .>= ex[1]) .& (altim["date"] .<= ex[2])
    altim_fit = curve_fit(model3, altim["date"][valid] .- date_intercept, altim["dm_gt"][valid], p3)
    valid = .!isnan.(altim["dm_gt"])
    CairoMakie.band!(axmain, altim["date"][valid], altim["dm_gt"][valid] .- altim["2sigma_gt"][valid], altim["dm_gt"][valid] .+ altim["2sigma_gt"][valid], color=(:black, 0.05))
    CairoMakie.lines!(altim["date"], altim["dm_gt"]; label="altimetry")

    valid = .!isnan.(grace["dm_gt"]) .& (grace["date"] .>= ex[1]) .& (grace["date"] .<= ex[2])
    grace_fit = curve_fit(model3, grace["date"][valid] .- date_intercept, grace["dm_gt"][valid], p3)

    # align grace with altim
    Δdm = altim_fit.param[1] - grace_fit.param[1]
    valid = .!isnan.(grace["dm_gt"])
    gdm = grace["dm_gt"] .+ Δdm
    CairoMakie.band!(axmain, grace["date"][valid], gdm[valid] .- grace["2sigma_gt"][valid], gdm[valid] .+ grace["2sigma_gt"][valid], color=(:black, 0.05))
    CairoMakie.lines!(grace["date"], gdm; label="gravimetry")

    zdm = [];
    if !isnothing(zemp)
        valid = .!isnan.(zemp["dm_gt"]) .& (zemp["date"] .>= ex[1]) .& (zemp["date"] .<= ex[2])
        zemp_fit = curve_fit(model2, zemp["date"][valid] .- date_intercept, zemp["dm_gt"][valid], p2)

        display(zemp_fit.param)
        # align zemp with altim
        Δdm = altim_fit.param[1] - zemp_fit.param[1]
        valid = .!isnan.(zemp["dm_gt"])
        zdm = zemp["dm_gt"] .+ Δdm
        CairoMakie.band!(axmain, zemp["date"][valid], zdm[valid] .- zemp["2sigma_gt"][valid], zdm[valid] .+ zemp["2sigma_gt"][valid], color=(:black, 0.05))
        CairoMakie.lines!(zemp["date"], zdm; label="in situ")

        textxoffset = 0.60;
        txt = "trend:            $(round(Int, altim_fit.param[2])) $(round(Int, grace_fit.param[2]))  $(round(Int, zemp_fit.param[2])) Gt yr⁻¹
acceleration: $(round(altim_fit.param[3], digits = 1)) $(round(grace_fit.param[3], digits = 1)) $(round(Int, zemp_fit.param[3])) Gt yr⁻²
amplitude:     $(abs(round(Int, altim_fit.param[4]))) $(abs(round(Int, grace_fit.param[4]))) Gt
area:             $(round(Int64,area)) km²"
    else
        textxoffset = 0.70
        txt = "trend:            $(round(Int, altim_fit.param[2])) $(round(Int, grace_fit.param[2])) Gt yr⁻¹
acceleration: $(round(altim_fit.param[3], digits = 1)) $(round(grace_fit.param[3], digits = 1)) Gt yr⁻²
amplitude:     $(abs(round(Int, altim_fit.param[4]))) $(abs(round(Int, grace_fit.param[4]))) Gt
area:             $(round(Int64,area)) km²"
    end

    text!(
        axmain, textxoffset, 0.95,
        text = txt,
        #font=:bold,
        color=(:black, 0.5),
        align=(:left, :top),
        space=:relative,
        fontsize=fontsize
    )

    axislegend(axmain; position=:lb, framevisible=false)

    # find range of visible data 
    valid = (altim["date"] .>= xlim[1]) .& (altim["date"] .<= xlim[2]) .& .!isnan.(altim["dm_gt"])
    aex = extrema(altim["dm_gt"][valid])
    
    valid = (grace["date"] .>= xlim[1]) .& (grace["date"] .<= xlim[2]) .& .!isnan.(grace["dm_gt"])
    gex = extrema(gdm[valid])

    if isnothing(zemp)
        ylim = (min(aex[1], gex[1]), max(aex[2], gex[2]))
    else
        valid = (zemp["date"] .>= xlim[1]) .& (zemp["date"] .<= xlim[2]) .& .!isnan.(zemp["dm_gt"])
        zex = extrema(zdm[valid])
        ylim = (min(aex[1], gex[1], zex[1]), max(aex[2], gex[2], zex[2]))
    end

    display(ylim)
    CairoMakie.ylims!(axmain, ylim)

    return f
end

function read_zemp2019(;datadir= setpaths().zemp_2019)
    # Read Zemp 2019 data
    fn_startswith = "Zemp_etal_results_region_"
    fn_endswith = ".csv"
    files = allfiles(datadir; fn_startswith, fn_endswith)

    all = Dict()
    for file in files
        # pull RGI number form file name
        foo = file[end-10:(end-7)]
        s = findfirst('_', foo)
        e = findlast('_', foo)
        rgi = "rgi$(foo[(s+1):(e-1)])"
        data = DataFrame((CSV.File(file; header=28, skipto=29, stripwhitespace=true)))
        push!(all, rgi => data)
    end

    date = DateTime.(1950:2016)
    rgi = "rgi" .* string.(1:19)

    drgi = Dim{:rgi}(rgi)
    ddate = Dim{:date}(date)

    dm_gt = DimArray(fill(NaN, (length(drgi), length(ddate))), (drgi, ddate))
    err_gt = DimArray(fill(NaN, (length(drgi), length(ddate))), (drgi, ddate))

    for r in rgi
        d = DateTime.(all[r].Year)
        dm_gt[At(r),At(d)] = cumsum(all[r].INT_Gt)
        err_gt[At(r),At(d)] = all[r].sig_Total_Gt
    end

    return (; dm_gt, err_gt, all)
end


function read_marzeion2020(;datadir= setpaths().marzeion_2020)

    foo = Dataset(datadir)
    rgi = "rgi" .* string.(Int.(foo["Region"]))
    date = DateTime.(foo["Time"])
    climate_model= getindex.(collect(foo["Climate_Model"].attrib), 2)
    glacier_model= getindex.(collect(foo["Glacier_Model"].attrib), 2)
    scenario = getindex.(collect(foo["Scenario"].attrib), 2)

    drgi = Dim{:rgi}(rgi)
    ddate = Dim{:date}(date)
    dclimate_model = Dim{:climate_model}(climate_model)
    dglacier_model = Dim{:glacier_model}(glacier_model)
    dscenario = Dim{:scenario}(scenario)

    foo = collect(foo["Mass"].var);
    dm_gt = DimArray(foo, (drgi, ddate, dclimate_model, dglacier_model, dscenario))

    return (;dm_gt)
end

function read_marzeion2012(; datadir=setpaths().marzeion_2012)
    foo = matread(datadir)

    scenario = collect(keys(foo))

    # for s in scenario
    s = first(scenario)
    date = vec(DateTime.(foo[s]["yr"]))
    climate_model = vec(foo[s]["model"])
    rgi = vec("rgi" .* string.(Int.(foo[s]["rgiid"])))
    gt2sle = foo[s]["gt2sle"] .* -1

    drgi = Dim{:rgi}(rgi)
    ddate = Dim{:date}(date)
    dclimate_model = Dim{:climate_model}(climate_model)
    dscenario = Dim{:scenario}(scenario)


    dm_gt = fill(NaN, (length(rgi), length(date), length(climate_model), length(scenario)))
    err_gt = fill(NaN, (length(rgi), length(date), length(climate_model), length(scenario)))

    dm_gt = DimArray(dm_gt, (drgi, ddate, dclimate_model, dscenario))
    err_gt = DimArray(err_gt, (drgi, ddate, dclimate_model, dscenario))

    for (i, s) in enumerate(scenario)
        climate_model0 = vec(foo[s]["model"])
        foo1 = foo[s]["mb"] ./ gt2sle
        foo1 = permutedims(foo1, (2, 1, 3))
        dm_gt[:, :, At(climate_model0), i] = foo1

        foo1 = foo[s]["mb_err"] ./ gt2sle
        foo1 = permutedims(foo1, (2, 1, 3))
        err_gt[:, :, At(climate_model0), i] = foo1
    end

    foo = (;dm_gt, err_gt)

    return foo
end


function plot_dm(
    dm;
    title=nothing,
    fontsize=18,
    colors=colorschemes[:hawaii100],
    scale_color = true,
    date_intercept=2012,
    date_alignall=(2000, 2005),
    area=NaN,
    stats_in_label = false, 
    xlim=(2000, 2024)
)

    f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98), size=(1000, 700), fontsize=fontsize)
    ga = f[1, 1] = GridLayout()

    # sample by number of 
    if scale_color
        colors = colors[1:floor(Int, length(colors) ./ length(dm)):end]
    end

    if !isnothing(title)
        Label(ga[1, 1, Top()],
            title, valign=:bottom,
            font=:bold,
            padding=(0, 0, 5, 0)
        )
    end

    kwargs = (; xminorgridvisible=true, xminorticks=IntervalsBetween(5))
    axmain = Axis(ga[1, 1]; ylabel="mass anomaly [Gt]", kwargs...)

    # fit model to average across all missions
    CairoMakie.xlims!(axmain, xlim)

    # make sure overlaping dates are used for model fit
    ex = Tuple{Float64,Float64}
    for (c, d) in enumerate(dm)
        valid = .!isnan.(d["dm_gt"])
        aex = extrema(d["date"][valid])
        if c == 1
            ex = aex
        else
            ex = (max(aex[1], ex[1]), min(aex[2], ex[2]))
        end
    end

    # allign all data to first mass change timeseries
    dm_fit = []
    dm_ref =
        for (c, d) in enumerate(dm)
            valid = .!isnan.(d["dm_gt"]) .& (d["date"] .>= ex[1]) .& (d["date"] .<= ex[2])
            push!(dm_fit, curve_fit(model2, d["date"][valid] .- date_intercept, d["dm_gt"][valid], p3))

            valid = .!isnan.(d["dm_gt"]) .& (d["date"] .>= date_alignall[1]) .& (d["date"] .<= date_alignall[2])
            Δdm = mean(d["dm_gt"][valid])

            valid = .!isnan.(d["dm_gt"])
            d["dm_gt"] .-= Δdm
            CairoMakie.band!(axmain, d["date"][valid], d["dm_gt"][valid] .- d["2sigma_gt"][valid], d["dm_gt"][valid] .+ d["2sigma_gt"][valid], color=(:black, 0.05))

            if stats_in_label
                label = "$(d["label"]): $(round(Int, dm_fit[c].param[2])) Gt yr⁻¹ $(round(dm_fit[c].param[3], digits=1)) Gt yr⁻²"
            else
                label = d["label"]
            end
            CairoMakie.lines!(d["date"], d["dm_gt"]; label, color=colors[c])

        end

    trend = []
    acc = []
    for f in dm_fit
        trend = push!(trend, round(Int, f.param[2]))
        acc = push!(acc, round(f.param[3], digits=1))
    end

    axislegend(axmain; position=:lb, framevisible=false)

    # find range of visible data 
    ylim = Tuple{Float64,Float64}
    for (c, d) in enumerate(dm)
        valid = (d["date"] .>= xlim[1]) .& (d["date"] .<= xlim[2]) .& .!isnan.(d["dm_gt"])
        aex = extrema(d["dm_gt"][valid])
        if c == 1
            ylim = aex
        else
            ylim = (min(aex[1], ylim[1]), max(aex[2], ylim[2]))
        end
    end

    CairoMakie.ylims!(axmain, ylim)
    
    d1 = Date(floor(Altim.decimalyear2datetime(ex[1]), Day))
    d2 = Date(floor(Altim.decimalyear2datetime(ex[2]), Day))
    txt = "trend fit to overlaping period: $(d1) to $(d2)"
    text!(
        axmain, 0.95, 0.00,
        text=txt,
        #font=:bold,
        color=(:black, 0.5),
        align=(:right, :bottom),
        space=:relative,
        fontsize=fontsize/1.5
    )
    return f
end

# read in hock_2019
function read_hock2019(; datadir=setpaths().hock_2019)
    #datadir= setpaths().hock_2019
    foo = Dataset(datadir)
    rgi = "rgi" .* string.(Int.(foo["region"]))
    date = DateTime.(foo["time"])

    # from metadata
    scenario = ["A1B", "RCP26", "RCP45", "RCP60", "RCP85"]
    glacier_model = ["BM", "GieOer", "GloGEM", "GloGEMebal", "HYOGA2", "RHB", "VASlangen"]
    climate_model = ["ACCESS1-0", "BCC-CSM1-1", "BNU-ESM", "CCSM4", "CGCM3-1(T63)", "CNRM-CM3", "CNRM-CM5", "CSIRO-Mk3-0", "CSIRO-Mk3-6-0", "CanESM2", "ECHAM5-MPI-OM", "GFDL-CM2-0", "GFDL-CM2-1", "GFDL-CM3", "GFDL-ESM2G", "GFDL-ESM2M", "GISS-E2-R", "HadGEM2-AO", "HadGEM2-CC", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-LR", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC-ESM-CHEM", "MIROC5", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-CGCM3", "NorESM1-M", "NorESM1-ME", "PCM", "UKMO-HadCM3"]

    drgi = Dim{:rgi}(rgi)
    ddate = Dim{:date}(date)
    dclimate_model = Dim{:climate_model}(climate_model)
    dglacier_model = Dim{:glacier_model}(glacier_model)
    dscenario = Dim{:scenario}(scenario)

    coord = (drgi, ddate, dclimate_model, dglacier_model, dscenario)
    dm_gt = DimArray(fill(NaN, length.(coord)), coord)

    for i in 1:length(dclimate_model)
        for j in 1:length(dglacier_model)
            for k in 1:length(dscenario)
                idx = (foo["scenario"] .== k) .& (foo["glaciermodel"] .== j) .& (foo["forcingmodel"] .== i)
                f = foo["volume"].var[:, :, findall(idx)]
                if isempty(f)
                    continue
                else
                    dm_gt[:, :, i, j, k][:] = f
                end
            end
        end
    end

    return (; dm_gt)
end

function read_grace_rgi(;datadir=setpaths()[:grace_rgi])

    grace0 = matread(datadir)

    # rgi leter to digit mapping
    old2new = Dict(
        "PAT" => "rgi18", 
        "GRE" => "rgi5", 
        "NAS" => "rgi10", 
        "ICE" => "rgi6", 
        "TRP" => "rgi16",
        "SAW" => "rgi14",
        "SAE" => "rgi15",
        "CDN" => "rgi3",
        "CAS" => "rgi13",
        "AND" => "rgi17",
        "CDS" => "rgi4",
        "ANT" => "rgi19",
        "CEU" => "rgi11",
        "ALA" => "rgi1",
        "SVB" => "rgi7",
        "WNA" => "rgi2", 
        "NEZ" => "rgi18",
        "RAI" => "rgi9",
        "SCA" => "rgi8"
    )


    grace = Dict(any(keys(old2new).== key) ? (old2new[key]) => val : (key) => val for (key, val) in grace0)
    return grace
end


function plot_height_time(dh1; geotile, fig_suffix, fig_folder, figure_suffix, mask = :glacier)
    for mission in keys(dh1)
        height = dims(dh1[mission], :height)
        valid = geotile["$(mask)_area_km2"] .> 0
        crange, = Altim.validrange(geotile["$(mask)_area_km2"] .> 0)

        # find valid range for entire mission so that all filled geotiles cover the same date range
        rrange, = Altim.validrange(vec(any(.!isnan.(dh1[mission]), dims=(1, 3))))

        # i need to plot twice, once to get colorbar and once to align x-axis.. makie does not yet support catigorical axis.. but coming soon
        color = cgrad(:balance, rev=true);
        clim = (-10, 10)

        fname = joinpath(fig_folder, "height_time_colormap.png")
        f = heatmap(rand(10, 10); clim, color, rightmargin=5Plots.mm)

        savefig(f, fname)
        display(f)

        for legend = [false]
            f = plot(layout=grid(2, 1, heights=[0.8, 0.2]), link=:x, framestyle=:box, legend=legend)

            heatmap!(dh1[mission][At(geotile.id), rrange, crange]; subplot=1, legend=false, framestyle=:box, clim, color)
            xlabel!("")
            ylabel!("")

            bar!(val(height[crange]), geotile["$(mask)_area_km2"][crange], subplot=2, legend=false, framestyle=:box, link=:x)
            ylabel!("area [km²]", subplot=2)
            xlabel!("height [m]")

            if legend
                fname = joinpath(fig_folder, "$(figure_suffix)_$(mission)_$(geotile.id)_binned_$(fig_suffix)_colorbar.png")
            else
                fname = joinpath(fig_folder, "$(figure_suffix)_$(mission)_$(geotile.id)_binned_$(fig_suffix).png")
            end
            savefig(f, fname)
            display(f)
        end
    end
end

# rsync -r devon:/mnt/bylot-r3/data/binned/2deg/figures /Users/gardnera/Research/20_01_GlobalGlacierChange/version\ 2/
function read_ipccar6(;datadir=setpaths()[:ipcc_ar6], start2007=false)
    # Data from IPCC AR6 Figure 9.21
    fn_endswith = ".csv"
    files = allfiles(datadir; fn_endswith)

    f2007 = contains.(files, Ref("2007"))
    if start2007
        files = files[f2007]
    else
        files = files[.!f2007]
    end

    err_files = contains.(files, Ref("error"));
    files_err = files[err_files]
    files = files[.!err_files]
    all = Dict()

    if start2007
        scenarios = [f[end-13:end-9] for f in files]
    else
        scenarios = [f[end-8:end-4] for f in files]
    end

    dscenario = Dim{:scenario}(scenarios)

    data = DataFrame((CSV.File(files[1]; header=1, skipto=2, stripwhitespace=true)))
    ridx = contains.(names(data), Ref("RGI"))

    rgi = "rgi" .* [n[5:end] for n in names(data)[ridx]]
    drgi = Dim{:rgi}(rgi)

    date = DateTime.(data[:, 1])
    ddate = Dim{:date}(date)

    dm_gt = DimArray(fill(NaN, length(drgi), length(ddate), length(dscenario)), (drgi, ddate, dscenario))
    err_gt = copy(dm_gt)

    for k in eachindex(files)
        # pull RGI number form file name
        data = DataFrame((CSV.File(files[k]; header=1, skipto=2, stripwhitespace=true)))
        for j = 1:19
            dm_gt[j,:,k] = data[:,j+1]
        end

        data = DataFrame((CSV.File(files_err[k]; header=1, skipto=2, stripwhitespace=true)))
        for j = 1:19
            err_gt[j, :, k] = data[:, j+1]
        end
    end

    foo = (; dm_gt, err_gt)

    return foo
end

function geotile_binarea!(geotile, ras, feature, bin_edges; invert=false, excludefeature=nothing)
        
    t1 = time();
    ras0 = read(Rasters.crop(ras, to=geotile.extent));

    # using rasterize is 15 times faster than reading its_live masks !!!!
    mask0 = Rasters.rasterize!(count, zeros(ras0.dims), feature; threaded=false, shape=:polygon, progress=false, verbose=false) .> 0;

    if invert
        mask0 = .!mask0
    end

    if .!isnothing(excludefeature)
        excludemask = Rasters.rasterize!(count, zeros(ras0.dims), excludefeature; threaded=false, shape=:polygon, progress=false, verbose=false) .> 0
        mask0 = mask0 .& .!excludemask
    end

    # calculate area per cell
    lon = lookup(ras0, X)
    lat = lookup(ras0, Y)
    d = Altim.meters2lonlat_distance.(Ref(1), lat)
    a = abs.((1 ./ getindex.(d, 2) * (lat[2] .- lat[1])) .* (1 / d[1][1] * (lon[2] - lon[1])))
    area_m2 = repeat(a', outer = [length(lon), 1])

    feature_area_km2 = mask0 .* area_m2 / (1000 * 1000)

    foo = geotile[out_var_name]

    df = DataFrame(ras = ras0[:], feature_area_km2 = feature_area_km2[:])
    dfbs = binstats(df, :ras, bin_edges, :feature_area_km2, col_function = sum, missing_bins = true)
    valid = .!(ismissing.(dfbs.feature_area_km2_sum))
    foo[valid] .= dfbs.feature_area_km2_sum[valid]

    t = round(time() - t1)
    printstyled("    -> $(geotile.id) hypsometry calculated: $(t)s\n"; color=:light_black)
    
    return geotile
end

function geotiles_mask_hyps(surface_mask, geotile_width)

    binned_folder = analysis_paths(; geotile_width).binned
    gt_file = joinpath(binned_folder, "geotile_$(surface_mask)_hyps.arrow");
    geotiles = DataFrame(Arrow.Table(gt_file));
    geotiles.extent = Extent.(getindex.(geotiles.extent, 1));
    geotiles = copy(geotiles);
    return geotiles
end




function plot_dh(dh_reg, w; title = "", xlims=(DateTime(2000), DateTime(2024)))
    dmission = dims(dh_reg, :mission)
    dmetric = Dim{:metric}(["mean", "trend", "acceleration", "amplitude", "date_intercept"])
    fit_param = fill(NaN, (dmission, dmetric))
    ddate = dims(dh_reg, :date)
    decdate = Altim.decimalyear.(ddate)
    for mission in dmission
        dh0 = dh_reg[At(mission), :]
        valid = .!isnan.(dh0)
        if any(valid)
            fit_param[At(mission), At("date_intercept")] = round(Int, mean(decdate[valid]))
            dh_fit = curve_fit(model3, decdate[valid] .- fit_param[At(mission), At("date_intercept")], dh0[valid], collect(w[At(mission), :][valid]), p3;)
            fit_param[At(mission), At("mean")] = dh_fit.param[1]
            fit_param[At(mission), At("trend")] = dh_fit.param[2]
            fit_param[At(mission), At("acceleration")] = dh_fit.param[3]
            fit_param[At(mission), At("amplitude")] = dh_fit.param[4]
        end
    end

    label = collect(dmission)
    prec = 2;
    label = replace.(label, "hugonnet" => "aster")
    
    dh_mean = round.(fit_param[:, At("mean")], digits=prec)
    trend = round.(fit_param[:, At("trend")], digits=prec)
    acc = round.(fit_param[:, At("acceleration")], digits=prec)
    amp = abs.(round.(fit_param[:, At("amplitude")], digits=prec))

    for i in eachindex(label)
        label[i] = "$(label[i]): $(dh_mean[i]) m $(trend[i]) m yr⁻¹ $(acc[i]) m yr⁻² $(amp[i]) m"
    end

    label = reverse(permutedims(label))
    p = plot(
        dh_reg';
        linewidth=2,
        ylabel="height anomaly [m]",
        label,
        legend=:outerright,
        margin=5Plots.mm,
        title,
        legendtitle="mission: mean trend acc amp",
        xlims
    )

    plot!(size=(1200, 400))
    return p, fit_param
end

function landoffset!(dh1, geotiles, landfit)

    t = Altim.decimalyear.(lookup(dh1[first(keys(dh1))], :date))
    height = lookup(dh1[first(keys(dh1))], :height)
    for rgi in lookup(landfit, :rgi)
        rgi_ind = findall(geotiles[:, rgi] .> 0)

        for mission in lookup(landfit, :mission)
            t0 = t .- landfit[At(mission), At(rgi), At("date_intercept")]
            land_dh = landfit[At(mission), At(rgi), At("mean")] ; #.+ landfit[At(mission), At(rgi), At("trend")] .* t0 .+ landfit[At(mission), At(rgi), At("acceleration")] .* t0.^2
            for i in rgi_ind
                for k in eachindex(height)
                    dh1[mission][i, :, k] .= dh1[mission][i, :, k] .- land_dh
                end
            end
        end
    end
    return dh1
end


function binned_filepath(binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct)
    if curvature_correct
        runid = "$(surface_mask)_dh_$(dem_id)_cc_$(binning_method)_$(project_id)"
    else
        runid = "$(surface_mask)_dh_$(dem_id)_$(binning_method)_$(project_id)"
    end

    binned_file = joinpath(binned_folder, "$(runid).jld2");

    return binned_file
end

function  binned_filled_filepath(binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)

    if curvature_correct
        runid = "$(surface_mask)_dh_$(dem_id)_cc_$(binning_method)_$(project_id)"
    else
        runid = "$(surface_mask)_dh_$(dem_id)_$(binning_method)_$(project_id)"
    end
    
    if amplitude_correct
        binned_filled_file = joinpath(binned_folder, "$(runid)_filled_ac_p$(paramater_set).jld2")
        out_id = replace(runid, "dh" => "dv")
        figure_suffix = "$(out_id)_ac"
    else
        binned_filled_file = joinpath(binned_folder, "$(runid)_filled_p$(paramater_set).jld2")
        figure_suffix = replace(runid, "dh" => "dv")
    end

    return binned_filled_file, figure_suffix
    
end
