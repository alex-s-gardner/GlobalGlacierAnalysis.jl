using Altim, DataFrames, CairoMakie, Arrow, Statistics, Dates, BinStatistics, Printf


outFigFolder = "/home/gardnera/Documents/Figures/Altim/offsets"
project_id = :v01
geotile_width = 2
domain = :all

paths = project_paths(project_id = project_id)
geotiles = project_geotiles(; geotile_width = geotile_width, domain = domain)

# geotiles = geotiles[ind,:]
dems = [:cop30_v2, :nasadem_v1, :arcticdem_v3_10m] #, :rema_v2_10m, :arcticdem_v3_10m];
sensors = [:icesat, :icesat2, :gedi]

sensor = sensors[3]
for sensor in sensors

    # plot offsets
    offsets = [DataFrame() for dem in dems]

    @time for i in eachindex(dems)
        suffix = ".$(dems[i])_offset"

        for geotile in eachrow(geotiles)
            fn = joinpath(paths[sensor].geotile, geotile.id * suffix)

            if isfile(fn)
                append!(offsets[i], DataFrame(Arrow.Table(fn)))
            end
        end
    end

    datetimes = [DateTime.(zeros(size(offset.forward))) for offset in offsets]
    detector_ids = [zeros(UInt8,size(offset.forward)) for  offset in offsets]

    function gedibeam2num(track)
        num = parse(UInt8, replace(track,"BEAM" => "" ); base=2)
        return num
    end

    @time for i in eachindex(dems)
        suffix = ".$(dems[i])_offset"
        start = 1;
        if sensor == :icesat2
            beamvar = :beam_id;
        elseif sensor == :gedi
            beamvar = :track;
        else
            beamvar = nothing;
        end

        b = copy(datetimes[i]);
        c = copy(detector_ids[i])

        for (i, geotile) in enumerate(eachrow(geotiles))
            fn = joinpath(paths[sensor].geotile, geotile.id * suffix)

            if isfile(fn)
                # extract date time 
                fn = joinpath(paths[sensor].geotile, geotile.id * ".arrow")
                foo = DataFrame(Arrow.Table(fn))
                valid = .!isempty.(foo.datetime)
                stop = start + length(valid) -1

                a = DateTime.(zeros(size(valid)))
                d = zeros(UInt8,size(valid))
                if any(valid)

                    if sensor == :icesat2
                        d[valid] = foo.detector_id[valid];
                    elseif sensor == :gedi
                        d[valid] = gedibeam2num.(foo.track[valid]);
                    else
                        d[valid] .= UInt8(1)
                    end

                    a[valid] = first.(foo.datetime[valid])
                end

                b[start:stop] = a;
                c[start:stop] = d;
                start = stop+1;
            end
        end
        datetimes[i] = b
        detector_ids[i] = c
    end

    mad(x) = median(abs.(x .- median(x)))

    for fun in [median, mad]
        for i in 1:length(dems)
            vnames = [:right, :forward, :above]
            for vname in vnames

                #vname = :right
                delta_bin = :daily
                
                if delta_bin == :weekly
                    db = 1/52;
                elseif delta_bin == :daily
                    db = 1/365;
                elseif delta_bin == :hourly12
                    db = 0.5/365;
                elseif delta_bin == :monthly
                    db = 1/12;
                end

                valid = .!isnan.(offsets[i][:, vname])
                y = offsets[i][!, vname][valid];
                weight = offsets[i].weight[valid];
                x = decimalyear.(datetimes[i][valid]);
                d = detector_ids[i][valid];
                
                deltax = 0.5
                minx = floor(minimum(x) ./ deltax) .* deltax
                maxx = ceil(maximum(x) ./ deltax) .* deltax

                p = sortperm(x)
                x = x[p];
                y = y[p]
                d = d[p]

                du = sort(unique(d))

                f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
                    resolution = (1000, 700), title = "test")

                ga = f[1, 1]
                ax = Axis(ga[1,1])

                bins = minimum(x):db:maximum(x)

                oset = 0;
                df0 = DataFrame(x = x, y = y)
                df3 = binstats(df0, :x, bins, :y; col_function = fun)
                delta = ceil(std(df3[!,Symbol("y_$(string(fun))")])*5 / .25) * .25

                for u in du
                    v = d .== u
                    w = weight[v]

                    #ry = rolling(sum, y[v].*w, ws) ./ rolling(sum, w, ws)
                        #y0 = vcat(zeros(eltype(y), padd), ry, zeros(eltype(y), padd)); 

                    df0 = DataFrame(x = x[v], y = y[v])
                    df3 = binstats(df0, :x, bins, :y; col_function = fun)
                    
                    x0 = bincenter.(df3[!,1])
                    y0 = df3[!,Symbol("y_$(string(fun))")]

                    med_offset = round(median(y0), digits = 2)
                    # inset nan's for large temporal gaps
                    dx0 = push!(x0[2:end] - x0[1:end-1], 0)
                    ind = findall(dx0 .> 1/52);
                    
                    c = 0
                    for in in ind
                        x0 = insert!(x0, in+c+1, x0[in+1])
                        y0 = insert!(y0, in+c+1, NaN)
                        c = c+1
                    end

                    bincenter.(df3[:,1])
                    lines!(ax, [2000, 2024], [oset, oset], color = :lightgray)

                    med_offset = @sprintf("%0.2f", med_offset)

                    lines!(ax, x0, y0 .+ oset, label = "d$u [$(med_offset)]")
                    #scatter!(x0, y0 .+ oset, label = "d$u [$(med_offset)]", markersize=5)

                    oset = oset - delta;
                end

                ax.yticks = (.6:-.2:-.6) .* delta

                Label(f[1, 1, Top()], "$sensor $delta_bin $vname of flight $(string(fun)) detector offset relative to $(dems[i]) [m]", valign = :bottom,
                    font = :bold,
                    padding = (0, 0, 5, 0),)

                leg = Legend(f[1,2], ax)

                ylims!(ax, oset, delta)
                xlims!(ax, minx, maxx)
                ax.xticks = minx:deltax:maxx

                #display(f); return
                
                fig_name =  "$(sensor)_$(delta_bin)_$(vname)_$(string(fun))_offset_to_$(dems[i]).pdf"
                save(joinpath(outFigFolder, fig_name), f)
            end
        end
    end
end
# rsync -r bylot:/home/gardnera/Documents/Figures/Altim/offsets/ /Users/gardnera/Documents/Figures/Altim/offsets/

f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1000, 700), title = "test")

ga = f[1, 1]# = GridLayout()
gb = f[2, 1]# = GridLayout()
gc = f[3, 1]# = GridLayout()
gd = f[1, 2]# = GridLayout()
axtop = Axis(ga[1,1])
axmiddle = Axis(gb)
axbottom = Axis(gc)
for i = eachindex(dems)
    if !isempty(offsets[i])
        println(i)
        valid = .!isnan.(offsets[i].forward)
        density!(axtop, offsets[i].forward[valid], boundary = (-50, 50), label = string(dems[i]))
        density!(axmiddle, offsets[i].right[valid], boundary = (-50, 50))
        density!(axbottom, offsets[i].above[valid], boundary = (-20, 20))
    end
end
leg = Legend(gd[1,1], axtop)

Label(ga[1, 3, Top()], "forward", valign = :bottom,
    font = :bold,
    padding = (0, 0, 5, 0))

Label(gb[1, 2, Top()], "right", valign = :bottom,
font = :bold,
padding = (0, 0, 5, 0))

Label(gc[1, 2, Top()], "above", valign = :bottom,
font = :bold,
padding = (0, 0, 5, 0))

f




