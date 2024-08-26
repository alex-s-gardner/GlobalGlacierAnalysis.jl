#using FileIO
#f = FileIO.load("/mnt/bylot-r3/data/binned/2deg/glacier_dh_best_cc_meanmadnorm3_v01_filled_ac_p1.jld2")
#f["dh_hyps"]["icesat2"]


# make regional plots
#if false

function plot_regional_dvdm(df ;
    # mission_featured = "synthesis"
    rgi = "rgi17",
    var0 = dm,
    mission_featured = "gedi",
    showmissions = true,
    )

    missions = unique(df.mission)

    # place mission_featured last for plotting 
    missions = missions[missions .!= mission_featured]
    missions = vcat(missions, mission_featured)

    title = Altim.rgi2label[rgi]

    index = findfirst((df.rgi.==rgi) .&  (df.var.==var0))
    ylabel = "$(Altim.units2label[df[index, "unit"]]) [$(df[index, "unit"])]"

    # make plots
    f = CairoMakie.Figure();
    ax = CairoMakie.Axis(f[1, 1];
        title,
        xticks = WilkinsonTicks(5),
        ylabel,
    )
    for mission in missions 
        index = (df.rgi .== rgi) .&  (df.var .== var0) .& (df.mission .== mission)
        if .!any(index)
            continue
        end
        mid = df_dmass[index, :mid][1] 
        low = df_dmass[index, :low][1] 
        high = df_dmass[index, :high][1]
        notnan = .!isnan.(mid)

        if any(notnan)
            if mission .== mission_featured
                CairoMakie.band!(decyear[notnan], low[notnan], high[notnan]; color=(Makie.wong_colors()[1], 0.6))
                CairoMakie.lines!(decyear[notnan], mid[notnan])
            elseif showmissions
                CairoMakie.band!(decyear[notnan], low[notnan], high[notnan]; color=(:lightgrey, 0.5))
                CairoMakie.lines!(decyear[notnan], mid[notnan], color=(:black, 0.3))
            end
        end
    end
    
    return f
end
 
#begin
    #df = deepcopy(df_dheight)
    df = deepcopy(df_dmass)
    rgis = reduce(vcat, (["rgi$i" for i in vcat(1:12, 16:19)], ["hma"]))

    color_pallet = Makie.Palette(:matter, [length(rgis)])

    showlines = false
    showmissions = true

    vars = unique(df_dmass.var)
    fontsize = 15
    #for var0 in vars
    begin
    var0 = first(vars)
                
        mission0 = "synthesis"

        # subset regions to include 
        index = (df.var .== var0) .& (df.mission .== mission0)

        reg_index = falses(size(index))
        for rgi in rgis
            reg_index = reg_index .| (df.rgi .== rgi)
        end
        index = index .& reg_index

        index0 = findfirst(index)
        title = "Glacier $(Altim.units2label[df[index0, "unit"]])"
        ylabel = "$(Altim.units2label[df[index0, "unit"]]) [$(df[index0, "unit"])]"

        if (df[index0, "unit"] .== "m") || (df[index0, "unit"] .== "m w.e.")
            delta_offset = -10
        else
            delta_offset = -200
        end

        dfX = df[index,:]

        sort!(dfX, :trend; rev=true)
        
        notnan= vec(any(.!isnan.(hcat(dfX.mid...)), dims = 2))

        # make plots
        xlims = extrema(decyear[notnan])
        xlims = (floor(xlims[1]), ceil(xlims[2]))

        f = CairoMakie.Figure(; backgroundcolor=RGBf(0.98, 0.98, 0.98), size=(500, 1000), fontsize);
        
        # initialize
        n = nrow(dfX)
        y2ticks = zeros(n)
        y2ticklabels = fill("", n)
        y2ticklablecolor = color_pallet;
        yticklables = [Altim.rgi2label[r] for r in dfX.rgi]

        ax1 = CairoMakie.Axis(f[1, 1];
            title,
            color = color_pallet,
            xticks = WilkinsonTicks(5),
            yticklabelrotation = 0 / (360/(2*pi)),
            ygridcolor  = :black,
            ygridwidth = 0.5,
            yticks = ((1:n)/n, string.((1:n)/n)),
            yticklabelcolor,
        )

        CairoMakie.xlims!(ax1, xlims...)

        offset = 0;
        for (i,r) in enumerate(eachrow(dfX))
        #i = 1; r = first(eachrow(dfX))

            println(r.rgi)
            o = offset - r.mid[notnan][1]
            CairoMakie.band!(ax1, decyear[notnan], r.low[notnan] .+ o, r.high[notnan] .+ o; color=(:lightgrey, 0.5))

            if showmissions
                for mission in missions
                    index = (df.var .== r.var) .& (df.mission .== mission) .& (df.rgi .== r.rgi) 
                    if any(index)
                    CairoMakie.lines!(ax1, decyear[notnan], df.mid[index][1][notnan] .+ o; color=(:black, 0.5))
                    end
                end
            end

            y = r.mid[notnan];
            ln = CairoMakie.lines!(ax1, decyear[notnan], y .+ o)
            
            y2ticklabels[i] = string(round(Int64, (y[end] - y[1]))) * " $(df[index0, "unit"])"
            y2ticks[i] = r.mid[notnan][end].+ o
            y2ticklablecolor[i] = ln.color.val;
            offset += delta_offset
        end
        
        CairoMakie.reset_limits!(ax1) # this is needed to be able to retrive limits
        minorgridspacing = -delta_offset/5
        ylim = round.(Int, ax1.yaxis.attributes.limits.val./minorgridspacing).*minorgridspacing
        yticks = 0:delta_offset:(offset-delta_offset)

        ax1.yticks = (yticks, yticklables)
        ax1.yticklabelcolor = y2ticklablecolor
        ax1.ygridvisible = showlines
 
        !showlines && CairoMakie.hidespines!(ax1)
        !showlines && CairoMakie.hidexdecorations!(ax1)


        ax2 = CairoMakie.Axis(f[1, 1]; 
            yaxisposition = :right,
            yticks = (y2ticks,y2ticklabels),
            yticklabelcolor,
        )

        foo = ax1.yaxis.attributes.limits.val
        ylims!(ax2,foo)

        #CairoMakie.linkyaxes!(ax1, ax2)
        minorgridspacing = -(delta_offset/5)

        a = ylim[1]:minorgridspacing:ylim[2]
        ax2.yminorticks = ylim[1]:minorgridspacing:ylim[2]

        ax2.yminorticksvisible = false
        ax2.yminorgridvisible = showlines
        ax2.yticksvisible = true
        ax2.ygridvisible = false
        CairoMakie.hidespines!(ax2)
        CairoMakie.hidexdecorations!(ax2)

        fname = joinpath("test.png")
        display(f)
        #save(fname, f)
    end
end