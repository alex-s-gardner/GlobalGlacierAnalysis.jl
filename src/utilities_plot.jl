
function plot_regional_dvdm(df;
    # featured_mission = "synthesis"
    rgi="rgi17",
    variable="dm",
    missions, # missions are listed in the oder that they will be plotted
    colors,
    xticks,
)

    dates = DataFrames.metadata(df, "date")
    decyear = Altim.decimalyear.(dates)

    title = Altim.rgi2label[rgi]
    
    # index into rows to plot
    index = falses(nrow(df))
    for mission in missions
         index =  index .| (df.mission .== mission)
    end
    index = index .& (df.rgi .== rgi) .& (df.var .== variable)
    if .!any(index)
        return nothing
    end

    dfX = df[index,:]

    ylabel = "$(Altim.units2label[df[findfirst(index), "unit"]]) [$(df[findfirst(index), "unit"])]"

    # make plots
    f = CairoMakie.Figure()
    ax = CairoMakie.Axis(f[1, 1];
        title,
        ylabel,
        xticks,
    )
    for mission in missions
        index_mission = findfirst(missions .== mission)
        if isnothing(findfirst(dfX.mission .== mission))
            continue
        end
        r = dfX[findfirst(dfX.mission .== mission),:]
        notnan = .!isnan.(r.mid)

        if any(notnan)
            CairoMakie.band!(ax,decyear[notnan], r.low[notnan], r.high[notnan]; color=(colors[index_mission], 0.1))
            CairoMakie.lines!(ax, decyear[notnan], r.mid[notnan]; color=colors[index_mission], label = Altim.mission2label[mission])
        end


    end
    axislegend(ax; framevisible = false)
    CairoMakie.xlims!(extrema(xticks)...)

    return f
end

function plot_multiregion_dvdm(df;
    variable="dm",
    featured_mission="synthesis",
    regions=reduce(vcat, (["rgi$i" for i in vcat(1:12, 16:19)], ["hma"])),
    showlines=false,
    showmissions=true,
    fontsize=15,
    cmap=:Dark2_4,
    regions_ordered=false,
    region_offsets=nothing,
    ylims=nothing,
    title = nothing,
    palette = nothing,
    delta_offset = nothing,
)

    missions = unique(df.mission)
    dates = DataFrames.metadata(df, "date")
    decyear = Altim.decimalyear.(dates)

    # place featured_mission last for plotting 
    missions = missions[missions.!=featured_mission]
    missions = vcat(missions, featured_mission)

    # subset regions to include 
    index = (df.var .== variable) .& (df.mission .== featured_mission)

    if !any(index)
        printstyled("NO (variable = $variable && mission = $featured_mission)\n"; color=:red)
        error()
    end

    reg_index = falses(size(index))
    for region in regions
        reg_index = reg_index .| (df.rgi .== region)
    end
    index = index .& reg_index

    index0 = findfirst(index)
    # title = "Glacier $(Altim.units2label[df[index0, "unit"]])"

    if isnothing(delta_offset)
        if (df[index0, "unit"] .== "m") || (df[index0, "unit"] .== "m w.e.")
            delta_offset = -10
        else
            delta_offset = -100
        end
    end

    dfX = df[index, :]
    notnan = vec(any(.!isnan.(hcat(dfX.mid...)), dims=2))

    #total mass change
    dfX[!, "total_change"] .= 0.0
    for r in eachrow(dfX)
        foo = r.mid[notnan]
        r.total_change = foo[end] - foo[1]
    end

    if !regions_ordered
        dfX = sort!(dfX, :total_change; rev=true)
    else
        # order same as input regions
        p = findfirst.(isequal.(regions), (dfX.rgi,))
        dfX = dfX[p, :]
    end

    # make plots
    xtick_delta = 4
    xticks = (floor(Int64, minimum(decyear[notnan]) ./ xtick_delta)*xtick_delta):xtick_delta:(ceil(Int64, maximum(decyear[notnan]) ./ xtick_delta)*xtick_delta)

    xlims = (xticks.start, xticks.stop)

    # backgroundcolor=RGBf(0.98, 0.98, 0.98)
    f = CairoMakie.Figure(; size=(500, 750), fontsize)

    # initialize
    n = nrow(dfX)
    
    if isnothing(palette)
        palette = (; color=Makie.resample_cmap(cmap, n))
    end

    y2ticks = zeros(n)
    y2ticklabels = fill("", n)
    y2ticklabelcolor = palette.color
    yticklabels = [Altim.rgi2label[rginum2txt[r]] for r in dfX.rgi]

    # this is a hack... axes need to be defined early or things break
    ax1 = CairoMakie.Axis(f[1, 1];
        #title,
        palette,
        xticks=WilkinsonTicks(5),
        yticklabelrotation=0 / (360 / (2 * pi)),
        ygridcolor=:black,
        ygridwidth=0.5,
        yticks=((1:n) / n, string.((1:n) / n)),
        yticklabelcolor=y2ticklabelcolor,
    )

    if .!isnothing(title)
        ax1.title = title
        ax1.titlecolor = RGBf(0.85, 0.85, 0.85)/2
    end

    ax2 = CairoMakie.Axis(f[1, 1];
        yaxisposition=:right,
        yticks=((1:n) / n, string.((1:n) / n)),
        yticklabelcolor=y2ticklabelcolor,
    )

    CairoMakie.xlims!(ax1, xlims...)

    ymin = Inf
    ymax = -Inf

    last = zeros(sum(notnan))

    if !isnothing(region_offsets)
        dfX[!, :offset] = region_offsets
    else
        dfX[!, :offset] .= 0.0
    end

    # plot error bounds first so all other lines are overlain
    first_iter = true
    for r in eachrow(dfX)
        #r = first(eachrow(dfX))
        ctr = r.mid[notnan][1]
        if first_iter || !isnothing(region_offsets)
            first_iter = false
        else
            r.offset = floor(min(minimum((last .- (r.mid[notnan] .- ctr))), 0) / delta_offset) * delta_offset .+ delta_offset
        end

        mid = r.mid[notnan] .+ r.offset .- ctr
        low = r.low[notnan] .+ r.offset .- ctr
        high = r.high[notnan] .+ r.offset .- ctr

        CairoMakie.band!(ax1, decyear[notnan], low, high; color=(:grey, 0.3))

        last = mid
        ymin = min(ymin, minimum(low))
        ymax = max(ymax, maximum(high))
    end

    # plot lines over error bounds
    for (i, r) in enumerate(eachrow(dfX))
        o = r.offset - r.mid[notnan][1]

        if showmissions
            for mission in missions
                index = (df.var .== r.var) .& (df.mission .== mission) .& (df.rgi .== r.rgi)
                if any(index)
                    CairoMakie.lines!(ax1, decyear[notnan], df.mid[index][1][notnan] .+ o; color=(:black, 0.5))
                end
            end
        end

        y = r.mid[notnan]
        ln = CairoMakie.lines!(ax1, decyear[notnan], y .+ o)

        y2ticklabels[i] = string(round(Int64, (y[end] - y[1]))) * " $(df[index0, "unit"])"
        y2ticks[i] = r.mid[notnan][end] .+ o
        y2ticklabelcolor[i] = ln.color.val
    end

    if isnothing(ylims)
        ylims = [ymin, ymax]
    end

    CairoMakie.ylims!(minimum(ylims), maximum(ylims))
    CairoMakie.reset_limits!(ax1) # this is needed to be able to retrive limits

    minorgridspacing = -delta_offset / 5
    ylim = round.(Int, ax1.yaxis.attributes.limits.val ./ minorgridspacing) .* minorgridspacing
    yticks = dfX.offset

    ax1.yticks = (yticks, yticklabels)
    ax1.yticklabelcolor = y2ticklabelcolor
    ax1.ygridvisible = showlines

    #!showlines && CairoMakie.hidexdecorations!(ax1,)
    try
        ax2.yticks = (y2ticks, y2ticklabels)
    catch e
        printstyled("NOTE: if you end up here there is a bug with CairoMakie in which a cryptic error is thrown when the ytick positions exceed the yaxis limits\n"; color=:red)
        rethrow(e)
    end


    ax2.yticklabelcolor = y2ticklabelcolor

    foo = ax1.yaxis.attributes.limits.val
    CairoMakie.ylims!(ax2, ylims)
    CairoMakie.ylims!(ax2, ylims)

    minorgridspacing = -(delta_offset / 5)

    ax2.yminorticks = ylim[1]:minorgridspacing:ylim[2]

    !showlines && CairoMakie.hidespines!(ax1, :t, :r, :l)
    CairoMakie.hidespines!(ax2)
    ax2.yminorticksvisible = false
    ax2.yminorgridvisible = showlines
    ax2.yticksvisible = false
    ax2.ygridvisible = false

    CairoMakie.hidexdecorations!(ax2)

    xtickcolor = RGBf(0.85, 0.85, 0.85)
    ax1.xaxisposition = :bottom
    ax1.xgridvisible = true
    ax1.xgridcolor = xtickcolor
    ax1.xticklabelcolor = xtickcolor ./ 2
    ax1.xticks = (Float64.(xticks), string.(xticks))
    ax1.xtickcolor = xtickcolor
    ax1.xlabelvisible = true
    ax1.yticksvisible = false
    ax1.xticksvisible = true
    ax1.bottomspinecolor = xtickcolor

    CairoMakie.ylims!(f.content[1], ylims)

    return f, dfX[:, :rgi], dfX[:, :offset], ylims
end





"""
    geotile_filled_dv_reg(dh, nobs0, geotiles, reg)

Calculate volume change per geotile and combine into regional estimates.

# Arguments
- `dh`: Dictionary mapping mission names to arrays of elevation changes
- `nobs0`: Dictionary mapping mission names to arrays of observation counts
- `geotiles`: DataFrame containing geotile information including area and region assignments
- `reg`: Vector of region identifiers

# Returns
- `dv_reg`: Array of volume changes per mission, region and date
- `nobs_reg`: Array of observation counts per mission, region and date  
- `area_reg`: Vector of total area per region

# Details
For each geotile, calculates volume change by multiplying elevation change by area.
Then aggregates geotile volume changes into regional totals based on region assignments.
Handles NaN values and zero-filling appropriately.
"""
function geotile_filled_dv_reg(dh, nobs0, geotiles, reg)

    vars = keys(dh)
    drgi = Dim{:rgi}(reg)
    dgeotile = dims(dh[first(vars)], :geotile)
    ddate = dims(dh[first(vars)], :date)
    dheight = dims(dh[first(vars)], :height)
    dmission = Dim{:mission}(collect(vars))

    # per geotile volume change
    dv2 = Dict()
    nobs2 = Dict()
    for var0 in vars
        dv2[var0] = fill(NaN, dgeotile, ddate)
        nobs2[var0] = fill(0, dgeotile, ddate)
    end

    for geotile in eachrow(geotiles)
        #geotile = geotiles[604,:]

        area = ones(length(ddate), 1) * hcat(geotile.area_km2)'

        for var0 in vars
            v0 = dh[var0][At(geotile.id), :, :] ./ 1000 .* area
            v0[isnan.(v0)] .= 0
            dv2[var0][At(geotile.id), :] = sum(v0, dims=:height)
            nobs2[var0][At(geotile.id), :] = sum(nobs0[var0][At(geotile.id), :, :], dims=:height)
        end
    end

    # combine RGI regions
    dv_reg = fill(NaN, dmission, drgi, ddate)
    nobs_reg = fill(0, dmission, drgi, ddate)
    area_reg = fill(0.0, drgi)

    for rgi in drgi
        #rgi = first(drgi)
        index = geotiles[:, rgi] .> 0.0
        area_reg[At(rgi)] = sum(reduce(vcat, geotiles.area_km2[index]))
        for mission in dmission
            v0 = dv2[mission][index, :]
            v0[isnan.(v0)] .= 0
            dv_reg[At(mission), At(rgi), :] = sum(v0, dims=:geotile)
            nobs_reg[At(mission), At(rgi), :] = sum(nobs2[mission][index, :], dims=:geotile)
        end
    end

    for mission in dmission
        v0 = dv_reg[At(mission), :, :]
        set2nan = all(v0 .== 0, dims=:rgi)
        v0[:, vec(set2nan)] .= NaN
        dv_reg[At(mission), :, :] = v0
    end

    return dv_reg, nobs_reg, area_reg
end