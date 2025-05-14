function plot_multiregion_dvdm(
    regions = regions;
    variables = ["dm", "dm_altim"], # last variable is plotted last
    units = "Gt",
    rgi_regions = setdiff(collect(dims(runs_rgi["dm_altim"], :rgi)), [13, 14, 15, 99]),
    showlines = false,
    fontsize = 15,
    cmap = :Dark2_4,
    region_order = nothing,
    ylims = nothing,
    title = nothing,
    palette = nothing,
    delta_offset = nothing,
    all_error_bounds = false,
    daterange = DateTime(2000, 1, 1):Month(1):DateTime(2025, 1, 1),
    numbersinylabel = false,

    )
 
    drgi = Dim{:rgi}(rgi_regions)
    dvarname = Dim{:varname}(variables)

    if isnothing(delta_offset)
        if (units .== "m") || (units .== "m w.e.")
            delta_offset = -10
        else
            delta_offset = -100
        end
    end

    total_mass_change = zeros(drgi)

    #total mass change of feature variable (i.e. last variable in variables)
    first_valid = findfirst(.!isnan.(regions[dvarname[end]][At(drgi[1]), minimum(daterange)..maximum(daterange), At(false)]))
    last_valid = findlast(.!isnan.(regions[dvarname[end]][At(drgi[1]), minimum(daterange)..maximum(daterange), At(false)]))

    if isempty(first_valid)
        error("No valid data for $dvarname[end] in $daterange")
    end

    for rgi in drgi
        a = regions[dvarname[end]][At(rgi), minimum(daterange)..maximum(daterange), At(false)][first_valid]
        b = regions[dvarname[end]][At(rgi), minimum(daterange)..maximum(daterange), At(false)][last_valid]
        total_mass_change[At(rgi)] = b - a
    end

    if isnothing(region_order)
        region_order = drgi[sortperm(collect(total_mass_change[:]); rev=true)]
    end

    total_mass_change = total_mass_change[At(rgi_regions)]

    # make plots
    xtick_delta = 5
    years = unique(year.(daterange))

    yticks = zeros(length(rgi_regions))
    xticks = (floor(Int64, minimum(years) ./ xtick_delta)*xtick_delta):xtick_delta:(ceil(Int64, maximum(years) ./ xtick_delta)*xtick_delta)

    xlims = (xticks.start, xticks.stop)

    # backgroundcolor=RGBf(0.98, 0.98, 0.98)
    f = CairoMakie.Figure(; size=(500, 750), fontsize)

    # initialize
    n = length(rgi_regions)
    (n == 1) ? (n = 2) : (n = n)

    if isnothing(palette)
        palette = (; color=Makie.resample_cmap(cmap, n))
    end

    y2ticks = zeros(n)
    y2ticklabels = fill("", n)
    y2ticklabelcolor = palette.color

    yticklabels = Altim.rginum2label.(region_order)

    if numbersinylabel
        yticklabels = ["$(Altim.rginum2label(id)) $(Altim.rginum2enclosed_alphanumerics[id])" for id in region_order]
    end

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
        ax1.titlecolor = RGBf(0.85, 0.85, 0.85) / 2
    end

    ax2 = CairoMakie.Axis(f[1, 1];
        yaxisposition=:right,
        yticks=((1:n) / n, string.((1:n) / n)),
        yticklabelcolor=y2ticklabelcolor,
    )

    CairoMakie.xlims!(ax1, xlims...)

    ymin = Inf
    ymax = -Inf

    region_offsets = zeros(region_order)
    
    # calculate offsets 
    varname = dvarname[end]
    valid_index = .!isnan.(regions[varname][At(region_order[1]), minimum(daterange)..maximum(daterange), At(false)])
    last = zeros(sum(valid_index))

    for rgi in region_order
        mid0 = regions[varname][At(rgi), minimum(daterange)..maximum(daterange), At(false)][valid_index]
        region_offsets[At(rgi)] = floor(min(minimum((last .- (mid0 .- mid0[1]))), 0) / delta_offset) * delta_offset .+ delta_offset .- mid0[1]

        last = mid0 .+ region_offsets[At(rgi)]
    end

    # plot error bounds for
    endi = length(dvarname)
    if all_error_bounds
        starti = 1
    else
        starti = length(dvarname)
    end
    
    for varname in dvarname[starti:endi]
        valid_index = .!isnan.(regions[varname][At(region_order[1]), minimum(daterange)..maximum(daterange), At(false)])
        last = zeros(sum(valid_index))

        for rgi in region_order
 
            mid0 = regions[varname][At(rgi), minimum(daterange)..maximum(daterange), At(false)][valid_index]
   
            if all(isnan.(mid0))
                continue
            end

            mid = mid0 .+ region_offsets[At(rgi)]
            low = mid .- regions[varname][At(rgi), minimum(daterange)..maximum(daterange), At(true)][valid_index]
            high = mid .+ regions[varname][At(rgi), minimum(daterange)..maximum(daterange), At(true)][valid_index]

            CairoMakie.band!(ax1, collect(Altim.decimalyear.(dims(mid0, :date))), collect(low), collect(high); color=(:grey, 0.4))

            ymin = min(ymin, minimum(low))
            ymax = max(ymax, maximum(high))
    
        end
    end

    # plot lines over error bounds
    for varname in dvarname
        #varname = "dm_altim"
        valid_index = .!isnan.(regions[varname][At(region_order[1]), minimum(daterange)..maximum(daterange), At(false)])
        last = zeros(sum(valid_index))

        for (i, rgi) in enumerate(region_order)
           
            #rgi = region_order[1]
            mid0 = regions[varname][At(rgi), minimum(daterange)..maximum(daterange), At(false)][valid_index] .+ region_offsets[At(rgi)]

            if all(isnan.(mid0))
                continue
            end
           
            if varname == dvarname[end]
                ln = CairoMakie.lines!(ax1, collect(Altim.decimalyear.(dims(mid0, :date))), collect(mid0))

                y2ticklabels[i] = string(round(Int64, (mid0[end] - mid0[1]))) * " $(units)"
                yticks[i] = mid0[1]
                y2ticks[i] = mid0[end] 
                y2ticklabelcolor[i] = ln.color.val
            else
                 ln = CairoMakie.lines!(ax1,collect(Altim.decimalyear.(dims(mid0, :date))), collect(mid0); color = (:black, 0.3))
            end

        end
    end

    if isnothing(ylims)
        ylims = [ymin, ymax]
    end

    CairoMakie.ylims!(ax1, minimum(ylims), maximum(ylims))
    CairoMakie.reset_limits!(ax1) # this is needed to be able to retrive limits

    minorgridspacing = -delta_offset / 5
    ylim = round.(Int, ax1.yaxis.attributes.limits.val ./ minorgridspacing) .* minorgridspacing
    #yticks = collect(region_offsets[At(collect(region_order))])

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

    return f, region_order, ylims
end


function show_error_bar_table(regional_results; cols2display = 3)
    # dump table to console
    colnames = names(regional_results)
    last_set = ceil(Int, (length(colnames) - 2) / (cols2display+1) - 1)
    for i = 0:last_set
        printstyled("\n------------------------------------------------------- $(colnames[3 + (4 * i)])-------------------------------------------------------\n", color=:yellow)
        if (i == last_set) && (length(colnames) > (last_set + 4))
            show(regional_results[:, [1:2; (3 .+ (4*i)):end]], allrows=true, allcols=true)
        else
            show(regional_results[:, [1:2; (3:(2+cols2display)) .+ (4 * i)]], allrows=true, allcols=true)
        end
    end
end

