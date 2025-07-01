"""
    plot_unbinned_height_anomalies(datetime, dh; title="")

Create a visualization of raw height anomalies and their monthly medians over time.

# Arguments
- `datetime`: Array of DateTime values representing measurement times
- `dh`: Array of height anomaly measurements [m]
- `title`: Optional title for the plot (default: "")

# Returns
- A Makie Figure object containing:
  - Top panel: Height anomalies with monthly medians, scaled to show detail
  - Bottom panel: Same data with expanded y-axis to show full range of values

# Description
Creates a two-panel plot showing height anomaly measurements over time. The top panel
uses a y-axis scaled to the 95th percentile of anomalies, while the bottom panel
shows the full range. Both panels include raw measurements and monthly median values.
"""
function plot_unbinned_height_anomalies(datetime, dh; title="")
    mm = 3.7795275590551176
    figure_width = 183mm

    sample_rate = ceil(Int, length(dh) / 1E6)

    df = binstats(DataFrame(Ti=decimalyear.(datetime), dh=dh), :Ti, decimalyear.(minimum(datetime):Month(1):maximum(datetime)), :dh, col_function=[median])

    f = Figure(size=(figure_width, figure_width * 2 / 3))

    ymax = quantile(abs.(dh[1:sample_rate:end]), 0.99)
    ax1 = CairoMakie.Axis(f[1:2, 1]; title, ylabel="height anomaly [m]")
    CairoMakie.scatter!(ax1, datetime[1:sample_rate:end], dh[1:sample_rate:end]; color=(:skyblue2, 0.2), label="raw dh")
    # add zero line
    CairoMakie.lines!(ax1, [minimum(datetime), maximum(datetime)], [0., 0.]; color=:snow4, linewidth=0.5)
    # add monthly medians
    CairoMakie.lines!(ax1, convert.(typeof(minimum(datetime)), decimalyear2datetime.(df.Ti)), df.dh_median; color=:black, label="monthly median")
    CairoMakie.ylims!(ax1, (-ymax, ymax))
    axislegend(ax1, position=:rt) # orientation=:horizontal, fontsize=10

    ax2 = CairoMakie.Axis(f[3, 1]; ylabel="height anomaly [m]")
    CairoMakie.scatter!(ax2, datetime[1:sample_rate:end], dh[1:sample_rate:end]; color=:skyblue2, label="raw")
    # add zero line
    CairoMakie.lines!(ax2, [minimum(datetime), maximum(datetime)], [0., 0.]; color=(:snow4, 0.2), linewidth=0.5)
    # add monthly medians
    CairoMakie.lines!(ax2, convert.(typeof(datetime[1]), decimalyear2datetime.(df.Ti)), df.dh_median; color=:black, label="monthly median")

    ymax = ceil(Int, maximum(abs.(df.dh_median)) / 5) * 5 + 5
    CairoMakie.ylims!(ax2, (-ymax, ymax))

    return f
end

"""
    plot_curvature(bin_center, dh_obs, dh_cor, nrow, mission, geotile)

Create a visualization showing the relationship between surface curvature and height anomalies.

# Arguments
- `bin_center`: Array of curvature bin centers [cm⁻¹]
- `dh_obs`: Array of observed height anomalies [m]
- `dh_cor`: Array of modeled height anomalies [m]
- `nrow`: Array of observation counts per bin
- `mission`: Mission identifier string
- `geotile`: Geotile object containing location information

# Returns
- A Makie Figure object containing:
  - Top panel: Height anomalies vs curvature showing observations, model fit, and corrected values
  - Bottom panel: Histogram of observation counts per curvature bin

# Description
Creates a two-panel plot showing how height anomalies vary with surface curvature. The top panel
displays raw observations, the model fit, and corrected values. The bottom panel shows the
distribution of observations across curvature bins.
"""
function plot_curvature(bin_center, dh_obs, dh_cor, nrow; title = "")
    mm = 3.7795275590551176
    figure_width = 183mm
    
    valid = .!ismissing.(dh_obs)
    valid_range, = validrange(valid)
   

    f = Figure(size=(figure_width * 1.5, figure_width))

    ax1 = CairoMakie.Axis(f[1:4, 1]; title, ylabel="height anomaly [m]")
    ax2 = CairoMakie.Axis(f[5:6, 1], ylabel="count [×1000]", xlabel="curvature [cm⁻¹]")

    hidexdecorations!(ax1; grid=false, minorgrid=false)

    plot!(ax1, bin_center[valid], dh_obs[valid]; label="observation")
    lines!(ax1, bin_center[valid_range], dh_cor[valid_range]; label="model")
    plot!(ax1, bin_center[valid], dh_obs[valid] .- dh_cor[valid]; label="corrected")
    barplot!(ax2, bin_center[valid], nrow[valid] / 1000)
    axislegend(ax1, framevisible=false, position=:lt)

    return f
end





## ---------------------------- elevation-time matrix plotting ----------------------------
"""
    plot_elevation_time(dh0; colorrange=(-20, 20))

Create a heatmap visualization of elevation changes over time.

# Arguments
- `dh0`: A dimensional array containing elevation change data with time and elevation dimensions
- `colorrange`: Tuple specifying the color scale range for the heatmap (default: (-20, 20))

# Returns
- A Makie Figure object containing the elevation-time heatmap

# Description
This function creates a heatmap where the x-axis represents time and y-axis represents elevation,
displaying elevation changes using a balanced color scheme. Only valid (non-NaN) data points are plotted.
The heatmap shows height anomalies [m] over time with elevation bins on the vertical axis.
"""
function plot_elevation_time(dh0; colorrange=(-20, 20))
    (valid_rows, valid_cols) = validrange(.!isnan.(dh0))

    f = Figure()
    ax = CairoMakie.Axis(f[1, 1], ylabel="height anomaly [m]")
    heatmap!(ax,
        dh0[valid_rows, valid_cols];
        colormap=Makie.Reverse(:balance), 
        colorrange,
        label=dh0.refdims[1][1]
    )

    #p.axis.xticklabel = Year.(p.axis.xaxis.tickvalues)
    return f
end

"""
    plot_elevation_time!(ax, dh0; colorrange=(-20, 20), xaxis_label=true, yaxis_label=true)

Create a heatmap visualization of elevation changes over time on an existing axis.

# Arguments
- `ax`: The axis object to plot on
- `dh0`: A dimensional array containing elevation change data with time and elevation dimensions
- `colorrange`: Tuple specifying the color scale range for the heatmap (default: (-20, 20))
- `xaxis_label`: Boolean to show x-axis label (default: true)
- `yaxis_label`: Boolean to show y-axis label (default: true)

# Returns
- The modified axis object with the heatmap plot

# Description
This function creates a heatmap on an existing axis where the x-axis represents time and y-axis represents
elevation, displaying elevation changes using a balanced color scheme. Only valid (non-NaN) data points are plotted.
This is the in-place version of `plot_elevation_time()` for use in multi-panel figures.
"""
function plot_elevation_time!(ax, dh0; colorrange = (-20,20), xaxis_label = true, yaxis_label = true)
    (valid_rows, valid_cols) = validrange(.!isnan.(dh0))

    if length(dh0.refdims) > 0
        label=dh0.refdims[1][1]
    else
        label=""
    end
    ax = heatmap!(ax,
        dh0[valid_rows, valid_cols];
        colormap=Makie.Reverse(:balance), 
        colorrange,
        label,
    )
    return ax
end

"""
    plot_elevation_time_multimission(dh; colorrange=(-20, 20), linkaxes=true, colorbar_label="height anomaly", mission_order=nothing, xtickspacing=5)

Create a multi-panel figure showing elevation changes over time for multiple altimetry missions.

# Arguments
- `dh`: Dictionary containing elevation change data for different missions/datasets
- `colorrange`: Tuple specifying the color scale range for all heatmaps (default: (-20, 20))
- `linkaxes`: Boolean to link x and y axes across subplots (default: true)
- `colorbar_label`: Label for the shared colorbar (default: "height anomaly")
- `mission_order`: Vector specifying the order of missions to display (default: nothing, uses all missions)
- `xtickspacing`: Spacing between x-axis ticks (default: 5)

# Returns
- A Figure object containing heatmap subplots for each mission with a shared colorbar

# Description
This function creates a multi-panel visualization where each subplot shows elevation changes
over time for a different altimetry mission or dataset. The subplots are arranged in a 2-column
grid layout, with each mission labeled at the top of its respective panel. A shared colorbar
on the right side shows the height anomaly scale in meters.
"""
function plot_elevation_time_multimission(dh; colorrange=(-20, 20), linkaxes=true, colorbar_label="height anomaly", mission_order = nothing, xtickspacing = 5)
   
    mm = 3.7795275590551176
    figure_width = 183mm

    missions = mission_proper_name.(collect(keys(dh)))
    dheight = dims(dh[first(keys(dh))], :height)
    Δheight = val(dheight)[2] - val(dheight)[1]
    ddate = dims(dh[first(keys(dh))], :date)
    
    valid = falses(size(dh[first(keys(dh))]))
    for mission in keys(dh)
        valid = valid .| .!isnan.(dh[mission])
    end
    date_range, = validrange(valid)

    xlims = (floor(Int, minimum(ddate[date_range] ./ xtickspacing)) * xtickspacing, ceil(Int, maximum(ddate[date_range] ./ xtickspacing)) * xtickspacing)

    f = Figure(size=(figure_width * 1.5, figure_width))

    if isnothing(mission_order)
        mission_order = collect(keys(dh))
    end

    height_min = Inf
    height_max = 0
    ax = []
    for (i, mission) in enumerate(mission_order)

        r = ceil(Int, i / 2)
        c = 2 - mod(i, 2)

        push!(ax, CairoMakie.Axis(f[r, c]))

        if !in(mission, keys(dh))
            continue
        end

        valid = .!isnan.(dh[mission])
        _, height_range = validrange(valid)
        height_min = min(height_min, minimum(parent(dheight[height_range])))
        height_max = max(height_max, maximum(parent(dheight[height_range])))

        plot_elevation_time!(ax[i], dh[mission]; colorrange)

        ax[i].title = "$(mission_proper_name(mission))"

        #Label(f[r, c, Top()], "$(missions[i]) - height anomaly [m]";
        #    valign=:bottom,
        #    padding=(0, 0, 5, 0))

        if (i != 1) && linkaxes
            linkyaxes!(ax[1], ax[i])
            linkxaxes!(ax[1], ax[i])
        end
        if r == 1
            ax[i].xticklabelsvisible = false
            ax[i].xticksvisible = false
        end

        if c == 2
            ax[i].ylabelvisible = false
            ax[i].yticklabelsvisible = false
            ax[i].yticksvisible = false
        end

        ax[i].ylabel = "elevation"
        ax[i].ytickformat = values -> ["$(round(Int,value))m" for value in values]


        ax[i].xticks = xlims[1]:xtickspacing:xlims[2]-xtickspacing
        xlims!(ax[i], xlims)
    end

    height_min = max(0.0, (floor(height_min / Δheight) * Δheight) - Δheight)
    height_max = ceil(height_max / Δheight) * Δheight .+ Δheight

    for (i, mission) in enumerate(mission_order)
        ylims!(ax[i], (height_min, height_max))
    end

    Colorbar(f[3, 1:2], f.content[1].scene.plots[1], vertical=false, label=colorbar_label, flipaxis=false, tickformat=values -> ["$(round(Int,value))m" for value in values])

    return f
end

"""
    plot_hypsometry!(ax, area_km2)

Plot glacier area-elevation distribution (hypsometry) on a given axis.

# Arguments
- `ax`: Makie axis to plot on
- `area_km2`: Area data for each elevation bin [km²]

# Returns
- Modified axis with hypsometry plot

# Description
Creates a polygon plot showing the distribution of glacier area across elevation bins.
The plot displays area on the x-axis and elevation on the y-axis with a skyblue fill
and black outline.
"""
function plot_hypsometry!(ax, area_km2)
    ax.xtickformat=values -> ["$(round(Int,value))km²" for value in values]
    ax.title="glacier area-elevation distribution"

    dheight = dims(area_km2, :height)
    Δheight = val(dheight)[2] - val(dheight)[1]

    pts0 = Point2f.(collect(area_km2), parent(val(dheight)))
    pts = similar(pts0, (length(pts0) * 2)+2)

    # set starting and ending points to zero area
    pts[1] = Point2f(0, dheight[1])
    pts[end] = Point2f(0, dheight[end])

    for i in eachindex(pts0)
        n = i * 2
        pts[n] = pts0[i]
        if i < length(pts0)
            pts[n+1] = (pts0[i+1][1], pts0[i][2])
        else
            pts[n+1] = pts0[i]
        end
    end

    poly!(ax, pts, color=:skyblue2, strokecolor=:black, strokewidth=1)
    xlims!(ax, (0, ceil(maximum(area_km2)*1.05 / Δheight) * Δheight))
    ax.ygridvisible = false
    ax.xgridvisible = false
    ax.yticklabelsvisible = false
    hidespines!(ax, :t, :r);
    ax.yticksvisible = false

    return ax
end

"""
    plot_area_average_height_anomaly!(ax, dh0, area_km2; cmap=:thermal, mission_color_width=nothing, mission_order=nothing, xtickspacing=5)

Plot area-averaged height anomalies for multiple missions on a given axis.

# Arguments
- `ax`: Makie axis to plot on
- `dh0`: Dictionary of dimensional arrays containing height anomaly data for each mission
- `area_km2`: Area data for each geotile [km²]
- `cmap`: Colormap for mission lines (default: :thermal)
- `mission_color_width`: Dictionary specifying custom colors and linewidths for missions
- `mission_order`: Array specifying order of missions to display
- `xtickspacing`: Spacing between x-axis ticks (default: 5)

# Returns
- Modified axis with area-averaged height anomaly plots

# Description
Plots time series of area-averaged height anomalies for multiple missions. Each mission is
displayed as a line with the median value shown in the legend. Supports custom styling
and mission ordering.
"""
function plot_area_average_height_anomaly!(ax, dh0, area_km2; cmap=:thermal, mission_color_width = nothing, mission_order = nothing, xtickspacing = 5)
    # mission_color_width = Dict("synthesis" => (color = :black, linewidth = 1))
    ax.yaxisposition=:right
    ax.ytickformat=values -> ["$(round(Int,value))m" for value in values]
    ax.title="area-averaged height anomaly"

    # ----------- area average height anomaly  ------------    
    ddate = dims(dh0[first(keys(dh0))], :date)
    x = collect(val(ddate))
    clrs = Makie.resample_cmap(cmap, length(keys(dh0))+1)

    valid = falses(size(dh0[first(keys(dh0))]))
    for mission in keys(dh0)
        valid = valid .| .!isnan.(dh0[mission])
    end
    date_range, = validrange(valid)

    if isnothing(mission_order)
        mission_order = collect(keys(dh0))
    end

    for (i, mission) in enumerate(mission_order)
        y =  collect(dh_area_average(dh0[mission], area_km2))
        y_median = round(median(y[.!isnan.(y)]), digits=1)

        if !in(mission, keys(dh0))
            continue
        end

        if isnothing(mission_color_width) || (!in(mission, keys(mission_color_width)))
            lines!(ax, x, y; label="$(mission_proper_name(mission)): $(y_median)m", color=clrs[i])
        else
            lines!(ax, x, y; label="$(mission_proper_name(mission)): $(y_median)m", color = mission_color_width[mission].color, linewidth = mission_color_width[mission].linewidth)
        end
    end

    xlims = (floor(Int, minimum(ddate[date_range] ./ xtickspacing)) * xtickspacing, ceil(Int, maximum(ddate[date_range] ./ xtickspacing)) * xtickspacing)
    ax.xticks = xlims[1]:xtickspacing:xlims[2]
    xlims!(ax, xlims)

    axislegend(ax, position=:rt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false

    #colgap!(f.layout, 10)
    #rowgap!(f.layout, 10)
    #hidespines!(ax, :t, :r)
    return ax
end

"""
    plot_elevation_time_multimission_geotiles(dh; geotiles2plot, area_km2, colorrange, label, colorbar_label, mission_order, hypsometry, area_averaged, plots_show, plots_save, plot_save_path_prefix, plot_save_format, cmap, xtickspacing)

Create multi-mission elevation-time plots for specified geotiles.

# Arguments
- `dh`: Dictionary of dimensional arrays containing elevation change data for each mission
- `geotiles2plot`: Array of geotile identifiers to plot (default: ["lat[+30+32]lon[+078+080]"])
- `area_km2`: Area data for each geotile [km²]
- `colorrange`: Color scale range for heatmaps (default: (-20, 20))
- `label`: Optional title label for the figure
- `colorbar_label`: Label for colorbar (default: "height anomaly")
- `mission_order`: Array specifying order of missions to display
- `hypsometry`: Whether to include hypsometry plot (default: false)
- `area_averaged`: Whether to include area-averaged height anomaly plot (default: false)
- `plots_show`: Whether to display plots (default: false)
- `plots_save`: Whether to save plots (default: false)
- `plot_save_path_prefix`: Prefix for saved plot filenames
- `plot_save_format`: File format for saved plots (default: ".png")
- `cmap`: Colormap for plots (default: :thermal)
- `xtickspacing`: Spacing between x-axis ticks (default: 5)

# Returns
- Makie Figure object containing the multi-mission elevation-time visualization

# Description
Creates comprehensive multi-mission elevation-time plots for specified geotiles. The function
can include hypsometry and area-averaged height anomaly panels. Supports both display and
file saving options.
"""
function plot_elevation_time_multimission_geotiles(
    dh;
    geotiles2plot=["lat[+30+32]lon[+078+080]"],
    area_km2 = nothing,
    colorrange=(-20, 20), 
    label = nothing, 
    colorbar_label = "height anomaly",
    mission_order = nothing,
    hypsometry = false, 
    area_averaged = false, 
    plots_show = false, 
    plots_save = false, 
    plot_save_path_prefix = "",
    plot_save_format = ".png",
    cmap=:thermal,
    xtickspacing = 5,
    )

    if !isnothing(geotiles2plot)
        for geotile2plot in geotiles2plot
            # need to use copy as the date lookup is modified
            dh0 = Dict() 
            for mission in keys(dh)
                dh0[mission] =dh[mission][geotile=At(geotile2plot)]
                if all(isnan.(dh0[mission]))
                    #println("all NaNs for $mission")
                    delete!(dh0, mission)
                end
            end
            
            if !isnothing(mission_order)
                mission_order = intersect(collect(keys(dh0)), mission_order)
            end

            # In Makie you are not able to speci
            (ddate, dheight) = dims(dh0[first(keys(dh0))])
            decyear = decimalyear.(collect(ddate))
            ddate = Dim{:date}(decyear)
            for mission in keys(dh0)
                dh0[mission] = DimArray(dh0[mission][:,:], (ddate, dheight))
            end

            f = plot_elevation_time_multimission(dh0; colorrange, colorbar_label, mission_order, xtickspacing)
            
            if hypsometry
                ax13 = CairoMakie.Axis(f[1, 3])
                plot_hypsometry!(ax13, area_km2[geotile=At(geotile2plot)])
                ylims!(ax13, f.content[1, 1].yaxis.attributes.limits.val)
            end

            if area_averaged
                ax23 = CairoMakie.Axis(f[2, 3])
                plot_area_average_height_anomaly!(ax23, dh0, area_km2[geotile=At(geotile2plot)]; mission_order, cmap, xtickspacing)
            end

            isnothing(label) ? nothing : Label(f[0, 1:2], label)

            plots_show ? display(f) : nothing

            if plots_save
                fname = plot_save_path_prefix * "_$(geotile2plot)$(plot_save_format)"
                CairoMakie.save(fname, f)
            end

            return f
        end
    end
end

"""
    plot_amplitude_correction(model_ref, model0, delta, mission, p0, p_ref)

Create a three-panel visualization comparing amplitude correction models.

# Arguments
- `model_ref`: Reference model height anomalies [m]
- `model0`: Original model height anomalies [m] 
- `delta`: Difference between models [m]
- `mission`: Mission identifier string
- `p0`: Original amplitude correction parameters
- `p_ref`: Reference amplitude correction parameters

# Returns
- A Makie Figure object with three heatmaps showing the comparison
"""
function plot_amplitude_correction(model_ref, model0, delta, mission, p0, p_ref)

    mission = mission_proper_names(mission)

    max_dh = ceil(Int, max(maximum(abs.(model_ref)), maximum(abs.(model0))))
    colorrange = (-max_dh, max_dh)
    f = Figure(size=(2000, 600))
    ax = map(1:2:6) do i
        CairoMakie.Axis(f[1, i])
    end

    hm = [heatmap!(ax[i], dh; colorrange) for (i, dh) in enumerate([model_ref, model0, delta])]

    CairoMakie.Axis(f[1, 1], title="ICESat-2 height anomaly [m]")
    CairoMakie.Axis(f[1, 3], title="$mission height anomaly [m]")
    CairoMakie.Axis(f[1, 5], title="delta height anomaly [m]")

    Colorbar(f[1, 2], hm[2])
    Colorbar(f[1, 4], hm[3])
    Colorbar(f[1, 6], hm[1])

    label_text = "$(mission) amplitude correction [parameters: old = $(round.(p0, sigdigits=2)), new = $(round.(p_ref, sigdigits=2))]"
    Label(f[0, :], label_text; fontsize=20)

    colgap!(f.layout, 10)
    rowgap!(f.layout, 10)
    return f
end



"""
    plot_area_average_height_anomaly_with_error(dh_area_average_median, dh_area_average_error; cmap=:thermal, mission_color_width=nothing, mission_order=nothing, xtickspacing=5, median_in_label=false)

Create a plot showing area-averaged height anomalies with error bands.

# Arguments
- `dh_area_average_median`: Dictionary of median height anomaly data by mission
- `dh_area_average_error`: Dictionary of error data by mission
- `cmap`: Color map for mission lines (default: :thermal)
- `mission_color_width`: Custom colors and line widths for missions
- `mission_order`: Order of missions to plot
- `xtickspacing`: Spacing between x-axis ticks (default: 5)
- `median_in_label`: Include median value in legend labels (default: false)

# Returns
- A Makie Figure object with the height anomaly plot
"""
function plot_area_average_height_anomaly_with_error(dh_area_average_median, dh_area_average_error; cmap=:thermal, mission_color_width=nothing, mission_order=nothing, xtickspacing=5, median_in_label=false)

    mm = 3.7795275590551176
    figure_width = 89mm

    f = Figure(size=(figure_width*1.5, figure_width));
    ax = CairoMakie.Axis(f[1, 1])
    plot_area_average_height_anomaly_with_error!(ax, dh_area_average_median, dh_area_average_error; cmap, mission_color_width, mission_order, xtickspacing, median_in_label)
    
    return f
end


"""
    plot_area_average_height_anomaly_with_error!(ax, dh_area_average_median, dh_area_average_error; cmap=:thermal, mission_color_width=nothing, mission_order=nothing, xtickspacing=5, median_in_label=false)

Plot area-averaged height anomalies with error bands on an existing axis.

# Arguments
- `ax`: Makie axis to plot on
- `dh_area_average_median`: Dictionary of median height anomaly data by mission
- `dh_area_average_error`: Dictionary of error data by mission

# Keyword Arguments
- `cmap`: Color map for mission lines (default: :thermal)
- `mission_color_width`: Custom colors and line widths for missions
- `mission_order`: Order of missions to plot
- `xtickspacing`: Spacing between x-axis ticks (default: 5)
- `median_in_label`: Include median value in legend labels (default: false)

# Returns
- The modified axis object

# Description
Unlike plot_area_average_height_anomaly!(), area_average is precomputed before passing to this function.
Plots error bands first, then lines on top with customizable styling and legend.
"""
function plot_area_average_height_anomaly_with_error!(ax, dh_area_average_median, dh_area_average_error; cmap=:thermal, mission_color_width=nothing, mission_order=nothing, xtickspacing=5, median_in_label=false)
    # mission_color_width = Dict("synthesis" => (color = :black, linewidth = 1))
    ax.ytickformat = values -> ["$(round(Int,value))m" for value in values]
    #ax.title = "median area-averaged height anomaly ± 2σ"

    # ----------- area average height anomaly  ------------    
    clrs = Makie.resample_cmap(cmap, length(keys(dh_area_average_median)) + 1)

    if isnothing(mission_order)
        mission_order = collect(keys(dh_area_average_median))
    end

    date_min = Inf
    date_max = -Inf

    for (i, mission) in enumerate(mission_order)
        ddate = dims(dh_area_average_median[mission], :date)
        x = collect(val(ddate))
        date_min = min(date_min, minimum(x))
        date_max = max(date_max, maximum(x))

        mid0 = dh_area_average_median[mission]
        valid_index = .!isnan.(mid0)

        if !any(valid_index)
            continue
        end

        mid0 = mid0
        low = mid0 .- dh_area_average_error[mission]
        high = mid0 .+ dh_area_average_error[mission]

        CairoMakie.band!(ax, x, collect(low), collect(high); color=(clrs[i], 0.4))
    end

    for (i, mission) in enumerate(mission_order)
        ddate = dims(dh_area_average_median[mission], :date)
        x = collect(val(ddate))
        y = collect(dh_area_average_median[mission])
        y_median = round(median(y[.!isnan.(y)]), digits=1)

        if !in(mission, keys(dh_area_average_median))
            continue
        end

        if median_in_label
            label="$(mission_proper_name(mission)): $(y_median)m"
        else
            label="$(mission_proper_name(mission))"
        end

        if isnothing(mission_color_width) || (!in(mission, keys(mission_color_width)))
            lines!(ax, x, y; label, color=clrs[i])
        else
            lines!(ax, x, y; label, color=mission_color_width[mission].color, linewidth=mission_color_width[mission].linewidth)
        end
    end

    xlims = (floor(Int, date_min / xtickspacing) * xtickspacing, ceil(Int, date_max / xtickspacing) * xtickspacing)
    ax.xticks = xlims[1]:xtickspacing:xlims[2]
    ax.ylabel = "height anomaly"
    xlims!(ax, xlims)

    axislegend(ax, position=:rt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false

    #colgap!(f.layout, 10)
    #rowgap!(f.layout, 10)
    #hidespines!(ax, :t, :r)
    return ax
end

"""
    plot_multiregion_dvdm(regions; kwargs...) -> Figure, Vector, Vector

Create a multi-region plot of glacier mass change time series with stacked regions.

This function visualizes mass change time series for multiple glacier regions, stacking them
vertically with appropriate offsets to avoid overlap. It supports plotting multiple variables
with error bounds and customizable styling.

# Arguments
- `regions`: Dictionary containing mass change data for different regions and variables

# Keyword Arguments
- `variables=["dm", "dm_altim"]`: Variables to plot (last variable is plotted on top)
- `units="Gt"`: Units for mass change (e.g., "Gt", "m", "m w.e.")
- `rgi_regions=setdiff(collect(dims(runs_rgi["dm_altim"], :rgi)), [13, 14, 15, 99])`: RGI regions to include
- `showlines=false`: Whether to show grid lines
- `fontsize=15`: Font size for plot text
- `cmap=:Dark2_4`: Color map for regions
- `region_order=nothing`: Custom ordering of regions (default: sorted by total mass change)
- `ylims=nothing`: Custom y-axis limits
- `title=nothing`: Plot title
- `palette=nothing`: Custom color palette
- `delta_offset=nothing`: Vertical offset between regions (default depends on units)
- `all_error_bounds=false`: Whether to show error bounds for all variables
- `daterange=DateTime(2000,1,1):Month(1):DateTime(2025,1,1)`: Date range to plot
- `numbersinylabel=false`: Whether to include numeric values in y-axis labels

# Returns
- `f`: The Makie Figure object
- `region_order`: Vector of region IDs in the order they were plotted
- `ylims`: The y-axis limits used in the plot
"""
function plot_multiregion_dvdm(
    regions=regions;
    variables=["dm", "dm_altim"], # last variable is plotted last
    units="Gt",
    rgi_regions=setdiff(collect(dims(runs_rgi["dm_altim"], :rgi)), [13, 14, 15, 99]),
    showlines=false,
    fontsize=15,
    cmap=:Dark2_4,
    region_order=nothing,
    ylims=nothing,
    title=nothing,
    palette=nothing,
    delta_offset=nothing,
    all_error_bounds=false,
    daterange=DateTime(2000, 1, 1):Month(1):DateTime(2025, 1, 1),
    numbersinylabel=false,)

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

    yticklabels = rginum2label.(region_order)

    if numbersinylabel
        yticklabels = ["$(rginum2label(id)) $(rginum2enclosed_alphanumerics[id])" for id in region_order]
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

            CairoMakie.band!(ax1, collect(decimalyear.(dims(mid0, :date))), collect(low), collect(high); color=(:grey, 0.4))

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
                ln = CairoMakie.lines!(ax1, collect(decimalyear.(dims(mid0, :date))), collect(mid0))

                y2ticklabels[i] = string(round(Int64, (mid0[end] - mid0[1]))) * " $(units)"
                yticks[i] = mid0[1]
                y2ticks[i] = mid0[end]
                y2ticklabelcolor[i] = ln.color.val
            else
                ln = CairoMakie.lines!(ax1, collect(decimalyear.(dims(mid0, :date))), collect(mid0); color=(:black, 0.3))
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


"""
    show_error_bar_table(regional_results; cols2display = 3)

Display a formatted table of regional results in the console, organized in sections.

Prints the regional results DataFrame in multiple sections, with each section showing a subset
of columns. The first two columns (typically region identifiers) are always displayed, followed
by `cols2display` data columns in each section. Section headers are printed in yellow and
contain the name of the first column in that section.

# Arguments
- `regional_results`: DataFrame containing regional data to be displayed
- `cols2display`: Number of data columns to show in each section (default: 3)
"""
function show_error_bar_table(regional_results; cols2display=3)
    # dump table to console
    colnames = names(regional_results)
    last_set = ceil(Int, (length(colnames) - 2) / (cols2display + 1) - 1)
    for i = 0:last_set
        printstyled("\n------------------------------------------------------- $(colnames[3 + (4 * i)])-------------------------------------------------------\n", color=:yellow)
        if (i == last_set) && (length(colnames) > (last_set + 4))
            show(regional_results[:, [1:2; (3 .+ (4*i)):end]], allrows=true, allcols=true)
        else
            show(regional_results[:, [1:2; (3:(2+cols2display)) .+ (4 * i)]], allrows=true, allcols=true)
        end
    end
end