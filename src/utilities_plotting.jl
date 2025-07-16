"""
    _publication_figure(; columns=1, rows=1)

Create a Figure with publication-standard dimensions for scientific plots.

# Arguments
- `columns`: Number of columns in the figure layout (1 or 2, default: 1)
- `rows`: Number of rows in the figure layout (default: 1)

# Returns
- A CairoMakie Figure object with standardized dimensions

# Description
Creates a Figure with dimensions optimized for publication layouts. The figure width
is determined by the number of columns: 89mm for single-column and 183mm for 
two-column layouts. The height is calculated proportionally based on the aspect ratio
of rows to columns. Uses a scaling factor of 1.5 for the width to accommodate
proper rendering and margins.

# Throws
- `ErrorException`: If columns is not 1 or 2
"""
function _publication_figure(; columns=1, rows=1)
    mm = 3.7795275590551176
    if columns == 1
        figure_width = 89mm
    elseif columns == 2
        figure_width = 183mm
    else
        error("columns must be 1 or 2: input was $columns")
    end
    
    f = CairoMakie.Figure(size=(figure_width * 1.5, figure_width * rows / columns))

    return f
end

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

    f = _publication_figure(columns=2, rows=2)

    sample_rate = ceil(Int, length(dh) / 1E6)

    df = binstats(DataFrame(Ti=decimalyear.(datetime), dh=dh), :Ti, decimalyear.(minimum(datetime):Month(1):maximum(datetime)), :dh, col_function=[median])

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
    f = _publication_figure(columns=2, rows=2)
    
    valid = .!ismissing.(dh_obs)
    valid_range, = validrange(valid)
   
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
   
    f = _publication_figure(columns=2, rows=2)

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
    plot_area_average_height_anomaly!(ax, dh0, area_km2; colormap=:thermal, mission_color_width=nothing, mission_order=nothing, xtickspacing=5)

Plot area-averaged height anomalies for multiple missions on a given axis.

# Arguments
- `ax`: Makie axis to plot on
- `dh0`: Dictionary of dimensional arrays containing height anomaly data for each mission
- `area_km2`: Area data for each geotile [km²]
- `colormap`: Colormap for mission lines (default: :thermal)
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
function plot_area_average_height_anomaly!(ax, dh0, area_km2; colormap=:thermal, mission_color_width = nothing, mission_order = nothing, xtickspacing = 5, median_in_label = false)
    # mission_color_width = Dict("synthesis" => (color = :black, linewidth = 1))
    ax.yaxisposition=:right
    ax.ytickformat=values -> ["$(round(Int,value))m" for value in values]
    ax.title="area-averaged height anomaly"

    # ----------- area average height anomaly  ------------    
    ddate = dims(dh0[first(keys(dh0))], :date)
    x = collect(val(ddate))
    clrs = Makie.resample_cmap(colormap, length(keys(dh0))+1)

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

        if median_in_label
            label = "$(mission_proper_name(mission)): $(y_median)m"
        else
            label = "$(mission_proper_name(mission))"
        end
        if isnothing(mission_color_width) || (!in(mission, keys(mission_color_width)))
            lines!(ax, x, y; label, color=clrs[i])
        else
            lines!(ax, x, y; label, color = mission_color_width[mission].color, linewidth = mission_color_width[mission].linewidth)
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
    plot_elevation_time_multimission_geotiles(dh; geotiles2plot, area_km2, colorrange, label, colorbar_label, mission_order, hypsometry, area_averaged, plots_show, plots_save, plot_save_path_prefix, plot_save_format, colormap, xtickspacing)

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
- `colormap`: Colormap for plots (default: :thermal)
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
    colormap=:thermal,
    xtickspacing = 5,
    median_in_label = true,
    )

    if !isnothing(geotiles2plot)
        figures = Dict()
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
                plot_area_average_height_anomaly!(ax23, dh0, area_km2[geotile=At(geotile2plot)]; mission_order, colormap, xtickspacing, median_in_label)
            end

            isnothing(label) ? nothing : Label(f[0, 1:2], label)

            plots_show ? display(f) : nothing

            if plots_save
                fname = plot_save_path_prefix * "_$(geotile2plot)$(plot_save_format)"
                CairoMakie.save(fname, f)
            end
            figures[geotile2plot] = f
        end
        return figures
    else
        return nothing
    end
end


"""
    load_and_plot_elevation_time_multimission_geotiles(; path2file, geotiles2plot=["lat[+30+32]lon[+078+080]"], colorrange=(-20, 20))

Load hypsometric height anomaly data from a file and plot multi-mission elevation time series for specified geotiles.

# Arguments
- `path2file`: Path to the file containing binned and filled elevation data.
- `geotiles2plot`: Array of geotile identifiers to plot (default: one example geotile).
- `colorrange`: Tuple specifying the color range for the plot (default: (-20, 20)).

# Returns
- `f`: Dictionary of CairoMakie Figure objects, keyed by geotile.
- `dh`: Loaded hypsometric height anomaly data.
- `param`: Parameters parsed from the file name.

# Description
This function loads hypsometric height anomaly data from the specified file, extracts plotting parameters,
and generates multi-mission elevation time series plots for each geotile in `geotiles2plot`. The plots include
hypsometry and area-averaged panels, and are configured for publication-quality output.
"""
function load_and_plot_elevation_time_multimission_geotiles(;
    path2file,
    geotiles2plot = ["lat[+30+32]lon[+078+080]"],
    colorrange = (-20, 20),
)

    param = binned_filled_fileparts(path2file)

    dh = FileIO.load(path2file, "dh_hyps")

    f = plot_elevation_time_multimission_geotiles(
        dh;
        geotiles2plot,
        area_km2 = _geotile_area_km2(; param.surface_mask, param.geotile_width),
        colorrange, 
        label = nothing,
        colorbar_label = "height anomaly",
        mission_order = vcat(plot_order["missions"], "Synthesis"),
        hypsometry = true, 
        area_averaged = true, 
        plots_show = false, 
        plots_save = false, 
        plot_save_path_prefix = "",
        plot_save_format = ".png",
        colormap = :thermal,
        xtickspacing = 5,
        median_in_label = false,
    )

    for geotile2plot in geotiles2plot
        Label(f[geotile2plot][0, :], geotile2plot, fontsize=20)
    end

    return f, dh, param
end



"""
    plot_area_average_height_anomaly_with_ensemble_spread(
        dh_area_averaged; 
        path2runs_synthesized_reference=nothing, 
        geotiles2plot=nothing, 
        ref_period=(Date(2000, 1, 1), Date(2001, 12, 31)), 
        xtickspacing=5, 
        p=0.95
    )

Plot area-averaged height anomaly time series for one or more geotiles, showing the ensemble spread.

For each geotile in `geotiles2plot`, this function:
- Extracts the area-averaged height anomaly ensemble,
- Converts the time axis to decimal years,
- Centers the anomalies relative to the mean over a reference period,
- Plots the ensemble spread and reference run using `plot_ensemble_spread!`.

# Arguments
- `dh_area_averaged`: DimArray or similar, containing area-averaged height anomaly ensemble data, with dimensions (run, date, geotile).
- `path2runs_synthesized_reference`: The run identifier (or index) to use as the reference ensemble member.
- `geotiles2plot`: Array of geotile identifiers to plot. If `nothing`, returns `nothing`.
- `ref_period`: Tuple of two `Date` objects specifying the reference period for centering (default: (2000-01-01, 2001-12-31)).
- `xtickspacing`: Spacing between x-axis ticks (in decimal years, default: 5).
- `p`: Quantile for ensemble spread (default: 0.95).

# Returns
- `figures`: Dictionary mapping geotile identifier to CairoMakie Figure object, or `nothing` if `geotiles2plot` is `nothing`.
"""
function plot_area_average_height_anomaly_with_ensemble_spread(
    dh_area_averaged; 
    path2runs_synthesized_reference = nothing, 
    geotiles2plot = nothing, 
    ref_period = (Date(2000, 1, 1), Date(2001, 12, 31)), 
    xtickspacing = 5, 
    p = 0.95
)
    if !isnothing(geotiles2plot)
        figures = Dict()
        for geotile2plot in geotiles2plot
            dh_area_averaged0 = dh_area_averaged[geotile=At(geotile2plot)]
            valid = .!isnan.(dh_area_averaged0)
            _, date_range, = validrange(valid)

            # convert to decimal year
            drun, ddate = dims(dh_area_averaged0)
            decyear = decimalyear.(collect(ddate))[date_range]
            ddate = Dim{:date}(decyear)

            # recreate DimArray
            dh_area_averaged0 = DimArray(dh_area_averaged0[:, date_range], (drun, ddate))

            # reference period
            ref_period = decimalyear.(ref_period)

            # center data around reference period
            ref_mean = mean(dh_area_averaged0[date=ref_period[1] .. ref_period[2]], dims=:date)
            @d dh_area_averaged0 .-= ref_mean

            figures[geotile2plot] = f = _publication_figure(columns=1, rows=1)
            ax = CairoMakie.Axis(figures[geotile2plot][1, 1])

            plot_ensemble_spread!(ax, dh_area_averaged0; path2runs_synthesized_reference, xtickspacing, p)
        end
        return figures
    else
        return nothing
    end
end

"""
    plot_ensemble_spread!(ax, dh_area_averaged0; path2runs_synthesized_reference=nothing, xtickspacing=5, p=0.95)

Plot the ensemble spread and reference time series of area-averaged height anomaly on a given axis.

# Arguments
- `ax`: A CairoMakie Axis object to plot on.
- `dh_area_averaged0`: A DimArray containing area-averaged height anomaly for each ensemble member and date.
- `path2runs_synthesized_reference`: The run identifier (or index) to use as the reference ensemble member.
- `xtickspacing`: Spacing between x-axis ticks (in decimal years, default: 5).
- `p`: Quantile for ensemble spread (default: 0.95).

# Description
Plots the following on the provided axis:
- The reference ensemble member time series (in black).
- The ensemble spread (quantile `p` of the absolute difference from the reference) as a shaded band.
- All ensemble members as faint lines.
- A legend indicating the reference, ensemble member, and ensemble uncertainty.

# Returns
- The modified axis object with the plot.
"""
function plot_ensemble_spread!(ax, dh_area_averaged0; path2runs_synthesized_reference=nothing, xtickspacing=5, p=0.95)

    # ensemble spread
    dh_area_averaged_spread = @d dh_area_averaged0 .- dh_area_averaged0[run=At(path2runs_synthesized_reference)]
    fun = x -> quantile(x, p)
    dh_area_averaged_spread = dropdims(mapslices(fun, abs.(parent(dh_area_averaged_spread)), dims=1), dims=1)

    ax.ytickformat = values -> ["$(round(Int,value))m" for value in values]
    drun, ddate = dims(dh_area_averaged0)
    x = collect(val(ddate))

    mid0 = collect(dh_area_averaged0[run=At(path2runs_synthesized_reference)])
    low = mid0 .- dh_area_averaged_spread
    high = mid0 .+ dh_area_averaged_spread

    CairoMakie.band!(ax, x, low, high; color=(:black, 0.3), label="2σ ensemble uncertainty")

    for run in drun
        lines!(ax, x, vec(dh_area_averaged0[run=At(run)]); color=(:mediumpurple2, 0.5), linewidth=0.1)
    end
    # this is just to get the legend entry for the ensemble member
    lines!(ax, x, vec(dh_area_averaged0[run=At(drun[1])]); color=:mediumpurple2, linewidth=1, label="ensemble member", visible=false)

    #plot reference
    lines!(ax, x, vec(dh_area_averaged0[run=At(path2runs_synthesized_reference)]); color=:black, linewidth=1, label="reference")

    date_min = minimum(collect(ddate))
    date_max = maximum(collect(ddate))
    xlims = (floor(Int, date_min / xtickspacing) * xtickspacing, ceil(Int, date_max / xtickspacing) * xtickspacing)
    ax.xticks = xlims[1]:xtickspacing:xlims[2]
    ax.ylabel = "height anomaly"
    xlims!(ax, xlims)

    axislegend(ax, position=:rt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1, visible=true) # orientation=:horizontal, framevisible=false

    return ax
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
    plot_area_average_height_anomaly_with_error(dh_area_average_median, dh_area_average_error; colormap=:thermal, mission_color_width=nothing, mission_order=nothing, xtickspacing=5, median_in_label=false)

Create a plot showing area-averaged height anomalies with error bands.

# Arguments
- `dh_area_average_median`: Dictionary of median height anomaly data by mission
- `dh_area_average_error`: Dictionary of error data by mission
- `colormap`: Color map for mission lines (default: :thermal)
- `mission_color_width`: Custom colors and line widths for missions
- `mission_order`: Order of missions to plot
- `xtickspacing`: Spacing between x-axis ticks (default: 5)
- `median_in_label`: Include median value in legend labels (default: false)

# Returns
- A Makie Figure object with the height anomaly plot
"""
function plot_area_average_height_anomaly_with_error(dh_area_average_median, dh_area_average_error; colormap=:thermal, mission_color_width=nothing, mission_order=nothing, xtickspacing=5, median_in_label=false)

    f = _publication_figure(columns=1, rows=1)
    ax = CairoMakie.Axis(f[1, 1])
    plot_area_average_height_anomaly_with_error!(ax, dh_area_average_median, dh_area_average_error; colormap, mission_color_width, mission_order, xtickspacing, median_in_label)
    
    return f
end


"""
    plot_area_average_height_anomaly_with_error!(ax, dh_area_average_median, dh_area_average_error; colormap=:thermal, mission_color_width=nothing, mission_order=nothing, xtickspacing=5, median_in_label=false)

Plot area-averaged height anomalies with error bands on an existing axis.

# Arguments
- `ax`: Makie axis to plot on
- `dh_area_average_median`: Dictionary of median height anomaly data by mission
- `dh_area_average_error`: Dictionary of error data by mission

# Keyword Arguments
- `colormap`: Color map for mission lines (default: :thermal)
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
function plot_area_average_height_anomaly_with_error!(ax, dh_area_average_median, dh_area_average_error; colormap=:thermal, mission_color_width=nothing, mission_order=nothing, xtickspacing=5, median_in_label=false)
    # mission_color_width = Dict("synthesis" => (color = :black, linewidth = 1))
    ax.ytickformat = values -> ["$(round(Int,value))m" for value in values]
    #ax.title = "median area-averaged height anomaly ± 2σ"

    # ----------- area average height anomaly  ------------    
    clrs = Makie.resample_cmap(colormap, length(keys(dh_area_average_median)) + 1)

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
    plot_model_fit(best_smb, best_fac, best_discharge, best_fit, dv0, cost_metric; xlims, xtickspacing, colormap=:thermal)

Create a two-panel figure visualizing the model fit and cost metric for glacier mass balance modeling.

# Arguments
- `best_smb`: Array of best-fit Surface Mass Balance (SMB) values over time.
- `best_fac`: Array of best-fit Firn Air Content (FAC) values over time.
- `best_discharge`: Array of best-fit discharge values over time.
- `best_fit`: Array of best-fit modeled height anomaly values over time.
- `dv0`: Array of observed height anomaly values over time.
- `cost_metric`: 2D array of cost metric values (e.g., RMSE or MAD) as a function of model parameters.
- `xlims`: Tuple specifying the x-axis limits for the time series plot.
- `xtickspacing`: Spacing between x-axis ticks for the time series plot.
- `colormap`: (Optional) Colormap to use for plotting (default: `:thermal`).

# Returns
- `f`: A CairoMakie Figure object with two subplots:
    - Left: Time series of SMB, FAC, discharge, observed, and modeled height anomalies.
    - Right: Contour plot of the cost metric as a function of model parameters.

"""
function plot_model_fit(best_smb, best_fac, best_discharge, best_fit, dv0, cost_metric; xlims, xtickspacing, colormap=:thermal)
    f = _publication_figure(columns=2, rows=1)
    ax1 = CairoMakie.Axis(f[1, 1])
    ax2 = CairoMakie.Axis(f[1, 2])

   

    plot_best_fit!(ax1, best_smb, best_fac, best_discharge, best_fit, dv0; xlims, xtickspacing, colormap)
    ax2, crange = plot_cost_metric!(ax2, cost_metric; colormap)
    Colorbar(f[1, 3]; limits=crange, colormap)#label="cost"

    return f
end

"""
    plot_best_fit!(ax, best_smb, best_fac, best_discharge, best_fit, dv0; xlims, xtickspacing, colormap=:thermal)

Plot the best-fit glacier mass balance model components and observations on a given axis.

# Arguments
- `ax`: A CairoMakie Axis object to plot on.
- `best_smb`: Array of best-fit Surface Mass Balance (SMB) values over time.
- `best_fac`: Array of best-fit Firn Air Content (FAC) values over time.
- `best_discharge`: Array of best-fit discharge values over time.
- `best_fit`: Array of best-fit modeled height anomaly values over time.
- `dv0`: Array of observed height anomaly values over time.

# Keyword Arguments
- `xlims`: Tuple specifying the x-axis limits for the time series plot.
- `xtickspacing`: Spacing between x-axis ticks for the time series plot.
- `colormap`: (Optional) Colormap to use for plotting (default: `:thermal`).

# Returns
- The modified axis object with the plotted lines and legend.

This function normalizes all time series to a common reference period (the first year), then plots
the SMB, FAC, discharge, observed, and modeled height anomalies as lines on the provided axis.
A legend is added for clarity.
"""
function plot_best_fit!(ax, best_smb, best_fac, best_discharge, best_fit, dv0; xlims, xtickspacing, colormap=:thermal)

    ddate = dims(best_smb, :date)
    decyear = decimalyear.(collect(val(ddate)))
    clrs = Makie.resample_cmap(colormap, 6)

    best_smb = best_smb .- mean(best_smb[ddate[1]..(ddate[1]+Year(1))])
    best_fac = best_fac .- mean(best_fac[ddate[1]..(ddate[1]+Year(1))])
    best_discharge = best_discharge .- mean(best_discharge[ddate[1]..(ddate[1]+Year(1))])
    best_fit = best_fit .- mean(best_fit[ddate[1]..(ddate[1]+Year(1))])
    dv0 = dv0 .- mean(dv0[ddate[1]..(ddate[1]+Year(1))])

    ax.ylabel="height anomaly"
    ax.ytickformat=values -> ["$(round(Int,value))m" for value in values]
    ax.xticks = xlims[1]:xtickspacing:xlims[2]
    xlims!(ax, xlims)
    #title="best fit for $single_geotile_test [pscale = $pscale0, Δheight = $Δheight0, mad = $(round(cost_metric_minimum; digits=2))]"

    CairoMakie.lines!(ax, decyear, parent(dv0); label="observed", color=clrs[1])
    CairoMakie.lines!(ax, decyear, parent(best_fac); label="FAC", color=clrs[2])
    CairoMakie.lines!(ax, decyear, parent(best_discharge); label="discharge", color=clrs[3])
    CairoMakie.lines!(ax,decyear, parent(best_smb); label="SMB", color=clrs[4])
    CairoMakie.lines!(ax,decyear, parent(best_fit); label="modeled", color=clrs[5])

    axislegend(ax, position=:lb, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false
    return ax
end

"""
    plot_cost_metric!(ax, cost_metric; colormap=:thermal)

Plot a contour map of the cost metric as a function of precipitation scale (pscale) and Δheight.

# Arguments
- `ax`: A CairoMakie Axis object to plot on.
- `cost_metric`: 2D array (or NamedDimsArray) of cost metric values, with dimensions :pscale and :Δheight.

# Keyword Arguments
- `colormap`: Colormap to use for the contour plot (default: `:thermal`).

# Returns
- `(ax, crange)`: The modified axis object and the color range tuple used for the plot.

This function creates a contour plot of the cost metric over the parameter space of precipitation scale and Δheight.
It highlights the minimum cost location with a marker and label, and adds a colorbar for reference.
"""
function plot_cost_metric!(ax, cost_metric; colormap=:thermal)
    
    ax.xlabel = "precipitation scaling"
    ax.ylabel = "height offset / melt scaling"
    ax.ytickformat = values -> ["$(round(Int,value))m" for value in values]
    dpscale = dims(cost_metric, :pscale)
    dΔheight = dims(cost_metric, :Δheight)


    step = ceil(Int, minimum(cost_metric)/4)
    crange = (0, step*6)
    contour!(ax, collect(val(dpscale)), collect(val(dΔheight)), parent(cost_metric); colorrange=crange, levels=0:step:maximum(cost_metric), colormap)
    #heatmap!(ax, cost_metric; colorrange=(0, 50))

    xlims!(ax, extrema(collect(val(dpscale))))
    ylims!(ax, extrema(collect(val(dΔheight))))

    index_minimum = argmin(cost_metric)
    pscale_best = DimPoints(cost_metric)[index_minimum][1]
    Δheight_best = DimPoints(cost_metric)[index_minimum][2]

    scatter!(ax, pscale_best, Δheight_best; color=:black, markersize=15, marker=:xcross)
    text!(ax, pscale_best, Δheight_best,
        text="  cost = $(round(cost_metric[index_minimum]; digits=1))",
        align=(:left, :center)
    )
    return ax, crange
end

"""
    plot_hist_gemb_altim_trend_amplitude(geotiles0)

Plot histograms comparing model and altimetry-derived trends and amplitudes.

# Arguments
- `geotiles0`: A DataFrame or table containing columns for model and altimetry trends and amplitudes.
    Expected columns: `"dv_trend"`, `"dv_altim_trend"`, `"dv_amplitude"`, `"dv_altim_amplitude"`.

# Returns
- `f`: A CairoMakie Figure object with two panels:
    - Left: Histogram of trends (model vs. altimetry)
    - Right: Histogram of amplitudes (model vs. altimetry)

# Description
Creates a two-panel figure. The first panel shows step histograms of the model and altimetry-derived trends.
The second panel shows step histograms of the model and altimetry-derived amplitudes. Legends are included for clarity.
"""
function plot_hist_gemb_altim_trend_amplitude(geotiles0)
    f = _publication_figure(columns=2, rows=1)
    var1 = "dh_trend"
    var2 = "dh_altim_trend"

    ax1 = CairoMakie.Axis(f[1, 1]; xlabel="trend", ylabel="count")
    ax1.xtickformat = values -> ["$(round(Int,value))m yr⁻¹" for value in values]
    CairoMakie.stephist!(geotiles0[:, var1], bins=-5:0.25:5; label="model")
    CairoMakie.stephist!(geotiles0[:, var2], bins=-5:0.25:5; label="synthesis")
    axislegend(ax1, framevisible=false)

    var1 = "dh_amplitude"
    var2 = "dh_altim_amplitude"
    ax2 = CairoMakie.Axis(f[1, 2]; xlabel="amplitude")
    ax2.xtickformat = values -> ["$(round(Int,value))m" for value in values]

    CairoMakie.stephist!(geotiles0[:, var1], bins=0:0.25:5; label="model")
    CairoMakie.stephist!(geotiles0[:, var2], bins=0:0.25:5; label="synthesis")
    axislegend(ax2, framevisible=false)

    return f
end



"""
    plot_hist_pscale_Δheight(gemb_fit)

Plot histograms of precipitation scaling and height offset parameters from GEMB fit results.

# Arguments
- `gemb_fit`: A DataFrame or table containing columns `:pscale` (precipitation scaling) and `:Δheight` (height offset).

# Returns
- `f`: A CairoMakie Figure object with two panels:
    - Left: Histogram of precipitation scaling values.
    - Right: Histogram of height offset values.

# Description
Creates a two-panel figure. The first panel shows a histogram of precipitation scaling factors (`pscale`).
The second panel shows a histogram of height offset values (`Δheight`). Both histograms have counts on the y-axis,
and the y-axis is set to start at zero.
"""
function plot_hist_pscale_Δheight(gemb_fit)
    f = _publication_figure(columns=2, rows=1)
    ax1 = CairoMakie.Axis(f[1, 1]; xlabel="precipitation scaling", ylabel="count")
    CairoMakie.hist!(ax1, gemb_fit[:, :pscale], bins=0.25:0.25:4)
    ylims!(low=0)
    ax2 = CairoMakie.Axis(f[1, 2]; xlabel="height offset / melt scaling")
    ax2.xtickformat = values -> ["$(round(Int,value))m" for value in values]
    CairoMakie.hist!(ax2, gemb_fit[:, :Δheight], bins=-3000:200:3000)
    ylims!(low=0)
    #axislegend(ax1, framevisible = false)
    return f
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
- `colormap=:Dark2_4`: Color map for regions
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
    colormap=:Dark2_4,
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
        palette = (; color=Makie.resample_cmap(colormap, n))
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

                if !occursin("yr⁻¹", units)
                    y2ticklabels[i] = string(round(Int64, (mid0[end] - mid0[1]))) * " $(units)"
                else
                    fit = ts_seasonal_model(mid0)
                    y2ticklabels[i] = string(round(Int64, fit.trend)) * " $(units)"
                end

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


"""
    plot_area_average_height_gemb_ensemble(gemb0, area_km2; vars2plot=["melt", "refreeze", "acc", "fac", "rain", "ec"], geotiles2plot=[geotiles_golden_test[1]], title_prefix="")

Plot area-averaged height change time series for GEMB ensemble variables and geotiles.

# Arguments
- `gemb0`: Dictionary or NamedTuple of GEMB output variables, each as a DimArray.
- `area_km2`: DimArray or array of glacier area (km²) for each geotile.
- `vars2plot`: Array of variable names to plot (default: ["melt", "refreeze", "acc", "fac", "rain", "ec"]).
- `geotiles2plot`: Array of geotile identifiers to plot (default: [geotiles_golden_test[1]]).
- `title_prefix`: String prefix for plot titles.

# Returns
- `figures`: Nested dictionary mapping geotile and variable names to Makie figures.
"""
function plot_area_average_height_gemb_ensemble(
    gemb0, area_km2;
    vars2plot = ["melt", "refreeze", "acc", "fac", "rain", "ec"],
    geotiles2plot = [geotiles_golden_test[1]],
    title_prefix = ""
)
    dpscale = dims(gemb0[first(vars2plot)], :pscale)
    if hasdim(gemb0[first(vars2plot)], :Δheight)
        dΔheight = dims(gemb0[first(vars2plot)], :Δheight)
        clrs = Makie.resample_cmap(:thermal, length(dpscale) * length(dΔheight) + 1)
    else
        dΔheight = nothing
        clrs = Makie.resample_cmap(:thermal, length(dpscale) + 1)
    end

    figures = Dict()

    for geotile2plot in geotiles2plot
        figures[geotile2plot] = Dict()
        for var0 in vars2plot
            figures[geotile2plot][var0] = _publication_figure(columns=1, rows=1)
            title = "$(title_prefix) $(var0)-$(geotile2plot)"
            ax = CairoMakie.Axis(figures[geotile2plot][var0][1, 1]; title)
            ax.ytickformat = values -> ["$(round(Int, value))m" for value in values]

            cnt = 1
            for pscale in dpscale
                if isnothing(dΔheight)
                    dh = dh_area_average(
                        gemb0[var0][geotile=At(geotile2plot), pscale=At(pscale)],
                        area_km2[geotile=At(geotile2plot)]
                    )
                    height_range, = validrange(.!isnan.(dh))
                    lines!(ax, dh[height_range]; label="p:$(pscale)", color=clrs[cnt])
                    cnt += 1
                else
                    for Δheight in dΔheight
                        dh = dh_area_average(
                            gemb0[var0][pscale=At(pscale), Δheight=At(Δheight)],
                            area_km2[geotile=At(geotile2plot)]
                        )
                        height_range, = validrange(.!isnan.(dh))
                        lines!(ax, dh[height_range]; label="p:$(pscale) Δh:$(Δheight)", color=clrs[cnt])
                        cnt += 1
                    end
                end
            end
            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1)
        end
    end
    return figures
end


"""
    plot_dh_gemb_ensemble(gemb_dv, area_km2; geotile = "", title_prefix = "")

Plot area-normalized height change ensemble results from GEMB model diagnostics.

# Arguments
- `gemb_dv`: Dictionary or NamedTuple of GEMB diagnostic variables, each as a DimArray with `:pscale` and `:Δheight` dimensions.
- `area_km2`: DimArray or array representing the area (in km²) for normalization.
- `geotile`: (Optional) String identifier for the geotile being plotted (used in plot titles).
- `title_prefix`: (Optional) String prefix for plot titles.

# Returns
- `figures`: Dictionary mapping each variable name to its corresponding Makie figure.

Each figure shows the area-normalized height change for all parameter scale (`pscale`) and height change (`Δheight`) ensemble members, with a legend indicating the parameter combinations.
"""
function plot_dh_gemb_ensemble(gemb_dv, area_km2; geotile = "", title_prefix="")
    area_total = sum(parent(area_km2))
    vars2plot = keys(gemb_dv)
    dpscale = dims(gemb_dv[first(vars2plot)], :pscale)
    dΔheight = dims(gemb_dv[first(vars2plot)], :Δheight)
    clrs = Makie.resample_cmap(:thermal, length(dpscale) * length(dΔheight) + 1)
    figures = Dict()

    for var0 in vars2plot
        figures[var0] = _publication_figure(columns=1, rows=1)
        title = "$(title_prefix) $(var0)-$(geotile)"
        ax = CairoMakie.Axis(figures[var0][1, 1]; title)
        ax.ytickformat = values -> ["$(round(Int,value))m" for value in values]

        cnt = 1
        for pscale in dpscale
            for Δheight in dΔheight
                dh = gemb_dv[var0][pscale=At(pscale), Δheight=At(Δheight)] ./ area_total
                height_range, = validrange(.!isnan.(dh))
                lines!(ax, dh[height_range]; label="p:$(pscale) Δh:$(Δheight)", color=clrs[cnt])
                cnt += 1
            end
        end

        axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
    end

    return figures
end



"""
    plot_point_location_river_flux(land_flux, glacier_flux; dates4plot, COMID, name)

Plot river flux components (land and glacier) for a specific location, without snow flux.

This is a convenience wrapper for `plot_point_location_river_flux` that omits the snow flux argument.

# Arguments
- `land_flux`: DimArray of land river flux (time × COMID).
- `glacier_flux`: DimArray of glacier river flux (time × COMID).
- `dates4plot`: Tuple of (start_date, end_date) for plotting.
- `COMID`: Integer river reach identifier.
- `name`: String for plot title.

# Returns
- `f`: Makie Figure object.
"""
function plot_point_location_river_flux(land_flux, glacier_flux; date_range, COMID, name, show_title=true)
    f = plot_point_location_river_flux(land_flux, glacier_flux, nothing; date_range, COMID, name, show_title)
    return f
end

"""
    plot_point_location_river_flux(land_flux, glacier_flux, snow_flux; dates4plot, COMID, name)

Plot river flux components and glacier contribution for a specific location.

Creates a figure with:
- Stacked area plots of land, glacier, and (optionally) snow river fluxes over time for a given COMID.
- Glacier fraction of total river flux, with seasonal maximum (gmax) highlighted.

# Arguments
- `land_flux`: DimArray of land river flux (time × COMID).
- `glacier_flux`: DimArray of glacier river flux (time × COMID).
- `snow_flux`: DimArray of snow river flux (time × COMID), or `nothing` if not available.
- `dates4plot`: Tuple of (start_date, end_date) for plotting.
- `COMID`: Integer river reach identifier.
- `name`: String for plot title.

# Returns
- `f`: Makie Figure object.
"""
function plot_point_location_river_flux(land_flux, glacier_flux, snow_flux; date_range, COMID, name, show_title=true)
   

    f = _publication_figure(columns=1, rows=1)

    date_interval = date_range[1]..date_range[2]

    land_flux0 = land_flux[Ti=date_interval, COMID=At(COMID)]
    glacier_flux0 = glacier_flux[Ti=date_interval, COMID=At(COMID)]

    if !isnothing(snow_flux)
         seperate_out_snow = true
         snow_flux0 = snow_flux[Ti = date_interval, COMID=At(COMID)]
    else
        seperate_out_snow = false
    end

    x = decimalyear.(val(dims(land_flux0, :Ti)))
    xticks = collect(floor(minimum(x)/5)*5:5:ceil(maximum(x)/5)*5)

    seperate_out_snow && (y_snow = parent(snow_flux0) ./ 1000)

    ax1 = CairoMakie.Axis(f[1:2, 1]; xticklabelsvisible=false, xticksvisible=false, xticks=xticks)
    if show_title
        ax1.title = name
    end

    # plot land flux
    y = parent(land_flux0) ./ 1000
    seperate_out_snow && (y .-= y_snow)

    x0 = vcat(x[1], x, x[end])
    yunits = Unitful.unit(land_flux[1])
    y = vcat(0 * yunits, y, 0 * yunits)
    max_y = maximum(y)
    poly!(ax1, GI.Point.(ustrip(x0), ustrip(y)); color=(:peru, 1), label="land")

    # plot snow flux
    if seperate_out_snow
        y = y_snow
        yunits = Unitful.unit(y[1])
        y = vcat(0.0 * yunits, y, 0.0 * yunits)
        max_y = max(maximum(y), max_y)
        poly!(ax1, GI.Point.(x0, ustrip(y)); color=(:gray, 0.5), label="snow")
    end

    y = parent(glacier_flux0) ./ 1000
    yunits = Unitful.unit(y[1])
    y = vcat(0 * yunits, y, 0 * yunits)
    poly!(ax1, GI.Point.(ustrip(x0), ustrip(y)); color=(:blue, 0.5), label="glacier")

    axislegend(ax1; merge=true, backgroundcolor=(:white, 0), framevisible=false, position=:lt)
    ax1.ylabel = "river flux (m³/s) x 10⁻³"
    xlims!(ax1, minimum(xticks), maximum(xticks))
    ylims!(ax1, 0, ustrip(max_y))

    ax2 = CairoMakie.Axis(f[3, 1], xticks=xticks)
    glacier_frac = glacier_flux0 ./ (land_flux0 .+ glacier_flux0) * 100

    gmax = groupby(glacier_frac, Ti => Bins(month, 1:12))
    gmax = mean.(gmax; dims=:Ti)
    gmax = cat(gmax...; dims=dims(gmax, :Ti))
    dTi = dims(gmax, :Ti)
    gmax_month = dTi[argmax(gmax)]
    gmax = round(Int8, maximum(gmax))

    foo = glacier_frac[month.(dims(glacier_frac, :Ti)).==gmax_month]
    gmax_min = round(Int8, minimum(foo))
    gmax_max = round(Int8, maximum(foo))

    lines!(ax2, x, parent(glacier_frac); color=(:blue, 0.5))
    lines!(ax2, [minimum(x), maximum(x)], [gmax_min, gmax_min]; color=(:gray, 0.5), linestyle=:dash)
    lines!(ax2, [minimum(x), maximum(x)], [gmax_max, gmax_max]; color=(:gray, 0.5), linestyle=:dash)
    lines!(ax2, [minimum(x), maximum(x)], [gmax, gmax]; color=(:black, 1))
    ax2.ylabel = "glacier fraction [%]"

    xlims!(ax2, minimum(xticks), maximum(xticks))
    text!(maximum(x), mean((gmax_max, gmax)), text="gmax = $gmax [$(Dates.monthname(gmax_month))] ", align=[:right, :center], color=(:black, 0.8), overdraw=true)

    return f
end