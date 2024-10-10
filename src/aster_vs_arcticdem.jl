begin
    using Altim
    using Arrow
    using DataFrames
    using NearestNeighbors
    using Plots
    using Statistics
    using DimensionalData
    using Dates

    project_id = :v01
    geotile_width = 2;
    paths = project_paths(; project_id)
    mission = "hugonnet"
    binning_method = "meanmadnorm3"

    min_frac = 0.2; 
    mask = :landice

    for rgi = ["rgi1", "rgi3", "rgi5", "rgi7", "rgi9"]
    for mask = [:landice, :land]
    
    mission_geotile_folder = paths[Symbol(mission)].geotile 
    geotiles = Altim.geotiles_w_mask(geotile_width)

    geotile = geotiles[findfirst((geotiles[:, rgi].>0).&(geotiles.glacier_frac.>min_frac)), :]

    path2altim = joinpath(mission_geotile_folder, geotile.id * ".arrow");
    path2masks = joinpath(mission_geotile_folder, geotile.id * ".masks");

    altim = DataFrame(Arrow.Table(path2altim));
    masks0 = select!(DataFrame(Arrow.Table(path2masks)), [mask])

    altim = copy(altim[masks0[:, mask],:])

    altim[!, :isaster] .= getindex.(altim.id, 1) .== 'A';

    altim, epsg = itslive_proj!(altim; height=nothing)

    # find non-aster data over ice and compare to nearby observations
    points = hcat(altim.X, altim.Y)
    #tree = KDTree(points[.!altim.isaster, :]; leafsize = 25)
    #tree = KDTree(points[1:1000, :]; leafsize = 25)
    altim[!, :aster_neighbor] .= 0.

    #for r in eachrow(df_notaster)
    altim[!, :decimalyear] = Altim.decimalyear.(altim.datetime)

    altim = altim[(abs.(altim.height) .< 9999) .& (altim.quality .| (.!altim.isaster)) .& (.!isnan.(altim.height)), :]
    altim[!, :dh] = altim.height .- altim.height_reference

    df_aster = altim[altim.isaster,:]
    df_arcticdem = altim[.!altim.isaster,:]

    # define date bin ranges... this should match what's in gemb_classes_binning.jl
    # define date and hight binning ranges 
    date_range, date_center = Altim.project_date_bins()
    Δd = abs(date_center[2] - date_center[1])
    height_range, height_center = Altim.project_height_bins()
    Δh = abs(height_center[2] - height_center[1])

    binningfun = Altim.binningfun_define(binning_method)

    aster, aster_nobs = Altim.geotile_bin2d(
    df_aster;
    var2bin="dh",
    dims_edges=("decimalyear" => Altim.decimalyear.(date_range), "height_reference" => height_range),
    binfunction=Altim.binningfun_define(binning_method))

    arcticdem, arcticdem_nobs = Altim.geotile_bin2d(
    df_arcticdem;
    var2bin="dh",
    dims_edges=("decimalyear" => Altim.decimalyear.(date_range), "height_reference" => height_range),
    binfunction=Altim.binningfun_define(binning_method))


    title = "wv - aster: $(mask) $(geotile.id)"

    p = Plots.heatmap(arcticdem .- aster; clims=(-5, 5), color=:redblue, xlims=(0, 2000), ylims=(2008, 2020), title)
    display(p)

    p = Plots.histogram((arcticdem.-aster)[:]; title, bins=-10:10, legend = nothing)
    vline!([Altim.nanmedian((arcticdem.-aster)[:])]; linewidth=2)
            display(p)




    for (i, ts) = enumerate([Altim.nanmedian.(eachrow(aster)), Altim.nanmedian.(eachrow(arcticdem)), Altim.nanmedian.(eachrow(arcticdem .- aster))])

        if i == 1
            title = "aster: $(mask) $(geotile.id)"
        elseif i == 2
            title = "wv: $(mask) $(geotile.id)"
        elseif i == 3
            title = "wv - aster: $(mask) $(geotile.id)"
        end

        p = Plots.plot(ts; xlims=(2008, 2021), title);

        dates = Altim.decimalyear2datetime.(dims(ts, :decimalyear))
        months = month.(dates)

        colors = cgrad(:phase, 12, categorical = true)

        for (c, m) in enumerate(unique(months));
            index = months .== m
            Plots.plot!(ts[index];xlims=(2008, 2021), seriestype=:scatter, color = colors[c], label = "$(Dates.monthname(m))", title, ylims = (-5, 5), foreground_color_legend = nothing, background_color_legend = nothing)
        end
        
        display(p)
    end




    end
end
    
end
