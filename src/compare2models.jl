"""
    compare2models.jl

Compare glacier runoff estimates between GEMB and models from Wimberly et al. 2024.

This script:
1. Loads glacier geometry and matches glaciers to major river basins
2. Processes glacier runoff data from GEMB and compares with Wimberly 2024 models
3. Calculates and visualizes differences in trends, means, and monthly patterns
4. Generates comparison plots for runoff trends, means, and seasonal patterns

The analysis focuses on the period defined by `dates2average` and produces
multiple visualization outputs in the figures directory.
"""

begin
    using Altim
    using NCDatasets
    using GeoDataFrames
    import GeometryOps as GO
    using DimensionalData
    using SortTileRecursiveTree
    using CairoMakie
    import GeoInterface as GI
    using FileIO
    using DataFrames
    using Dates
    using Statistics
    using LsqFit
    using Unitful
    using Altim.MyUnits

    dates2average = [Date("2000-01-01"), Date("2025-01-01")]
    #glacier_runoff_reference = "binned_unfiltered/2deg/glacier_dh_best_cc_meanmadnorm5_v01_filled_ac_p2_synthesized_perglacier.jld2"
    glacier_summary_file = joinpath("/mnt/bylot-r3/data/project_data", "gardner2025_glacier_summary.nc")
   
    paths = Altim.pathlocal
    km2Gt = 910/1000
end

"""
Match individual glaciers to major river basins using spatial intersection.

Loads glacier outlines and river basin boundaries, then identifies which glaciers
fall within the river basins analyzed in Wimberly 2024 study.
"""
begin #[70s cold]
    
    # load glacier geometry
    glaciers = GeoDataFrames.read(paths[:glacier_individual_outlines])
    glaciers[!, :center_pts] = GI.Point.(glaciers.CenLon, glaciers.CenLat)

    # load large river basins and subset to those in wimberly 2024
    basins = GeoDataFrames.read(paths[:river_major_basins])

    # load wimberly 2024 data to subset basins [is in units of Gt]
    runoff_modeled = Altim.wimberly2024()
    dbasin = collect(dims(runoff_modeled, :Basin))

    # this is done in 2 steps to minimize risk of `occursin` returning an incorrect match
    basin_indices = [findfirst(basins.RIVERBASIN .== basin) for basin in dbasin]
    missing_index = isnothing.(basin_indices)
    basin_indices[missing_index] = [findfirst(occursin.(Ref(basin), basins.RIVERBASIN)) for basin in dbasin[missing_index]]
    basins_subset_df = basins[basin_indices, :]

    # identify those glaciers that fall within select river basins
    center_pts_tree = STRtree(glaciers.center_pts)
    glaciers[!, :basin] .= "";

    Threads.@threads for i in eachindex(basins_subset_df.geometry)
        query_index = SortTileRecursiveTree.query(center_pts_tree, basins_subset_df.geometry[i])
        index = query_index[GO.intersects.(Ref(basins_subset_df.geometry[i]), glaciers.center_pts[query_index])]
        glaciers[index, :basin] .= dbasin[i]
    end
end

"""
Process glacier runoff data and calculate basin-level statistics.

Loads glacier runoff data, groups by basin, and calculates monthly patterns,
trends, and means for comparison with Wimberly 2024 models.
"""
begin
    # load glacier runoff data
    glacier_runoff = NCDataset(glacier_summary_file)

    # convert to DimensionalData
    glacier_runoff = Altim.nc2dd(glacier_runoff["runoff"])

    
    dTi = dims(glacier_runoff, :Ti)

    # convert form DateTime to Date to be compatible with wimberly 2024
    ddate = Dim{:date}(Date.(dTi))
    drgiid = dims(glacier_runoff, :rgiid)

    sindex = sortperm(val(drgiid))
    glacier_runoff = glacier_runoff[rgiid=sindex]

    drgiid = dims(glacier_runoff, :rgiid)
    
    # ensure that glaciers0 RGIId order matches glaciers order
    @assert collect(val(drgiid)) == glaciers.RGIId

    # add glacier runoff
    glaciers[!, :runoff] = eachrow(parent(glacier_runoff))
    glaciers0 = deepcopy(glaciers)

    # drop glaciers not in wimberly 2024 basins
    glaciers0 = glaciers0[glaciers0.basin .!= "", :]

    # group glaciers0 by basin
    gdf = groupby(glaciers0, :basin)

    # sum runoff by basin
    dbasin = dims(runoff_modeled, :Basin)

    index_date = (ddate .> dates2average[1]) .& (ddate .< dates2average[2])

    basin_runoff = zeros(Ti(collect(ddate[findall(index_date)]) .- Day(15)), dbasin)* u"Gt"

    for basin in dbasin
        foo = reduce(hcat, gdf[(basin,)].runoff);
        foo[isnan.(foo)] .= 0 * u"Gt"
        cummulative_runoff = vec(sum(foo, dims=2))
        runoff = diff(cummulative_runoff)

        # `diff` shifts time to the left ... instead we want to shift time to the right
        basin_runoff[:, At(basin)] = runoff[index_date[2:end]]
    end

    # compare to wimberly 2024

    # trim data to overlaping period
    runoff_modeled = runoff_modeled[dates2average[1]..dates2average[2], :, :, :, :]

    #for basin in dbasin
    dssp = dims(runoff_modeled, :SSP)
    dmodel = dims(runoff_modeled, :Model)

    dmodel = Dim{:model}(vcat("GEMB", collect(dmodel)))
    dmonth = Dim{:month}(1:12)
    basin_runoff_monthly = zeros(dmonth, dbasin, dssp, dmodel)
    basin_runoff_trend = zeros(dbasin, dssp, dmodel)
    basin_runoff_mean = zeros(dbasin, dssp, dmodel)

    for basin in dbasin
        for ssp in dssp
            for model in dmodel
                if model == "GEMB"
                    ts = basin_runoff[:, At(basin)]
                else
                    # take the median across GCMs
                    ts = dropdims(median(runoff_modeled[:, :, At(ssp), At(basin), At(model)], dims=:GCM), dims=:GCM);   
                end

                basin_runoff_monthly[:, At(basin), At(ssp), At(model)] .= ustrip(collect(mean.(groupby(ts, Ti => Bins(month, 1:12)))))
                
                gts = groupby(ts, Ti => year)
                ts_annual = ustrip(sum.(gts[length.(gts) .== 12]))
                fit = curve_fit(Altim.model_trend, 1:length(ts_annual),ts_annual,  Altim.p_trend)
                basin_runoff_trend[At(basin), At(ssp), At(model)] = fit.param[2] 

                basin_runoff_mean[At(basin), At(ssp), At(model)] = mean(ts_annual);
            end
        end
    end
end

"""
Generate histograms comparing differences in trends and means between GEMB and other models.

Creates histograms showing percentage differences between GEMB and Wimberly 2024 models
for both runoff trends and mean annual runoff values.
"""		
palette = (; color=reverse(Makie.resample_cmap(:Paired_4, length(dmodel))))

ssp = "ssp585"
scale = identity
for (k, da) in enumerate([basin_runoff_trend, basin_runoff_mean])
    f = Figure()
    if k == 1
        xlabel = "difference in trend [%]"
        outfile = joinpath(paths[:figures], "wimberly2024_trend_hist.png")
    else
        xlabel = "difference in mean annual runoff [%]"
        outfile = joinpath(paths[:figures], "wimberly2024_mean_hist.png")
    end
    print("\n#### $xlabel ###\n")
    ax = Axis(f[1, 1]; yscale=scale, xscale=scale, xlabel=xlabel, ylabel="count")

    a = da[:, At(ssp), At("GEMB")]
    for (i, model) in enumerate(dmodel[2:end])
        b = da[:, At(ssp), At(model)]
        delta = ((b-a) ./ a *100)
        println("standard deviation between GEMB and $model = $(round(Int,std(delta)))%")
        println("mean absolute differnce between GEMB and $model = $(round(Int,mean(abs.(delta))))%")
        stephist!(ax, delta; label=model, color=palette[1][i+1], linewidth=3)
    end
    axislegend(ax; position=:rt, framevisible=false)    
    display(f)

    save(outfile, f)
end

"""
Generate scatter plots comparing GEMB with other models for trends and means.

Creates scatter plots showing the relationship between GEMB estimates and 
Wimberly 2024 model estimates for both runoff trends and mean annual runoff.
"""
scale = identity
for (k, da) in enumerate([basin_runoff_trend, basin_runoff_mean])

    f = Figure()

    if k == 1
        limits = [-0.5, 3]*0.1
        units = "[Gt/yr²]"
        outfile = joinpath(paths[:figures], "wimberly2024_trend_scatter.png")
    else
        limits = [0, 6]
        units = "[Gt/yr]"
        outfile = joinpath(paths[:figures], "wimberly2024_rate_scatter.png")
    end

    r = 1
    c = 1
    for ssp in dssp
        
        xlabel = r == 2 ? "GEMB $units" : ""
        ylabel = c == 1 ? "Other $units" : ""

        ax = Axis(f[r, c]; title=ssp, xlabel=xlabel, ylabel=ylabel, yscale=scale, xscale=scale)

        x = collect(da[:, At(ssp), At("GEMB")])
        for (i, model) in enumerate(dmodel[2:end])
            y = collect(da[:, At(ssp), At(model)])
            scatter!(ax, x, y; label= model, color=palette[1][i+1])
        end

        lines!(ax, limits, limits; label=nothing, color=(:black, 0.5))

        if (r==1) && (c ==1)
            axislegend(ax; position=:lt, framevisible=false)
        end

        c += 1
        if c > 2
            c = 1
            r += 1
        end
        xlims!(ax, limits)
        ylims!(ax, limits)

    end

    display(f)
    save(outfile, f)
end

"""
Generate monthly runoff comparison plot between GEMB and other models.

Creates a line plot showing the seasonal pattern of glacier runoff for
GEMB and Wimberly 2024 models.
"""
begin
    total_runoff_monthly = dropdims(sum(basin_runoff_monthly, dims=:Basin), dims=:Basin)
    outfile = joinpath(paths[:figures], "wimberly2024_monthly_runoff.png")
    f = Figure()
    limits = [1,12]
    ax = Axis(f[1, 1]; xticks=1:12)

    ssp = dssp[1]
    for (i, model) in enumerate(dmodel)
        if model == "GEMB"
            lines!(ax, 1:12,collect(total_runoff_monthly[:, At(ssp), At(model)]), label=model, color=palette[1][i], linewidth=5)
        else
            lines!(ax, 1:12,collect(total_runoff_monthly[:, At(ssp), At(model)]), label=model, color=palette[1][i], linewidth=2)
        end
    end
    xlims!(ax, limits)

    axislegend(ax; position=:lt, framevisible=false)
    display(f)

    save(outfile, f)
end
