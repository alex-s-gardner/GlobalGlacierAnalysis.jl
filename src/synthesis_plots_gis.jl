"""
    synthesis_plots_gis.jl

Generate plots and GIS files from glacier mass change synthesis data.

This script processes the reference synthesis run to extract trends and amplitudes
of glacier elevation and volume changes. It creates visualization outputs and
exports the data to GIS formats for further analysis.

The workflow:
1. Loads the reference synthesis run data
2. Calculates temporal trends for all variables
3. Exports results to Arrow and GeoPackage formats
4. Creates histograms comparing GEMB and altimetry-derived trends and amplitudes
5. Visualizes parameter distributions (pscale and Δheight)

Key outputs:
- GIS files with geotile-level rates in km³/yr
- Histograms comparing different measurement approaches
- Visualizations of model parameter distributions
"""

begin
    using GlobalGlacierAnalysis
    using FileIO
    using Dates
    using GeoDataFrames
    using CairoMakie
    import GeoFormatTypes as GFT
    using DataFrames


    reference_run = "binned/2deg/glacier_dh_best_meanmadnorm5_v01_filled_ac_p2_synthesized.jld2"

    paths = GlobalGlacierAnalysis.pathlocal
end

# export trends and amplitudes for plotting of reference_run only
begin     
    path2reference = joinpath(paths.data_dir, reference_run)
    binned_synthesized_dv_file = replace(path2reference, ".jld2" => "_gembfit_dv.jld2")
   
    geotiles0 = FileIO.load(binned_synthesized_dv_file, "geotiles")

    vars_no_write = setdiff(names(geotiles0), ["id", "glacier_frac", "landice_frac", "floating_frac", "geometry", "group", "pscale", "Δheight", "mie2cubickm", "rgiid"])
    vars_ts = setdiff(vars_no_write, ["extent", "area_km2"])
    
    # Fit temporal trends to all variables
    source_crs1 = GFT.EPSG(4326)
    geotiles0 = GlobalGlacierAnalysis.df_tsfit!(geotiles0, vars_ts; datelimits = (DateTime(2000,1,1), DateTime(2023,1,1)))
    geotiles0[!,:area_km2] = sum.(geotiles0[!,:area_km2])
    outfile = joinpath(paths.data_dir, "project_data", "geotiles_rates_km3yr.arrow");
    GeoDataFrames.write(outfile, geotiles0[:, Not(vars_no_write)]; crs=source_crs1)

    # Export results, excluding raw time series
    isvec = []
    for i = 1:ncol(geotiles0)
        push!(isvec, typeof(geotiles0[1,12]) >: Vector)
    end

    GeoDataFrames.write(joinpath(paths.data_dir, "project_data", "geotiles_rates_km3yr.gpkg"), geotiles0[:, Not(vars_no_write)]; crs=source_crs1)

    # now do the same but for area averaged rates of change 
    for varname in vars_ts
        geotiles0[:, varname] ./= geotiles0[:, :mie2cubickm]
    end

    # Fit temporal trends to all variables (easier than finding all of the fits and then multiplying by mie2cubickm)
    geotiles0 = GlobalGlacierAnalysis.df_tsfit!(geotiles0, vars_ts; datelimits = (DateTime(2000,1,1), DateTime(2025,1,1)))

    # plot a histogram of the rates
    f = Figure();
    var1 = "dv_trend"
    var2 = "dv_altim_trend"
    linewidth = 3
    ax1 = f[1, 1] = Axis(f; xlabel = "trend [m i.e. yr⁻¹]", ylabel = "count")
    CairoMakie.stephist!(geotiles0[:, var1], bins = -5:0.25:5; label = "gemb", linewidth)
    CairoMakie.stephist!(geotiles0[:, var2], bins = -5:0.25:5; label = "GlobalGlacierAnalysis", linewidth); 
    axislegend(ax1, framevisible = false); 
    
    var1 = "dv_amplitude"
    var2 = "dv_altim_amplitude"
    ax2 = f[1, 2] = Axis(f; xlabel = "amplitude [m i.e.]", ylabel = "count")
    CairoMakie.stephist!(geotiles0[:, var1], bins = 0:0.25:5; label = "gemb", linewidth)
    CairoMakie.stephist!(geotiles0[:, var2], bins = 0:0.25:5; label = "GlobalGlacierAnalysis", linewidth); 
    axislegend(ax2, framevisible = false); 
    display(f)

    outfile = joinpath(paths[:figures], "hist_gemb_altim_trend_amplitude.png")
    CairoMakie.save(outfile, f)
 
    # plot a histogram of pscale and Δheight
    synthesized_gemb_fit = replace(path2reference, ".jld2" => "_gemb_fit.arrow")
    gemb_fit = GeoDataFrames.read(synthesized_gemb_fit);

    f = Figure();
    ax1 = f[1, 1] = Axis(f; xlabel = "pscale", ylabel = "count")
    CairoMakie.hist!(gemb_fit[:, :pscale], bins = 0.25:0.25:4, linewidth)
    ylims!(low = 0)
    ax2 = f[1, 2] = Axis(f; xlabel = "Δheight [m]", ylabel = "count")
    CairoMakie.hist!(gemb_fit[:, :Δheight], bins = -3000:100:3000, linewidth); 
    ylims!(low = 0)
    #axislegend(ax1, framevisible = false); 
    display(f)

    outfile = joinpath(paths[:figures], "hist_pscale_dheight.png")
    CairoMakie.save(outfile, f)

    # rsync -rav devon:/mnt/bylot-r3/data/project_data/ ~/data/Altim/project_data/
    # Sync command for reference
end
