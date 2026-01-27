# Load packages and data paths
begin
    using CSV
    using DataFrames
    import GlobalGlacierAnalysis as GGA
    using DimensionalData
    using Dates
    using CairoMakie
    using Statistics
    using ProgressMeter
    using NCDatasets
    using FileIO
    using CFTime
    using GeoDataFrames
    
    # Import functions from utilities
    include("utilities_response2reviewers.jl")

    # data paths
    paths = (
        wgms_mb = "/mnt/bylot-r3/data/glacier_mb/WGMS/DOI-WGMS-FoG-2025-02b/data/mass_balance.csv",
        wgms_glacier = "/mnt/bylot-r3/data/glacier_mb/WGMS/DOI-WGMS-FoG-2025-02b/data/glacier.csv",

        # GEMB dataset paths
        gemb_corrected = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_NH_1979to2024_820_40_racmo_grid_lwt_e97_0_corrected_geotile_dv.jld2",
        gemb_original = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2024_820_40_racmo_grid_lwt_e97_0_geotile_dv.jld2",

        path2glacier = GGA.pathlocal[Symbol("glacier_individual")],
    )

    # Analysis parameters
    study_start_year = 2000
    min_neighbors_variogram = 3
    geotile_thresholds = [1, 20, 65]
    histogram_bounds = (-5, 3)
    historic_analysis_mincount_range = 5:15

    # RGI glacier counts for context
    rgi6_count = 215547
    rgi7_count = 274531
end



# --- Example: compare modeled runoff to external (Rounce 2022) dataset (disabled/block comment) ---
begin
    date_range = DateTime(2000, 1, 1).. DateTime(2025, 1, 1)
    scenario = "ssp126"
    var2extract = "glac_runoff_monthly"
    path2regional = joinpath(GGA.pathlocal[:project_dir], "Gardner2025_regional_timseries.nc")
    println("Gardner2025_regional_timseries.nc last modified: $(Dates.unix2datetime(mtime(path2regional)))")
    ds_gardner = GGA.netcdf2dimstack(path2regional)
    rain_gardner_include = true
    drgi = dims(ds_gardner[:fac],:rgi);
    dstudy = Dim{:study}(["gardner", "rounce"])
    runoff = zeros(drgi, dstudy)


    for rgi_region in drgi
    #rgi_region = 1
    # Load Rounce data and Gardner output for comparison
        begin
            folder = "/mnt/bylot-r3/data/glacier_mb/PyGEM/CMIP_global_2000_2100/$(lpad(string(rgi_region), 2, '0'))/"

            if !isdir(folder)
                continue
            end

            path2rounce = GGA.allfiles(folder; 
                fn_startswith="R$(lpad(string(rgi_region), 2, '0'))_glac_runoff_monthly_$(scenario)",
                fn_endswith="2000_2100_all.nc"
            )

            ds_example = NCDataset(path2rounce[1])
            dTi = Dim{:Ti}(reinterpret.(DateTime, ds_example[:time]))
            dmodel = Dim{:model}(ds_example[:Climate_Model])

            runoff_rounce =zeros(dTi)

            for path in path2rounce
                NCDataset(path) do ds_rounce0
                    drgi = Dim{:rgi}(ds_rounce0[:RGIId])
                    runoff_rounce0 = ds_rounce0[var2extract][:, :, :]::Array{Union{Missing,Float64},3}
                    runoff_rounce0 = DimArray(runoff_rounce0, (dTi, drgi, dmodel)) * 1E-9 # m^3/day to km^3
                    runoff_rounce0 = mean(runoff_rounce0, dims=:model)
                    runoff_rounce0 = coalesce.(runoff_rounce0, 0.0)
                    runoff_rounce0 = dropdims(sum(runoff_rounce0, dims=:rgi), dims=(:rgi, :model))
                    runoff_rounce .+= runoff_rounce0
                end
            end

            runoff_rounce = runoff_rounce[Ti = date_range]

            if rain_gardner_include
                runoff_gardner = ds_gardner[:runoff][date=date_range, rgi=At(rgi_region)] .+ ds_gardner[:rain][date=date_range, rgi=At(rgi_region)]
            else
                runoff_gardner = ds_gardner[:runoff][date=date_range, rgi=At(rgi_region)]
            end
            runoff_gardner .-= runoff_gardner[1]

            runoff_gardner = diff(runoff_gardner)

            # Plot comparison
            f = Figure();
            ax = f[1, 1] = Axis(f, ylabel="runoff [Gt]", title="RGI $(rgi_region)");
            lines!(ax, val(dims(runoff_rounce, :Ti)), collect(runoff_rounce); label="Rounce 2022")
            lines!(ax, val(dims(runoff_gardner, :date)), collect(runoff_gardner); label="Gardner 2025")
            f[1, 2] = Legend(f, ax)

            println("mean fractional difference rgi $(rgi_region): Rounce 2022 minus Gardner 2025: $(round(Int, ((sum(runoff_rounce) .- sum(runoff_gardner)) ./ sum(runoff_gardner)) * 100))%")

            runoff[rgi = At(rgi_region), study = At("gardner")] = sum(runoff_gardner)
            runoff[rgi = At(rgi_region), study = At("rounce")] = sum(runoff_rounce)
        end; 
        display(f)
    end

    runoff0 = deepcopy(runoff);
    # convert from total runoff in Gt to Gt/yr
    years = Dates.Day(maximum(date_range) - minimum(date_range)).value / 365.25
    runoff[study=At("gardner")] ./= years
    runoff[study=At("rounce")] ./= years

    # polulate HMA and Global runoff
    runoff[rgi = At(99), study = At("gardner")] = sum(runoff[study=At("gardner")])
    runoff[rgi = At(99), study = At("rounce")] = sum(runoff[study=At("rounce")])

    runoff[rgi = At(98), study = At("gardner")] = sum(runoff[study=At("gardner")][rgi = At(13:15)])
    runoff[rgi = At(98), study = At("rounce")] = sum(runoff[study=At("rounce")][rgi = At(13:15)])


    f = Figure();
    ax = Axis(f[1, 1]; ylabel="PyGEMv1.0.2", xlabel="This study");
    scatter!(ax, vec(runoff[study = At("gardner")]), vec(runoff[study = At("rounce")]));
    xmax = ceil(maximum(runoff[study=At("gardner")]), digits=-2)
    xlims!(ax, 0, xmax)
    ylims!(ax, 0, xmax)
    lines!(ax, [0, xmax], [0, xmax]);
    f

    delta = (runoff[study=At("rounce")] .- runoff[study=At("gardner")]) ./ runoff[study=At("gardner")] * 100

    rgis2plot = [19, 18, 17, 16, 98, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]

    f = GGA._publication_figure(columns=1, rows=1);
    ax = Axis(f[1, 1]; yticks = (1:length(rgis2plot), GGA.rginum2label.(rgis2plot)), xtickformat = "{:.0f}%", title = "difference in runoff between PyGEMv1.0.2 and this study")

    barplot!(ax, delta[rgi = At(rgis2plot)].data; direction=:x)
    display(f)
    save(joinpath(GGA.pathlocal.figures, "runoff_comparison.png"), f)

    println("PyGEMv1.0.2 produces $(round(Int,delta[rgi = At(99)]))% runoff compared to this study")
    gardner_minus_rounce = (runoff[study=At("gardner"), rgi = At(99)].-runoff[study=At("rounce"), rgi = At(99)]) ./ runoff[study=At("rounce"), rgi = At(99)] * 100
    println("This study produces $(round(Int,gardner_minus_rounce))% runoff compared to PyGEMv1.0.2")

    rounce_gardner_total_flux = (runoff[study=At("rounce"), rgi=At(99)] .- runoff[study=At("gardner"), rgi=At(99)]) ./ (runoff[study=At("gardner"), rgi=At(99)] + 155) * 100
    println("PyGEMv1.0.2 produces $(round(Int,rounce_gardner_total_flux))% total freshwater flux ompared to this study")
    gardner_rounce_total_flux = (runoff[study=At("gardner"), rgi=At(99)] .- runoff[study=At("rounce"), rgi=At(99)]) ./ (runoff[study=At("rounce"), rgi=At(99)] + 155) * 100
    println("This study produces $(round(Int,gardner_rounce_total_flux))% total freshwater flux ompared to PyGEMv1.0.2")
end


# Block for manually examining differences between corrected and original GEMB data (disabled by default)
if false
    # Load data volumes from GEMB datasets
    gemb_corrected = load(paths.gemb_corrected)
    gemb_original = load(paths.gemb_original)

    # Specify variable, spatial, temporal, and physical criteria to examine
    var2examine = "runoff"
    geotile2examine = "lat[+30+32]lon[+078+080]"
    dates2examine = DateTime(2000, 1, 1)..DateTime(2024, 12, 31)
    deltaheight2examine = -2000..0
    pscale2examine = 2.5

    # Slice data accordingly
    a = gemb_original[var2examine][geotile = At(geotile2examine), date = dates2examine, pscale = At(pscale2examine), ΔT = deltaheight2examine]
    b = gemb_corrected[var2examine][geotile = At(geotile2examine), date = dates2examine, pscale = At(pscale2examine), ΔT = deltaheight2examine]

    # Visualize raw vs corrected data and their fractional difference
    begin
        f = Figure(size=(1200, 400));
        ax1 = Axis(f[1, 1], title="original", ylabel="$(DimensionalData.name(dims(a)[2]))");
        heatmap!(ax1, a.data, colorrange=extrema(a));

        ax2 = Axis(f[1, 2], title="corrected");
        hm2 = heatmap!(ax2, b.data, colorrange=extrema(a));
        Colorbar(f[1, 3], hm2)

        ax4 = Axis(f[1, 4], title="fractional difference");
        frac_diff = (a .- b) ./ a
        hm4 = heatmap!(ax4, frac_diff.data, colorrange=extrema(frac_diff));
        Colorbar(f[1, 5], hm4)

        titlelayout = GridLayout(f[0, 1], halign = :left, tellwidth = false)
        Label(titlelayout[1, 1], "$geotile2examine - $var2examine: pscale = $pscale2examine", halign = :left, fontsize = 30, font = "TeX Gyre Heros Bold Makie")
    end; f

    # Alternate/quick visualization
    diff_frac = (a .- b) ./ a
    heatmap(diff_frac)    
end


# Plot global firn air content anomaly over time and save figure
begin
    filename0 = joinpath(GGA.pathlocal[:project_dir], "Gardner2025_regional_timseries.nc")
    regional_ts = GGA.netcdf2dimstack(filename0)

    rgi = 99
    var = :fac
    date_range = DateTime(2000, 1, 1) .. DateTime(2025, 1, 1)
    ts = regional_ts[var][rgi = At(rgi), date = date_range]

    f = GGA._publication_figure(columns=1, rows=1);
    ax = Axis(f[1, 1]; ylabel="volume [km³]", title="Global firn air content anomaly")
    lines!(ax, val(dims(ts, :date)), ts.data);
    xlims!(ax, minimum(val(dims(ts, :date)))-Month(1), maximum(val(dims(ts, :date)))+Month(1))
    display(f)
    save(joinpath(GGA.pathlocal.figures, "global_firn_air_content_anomaly.png"), f);
end

begin
    # Load WGMS data
    wgms_mb, wgms_glacier, ds, glaciers = load_wgms_data(paths);

    # Apply function to fill WGMS table with synthetic results
    add_runoff_and_dm!(wgms_mb, ds);

    # Select valid observations within the study period (after 2000)
    index = wgms_mb.both_valid .& (wgms_mb.begin_date .> DateTime(study_start_year, 1, 1));

    # Fill in missing glacier areas by joining to glacier polygons based on intersection
    @showprogress desc = "Adding missing areas..." Threads.@threads for r in eachrow(wgms_mb)[index .& ismissing.(wgms_mb.area)]
        ind = findfirst(GGA.GO.intersects.(glaciers.geometry, Ref(r.geometry)))
        if !isnothing(ind)
            r.area = round(glaciers[ind, :Area] * 1e6)
        end
    end;

    # Write one validation point per unique glacier to geopackage (point locations only)
    begin
        gdf = groupby(wgms_mb[index, :], :glacier_name)
        df2 = DataFrames.combine(gdf, :geometry => first => :geometry, :geometry => length => :count)
        GGA.GeoDataFrames.write(joinpath(GGA.pathlocal[:project_dir], "wgms_mb_validation_points.gpkg"),
                            df2; crs=GGA.EPSG(4326))
    end;

    # Attach glacier geometry to validation points and write to geopackage
    begin
        rename!(df2, :geometry => :point_location)
        df2[!, :geometry] .= repeat([glaciers.geometry[1]], nrow(df2))

        @showprogress desc = "Adding glacier geometry to validation points..." Threads.@threads for r in eachrow(df2)
            ind = findfirst(GGA.GO.intersects.(glaciers.geometry, Ref(r.point_location)))
            if !isnothing(ind)
                r.geometry = glaciers.geometry[ind]
            end
        end

        # Exclude points with no valid geometry and write result
        no_intersection = df2[:, :geometry] .== Ref(glaciers.geometry[1])
        GGA.GeoDataFrames.write(joinpath(GGA.pathlocal[:project_dir], "wgms_mb_validation.gpkg"),
            df2[.!no_intersection, Not(:point_location)]; crs=GGA.EPSG(4326))
    end;

    # --- Peyto Glacier investigation (example, disabled) ---
    #=
    # Identify rows corresponding to Peyto Glacier after 1999
    index_peyto = (wgms_mb.glacier_id .== 57) .& (wgms_mb.year .> 1999) 

    # Compute modeled and observed cumulative mass balance and visualize
    area = ds.area[X=Near(wgms_mb.longitude[findlast(index_peyto)]), Y=Near(wgms_mb.latitude[findlast(index_peyto)])]::Float64
    dm = ds.dm[X=Near(wgms_mb.longitude[findlast(index_peyto)]), Y=Near(wgms_mb.latitude[findlast(index_peyto)]), Ti=Touches(DateTime(2000, 1, 1), DateTime(2025, 1, 1))]::DimVector{Float64}
    dm_cumulative = cumsum(dm) ./ area / 1000
    f = Figure();
    ax1 = Axis(f[1, 1]);
    scatterlines!(wgms_mb[index_peyto, :year], cumsum(wgms_mb[index_peyto, :annual_balance]); label="in situ");
    scatterlines!(GGA.decimalyear.(dims(dm_cumulative, :Ti).val), dm_cumulative.data; label="modeled");
    f
    =#

    wgms_mean = mean(wgms_mb[index, :].annual_balance)
    model_mean = mean(wgms_mb[index, :].dm)

    println("FoG glaciers have a mean of $(round(wgms_mean, digits=2)) m w.e. and the extracted model values have a mean of $(round(model_mean, digits=2)) m w.e. ")

    # ---- Summary statistics and coverage for WGMS mass balance data ----

    # Calculate summary statistics
    index_all, index_area, df2 = calculate_wgms_summary_statistics(wgms_mb, index);

    # Calculate and plot area coverage
    wgms_annual_coverage = calculate_annual_coverage(df2, ds);
    f_area = plot_area_coverage(df2, wgms_annual_coverage)
    save(joinpath(GGA.pathlocal.figures, "wgms_mb_area_coverage.png"), f_area)

    # Calculate and plot glacier count statistics
    df2_glacier = calculate_glacier_count_statistics(wgms_mb, index_area);
    f_count = plot_glacier_count_timeseries(df2_glacier);
    save(joinpath(GGA.pathlocal.figures, "wgms_mb_glacier_count.png"), f_count)

    # ---- Calculate spatial variogram of annual balance ----
    variogram, ddistance = calculate_spatial_variogram(wgms_mb; mincount=min_neighbors_variogram);

    # Compute and summarize geotile misfit statistics
    df2 = calculate_geotile_misfit_statistics(wgms_mb, index)
    radius, synth_minus_obs_mad = calculate_misfit_summary_statistics(wgms_mb, index, df2);

    # Plot variogram
    f_var = plot_variogram(ddistance, variogram, radius, synth_minus_obs_mad)
    display(f_var)
    save(joinpath(GGA.pathlocal.figures, "wgms_mb_variogram.png"), f_var)

    println("at 100km observations have a mean absolute difference of $(round(variogram[distance = At(100_000)], digits=2)) m w.e. yr⁻¹ and the mean absolute difference between in situ and this study is $(round(synth_minus_obs_mad, digits=2)) m w.e. yr⁻¹")

    # Coverage calculation for each year (repeat, probably for reporting)
    index_area = index .& .!ismissing.(wgms_mb.area)
    gdf = groupby(wgms_mb[index_area, :], :year)
    df2 = DataFrames.combine(gdf, :area => sum => :area_sum, :annual_balance => std => :annual_balance_std, :synth_minus_obs => std => :synth_minus_obs_std)
    wgms_annual_coverage = df2.area_sum ./ sum(ds.area[:]) * 100
    println("average annual coverage: $(round(mean(wgms_annual_coverage), digits=2))%")

    # --- 2D histogram (heatmap) of observed vs modeled mass balance ---
    f_hist, density = create_2d_histogram(wgms_mb[index, :]; bound_lower=histogram_bounds[1], bound_upper=histogram_bounds[2], colormap=:tempo);
    save(joinpath(GGA.pathlocal.figures, "wgms_mb_2d_histogram.png"), f_hist)

    # --- Histograms for misfit statistics by season/period ---
    create_misfit_histograms(wgms_mb, index)
    save(joinpath(GGA.pathlocal.figures, "wgms_mb_misfit_histograms.png"), f_hist)
end

# load hugonnet data
# load in all regions
begin
    hugonnet = DataFrame()
    for rgi = 1:19
        rgi_str = lpad(string(rgi), 2, '0')
        path2hugonnet = "/mnt/bylot-r3/data/glacier_mb/Hugonnet2021/time_series/time_series_$(rgi_str)/dh_$(rgi_str)_rgi60_pergla_rates.csv"
        hugonnet = vcat(hugonnet, CSV.read(path2hugonnet, DataFrame))
    end

    # find annual balance and format similar to wgms_mb for spatial variogram
    hugonnet[!, "start_date"] = DateTime.(getindex.(split.(hugonnet.period, "_"), 1))
    hugonnet[!, "end_date"] = DateTime.(getindex.(split.(hugonnet.period, "_"), 2))

    hugonnet[!, "end_year"] = year.(hugonnet[!, "end_date"])
    hugonnet[!, "start_year"] = year.(hugonnet[!, "start_date"])
    hugonnet[!, "dt_years"] = hugonnet[!, "end_year"] - hugonnet[!, "start_year"]
    hugonnet[!, "year"] = hugonnet[!, "end_year"]
    index = hugonnet[!, "dt_years"] .== 1
    hugonnet = hugonnet[index, :]
    hugonnet[!, :annual_balance] = hugonnet.dhdt * 0.85


    rgi6 = GeoDataFrames.read(GGA.pathlocal.glacier_individual)
    index = indexin(hugonnet.rgiid, rgi6.RGIId)

    # it seems some ids do not match... I think this is because hugonnet updated glaciers that are not in rgi6
    index_not_match = isnothing.(index)
    hugonnet = hugonnet[.!index_not_match, :]

    hugonnet[!, :longitude] = rgi6.CenLon[index[.!index_not_match]]
    hugonnet[!, :latitude] = rgi6.CenLat[index[.!index_not_match]]

subset_random = rand(1:nrow(hugonnet), round(Int, nrow(hugonnet)/100))
variogram, ddistance = calculate_spatial_variogram(deepcopy(hugonnet[subset_random, [:annual_balance, :longitude, :latitude, :year]]); mincount=min_neighbors_variogram);

f_var_hugonnet = plot_variogram(ddistance, variogram, 0, 0)
display(f_var_hugonnet)
save(joinpath(GGA.pathlocal.figures, "wgms_mb_variogram_hugonnet.png"), f_var_hugonnet)

end


# link hugonnet glaciers to wgms_mb
index_into_wgms = indexin(hugonnet.rgiid, wgms_glacier.rgi60_ids)
index_valid = .!isnothing.(index_into_wgms)
hugonnet[!, :wgms_id] .= 0
hugonnet[index_valid, :wgms_id] = wgms_glacier.id[index_into_wgms[index_valid]]


wgms_mb[!, :annual_balance_hugonnet] .= NaN
for r in eachrow(wgms_mb)
    index = findfirst((hugonnet.dt_years .== 1) .& (hugonnet.wgms_id .== r.glacier_id) .& (hugonnet.year .== r.year))
    if !isnothing(index)
        r.annual_balance_hugonnet = hugonnet.dhdt[index] .* 0.85
    end
end

delta = wgms_mb.annual_balance_hugonnet .- wgms_mb.annual_balance
index_valid = (.!ismissing.(delta)) .& (.!isnan.(delta))
f = Figure();
ax = Axis(f[1, 1]; xlabel="wgms minus hugonnet [m w.e. yr⁻¹]", ylabel="count");
hist!(ax, delta[index_valid])
f



# Create precipitation raster of maximum monthly precipitation
begin
    path2precip = "/mnt/bylot-r3/data/glacier_mb/GPCP/precip.mon.ltm.1991-2020.nc"

    precip = GGA.Rasters.Raster(path2precip, name=:precip);
    precip_max = dropdims(getindex.(argmax(precip, dims=3), 3), dims=:Ti)

    output_path = joinpath(splitpath(path2precip)[1:end-1]..., "precip.mon_of_max.ltm.1991-2020.tif")
    precip_max = reverse(precip_max, dims=:Y)

    precip_max = rebuild(precip_max, data = UInt8.(precip_max.data))
    precip_max = set(precip_max, X => dims(precip_max, X).val .-180)

    GGA.Rasters.write(output_path, precip_max, force=true)
end;



