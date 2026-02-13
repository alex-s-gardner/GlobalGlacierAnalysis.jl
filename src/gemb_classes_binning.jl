# =============================================================================
# Process GEMB (Glacier Energy and Mass Balance) model output for global glacier analysis.
# =============================================================================
#
# This script:
# 1. Loads and combines raw GEMB output from multiple simulations
# 2. Organizes data into geotiles with consistent spatial and temporal dimensions  
# 3. Fills data gaps through interpolation and extrapolation across elevation bands
# 4. Extends model results to cover additional precipitation scaling factors
#
# The processed data is saved at multiple stages to enable efficient reuse and analysis.
# GEMB output is in units of meters ice equivalent (m i.e.) assuming an ice density of 910 kg/m³.

begin
    import GlobalGlacierAnalysis as GGA
    using Dates
    using FileIO
    using ProgressMeter
    using DimensionalData
    using CairoMakie
    using Statistics
   
    # mirror files from Nicole's home directory to RAID storage
    mirror_raw_files = false;

    single_geotile_test = nothing #GGA.geotiles_golden_test[1] #"lat[+60+62]lon[-142-140]"; GGA.geotiles_golden_test[1]

    # run parameters 
    project_id = :v01;
    geotile_width = 2;
    surface_mask = :glacier
    height_bins_extrapolate = 1; # number bins above and below the valid height range to extrapolate with each search iteration
    maximum_extrap_fraction = 0.1; # if there is < maximum_extrap_fraction of the fraction of area below/above the valid height range, then extrapolate the height range to the bottom/top of the glacier
    minimum_extrap_area = 10; # if there is <  minimum_extrap_area of the area below/above the valid height range, then extrapolate the height range to the bottom/top of the glacier
    maximum_search_iterations = 100; # maximum number of search iterations to populate the elevation profile
    minimum_land_coverage_fraction = 0.70; # minimum land coverage fraction for acceptance of GEMB data
    search_buffer_increment = 100_000; # increment of search buffer in meters with each search iteration [m]
    show_interp_extrap_plots = false;
    show_interp_extrap_stats = false;
    gemb_run_id = 5;

    elevation_classes_method = :mscale #:mscale # [:none, :Δelevation, :mscale]

    # exclude derived variables of smb and runoff (these are calculated later to ensure mass conservation after interpolation)
    vars2extract = ["fac", "acc", "refreeze", "melt", "rain", "ec"]
    dims2extract = ["latitude", "longitude", "date", "height"]

    gembinfo = GGA.gemb_info(; gemb_run_id);
 
    # define date and hight binning ranges 
    date_range, date_center = GGA.project_date_bins()

    # expand daterange to 1940 by make sure to match exisiting project ranges 
    date_end_new = Date(1970,1,1)
    Δd = 30
    date_range = reverse(last(date_range):-Day(Δd):date_end_new)
    date_center = date_range[1:end-1] .+ Day(Δd / 2)
    ddate = Dim{:date}(date_center)

    length(date_range) == (length(date_center) +1) ? nothing : error("date_range and date_center do not have the same length")

    # load geotile definitions with corresponding hypsometry
    geotiles = GGA._geotile_load_align(; surface_mask, geotile_order=nothing, only_geotiles_w_area_gt_0=true)
    dgeotile = Dim{:geotile}(geotiles.id)

    area_km2 = GGA._geotile_area_km2(; surface_mask, geotile_width)
  
    if .!isnothing(single_geotile_test)
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)], OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!"
    end
end;

# move files from Nicole's home directory to RAID storage
if mirror_raw_files
    run(`bash -c "rsync -av /home/schlegel/Share/GEMBv1/SH_sample/*_corrected_* /mnt/bylot-r3/data/gemb/mat/lw_correction/"`)
    run(`bash -c "rsync -av /home/schlegel/Share/GEMBv1/NH_sample/*_corrected_* /mnt/bylot-r3/data/gemb/mat/lw_correction/"`)
    run(`bash -c "rsync -av --exclude='*_corrected_*' /home/schlegel/Share/GEMBv1/SH_sample/ /mnt/bylot-r3/data/gemb/mat/no_lw_correction/"`)
    run(`bash -c "rsync -av --exclude='*_corrected_*' /home/schlegel/Share/GEMBv1/NH_sample/ /mnt/bylot-r3/data/gemb/mat/no_lw_correction/"`)
end

# Sanity check block to visualize and inspect GEMB .mat data for typical errors/anomalies.
sanity_check = false;
if sanity_check

    # Gather all GEMB .mat file paths according to filtering criteria
    gemb_files = vcat(GGA.allfiles.(gembinfo.gemb_folder; subfolders=false, fn_endswith=".mat", fn_contains=gembinfo.file_uniqueid)...)

    # Pick the first file whose name contains "p2_t2" (a particular experiment)
    fn = gemb_files[findfirst(occursin.("p2_t2", gemb_files))]
    file = GGA.matopen(fn)
    foo = GGA.MAT.read(file)

    # plot land sea mask with accepted and rejected points 
    land_sea_mask = GGA.Raster(GGA.pathlocal.era5_land_sea_mask)[:,:,1];
    d = 0.25;
    y = X(90:-d:-90; sampling=DimensionalData.Intervals(DimensionalData.Center()));
    x = Y(0:d:359.75); sampling=DimensionalData.Intervals(DimensionalData.Center());

    xpt =foo["lon"][:];
    index = xpt .< 0;
    xpt[index] .+= 360;
    ypt = foo["lat"][:]
  
    pt = tuple.(X.(Near.(xpt)), Y.(Near.(ypt)));
    val = land_sea_mask[pt];

    land_mask = DimArray(land_sea_mask.data', (y,x))
    figure = Figure();
    markersize=3;
    ax = Axis(figure[1, 1]);
    rotate!(ax.scene, +pi/2);
    p = heatmap!(ax, land_mask);
    scatter!(ax,ypt, xpt;color=:red, markersize);

    index = val .>=minimum_land_coverage_fraction;
    scatter!(ax,ypt[index], xpt[index] ; color=:green, markersize);
    display(figure);

    println("######### plotting sanity check for a singe file ######################## ")
    println(fn)
    println("######################################################## ")

    # Scatter plot of longitude vs. latitude of all points in foo
    p = plot(vec(foo["lon"]), vec(foo["lat"]));

    # Identify valid rows: any non-NaN and nonzero 'Melt' values across time (dims=2)
    index_valid = vec(any(.!isnan.(foo["Melt"]) .& (foo["Melt"] .!= 0), dims=2))
    # Overlay in red those locations with any valid 'Melt' values
    plot!(vec(foo["lon"][index_valid]), vec(foo["lat"][index_valid]); color=:red)
    display(p)

    # Loop through all loaded keys (variables) in the mat file
    for k in keys(foo)
        # Only plot if the variable matches the size of "Rain" (assume same grid & time structure)
        if length(foo[k]) == length(foo["Rain"])
            f = Figure()
            ax = Axis(f[1, 1]; title="$k")
            # Plot the first valid location's time series for this variable
            lines!(ax, vec(foo["time"]), vec(foo[k][index_valid, :][1, :]))
            display(f)
        end
    end

    # plot dv
    begin
        point_2_select = 7;
        ind = findall(GGA.within.(Ref(GGA.geotile_extent(GGA.geotiles_golden_test[1])), foo["lon"][:], foo["lat"][:]))[point_2_select]
        dv = vec(cumsum(foo["Accumulation"][ind, :] .- foo["Melt"][ind, :] .+ foo["Refreeze"][ind, :] .- foo["EC"][ind, :]) ./ 1000 .+ foo["FAC"][ind, :])
        figure = Figure();
        ax = Axis(figure[1, 1],title="dv: lat = $(round(foo["lat"][ind], digits=3)), lon = $(round(foo["lon"][ind], digits=3))");
        lines!(ax, vec(foo["time"]), dv; )
        display(figure);
    end

    begin
        gemb0 = GGA.gemb_read2(fn; datebin_edges = GGA.decimalyear.(date_range))

        point_2_select = 7;
        ind = findall(GGA.within.(Ref(GGA.geotile_extent(GGA.geotiles_golden_test[1])), gemb0["longitude"][:], gemb0["latitude"][:]))[point_2_select]
        dv = vec(gemb0["acc"][ind, :] .- gemb0["melt"][ind,:] .+ gemb0["refreeze"][ind,:] .- gemb0["ec"][ind,:] .+ gemb0["fac"][ind,:])
        figure = Figure();
        ax = Axis(figure[1, 1],title="dv: lat = $(round(gemb0["latitude"][ind], digits=3)), lon = $(round(gemb0["longitude"][ind], digits=3))");
        lines!(ax, vec(gemb0["date"]), dv; )
        display(figure);
    end

end

# =============================================================================
# GEMB DATA MERGING AND GEOTILE PROCESSING
# =============================================================================
# Merges GEMB model data from multiple .mat files into a single combined file
# and organizes data into geotiles with elevation and precipitation classes.
# 
# Process:
# 1. Collect and read all GEMB .mat files (~3 min)
# 2. Standardize coordinates and merge simulations  
# 3. Organize data into geotiles with elevation/precipitation bins (~3.5 min)
# 4. Apply spatial buffering to ensure minimum coverage requirements
# 5. Save processed data to JLD2 files for each geotile
# 
# Output: Combined GEMB data saved to JLD2 file with geotile structure
# 
# Parameters:
# - gemb_files: Vector of GEMB .mat file paths
# - gembinfo: GEMB information structure containing file metadata
# - geotiles: DataFrame with geotile definitions and hypsometry
# - search_buffer: Buffer distance in meters for spatial coverage
# - min_gemb_coverage: Minimum required GEMB data coverage (default: 0.75)
# - force_remake_before: DateTime threshold for forcing file regeneration
# - single_geotile_test: Optional single geotile ID for testing

begin
    gemb_files = vcat(GGA.allfiles.(gembinfo.gemb_folder; subfolders=false, fn_endswith=".mat", fn_contains=gembinfo.file_uniqueid)...)

    # ensure expected number of files found

    if any(occursin.("schlegel", gembinfo.gemb_folder))
        expected_number_of_files = length(gembinfo.elevation_delta) .* length(gembinfo.precipitation_scale) * length(gembinfo.gemb_folder); # multiply by 2 for NH and SH
    else
        expected_number_of_files = length(gembinfo.elevation_delta) .* length(gembinfo.precipitation_scale) * 2; # multiply by 2 for NH and SH
    end

    if length(gemb_files) != expected_number_of_files
        error("Expected $(expected_number_of_files) files but found $(length(gemb_files)): check that file_uniqueid in utilites_project.jl is not specific to a single hemisphere")
    end

    # units of m i.e. [m of air for fac]
    gembX = GGA.read_gemb_files(gemb_files, gembinfo; vars2extract=vcat(dims2extract, vars2extract), date_range, date_center, path2land_sea_mask=GGA.pathlocal.era5_land_sea_mask, minimum_land_coverage_fraction)

    # check that elvation and precipitaiton classes of the raw data look correct
    if !isnothing(single_geotile_test)

        point_2_select = 7;
        index_geotile = GGA.within.(Ref(GGA.geotile_extent(single_geotile_test)), gembX["longitude"], gembX["latitude"])

        index = findall(index_geotile .& (gembX["precipitation_scale"] .== 1.0) .& (gembX["elevation_delta"] .== 0))[point_2_select]
     
        index_point = (gembX["latitude"] .== gembX["latitude"][findall(index_geotile)[point_2_select]]) .& (gembX["longitude"] .== gembX["longitude"][findfirst(index_geotile)])

        pscale = gembX["precipitation_scale"]
        Δelevation = gembX["elevation_delta"]

        index = (pscale .== 1.0) .& index_point

        cmap = Makie.resample_cmap(:thermal, length(elevation_delta[index])+1);
        for k in keys(gembX)
            # Only plot if the variable matches the size of "Rain" (assume same grid & time structure)
            if length(gembX[k]) == length(gembX["fac"])
                f = Figure();
                ax = Axis(f[1, 1], title="gemb (lat = $(round(gembX["latitude"][findfirst(index_geotile)], digits=3)), lon = $(round(gembX["longitude"][findfirst(index_geotile)], digits=3))): $k");
                for edelta in [0] #sort(elevation_delta[index])
                    index0 = findfirst((elevation_delta .== edelta) .& index)
                    lines!(ax, gembX["date"], gembX[k][index0,:]; label="$edelta", color=cmap[findfirst(edelta .== sort(elevation_delta[index]))])
                end
            f[1, 2] = Legend(f, ax, "Δelevation [m]", framevisible = false)
            display(f)
            end
        end

        index = (Δelevation .== 0) .& index_point
        for k in keys(gembX)
            # Only plot if the variable matches the size of "Rain" (assume same grid & time structure)
            if length(gembX[k]) == length(gembX["fac"])
                f = Figure();
                ax = Axis(f[1, 1], title="gemb (lat = $(round(gembX["latitude"][findfirst(index_geotile)], digits=3)), lon = $(round(gembX["longitude"][findfirst(index_geotile)], digits=3))): $k");
                for pscale0 in [1] #sort(pscale[index])
                    index0 = findfirst((pscale .== pscale0) .& index)
                    lines!(ax, gembX["date"], gembX[k][index0,:]; label="$pscale0", color=cmap[findfirst(pscale0 .== sort(pscale[index]))])
                end
            f[1, 2] = Legend(f, ax, "pscale", framevisible=false)
            display(f)
            end
        end

        f = Figure();
        ax = Axis(f[1, 1], title="gemb (lat = $(round(gembX["latitude"][findfirst(index_geotile)], digits=3)), lon = $(round(gembX["longitude"][findfirst(index_geotile)], digits=3))): dv");
        for pscale0 in [1] #sort(pscale[index])
            index0 = findfirst((pscale .== pscale0) .& index)
            dv =  gembX["acc"][index0,:] .- gembX["fac"][index0,:] .- gembX["melt"][index0,:] .+ gembX["refreeze"][index0,:] .- gembX["ec"][index0,:]
            lines!(ax, gembX["date"], dv; label="$pscale0", color=cmap[findfirst(pscale0 .== sort(pscale[index]))])
        end
        display(f)



        index0 = findall(index_geotile .& (gembX["precipitation_scale"] .== 1.0) .& (gembX["elevation_delta"] .== 0))[point_2_select]
        f = Figure();
        ax = Axis(f[1, 1], title="gemb (lat = $(round(gembX["latitude"][index0], digits=3)), lon = $(round(gembX["longitude"][index0], digits=3))): dv");
        dv =  gembX["acc"][index0,:] .- gembX["fac"][index0,:] .- gembX["melt"][index0,:] .+ gembX["refreeze"][index0,:] .- gembX["ec"][index0,:]
        lines!(ax, gembX["date"], dv;)
        display(f)

    end

    if !isnothing(single_geotile_test)
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)], OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!"
        geotiles = geotiles[geotiles.id .== single_geotile_test, :]
    end

    # this takes about 3 to 8 minutes depending on method used
    gemb_dv0 = GGA.process_gemb_geotiles(
        gembX,
        geotiles;
        maximum_search_iterations,
        height_bins_extrapolate,
        maximum_extrap_fraction,
        search_buffer_increment,
        show_interp_extrap_plots,
        show_interp_extrap_stats,
        single_geotile_test,
        elevation_classes_method,
    );

    # there should be no NaN in any geotiles for Date(2000,1,1)
    for k in keys(gemb_dv0)
        if any(isnan.(gemb_dv0[k][date = Near(Date(2000,1,1))]))
            error("NaN found in $(k) for nearest Date(2000,1,1)")
        end
    end

    # ensure that there are not negative values for melt rates.
    for k in [:refreeze, :rain, :acc, :melt, :runoff]
        diff_melt = diff(gemb_dv0[k], dims=:date)
        if any(diff_melt .< -1E-9)
            error("!!! Negative values found in $(k) !!!!")
        end
    end

    if any(diff(gemb_dv0[:refreeze], dims=:date) .- diff(gemb_dv0[:melt], dims=:date) .> 1E-9)
        error("!!! refreeze is greater than melt !!!!")
    end
    
    if !isnothing(single_geotile_test)
        # plot the first variable
        dpscale =  dims(gemb_dv0, :pscale)
        dmscale = dims(gemb_dv0, :mscale)
        cmap = Makie.resample_cmap(:thermal, max(length(dpscale), length(dmscale)))

        if elevation_classes_method == :mscale
                dmscale_center = 1
        else
                dmscale_center = 0
        end

        for k in keys(gemb_dv0)
            f = GGA._publication_figure(; columns=1, rows=1)
            ax = Axis(f[1, 1]; title="$single_geotile_test [pscale = 1]: $k")
        
            for pscale in [1] #dpscale
                for mscale in dmscale
                    lines!(ax, gemb_dv0[k][geotile = At(single_geotile_test), pscale = At(pscale), mscale = At(mscale)]; label="mscale: $mscale", color=cmap[findfirst(mscale .== dmscale)])
                end
            end

            axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
            display(f)

            for k in keys(gemb_dv0)
                f = GGA._publication_figure(; columns=1, rows=1)
                ax = Axis(f[1, 1]; title="$single_geotile_test [mscale = $dmscale_center]: $k")

                for pscale in dpscale
                    for mscale in [dmscale_center]
                        lines!(ax, gemb_dv0[k][geotile=At(single_geotile_test), pscale=At(pscale), mscale=At(mscale)]; label="pscale: $pscale", color=cmap[findfirst(pscale .== dpscale)])
                    end
                end
                axislegend(ax, position=:lt, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false)
                display(f)
            end
        end
    end

    if isnothing(single_geotile_test)
        # construct output filename
        gemb_geotile_filename_dv = replace(gembinfo.filename_gemb_combined, ".jld2" => "_geotile_dv.jld2");
        save(gemb_geotile_filename_dv, "gemb_dv", gemb_dv0)
    end
end