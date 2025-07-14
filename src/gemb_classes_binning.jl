# =============================================================================
# Process GEMB (Glacier Energy and Mass Balance) model output for global glacier analysis.
# =============================================================================
#
# This script:
# 1. Loads and combines raw GEMB output from multiple simulations
# 2. Organizes data into geotiles with consistent spatial and temporal dimensions  
# 3. Fills data gaps through interpolation and extrapolation across elevation bands
# 4. Extends model results to cover additional precipitation scaling factors
# 5. Extrapolates data temporally using climatological means
#
# The processed data is saved at multiple stages to enable efficient reuse and analysis.
# GEMB output is in units of meters ice equivalent (m i.e.) assuming an ice density of 910 kg/m³.

#begin
begin
    import GlobalGlacierAnalysis as GGA
    using Dates
    using FileIO
    using ProgressMeter
    using DimensionalData
   
    single_geotile_test = GGA.geotiles_golden_test[1]
    geotiles2plot = nothing; #[GGA.geotiles_golden_test[1]]

    # run parameters 
    force_remake_before = DateTime(2025, 7, 4, 10, 0, 0) + GGA.local2utc
    project_id = :v01;
    geotile_width = 2;
    surface_mask = :glacier
    geotile_buffer = 50000 # distance in meters outside of geotiles to look for data
    gemb_run_id=4;

    # exclude derived variables of smb and runoff (these are calculated later to ensure mass conservation after interpolation)
    vars2extract = ["fac", "acc", "refreeze", "melt", "rain", "ec"]
    dims2extract = ["latitude", "longitude", "date", "height"]

    # min gemb coverage 
    min_gemb_coverage =  0.75

    gembinfo = GGA.gemb_info(; gemb_run_id)
 
    #Δheight simulates changing model elevation to increase / decrease melt, this is done in the regional Δvolume calculation
    height_center = GGA.project_height_bins()[2]
    height_bin_interval = height_center[2] - height_center[1]
    dΔheight = Dim{:Δheight}(-3000:height_bin_interval:3000) 

    # define date and hight binning ranges 
    date_range, date_center = GGA.project_date_bins()

    # expand daterange to 1940 by make sure to match exisiting project ranges 
    foo = collect(date_range[1]:-Day(30):DateTime(1940))
    date_range = foo[end]:Day(30):date_range[end]
    date_center = date_range[1:end-1] .+ Day(Day(30) / 2)
    ddate = Dim{:date}(date_center)

    # load geotile definitions with corresponding hypsometry
    geotiles = GGA._geotile_load_align(; surface_mask, geotile_order=nothing, only_geotiles_w_area_gt_0=true)
    dgeotile = Dim{:geotile}(geotiles.id)

    area_km2 = GGA._geotile_area_km2(; surface_mask, geotile_width)
    dpscale_new = Dim{:pscale}(0.25:0.25:4.0)
  
    if .!isnothing(single_geotile_test)
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)], OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!"
        geotiles = geotiles[geotiles.id .== single_geotile_test, :]
        geotiles2plot = [single_geotile_test]
    end
end;

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
# - geotile_buffer: Buffer distance in meters for spatial coverage
# - min_gemb_coverage: Minimum required GEMB data coverage (default: 0.75)
# - force_remake_before: DateTime threshold for forcing file regeneration
# - single_geotile_test: Optional single geotile ID for testing

begin 
    if .!isnothing(single_geotile_test)
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)], OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!"
    end
    
    gemb_files = vcat(GGA.allfiles.(gembinfo.gemb_folder; subfolders=false, fn_endswith=".mat", fn_contains=gembinfo.file_uniqueid)...)

    # ensure expected number of files found
    expected_number_of_files = length(gembinfo.elevation_delta) .* length(gembinfo.precipitation_scale) * 2; # multiply by 2 for NH and SH
    if length(gemb_files) != expected_number_of_files
        error("Expected $(expected_number_of_files) files but found $(length(gemb_files))")
    end

    gemb = GGA.read_gemb_files(gemb_files, gembinfo; vars2extract=vcat(dims2extract, vars2extract), date_range)

    @showprogress desc="Populating geotiles with GEMB data, this will take ~7 min on 128 threads [90GB of memory]" Threads.@threads for geotile_row in eachrow(geotiles)

        gemb_geotile_filename = replace(gembinfo.filename_gemb_combined, ".jld2" => "_$(geotile_row.id).jld2")

        if isfile(gemb_geotile_filename) && (isnothing(force_remake_before) || (Dates.unix2datetime(mtime(gemb_geotile_filename)) > force_remake_before)) && isnothing(single_geotile_test)
            printstyled("    -> Skipping $(gemb_geotile_filename) as it was created after the force_remake_before date: $force_remake_before \n"; color=:light_green)
        else

            # TODO: explore interpolating locally to fill gaps in the raw data instead of seacrhing for data in the buffer
            gemb0 = GGA.gemb2geotile(gemb, geotile_row; geotile_buffer, min_gemb_coverage, ddate)

            # plot after grouping into geotiles
            if !isnothing(geotiles2plot) && in(geotile_row.id, geotiles2plot)

                geotiles2plot0 = geotiles2plot[geotiles2plot .== geotile_row.id]
                figures = GGA.plot_area_average_height_gemb_ensemble(gemb0, area_km2; vars2plot=vars2extract, geotiles2plot=geotiles2plot0, title_prefix= "[1] binned raw data:")
                
                for geotile2plot in keys(figures)
                    for var0 in keys(figures[geotile2plot])
                        display(figures[geotile2plot][var0])
                    end
                end
            end

            # save here so that further processing does not need to carry `gemb` raw input in memory
            if isnothing(single_geotile_test)
                save(gemb_geotile_filename, gemb0)
            end
        end
    end
end

# =============================================================================
# GEMB DATA FILLING, ΔHEIGHT CLASS ADDITION, ADDITIONAL PSCALE CLASSES, AND VOLUME CHANGE
# =============================================================================
# Fills gaps in GEMB data and adds elevation change classes to the data.
# 
# Process:
# 1. Determine force remake threshold based on existing file timestamps
# 2. Initialize volume change data structure for all geotiles
# 3. Process each geotile in parallel:
#    - Load geotile-specific GEMB data
#    - Fill gaps and add elevation change classes
#    - Generate plots for visualization (if requested)
#    - Calculate volume changes with extended precipitation scaling
#    - Store results in combined data structure
# 4. Save combined volume change data to disk
#
# Performance: ~24 hours on 128 threads, consumes ~900GB memory
# Input: Individual geotile GEMB files
# Output: Combined volume change data with extended precipitation scaling classes
begin
    # Determine timestamp threshold for forcing file regeneration
    force_remake_before = maximum(Dates.unix2datetime.(mtime.(GGA.allfiles("/mnt/bylot-r3/data/gemb/raw/"; subfolders=false, fn_endswith = "].jld2"))))

    # Handle single geotile test mode
    if .!isnothing(single_geotile_test)
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)], OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!"
    end   
  
    # Initialize volume change data structure for all geotiles
    gemb_geotile_filename_dv = replace(gembinfo.filename_gemb_combined, ".jld2" => "_geotile_dv.jld2")

    # Check if output file exists and is newer than force remake threshold
    if isfile(gemb_geotile_filename_dv) && (isnothing(force_remake_before) || (Dates.unix2datetime(mtime(gemb_geotile_filename_dv)) > force_remake_before)) && .!isnothing(single_geotile_test)
        printstyled("    -> Skipping $(gemb_geotile_filename_dv) as it was created after the force_remake_before date: $force_remake_before \n"; color=:light_green)
    else

        # Initialize dictionary to store volume change data for all variables
        gemb_dv = Dict()
       
        # Pre-allocate arrays for each variable with NaN values
        for k in vcat(vars2extract, ["smb", "runoff"])
            gemb_dv[k] = fill(NaN, (dgeotile, ddate, dpscale_new, dΔheight))
        end

        #TODO: there is a memory leak here, it looks like the memory usage is not being freed after the loop... move more into funcitons

        #TODO: Add check-point saving of gemb_dv
        
        # Process all geotiles .. it looks like lower @threads is 7x faster than using @threads here
        @showprogress desc="filling gemb geotiles, adding Δheight classes, extending pscale classes, and saving volume change data to disk [~8 hrs on 128 threads] [consumes 250 GB of memory!!]... " for geotile_row in eachrow(geotiles) 
        
            # Load geotile-specific GEMB data from disk
            gemb_geotile_filename = replace(gembinfo.filename_gemb_combined, ".jld2" => "_$(geotile_row.id).jld2")
            gemb0 = load(gemb_geotile_filename) #[4s]
                
            # Fill data gaps and add elevation change classes (~27 seconds per geotile)
            gemb0 = GGA.gemb_fill_gaps_and_add_Δheight_classes(gemb0, dΔheight) #[27s]
    
            # Generate visualization plots after adding elevation change classes
            if !isnothing(geotiles2plot) && in(geotile_row.id, geotiles2plot)
                vars2plot = collect(keys(gemb0))
                geotiles2plot0 = geotiles2plot[geotiles2plot .== geotile_row.id]
                figures = GGA.plot_area_average_height_gemb_ensemble(gemb0, area_km2; vars2plot, geotiles2plot=geotiles2plot0, title_prefix= "[2] Δheight classes added:")
                
                # Display all generated figures
                for geotile2plot in keys(figures) 
                    for var0 in keys(figures[geotile2plot])
                        display(figures[geotile2plot][var0])
                    end
                end
            end
            
            # Calculate volume changes with extended precipitation scaling classes (~8 seconds per geotile)
            gemb_dv0 = GGA.gemb_dv(gemb0, area_km2[geotile=At(geotile_row.id)], dpscale_new) #[3s]

            # Store results in combined data structure for all geotiles
            for k in keys(gemb_dv) #[0.2s]
                gemb_dv[k][geotile = At(geotile_row.id)] = gemb_dv0[k]
            end

            # Generate visualization plots after adding precipitation scaling classes
            if !isnothing(geotiles2plot) && in(geotile_row.id, geotiles2plot)
                # Convert volume changes to height changes for plotting
                gemb_dh0 = deepcopy(gemb_dv0)
                for k in keys(gemb_dh0)
                    gemb_dh0[k] = gemb_dv0[k] ./ sum(area_km2[geotile=At(geotile_row.id)]) * 1000
                end
                figures = GGA.plot_dh_gemb_ensemble(gemb_dh0, area_km2; geotile=geotile_row.id, title_prefix= "[3] pscale classes added:")
                for var0 in keys(figures)
                    display(figures[var0])
                end
            end        
        end

        # Save combined volume change data to disk (skip if in test mode)
        if isnothing(single_geotile_test)
            println("Saving $(gemb_geotile_filename_dv)")
            save(gemb_geotile_filename_dv, gemb_dv)
        end
    end
end

mv("/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2024_820_40_racmo_grid_lwt_e97_0_geotile_dv_2.jld2", "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2024_820_40_racmo_grid_lwt_e97_0_geotile_dv.jld2")