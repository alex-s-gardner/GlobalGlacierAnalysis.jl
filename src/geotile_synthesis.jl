"""
    geotile_synthesis.jl

Main script for synthesizing glacier elevation change data across geotiles.

This script performs several key operations:
1. Calculates glacier hypsometry for different glacier inventories
2. Estimates global glacier discharge
3. Synthesizes elevation change data from multiple altimetry missions
4. Calibrates GEMB (Glacier Energy and Mass Balance) model to altimetry data
5. Computes volume and mass change time series for each geotile

The workflow handles multiple surface masks (glacier inventories), combines
measured and modeled discharge data, and applies calibration parameters to
GEMB model outputs to match altimetry observations.
"""

begin
    using Dates
    using ProgressMeter
    force_remake_before = Date(2025,6,17);

    # set force_remake == true to redo all steps from scratch
    force_remake_hypsometry = false; # these files are not altimetry dependent
    force_remake_discharge = false; # these files are not altimetry dependent
    force_remake_geotile_groups = false; # these files are not altimetry dependent
end

# LOAD PACKAGES,SET PATHS AND PARAMETERS
#begin
    import GlobalGlacierAnalysis as GGA
    using FileIO
    using DimensionalData
    using Statistics
    using Dates
    using LsqFit
    using Distances
    using ProgressMeter
    using DataFrames
    using Rasters
    using GeoDataFrames
    import GeometryOps as GO
    import GeoInterface as GI
    import GeoFormatTypes as GFT
    using BinStatistics
    using GeoFormatTypes
    using StatsBase
    using DataInterpolations
    using SortTileRecursiveTree
    using CairoMakie
    using GeometryOps
    using Arrow
    using Loess
    using DataInterpolations

    paths = GGA.pathlocal
    geotile_width = 2;
    downscale_to_glacier_method = "area"

    # to include in uncertainty
    project_id = ["v01"]
    surface_mask=["glacier", "glacier_rgi7"]
    dem_id=["best", "cop30_v2"]
    curvature_correct=[false, true]
    amplitude_correct=[true]
    binning_method = ["median","nmad5", "nmad3"]
    paramater_set=[1, 2, 3, 4]
    binned_folder=["/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"]

    filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0.jld2"
    filename_gemb_geotile = replace(filename_gemb_combined, ".jld2" => "_geotile.jld2")
    filename_gemb_geotile_filled_dv = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0_geotile_filled_extra_extrap_dv.jld2"

    geotile_groups_fn = joinpath(paths.data_dir, "project_data", "geotile_groups.arrow")

    globaldischarge_fn = joinpath(paths.data_dir, "GlacierOutlines/GlacierDischarge/global_glacier_discharge.jld2")

    # for latitudes below this set discharge2smb
    discharge2smb_max_latitude = -60;
    discharge2smb_equilibrium_period = (Date(1979), Date(2000))

    param_nt = (;project_id, surface_mask, dem_id, curvature_correct, amplitude_correct, binning_method, paramater_set, binned_folder)
    params = GGA.ntpermutations(param_nt)

    # only include files that exist
    path2runs = String[]
    for param in params
    #param = params[1]
        binned_aligned_file, _ = GGA.binned_filled_filepath(; param...)
        if isfile(binned_aligned_file)
            push!(path2runs, binned_aligned_file)
        end
    end

    ## calculate individual glacier hypsometry
    geomfile = Dict()
    for sm in surface_mask
        if sm == "glacier"
            geomfile[sm] = joinpath(paths.data_dir, "GlacierOutlines/rgi60/rgi60_Global.gpkg")
        elseif sm == "glacier_rgi7"
            geomfile[sm] = joinpath(paths.data_dir, "GlacierOutlines/RGI2000-v7.0-G-global-fix/rgi70_Global.gpkg")
        end
    end

    glacier_geotile_hyps_fn = Dict()
    for sm in surface_mask
        glacier_geotile_hyps_fn[sm] = replace(geomfile[sm], ".gpkg" => "geotile_hyps.jld2")
    end

    #geomfile_rgi6 = joinpath(paths.data_dir, "GlacierOutlines/rgi60/rgi60_Global.gpkg")
    #geomfile_rgi7 = joinpath(paths.data_dir, "GlacierOutlines/RGI2000-v7.0-G-global-fix/rgi70_Global.gpkg")

    # Manual override for specific geotile groups
    # This handles cases where large glaciers cross multiple tiles but should be treated separately
    # NOTE: groupings are only updated if force_remake == true
    geotile_groups_manual = [
        ["lat[+62+64]lon[-148-146]"],
        ["lat[+62+64]lon[-146-144]"],
        ["lat[+60+62]lon[-138-136]"],
        ["lat[+58+60]lon[-138-136]"],
        ["lat[+58+60]lon[-136-134]"],
        ["lat[+58+60]lon[-134-132]"],
        ["lat[+56+58]lon[-134-132]"],
        ["lat[+78+80]lon[-084-082]"],
        ["lat[-52-50]lon[-074-072]"],
        ["lat[+56+58]lon[-130-128]"],
        ["lat[-74-72]lon[-080-078]", "lat[-74-72]lon[-078-076]"],
        ["lat[+78+80]lon[+010+012]", "lat[+78+80]lon[+012+014]", "lat[+78+80]lon[+014+016]"],
        ["lat[+76+78]lon[+062+064]", "lat[+76+78]lon[+064+066]", "lat[+76+78]lon[+066+068]", "lat[+76+78]lon[+068+070]", "lat[+74+76]lon[+062+064]", "lat[+74+76]lon[+064+066]", "lat[+74+76]lon[+066+068]", "lat[+74+76]lon[+066+068]"],
    ]

end


"""
    Calculate glacier hypsometry for each surface mask inventory

Processes each glacier inventory to calculate area-elevation distributions:
1. Determines appropriate glacier ID field based on RGI version
2. Reads glacier outlines from the specified inventory
3. Loads COP30 DEM and defines elevation bins
4. Identifies geotiles containing glaciers
5. Calculates hypsometric area distribution for each glacier-geotile combination
6. Saves results to JLD2 file for later use

Note: This calculation is computationally intensive (~37 minutes per inventory)
"""

for sm in surface_mask 
    begin
        # Set the glacier ID attribute name based on RGI version
        persistent_attribute = if Base.contains(geomfile[sm], "rgi6")
            :RGIId
        else
            :rgi_id
        end

        # Only process if output file doesn't exist or force_remake is true
        if !isfile(glacier_geotile_hyps_fn[sm]) || force_remake_hypsometry

            # Read glacier outlines
            glacier_geom = GeoDataFrames.read(geomfile[sm])
            
            # Load elevation data and define height bins
            h = Raster(paths.cop30_v2, lazy=true)
            height_range, height_center = GGA.project_height_bins()
            
            # Get geotiles that contain glaciers
            geotiles = GGA.geotiles_w_mask(geotile_width)
            geotiles = geotiles[geotiles.glacier_frac .> 0.0, :]

            
            # NOTE: glaciers in this geotile have elevations < 0 in COP30 DEM]
            # for testing
            # gtid = "lat[-76-74]lon[-102-100]"
            # geotiles = geotiles[geotiles.id .== gtid, :]

            # Calculate hypsometry for each glacier-geotile combination
            # This is computationally intensive (~37 min)
            glacier_geom = geotile_zonal_area_hyps(h, height_range, glacier_geom, geotiles.id; persistent_attribute)
            
            # Save results
            FileIO.save(glacier_geotile_hyps_fn[sm], Dict("glaciers" => glacier_geom))
        end
    end 
end

"""
Process and calculate global glacier discharge [aka frontal ablation] data.

This section:
1. Loads glacier hypsometry data from previously generated files
2. Downscales geotile-level GEMB model data to individual glaciers using either hypsometry or area-based methods
3. Estimates discharge for unmeasured glaciers based on surface mass balance
4. Combines estimated discharge with measured discharge data
5. Saves the combined dataset to a file

Parameters:
- Uses global variables for file paths and configuration settings
- Supports "hypsometry" or "area" downscaling methods
- Applies discharge estimation only for glaciers with non-zero area
"""

if .!isfile(globaldischarge_fn) || force_remake_discharge # Load RGI6 glacier hypsometry data from previously generated file glacier_geotile_hyps_fn = joinpath(paths.data_dir, "GlacierOutlines/rgi60/rgi60_Global.gpkg")
    sm = "glacier"
    
    glaciers0 = load(glacier_geotile_hyps_fn[sm], "glaciers")

    # Get geotiles containing glaciers and initialize parameters:
    # - Filter to only tiles with glaciers
    # - Set precipitation scaling to 1.0 (no scaling)
    # - Set height adjustment to 0
    geotiles = GGA.geotiles_w_mask(geotile_width)
    geotiles = geotiles[(geotiles.glacier_frac.>0.0), :]

    pscale = 1;
    Δheight = 0;

    # load example geotile to get dimensions
    if downscale_to_glacier_method == "hypsometry"
        filename_gemb_geotile_filled = replace(filename_gemb_geotile, ".jld2" => "_$(glaciers0.geotile[1])_filled.jld2")
        smb = load(filename_gemb_geotile_filled, "smb")

        dgeotile = Dim{:geotile}(geotiles.id)
        ddate = dims(smb, :date)
        dheight = dims(smb, :height)

        var1 = fill(NaN, dgeotile, ddate, dheight)

        @time for varname = ["smb"]
            @showprogress desc="Populate $(varname)..." Threads.@threads for geotile in geotiles.id
                filename_gemb_geotile_filled = replace(filename_gemb_geotile, ".jld2" => "_$(geotile)_filled.jld2")
                var0 = load(filename_gemb_geotile_filled, varname)
                var1[At(geotile), :, :] = var0[:,:,At(pscale), At(Δheight)]
            end

            # Convert geotile-level data to per-glacier values
            glaciers0 = GGA.geotile2glacier!(glaciers0, var1; varname)
        end
    elseif downscale_to_glacier_method == "area"
        glacier_area_km2 = sum.(glaciers0.area_km2)
        smb = load(filename_gemb_geotile_filled_dv, "smb")
        smb =   smb[:,:,At(pscale), At(Δheight)]
            
        ddate = dims(smb, :date)
        glaciers0[!,"smb"] =  [fill(NaN, ddate) for _ in 1:nrow(glaciers0)]

        for geotile in unique(glaciers0.geotile)
        #geotile = "lat[+28+30]lon[+082+084]"
            gindex = findall(geotile .== glaciers0.geotile)
            gt_dv = smb[At(geotile), :]
            garea = sum.(glaciers[gindex, :area_km2])
            gweighting = garea ./ sum(garea)
            for (i,ig) in enumerate(gindex)
                glaciers[ig,"smb"][:] = gt_dv * gweighting[i] ./ garea[i] * 1000 # convert from km3 to mwe
            end
        end
    else
        error("downscale_to_glacier_method must be either \"hypsometry\" or \"area\"")
    end

    # For unmeasured glaciers, estimate discharge using surface mass balance
    # Only process glaciers with non-zero area
    discharge0 = GGA.discharge2smb(
        glaciers0[sum.(glaciers0.area_km2).>0, :];
        discharge2smb_max_latitude,
        discharge2smb_equilibrium_period
    )

    # Combine estimated discharge with measured discharge data
    discharge = GGA.glacier_discharge(; datadir=GGA.pathlocal[:data_dir])
    discharge = vcat(discharge, discharge0)
    
    # Set any negative discharge values to zero
    discharge[discharge.discharge_gtyr.<0, :discharge_gtyr] .= 0

    # Save the combined measured and estimated discharge data
    save(globaldischarge_fn, Dict("discharge"=>discharge))
else
    # If file exists, just load the discharge data
    discharge = load(globaldischarge_fn, "discharge")
end



path2geotile_synthesis_error = GGA.geotile_synthesis_error(; 
    path2runs, 
    outfile="/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2", 
    force_remake_before,
)


"""
    geotile_synthesize_runs(;
        path2runs,
        path2geotile_synthesis_error,
        missions2include=["hugonnet" "gedi" "icesat" "icesat2"],
        force_remake_before
    )

Synthesize elevation change data from multiple altimetry missions into a unified dataset.

# Arguments
- `path2runs`: Paths to the aligned altimetry data files
- `path2geotile_synthesis_error`: Path to the geotile synthesis error file
- `missions2include`: Array of mission names to include in the synthesis
- `force_remake_before`: Date threshold for forcing remake of existing files
"""

GGA.geotile_synthesize_runs(;
    path2runs,
    path2geotile_synthesis_error, 
    missions2include=["hugonnet" "gedi" "icesat" "icesat2"],
    force_remake_before
)


"""
Load and prepare geotile data for multiple surface masks.

This function:
1. Updates file paths to point to synthesized data
2. Loads elevation change data to extract geotile dimensions
3. Processes geotiles for each surface mask by:
   - Generating geotiles with hypsometry information
   - Making geotiles mutually exclusive by RGI region
   - Aligning geotile indices with elevation change dimensions
   - Standardizing area column names
4. Returns dictionaries of geotiles and glaciers for each surface mask

Returns dictionaries of geotiles and glaciers indexed by surface mask.
"""

begin
    # Update paths to point to synthesized data
    path2runs = replace.(path2runs, "aligned.jld2" => "synthesized.jld2")

    # Load elevation change data to get geotile dimensions
    dh = FileIO.load(path2runs[1], "dh_hyps")

    # Get geotile dimensions
    dgeotile = dims(dh, :geotile)

    # Process geotiles for each surface mask
    geotiles = Dict()
    glaciers = Dict()
    
    for sm in surface_mask
    #sm = surface_mask[1]
        # Generate geotiles for this mask
        geotiles0 = GGA.geotiles_mask_hyps(sm, geotile_width)
        
        # Make geotiles mutually exclusive by RGI region
        geotiles0, reg = GGA.geotiles_mutually_exclusive_rgi!(copy(geotiles0))
        
        # Align geotile indices with elevation change dimensions
        gt_ind = [findfirst(geotiles0.id .== g0) for g0 in collect(dgeotile)]
        
        # Rename area column to include surface mask
        rename!(geotiles0,"$(sm)_area_km2" => "area_km2")
        
        # Store in dictionary
        geotiles[sm] = geotiles0[gt_ind, :]

        glaciers[sm] = FileIO.load(glacier_geotile_hyps_fn[sm], "glaciers")
    end
end


"""
    Group geotiles based on glacier overlap

This code identifies and groups geotiles that are connected by large glaciers
that cross tile boundaries. This ensures that glaciers spanning multiple tiles
are treated consistently in subsequent analysis.

The function:
1. Uses a minimum glacier area threshold (default 100 km²) to identify large glaciers
2. Identifies which geotiles these large glaciers overlap with
3. Creates groups of connected geotiles based on shared glaciers
4. Applies manual grouping overrides for specific regions
5. Saves the grouping information to a file for later use

The grouping is only recalculated if the output file doesn't exist or if forced.
"""

if !isfile(geotile_groups_fn) || force_remake_geotile_groups
#@time begin
    # Filter for large glaciers (>100 km2) that may cross tile boundaries
    min_area_km2 = 100
    sm = surface_mask[1]

    # Find intersection between geotiles and glaciers
    geotiles0 = deepcopy(geotiles[sm])
    
    # Load glacier outlines
    glaciers0 = GeoDataFrames.read(geomfile[sm])

    # identify groups of geotiles that are connected by large overlapping glaciers
    geotiles0 = geotile_grouping!(geotiles0, glaciers0, min_area_km2; geotile_groups_manual)
    
    # Write geotile groups to file, excluding some columns
    GeoDataFrames.write(geotile_groups_fn, geotiles0[:,[:geometry, :id, :group]]; crs=GFT.EPSG(4326))
end

"""
    Calibrate GEMB model to altimetry data

This function calibrates the Glacier Energy and Mass Balance (GEMB) model to match altimetry observations:
1. Loads GEMB surface mass balance (SMB) and firn air content (FAC) data
2. Assigns geotile groups based on glacier connectivity
3. For each altimetry dataset:
   - Converts elevation change to volume change
   - Finds optimal GEMB parameters by fitting to altimetry data
   - Accounts for glacier groups to handle glaciers spanning multiple tiles
   - Calculates mean absolute deviation (MAD) as a quality metric
   - Saves calibration parameters to Arrow files for later use

The calibration is only performed if output files don't exist or if forced by date.
"""

begin
    # load gemb data
    # Note: GEMB volume change is for a single surface mask... this is a hack has making gemb dv for multiple surface masks is onerious... that said, because GEMB data is calibrated to altimetry that uses different surface masks.. this is mostly a non-issue

    # NOTE: Variables are in units of mie
    smb, fac = load(filename_gemb_geotile_filled_dv, ("smb", "fac"))
 
    # add in groups
    sm = surface_mask[1]
    geotiles0 = deepcopy(geotiles[sm])
    geotiles0[!,:group] .= 0
    geotiles_groups = GeoDataFrames.read(geotile_groups_fn)
    for row in eachrow(geotiles0)
        row.group = geotiles_groups[findfirst(isequal(row.id), geotiles_groups.id), :group]
    end

    @showprogress desc="Finding optimal GEMB fits to altimetry data..." Threads.@threads for binned_synthesized_file in path2runs

        synthesized_gemb_fit = replace(binned_synthesized_file,  ".jld2" => "_gemb_fit.arrow")

        if !(isfile(synthesized_gemb_fit)) || !isnothing(force_remake_before)

            if isfile(synthesized_gemb_fit) && !isnothing(force_remake_before)
                if Dates.unix2datetime(mtime(synthesized_gemb_fit)) > force_remake_before
                    @warn "Skipping $(synthesized_gemb_fit) because it was created after $force_remake_before"
                    continue
                end
            end
            
            t1 = time()

            run_parameters = GGA.binned_filled_fileparts(synthesized_gemb_fit)
            sm = run_parameters.surface_mask;

            dh = FileIO.load(binned_synthesized_file, "dh_hyps")
    
            # convert elevation change to volume change
            dv_altim = GGA.dh2dv(dh, geotiles[sm]);
  
            if any(all(isnan.(dv_altim), dims = :date))
                # keep this error message
                dgeotile = dims(dv_altim, :geotile)
                all_nan_geotiles = dgeotile[findall(dropdims(all(isnan.(dv_altim), dims = :date), dims = :date))]
                printstyled("Geotiles that only contain NaNs:\n"; color = :red)
                println("$(collect(all_nan_geotiles)))")
                printstyled("Possible sources of error include:\n"; color = :red)
                printstyled("  [1] synthesis_error_file was run on a subset of files that do not include the full error range of the data, to fix this you need to delete the exisiting error file and any synthesis files that were created :\n [2] the wrong surface mask was passed to GGA.dh2dv\n"; color = :red)
                error("geotile volume change contains all NaNs")
            end


            # find optimal fit to gemb data
            # there are issues with calibrating the smb model to individual geotiles glaciers 
            #can cross multiple geotiles therefore we calibrate the model for groups of 
            # distinct geotiles.

            # using geotiles0 instead of geotiles[sm] is fine here
            examine_model_fits = nothing; #"lat[+76+78]lon[-094-092]"
            #examine_model_fits = "lat[+30+32]lon[+096+098]"
            #examine_model_fits = "lat[-50-48]lon[-074-072]"
            #examine_model_fits = "lat[+58+60]lon[-132-130]"
            #examine_model_fits = "lat[+52+54]lon[-128-126]"
        
            df = GGA.gemb_bestfit_grouped(dv_altim, smb, fac, discharge, geotiles0; examine_model_fits)

            df[!,:area_km2] = sum.(geotiles[sm].area_km2)
            df.mad = df.mad ./ df.area_km2
            rename!(df, "mad"=>"mad_m")
            df.geometry = GGA.extent2rectangle.(df.extent)
            df = df[:,Not(:extent)]
            GeoDataFrames.write(synthesized_gemb_fit, df)
            println("\n $binned_synthesized_file optimal GEMB fit found: $(round(Int,time() -t1))s")
        end
    end
end


"""
Generate calibrated geotile-level timeseries by combining altimetry and GEMB model data.

This function:
1. Loads synthesized altimetry data and GEMB fit parameters
2. Processes each geotile containing glaciers
3. Applies calibration parameters from GEMB model fits
4. Converts elevation changes to volume and mass changes
5. Handles discharge data and firn air content corrections
6. Aggregates data for glacier groups that span multiple geotiles
7. Computes regional estimates by RGI region

The function creates comprehensive timeseries including:
- Altimetry-derived volume changes (dv_altim)
- Surface mass balance components (runoff, smb, rain, etc.)
- Discharge estimates
- Combined model-based volume/mass changes (dv, dm)
- Altimetry-derived mass changes with firn correction (dm_altim)

Output is saved as JLD2 files with "_gembfit_dv.jld2" suffix.
"""
begin
    # Load synthesized data and GEMB fit parameters
    discharge = load(globaldischarge_fn, "discharge")
    volume2mass = GGA.δice / 1000

    has_glacier = Dict()
    for keys in keys(geotiles)
        has_glacier[keys] = [id in unique(glaciers[keys].geotile) for id in geotiles[keys].id]
    end

    # removed Threads.@threads due to memory issues
    @showprogress desc="Computing calibrated geotile level data timeseries for all runs..." for binned_synthesized_file in path2runs 

        binned_synthesized_dv_file = replace(binned_synthesized_file, ".jld2" => "_gembfit_dv.jld2")

        run_parameters = GGA.binned_filled_fileparts(binned_synthesized_file)
        sm = run_parameters.surface_mask

        glaciers0 = copy(glaciers[sm])

        if !isfile(binned_synthesized_dv_file) || !isnothing(force_remake_before)

            if isfile(binned_synthesized_dv_file) && !isnothing(force_remake_before)
                if Dates.unix2datetime(mtime(binned_synthesized_dv_file)) > force_remake_before
                    @warn "Skipping $(binned_synthesized_dv_file) because it was created after $force_remake_before"
                    continue
                end
            end

            synthesized_gemb_fit = replace(binned_synthesized_file, ".jld2" => "_gemb_fit.arrow")
        
            gemb_fit = GeoDataFrames.read(synthesized_gemb_fit);
        
            geotiles0 = copy(geotiles[sm])
            geotile_groups = GeoDataFrames.read(geotile_groups_fn)
        
            # Filter to only include geotiles containing glaciers
            geotiles0 = geotiles0[has_glacier[sm],:]

            # Define variables to process
            varnames = ["dv_altim", "runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze"]

            # Add group assignments from geotile_groups
            geotiles0[!, :group] .= 0
            for geotile in eachrow(geotiles0)
                geotile.group = geotile_groups[findfirst(isequal(geotile.id), geotile_groups.id), :group]
            end

            # Add GEMB fit parameters and mass conversion factors
            geotiles0[!, :pscale] .= 0.0
            geotiles0[!, :Δheight] .= 0.0
            geotiles0[!, :mie2cubickm] .= 0.0

            # load to get date vector
            var0 = load(filename_gemb_geotile_filled_dv, "smb")
            ddate = dims(var0, :date)
            Δdecyear = GGA.decimalyear.(ddate) .- GGA.decimalyear.(ddate[1])
            geotiles0[!, :discharge] .=[fill(NaN, length(ddate)) for _ in 1:nrow(geotiles0)]
            
            for geotile in eachrow(geotiles0)
                    fit_index = findfirst(isequal(geotile.id), gemb_fit.id)
                    geotile.pscale = gemb_fit[fit_index, :pscale]
                    geotile.Δheight = gemb_fit[fit_index, :Δheight]
                    area_km2 = sum(geotile.area_km2)
                    geotile.mie2cubickm = area_km2/1000  # Convert meters ice equivalent to cubic kilometers

                    # include discharge average over total glacier area [i.e. save in units of mie]
                    index = GGA.within.(Ref(geotile.extent), discharge.longitude, discharge.latitude)
                    geotile.discharge = sum(discharge.discharge_gtyr[index])/volume2mass * Δdecyear # Gt to units of mie
            end

            # add discharge date metadata
            colmetadata!(geotiles0, "discharge", "date", collect(ddate), style=:note)
            colmetadata!(geotiles0, "discharge", "units", "km3 [i.e.]", style=:note)

            # Process each variable
            for varname in varnames

                # Special handling for height change data
                if varname == "dv_altim"
                    dh = load(binned_synthesized_file, "dh_hyps")
                    var0 = GGA.dh2dv(dh, geotiles0);
                else
                    var0 = load(filename_gemb_geotile_filled_dv, varname)
                end
                
                # Initialize time series arrays
                ddate = dims(var0, :date)
                geotiles0[!, varname] = [fill(NaN, length(ddate)) for _ in 1:nrow(geotiles0)]

                # Add date metadata for time series
                colmetadata!(geotiles0, varname, "date",collect(ddate), style=:note)
                colmetadata!(geotiles0, varname, "units", "km3 [i.e.]", style=:note)
                
                # Apply GEMB scaling and convert units for each geotile
                for geotile in eachrow(geotiles0)
                    geotile[varname] = var0[At(geotile.id), :, At(geotile.pscale), At(geotile.Δheight)];
                end
            end

            # modeled dv
            geotiles0[!, :dv] = (geotiles0[:, :smb] .- geotiles0[:, :discharge] .+ geotiles0[:, :fac])
            colmetadata!(geotiles0, "dv", "date", collect(ddate), style=:note) # to units of km3 assuming an ice density of 910 kg/m3
            colmetadata!(geotiles0, "dv", "units", "km3 [i.e.]", style=:note)

            # modeled dm
            geotiles0[!, :dm] = (geotiles0[:, :smb] .- geotiles0[:, :discharge]) .* volume2mass # to units of Gt
            colmetadata!(geotiles0, "dm", "date", collect(ddate), style=:note)
            colmetadata!(geotiles0, "dm", "units", "Gt", style=:note)

            # altimetry dm
            geotiles0[!, :dm_altim] = copy(geotiles0[!, :dv_altim])
            for r in eachrow(geotiles0)
                model = LinearInterpolation(r.fac, GGA.decimalyear.(colmetadata(geotiles0, "fac", "date")))
                fac = model(GGA.decimalyear.(colmetadata(geotiles0, "dv_altim", "date")))
                r.dm_altim = (r.dv_altim .- fac) .* volume2mass # in units of Gt
            end

            colmetadata!(geotiles0, "dm_altim", "date", colmetadata(geotiles0, "dv_altim", "date"), style=:note)
            colmetadata!(geotiles0, "dm_altim", "units", "Gt", style=:note)

            varnames_all = vcat(varnames, "discharge", "dv", "dm", "dm_altim")
            
            # `dichage`, `dm` and `dv` must be average by geotile groups as they are not valid for single geotiles
            vars2average = ["discharge", "dv", "dm", "dv_altim", "dm_altim"]
            gdf = groupby(geotiles0, :group)

            for g in gdf
                #g = gdf[4]
                if nrow(g) > 1
                    for varn in vars2average
                    #varn = "dv_altim"
                        foo = zeros(eltype(g[1,varn]),length(g[1,varn]))
                        for r in eachrow(g)
                            foo .+= r[varn]
                        end
                        foo = foo ./ sum(g[:, :mie2cubickm])

                        for r in eachrow(g)
                            r[varn] = foo .* r[:mie2cubickm]
                        end
                    end
                end
            end

            #### now compute regional estimates
            geotiles0[!, :rgiid] .= Int8(0)
            for i = 1:19
                rgi = "rgi$i"
                geotiles0[geotiles0[:, rgi] .> 0, :rgiid] .= i
                geotiles0 = geotiles0[:, Not(rgi)]
            end

            # sanity check
            for varname in varnames_all
                if any([all(isnan.(v)) for v in geotiles0[:,varname]])
                    println("$(varname) has all NaNs")
                end
            end

            FileIO.save(binned_synthesized_dv_file, Dict("geotiles" => geotiles0))
        end
    end
end