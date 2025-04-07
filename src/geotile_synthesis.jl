# This script performs synthesis and analysis of glacier elevation change data.
# Key processing steps:
#
# 1. Initialization and Setup (~1 min):
#    - Imports required packages for data processing and analysis
#    - Sets up paths and parameters for geotile synthesis
#    - Configures project ID, surface masks, DEM sources
#    - Builds parameter combinations for uncertainty analysis
#
# 2. Glacier Hypsometry Processing (~37 min per inventory):
#    - Processes both RGI6 and RGI7 glacier inventories
#    - Calculates area-elevation distributions for each glacier
#    - Saves hypsometry data to JLD2 files
#
# 3. Discharge Data Processing:
#    - Handles measured discharge data
#    - Estimates discharge for unmeasured areas using SMB
#    - Combines measured and estimated discharge
#
# 4. Synthesis Error Calculation (~2 hrs):
#    - Calculates uncertainties for combining multiple altimetry missions
#    - Saves error estimates for later weighting
#
# 5. Multi-mission Data Synthesis (~2 hrs):
#    - Combines data from multiple altimetry missions
#    - Uses synthesis errors to weight different missions
#    - Processes 252 run files
#
# 6. Geotile to Glacier Downscaling:
#    - Converts geotile-level data to individual glaciers
#    - Uses either hypsometry or area-based methods
#    - Handles GEMB variables (SMB, FAC, runoff)
#
# 7. Time Series Analysis and Export:
#    - Fits models to glacier time series (2000-2023)
#    - Calculates glacier centroids
#    - Exports rates to geospatial formats
#    - Creates gridded SMB averages
#
# Total runtime is approximately 5-6 hours, with most time spent on
# synthesis error calculation and multi-mission synthesis steps.


using Dates
using ProgressMeter
force_remake_before = Date(2025,3,21);

# set force_remake == true to redo all steps from scratch
force_remake_hypsometry = false; # these files are not altimetry dependent
force_remake_discharge = false; # these files are not altimetry dependent
force_remake_geotile_groups = false; # these files are not altimetry dependent

# LOAD PACKAGES,SET PATHS AND PARAMETERS
begin
    using Altim
    using FileIO
    using DimensionalData
    using Statistics
    using Dates
    using LsqFit
    using Distances
    using ProgressMeter
    using Plots
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
    include("utilities_synthesis.jl")

    # 1. Loads local configuration paths
    # 2. Sets key parameters:
    #    - Project ID, surface masks, DEM sources, and processing options for uncertainty analysis
    #    - Latitude threshold for discharge-to-SMB conversion
    #    - Geotile width of 2 degrees
    #    - Time period for equilibrium calculations
    # 3. Builds parameter combinations and filters to only include existing files
    # 4. Sets up paths for glacier geometry files (RGI v6 and v7)

    paths = Altim.pathlocal
    geotile_width = 2;
    downscale_to_glacier_method = "area"
    reference_run = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_dh_best_cc_meanmadnorm3_v01_filled_ac_p1_synthesized.jld2"
    path2runs_override = nothing #[reference_run]

    # to include in uncertainty
    project_id = ["v01"]
    surface_mask=["glacier", "glacier_rgi7"]
    dem_id=["best", "cop30_v2"]
    curvature_correct=[false, true]
    amplitude_correct=[true]
    binning_method = ["median","meanmadnorm5", "meanmadnorm3"]
    paramater_set=[1, 2, 3, 4]
    binned_folder=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg")

    filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0.jld2"
    filename_gemb_geotile = replace(filename_gemb_combined, ".jld2" => "_geotile.jld2")
    filename_gemb_geotile_filled_dv = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0_geotile_filled_extra_extrap_dv.jld2"

    geotile_groups_fn = joinpath(paths.data_dir, "project_data", "geotile_groups.arrow")

    globaldischarge_fn = joinpath(paths.data_dir, "GlacierOutlines/GlacierDischarge/global_glacier_discharge.jld2")

    # for latitudes below this set discharge2smb
    discharge2smb_max_latitude = -60;
    discharge2smb_equilibrium_period = (Date(1979), Date(2000))

    param_nt = (;project_id, surface_mask, dem_id, curvature_correct, amplitude_correct, binning_method, paramater_set, binned_folder)
    params = Altim.ntpermutations(param_nt)

    # only include files that exist
    path2runs = String[]
    for param in params
        binned_aligned_file = Altim.binned_aligned_filepath(; param...)
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

# GLACIER HYPSOMETRY BY GEOTILE
# For each inventory:
# 1. Generates output filename for storing hypsometry data
# 2. Determines glacier ID field name based on RGI version (RGIId for v6, rgi_id for v7)
# 3. If output file doesn't exist:
#    - Reads glacier outlines from GeoPackage file
#    - Loads Copernicus DEM elevation data 
#    - Defines elevation bins for hypsometry calculation
#    - Gets geotiles (2-degree grid cells) containing glaciers
#    - Calculates area-elevation distribution (hypsometry) for each glacier-geotile combination
#    - Saves results to JLD2 file
# Note: Hypsometry calculation takes ~37 minutes per inventory
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
            height_range, height_center = Altim.project_height_bins()
            
            # Get geotiles that contain glaciers
            geotiles = Altim.geotiles_w_mask(geotile_width)
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


# GLACIER DISCHARGE
# This section processes glacier discharge data by:
# 1. Loading RGI6 glacier hypsometry data
# 2. Getting geotiles containing glaciers and initializing parameters:
#    - Filters to tiles with glaciers
#    - Sets precipitation scaling to 1.0 (no scaling) 
#    - Sets height adjustment to 0
# 3. Loading and processing SMB data:
#    - Loads SMB data for each geotile
#    - Converts geotile-level data to per-glacier values
# 4. Handling discharge data:
#    - Estimates discharge for unmeasured glaciers using SMB
#    - Combines with measured discharge data
#    - Sets negative discharge values to zero
# 5. Saves the final combined discharge data to file

if .!isfile(globaldischarge_fn) || force_remake_discharge # Load RGI6 glacier hypsometry data from previously generated file glacier_geotile_hyps_fn = joinpath(paths.data_dir, "GlacierOutlines/rgi60/rgi60_Global.gpkg")
    sm = "glacier"
    
    glaciers0 = load(glacier_geotile_hyps_fn[sm], "glaciers")

    # Get geotiles containing glaciers and initialize parameters:
    # - Filter to only tiles with glaciers
    # - Set precipitation scaling to 1.0 (no scaling)
    # - Set height adjustment to 0
    geotiles = Altim.geotiles_w_mask(geotile_width)
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
            glaciers0 = Altim.geotile2glacier!(glaciers0, var1; varname)
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
    discharge0 = Altim.discharge2smb(
        glaciers0[sum.(glaciers0.area_km2).>0, :];
        discharge2smb_max_latitude,
        discharge2smb_equilibrium_period
    )

    # Combine estimated discharge with measured discharge data
    discharge = Altim.glacier_discharge(; datadir=Altim.pathlocal[:data_dir])
    discharge = vcat(discharge, discharge0)
    
    # Set any negative discharge values to zero
    discharge[discharge.discharge_gtyr.<0, :discharge_gtyr] .= 0

    # Save the combined measured and estimated discharge data
    save(globaldischarge_fn, Dict("discharge"=>discharge))
else
    # If file exists, just load the discharge data
    discharge = load(globaldischarge_fn, "discharge")
end


# MISSION ERROR AS SPREAD OF ALL MODEL RUNS
# 1. Calls geotile_synthesis_error() function which:
#    - Takes path to run files as input
#    - Processes ~252 run files (takes ~2 hours)
#    - Calculates errors/uncertainties for synthesizing different missions
#    - Saves results to a JLD2 file at specified output path
#    - Has option to force remaking the error file
#
# 2. The synthesis errors will be used later to:
#    - Weight different missions when combining their measurements
#    - Account for varying uncertainties between missions
#    - Optimize the final synthesis of elevation changes
@time path2geotile_synthesis_error = geotile_synthesis_error(; 
    path2runs, 
    outfile="/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2", 
    force_remake_before,
)

# SYNTHESIS OF ALL MISSIONS
# 1. Takes ~13 min to process 252 run files
# 2. Combines data from specified missions:
#    - Hugonnet glacier mass balance
#    - GEDI laser altimetry 
#    - ICESat laser altimetry
#    - ICESat-2 laser altimetry
# 3. Uses synthesis error calculations from previous step to:
#    - Weight contributions from each mission
#    - Account for varying uncertainties
# 4. Can force reprocessing of synthesis with force_remake_before
# rm.(allfiles("/mnt/bylot-r3/data/binned/2deg/"; fn_endswith = "synthesized.jld2"))
geotile_synthesize_runs(;
    path2runs,
    path2geotile_synthesis_error, missions2include=["hugonnet" "gedi" "icesat" "icesat2"],
    force_remake_before
)

# Load geotile info for mutiple surface masks
begin
    # Update paths to point to synthesized data
    if isnothing(path2runs_override)
        path2runs = replace.(path2runs, "aligned.jld2" => "synthesized.jld2")
    else
        path2runs = path2runs_override
    end
    
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
        geotiles0 = Altim.geotiles_mask_hyps(sm, geotile_width)
        
        # Make geotiles mutually exclusive by RGI region
        geotiles0, reg = Altim.geotiles_mutually_exclusive_rgi!(copy(geotiles0))
        
        # Align geotile indices with elevation change dimensions
        gt_ind = [findfirst(geotiles0.id .== g0) for g0 in collect(dgeotile)]
        
        # Rename area column to include surface mask
        rename!(geotiles0,"$(sm)_area_km2" => "area_km2")
        
        # Store in dictionary
        geotiles[sm] = geotiles0[gt_ind, :]

        glaciers[sm] = FileIO.load(glacier_geotile_hyps_fn[sm], "glaciers")
    end
end

# IDENTIFY GEOTILE GROUPINGS FOR MODEL CALLIBRATION
# 1. Initializes columns in geotiles dataframe to store:
#    - Indices of intersecting glaciers (discharge_ind)
#    - Indices of connected geotiles (geotile_intersect)
#
# 2. Finds intersections between large glaciers (>100 km2) and geotiles:
#    - Builds spatial index tree for efficient intersection testing
#    - For each geotile, finds intersecting glaciers using tree queries
#
# 3. Groups connected geotiles:
#    - Identifies geotiles that share glaciers
#    - Creates groups of connected geotiles
#    - Applies manual overrides for specific cases where large glaciers 
#      cross tiles but should be treated separately
#
# 4. Saves results:
#    - Converts geotile extents to rectangles
#    - Writes geotile groups to file, excluding temporary columns
# This takes ~10 seconds

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


# FIND OPTIMAL GEMB MODEL FIT TO ALTIMETRY DATA
# 1. Loads GEMB data:
#    - Reads SMB (Surface Mass Balance) and FAC (Firn Air Content) from filled geotile data
#
# 2. For each binned/synthesized altimetry file:
#    - Converts elevation change (dh) to volume change (dv) for each geotile
#    - Performs sanity checks to ensure no geotiles contain only NaN values
#
#    - Handles troubleshooting for specific geotiles if needed
# 3. Fits GEMB model to altimetry data:
#    - Calibrates model at regional rather than individual geotile scale
#    - This approach handles cases where glacier discharge crosses multiple geotiles
#    - Future improvements could include:
#      * Identifying geotile groups that enclose glacier groups
#      * Calculating flux across geotiles
#
# 4. Saves results:
#    - Stores fit parameters and statistics for each geotile group
#    - Includes area normalization of mean absolute deviation (MAD)

# Takes ~13 seconds
begin
    # load gemb data
    # Note: GEMB volume change is for a single surface mask... this is a hack has making gemb dv for multiple surface masks is onerious... that said, because GEMB data is calibrated to altimetry that uses different surface masks.. this is mostly a non-issue
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
    #for binned_synthesized_file in [reference_run]
    #binned_synthesized_file = reference_run
        synthesized_gemb_fit = replace(binned_synthesized_file,  ".jld2" => "_gemb_fit.arrow")


        if !(isfile(synthesized_gemb_fit)) || !isnothing(force_remake_before)

            if isfile(synthesized_gemb_fit) && !isnothing(force_remake_before)
                if Dates.unix2datetime(mtime(synthesized_gemb_fit)) > force_remake_before
                    @warn "Skipping $(synthesized_gemb_fit) because it was created after $force_remake_before"
                    continue
                end
            end
            
            t1 = time()

            run_parameters = Altim.binned_filled_fileparts(synthesized_gemb_fit)
            sm = run_parameters.surface_mask;

            dh = FileIO.load(binned_synthesized_file, "dh_hyps")
    
            # convert elevation change to volume change
            dv_altim = Altim.dh2dv(dh, geotiles[sm]);
  
            if any(all(isnan.(dv_altim), dims = :date))
                # keep this error message
                dgeotile = dims(dv_altim, :geotile)
                all_nan_geotiles = dgeotile[findall(dropdims(all(isnan.(dv_altim), dims = :date), dims = :date))]
                printstyled("Geotiles that only contain NaNs:\n"; color = :red)
                println("$(collect(all_nan_geotiles)))")
                printstyled("Possible sources of error include:\n"; color = :red)
                printstyled("  [1] synthesis_error_file was run on a subset of files that do not include the full error range of the data, to fix this you need to delete the exisiting error file and any synthesis files that were created :\n [2] the wrong surface mask was passed to Altim.dh2dv\n"; color = :red)
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
        
            df = Altim.gemb_bestfit_grouped(dv_altim, smb, fac, discharge, geotiles0; examine_model_fits)

            df[!,:area_km2] = sum.(geotiles[sm].area_km2)
            df.mad = df.mad ./ df.area_km2
            rename!(df, "mad"=>"mad_m")
            df.geometry = Altim.extent2rectangle.(df.extent)
            df = df[:,Not(:extent)]
            GeoDataFrames.write(synthesized_gemb_fit, df)
            println("\n $binned_synthesized_file optimal GEMB fit found: $(round(Int,time() -t1))s")
        end
    end
end

# DOWNSCALE OPTIMAL GEMB MODEL TO INDIVIDUAL GLACIERS
# 1. Loading the GEMB fit parameters and height change data for each geotile
# 2. Converting geotile height changes to per-glacier values
# 3. Scaling GEMB variables (SMB, FAC, runoff) using fit parameters and converting to per-glacier values
# This takes ~8 min minutes per input file [5.5 hours for all]
if true

    @showprogress desc="Downscaling optimal GEMB data to individual glaciers..." for binned_synthesized_file in path2runs
        #binned_synthesized_file = reference_run

        run_parameters = Altim.binned_filled_fileparts(binned_synthesized_file)
        sm = run_parameters.surface_mask;

        glaciers0 = copy(glaciers[sm])
        geotiles0 = copy(geotiles[sm])
        
        synthesized_gemb_fit = replace(binned_synthesized_file,  ".jld2" => "_gemb_fit.arrow")
        perglacier_synthesized_file = replace(binned_synthesized_file, ".jld2" => "_perglacier.jld2")
        gemb_fit = GeoDataFrames.read(synthesized_gemb_fit)

        if !isfile(perglacier_synthesized_file) || !isnothing(force_remake_before)

            if isfile(perglacier_synthesized_file) && !isnothing(force_remake_before)
                if Dates.unix2datetime(mtime(perglacier_synthesized_file)) > force_remake_before
                    @warn "Skipping $(perglacier_synthesized_file) because it was created after $force_remake_before"
                    return
                end
            end

            ### DOWNSCALE TO EACH GLACIER BY HYPSOMETRY 
            # NOTE: you get strange spatial gradients in smb due to a lack of spatial gradents in 
            # climate forcing [like leeward side of mountain ranges]
            # therefore it is currently recommended to use the downscale by "area" method

            #println("start downscaling from geotile to glacier which will take ~8 minutes per file...")
            t1 = time()

            if downscale_to_glacier_method == "hypsometry"

                dh = load(binned_synthesized_file, "dh_hyps")
            
                glaciers0 = Altim.geotile2glacier!(glaciers0, dh; varname = "dh");

                # load example geotile to get dimensions
                filename_gemb_geotile_filled = replace(filename_gemb_geotile, ".jld2" => "_$(gemb_fit.id[1])_filled.jld2")
                smb = load(filename_gemb_geotile_filled, "smb")
            
                dgeotile = dims(dh, :geotile)
                ddate = dims(smb, :date)
                dheight = dims(smb, :height)
                dpscale = dims(smb, :pscale)
                dΔheight = dims(smb, :Δheight)

                var = fill(NaN, dgeotile, ddate, dheight)

                # scale smb, fac, runoff [NOTE: all are in unit of m i.e.]
                for varname = ["smb", "fac", "runoff"]
                    @showprogress desc="Populate calibrated $(varname)..." Threads.@threads for gt in eachrow(gemb_fit)
                        filename_gemb_geotile_filled = replace(filename_gemb_geotile, ".jld2" => "_$(gt.id)_filled.jld2")
                        var0 = load(filename_gemb_geotile_filled, varname)
                        var[At(gt.id), :, :] = var0[:,:,At(gt.pscale), At(gt.Δheight)]
                    end

                    glaciers0 = Altim.geotile2glacier!(glaciers0, var; varname)
                end
            elseif downscale_to_glacier_method == "area"
                ### DOWNSCALE TO EACH GLACIER FRACTIONAL AREA WITHIN EACH GEOTILE
                glacier_area_km2 = sum.(glaciers0.area_km2)

                for varname = ["dh", "smb", "fac", "runoff"]

                    if varname == "dh"
                        dh = load(binned_synthesized_file, "dh_hyps")
                        var0 = Altim.dh2dv(dh, geotiles0);
                    else
                        var0 = load(filename_gemb_geotile_filled_dv, varname )
                    end
                    
                    ddate = dims(var0, :date)
                    glaciers0[!,varname] =  [fill(NaN, ddate) for _ in 1:nrow(glaciers0)]

                    for geotile in unique(glaciers0.geotile)
                    #geotile = "lat[+28+30]lon[+082+084]"
                        gindex = findall(geotile .== glaciers0.geotile)
                        if varname != "dh"
                            gfit = gemb_fit[findfirst(isequal(geotile), gemb_fit.id),:]
                            gt_dv = var0[At(geotile), :, At(gfit.pscale), At(gfit.Δheight)]
                        else
                            gt_dv = var0[At(geotile), :, :]
                        end

                        garea = sum.(glaciers0[gindex, :area_km2])
                        gweighting = garea ./ sum(garea)
                        for (i,ig) in enumerate(gindex)
                            glaciers0[ig,varname][:] = gt_dv * gweighting[i] ./ garea[i] * 1000 # convert from km3 to m i.e.
                        end
                    end
                end 
            else
                error("downscale_to_glacier_method must be either \"hypsometry\" or \"area\"")
            end

            FileIO.save(perglacier_synthesized_file, Dict("glaciers" => glaciers0))

            fn = splitpath(perglacier_synthesized_file)
            #println("downscaling from geotile to glacier complete $(round(Int,time() -t1))s: $(fn[end])")
        end
    end
end


# This section processes synthesized altimetry and GEMB data at the geotile level
# Key steps:
#
# 1. Data Loading and Setup:
#    - Loads discharge data and mass conversion factors
#    - Processes each run file from path2runs
#    - Copies glacier and geotile data for current surface mask
#
# 2. Geotile Processing:
#    - Filters geotiles to only those containing glaciers
#    - Assigns group IDs from geotile_groups
#    - Adds GEMB fit parameters and mass conversion factors
#    - Processes discharge data for each geotile
#
# 3. Variable Processing:
#    - Handles multiple variables: dv_altim, runoff, fac, smb, etc.
#    - Converts units and applies GEMB scaling
#    - Calculates derived variables like dv, dm, dm_altim
#
# 4. Assign rgi ids and save:
#    - Averages discharge, dv, dm by geotile groups
#    - Assigns RGI region IDs
#    - Saves processed data to JLD2 files
#

begin
    # Load synthesized data and GEMB fit parameters
    discharge = load(globaldischarge_fn, "discharge")
    volume2mass = Altim.δice / 1000

    has_glacier = Dict()
    for keys in keys(geotiles)
        has_glacier[keys] = [id in unique(glaciers[keys].geotile) for id in geotiles[keys].id]
    end

    # removed Threads.@threads due to memory issues
    @showprogress desc="Computing calibrated geotile level data timeseries for all runs..." for binned_synthesized_file in path2runs 

    #for binned_synthesized_file in [reference_run]
    #binned_synthesized_file = reference_run
        
        binned_synthesized_dv_file = replace(binned_synthesized_file, ".jld2" => "_gembfit_dv.jld2")

        run_parameters = Altim.binned_filled_fileparts(binned_synthesized_file)
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
            Δdecyear = Altim.decimalyear.(ddate) .- Altim.decimalyear.(ddate[1])
            geotiles0[!, :discharge] .=[fill(NaN, length(ddate)) for _ in 1:nrow(geotiles0)]
            
            for geotile in eachrow(geotiles0)
                    fit_index = findfirst(isequal(geotile.id), gemb_fit.id)
                    geotile.pscale = gemb_fit[fit_index, :pscale]
                    geotile.Δheight = gemb_fit[fit_index, :Δheight]
                    area_km2 = sum(geotile.area_km2)
                    geotile.mie2cubickm = area_km2/1000  # Convert meters ice equivalent to cubic kilometers

                    # include discharge average over total glacier area [i.e. save in units of mie]
                    index = Altim.within.(Ref(geotile.extent), discharge.longitude, discharge.latitude)
                    geotile.discharge = sum(discharge.discharge_gtyr[index])/volume2mass * Δdecyear # to units of mie
            end

            # add discharge date metadata
            colmetadata!(geotiles0, "discharge", "date", collect(ddate), style=:note)

            # Process each variable
            for varname in varnames

                # Special handling for height change data
                if varname == "dv_altim"
                    dh = load(binned_synthesized_file, "dh_hyps")
                    var0 = Altim.dh2dv(dh, geotiles0);
                else
                    var0 = load(filename_gemb_geotile_filled_dv, varname)
                end
                
                # Initialize time series arrays
                ddate = dims(var0, :date)
                geotiles0[!, varname] = [fill(NaN, length(ddate)) for _ in 1:nrow(geotiles0)]

                # Add date metadata for time series
                colmetadata!(geotiles0, varname, "date",collect(ddate), style=:note)
                
                # Apply GEMB scaling and convert units for each geotile
                for geotile in eachrow(geotiles0)
                    geotile[varname] = var0[At(geotile.id), :, At(geotile.pscale),At(geotile.Δheight)];
                end
            end

            # modeled dv
            geotiles0[!, :dv] = (geotiles0[:, :smb] .- geotiles0[:, :discharge] .+ geotiles0[:, :fac])
            colmetadata!(geotiles0, "dv", "date", collect(ddate), style=:note)

            # modeled dm
            geotiles0[!, :dm] = (geotiles0[:, :smb] .- geotiles0[:, :discharge]) .* volume2mass # to units of mwe [kg/m²]
            colmetadata!(geotiles0, "dm", "date", collect(ddate), style=:note)

            # altimetry dm
            geotiles0[!, :dm_altim] = copy(geotiles0[!, :dv_altim])
            for r in eachrow(geotiles0)
                model = LinearInterpolation(r.fac, Altim.decimalyear.(colmetadata(geotiles0, "fac", "date")))
                fac = model(Altim.decimalyear.(colmetadata(geotiles0, "dv_altim", "date")))
                r.dm_altim = (r.dv_altim .- fac) .* volume2mass # in units of mwe [kg/m²]
            end

            colmetadata!(geotiles0, "dm_altim", "date", colmetadata(geotiles0, "dv_altim", "date"), style=:note)

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

# export trends and amplitudes for plotting of reference_run only
begin     
    binned_synthesized_file = reference_run
    binned_synthesized_dv_file = replace(binned_synthesized_file, ".jld2" => "_gembfit_dv.jld2")

    outfile = joinpath(paths.data_dir, "project_data", "geotiles_rates.fgb");

    if !isfile(outfile) || !isnothing(force_remake_before)

        if isfile(outfile) && !isnothing(force_remake_before)
            if Dates.unix2datetime(mtime(outfile)) > force_remake_before
                @warn "Skipping $(outfile) because it was created after $force_remake_before"
                return
            end
        end

        geotiles0 = FileIO.load(binned_synthesized_dv_file, "geotiles")

        vars_no_write = setdiff(names(geotiles0), ["id", "glacier_frac", "landice_frac", "floating_frac", "geometry", "group", "pscale", "Δheight", "mie2cubickm", "rgiid"])
        vars_ts = setdiff(vars_no_write, ["extent", "area_km2"])
        

        # Fit temporal trends to all variables
        geotiles0 = Altim.df_tsfit!(geotiles0, vars_ts; datelimits = (DateTime(2000,1,1), DateTime(2023,1,1)))
        
        # Export results, excluding raw time series
        source_crs1 = GFT.EPSG(4326)
        isvec = []
        for i = 1:ncol(geotiles0)
            push!(isvec, typeof(geotiles0[1,12]) >: Vector)
        end

        GeoDataFrames.write(joinpath(paths.data_dir, "project_data", "geotiles_rates_km3yr.fgb"), geotiles0[:, Not(vars_no_write)]; crs=source_crs1)


        # now do the same but for area averaged rates of change 
        for varname in vars_ts
            geotiles0[:, varname] ./= geotiles0[:, :mie2cubickm]
        end

        # Fit temporal trends to all variables (easier than finding all of the fits and then multiplying by mie2cubickm)
        geotiles0 = Altim.df_tsfit!(geotiles0, vars_ts; datelimits = (DateTime(2000,1,1), DateTime(2025,1,1)))


        # plot a histogram of the rates
        f = Figure();

        ax1 = f[1, 1] = Axis(f; xlabel = "trend [m i.e. yr⁻¹]", ylabel = "count")
        CairoMakie.stephist!(geotiles0[:, :dv_trend], bins = -5:0.25:5; label = "dh")
        CairoMakie.stephist!(geotiles0[:, :dv_altim_trend], bins = -5:0.25:5; label = "dh altim"); 
        axislegend(ax1, framevisible = false); 

        ax2 = f[1, 2] = Axis(f; xlabel = "amplitude [m i.e.]", ylabel = "count")
        CairoMakie.stephist!(geotiles0[:, :dv_amplitude], bins = 0:0.25:5; label = "dh")
        CairoMakie.stephist!(geotiles0[:, :dv_altim_amplitude], bins = 0:0.25:5; label = "dh altim"); 
        axislegend(ax2, framevisible = false); 
       
        display(f)

        # plot a histogram of pscale and Δheight
        synthesized_gemb_fit = replace(binned_synthesized_file, ".jld2" => "_gemb_fit.arrow")
        gemb_fit = GeoDataFrames.read(synthesized_gemb_fit);

        f = Figure();
        ax1 = f[1, 1] = Axis(f; xlabel = "pscale", ylabel = "count")
        CairoMakie.stephist!(gemb_fit[:, :pscale], bins = 0.25:0.25:4)
        ax2 = f[1, 2] = Axis(f; xlabel = "Δheight [m]", ylabel = "count")
        CairoMakie.stephist!(gemb_fit[:, :Δheight], bins = -3000:200:3000); 
        #axislegend(ax1, framevisible = false); 
        display(f)

        CairoMakie.save(replace(outfile, ".fgb" => ".png"), f)

        GeoDataFrames.write(outfile, geotiles0[:, Not(vars_no_write)]; crs=source_crs1)
    end

    # rsync -rav devon:/mnt/bylot-r3/data/project_data/ ~/data/Altim/project_data/
    # Sync command for reference
end
