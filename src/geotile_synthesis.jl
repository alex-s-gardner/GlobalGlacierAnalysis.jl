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

# set force_remake == true to redo all steps from scratch
force_remake = false;

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
    #    - Geotile width of 2 degrees
    #    - Project ID, surface masks, DEM sources, and processing options for uncertainty analysis
    #    - Latitude threshold for discharge-to-SMB conversion
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
    dem_id=["best" "cop30_v2"]
    curvature_correct=[false true]
    amplitude_correct=[true]
    binning_method=["median" "meanmadnorm5" "meanmadnorm3"] # ["median" "meanmadnorm10" "meanmadnorm5" "meanmadnorm3"]
    paramater_set=[1, 2, 3, 4]
    binned_folder=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg")

    filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0.jld2"
    filename_gemb_geotile = replace(filename_gemb_combined, ".jld2" => "_geotile.jld2")
    #filename_gemb_geotile_filled_dv = replace(filename_gemb_combined, ".jld2" => "_geotile_filled_dv.jld2")
    filename_gemb_geotile_filled_dv = replace(filename_gemb_combined, ".jld2" => "_geotile_filled_extra_dv.jld2")

    geotile_groups_fn = joinpath(paths.data_dir, "project_data", "geotile_groups.arrow")

    globaldischage_fn = joinpath(paths.data_dir, "GlacierOutlines/GlacierDischarge/global_glacier_dischage.jld2")

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
begin
    for sm in surface_mask 
        # Set the glacier ID attribute name based on RGI version
        persistent_attribute = if Base.contains(geomfile[sm], "rgi6")
            :RGIId
        else
            :rgi_id
        end

        # Only process if output file doesn't exist or force_remake is true
        if !isfile(glacier_geotile_hyps_fn[sm]) || force_remake

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

if .!isfile(globaldischage_fn) || force_remake # Load RGI6 glacier hypsometry data from previously generated file glacier_geotile_hyps_fn = joinpath(paths.data_dir, "GlacierOutlines/rgi60/rgi60_Global.gpkg")
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
    save(globaldischage_fn, Dict("discharge"=>discharge))
else
    # If file exists, just load the discharge data
    discharge = load(globaldischage_fn, "discharge")
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
path2geotile_synthesis_error = geotile_synthesis_error(; 
    path2runs, 
    outfile="/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2", 
    force_remake
)

# SYNTHESIS OF ALL MISSIONS
# 1. Takes ~2 hours to process 252 run files
# 2. Combines data from specified missions:
#    - Hugonnet glacier mass balance
#    - GEDI laser altimetry 
#    - ICESat laser altimetry
#    - ICESat-2 laser altimetry
# 3. Uses synthesis error calculations from previous step to:
#    - Weight contributions from each mission
#    - Account for varying uncertainties
# 4. Can force reprocessing of synthesis with force_remake flag
geotile_synthesize_runs(;
    path2runs,
    path2geotile_synthesis_error, 
    missions2include=["hugonnet" "gedi" "icesat" "icesat2"],
    force_remake
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
#    - Indices of intersecting glaciers (dischage_ind)
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

if !isfile(geotile_groups_fn) || force_remake
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
#    - Handles troubleshooting for specific geotiles if needed
#
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
@time begin

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
    #for binned_synthesized_file in path2runs

        synthesized_gemb_fit = replace(binned_synthesized_file,  ".jld2" => "_gemb_fit.arrow")

        if !(isfile(synthesized_gemb_fit)) || force_remake
            t1 = time()

            run_parameters = Altim.binned_filled_fileparts(synthesized_gemb_fit)
            sm = run_parameters.surface_mask;

            dh = FileIO.load(binned_synthesized_file, "dh_hyps")
    
            # convert elevation change to volume change
            dv_altim = Altim.dh2dv(dh, geotiles[sm]);

            # sanity check
            if false
                gtids = ["lat[-14-12]lon[-076-074]", "lat[-10-08]lon[-078-076]", "lat[+28+30]lon[+082+084]", "lat[+58+60]lon[-136-134]", "lat[+68+70]lon[-070-068]", "lat[+62+64]lon[-044-042]", "lat[+32+34]lon[+076+078]", "lat[+28+30]lon[+092+094]", "lat[-46-44]lon[+168+170]", "lat[-70-68]lon[-068-066]", "lat[-52-50]lon[-074-072]", "lat[+58+60]lon[-136-134]"]

                gtids = ["lat[+38+40]lon[+070+072]"]
                for gtid in gtids
                    gemb_fit = Altim.gemb_bestfit_grouped(dv_altim, smb, fac, discharge, geotiles0; examine_model_fits = gtid) 
                end
            end
            
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
            # df = Altim.gemb_bestfit(dv_altim, smb, fac, discharge, geotiles[sm])

            # there are issues with calibrating the smb model to individual geotiles glaciers 
            #can cross multiple geotiles therefore we calibrate the model for groups of 
            # distinct geotiles.

            # using geotiles0 instead of geotiles[sm] is fine here
            df = Altim.gemb_bestfit_grouped(dv_altim, smb, fac, discharge, geotiles0)

            df[!,:area_km2] = sum.(geotiles[sm].area_km2)
            df.mad = df.mad ./ df.area_km2
            rename!(df, "mad"=>"mad_m")
            df.geometry = Altim.extent2rectangle.(df.extent)
            df = df[:,Not(:extent)]
            GeoDataFrames.write(synthesized_gemb_fit, df)
            println("$binned_synthesized_file optimal GEMB fit found: $(round(Int,time() -t1))s")
        end
    end
end

# DOWNSCALE OPTIMAL GEMB MODEL TO INDIVIDUAL GLACIERS
# 1. Loading the GEMB fit parameters and height change data for each geotile
# 2. Converting geotile height changes to per-glacier values
# 3. Scaling GEMB variables (SMB, FAC, runoff) using fit parameters and converting to per-glacier values
# This takes ~8 min minutes per input file
@time begin

    #for binned_synthesized_file in path2runs
    binned_synthesized_file = reference_run

    run_parameters = Altim.binned_filled_fileparts(binned_synthesized_file)
    sm = run_parameters.surface_mask;

    glaciers0 = copy(glaciers[sm])
    
    synthesized_gemb_fit = replace(binned_synthesized_file,  ".jld2" => "_gemb_fit.arrow")
    perglacier_synthesized_file = replace(binned_synthesized_file, ".jld2" => "_perglacier.jld2")
    gemb_fit = GeoDataFrames.read(synthesized_gemb_fit)

    if !isfile(perglacier_synthesized_file) || force_remake
        ### DOWNSCALE TO EACH GLACIER BY HYPSOMETRY 
        # NOTE: you get strange spatial gradients in smb due to a lack of spatial gradents in 
        # climate forcing [like leeward side of mountain ranges]
        # therefore it is currently recommended to use the downscale by "area" method

        println("start downscaling from geotile to glacier which will take ~8 minutes per file...")
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
        println("downscaling from geotile to glacier complete $(round(Int,time() -t1))s: $(fn[end])")
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
    discharge = load(globaldischage_fn, "discharge")
    volume2mass = Altim.δice / 1000

    has_glacier = Dict()
    for keys in keys(geotiles)
        has_glacier[keys] = [id in unique(glaciers[keys].geotile) for id in geotiles[keys].id]
    end

    # removed Threads.@threads due to memory issues
    @showprogress desc="Computing calibrated geotile level data timeseries for all runs..." for binned_synthesized_file in path2runs 
    #binned_synthesized_file = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_dh_cop30_v2_median_v01_filled_ac_p1_synthesized.jld2"

        binned_synthesized_dv_file = replace(binned_synthesized_file, ".jld2" => "_gembfit_dv.jld2")

        println(binned_synthesized_file)
        run_parameters = Altim.binned_filled_fileparts(binned_synthesized_file)
        sm = run_parameters.surface_mask

        glaciers0 = copy(glaciers[sm])

        if !isfile(binned_synthesized_dv_file) || force_remake

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
                r.dm_altim = (r.dm_altim .- fac) .* volume2mass # in units of mwe [kg/m²]
            end

            colmetadata!(geotiles0, "dm_altim", "date", colmetadata(geotiles0, "dv_altim", "date"), style=:note)

            varnames_all = vcat(varnames, "discharge", "dv", "dm", "dm_altim")
            
            # `dichage`, `dm` and `dv` must be average by geotile groups as they are not valid for single geotiles
            vars2average = ["discharge", "dv", "dm"]
            gdf = groupby(geotiles0, :group)
            for g in gdf
                #g = gdf[4]
                if nrow(g) > 1
                    for varn in vars2average
                        
                        foo = g[1,varn] .* g[1, :mie2cubickm]
                        for i = 2:nrow(g)
                            foo .+= g[i, varn] .* g[i, :mie2cubickm]
                        end
                        foo = foo ./ sum(g[!, :mie2cubickm])
                        
                        for i = 1:nrow(g)
                            g[i, varn] = foo
                        end           
                    end
                end
            end

            #### now compute regional estimates
            geotiles0[!, :rgiid] .= Int8(0)
            for i = 1:19
                rgi = "rgi$i"
                geotiles0[geotiles0[:, rgi].>0, :rgiid] .= i
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
    geotiles0 = Altim.df_tsfit!(geotiles0, vars_ts; datelimits = (DateTime(2000,1,1), DateTime(2023,1,1)))

    GeoDataFrames.write(joinpath(paths.data_dir, "project_data", "geotiles_rates.fgb"), geotiles0[:, Not(vars_no_write)]; crs=source_crs1)

    # Sync command for reference
    # rsync -rav devon:/mnt/bylot-r3/data/project_data/ ~/data/Altim/project_data/
end

# make multi-region plot
#begin
    

    

    # read in example file
    regions = FileIO.load(replace(reference_run, ".jld2" => "_gembfit_dv.jld2"), "geotiles")
    varnames = setdiff(names(geotiles0), ["id", "extent", "glacier_frac", "area_km2", "discharge", "landice_frac", "floating_frac", "geometry", "group", "pscale", "Δheight", "mie2cubickm", "rgiid"])

    drun = Dim{:run}([last(splitpath.(path2run)) for path2run in path2runs])
    drgi = Dim{:rgi}([1:19;98])
    regions0 = Dict()
    for varname in varnames
        ddate = Dim{:date}(colmetadata(regions, varname, "date"))
        regions0[varname] = fill(NaN, drun, drgi, ddate)
    end


    @time for binned_synthesized_file in path2runs
    #binned_synthesized_file = path2runs[1]
        run = splitpath(binned_synthesized_file)[end]
        binned_synthesized_dv_file = replace(binned_synthesized_file, ".jld2" => "_gembfit_dv.jld2")
        
        println(run)

        geotiles0 = FileIO.load(binned_synthesized_dv_file, "geotiles")

        # group by rgi and sum
        geotiles_reg = groupby(geotiles0, :rgiid)

        # sum across timeseries
        regions = DataFrames.combine(geotiles_reg, varnames .=> Ref ∘ sum; renamecols=false)

        # add an HMA region 
        index = (regions[:, :rgiid] .<= 15) .& (regions[:, :rgiid] .>= 13)
        regions = append!(regions, DataFrames.combine(regions[index, :], vcat("rgiid", varnames) .=> Ref ∘ sum; renamecols=false))
        regions[end, :rgiid] = 98

        for varname in varnames
            regions0[varname][At(run),At(regions.rgiid),:] = reduce(hcat,regions[!,varname])'
        end
    end

    # create df for reference run
    binned_synthesized_file = reference_run
    binned_synthesized_dv_file = replace(binned_synthesized_file, ".jld2" => "_gembfit_dv.jld2")
    
    geotiles0 = FileIO.load(binned_synthesized_dv_file, "geotiles")

    # group by rgi and sum
    geotiles_reg = groupby(geotiles0, :rgiid)

    # sum across timeseries
    regions = DataFrames.combine(geotiles_reg, varnames .=> Ref ∘ sum; renamecols=false)

    # add an HMA region 
    index = (regions[:, :rgiid] .<= 15) .& (regions[:, :rgiid] .>= 13)
    regions = append!(regions, DataFrames.combine(regions[index, :], vcat("rgiid", varnames) .=> Ref ∘ sum; renamecols=false))
    regions[end, :rgiid] = 98

    # copy metadata
    for varname in names(regions)[2:end]
        cmetadata = colmetadata(geotiles0, varname, "date")
        regions = colmetadata!(regions, varname, "date", cmetadata; style=:note)
    end

    # add error columns to regions
    #for varname in varnames
    varname = "dv"
        varname_err = "$(varname)_err"

        foo = std(regions0[varname], dims = :run)

        ### There are lots of regaions with all NaNs.. and volume change needs to be centered around zero
        CairoMakie.heatmap(regions0[varname][:,9,:])

        regions[!, varname_err] = eachrow(collect(dropdims(std(regions0[varname], dims = :run), dims = :run)))
        regions = colmetadata!(regions, varname_err, "date", colmetadata(regions, varname, "date"); style=:note)
    end
    
    # read in GRACE data
    grace = Altim.read_grace_rgi(; datadir=setpaths()[:grace_rgi])

    regions[!, :dm_grace] = [fill(NaN, length(grace["rgi1"]["dM_gt_mdl_fill"])) for _ in 1:nrow(regions)]
    gdate = vec(Altim.datenum2date.(grace["rgi1"]["dM_gt_mdl_fill_date"]));
    regions = colmetadata!(regions, :dm_grace, "date", gdate; style=:note)
    regions[!, :dm_grace_err] = [fill(NaN, length(grace["rgi1"]["dM_gt_mdl_fill"])) for _ in 1:nrow(regions)]
    regions = colmetadata!(regions, :dm_grace_err, "date", gdate; style=:note)

    for r in eachrow(regions)
        #r = first(eachrow(regions))
        if r.rgiid == 98
            rgi = "HMA"
        else
            rgi = "rgi$(r.rgiid)"
        end

        haskey(grace, rgi) || continue

        r.dm_grace = vec(grace[rgi]["dM_gt_mdl_fill"])
        r.dm_grace_err = vec(grace[rgi]["dM_sigma_gt_mdl_fill"])
    end

    # reformat for plotting
    varnames = setdiff(names(regions), ["rgiid"])
    dates = DateTime(2000, 1, 15):Month(1):DateTime(2023, 1, 15)
    dates_decyear = Altim.decimalyear.(dates)
    df = DataFrame()

    for varname in varnames
        #varname = "discharge"

        date0 = Altim.decimalyear.(colmetadata(regions, varname, "date"))
        v = fill(NaN, length(dates_decyear))

        if Base.contains(varname, "altim")
            mission = "synthesis"
        elseif Base.contains(varname, "grace")
            mission = "grace"
        else
            mission = "gemb"
        end

        if Base.contains(varname,"dm")
            if !Base.contains(varname,"err")
                varname_out = "dm"
            else
                varname_out = "dm_err"
            end
        elseif Base.contains(varname,"dv")
            if !Base.contains(varname,"err")
                varname_out = "dv"
            else
                varname_out = "dv_err"
            end
        else
            varname_out = varname
        end

        for r in eachrow(regions)
            #r = first(eachrow(regions))
            model = LinearInterpolation(r[varname], date0)
            index = (dates_decyear .> minimum(date0)) .& (dates_decyear .< maximum(date0))
            v[index] = model(dates_decyear[index])
            varname_mid = "$(varname)"
            df = vcat(df, DataFrame("rgi" => r[:rgiid], "mission" => mission, "var" => varname_out, "mid" => Ref(copy(v))))
        end
    end

    # set GRACE to NaN for Greenland and Antarctica
    index = findfirst((df.rgi .== 5) .& (df.mission .== "grace"))
    df[index, :mid] .= NaN
    index = findfirst((df.rgi .== 19) .& (df.mission .== "grace"))
    df[index, :mid] .= NaN

    df[!, :low] = deepcopy(df[!, :mid])
    df[!, :high] = deepcopy(df[!, :mid])
    df[!, :unit] .= "gt"
    metadata!(df, "date", dates; style=:note)

    gdf = groupby(df, "rgi")

    #for g in gdf
    g = gdf[3]
        #for mission in unique(g.mission)
        mission = unique(g.mission)[1]
            mission_index = g.mission .== mission

            vars = unique(g[mission_index, :var])
            vars = vars[.!Base.contains.(vars, Ref("err"))]

            #for var0 in vars
            var0 = vars[2]
                var_index = g.var .== var0
                ind_mid = findfirst(var_index .& mission_index)
                ind_err = findfirst((g.var .== "$(var0)_err") .& mission_index)
                if !isnothing(ind_err)
                    g[ind_mid, :low] .-= g[ind_err, :mid]
                    g[ind_mid, :high] .+= g[ind_err, :mid]
                end
            end
        end
    end

    # remove "dm_grace_err"
    df = df[.!Base.contains.(df.var, Ref("err")), :]

    # align to altim
    gdf = groupby(df, "rgi")
    for g in gdf
        #g = gdf[1]
        index_var = g.var .== "dm"

        grace0 = g[findfirst((g.mission .== "grace") .& index_var), :mid]
        altim0 = g[findfirst((g.mission .== "synthesis") .& index_var), :mid]
        gemb0 = g[findfirst((g.mission .== "gemb") .& index_var), :mid]

        index = .!isnan.(altim0)
        for ts in (grace0, gemb0)
            index2 = .!isnan.(ts)
            if any(index2)
                index = index .& index2
            end
        end

        delta = mean(altim0[index])
        g[findfirst((g.mission .== "synthesis") .& index_var), :mid] .-= delta
        g[findfirst((g.mission .== "synthesis") .& index_var), :low] .-= delta
        g[findfirst((g.mission .== "synthesis") .& index_var), :high] .-= delta

        delta = mean(grace0[index])
        g[findfirst((g.mission .== "grace") .& index_var), :mid] .-= delta
        g[findfirst((g.mission .== "grace") .& index_var), :low] .-= delta
        g[findfirst((g.mission .== "grace") .& index_var), :high] .-= delta

        delta = mean(gemb0[index])
        g[findfirst((g.mission .== "gemb") .& index_var), :mid] .-= delta
        g[findfirst((g.mission .== "gemb") .& index_var), :low] .-= delta
        g[findfirst((g.mission .== "gemb") .& index_var), :high] .-= delta
    end

    # exclude rgi 13, 14, 15
    index = .!((df.rgi .== 13) .| (df.rgi .== 14) .| (df.rgi .== 15))
    df = df[index, :]


exclude_mission = df.mission .== "gemb"
f, rgi, offset, ylims = Altim.plot_multiregion_dvdm(df[.!exclude_mission, :];
    variable="dm",
    featured_mission="synthesis",
    regions=unique(df.rgi),
    showlines=false,
    showmissions=true,
    fontsize=15,
    cmap=:Dark2_4,
    regions_ordered=false,
    region_offsets=nothing,
    ylims=nothing,
    title=nothing,
    palette=nothing,
    delta_offset=nothing,
)
f

end
