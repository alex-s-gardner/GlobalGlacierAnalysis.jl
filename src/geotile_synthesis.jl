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

# This section imports required packages for geotile synthesis:
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
    include("utilities_synthesis.jl")
end


# This section sets up paths and parameters for geotile synthesis:
#
# 1. Loads local configuration paths
# 2. Sets key parameters:
#    - Geotile width of 2 degrees
#    - Project ID, surface masks, DEM sources, and processing options for uncertainty analysis
#    - Latitude threshold for discharge-to-SMB conversion
#    - Time period for equilibrium calculations
# 3. Builds parameter combinations and filters to only include existing files
# 4. Sets up paths for glacier geometry files (RGI v6 and v7)
begin
    paths = Altim.pathlocal
    geotile_width = 2;
    downscale_to_glacier_method = "area"

    # to include in uncertainty
    project_id = ["v01"]
    surface_mask=["glacier", "glacier_rgi7", "glacier_b1km"]
    dem_id=["best" "cop30_v2"]
    curvature_correct=[false true]
    amplitude_correct=[true]
    binning_method=["median" "meanmadnorm5" "meanmadnorm3" "meanmadnorm10"] # ["median" "meanmadnorm10" "meanmadnorm5" "meanmadnorm3"]
    paramater_set=[1, 2, 3, 4]
    binned_folder=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg")

    filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0.jld2"
    filename_gemb_geotile = replace(filename_gemb_combined, ".jld2" => "_geotile.jld2")

    geotile_groups_fn = joinpath(paths.data_dir, "project_data", "geotile_groups.gpkg")

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
    glacier_geom = Dict();
    geomfile_rgi6 = joinpath(paths.data_dir, "GlacierOutlines/rgi60/rgi60_Global.gpkg")
    geomfile_rgi7 = joinpath(paths.data_dir, "GlacierOutlines/RGI2000-v7.0-G-global-fix/rgi70_Global.gpkg")
end
# Process both RGI6 and RGI7 glacier outlines
# This section processes glacier hypsometry data for both RGI6 and RGI7 glacier inventories
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
    for geomfile_rgi = [geomfile_rgi6, geomfile_rgi7]
        # Generate output filename for hypsometry data
        glacier_geotile_hyps_fn = replace(geomfile_rgi, ".gpkg" => "geotile_hyps.jld2")

        # Set the glacier ID attribute name based on RGI version
        persistent_attribute = if Base.contains(geomfile_rgi, "rgi6")
            :RGIId
        else
            :rgi_id
        end

        # Only process if output file doesn't exist
        if !isfile(glacier_geotile_hyps_fn) || force_remake
            # Read glacier outlines
            glacier_geom = GeoDataFrames.read(geomfile_rgi)
            
            # Load elevation data and define height bins
            h = Raster(paths.cop30_v2, lazy=true)
            height_range, height_center = Altim.project_height_bins()
            
            # Get geotiles that contain glaciers
            geotiles = Altim.geotiles_w_mask(geotile_width)
            geotiles = geotiles[geotiles.glacier_frac .> 0.0, :]

            # Calculate hypsometry for each glacier-geotile combination
            # This is computationally intensive (~37 min)
            @time geotile_zonal_area_hyps(h, height_range, glacier_geom, geotiles.id; persistent_attribute)
            
            # Save results
            FileIO.save(glacier_geotile_hyps_fn, Dict("glaciers" => glacier_geom))
        end
    end 
end

## measured dischage and dischage for unmeasured areas by setting dischage == unscaled smb
# This section processes glacier discharge data by:
# 1. Setting up paths and checking if global discharge file exists
begin
    globaldischage_fn = joinpath(paths.data_dir, "GlacierOutlines/GlacierDischarge/global_glacier_dischage.jld2")
    if .!isfile(globaldischage_fn) || force_remake
        # Load RGI6 glacier hypsometry data from previously generated file
        glacier_geotile_hyps_fn = joinpath(paths.data_dir, "GlacierOutlines/rgi60/rgi60_Global.gpkg")
        glacier_geotile_hyps_fn = replace(glacier_geotile_hyps_fn, ".gpkg" => "geotile_hyps.jld2")
        glaciers = load(glacier_geotile_hyps_fn, "glaciers")

        # Load GEMB surface mass balance data that has been filled at geotile level
        filename_gemb_geotile_filled = replace(filename_gemb_combined, ".jld2" => "_geotile_filled.jld2")
        gemb = Dict("smb" => load(filename_gemb_geotile_filled, "smb"))

        # Get geotiles containing glaciers and initialize parameters:
        # - Filter to only tiles with glaciers
        # - Set precipitation scaling to 1.0 (no scaling)
        # - Set height adjustment to 0
        geotiles = Altim.geotiles_w_mask(geotile_width)
        geotiles = geotiles[(geotiles.glacier_frac.>0.0), :]
        geotiles[!, :pscale] .= 1.0 # no precipitation scaling
        geotiles[!, :Δheight] .= 0.0

        # Convert geotile-level data to per-glacier values
        @time glaciers = Altim.geotile2glacier!(glaciers, gemb, geotiles)

        # For unmeasured glaciers, estimate discharge using surface mass balance
        # Only process glaciers with non-zero area
        discharge0 = Altim.discharge2smb(
            glaciers[sum.(glaciers.area_km2).>0, :];
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
end

# This section calculates synthesis errors for combining multiple altimetry missions:
#
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
    outfile="/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2"; 
    force_remake
)

# This section synthesizes elevation change data from multiple altimetry missions:
#
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
    missions2include=["hugonnet" "gedi" "icesat" "icesat2"];
    force_remake
)

# Load geotile info for mutiple surface masks
begin
    # Update paths to point to synthesized data
    path2runs = replace.(path2runs, "aligned.jld2" => "synthesized.jld2")

    # Set parameters
    geotile_width = 2

    # Load elevation change data
    dh = FileIO.load(path2runs[1], "dh_hyps")

    # Get geotile dimensions
    dgeotile = dims(dh, :geotile)

    # Extract run parameters and surface masks
    run_parameters_all = Altim.binned_filled_fileparts.(path2runs)
    surface_masks = unique(getindex.(run_parameters_all, :surface_mask))

    # Process geotiles for each surface mask
    geotiles = Dict()
    for surface_mask in surface_masks
        # Generate geotiles for this mask
        geotiles0 = Altim.geotiles_mask_hyps(surface_mask, geotile_width)
        
        # Make geotiles mutually exclusive by RGI region
        geotiles0, reg = Altim.geotiles_mutually_exclusive_rgi!(copy(geotiles0))
        
        # Align geotile indices with elevation change dimensions
        gt_ind = [findfirst(geotiles0.id .== g0) for g0 in collect(dgeotile)]
        
        # Rename area column to include surface mask
        rename!(geotiles0,"$(surface_mask)_area_km2" => "area_km2")
        
        # Store in dictionary
        geotiles[surface_mask] = geotiles0[gt_ind, :]
    end
end

# This section matches glaciers with geotiles and groups connected geotiles:
#
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
if !isfile(geotile_groups_fn) || force_remake
    # Find intersection between geotiles and glaciers
    geotiles0 = geotiles[surface_mask[1]] 
    
    # Initialize columns to store discharge indices and geotile intersections
    geotiles0[!,:dischage_ind] .= [Int64[]]
    geotiles0[!,:geotile_intersect] .= [Int64[]]

    # Load glacier outlines
    glaciers = GeoDataFrames.read(geomfile_rgi6)

    # Filter for large glaciers (>100 km2) that may cross tile boundaries
    min_area_km2 = 100;
    glaciers_large = glaciers[glaciers.Area .> min_area_km2, :]
    geometry_column = first(GI.geometrycolumns(glaciers_large))

    # Build spatial index tree for efficient intersection testing
    tree = STRtree(glaciers_large[:, geometry_column]; nodecapacity=3)
    
    # For each geotile, find intersecting glaciers
    @showprogress dt=1 desc="Match glaciers_large with geotiles ..." Threads.@threads for gt in eachrow(geotiles0)
        # Query tree for potential intersecting glaciers
        potential_idxs = SortTileRecursiveTree.query(tree, gt.extent)

        # Find glaciers that actually intersect the geotile
        intersecting = findall(Base.Fix1(GO.intersects, Altim.extent2rectangle(gt.extent)), 
                             view(glaciers_large[:, geometry_column], potential_idxs))
        
        gt.dischage_ind = potential_idxs[intersecting]
    end

    # Group geotiles that share glaciers
    @showprogress dt=1 desc="Group geotiles by glaciers_large ..." for (i, gt) in enumerate(eachrow(geotiles0))
        if .!isempty(gt.dischage_ind)
            # Find other geotiles that share glaciers with this one
            intersecting_geotiles = .!isdisjoint.(Ref(gt.dischage_ind), geotiles0.dischage_ind)
            intersecting_geotiles[i] = false
            
            if any(intersecting_geotiles)
                gt.geotile_intersect = findall(intersecting_geotiles)
            end
        end
    end

    # Identify connected groups of geotiles
    connectivity = geotiles0.geotile_intersect
    geotiles0[!,:group] = Altim.connected_groups(geotiles0.geotile_intersect)

    # Manual override for specific geotile groups
    # This handles cases where large glaciers cross multiple tiles but should be treated separately
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
        ["lat[-74-72]lon[-080-078]", "lat[-74-72]lon[-078-076]"]
    ]
    
    # Assign new group numbers to manual overrides
    group0 = maximum(geotiles0.group)
    for grp in geotile_groups_manual
        for g in grp
            geotiles0[findfirst(geotiles0.id .== g), :group] = group0
        end
        group0 += 1
    end

    # Convert geotile extents to rectangles and save results
    geometry_column = first(GI.geometrycolumns(geotiles0))
    geotiles0[!,geometry_column] = Altim.extent2rectangle.(geotiles0.extent)
    
    # Write geotile groups to file, excluding some columns
    GeoDataFrames.write(geotile_groups_fn,
                       geotiles0[:, Not(:extent, :area_km2, :dischage_ind, :geotile_intersect)]; 
                       crs=GFT.EPSG(4326))
end


# This section finds the optimal fit between GEMB (Glacier Energy and Mass Balance) data and altimetry derived volume change:
#
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
if !isfile(synthesized_gemb_fit) || force_remake
    # load gemb data
    filename_gemb_geotile_filled_dv = replace(filename_gemb_combined, ".jld2" => "_geotile_filled_dv.jld2")
    dv_gemb = Dict()
    smb, fac = load(filename_gemb_geotile_filled_dv, ("smb", "fac"))
    geotiles0 = GeoDataFrames.read(geotile_groups_fn)

    for binned_synthesized_file in path2runs
    #binned_synthesized_file = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_dh_best_cc_meanmadnorm3_v01_filled_ac_p1_synthesized.jld2"
        println(binned_synthesized_file)

        synthesized_gemb_fit = replace(binned_synthesized_file,  ".jld2" => "_gemb_fit.jld2")

        if !(isfile(synthesized_gemb_fit)) || force_remake
            t1 = time()

            run_parameters = Altim.binned_filled_fileparts(synthesized_gemb_fit)
            surface_mask = run_parameters.surface_mask;
            dh = FileIO.load(binned_synthesized_file, "dh_hyps")
    
            # convert elevation change to volume change
            dv_altim = Altim.dh2dv(dh, geotiles0);

            # sanity check
            if false
                gtids = ["lat[-14-12]lon[-076-074]", "lat[-10-08]lon[-078-076]", "lat[+28+30]lon[+082+084]", "lat[+58+60]lon[-136-134]", "lat[+68+70]lon[-070-068]", "lat[+62+64]lon[-044-042]", "lat[+32+34]lon[+076+078]", "lat[+28+30]lon[+092+094]", "lat[-46-44]lon[+168+170]", "lat[-70-68]lon[-068-066]", "lat[-52-50]lon[-074-072]", "lat[+58+60]lon[-136-134]"]

                gtids = ["lat[+58+60]lon[-136-134]"]
                for gtid in gtids
                    gemb_fit = Altim.gemb_bestfit_grouped(dv_altim, smb, fac, discharge, geotiles0; troubleshoot_geotile = gtid)
                end
            end

            if any(all(isnan.(dv_altim), dims = :date))
                dgeotile = dims(dv_altim, :geotile)
                all_nan_geotiles = dgeotile[findall(dropdims(all(isnan.(dv_altim), dims = :date), dims = :date))]
                printstyled("Geotiles that only contain NaNs:\n"; color = :red)
                println("$(collect(all_nan_geotiles)))")
                printstyled("Possible sources of error include:\n"; color = :red)
                printstyled("  [1] synthesis_error_file was run on a subset of files that do not include the full error range of the data, to fix this you need to delete the exisiting error file and any synthesis files that were create :\n"; color = :red)
                error("geotile volume change contains all NaNs")
            end

            # find optimal fit to gemb data
            # df = Altim.gemb_bestfit(dv_altim, smb, fac, discharge, geotiles[surface_mask])

            # there are issues with calibrating the smb model to individual geotiles glaciers 
            #can cross multiple geotiles therefore we calibrate the model for groups of 
            # distinct geotiles.

            df = Altim.gemb_bestfit_grouped(dv_altim, smb, fac, discharge, geotiles0)

            df[!,:area_km2] = sum.(geotiles0.area_km2)
            df.mad = df.mad ./ df.area_km2
            rename!(df, "mad"=>"mad_m")
            save(synthesized_gemb_fit, Dict("gemb_fit" => df))
            println("$binned_synthesized_file optimal GEMB fit found: $(round(Int,time() -t1))s")
        end
    end
end





# Downscale geotile-level data to individual glaciers by:
# 1. Loading the GEMB fit parameters and height change data for each geotile
# 2. Converting geotile height changes to per-glacier values
# 3. Scaling GEMB variables (SMB, FAC, runoff) using fit parameters and converting to per-glacier values


for binned_synthesized_file in path2runs
    #binned_synthesized_file = "/mnt/bylot-r3/data/binned/2deg/glacier_dh_best_cc_meanmadnorm3_v01_filled_ac_p2_synthesized.jld2"

    glacier_geotile_hyps_fn = replace(geomfile_rgi, ".gpkg" => "geotile_hyps.jld2")
    glaciers0 = FileIO.load(glacier_geotile_hyps_fn, "glaciers")
    gemb_fit = load(synthesized_gemb_fit, "gemb_fit")
    glaciers = copy(glaciers0)


  
    synthesized_gemb_fit = replace(binned_synthesized_file, ".jld2" => "_gemb_fit.jld2")
    perglacier_synthesized_file = replace(binned_synthesized_file, ".jld2" => "_perglacier.jld2")
    
    if !isfile(perglacier_synthesized_file) || force_remake
        ### DOWNSCALE TO EACH GLACIER BY HYPSOMETRY 
        # NOTE: you get strange spatial gradients in smb due to a lack of spatial gradents in 
        # climate forcing [like leeward side of mountain ranges]
        # therefore it is currently recommended to use the downscale by "area" method

        if downscale_to_glacier_method == "hypsometry"
            println("start downscaling from geotile to glacier, each iteration takes about 8 min so be patient")
            t1 = time()
           
            dh = load(binned_synthesized_file, "dh_hyps")
           
            glaciers = Altim.geotile2glacier!(glaciers, dh; varname = "dh");

            # load example geotile to get dimensions
            filename_gemb_geotile_filled = replace(filename_gemb_geotile, ".jld2" => "_$(gemb_fit.id[1])_filled.jld2")
            smb = load(filename_gemb_geotile_filled, "smb")
        
            dgeotile = dims(dh, :geotile)
            ddate = dims(smb, :date)
            dheight = dims(smb, :height)
            dpscale = dims(smb, :pscale)
            dΔheight = dims(smb, :Δheight)

            var = fill(NaN, dgeotile, ddate, dheight)

            # scale smb, fac, runoff
            for varname = ["smb", "fac", "runoff"]
                @showprogress desc="Populate calibrated $(varname)..." Threads.@threads for gt in eachrow(gemb_fit)
                    filename_gemb_geotile_filled = replace(filename_gemb_geotile, ".jld2" => "_$(gt.id)_filled.jld2")
                    var0 = load(filename_gemb_geotile_filled, varname)
                    var[At(gt.id), :, :] = var0[:,:,At(gt.pscale), At(gt.Δheight)]
                end

                glaciers = Altim.geotile2glacier!(glaciers, var; varname)
            end
            
            FileIO.save(perglacier_synthesized_file, Dict("glaciers" => glaciers))

            fn = splitpath(perglacier_synthesized_file)
            println("downscaling from geotile to glacier complete $(round(Int,time() -t1))s: $(fn[end])")
        end
    elseif downscale_to_glacier_method == "area"
        ### DOWNSCALE TO EACH GLACIER FRACTIONAL AREA WITHIN EACH GEOTILE
        glacier_area_km2 = sum.(glaciers.area_km2)

        @time for varname = ["dh", "smb", "fac", "runoff"]

            if varname == "dh"
                dh = load(binned_synthesized_file, "dh_hyps")
                var0 = Altim.dh2dv(dh, geotiles0);
            else
                var0 = load(filename_gemb_geotile_filled_dv, varname )
            end
            
            ddate = dims(var0, :date)
            glaciers[!,varname] =  [fill(NaN, ddate) for _ in 1:nrow(glaciers)]
            glaciers_grouped = groupby(glaciers, :geotile)

            for geotile in unique(glaciers.geotile)
            #geotile = "lat[+28+30]lon[+082+084]"
                gindex = findall(geotile .== glaciers.geotile)
                if varname != "dh"
                    gfit = gemb_fit[findfirst(isequal(geotile), gemb_fit.id),:]
                    gt_dv = var0[At(geotile), :, At(gfit.pscale), At(gfit.Δheight)]
                else
                    gt_dv = var0[At(geotile), :, :]
                end

                garea = sum.(glaciers[gindex, :area_km2])
                gweighting = garea ./ sum(garea)
                for (i,ig) in enumerate(gindex)
                    glaciers[ig,varname][:] = gt_dv * gweighting[i] ./ garea[i] * 1000 # convert from km3 to mwe
                end
            end
        end

        FileIO.save(perglacier_synthesized_file, Dict("glaciers" => glaciers))
    else
        error("downscale_to_glacier_method must be either \"hypsometry\" or \"area\"")
    end
end


#TODO: need to constuct proper file names for geospatial files. 

# This section processes each binned synthesized file in path2runs to:
# 1. Load glacier data from the per-glacier synthesized file
# 2. Add metadata for time series variables (dh, smb, fac, runoff)
# 3. Fit time series models to each variable for 2000-2023 period
# 4. Calculate glacier centroids and export rates to a FlatGeobuf file
# 5. Create a 2-degree grid and calculate area-weighted SMB means:
#    - Makes DataFrame with extent and geometry for each grid cell
#    - Gets glacier coordinates and areas
#    - For each grid cell:
#      * Finds glaciers within the cell
#      * Calculates area-weighted mean SMB trend
#    - Filters to cells with valid SMB values
#    - Exports gridded results to FlatGeobuf file
for binned_synthesized_file in path2runs
    #binned_synthesized_file = "/mnt/bylot-r3/data/binned/2deg/glacier_dh_best_cc_meanmadnorm3_v01_filled_ac_p2_synthesized.jld2"

    # for sanity checking in QGIS
    perglacier_synthesized_file = replace(binned_synthesized_file, ".jld2" => "_perglacier.jld2")
    glaciers = FileIO.load(perglacier_synthesized_file, "glaciers")

    ## this is a hack until arrow can support metadata
    tsvars= ["dh",  "smb", "fac", "runoff"]
    for tsvar in tsvars
        t = collect(dims(glaciers[1, tsvar], :date))
        # add back  metadata
        colmetadata!(glaciers, tsvar, "date", t, style=:note)
    end

    # this takes 5 min
    @time Altim.df_tsfit!(glaciers, [:dh, :smb, :fac, :runoff]; datelimits = (DateTime(2000,1,1), DateTime(2023,1,1)))

    outvars = ["geometry", "RGIId"]
    for tsvar in tsvars
        push!(outvars, "$(tsvar)_offset")
        push!(outvars, "$(tsvar)_trend")
        push!(outvars, "$(tsvar)_amplitude")
        push!(outvars, "$(tsvar)_phase")
    end

    # find glacier centroids
    perglacier_rates_path = replace(perglacier_synthesized_file, ".jld2" => "_rates.fgb")
    glaciers[!, :geom] = GI.Point.(GO.centroid.(glaciers.geom))
    valid = .!isnan.(getindex.(GI.coordinates.(glaciers.geom),1))
    rename!(glaciers, "geom" => "geometry")
    GeoDataFrames.write(perglacier_rates_path, glaciers[valid, outvars]; crs=GFT.EPSG(4326))
    # rsync devon:/mnt/bylot-r3/data/binned/2deg/glacier_dh_best_cc_meanmadnorm3_v01_filled_ac_p2_synthesized_perglacier_rates.fgb ~/data/Altim/project_data/

    # deteremine average rate per 1 degree grid
    Δdeg = 2
    lat = (-90+Δdeg/2):Δdeg:(90-Δdeg/2)
    lon = (-180+Δdeg/2):Δdeg:(180-Δdeg/2)
    source_crs1 = GFT.EPSG(4326)

    glaciers_grouped = DataFrame(extent=[extent = Extent(X=(x - Δdeg / 2, x + Δdeg / 2), Y=(y - Δdeg / 2, y + Δdeg / 2)) for y in lat, x in lon][:])
    glaciers_grouped[!, :geometry] .= [GI.Point(mean(gextent.X), mean(gextent.Y), crs=source_crs1) for gextent in glaciers_grouped.extent]

    # group for each lat/lon
    glacier_coords = GI.coordinates.(glaciers.geometry)
    glaciers_x = getindex.(glacier_coords,1)
    glaciers_y= getindex.(glacier_coords,2)

    glaciers_grouped[!, :smb] .= 0.0
    area = sum.(glaciers.area_km2)

    for gg in eachrow(glaciers_grouped)
    #gg = eachrow(glaciers_grouped)[1875]
        index = Altim.within.(Ref(gg.extent), glacier_x, glacier_y) .& (.!isnan.(glaciers.smb_trend))

        if any(index)
            # area weighted mean
            gg.smb = sum(glaciers[index, :smb_trend] .* area[index]) ./ sum(area[index])
        end
    end

    glaciers_grouped = glaciers_grouped[(glaciers_grouped.smb .!= 0.0) .& (.!isnan.(glaciers_grouped.smb)), :]
    GeoDataFrames.write(joinpath(paths.data_dir, "project_data","glaciers_grouped.fgb"), glaciers_grouped[:, [:geometry, :smb]])
end








binned_synthesized_file  = "/mnt/bylot-r3/data/binned/2deg/glacier_dh_best_meanmadnorm3_v01_filled_ac_p1_synthesized.jld2"
synthesized_gemb_fit = replace(binned_synthesized_file, ".jld2" => "_gemb_fit.jld2")
perglacier_synthesized_file = replace(binned_synthesized_file, ".jld2" => "_perglacier.jld2")
glaciers = FileIO.load(perglacier_synthesized_file, "glaciers")

# remove glaciers without area
zero_area_index  = sum.(glaciers.area_km2) .== 0
glaciers = glaciers[.!zero_area_index,:]

# remove glaciers where the elevation profile fell outside of valid range (only 12 glaciers)
glaciers = glaciers[.!map(x -> all(isnan.(x)), glaciers.dh),:]

# read in regions
rgi_regions = GeoDataFrames.read(paths.rgi6_regions_shp)

# read in GRACE data
grace = Altim.read_grace_rgi(; datadir=setpaths()[:grace_rgi])

# select a region
begin
# TODO: rgi 13/14/15 gemb and obs do not match.

if false
    geotile = "lat[+28+30]lon[+082+084]"# "lat[+56+58]lon[-134-132]"
    extent = geotiles["glacier"][findfirst(geotiles["glacier"].id .== geotile), :extent]
    region_geom = Altim.extent2rectangle(extent)
    in_region = GO.within.(tuple.(glaciers.CenLon, glaciers.CenLat), Ref(region_geom))
    in_region_dischage = GO.within.(tuple.(discharge.longitude, discharge.latitude), Ref(region_geom))
    
    title = geotile
else
    #rgi = 13
    rgi = 19
    in_region = falses(nrow(glaciers))
    in_region_dischage = falses(nrow(discharge))
    for region_geom = rgi_regions[rgi_regions.RGI_CODE.==rgi,:geometry]
        in_region .|= GO.within.(tuple.(glaciers.CenLon, glaciers.CenLat), Ref(region_geom))
        in_region_dischage .|= GO.within.(tuple.(discharge.longitude, discharge.latitude), Ref(region_geom))
    end
    title = Altim.rgi2label["rgi$rgi"]
end

# export for plotting
if false
    var_name = "runoff"
    v0 = copy(glaciers[in_region, :])
    start_date = DateTime(2000,04,1)
    end_date= DateTime(2023,04,1)

    for var_name in ["runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze", "dh"]
        sindex = findfirst(dims(v0[1, var_name], :date) .> start_date)
        eindex = findlast(dims(v0[1, var_name], :date) .< end_date)
        v1 = getindex.(v0[:, var_name],  eindex ) .- getindex.(v0[:, var_name],  sindex)
        v1 = v1.*sum.(glaciers[in_region, :area_km2])/1000
        v0[!, var_name] = v1
    end

    rename!(v0, "geom" => "geometry")
    #v0[!, :geometry] = tuple.(v0.CenLon,v0.CenLat)
    GeoDataFrames.write("glacier_2000_2023_anomaly.gpkg", v0[:, ["geometry", "RGIId", "runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze", "dh"]])
end

# calculate retional volume sums all sampled to the dh data
reg = Dict()
ddate = dims(glaciers[1, :dh], :date)
decyear = Altim.decimalyear.(ddate)
Δdecyear = decyear .- decyear[1]

for var_name in ["dh", "runoff", "fac", "smb"]; #, "rain", "acc", "melt", "ec", "refreeze", "dh"]
#var_name = "dh"
    
    v0 = glaciers[in_region , var_name]

    decyear0 = Altim.decimalyear.(dims(v0[1], :date))
    area_km2 = sum.(glaciers[in_region, :area_km2])
    v0 = v0 ./ 1000 .* area_km2
    v0 = vec(sum(reduce(hcat, v0), dims=2))
   
    notnan = .!isnan.(v0)
    decyear0 = decyear0[notnan]
    v0 = v0[notnan]
    index = (decyear .>= minimum(decyear0)) .& (decyear .<= maximum(decyear0))

    interp = DataInterpolations.LinearInterpolation(v0, decyear0)
    out = fill(NaN, ddate)
    out[index] = interp(decyear[index])

    reg[var_name] = out
end

reg["dischage"] = sum(discharge[in_region_dischage, :discharge_gtyr])

# add grace data
begin
    var_name = "grace"

    grace_date = Dim{:date}(Altim.datenum2date.(vec(grace["rgi$(rgi)"]["dM_gt_mdl_fill_date"])))
    grace_data = DimArray(vec(grace["rgi$(rgi)"]["dM_gt_mdl_fill"]), grace_date)

    grace_notnan = .!isnan.(grace_data)
    grace_decyear = Altim.decimalyear.(grace_date[grace_notnan])
    index = (decyear .>= minimum(grace_decyear)) .& (decyear .<= maximum(grace_decyear))

    grace_interp = fill(NaN, ddate)
    interp = DataInterpolations.LinearInterpolation(grace_data[grace_notnan], Altim.decimalyear.(grace_date[grace_notnan]))
    grace_interp[index] = interp(decyear[index])
    reg[var_name] = grace_interp
end



begin
    volume2mass = Altim.δice / 1000

    if false
        A = reg["dh"] .- reg["fac"]
        A_label = "this study"
        B = reg["grace"]
        B_label = "GRACE/-FO"
        y_label = "mass anomaly [Gt]"
    else
        A = reg["dh"] 
        A_label = "observed"
        B = reg["smb"] ./ volume2mass .+ reg["fac"] .- (reg["dischage"] .* Δdecyear) ./ volume2mass
        B_label = "modeled"
        y_label = "volume anomaly [km³]"
    end

    # align B to A
    index = .!isnan.(A) .& .!isnan.(B)
    B = B .- mean(B[index].-A[index])

    p = plot(A; title, label=A_label)
    plot!(B;title, label=B_label)
    ylabel!(p,y_label)
    xlims!((DateTime(2000, 1, 1), DateTime(2024, 1, 1)))
    display(p)
end

end
