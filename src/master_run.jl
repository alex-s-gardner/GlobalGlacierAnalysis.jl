# consider running `geotile_remove_files.jl` to clear files from folder... if there is an 
# error this will allow you to restart processing where you left off, unlike 
# force_remake = true


# geotile_build_archive.jl processes altimetry data from ICESat, ICESat-2, and GEDI missions into geotiles.
#
# The script:
# 1. Sets up project configuration:
#    - Defines project ID and geotile width (2 degrees)
#    - Loads project paths and product definitions
#    - Filters geotiles to only include those with land ice
#
# 2. For each altimetry product (ICESat/ICESat-2/GEDI):
#    - Creates required directories for raw data and geotiles
#    - Searches for granules matching geotile extents using parallel 'find'
#    - Downloads granules from remote source using aria2c with retry logic
#    - Loads and processes granule data
#    - Sorts granules by longitude/latitude for improved processing speed
#    - Builds geotiles from the granule data
#
# Processing times (on old RAID @ 100 MB/s):
# - GEDI: ~4 days
# - ICESat-2: ~1 week  
# - ICESat: ~3 hours
include("/home/gardnera/Documents/GitHub/Altim.jl/src/geotile_build_archive.jl")


# geotile_build_hugonnet.jl processes Hugonnet glacier elevation change data into geotiles
#
# Processing steps:
# 1. Set up configuration:
#    - Define geotile width (2 degrees)
#    - Load project paths and create geotile grid
#    - Optional: Filter geotiles to specific regions/extents
#
# 2. Create output directory for Hugonnet geotiles if needed
#
# 3. Load Hugonnet data:
#    - Get catalogue of Hugonnet elevation change stacks
#    - Detect if using old/new data format
#
# 4. Build geotiles:
#    - Iterate through each geotile
#    - Process Hugonnet data into geotile format
#    - Save geotiles to output directory
include("/home/gardnera/Documents/GitHub/Altim.jl/src/geotile_build_hugonnet.jl")


# geotile_dem_extract.jl extracts elevation data and derivatives from global DEMs:
#
# Configuration:
# - Uses project ID v01 and 2-degree geotile width
# - Extracts slope and curvature in addition to elevation
# - Processes 4 DEMs: REMA v2, COP30 v2, ArcticDEM v4, NASADEM v1
#
# Processing:
# 1. Sets up project paths and loads product definitions
# 2. Filters geotiles to only include those with land ice
# 3. Extracts DEM data for each geotile using Altim.geotile_extract_dem()
#    - Includes elevation, slope and curvature
#    - Only processes if force_remake=true or output doesn't exist
include("/home/gardnera/Documents/GitHub/Altim.jl/src/geotile_dem_extract.jl")

# geotile_mask_extract.jl extracts surface type masks for geotiles:
#
# Configuration:
# - Uses project ID v01 and 2-degree geotile width
# - Filters to only include geotiles with land ice
# - Extracts masks for: floating ice, glacier ice, inland water, 
#   land, land ice, and ocean
#
# Processing:
# 1. Sets up project paths and loads product definitions
# 2. Filters geotiles to only include those with land ice
# 3. For each product:
#    - Extracts mask data using geotile_extract_mask()
#    - Only processes if force_remake=true or output doesn't exist
include("/home/gardnera/Documents/GitHub/Altim.jl/src/geotile_mask_extract.jl")

# geotile_canopyh_extract.jl extracts canopy height data for geotiles:
#
# Configuration:
# - Uses project ID v01 and 2-degree geotile width
# - Processes only geotiles containing land ice
# - Uses ETH Global Canopy Height 10m 2020 dataset
#
# Processing:
# 1. Sets up project paths and loads product definitions
# 2. Filters geotiles to only include those with land ice
# 3. Reads canopy height data from source GeoTIFF
# 4. Extracts canopy height values for each geotile using geotile_pointextract()
#    - Processes all altimetry missions
#    - Only processes if force_remake=true or output doesn't exist
#    - Uses nodata value of 255
include("/home/gardnera/Documents/GitHub/Altim.jl/src/geotile_canopyh_extract.jl")

# geotile_hyps.jl bins raw data into time-elevation geotile datacubes
#
# Configuration:
# - Uses project ID v01 and 2-degree geotile width
# - Processes only geotiles containing land ice
# - Uses ETH Global Canopy Height 10m 2020 dataset
#
# Processing:
# 1. Sets up project paths and loads product definitions
# 2. Filters geotiles to only include those with land ice
# 3. Reads canopy height data from source GeoTIFF
# 4. Extracts canopy height values for each geotile using geotile_pointextract()
#    - Processes all altimetry missions
#    - Only processes if force_remake=true or output doesn't exist
#    - Uses nodata value of 255
include("/home/gardnera/Documents/GitHub/Altim.jl/src/geotile_hyps.jl")

# gemb_classes_binning.jl processes GEMB (Glacier Energy Mass Balance) model data into geotiles and regions.
# Key processing steps:
#
# 1. Initialization (~1 min):
#    - Sets up project parameters (dates, heights, precipitation scales)
#    - Loads required packages and geotile definitions
#    - Configures input/output paths
#
# 2. GEMB Data Merging (~3 min):
#    - Combines multiple GEMB .mat files into single dataset
#    - Handles precipitation scaling and elevation delta variations
#    - Processes coordinates and saves merged data
#
# 3. Geotile Processing (~3.5 min):
#    - Bins GEMB data into time-height geotiles
#    - Buffers geotile extents to ensure coverage
#    - Tracks observation counts
#    - Handles multiple precipitation scales
#
# 4. Gap Filling and Height Classes (~9 hrs):
#    - Fills missing values across elevations using LOESS smoothing
#    - Creates new elevation classes through height shifts
#    - Applies physical constraints on variables
#    - Calculates volume changes using glacier areas
#    - Saves filled data for each geotile
#
# 5. Regional Aggregation (~11 sec):
#    - Aggregates geotile data into RGI regions
#    - Calculates regional glacier areas
#    - Sums volume changes by region
#    - Saves regional metrics
#
# Total runtime is approximately 9.5 hours, with most time spent on gap filling
# and creating height classes in step 4.
include("/home/gardnera/Documents/GitHub/Altim.jl/src/gemb_classes_binning.jl")

# geotile_binning.jl performs elevation data processing and analysis for glacier studies.
# Key processing steps:
#
# 1. Configuration (~1 min):
#    - Sets up project parameters (ID, geotile width, paths)
#    - Configures processing flags and thresholds
#    - Defines surface masks, DEM sources, and statistical methods
#    - Sets parameters for gap filling and mission alignment
#
# 2. Geotile Binning (~30 min):
#    - Bins elevation data into geotiles using Altim.geotile_binning()
#    - Processes multiple surface types (glacier, land etc)
#    - Applies different binning methods and corrections
#    - Handles filtered and unfiltered data streams
#
# 3. Gap Filling (~2 hrs):
#    - Fills data gaps using Altim.geotile_binned_fill()
#    - Applies multiple parameter sets for robust filling
#    - Handles amplitude corrections and mission normalization
#    - Generates diagnostic plots if requested
#
# 4. Mission Alignment (~1 hr):
#    - Aligns data between ICESat-2 and ICESat missions
#    - Corrects for land surface trends
#    - Replaces data with models in specific regions
#    - Ensures consistent measurements across missions
#
# 5. Regional Aggregation (~3 min):
#    - Calculates regional volume changes
#    - Aggregates geotile data to larger regions
#    - Processes multiple surface types and parameter sets
#
# Total runtime is approximately 3.5 hours, with most time spent on
# gap filling and mission alignment steps.
include("/home/gardnera/Documents/GitHub/Altim.jl/src/geotile_binning.jl")

# geotile_synthesis.jl performs synthesis and analysis of glacier elevation change data.
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
include("/home/gardnera/Documents/GitHub/Altim.jl/src/geotile_synthesis.jl")


# glacier_routing.jl analyzes glacier meltwater routing and discharge by:
#
# 1. Finding Lowest Points:
#    - Identifies lowest elevation point within each glacier boundary
#    - Uses 30m DEM to determine where meltwater enters river network
#
# 2. Mapping Drainage Networks:
#    - Maps glaciers to drainage basins and closest rivers
#    - Traces downstream river paths for each glacier
#    - Identifies rivers receiving glacier meltwater
#    - Determines ocean vs inland termination points
#
# 3. Processing Glacier Variables:
#    - Converts per-glacier geotile data to single values
#    - Interpolates variables (elevation change, mass balance, etc) to monthly values
#    - Converts units from m/year to m³/s where needed
#
# 4. Routing Glacier Runoff:
#    - Maps glacier runoff to river segments
#    - Calculates monthly mean runoff rates
#    - Computes fraction of river flow from glacier melt
#
# 5. Analyzing Discharge Points:
#    - Identifies terminal points where rivers/glaciers end
#    - Maps direct ocean drainage points
#    - Incorporates discharge measurement points
#    - Aggregates data into gridded format for visualization
#
# 6. Calculating Statistics:
#    - Converts runoff to water equivalent units
#    - Fits trends and seasonal cycles to time series
#    - Calculates spatial averages on lat/lon grid
#    - Produces final datasets for analysis and visualization
#
# The script processes global glacier and river datasets to understand:
# - How glacier meltwater enters and moves through river networks
# - Relative contributions of glacier melt to river discharge
# - Spatial and temporal patterns in glacier changes and runoff
include("/home/gardnera/Documents/GitHub/Altim.jl/src/glacier_routing.jl")