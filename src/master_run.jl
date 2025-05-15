"""
    master_run.jl

Main workflow script for the GlobalGlacierAnalysis glacier and river analysis pipeline.

Orchestrates the entire data processing workflow by executing component scripts in sequence:
1. Builds satellite altimetry archives (GEDI, ICESat-2, ICESat)
2. Processes Hugonnet glacier elevation change data
3. Validates existing ancillary data (DEMs, masks) for geotiles
4. Extracts DEM data for each geotile
5. Extracts mask data (glacier, land ice) for each geotile
6. Extracts canopy height data for each geotile
7. Performs hypsometric analysis of elevation data
8. Creates glacier model (GEMB) classes for each geotile
9. Performs statistical binning of elevation data
10. Synthesizes processed geotile data into global/regional datasets
11. Generates summary files with key statistics and metrics
12. Creates visualization plots and GIS outputs
13. Routes land surface model runoff through river networks
14. Routes glacier runoff through river networks
15. Calculates gmax (maximum glacier contribution to river flux)
16. Analyzes population affected by glacier-fed river changes
17. Generates point-based figures for gmax visualization
18. Produces regional results for analysis and sharing

Implements checkpoint logic to avoid redundant processing, allowing efficient restarts.

Processing times (on old RAID @ 100 MB/s):
- GEDI: ~4 days
- ICESat-2: ~1 week  
- ICESat: ~3 hours
"""


# [1] Build archives from satellite altimetry data (GEDI, ICESat-2, ICESat)
include("geotile_build_archive.jl")

# [2] Process Hugonnet glacier elevation change data
include("geotile_build_hugonnet.jl")

# [3] Validate ancillary data that might already exist (DEMs, masks, etc.) for geotiles
include("geotile_ancillary_check.jl")

# [4] Extract DEM data for each geotile
include("geotile_dem_extract.jl")

# [5] Extract mask data (glacier, land ice, etc.) for each geotile
include("geotile_mask_extract.jl")

# [6] Extract canopy height data for each geotile
include("geotile_canopyh_extract.jl")

# [7] Perform hypsometric analysis of elevation data for each geotile
include("geotile_hyps.jl")

# [8] Create glacier model (GEMB) classes for each geotile
include("gemb_classes_binning.jl")

# [9] Perform statistical binning of elevation data for each geotile
include("geotile_binning.jl")

# [10] Synthesize processed geotile data into global and regional datasets
include("geotile_synthesis.jl")

# [11] Generate summary files with key statistics and metrics for further analysis and sharing
include("summary_files.jl")

# [12] Generate visualization plots and GIS outputs for data analysis and presentation
include("synthesis_plots_gis.jl")

# [13] Route land surface model runoff through river networks for river flux calculation
include("land_surface_model_routing.jl")

# [14] Route glacier runoff through river networks for glacier flux calculation
include("glacier_routing.jl")

# [15] Calculate gmax (maximum glacier contribution to river flux)
include("gmax_global.jl")

# [16] Analyze population affected by glacier-fed river changes using buffer analysis
include("river_buffer_population.jl")

# [17] Generate point-based figures for gmax visualization
include("gmax_point_figure.jl")

# [18] Generate regional results for further analysis and sharing
include("regional_results.jl")
