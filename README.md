# Altim.jl

A Julia package for comprehensive glacier elevation change analysis using multi-mission satellite altimetry data.

## Overview

Altim.jl provides an end-to-end workflow for processing satellite altimetry data to understand glacier mass changes and their hydrological impacts. The package handles data from multiple satellite missions, combines them into a consistent framework, and analyzes glacier mass balance and meltwater routing.

## Installation

```julia
using Pkg
Pkg.add("https://github.com/alex-s-gardner/Altim.jl")
```

## Processing Pipeline

### 1. Building Altimetry Archive (geotile_build_archive.jl)

**Purpose:** Creates the foundational dataset by processing raw satellite measurements into a standardized format.

This step:
- Downloads and processes data from ICESat, ICESat-2, and GEDI missions
- Organizes measurements into 2-degree geotiles for efficient processing
- Filters data to focus on land ice areas
- Ensures consistent spatial and temporal organization of measurements

Processing times (100 MB/s IO):
- GEDI: ~4 days
- ICESat-2: ~1 week  
- ICESat: ~3 hours

### 2. Processing Hugonnet Data (geotile_build_hugonnet.jl)

**Purpose:** Incorporates independent glacier elevation change data as a validation source and gap-filler.

This step:
- Catalogs existing elevation change stacks
- Converts data into the same geotile format as satellite measurements
- Provides an independent dataset for cross-validation
- Helps fill gaps in satellite coverage

### 3. DEM Data Extraction (geotile_dem_extract.jl)

**Purpose:** Provides baseline elevation data and topographic parameters needed for accurate change detection.

Processes multiple DEMs:
- REMA (Antarctica)
- COP30 (global)
- ArcticDEM (Arctic)
- NASADEM (global)

Extracts critical parameters:
- Base elevations
- Slope measurements
- Surface curvature
- Topographic derivatives

### 4. Surface Mask Processing (geotile_mask_extract.jl)

**Purpose:** Identifies different surface types to ensure accurate data processing and interpretation.

Creates masks for:
- Floating ice (requires different processing)
- Glacier ice (primary analysis target)
- Inland water (affects routing)
- Land/ocean boundaries (affects discharge calculations)

### 5. Canopy Height Processing (geotile_canopyh_extract.jl)

**Purpose:** Accounts for vegetation effects on elevation measurements.

This step:
- Extracts vegetation height data from ETH Global Canopy Height dataset
- Helps correct elevation measurements in vegetated areas
- Improves accuracy of glacier boundary detection

### 6. Data Synthesis (geotile_synthesis.jl)

**Purpose:** Combines all data sources into a coherent analysis framework.

Key components (~5-6 hours total):
1. Glacier Hypsometry Processing (~37 min per inventory)
   - Calculates area-elevation distributions
   - Essential for mass balance calculations

2. Synthesis Error Calculation (~2 hrs)
   - Determines uncertainties between missions
   - Enables proper weighting of different data sources

3. Multi-mission Synthesis (~2 hrs)
   - Combines data from all missions
   - Accounts for different uncertainties and biases

4. Geotile to Glacier Downscaling
   - Converts gridded data to individual glacier measurements
   - Handles GEMB variables (SMB, FAC, runoff)

5. Time Series Analysis
   - Fits models to glacier changes (2000-2023)
   - Calculates rates and trends
   - Exports results in standard formats

### 7. Glacier Routing Analysis (glacier_routing.jl)

**Purpose:** Tracks how glacier meltwater moves through the landscape.

Analysis components:
1. Lowest Point Identification
   - Finds where meltwater enters river networks
   - Uses high-resolution 30m DEM

2. Drainage Network Mapping
   - Traces water flow paths
   - Identifies receiving rivers
   - Maps ocean vs inland termination

3. Variable Processing
   - Converts measurements to discharge units
   - Calculates monthly values
   - Handles unit conversions

4. Runoff Routing
   - Maps glacier contributions to rivers
   - Calculates flow rates
   - Determines relative glacier contributions

5. Discharge Analysis
   - Identifies terminal points
   - Maps ocean drainage
   - Incorporates measurement stations

## Dependencies

Core dependencies:
- DimensionalData: For labeled array operations
- Rasters: For spatial data handling
- GeoInterface: For geometric operations
- CairoMakie: For visualization
- Statistics: For numerical analysis
- DataFrames: For data organization
- FileIO: For file operations

## License

MIT License

## Contributing

Contributions welcome! Please submit issues and pull requests on GitHub.

## Citation

If you use Altim.jl in your research, please cite:
[Citation information to be added]
