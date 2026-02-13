# GlobalGlacierAnalysis.jl

A Julia package for comprehensive glacier elevation change and hydrological analysis using multi-mission satellite altimetry and model data.

## Overview

GlobalGlacierAnalysis.jl provides an end-to-end workflow for processing satellite altimetry, glacier elevation change, and model data to quantify glacier mass changes and their impacts on river systems and populations.

## Installation

```julia
using Pkg
Pkg.add("https://github.com/alex-s-gardner/GlobalGlacierAnalysis.jl")
```

## Quick Start

### 1. Data Requirements

**Note:** This package is designed for use on JPL servers (bylot, devon, baffin) and expects data at specific locations. If running elsewhere, adapt `src/local_paths.jl` to point to your data directories.

**Required datasets:**
- Satellite altimetry data (GEDI, ICESat, ICESat-2)
- Glacier elevation change stacks (Hugonnet et al.)
- Global DEMs (REMA, COP30, ArcticDEM, NASADEM)
- Glacier outlines (RGI v6/v7)
- Canopy height data (ETH Global Canopy Height)
- River network and basin data (MERIT Hydro, BasinATLAS, GRDC)
- Validation datasets (GRACE, Zemp2019, etc.)

### 2. Run the Workflow

The main workflow is orchestrated by `src/run_all.jl`:

```bash
julia --project src/run_all.jl
```

This executes the complete pipeline in the following sequence:

1. **Build satellite altimetry archives** - Processes GEDI, ICESat-2, ICESat data
2. **Process Hugonnet glacier elevation change data** - Extracts and validates glacier change datasets
3. **Validate existing ancillary data** - Checks existing DEMs, masks for geotiles
4. **Extract DEM data** - Processes elevation data for each geotile with slope/curvature calculations
5. **Extract mask data** - Processes glacier, land ice, and other surface masks
6. **Extract canopy height data** - Processes vegetation height information
7. **Perform hypsometric analysis** - Analyzes elevation distribution patterns
8. **Perform statistical binning** - Applies statistical methods to elevation change data
9. **Fill, extrapolate, and adjust binned data** - Completes missing data using interpolation methods
10. **Synthesize processed geotile data** - Combines multiple altimetry missions with error corrections
11. **Create glacier model (GEMB) classes** - Generates glacier energy and mass balance classifications
12. **Calculate global glacier discharge** - Computes ice discharge with filled Antarctic data using modeled SMB
13. **Calibrate GEMB model to altimetry data** - Optimizes glacier model parameters for grouped geotiles
14. **Calibrate GEMB model for each synthesized geotile dataset** - Applies model calibration to individual datasets
15. **Generate glacier-level summary files** - Creates key statistics and metrics for further analysis
16. **Export geotile-level glacier change trends** - Produces GIS-compatible files with trends and amplitudes
17. **Route land surface model runoff** - Processes terrestrial water fluxes through river networks
18. **Route glacier runoff** - Processes glacier meltwater through river networks
19. **Calculate gmax** - Determines maximum glacier contribution to river flux
20. **Generate point-based figures** - Creates gmax visualization plots
21. **Produce regional results** - Generates regional analysis outputs for sharing
22. **Generate extended data figures** - Produces manuscript-quality visualizations
23. **Analyze population affected by glacier-fed river changes** - Quantifies population impact using buffer analysis
24. **Create output files for the manuscript** - Generates submission-ready data and figures

**Processing times (on RAID @ 100 MB/s):**
- GEDI: ~4 days
- ICESat-2: ~1 week  
- ICESat: ~3 hours
- Full workflow: several days depending on data and hardware

### 3. Configuration

The workflow can be customized by modifying parameters at the top of `src/run_all.jl` and parameters within scripts called by `src/run_all.jl`


### 4. Checkpoint System

The workflow implements checkpoint logic to avoid redundant processing:
- Set `force_remake = false` for normal operation
- Use `force_remake_before` parameters for selective reprocessing
- Comment/uncomment specific steps in `run_all.jl` for partial runs

## Key Components

### Core Processing Scripts
- `run_all.jl` — Main workflow orchestrator
- `gemb_classes_binning.jl` — Glacier model classification
- `land_surface_model_routing.jl` — Terrestrial runoff routing
- `glacier_routing.jl` — Glacier runoff routing (called from `run_all.jl`)
- `gmax_global.jl` — Maximum glacier contribution analysis
- `gmax_point_figure.jl` — Visualization generation
- `regional_results.jl` — Regional analysis outputs
- `manuscript_extended_data_figures.jl` — Publication figures
- `river_buffer_population.jl` — Population impact analysis
- `create_output_files.jl` — Manuscript submission outputs

### Utility Modules

The package includes 15 utility modules that provide specialized functionality:

**Core Utilities:**
- `utilities.jl` — General utility functions and helpers
- `utilities_project.jl` — Project management and configuration utilities
- `utilities_main.jl` — Main workflow and orchestration functions

**Data Processing:**
- `utilities_build_archive.jl` — Satellite altimetry archive building functions
- `utilities_hugonnet.jl` — Hugonnet glacier elevation change data processing
- `utilities_readers.jl` — Data reading and input/output functions
- `utilities_binning.jl` — Statistical binning and data aggregation utilities
- `utilities_binning_lowlevel.jl` — Low-level binning algorithms and functions

**Analysis & Synthesis:**
- `utilities_gemb.jl` — Glacier Energy and Mass Balance (GEMB) model utilities
- `utilities_synthesis.jl` — Data synthesis and combination functions
- `utilities_routing.jl` — River routing and hydrological analysis utilities

**Output & Visualization:**
- `utilities_postprocessing.jl` — Post-processing and data export functions
- `utilities_plotting.jl` — Visualization and plotting utilities
- `utilities_manuscript.jl` — Manuscript figure and table utilities
- `utilities_response2reviewers.jl` — Response-to-reviewers analysis and outputs

**Geospatial:**
- `mapzonal.jl` — Zonal statistics and spatial analysis functions

## Dependencies

Core dependencies include:
- **Geospatial:** Proj, GeoArrays, SpaceLiDAR, Geodesy, Rasters, Shapefile, LibGEOS, GeoInterface
- **Data Processing:** DataFrames, Statistics, FileIO, Arrow, HTTP, MAT, CSV, NCDatasets, JLD2
- **Visualization:** CairoMakie, ColorSchemes
- **Scientific Computing:** NonlinearSolve, LsqFit, Distributions, LinearAlgebra
- **Custom Packages:** RangeExtractor, GeoTiles, BinStatistics

Install all dependencies with:
```julia
Pkg.instantiate()
```

## Data Outputs

The workflow produces:
- **Geotile-level data** - Processed elevation change and statistics
- **Global datasets** - Synthesized glacier mass change estimates
- **River routing results** - Glacier contribution to river systems
- **Population analysis** - Affected population statistics
- **GIS outputs** - Geospatial files for further analysis
- **Visualization figures** - Publication-ready plots and maps

## License

MIT License

## Contributing

Contributions are welcome! Please submit issues and pull requests on GitHub.

---

**For questions or help, please open an issue on GitHub.**
