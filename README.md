# Altim.jl

A Julia package for comprehensive glacier elevation change and hydrological analysis using multi-mission satellite altimetry and model data.

## Overview

Altim.jl provides an end-to-end workflow for processing satellite altimetry, glacier elevation change, and model data to quantify glacier mass changes and their impacts on river systems and populations. The package orchestrates a sequence of processing steps, from raw data ingestion to advanced analysis and visualization.

## Quick Start

### 1. Installation

```julia
using Pkg
Pkg.add("https://github.com/alex-s-gardner/Altim.jl")
```

### 2. Data Requirements

**Note:** This package is designed for use on JPL servers (bylot, devon, baffin) and expects data to be available at specific locations. If you are running elsewhere, you must adapt `src/local_paths.jl` to point to your data directories.

**You must download or have access to:**
- Satellite altimetry data (GEDI, ICESat, ICESat-2)
- Glacier elevation change stacks (Hugonnet et al.)
- Global DEMs (REMA, COP30, ArcticDEM, NASADEM)
- Glacier outlines (RGI v6/v7)
- Canopy height data (ETH Global Canopy Height)
- River network and basin data (MERIT Hydro, BasinATLAS, GRDC)
- Validation datasets (GRACE, Zemp2019, etc.)

Paths to these datasets are set in `src/local_paths.jl`. If you are not on a supported machine, you must edit this file to match your environment.

### 3. Running the Workflow

The main workflow is orchestrated by `src/master_run.jl`. This script sequentially runs all major processing steps, including data preparation, analysis, and output generation.

**To run the full pipeline:**
```bash
julia --project src/master_run.jl
```

This will:
1. Build satellite altimetry archives
2. Process glacier elevation change data
3. Validate and extract DEMs, masks, and canopy height
4. Perform hypsometric and statistical analyses
5. Synthesize data into global/regional datasets
6. Generate summary files and GIS outputs
7. Route runoff through river networks
8. Analyze glacier contributions to rivers and affected populations
9. Produce figures and regional results

**Processing times:**  
- GEDI: ~4 days  
- ICESat-2: ~1 week  
- ICESat: ~3 hours  
(Other steps are faster, but the full workflow may take several days depending on data and hardware.)

### 4. Script Descriptions

Each step in the workflow is implemented as a separate script in `src/` and included by `master_run.jl`. Advanced users can run these scripts individually for modular processing or debugging.

Key scripts (in order of execution):
- `geotile_build_archive.jl` — Build altimetry archives
- `geotile_build_hugonnet.jl` — Process Hugonnet glacier change data
- `geotile_ancillary_check.jl` — Validate ancillary data
- `geotile_dem_extract.jl` — Extract DEMs
- `geotile_mask_extract.jl` — Extract surface masks
- `geotile_canopyh_extract.jl` — Extract canopy height
- `geotile_hyps.jl` — Hypsometric analysis
- `gemb_classes_binning.jl` — Create glacier model classes
- `geotile_binning.jl` — Statistical binning
- `geotile_synthesis.jl` — Synthesize data
- `summary_files.jl` — Generate summary files
- `synthesis_plots_gis.jl` — Create plots and GIS outputs
- `land_surface_model_routing.jl` — Route land surface model runoff
- `glacier_routing.jl` — Route glacier runoff
- `gmax_global.jl` — Calculate maximum glacier river flux
- `river_buffer_population.jl` — Analyze population affected by glacier-fed rivers
- `gmax_point_figure.jl` — Generate gmax figures
- `regional_results.jl` — Produce regional results

### 5. Customization

- **Data paths:** Edit `src/local_paths.jl` to match your data storage locations if not running on a supported JPL server.
- **Partial runs:** You may comment/uncomment `include(...)` lines in `master_run.jl` to run only specific steps.
- **Checkpoints:** The workflow implements checkpoint logic to avoid redundant processing and allow efficient restarts.

### 6. Dependencies

Core dependencies (see `Project.toml` for full list):
- DimensionalData
- Rasters
- GeoInterface
- CairoMakie
- Statistics
- DataFrames
- FileIO
- NCDatasets
- Unitful

Install all dependencies with:
```julia
Pkg.instantiate()
```

## License

MIT License

## Contributing

Contributions are welcome! Please submit issues and pull requests on GitHub.


---

**For questions or help, please open an issue on GitHub.**
