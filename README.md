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

The main workflow is orchestrated by `src/master_run.jl`:

```bash
julia --project src/master_run.jl
```

This executes the complete pipeline:
1. Build satellite altimetry archives
2. Process glacier elevation change data
3. Extract and validate DEMs, masks, and canopy height
4. Perform hypsometric and statistical analyses
5. Synthesize data into global/regional datasets
6. Route runoff through river networks
7. Analyze glacier contributions to rivers and affected populations
8. Generate figures and regional results

**Processing times:**
- GEDI: ~4 days
- ICESat-2: ~1 week
- ICESat: ~3 hours
- Full workflow: several days depending on data and hardware

### 3. Customization

- **Data paths:** Edit `src/local_paths.jl` to match your data storage locations
- **Partial runs:** Comment/uncomment `include(...)` lines in `master_run.jl` to run specific steps
- **Checkpoints:** The workflow implements checkpoint logic for efficient restarts

## Key Scripts

Each processing step is implemented as a separate script in `src/`:

- `geotile_build_archive.jl` — Build altimetry archives
- `geotile_build_hugonnet.jl` — Process Hugonnet glacier change data
- `geotile_hyps.jl` — Hypsometric analysis
- `geotile_binning.jl` — Statistical binning
- `geotile_synthesis.jl` — Synthesize data
- `glacier_routing.jl` — Route glacier runoff
- `river_buffer_population.jl` — Analyze population affected by glacier-fed rivers
- `regional_results.jl` — Produce regional results

## Dependencies

Core dependencies include:
- DimensionalData, Rasters, GeoInterface
- CairoMakie, DataFrames, Statistics
- FileIO, NCDatasets

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
