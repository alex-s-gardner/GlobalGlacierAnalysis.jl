# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

GlobalGlacierAnalysis.jl is a comprehensive Julia package for analyzing global glacier mass changes using multi-mission satellite altimetry (ICESat, ICESat-2, GEDI) combined with Hugonnet et al. glacier elevation change stacks. The workflow processes raw satellite data through to publication-ready figures and datasets, including routing glacier runoff through river networks and quantifying population impacts.

**Citation:** Gardner, A. S., Schlegel, N.-J., Greene, C. A., Hugonnet, R., Menounos, B., Wiese, D. N., & Berthier, E. (in review). Glacier contributions to rivers and oceans in the early twenty first century.

## System Requirements

- **Julia version:** 1.12+ (currently using 1.12.6)
- **Platform:** Designed for JPL servers (bylot, devon, baffin) with access to shared data directories
- **Data location:** Expects data at `/mnt/bylot-r3/data/` and `/mnt/devon-r2/shared_data/`
  - If running elsewhere, adapt `src/local_paths.jl` to point to your data locations
- **Storage:** Several TB of disk space for raw altimetry and processed geotiles

## Quick Start

```bash
# Activate the project environment
julia --project

# Install/update dependencies
julia --project -e 'using Pkg; Pkg.instantiate()'

# Run the complete workflow (takes several days)
julia --project src/run_all.jl

# Or run individual processing steps (see src/run_all.jl for available scripts)
julia --project src/regional_results.jl
```

## Development Commands

### Running Individual Scripts
```bash
# Run regional analysis
julia --project src/regional_results.jl

# Generate manuscript figures
julia --project src/manuscript_extended_data_figures.jl

# Run routing analysis
julia --project src/glacier_routing.jl
julia --project src/land_surface_model_routing.jl

# Create final outputs
julia --project src/create_output_files.jl
```

### Interactive Development
```bash
# Start Julia REPL with project
julia --project

# Load the module
julia> import GlobalGlacierAnalysis as GGA
julia> using Unitful
julia> Unitful.register(GGA.MyUnits)

# Access functions
julia> GGA.geotile_build_archive(...)
```

### Testing Single Geotiles
Edit `src/run_all.jl` and set:
```julia
single_geotile_test = "lat+45+47lon-123-121"  # Test one geotile
# or
single_geotile_test = GGA.geotiles_golden_test[1]  # Use predefined test
```

### Reprocessing Data
```julia
# Selective reprocessing by date
force_remake_before = DateTime(2025, 7, 1)

# Force remake all (use with caution)
force_remake = true

# Or delete specific output files and re-run
```

## Key Architecture Concepts

### Geotile-Based Processing

The core spatial unit is a "geotile" - a geographic grid cell (typically 2° × 2°) identified by strings like `"lat+45+47lon-123-121"`. All satellite altimetry and analysis results are organized by geotiles, enabling:
- Parallel processing of independent geographic regions
- Efficient spatial subsetting and data access
- Hierarchical aggregation from local to regional to global scales

Geotile IDs encode their bounding boxes: `lat[min_lat][max_lat]lon[min_lon][max_lon]` with latitude [-90, 90] and longitude [-180, 180].

### Multi-Mission Data Synthesis

The workflow combines four elevation change sources:
1. **ICESat-2** (ATL06 v6): Dense laser altimetry, 2018-present
2. **ICESat** (GLAH06 v34): Laser altimetry, 2003-2009
3. **GEDI** (GEDI02_A v2): Lidar, 2019-present, ±52° latitude
4. **Hugonnet** (HSTACK v1): Optical stereo DEMs, 2000-2020

Each mission has unique spatial coverage, temporal sampling, and error characteristics. The synthesis process:
- Co-registers missions to a common reference frame
- Applies amplitude and curvature corrections
- Fills gaps using spatial/temporal interpolation
- Propagates uncertainties through the ensemble

### Three-Stage Data Flow

1. **Archive Building** (`utilities_build_archive.jl`)
   - Downloads/ingests raw satellite granules
   - Extracts elevation points with geolocation
   - Organizes by geotile for efficient access

2. **Statistical Binning** (`utilities_binning.jl`, `utilities_binning_lowlevel.jl`)
   - Groups elevation changes by surface type (glacier, land ice)
   - Applies robust statistics (median, NMAD) to filter outliers
   - Extracts DEMs, masks, and canopy height for each altimetry point
   - Creates hypsometric (elevation-binned) summaries

3. **Synthesis & Analysis** (`utilities_synthesis.jl`, `utilities_postprocessing.jl`)
   - Fills missing data using spatial/temporal interpolation (4 fill parameters)
   - Calibrates GEMB (Glacier Energy and Mass Balance) model to observations
   - Routes glacier meltwater through river networks
   - Generates regional summaries and figures

### Configuration System

Processing is controlled by parameters in `src/run_all.jl`:
- `project_id = :v01` - version identifier for output paths
- `geotile_width = 2` - degrees (use 2 for full runs, smaller for testing)
- `missions = (:icesat2, :icesat, :gedi, :hugonnet)` - which datasets to process
- `force_remake = false` - checkpoint system to skip completed steps
- `surface_masks = [:glacier, :glacier_rgi7]` - RGI versions to analyze
- `binning_methods = ["nmad3", "nmad5", "median"]` - outlier filtering approaches
- `fill_params = [1, 2, 3, 4]` - gap-filling aggressiveness levels

The `binned_filled_filepaths()` function generates Cartesian products of processing options to create ensemble members.

## File Organization

```
src/
├── run_all.jl                          # Main orchestrator script
├── GlobalGlacierAnalysis.jl            # Module definition and exports
├── local_paths.jl                      # Machine-specific data paths
│
├── utilities_project.jl                # Project config (products, paths, constants)
├── utilities_main.jl                   # High-level workflow functions
├── utilities.jl                        # General helpers and utilities
│
├── utilities_build_archive.jl          # Satellite data ingestion
├── utilities_hugonnet.jl               # Hugonnet stack processing
├── utilities_binning.jl                # Statistical binning (high-level)
├── utilities_binning_lowlevel.jl       # Binning algorithms (low-level)
├── utilities_synthesis.jl              # Multi-mission data synthesis
├── utilities_gemb.jl                   # GEMB model calibration
├── utilities_routing.jl                # River network routing
├── utilities_postprocessing.jl         # Export and aggregation
├── utilities_plotting.jl               # Visualization functions
├── utilities_manuscript.jl             # Publication figures
├── utilities_readers.jl                # Data I/O functions
├── utilities_response2reviewers.jl     # Review response analysis
│
├── gemb_classes_binning.jl             # GEMB classification workflow
├── glacier_routing.jl                  # Glacier runoff routing
├── land_surface_model_routing.jl       # LSM runoff routing
├── gmax_global.jl                      # Maximum glacier contribution
├── gmax_point_figure.jl                # Gmax visualization
├── regional_results.jl                 # Regional summaries
├── river_buffer_population.jl          # Population impact analysis
├── manuscript_extended_data_figures.jl # Extended data figures
├── create_output_files.jl              # Final NetCDF outputs
└── mapzonal.jl                         # Zonal statistics functions
```

Output data is organized in `pathlocal.data_dir` with this structure:
```
data_dir/
├── icesat2/ATL06/006/geotile/2deg/     # ICESat-2 geotiles
├── icesat/GLAH06/034/geotile/2deg/     # ICESat geotiles
├── gedi/GEDI02_A/002/geotile/2deg/     # GEDI geotiles
├── hugonnet/geotile/2deg/              # Hugonnet geotiles
├── binned/2deg/                        # Binned elevation changes
├── binned_unfiltered/2deg/             # Binned without outlier filtering
└── filled/2deg/                        # Gap-filled synthesis products
```

## Important Constants and Settings

From `utilities_project.jl`:
- `δice = 910` - glacier ice density [kg/m³]
- `ocean_area_km2 = 362.5 * 1E6` - for global flux calculations
- `reference_ensemble_file` - the "best" synthesis used as reference for error calculations

Error model parameters:
- `error_quantile = 0.95` - percentile for uncertainty estimates
- `error_scaling = 1.5` - conservative scaling factor for 2σ errors

## Data Dependencies

**Required datasets** (see README.md and `local_paths.jl` for sources):
- Satellite altimetry: GEDI (4 days to process), ICESat-2 (1 week), ICESat (3 hours)
- Hugonnet glacier elevation change stacks (request from authors)
- Global DEMs: REMA v2, ArcticDEM v4, COP30 v2, NASADEM v1
- Glacier outlines: RGI v6/v7 from GLIMS
- ETH Global Canopy Height 10m 2020 (download via aria2c from ETH)
- River networks: MERIT Hydro, BasinATLAS, GRDC
- Validation: GRACE, Zemp2019 glacier mass balance

## Common Development Patterns

### Running Partial Workflows

The workflow uses checkpoint logic - set `force_remake = false` and it will skip steps where output files already exist. To rerun specific steps:
1. Delete output files for those steps, or
2. Set `force_remake_before = DateTime(2025, 7, 1)` to reprocess files older than that date, or
3. Comment out unneeded steps in `run_all.jl`

### Testing on Single Geotiles

For rapid iteration, set `single_geotile_test = "lat+45+47lon-123-121"` in `run_all.jl` to process only one geotile. Use `GGA.geotiles_golden_test` for a predefined test set.

### Working with DimArrays and Rasters

The package heavily uses DimensionalData.jl and Rasters.jl for labeled arrays:
- Dimensions: `Ti` (time), `X` (longitude), `Y` (latitude), `:geotile`, `:varname`, `:error`
- Access via named dimensions: `data[varname=At("runoff"), geotile=At(gt_id), error=At(false)]`
- Units via Unitful.jl: most elevation values have units attached (e.g., `u"mm"`)

### Parallel Processing Patterns

The workflow is designed for parallel execution at the geotile level but typically runs serially due to I/O bottlenecks. When modifying, maintain geotile independence - each should process without needing data from other geotiles.

### Plotting and Visualization

Use `CairoMakie.jl` for all figures. Set `plots_show = true` and `plots_save = true` in scripts to enable visualization during development. Figures save to `pathlocal.figures`.

### Error Handling and Debugging

When debugging issues:
1. Check file modification times match expectations (checkpoint system compares these)
2. Enable warnings: `warnings = true` in binning parameters
3. Use `@info` or `println` for progress tracking - the workflow has minimal logging by default
4. Check intermediate `.jld2` files with `JLD2.load(filename)` to inspect data
5. Look for `missing` values in DimArrays - these indicate gaps in processing
6. Verify geotile IDs are formatted correctly: `"lat[min][max]lon[min][max]"` with signs

### Ensemble Processing

The workflow generates many ensemble members through parameter combinations:
- Different surface masks (`:glacier`, `:glacier_rgi7`)
- Different DEMs (`:best`, `:cop30_v2`)
- Curvature correction (true/false)
- Amplitude correction (true/false)
- Binning methods (`"median"`, `"nmad3"`, `"nmad5"`)
- Fill parameters (1, 2, 3, 4)
- Binned vs unfiltered data

Use `binned_filled_filepaths()` to generate paths to all ensemble members. The reference ensemble is specified in `utilities_project.jl` as `reference_ensemble_file`.

## Module Structure

The main module `GlobalGlacierAnalysis` (abbreviated `GGA` in scripts) includes 15 utility files that export hundreds of functions. Key functions by module:

**utilities_main.jl:** `geotile_build_archive`, `geotile_build_hugonnet`, `geotile_dem_extract`, `geotile_mask_extract`, `geotile_hyps_extract`

**utilities_binning.jl:** `geotile_binning`, `geotile_filled_extrapolated`, `geotiles_mean_error`

**utilities_synthesis.jl:** `synthesize_geotile_data`, `gemb_fit_to_altimetry`

**utilities_routing.jl:** `route_glacier_runoff`, `route_lsm_runoff`, `calculate_gmax`

**utilities_postprocessing.jl:** `create_regional_summaries`, `export_geotile_trends`, `dimstack2netcdf`, `netcdf2dimstack`

## Testing

`test/runtests.jl` runs the suite (~2350 tests, ~4 minutes). It uses synthetic data with known
ground truth and mocks external files, so it needs no access to the JPL data directories:

```bash
julia --project test/runtests.jl
```

Tests are organized as `test/unit/` (pure functions), `test/algorithms/` (binning, GEMB, routing,
synthesis with synthetic inputs), `test/io/` (round-trips), and `test/integration/` (end-to-end
workflows and statistical invariants). Generators live in `test/fixtures/`. See `test/README.md`.

Scientific validation is separate from the test suite and happens through:
- Comparison with published datasets (GRACE, Zemp2019)
- Visual inspection of figures generated in `manuscript_extended_data_figures.jl`
- Running single-geotile tests before full workflow execution
- Checkpoint system that validates output files exist and are readable

When adding new functionality:
1. Test on a single geotile first using `single_geotile_test`
2. Verify outputs visually with plots (`plots_show = true`)
3. Check intermediate files are created correctly in expected output directories
4. Run on the golden test set: `GGA.geotiles_golden_test`

## Known Gotchas

- **Path assumptions:** Code assumes JPL server environment. Running elsewhere requires editing `local_paths.jl`
- **Memory usage:** Processing all missions simultaneously can consume >64GB RAM for large geotiles
- **Float precision:** Some operations mix Float32 and Float64 - be careful with type stability in hot loops
- **Coordinate conventions:** Longitude is [-180, 180], latitude is [-90, 90]
- **Time zones:** Timestamps are UTC; convert with `local2utc = Hour(7)` constant
- **File formats:** JLD2 for intermediate data, NetCDF for final outputs, Arrow for tabular data
- **Custom units:** The package defines `Gt` (gigatons) via `MyUnits` module - must call `Unitful.register(GGA.MyUnits)`
- **Two kinds of validation:** `test/runtests.jl` covers algorithms with synthetic data; scientific correctness is established separately by comparison against published datasets
- **Geotile independence:** Each geotile processes independently - don't create cross-geotile dependencies when modifying code
