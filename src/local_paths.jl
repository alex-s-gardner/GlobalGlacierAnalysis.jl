"""
    setpaths()

Returns a named tuple containing local file paths for data, DEMs, masks, and other resources
based on the current machine's hostname. Currently configured for JPL servers (bylot, devon, baffin).

The paths include:
- Data directories for altimetry and project data
- Figure output directories
- Global DEM VRTs (Copernicus, NASADEM, ArcticDEM, REMA)
- Vector files (glacier outlines, ice masks, etc.)
- Validation datasets (GRACE, glacier mass balance studies)
- River network data

Throws an error if run on an unrecognized machine.

NOTE:
If the GlobalGlacierAnalysis.jlacierAnalysis.jlacierAnalysis.jl package is installed from a repository, the local paths will need to be set 
manually and all dependent files will need to be downloaded and stored in the local paths.
"""
function setpaths() 
    hostname = gethostname()
    project_dir = "/mnt/bylot-r3/data/project_data/"
    data_dir = "/mnt/bylot-r3/data/"

    if (hostname == "bylot.jpl.nasa.gov") || (hostname == "devon.jpl.nasa.gov") || (hostname == "baffin.jpl.nasa.gov")

        pathlocal = (
            
            # all altimetry and geotile data will be stored in a stucture within data_dir
            data_dir = data_dir,

            # folder where random project data is stored
            project_dir=project_dir,

            # figure output directory
            figures = "/mnt/bylot-r3/altim_figs/",

            # comp30 DEM - see cop30_vrt_build.jl for more details
            cop30_v1 = "/mnt/devon-r2/shared_data/COP-DEM_GLO-30-DGED_PUBLIC/GLO30_DGED_hgt.vrt",
            cop30_v2 = "/mnt/devon-r2/shared_data/copernicus-dem-30m/DEM.vrt",

            # NASADEM - see nasadem_vrt_build.jl for more details
            nasadem_v1 = "/mnt/devon-r2/shared_data/NASADEM/mosaic/NASADEM_hgt.vrt",

            # ArcticDEM - see https://www.pgc.umn.edu/data/arcticdem/
            arcticdem_v3_10m = "/mnt/devon-r2/shared_data/arcticdem/mosaics/v3.0/10m/arcticdem_v3.0_10m.vrt",
            arcticdem_v4_10m = "/mnt/devon-r2/shared_data/arcticdem/mosaics/v4.1/10m/arcticdem_v4.1_10m.vrt",

            # REMA - see https://www.pgc.umn.edu/data/rema/
            rema_v2_10m = "/mnt/devon-r2/shared_data/rema/mosaics/v2.0/10m/rema_v2.0_10m.vrt",

            # ASTER and WorldView elevation data: must be requested from the author: Romain Hugonnet <hugonnet@uw.edu>
            hugonnet_v1_stacks = "/mnt/bylot-r3/data/hugonnet/",

            # Global Canopy Height data: must be downloaded from ETH
            # Download canopy height data in teminal
            # aria2c -x10 https://share.phys.ethz.ch/~pf/nlangdata/ETH_GlobalCanopyHeight_10m_2020_version1/ETH_GlobalCanopyHeight_10m_2020_version1.tar.gz 
            # tar -xvzf ETH_GlobalCanopyHeight_10m_2020_version1.tar.gz
            canopyheight_10m_v1 = "/mnt/devon-r0/shared_data/canopy_height/ETH_GlobalCanopyHeight_10m_2020_version1/ETH_GlobalCanopyHeight_10m_2020_mosaic_Map.vrt",

            # Geoid files
            # these should automatically download from https://www.agisoft.com/downloads/geoids/
            # when the utilities.jl/geoid() function is called
            geoid_dir = "/mnt/devon-r2/shared_data/geoids",

            # its_live parameter files: can be downloaded from https://its-live-data.s3.amazonaws.com/index.html#autorift_parameters/v001/
            itslive_parameters = "/mnt/devon-r2/data/its-live-data/autorift_parameters/v001/",   

            # Landice mask: This is a custom mask created by Gardner et al. 2025, please request from the authors: Alex Gardner <alex.s.gardner@jpl.nasa.gov>
            icemask = "/mnt/devon2-r1/devon0/gardnera/Documents/GitHub/GlobalGlacierAnalysis/data/land_ice_0.0083_cog.tif",

            # RGI outlines: can be downloaded from https://www.glims.org/RGI/
            # RGI 7 region outlines
            RGI_regions = "/mnt/devon-r2/data/GlacierOutlines/RGI2000-v7.0-regions",
            # RGI 7 meged with only subset of Greenland (cl0cl1) glaciers selected
            glacier_shp="/mnt/bylot-r3/data/vector_files/glacier_cl0cl1.shp",
            # RGI 7 merged outlines
            glacier_rgi7_shp="/mnt/bylot-r3/data/vector_files/RGI2000-v7.0-C-01_global.shp",
            # RGI 6 merged outlines buffered outward by 1km [Optional]
            glacier_b1km_shp="/mnt/bylot-r3/data/vector_files/glacier_b1km.shp",
            # RGI 6 merged outlines buffered outward by 10km [Optional]
            glacier_b10km_shp="/mnt/bylot-r3/data/vector_files/glacier_b10km.shp",
            # RGI 6 individual outlines combined into a geopackage
            glacier_individual="/mnt/bylot-r3/data/GlacierOutlines/rgi60/rgi60_Global.gpkg",
            # RGI 7 individual outlines combined into a geopackage
            glacier_rgi7_individual="/mnt/bylot-r3/data/GlacierOutlines/RGI2000-v7.0-G-global-fix/rgi70_Global.gpkg",
            # landice mask: this is a custom mask created by Gardner et al. 2025, please request from the authors: Alex Gardner <alex.s.gardner@jpl.nasa.gov>
            landice_shp="/mnt/bylot-r3/data/vector_files/land_ice.shp",
            # landice mask buffered outward by 1km: this is a custom mask created by Gardner et al. 2025, please request from the authors: Alex Gardner <alex.s.gardner@jpl.nasa.gov>
            landice_b1km_shp="/mnt/bylot-r3/data/vector_files/land_ice_b1km.shp",
            # landice mask buffered outward by 10km: this is a custom mask created by Gardner et al. 2025, please request from the authors: Alex Gardner <alex.s.gardner@jpl.nasa.gov>
            landice_b10km_shp="/mnt/bylot-r3/data/vector_files/land_ice_b10km.shp",
            # floating ice mask: this is a custom mask created by Gardner et al. 2025, please request from the authors: Alex Gardner <alex.s.gardner@jpl.nasa.gov>
            floating_shp="/mnt/bylot-r3/data/vector_files/floating_ice.shp",
            # landice mask: this is a custom mask created by Gardner et al. 2025, please request from the authors: Alex Gardner <alex.s.gardner@jpl.nasa.gov>
            not_ocean_shp="/mnt/bylot-r3/data/vector_files/not_ocean.shp", # GGSHHS_f_L1 +  scripps_antarctica_grounded
            # RGI 6 regions file
            rgi6_regions_shp="/mnt/bylot-r3/data/vector_files/00_rgi60_O1Regions_fix.shp",
            # water polygons from Global Self-consistent, Hierarchical, High-resolution Shorelines: https://www.soest.hawaii.edu/pwessel/gshhg/
            water_shp="/mnt/bylot-r3/data/vector_files/GSHHS_f_L2.shp", #GSHHG lakes
            # ESRI Continents: https://hub.arcgis.com/datasets/esri::world-continents/about
            continents="/mnt/bylot-r3/data/vector_files/World_Continents.geojson",
            # Natural Earth Countries: https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-0-countries/
            countries="/mnt/bylot-r3/data/vector_files/ne_10m_admin_0_countries.geojson",
            # grace timeseries data for validation [created in MATLAB using jplMasconTS.m]
            # this is a custom mask created by Gardner et al. 2025, please request from the authors: Alex Gardner <alex.s.gardner@jpl.nasa.gov>
            grace_rgi = "/mnt/bylot-r3/data/GRACE/dM_grace_rgi.mat",

            # glacier mass balance data for comparision [all optional and from original papers]
            zemp_2019 = "/mnt/bylot-r3/data/glacier_mb/Zemp2019/Zemp_etal_DataTables2a-t_results_regions_global/",
            marzeion_2012 = "/mnt/bylot-r3/data/glacier_mb/Marzeion2012/Marzeion2012.mat",
            marzeion_2020 = "/mnt/bylot-r3/data/glacier_mb/Marzeion2020/suppl_GlacierMIP_results.nc",
            hock_2019 = "/mnt/bylot-r3/data/glacier_mb/Hock2019/S0022143019000224sup001.nc",
            ipcc_ar6 = "/mnt/bylot-r3/data/glacier_mb/IPCCAR6/",

            # glacier runoff data for comparision: https://datadryad.org/dataset/doi:10.5061/dryad.pk0p2ngxf
            wimberly_2024 = "/mnt/bylot-r3/data/glacier_mb/Wimberly2024/all_rf_aligned_data.csv",

            #GlaMBIE annula mass balance data: https://wgms.ch/data_glambie/
            glambie_2024 = "/mnt/bylot-r3/data/glacier_mb/glambie/41586_2024_8545_MOESM4_ESM_Gt.csv",

            # path to river files
            # River reaches: https://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_Hydro/
            river = "/mnt/bylot-r3/data/rivers/MERIT_Hydro_v07_Basins_v01",
            # River basins: https://www.reachhydro.org/home/params/merit-basins
            river_basins = "/mnt/bylot-r3/data/rivers/BasinATLAS_Data_v10.gdb/BasinATLAS_v10_lev02.geojson",
            # Major river basins: https://grdc.bafg.de/products/basin_layers/major_rivers/
            river_major_basins = "/mnt/bylot-r3/data/rivers/GRDC_Major_River_Basins/mrb_basins.json",

            # Kochtitzky NH discharge and terminus retreate 
            discharge_nh="/mnt/bylot-r3/data/GlacierOutlines/GlacierDischarge/Kochtitzky2022/41467_2022_33231_MOESM4_ESM.csv",
            discharge_npi="/mnt/bylot-r3/data/GlacierOutlines/GlacierDischarge/Fuerst2023/fuerst_2023_npi_comparison_ice_discharge_v1.0.0.txt",
            discharge_spi="/mnt/bylot-r3/data/GlacierOutlines/GlacierDischarge/Fuerst2023/fuerst_2023_spi_comparison_ice_discharge_v1.0.0.txt",
            rgi6_southern_andes = "/mnt/bylot-r3/data/GlacierOutlines/rgi60/17_rgi60_SouthernAndes/17_rgi60_SouthernAndes.shp",

            # discharge data for global discharge analysis [created as part of processing pipeline]
            discharge_global = "/mnt/bylot-r3/data/GlacierOutlines/GlacierDischarge/global_glacier_discharge.jld2",

            glacier_summary = joinpath(project_dir, "gardner2025_glacier_summary.nc")
        )

        return pathlocal
    else
        error("unrecognized machine, Atlim.jl requires local paths to be set in setpaths() in order to run")
    end
end