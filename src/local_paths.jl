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
"""
function setpaths() 
    hostname = gethostname()
    if (hostname == "bylot.jpl.nasa.gov") || (hostname == "devon.jpl.nasa.gov") || (hostname == "baffin.jpl.nasa.gov")

        pathlocal = (
            # all altimetry and geotile data will be stored in a stucture within data_dir
            data_dir = "/mnt/bylot-r3/data/",
            project_dir = "/mnt/bylot-r3/data/project_data/",

            # figure output directory
            figures = "/mnt/bylot-r3/altim_figs/",

            # location of global DEM vrts
            cop30_v1 = "/mnt/devon-r2/shared_data/COP-DEM_GLO-30-DGED_PUBLIC/GLO30_DGED_hgt.vrt",
            cop30_v2 = "/mnt/devon-r2/shared_data/copernicus-dem-30m/DEM.vrt",
            nasadem_v1 = "/mnt/devon-r2/shared_data/NASADEM/mosaic/NASADEM_hgt.vrt",
            arcticdem_v3_10m = "/mnt/devon-r2/shared_data/arcticdem/mosaics/v3.0/10m/arcticdem_v3.0_10m.vrt",
            arcticdem_v4_10m = "/mnt/devon-r2/shared_data/arcticdem/mosaics/v4.1/10m/arcticdem_v4.1_10m.vrt",
            #arcticdem_v4_10m = "/mnt/baffin-r1/shared_data/arcticdem/mosaics/v4.1/10m/arcticdem_v4.1_10m.vrt",
            rema_v2_10m = "/mnt/devon-r2/shared_data/rema/mosaics/v2.0/10m/rema_v2.0_10m.vrt",
            #rema_v2_10m = "/mnt/baffin-r1/shared_data/rema/mosaics/v2.0/10m/rema_v2.0_10m.vrt",
            hugonnet_v1_stacks = "/mnt/bylot-r3/data/hugonnet/",
            canopyheight_10m_v1 = "/mnt/devon-r0/shared_data/canopy_height/ETH_GlobalCanopyHeight_10m_2020_version1/ETH_GlobalCanopyHeight_10m_2020_mosaic_Map.vrt",

            # location of goid files
            geoid_dir = "/mnt/devon-r2/shared_data/geoids",

            # locaiton of its_live parameter files
            itslive_parameters = "/mnt/devon-r2/data/its-live-data/autorift_parameters/v001/",   

            # masks
            icemask = "/mnt/devon2-r1/devon0/gardnera/Documents/GitHub/Altim/data/land_ice_0.0083_cog.tif",

            # vector files
            RGI_dissolved = "/mnt/devon-r2/data/GlacierOutlines/RGI2000-v7.0-C-global-fix",
            RGI_regions = "/mnt/devon-r2/data/GlacierOutlines/RGI2000-v7.0-regions",
            glacier_shp="/mnt/bylot-r3/data/vector_files/glacier_cl0cl1.shp",
            glacier_rgi7_shp="/mnt/bylot-r3/data/vector_files/RGI2000-v7.0-C-01_global.shp",
            glacier_b1km_shp="/mnt/bylot-r3/data/vector_files/glacier_b1km.shp",
            glacier_b10km_shp="/mnt/bylot-r3/data/vector_files/glacier_b10km.shp",
            glacier_individual_outlines="/mnt/bylot-r3/data/GlacierOutlines/rgi60/rgi60_Global.gpkg",
            landice_shp="/mnt/bylot-r3/data/vector_files/land_ice.shp",
            landice_b1km_shp="/mnt/bylot-r3/data/vector_files/land_ice_b1km.shp",
            landice_b10km_shp="/mnt/bylot-r3/data/vector_files/land_ice_b10km.shp",
            floating_shp="/mnt/bylot-r3/data/vector_files/floating_ice.shp",
            not_ocean_shp="/mnt/bylot-r3/data/vector_files/not_ocean.shp", # GGSHHS_f_L1 +  scripps_antarctica_grounded
            rgi6_regions_shp="/mnt/bylot-r3/data/vector_files/00_rgi60_O1Regions_fix.shp",
            water_shp="/mnt/bylot-r3/data/vector_files/GSHHS_f_L2.shp", #GSHHG lakes
            continents="/mnt/bylot-r3/data/vector_files/World_Continents.geojson",
            countries="/mnt/bylot-r3/data/vector_files/ne_10m_admin_0_countries.geojson",

            # grace timeseries data for validation [created in MATLAB using jplMasconTS.m]
            grace_rgi = "/mnt/bylot-r3/data/GRACE/dM_grace_rgi.mat",

            # glacier mass balance data for comparision
            zemp_2019 = "/mnt/bylot-r3/data/glacier_mb/Zemp2019/Zemp_etal_DataTables2a-t_results_regions_global/",
            marzeion_2012 = "/mnt/bylot-r3/data/glacier_mb/Marzeion2012/Marzeion2012.mat",
            marzeion_2020 = "/mnt/bylot-r3/data/glacier_mb/Marzeion2020/suppl_GlacierMIP_results.nc",
            hock_2019 = "/mnt/bylot-r3/data/glacier_mb/Hock2019/S0022143019000224sup001.nc",
            ipcc_ar6 = "/mnt/bylot-r3/data/glacier_mb/IPCCAR6/",
            wimberly_2024 = "/mnt/bylot-r3/data/glacier_mb/Wimberly2024/all_rf_aligned_data.csv",
            glambie_2024 = "/mnt/bylot-r3/data/glacier_mb/glambie/41586_2024_8545_MOESM4_ESM_Gt.csv",

            # path to river files
            river = "/mnt/bylot-r3/data/rivers/MERIT_Hydro_v07_Basins_v01",
            river_basins = "/mnt/bylot-r3/data/rivers/BasinATLAS_Data_v10.gdb/BasinATLAS_v10_lev02.geojson",
            river_major_basins = "/mnt/bylot-r3/data/rivers/GRDC_Major_River_Basins/mrb_basins.json",
        )

        return pathlocal
    else
        error("unrecognized machine, Atlim.jl requires local paths to be set in setpaths() in order to run")
    end
end