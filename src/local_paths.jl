## MODIFY LOCAL PATHS FOR EACH MACHINE
# TODO:  maybe use "Preferences.jl" @set_preferences!(pref_name => key) and @load_preference(pref_name) macros in future

"""
    setpaths()

Returns local paths
"""

function setpaths() 
    hostname = gethostname()
    if hostname == "MT-502549"

        data_dir = "$(userdir)/data/"
        pathlocal = (

            # all altimetry and geotile data will be stored in a stucture within data_dir
            data_dir=data_dir,

            # location of global DEM vrts
            cop30_v1="$(data_dir)/COP-DEM_GLO-30-DGED_PUBLIC/GLO30_DGED_hgt.vrt",
            cop30_v2="$(data_dir)/copernicus-dem-30m/DEM.vrt",
            nasadem_v1="$(data_dir)/NASADEM/mosaic/NASADEM_hgt.vrt",
            arcticdem_v3_10m="$(data_dir)/arcticdem/mosaics/v3.0/10m/arcticdem_v3.0_10m.vrt",
            rema_v2_10m="$(data_dir)/rema/mosaics/v2.0/10m/rema_v2.0_10m.vrt",
            hugonnet_v1_stacks="$(data_dir)/hugonnet/",
            canopyheight_10m_v1="$(data_dir)/canopy_height/ETH_GlobalCanopyHeight_10m_2020_version1/ETH_GlobalCanopyHeight_10m_2020_mosaic_Map.vrt",
            
            # location of goid files
            geoid_dir="$(data_dir)/geoids",

            # locaiton of its_live parameter files
            itslive_parameters="$(data_dir)/its-live-data/autorift_parameters/v001/",

            # masks
            icemask = "/Users/gardnera/data/GlacierOutlines/WorldMask_20190513/land_ice_0.0083_cog.tif"
        )

        return pathlocal

    elseif (hostname == "bylot.jpl.nasa.gov") || (hostname == "devon.jpl.nasa.gov") || (hostname == "baffin.jpl.nasa.gov")

        pathlocal = (
            # all altimetry and geotile data will be stored in a stucture within data_dir
            data_dir = "/mnt/bylot-r3/data/",

            # location of global DEM vrts
            cop30_v1 = "/mnt/devon-r2/shared_data/COP-DEM_GLO-30-DGED_PUBLIC/GLO30_DGED_hgt.vrt",
            cop30_v2 = "/mnt/devon-r2/shared_data/copernicus-dem-30m/DEM.vrt",
            nasadem_v1 = "/mnt/devon-r2/shared_data/NASADEM/mosaic/NASADEM_hgt.vrt",
            arcticdem_v3_10m = "/mnt/devon-r2/shared_data/arcticdem/mosaics/v3.0/10m/arcticdem_v3.0_10m.vrt",
            #arcticdem_v4_10m = "/mnt/devon-r2/shared_data/arcticdem/mosaics/v4.1/10m/arcticdem_v4.1_10m.vrt",
            arcticdem_v4_10m = "/mnt/baffin-r1/shared_data/arcticdem/mosaics/v4.1/10m/arcticdem_v4.1_10m.vrt",
            #rema_v2_10m = "/mnt/devon-r2/shared_data/rema/mosaics/v2.0/10m/rema_v2.0_10m.vrt",
            rema_v2_10m = "/mnt/baffin-r1/shared_data/rema/mosaics/v2.0/10m/rema_v2.0_10m.vrt",
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
            glacier_shp="/mnt/bylot-r3/data/vector_files/glacier.shp",
            glacier_b1km_shp="/mnt/bylot-r3/data/vector_files/glacier_b1km.shp",
            glacier_b10km_shp="/mnt/bylot-r3/data/vector_files/glacier_b10km.shp",
            landice_shp="/mnt/bylot-r3/data/vector_files/land_ice.shp",
            landice_b1km_shp="/mnt/bylot-r3/data/vector_files/land_ice_b1km.shp",
            landice_b10km_shp="/mnt/bylot-r3/data/vector_files/land_ice_b10km.shp",
            floating_shp="/mnt/bylot-r3/data/vector_files/floating_ice.shp",
            rgi6_regions_shp="/mnt/bylot-r3/data/vector_files/00_rgi60_O1Regions_fix.shp",
        )

        return pathlocal
    else
        error("unrecognized machine, Atlim.jl requires local paths to be set in setpaths() in order to run")
    end
end