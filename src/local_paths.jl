"""
    setpaths()

Local file paths on the JPL servers (bylot, devon, baffin): the paths of the volume chain from
`GlacierVolumeChange.setpaths()`, plus the validation, river, discharge and output files this package
uses. Throws an error on any other host.
"""
function setpaths()
    gvc_paths = GVC.setpaths()
    project_dir = gvc_paths.project_dir

    pathlocal = (;
        gvc_paths...,

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
        # HydroBASINS_v1c major river basins:https://www.hydrosheds.org/products/hydrobasins
        river_level03_as_basins = "/mnt/bylot-r3/data/rivers/HydroBASINS_v1c/hybas_as_lev01-12_v1c/hybas_as_lev03_v1c.shp",
        river_level03_si_basins="/mnt/bylot-r3/data/rivers/HydroBASINS_v1c/hybas_si_lev01-12_v1c/hybas_si_lev03_v1c.shp",

        # Kochtitzky NH discharge and terminus retreate 
        discharge_nh="/mnt/bylot-r3/data/GlacierOutlines/GlacierDischarge/Kochtitzky2022/41467_2022_33231_MOESM4_ESM.csv",
        discharge_npi="/mnt/bylot-r3/data/GlacierOutlines/GlacierDischarge/Fuerst2023/fuerst_2023_npi_comparison_ice_discharge_v1.0.0.txt",
        discharge_spi="/mnt/bylot-r3/data/GlacierOutlines/GlacierDischarge/Fuerst2023/fuerst_2023_spi_comparison_ice_discharge_v1.0.0.txt",
        rgi6_southern_andes = "/mnt/bylot-r3/data/GlacierOutlines/rgi60/17_rgi60_SouthernAndes/17_rgi60_SouthernAndes.shp",

        # discharge data for global discharge analysis [created as part of processing pipeline]
        discharge_global = "/mnt/bylot-r3/data/GlacierOutlines/GlacierDischarge/global_glacier_discharge.jld2",

        glacier_summary = joinpath(project_dir, "Gardner2025_glacier_summary.nc"),

        era5_land_sea_mask = "/mnt/bylot-r3/data/era5/LandSeaMask.nc",
    )

    return pathlocal
end

pathlocal = setpaths()
