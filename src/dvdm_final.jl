
"""
Synthesize and process glacier mass change data from multiple sources.

This script combines and processes glacier mass change data from various sources 
(altimetry, GRACE, GEMB) to create a comprehensive synthesis dataset. It handles
data from different regions, missions and processing methods to generate uncertainty
estimates.

Key operations:
1. Synthesizes mass change data using multiple processing configurations:
   - Different surface masks (glacier, glacier_rgi7, glacier_b1km)
   - Different DEMs (best, COP30)
   - Various binning methods and parameters
   - Curvature and amplitude corrections
2. Combines regions into larger groupings (e.g. High Mountain Asia)
3. Adds GRACE satellite data and global estimates
4. Calculates trends and uncertainties
5. Saves final processed dataset

Parameters:
- geotile_width: Width of geotiles in degrees (default: 2)
- project_id: Project identifier (default: :v01) 
- dvdm_synthesis_id: Identifier for synthesis dataset
- Various processing parameters for uncertainty estimation

The script outputs:
- A JLD2 file containing the processed DataFrame with:
  - Mass changes for each region/mission
  - Uncertainty estimates
  - Trends and accelerations
  - Combined regional estimates

Dependencies:
Altim, FileIO
"""

begin
    using Altim
    using FileIO

    geotile_width = 2
    project_id = :v01

    dvdm_synthesis_id = "glacier_filtered"
    
    df = Altim.geotile_dvdm_synthesize(;
        # best estimate,
        project_id,
        surface_mask_best = "glacier",
        dem_best="best",
        curvature_correct_best=true,
        amplitude_correct_best=true,
        binning_method_best="meanmadnorm3",
        fill_param_best = 1,
        binned_folder_best="/mnt/bylot-r3/data/binned/2deg",

        # to include in uncertainty
        surface_masks=["glacier" "glacier_rgi7" "glacier_b1km"],
        dems=["best" "cop30_v2"],
        curvature_corrects=[false true],
        amplitude_corrects=[true],
        binning_methods=["median" "meanmadnorm5" "meanmadnorm3" "meanmadnorm10"],# ["median" "meanmadnorm10" "meanmadnorm5" "meanmadnorm3"]
        fill_params=[1, 2, 3, 4],
        binned_folders=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"),

        # manual adjustments
        regions_to_overwrite_hugonnet_data_with_model_fit=["rgi19"],

        ## combine regions at the per-sensor level
        # this gets a bit complicated a only regions with consistnet mission inclusion can be grouped for proper 
        region_col="rgi",

        region_combines_before_fac=(
            ("rgi13", "rgi14", "rgi15") => "hma",
        ),

        region_combines_after_fac=(
            ("rgi1", "rgi3", "rgi4", "rgi5", "rgi6", "rgi7", "rgi8", "rgi9", "rgi19") => "hll",
            ("rgi2", "rgi10", "rgi11", "rgi12", "hma", "rgi17", "rgi18") => "ghll",
            ("rgi1", "rgi3", "rgi4", "rgi6", "rgi7", "rgi8", "rgi9") => "hll_ep",
        ),

        combine_vars=["area_km2", "val", "nobs"],

        missions_to_filter_on_nobs_km2=["icesat", "gedi"],

        path2gemb = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0_geotile_filled_d_reg.jld2",
        )

    df = Altim.geotile_dvdm_addgrace!(df)
    df = Altim.geotile_dvdm_global!(df)

    # remove GRACE rgi12 wich is a syntheis transplant for generating global extimate
    index = (df.rgi .== "rgi12") .& (df.mission .== "grace")
    df = df[.!index, :]
    
    df = Altim.geotile_dvdm_add_trend!(df; iterations=1000)

    final_data_dir = joinpath(Altim.pathlocal.data_dir, "altim_final", "$(geotile_width)_$(project_id)")

    final_filename = "dvdm_$(dvdm_synthesis_id).jld2"

    save(joinpath(final_data_dir, final_filename), Dict("df" => df))
end