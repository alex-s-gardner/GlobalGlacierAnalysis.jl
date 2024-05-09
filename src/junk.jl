using Altim

project_id = :v01;
geotile_width = 2;
binned_folder = analysis_paths(; geotile_width).binned
binning_method_gemb = "mean"; 
runid_gemb = "glacier_gemb_$(binning_method_gemb)_$(project_id)"
gemb_file = joinpath(binned_folder, "$(runid_gemb).jld2");


 fac = load(gemb_file, "fac_hyps");
    nobs_gemb = load(gemb_file, "nobs_hyps");
    smb = load(gemb_file, "smb_hyps");
    
