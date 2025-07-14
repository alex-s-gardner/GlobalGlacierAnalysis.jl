import GlobalGlacierAnalysis as GGA

fn = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_rgi7_dh_cop30_v2_median_v01_filled_ac_p4_aligned.jld2"
fn = "/mnt/bylot-r3/data/binned/2deg/glacier_dh_best_cc_nmad5_v01_filled_ac_p4_aligned.jld2"
geotiles2plot = ["lat[+54+56]lon[-126-124]"] #GGA.geotiles_golden_test

#TODO: add plot of location with glacier coverage
#TODO: merge icesat and icesat2 if Synthesis exists
f, dh, param = GGA.load_and_plot_elevation_time_multimission_geotiles(;
    path2file=fn,
    geotiles2plot
);

for geotile in geotiles2plot
    display(f[geotile])
end






