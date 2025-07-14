import GlobalGlacierAnalysis as GGA

fn = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_rgi7_dh_best_nmad3_v01_filled_ac_p3_aligned.jld2"

f, dh, param = GGA.load_and_plot_elevation_time_multimission_geotiles(;
    path2file=fn,
    geotiles2plot=GGA.geotiles_golden_test,
    colorrange=(-20, 20),
);

f[GGA.geotiles_golden_test[1]]
f[GGA.geotiles_golden_test[2]]
f[GGA.geotiles_golden_test[3]]
f[GGA.geotiles_golden_test[4]]