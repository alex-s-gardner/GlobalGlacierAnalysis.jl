begin
    import GlobalGlacierAnalysis as GGA
    using FileIO
    using DimensionalData

    # make figures for the paper
    geotile2plot = "lat[+28+30]lon[+082+084]"

    modelrun2plot_notfilled = "/mnt/bylot-r3/data/binned/2deg/glacier_dh_best_cc_nmad5_v01.jld2"
    fill_suffix = "ac_p2"
    #fill_suffix = "p2"
end

begin
    dh = load(modelrun2plot_notfilled, "dh_hyps")
    f = GGA.plot_elevation_time_multimission(dh, geotile2plot)
    display(f)


    path2run = replace(modelrun2plot_notfilled, ".jld2" => "_filled_$(fill_suffix).jld2")
    dh = load(path2run, "dh_hyps")
    f = GGA.plot_elevation_time_multimission(dh, geotile2plot)
    display(f)

    path2run = replace(modelrun2plot_notfilled, ".jld2" => "_filled_$(fill_suffix)_aligned.jld2")
    dh = load(path2run, "dh_hyps")
    f = GGA.plot_elevation_time_multimission(dh, geotile2plot)
    display(f)


    path2run = replace(modelrun2plot_notfilled, ".jld2" => "_filled_$(fill_suffix)_synthesized.jld2")
    dh = load(path2run, "dh_hyps")
    f = GGA.plot_elevation_time(dh[geotile=At(geotile2plot)])
    display(f)
end