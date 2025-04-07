 binned_synthesized_file = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_rgi7_dh_best_median_v01_filled_ac_p1.jld2"

            dh = FileIO.load(binned_synthesized_file, "dh_hyps")

            geotile_num = 802;
            p = CairoMakie.plot(dh["icesat2"][geotile_num,:,At(2550)])
            CairoMakie.plot!(dh["icesat"][geotile_num,:,At(2550)]); p
            CairoMakie.plot!(dh["gedi"][geotile_num,:,At(2550)]); p
            CairoMakie.plot!(dh["hugonnet"][geotile_num,:,At(2550)]); p

            geotile_num = 1103;
            p = CairoMakie.plot(dh["icesat2"][geotile_num,:,At(2550)])
            CairoMakie.plot!(dh["icesat"][geotile_num,:,At(2550)]); p
            CairoMakie.plot!(dh["gedi"][geotile_num,:,At(2550)]); p
            CairoMakie.plot!(dh["hugonnet"][geotile_num,:,At(2550)]); p