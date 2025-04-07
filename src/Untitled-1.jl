reference_run = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_dh_best_cc_meanmadnorm3_v01_filled_ac_p1.jld2"

gtid = "lat[+58+60]lon[-132-130]"

f = load(reference_run)

p = CairoMakie.lines(f["dh_hyps"]["icesat2"][At(gtid),:,At(2050)], label="ICESat-2");
CairoMakie.lines!(f["dh_hyps"]["icesat"][At(gtid), :, At(2050)], label="ICESat");
CairoMakie.lines!(f["dh_hyps"]["hugonnet"][At(gtid), :, At(2050)], label="Hugonnet");
CairoMakie.lines!(f["dh_hyps"]["gedi"][At(gtid), :, At(2050)], label="GEDI"); p


reference_run = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_dh_best_cc_meanmadnorm3_v01_filled_p1.jld2"
f = load(reference_run)

p = CairoMakie.lines(f["dh_hyps"]["icesat2"][At(gtid), :, At(2150)], label="ICESat-2");
CairoMakie.lines!(f["dh_hyps"]["icesat"][At(gtid), :, At(2150)], label="ICESat");
CairoMakie.lines!(f["dh_hyps"]["hugonnet"][At(gtid), :, At(2150)], label="Hugonnet");
CairoMakie.lines!(f["dh_hyps"]["gedi"][At(gtid), :, At(2150)], label="GEDI");
p






    project_id=:v01
    geotile_width=2
    force_remake_before=DateTime(2029, 1, 1)
    update_geotile=false # this will load in prevous results to update select geotiles or missions
    update_geotile_missions=["icesat2"]
    plot_dh_as_function_of_time_and_elevation=true
    mission_reference_for_amplitude_normalization="icesat2"
    all_permutations_for_glacier_only=true
    surface_masks=[:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km]
    binned_folders=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg")
    dem_ids=[:best, :cop30_v2]
    binning_methods=["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"]
    curvature_corrects=[true, false]
    paramater_sets=[1, 2]
    amplitude_corrects=[true, false]
    showplots = true


    params = NamedTuple[]
    for binned_folder in binned_folders
        for paramater_set in paramater_sets
            for surface_mask in surface_masks
                for dem_id in dem_ids
                    for binning_method in binning_methods
                        for curvature_correct in curvature_corrects
                            for amplitude_correct = amplitude_corrects
                                push!(params, (; project_id, binned_folder, paramater_set, surface_mask, dem_id, binning_method, curvature_correct, amplitude_correct))
                            end
                        end
                    end
                end
            end
        end
    end

    # Threads is throwing errors due to reading of JLD2 files, Threads is implimented at
    # lower level with reasonable performance
    #@showprogress desc = "Filling hypsometric elevation change data ..." for param in params
        param = first(params)

        fig_folder = joinpath(param.binned_folder, "figures")
        !isdir(fig_folder) && mkdir(fig_folder)

        param_filling = Altim.binned_filling_parameters[param.paramater_set]

        # load geotiles
        geotiles = Altim.geotiles_mask_hyps(param.surface_mask, geotile_width)

        # make geotile rgi regions mutually exexclusive 
        geotiles, reg = Altim.geotiles_mutually_exclusive_rgi!(geotiles)

        # skip permutations if all_permutations_for_glacier_only = true
        if all_permutations_for_glacier_only
            if ((!(param.surface_mask == :glacier) && Base.contains(param.binned_folder, "unfiltered")) &&
                ((param.dem_id != :best) && (param.binning_method != "meanmadnorm3") && param.curvature_correct))

                continue
            end
        end

        # paths to files
        binned_file = Altim.binned_filepath(; param.binned_folder, param.surface_mask, param.dem_id, param.binning_method, project_id, param.curvature_correct)
        binned_file_land = Altim.binned_filepath(; param.binned_folder, surface_mask=:land, param.dem_id, param.binning_method, project_id, param.curvature_correct)
        binned_filled_file, figure_suffix = Altim.binned_filled_filepath(; param.binned_folder, param.surface_mask, param.dem_id, param.binning_method, project_id, param.curvature_correct, param.amplitude_correct, param.paramater_set)

        #if (!isnothing(force_remake_before) || !isfile(binned_filled_file)) && isfile(binned_file)

            ##################################### HACK FOR UPDATE #####################################
            if isfile(binned_filled_file) && !isnothing(force_remake_before)
                if Dates.unix2datetime(mtime(binned_file)) > force_remake_before
                    @warn "Skipping $(binned_filled_file) because it was created after $force_remake_before"
                    continue
                end
            end

            println("\nfilling binned data: surface_mask = $(param.surface_mask); binning_method = $(param.binning_method); dem_id = $(param.dem_id); curvature_correct = $(param.curvature_correct); amplitude_correct = $(param.amplitude_correct)")

            dh1 = FileIO.load(binned_file, "dh_hyps")

   
            if .!any(.!isnan.(dh1["hugonnet"]))
                println("NOT DATA: $binned_file")
                continue
            end


            if update_geotile && isfile(binned_filled_file)
                old_keys = setdiff(keys(dh1), update_geotile_missions)
                for k in old_keys
                    delete!(dh1, k)
                end
            end

            nobs1 = FileIO.load(binned_file, "nobs_hyps")

            # align geotile dataframe with DimArrays
            gt = collect(dims(dh1[first(keys(dh1))], :geotile))
            gt_ind = [findfirst(geotiles.id .== g0) for g0 in gt]
            geotiles = geotiles[gt_ind, :]

            # create a data frame to store model parameters
            # initialize dimensional arrays
            params_fill = Dict()
            for mission in keys(dh1)
                n = length(dims(dh1[mission], :geotile))
                push!(params_fill, mission => DataFrame(geotile=val(dims(dh1[mission], :geotile)), nobs_raw=zeros(n), nbins_raw=zeros(n), nobs_final=zeros(n), nbins_filt1=zeros(n), param_m1=[fill(NaN, size(Altim.p1)) for i in 1:n], h0=fill(NaN, n), t0=fill(NaN, n), dh0=fill(NaN, n), bin_std=fill(NaN, n), bin_anom_std=fill(NaN, n)))
            end

            if plot_dh_as_function_of_time_and_elevation
                dh_time_elevation_idx = findfirst(geotiles.id .== "lat[+28+30]lon[+082+084]")
                Altim.plot_height_time(dh1; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="raw", fig_folder, figure_suffix, mask=param.surface_mask, showplots)
            end

            # I think that the model fitting can have a large influence on how off ice dh is treated. When lots of off ice area is included (e.g. param.surface_mask == glacier_b1km of glacier_b10km) there can be a very sharp transition at lower elevations from high rates to low rates of elevation change. 
            if param_filling.adjust2land && any(keys(dh1) .== "huggonet")
                dh_land = FileIO.load(binned_file_land, "dh_hyps")
                nobs_land = FileIO.load(binned_file_land, "nobs_hyps")

                dh1["hugonnet"] = Altim.geotile_adjust!(
                    dh1["hugonnet"],
                    dh_land["hugonnet"],
                    nobs_land["hugonnet"],
                    dh_land["icesat2"];
                    ref_madnorm_max=10,
                    minnobs=45,
                    ref_minvalid=12,
                )
            end

            Altim.hyps_model_fill!(dh1, nobs1, params_fill; bincount_min=param_filling.bincount_min,
                model1_madnorm_max=param_filling.model1_madnorm_max, smooth_n=param_filling.smooth_n,
                smooth_h2t_length_scale=param_filling.smooth_h2t_length_scale)

            variogram_range_ratio = false

            # make plot of the height to time variogram range ratio
            if variogram_range_ratio
                range_ratio = Altim.hyps_model_fill!(dh1, nobs1, params_fill; bincount_min=param_filling.bincount_min,
                    model1_madnorm_max=param_filling.model1_madnorm_max, smooth_n=param_filling.smooth_n,
                    smooth_h2t_length_scale=param_filling.smooth_h2t_length_scale, variogram_range_ratio)

                fontsize = 18
                f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98), size=(700, 700), fontsize=fontsize)

                for (i, mission) in enumerate(keys(dh1))
                    valid = .!isnan.(range_ratio[mission])
                    title = replace.(mission, "hugonnet" => "aster")

                    title = "$title  mean = $(round(Int, mean(range_ratio[mission][valid])))"
                    Axis(f[i, 1], title=title)
                    CairoMakie.hist!(collect(range_ratio[mission][valid]); title=mission, bins=0:250:3000)
                end

                fname = joinpath(fig_folder, "$(figure_suffix)_variogram_range_ratio.png")
                save(fname, f)
                display(f)
            end

            # filter and fit model to geotiles [4 min for for all glacierized geotiles and 4 missions]
            if plot_dh_as_function_of_time_and_elevation
                Altim.plot_height_time(dh1; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="filled", fig_folder, figure_suffix, mask=param.surface_mask)
            end

            # apply seasonal amplitude normalization
            if param.amplitude_correct

                # check if reference mission is present, if not, add
                ref_added = false
                if update_geotile && (!any(keys(dh1) .== mission_reference_for_amplitude_normalization)) && isfile(binned_filled_file)
                    ref_added = true
                    (dh_hyps, model_param) = FileIO.load(binned_filled_file, ("dh_hyps", "model_param"))

                    for k in setdiff(keys(dh_hyps), [mission_reference_for_amplitude_normalization])
                        delete!(dh_hyps, k)
                        delete!(model_param, k)
                    end

                    dh1 = merge(dh1, dh_hyps)
                    params_fill = merge(params_fill, model_param)
                end

                Altim.hyps_amplitude_normalize!(dh1, params_fill; mission_reference=mission_reference_for_amplitude_normalization)

                # remove mission_reference_for_amplitude_normalization if it was added.
                if ref_added
                    delete!(dh1, mission_reference_for_amplitude_normalization)
                    delete!(params_fill, mission_reference_for_amplitude_normalization)
                end
            end

            if plot_dh_as_function_of_time_and_elevation
                Altim.plot_height_time(dh1; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="filled_ampnorm", fig_folder, figure_suffix, mask=param.surface_mask, showplots)
            end

            dh1 = Altim.hyps_fill_empty!(dh1, params_fill, geotiles; mask=param.surface_mask)

            dh1 = Altim.hyps_fill_updown!(dh1, geotiles; mask=param.surface_mask)

            if plot_dh_as_function_of_time_and_elevation
                Altim.plot_height_time(dh1; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="filled_updown", fig_folder, figure_suffix, mask=param.surface_mask, showplots)
            end

            if update_geotile && isfile(binned_filled_file)
                (dh_hyps, nobs_hyps, model_param) = FileIO.load(binned_filled_file, ("dh_hyps", "nobs_hyps", "model_param"))


                for k in update_geotile_missions
                    delete!(dh_hyps, k)
                    delete!(nobs_hyps, k)
                    delete!(model_param, k)
                end

                dh1 = merge(dh1, dh_hyps)
                nobs1 = merge(nobs1, nobs_hyps)
                params_fill = merge(params_fill, model_param)
            end

            # save filled geotiles
            save(binned_filled_file, Dict("dh_hyps" => dh1, "nobs_hyps" => nobs1, "model_param" => params_fill))
        end
    end
end
