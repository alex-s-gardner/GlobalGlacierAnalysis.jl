# add amplitude scaling and aplitude magnitude to image 
begin
    using Arrow
    using Altim
    using DataFrames
    using Extents
    using DimensionalData
    using DimensionalData.LookupArrays
    using Statistics
    using JLD2
    using FileIO
    using RollingFunctions
    using NearestNeighbors
    using Interpolations
    using ScatteredInterpolation
    using DataInterpolations
    using GeoTiles
    using Distances
    using LsqFit
    using Dates
    using DataInterpolations
    using Plots

    include("utilities_hyps.jl")

    project_id = :v01;
    geotile_width = 2;
    fill_dh = false;
    showplots = false;

    regional_offset_apply = true;
    fac_scale_apply = true

    mission_reference = "icesat2"
    
    binned_folder = analysis_paths(; geotile_width).binned
    fig_folder = joinpath(binned_folder, "figures")
    !isdir(fig_folder) && mkdir( fig_folder)

    showplots && (using Plots)

    paramater_sets = [1,2,3,4]
    dem_ids = [:best, :cop30_v2]
    binning_methods = ["median", "meanmadnorm3", "meanmadnorm5", "meanmadnorm10"];
    curvature_corrects = [true, false]
    surface_masks = [:glacier, :glacier_b1km, :land, :glacier_b10km]
    amplitude_corrects = [true, false]

    if false
        surface_masks = [:glacier, :glacier_b1km, :land, :glacier_b10km]
        dem_ids = [:best]
        binning_methods = ["meanmadnorm3"]
        curvature_corrects = [true]
        amplitude_corrects = [true]
    end

    model2::Function = model2(h, p) = p[1] .+ p[2] .* h .+ p[3] .* h .^ 2;
    p2 = zeros(3);
    lb2 = [-30.0, -0.1, -0.01];
    ub2 = [+30.0, +0.1, 0.01];
    #plotly()
end


for paramater_set in paramater_sets
    if paramater_set == 1
        param = (
            bincount_min=Dict("icesat" => 9,
                "icesat2" => 9,
                "gedi" => 9,
                "hugonnet" => 51,
            ), smooth_n=Dict("icesat" => 5,
                "icesat2" => 5,
                "gedi" => 5,
                "hugonnet" => 21,
            ), smooth_h2t_length_scale=800, # 800 m = 1 year in distance for anomaly from variogram analysis =
            model1_madnorm_max=5, # this is a sigma-equivelent threshold
        )

    elseif paramater_set == 2
        param = (
            bincount_min=Dict("icesat" => 9,
                "icesat2" => 9,
                "gedi" => 9,
                "hugonnet" => 9,
            ), smooth_n=Dict("icesat" => 5,
                "icesat2" => 5,
                "gedi" => 5,
                "hugonnet" => 5,
            ), smooth_h2t_length_scale=800, # 800 m = 1 year in distance for anomaly from variogram analysis =
            model1_madnorm_max=10, # this is a sigma-equivelent threshold
        )
    elseif paramater_set == 3
        param = (
            bincount_min=Dict("icesat" => 11,
                "icesat2" => 11,
                "gedi" => 11,
                "hugonnet" => 11,
            ), smooth_n=Dict("icesat" => 5,
                "icesat2" => 5,
                "gedi" => 5,
                "hugonnet" => 11,
            ), smooth_h2t_length_scale=800, # 800 m = 1 year in distance for anomaly from variogram analysis =
            model1_madnorm_max=15, # this is a sigma-equivelent threshold
        )

    elseif paramater_set == 4
        param = (
            bincount_min=Dict("icesat" => 5,
                "icesat2" => 5,
                "gedi" => 5,
                "hugonnet" => 5,
            ), smooth_n=Dict("icesat" => 7,
                "icesat2" => 7,
                "gedi" => 7,
                "hugonnet" => 7,
            ), smooth_h2t_length_scale=800, # 800 m = 1 year in distance for anomaly from variogram analysis =
            model1_madnorm_max=5, # this is a sigma-equivelent threshold
        )
    end

    for surface_mask in surface_masks
        # ~5 min for filling and plotting per iteration, 1.5 hours for all permutations and combinations
        geotiles = geotiles_mask_hyps(surface_mask, geotile_width)
        # make geotile rgi regions mutually exexclusive 
        geotiles, reg = geotiles_mutually_exclusive_rgi!(geotiles)    

        for binning_method = binning_methods
            for dem_id in dem_ids 
                for curvature_correct in curvature_corrects
                    for amplitude_correct = amplitude_corrects

                        if false
                            binning_method = "meanmadnorm10"
                            dem_id = :best 
                            curvature_correct = true
                            amplitude_correct = true
                        end

                        # paths to files
                        binned_file = binned_filepath(binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct)
                        binned_filled_file, figure_suffix = binned_filled_filepath(binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)

                        #fill_dh = true; showplots = true; (using Plots)
                        if (fill_dh || !isfile(binned_filled_file)) && isfile(binned_file)
                            dh1 = load(binned_file, "dh_hyps");

                            if .!any(.!isnan.(dh1["hugonnet"]))
                                println("NOT DATA: $binned_file")
                                continue
                            end

                            nobs1 = load(binned_file, "nobs_hyps");

                            # align geotile dataframe with DimArrays
                            gt = collect(dims(dh1[first(keys(dh1))], :geotile))
                            gt_ind = [findfirst(geotiles.id .== g0) for g0 in gt]
                            geotiles = geotiles[gt_ind,:]

                            # create a data frame to store model parameters
                            # initialize dimensional arrays
                            params = Dict();
                            for  mission in keys(dh1)
                                n = length(dims(dh1[mission], :geotile))
                                push!(params, mission => DataFrame(geotile = val(dims(dh1[mission], :geotile)), nobs_raw = zeros(n), nbins_raw = zeros(n), nobs_final = zeros(n), nbins_filt1 = zeros(n), param_m1 = [fill(NaN, size(p1)) for i in 1:n], h0 = fill(NaN,n), t0 = fill(NaN,n), dh0 = fill(NaN,n), bin_std = fill(NaN,n), bin_anom_std =  fill(NaN,n)));
                            end

                            showplots && (i = 652);
                            showplots && (i = findlast(Extents.intersects.(Ref(Extent(X=(23.6, 25.8), Y=(79.2, 79.7))), geotiles.extent)))

                            showplots && plot_height_time(dh1; geotile = geotiles[i,:], fig_suffix = "raw", fig_folder, figure_suffix, mask = surface_mask)

                            println("surface_mask = $surface_mask; binning_method = $binning_method; dem_id = $dem_id; curvature_correct = $curvature_correct; amplitude_correct = $amplitude_correct")

                            # I think that the model fitting can have a large influence on how off ice dh is treated. When lots of off ice area is included (e.g. surface_mask == glacier_b1km of glacier_b10km) there can be a very sharp transition at lower elevations from high rates to low rates of elevation change. 
                            hyps_model_fill!(dh1, nobs1, params; bincount_min = param.bincount_min, 
                                model1_madnorm_max = param.model1_madnorm_max, smooth_n = param.smooth_n,  
                                smooth_h2t_length_scale = param.smooth_h2t_length_scale)

                            variogram_range_ratio = false
                            # make plot of the height to time variogram range ratio
                            if variogram_range_ratio
                                range_ratio = hyps_model_fill!(dh1, nobs1, params; bincount_min = param.bincount_min, 
                                model1_madnorm_max = param.model1_madnorm_max, smooth_n = param.smooth_n,  
                                smooth_h2t_length_scale = param.smooth_h2t_length_scale, variogram_range_ratio)

                                using GeoStats
                                fontsize = 18
                                f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98), size=(700, 700), fontsize=fontsize);

                                for (i, mission) in enumerate(keys(dh1))
                                        valid = .!isnan.(range_ratio[mission])
                                        title = replace.(mission, "hugonnet" => "aster")

                                        title = "$title  mean = $(round(Int, mean(range_ratio[mission][valid])))"
                                        Axis(f[i, 1], title = title)
                                        CairoMakie.hist!(collect(range_ratio[mission][valid]); title = mission,  bins=0:250:3000)
                                end

                                fname = joinpath(fig_folder, "variogram_range_ratio_$(figure_suffix)_p$(paramater_set).png")
                                save(fname, f)
                                display(f)
                            end

                            # filter and fit model to geotiles [4 min for for all glacierized geotiles and 4 missions]        
                            showplots && plot_height_time(dh1; geotile = geotiles[i,:], fig_suffix = "filled", fig_folder, figure_suffix, mask = surface_mask)
                            
                            # apply seasonal amplitude normalization
                            if amplitude_correct
                                hyps_amplitude_normalize!(dh1, params; mission_reference)
                            end
                            
                            showplots && plot_height_time(dh1; geotile = geotiles[i,:], fig_suffix = "filled_ampnorm", fig_folder, figure_suffix, mask = surface_mask)

                            hyps_fill_empty!(dh1, params, geotiles; mask = surface_mask)

                            hyps_fill_updown!(dh1, geotiles; mask = surface_mask)

                            showplots && plot_height_time(dh1; geotile = geotiles[i,:], fig_suffix = "filled_updown", fig_folder, figure_suffix, mask = surface_mask)

                            # save filled geotiles
                            save(binned_filled_file, Dict("dh_hyps" => dh1, "nobs_hyps" => nobs1, "model_param" =>   params ));
                        end
                    end
                end
            end
        end
    end


    # find offsets to low vegitated land surfaces
    begin
        surface_mask = :land
        binning_method = "meanmadnorm3"
        dem_id = :best 
        curvature_correct = true
        amplitude_correct = false

        # ~5 min for filling and plotting per iteration, 1.5 hours for all permutations and combinations
        geotiles = geotiles_mask_hyps(surface_mask, geotile_width)

        # make geotile rgi regions mutually exexclusive 
        geotiles, reg = geotiles_mutually_exclusive_rgi!(geotiles)    

        # bin method
        if curvature_correct
            runid = "$(surface_mask)_dh_$(dem_id)_cc_$(binning_method)_$(project_id)"
        else
            runid = "$(surface_mask)_dh_$(dem_id)_$(binning_method)_$(project_id)"
        end

        landfit_param_file = joinpath(binned_folder, "$(runid)_fit_param_p$(paramater_set).jld2")
        if !isfile(landfit_param_file)

            # paths to files
            binned_file = binned_filepath(binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct)
            binned_filled_file, figure_suffix = binned_filled_filepath(binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)

            model_param = load(binned_filled_file, "model_param")
            dh1 = load(binned_filled_file, "dh_hyps");
            nobs1 = load(binned_filled_file, "nobs_hyps");

            # align geotiles
            gt = collect(dims(dh1[first(keys(dh1))], :geotile))
            gt_ind = [findfirst(geotiles.id .== g0) for g0 in gt]
            geotiles = geotiles[gt_ind,:]

            drgi = Dim{:rgi}(reg)
            ddate = dims(dh1[first(keys(dh1))], :date)
            dmission = Dim{:mission}(collect(keys(dh1)))
            dh_reg = fill(NaN, (dmission, drgi, ddate))
            nobs_reg = fill(NaN, (dmission, drgi, ddate))

            for mission in dmission
                nopbs_sum = sum(nobs1[mission], dims = :height);
                nopbs_sum = dropdims(nopbs_sum, dims = :height)
                foo = dh1[mission] .* nobs1[mission]
                foo[isnan.(foo)] .= 0
                dh_mean = sum(foo, dims = :height) ./ nopbs_sum
                dh_mean = dropdims(dh_mean, dims = :height)

                for rgi in drgi
                    rgi_ind = geotiles[:, rgi] .> 0
                    geotile_id = geotiles[rgi_ind, :id]

                    foo = dh_mean[At(geotile_id), :];
                    foo[isnan.(foo)] .= 0
                    dh_reg[At(mission), At(rgi), :] = sum(foo .* nopbs_sum[At(geotile_id), :], dims= :geotile) ./ sum(nopbs_sum[At(geotile_id), :], dims= :geotile)

                    nobs_reg[At(mission), At(rgi), :] = sum(nopbs_sum[At(geotile_id), :], dims= :geotile)
                end
            end

            decdate = Altim.decimalyear.(lookup(dh_reg, :date))
            dmission = dims(dh_reg, :mission)
            dmetric = Dim{:metric}(["mean", "trend", "acceleration", "amplitude", "date_intercept"])
            fit_param = fill(NaN, (dmission, drgi, dmetric))

            for rgi in drgi
                title ="Randolph Glacier Inventory: Region $(rgi[4:end])"
                f, fit_param[:,At(rgi),:] = plot_dh(dh_reg[:, At(rgi), :], (nobs_reg[:, At(rgi), :]); title, xlims=(DateTime(2000), DateTime(2024)))

                fname = joinpath(fig_folder, "$(figure_suffix)_$(rgi)_p$(paramater_set).png")
                save(fname, f)
                showplots && display(f)
            end

            save(landfit_param_file, Dict("landfit" => fit_param));
        end
        
        landfit = load(landfit_param_file, "landfit");
        
    end

    # firn correction [same for all iterations so simply calculate once]
    begin
        surface_mask = :glacier
        binning_method = "meanmadnorm3"
        dem_id = :best
        curvature_correct = true
        amplitude_correct = true
        mission_ref_fac = "hugonnet"

        # paths to files
        binned_file = binned_filepath(binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct)
        binned_filled_file, figure_suffix = binned_filled_filepath(binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)

        idx = findlast("_", binned_filled_file)
        foo = binned_filled_file[1:(idx[1]-1)]

        file_parts = splitpath(foo)
        file_parts[end] = "gemb_"* file_parts[end] * ".jld2"
        gemb_fit_file = joinpath(file_parts);

        if !isfile(gemb_fit_file)
            dh1 = load(binned_filled_file, "dh_hyps")

            fac = load(gemb_fit_file, "fac")

            # make FAC cube
            fac = copy(dh1[mission_ref_fac])
            fac[:] .= NaN
            smb = copy(fac)

            geotiles = geotiles_mask_hyps(surface_mask, geotile_width)

            # make geotile rgi regions mutually exexclusive 
            geotiles, reg = geotiles_mutually_exclusive_rgi!(geotiles)    

            geotiles = geotiles[(geotiles.glacier_frac .> 0), :]

            # find the best matched GEMB data within geotile
            geotile_buffer = 1;
            # this takes 1hr
            Threads.@threads for geotile in eachrow(geotiles)
                println(geotile.id)
                #geotiles0 = geotiles.id[(geotiles.rgi1 .> 0) .& (geotiles.glacier_frac .> 0.3), :]
                #geotile = first(geotiles0)

                dh0 = dh1[mission_ref_fac][At(geotile.id),:,:]
                fac0 = fac[At(geotile.id),:,:]
                smb0 = smb[At(geotile.id),:,:]

                g0 = gemb[1]

                # round to the nearest day
                tdh = round.(Altim.decimalyear.(lookup(dh0, :date)), digits = 3);

                valid = .!isnan.(dh0)
                if !any(valid)
                    continue 
                end
                trange_dh, hrang_dh = Altim.validrange(.!isnan.(dh0))
                e_dh = extrema(tdh[trange_dh])

                tg = round.(vec(g0["datetime"]), digits = 3);

                _, trange_g = Altim.validrange(.!isnan.(g0["fac"]))
                e_g = extrema(tg[trange_g])

                ext = (max(e_dh[1], e_g[1]), min(e_dh[2], e_g[2]))

                trange_dh, = Altim.validrange((tdh .>= ext[1]) .& (tdh .<= ext[2]))
                trange_g, = Altim.validrange((tg .>= ext[1]) .& (tg .<= ext[2]))
                (length(trange_dh) != length(trange_g)) && error("gemb and dh do not have the same length")

                dh_altim = dh0[trange_dh,hrang_dh]
                dh_altim0 = copy(dh_altim)
                tdh0 = vec(tdh[trange_dh] .- mean(tdh[trange_dh]))

                # remove 3rd polynomal
                for h in lookup(dh_altim, :height)
                    fit0 = curve_fit(model2, tdh0, dh_altim[:,At(h)], p2)
                    dh_altim0[:,At(h)] = fit0.resid
                end

                #in geotile
                ext = GeoTiles.extent(geotile.id);
                # buffer by 1 degree all arround
                ext = Extents.buffer(ext, (X=geotile_buffer, Y=geotile_buffer))
                in_geotile = findall(vec([Altim.within(ext,x,y) for (x,y) in zip(g0["longitude"], g0["latitude"])]))

                dgembts = Dim{:gembts}(1:sum(in_geotile))
                dheight = dims(dh_altim, :height)
                p_scale = getindex.(gemb,"precipitation_scale")
                e_delta = getindex.(gemb,"elevation_delta")

                dgembrun = Dim{:gembrun}(1:length(p_scale))
                rmse = fill(NaN, (dgembrun, dgembts, dheight))

                for r in lookup(rmse, :gembrun)
                    g0 = gemb[r]

                    dh_gem = ((g0["smb"][in_geotile, trange_g] .* (1/.91)) .+ g0["fac"][in_geotile, trange_g])'
                    dh_gem0 = copy(dh_gem)
                    tg0 = vec(tdh[trange_g] .- mean(tg[trange_g]))
                    for j in axes(dh_gem, 2)
                        fit0 = curve_fit(model2, tdh0,  dh_gem[:,j], p2)
                        dh_gem0[:,j] = fit0.resid
                    end

                    for h in lookup(dh_altim, :height)
                        alt = dh_altim0[:,At(h)]
                        for j in axes(dh_gem0, 2)
                            rmse[At(r),At(j), At(h)] = sqrt(mean((alt .- dh_gem0[:,j]).^2))
                        end
                    end
                end

                # loop for each height range and find best run and repective elevation
                gembrun_best_fit = fill(0, (dheight))
                gembts_best_fit = fill(0, (dheight))
                rmse_best_fit = fill(0., (dheight))

                for h in lookup(dh_altim, :height)
                    foo = rmse[:,:,At(h)]
                    valid = .!isnan.(foo)
                    
                    if !any(valid)
                        continue
                    end

                    mrmse = minimum(foo[valid])
                    idx = findfirst(foo .== minimum(mrmse))

                    rmse_best_fit[At(h)] = mrmse
                    gembrun_best_fit[At(h)] = lookup(rmse, :gembrun)[idx[1]]
                    gembts_best_fit[At(h)] =  lookup(rmse, :gembts)[idx[2]]

                    fac0[:, At(h)] =  gemb[gembrun_best_fit[At(h)]]["fac"][in_geotile[gembts_best_fit[At(h)]],:]
                    smb0[:, At(h)] =  gemb[gembrun_best_fit[At(h)]]["smb"][in_geotile[gembts_best_fit[At(h)]], :]
                end

                fac[At(geotile.id),:,:] = fac0
                smb[At(geotile.id),:,:] = smb0
            end

            save(gemb_fit_file, Dict("fac" => fac, "smb" => smb,));
        end
        
        surface_mask = :glacier
        geotiles = geotiles_mask_hyps(surface_mask, geotile_width)
        # make geotile rgi regions mutually exexclusive 
        geotiles, reg = geotiles_mutually_exclusive_rgi!(geotiles)    
        geotiles = geotiles[(geotiles.glacier_frac .> 0), :]

        fac = Dict(mission_ref_fac => load(gemb_fit_file, "fac"));
        fac[mission_ref_fac][isnan.(fac[mission_ref_fac])] .= 0

        fac_nobs = deepcopy(fac)
        fac_nobs[mission_ref_fac] .= 1

        facv, nobsv = hyps_volume_change(fac, fac_nobs, geotiles; mask=:glacier)
        facv_reg = hyps_geotile_aggrigate(facv, geotiles, reg; fun = sum)

        smb = Dict(mission_ref_fac => load(gemb_fit_file, "smb"));
        smb[mission_ref_fac][isnan.(smb[mission_ref_fac])] .= 0
        
        smbv, nobsv = hyps_volume_change(smb, fac_nobs, geotiles; mask=:glacier)
        smbv_reg = hyps_geotile_aggrigate(smbv, geotiles, reg; fun = sum)
    end

    # compute volume and mass change for each region
    ## check if file exists and contains all binned_filled_files (need to load outfile to check this condition) 
    for surface_mask in surface_masks[1:end]

        fName = "dvdm_reg_$(surface_mask)_$(project_id)_p$(paramater_set).jld2"
        outfile = joinpath(binned_folder, fName)

        compute_dvdm_reg = false
        
        if isfile(outfile)
            df = load(outfile,"df")

            for binning_method = binning_methods
                for dem_id in dem_ids 
                    for curvature_correct in curvature_corrects
                        for amplitude_correct = amplitude_corrects

                            binned_filled_file, figure_suffix = binned_filled_filepath(binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)

                            if !(isfile(binned_filled_file))
                                continue
                            end

                            if .!any((df.binning_method .== binning_method) .& (df.dem_id .== string(dem_id)) .& (df.curvature_correct .== curvature_correct) .& (df.amplitude_correct .==  amplitude_correct))
                                compute_dvdm_reg = true;
                            end
                        end
                    end
                end
            end
        else
            compute_dvdm_reg = true;
        end
            
        if compute_dvdm_reg
            # initialize
            df = DataFrame()
            dm_reg_masked_all = []
            dv_reg_masked_all = []
            nobs_reg_masked_all = []

            geotiles = geotiles_mask_hyps(surface_mask, geotile_width)

            # make geotile rgi regions mutually exexclusive 
            geotiles, reg = geotiles_mutually_exclusive_rgi!(geotiles)

            geotiles = geotiles[(geotiles.glacier_frac.>0), :]

            drgi = Dim{:rgi}(reg)
            area_reg = DimArray(fill(0.0, length(reg)), drgi)
            for rgi in reg
                rgi_ind = geotiles[:, rgi] .> 0
                area_reg[At(rgi)] = sum(sum(geotiles[rgi_ind, "$(surface_mask)_area_km2"]))
            end

            # calculate regional volume and mass changes
            # this takes about 10 min for all permutations and combinations
            @time for binning_method = binning_methods
                for dem_id in dem_ids 
                    for curvature_correct in curvature_corrects
                        for amplitude_correct = amplitude_corrects

                            if false
                                binning_method = "meanmadnorm3"
                                dem_id = :cop30_v2 
                                curvature_correct = true
                                amplitude_correct = true
                            end

                            # paths to files
                            binned_file = binned_filepath(binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct)
                            binned_filled_file, figure_suffix = binned_filled_filepath(binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)

                            if !(isfile(binned_filled_file))
                                continue
                            end

                            model_param = load(binned_filled_file, "model_param")
                            dh1 = load(binned_filled_file, "dh_hyps");
                            nobs1 = load(binned_filled_file, "nobs_hyps");

                            # apply offset to dh
                            apply_landoffset = true;
                            if apply_landoffset
                                # remove only GEDI trend
                                landfit[At(["icesat", "icesat2", "hugonnet"]),:,:] .= 0
                                dh1 = landoffset!(dh1, geotiles, landfit)
                            end

                            # calculate geotile volume change
                            dv, nobs2 = hyps_volume_change(dh1, nobs1, geotiles; mask = surface_mask)

                            dv_reg = hyps_geotile_aggrigate(dv, geotiles, reg; fun = sum)
                            nobs_reg = hyps_geotile_aggrigate(nobs2, geotiles, reg; fun = sum)

                            # allign to mission_reference
                            if regional_offset_apply 
                                dv_offset, dv_reg = hyps_offset2reference!(dv_reg, mission_reference)
                            end

                            # align region 9 to icesat
                            regional_offset_apply_rgi9icesat  = false
                            if regional_offset_apply_rgi9icesat 
                                dv_offsetX, dv_regX = hyps_offset2reference!(copy(dv_reg), "icesat")

                                dv_offset[:,At("rgi9")] .= dv_offsetX[:,At("rgi9")]
                                dv_reg[:,At("rgi9"),:] .= dv_regX[:,At("rgi9"),:]
                            end

                            ## scale fac
                            if fac_scale_apply 
                                fac_scale0, facv_reg0, smbv_reg0  = hyps_scale_fac!(copy(facv_reg[At(mission_ref_fac),:,:]), copy(smbv_reg[At(mission_ref_fac),:,:]), dv_reg)
                            end

                            # calculate regional mass change 
                            dm_reg = hyps_mass_change(dv_reg, facv_reg0)

                            # start at zero mass and volume change
                            ref_mission = "hugonnet"
                            ref_epoch = findfirst(vec(any(.!isnan.(dm_reg[At(ref_mission),:,:]), dims = (1))))
                            for rgi in reg
                                dm_reg[:, At(rgi), :] .-= dm_reg[At(ref_mission), At(rgi), ref_epoch]
                                dv_reg[:, At(rgi), :] .-= dv_reg[At(ref_mission), At(rgi), ref_epoch]
                            end
                
                            # exclude estimates if they have insuffient observations
                            dm_reg_masked = copy(dm_reg)
                            dv_reg_masked = copy(dv_reg)
                            nobs_reg_masked = copy(nobs_reg)

                            for rgi in dims(dm_reg, :rgi)
                                for mission in dims(dm_reg, :mission)
                                    obs_per_km2 = nobs_reg[At(mission),At(rgi), :] ./ area_reg[At(rgi)]
                                    
                                    valid = obs_per_km2 .> 0
                                    if !any(valid)
                                        dm_reg_masked[At(mission),At(rgi), :] .= NaN
                                        dv_reg_masked[At(mission),At(rgi), :] .= NaN
                                        nobs_reg_masked[At(mission),At(rgi), :] .= 0
                                        continue 
                                    end

                                    obs_per_km2 = median(obs_per_km2[valid])

                                    if obs_per_km2 < .01
                                        dm_reg_masked[At(mission),At(rgi), :] .= NaN
                                        dv_reg_masked[At(mission),At(rgi), :] .= NaN
                                        nobs_reg_masked[At(mission),At(rgi), :] .= 0
                                        continue 
                                    end
                                end
                            end
                                
                            # add trend, acceleration and amplitude to data frame
                            ddate = dims(dm_reg_masked, :date)
                            dm_global = DimArray(zeros(length(ddate)), (ddate,))
                            dv_global = copy(dm_global)

                            date_intercept = 2012.;
                
                            for rgi in reg
                                dm_gt, dm_fit = region_fit(dm_reg_masked[:, At(rgi), :], date_intercept)
                                dm_global .+= dm_gt
                            
                                dv_km3, dv_fit = region_fit(dv_reg_masked[:, At(rgi), :], date_intercept)
                                dv_global .+= dv_km3

                                trend_gtyr = dm_fit.param[2]
                                acc_gty2 = dm_fit.param[3]
                                amp_gt = dm_fit.param[4]
                                err_2sigma_gt = std(dm_fit.resid) *2;

                                trend_km3yr = dv_fit.param[2]
                                acc_km3y2 = dv_fit.param[3]
                                amp_km3 = dv_fit.param[4]

                                dm_gt = [dm_gt.data]
                                dv_km3 = [dv_km3.data]
                                err_2sigma_km3 = std(dv_fit.resid) *2;
                                    
                                # include icesat - hugonnet offset
                                dm1 = dm_reg_masked[At("icesat"), At(rgi), :]
                                dm2 = dm_reg_masked[At("hugonnet"), At(rgi), :]
                                valid = .!isnan.(dm1) .& .!isnan.(dm2)
                                if any(valid)
                                    icesat_hugonnet_offset_Gt = median(abs.(dm1[valid] .- dm2[valid]))
                                else
                                    icesat_hugonnet_offset_Gt = NaN
                                end

                                df = append!(df, DataFrame(; binning_method, dem_id, curvature_correct, amplitude_correct, rgi, trend_gtyr, acc_gty2, amp_gt, err_2sigma_gt, trend_km3yr, acc_km3y2, amp_km3, dm_gt, dv_km3, err_2sigma_km3, icesat_hugonnet_offset_Gt))
                            end

                            # add global 
                            _, dm_fit = region_fit(dm_global, date_intercept)
                            _, dv_fit = region_fit(dv_global, date_intercept)

                            trend_gtyr = dm_fit.param[2]
                            acc_gty2 = dm_fit.param[3]
                            amp_gt = dm_fit.param[4]
                            err_2sigma_gt = std(dm_fit.resid) *2;

                            trend_km3yr = dv_fit.param[2]
                            acc_km3y2 = dv_fit.param[3]
                            amp_km3 = dv_fit.param[4]
                            err_2sigma_km3 = std(dv_fit.resid) *2;

                            rgi = "global"
                            dm_gt = [dm_global.data]
                            dv_km3 = [dv_global.data]
                            
                            # include icesat - hugonnet offset
                            dm1 = dm_reg_masked[At("icesat"), :, :]
                            dm2 = dm_reg_masked[At("hugonnet"), :, :]
                            valid = .!isnan.(dm1) .& .!isnan.(dm2)
                            if any(valid)
                                icesat_hugonnet_offset_Gt = mean(abs.(dm1[valid] .- dm2[valid]))
                            else
                                icesat_hugonnet_offset_Gt = NaN
                            end

                            df = append!(df, DataFrame(; binning_method, dem_id, curvature_correct, amplitude_correct, rgi, trend_gtyr, acc_gty2, amp_gt, err_2sigma_gt, trend_km3yr, acc_km3y2, amp_km3, dm_gt, dv_km3, err_2sigma_km3, icesat_hugonnet_offset_Gt))

                            dm_reg_masked_all = push!(dm_reg_masked_all, figure_suffix => dm_reg_masked)
                            dv_reg_masked_all = push!(dv_reg_masked_all, figure_suffix => dv_reg_masked)
                            nobs_reg_masked_all = push!(nobs_reg_masked_all, figure_suffix => nobs_reg_masked)
                        end
                    end
                end
            end

            # save output
            save(outfile, Dict("dm" => Dict(dm_reg_masked_all...), "dv" => Dict(dv_reg_masked_all...), "nobs" => Dict(nobs_reg_masked_all...), "df" => df, "area" => area_reg));
        end
    end
end