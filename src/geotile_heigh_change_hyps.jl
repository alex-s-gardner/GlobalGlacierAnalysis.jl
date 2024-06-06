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

    include("utilities_hyps.jl")

    project_id = :v01;
    geotile_width = 2;
    fill_dh = false;
    showplots = false;



    regional_offset_apply = true;
    fac_scale_apply = true

    mission_reference = "icesat2"
    
    param = (
        # filter parameters 

        #=
        # scale bincount_min by the raw std of binned values
        bincount_min = Dict("icesat" => 11, 
            "icesat2" => 21, 
            "gedi" => 21,
            "hugonnet" => 51,
        ),

        smooth_n = Dict("icesat" => 3,
            "icesat2" => 9,
            "gedi" => 9,
            "hugonnet" => 13,
        ),
        =#

        bincount_min = Dict("icesat" => 9, 
            "icesat2" => 9, 
            "gedi" => 9,
            "hugonnet" => 51,
        ),

        smooth_n = Dict("icesat" => 5,
            "icesat2" => 5,
            "gedi" => 5,
            "hugonnet" => 15,
        ),

        smooth_h2t_length_scale = 800, # 800 m = 1 year in distance for anomaly from variogram analysis =
        model1_madnorm_max = 5, # this is a sigma-equivelent threshold
    )

    binned_folder = analysis_paths(; geotile_width).binned
    fig_folder = joinpath(binned_folder, "figures")
    !isdir(fig_folder) && mkdir( fig_folder)

    mask  = :glacier;
    gt_file = joinpath(binned_folder, "geotile_$(mask)_hyps.arrow");
    geotiles = DataFrame(Arrow.Table(gt_file));
    geotiles.extent = Extent.(getindex.(geotiles.extent, 1));
    geotiles = copy(geotiles);

    # make geotile rgi regions mutually exexclusive
    geotiles, reg = geotiles_mutually_exclusive_rgi!(geotiles)    

    showplots && (using Plots)

    dem_ids = [:best, :cop30_v2]
    binning_methods = ["median", "meanmadnorm3", "meanmadnorm5"]
    curvature_corrects = [true, false]
    amplitude_corrects = [true, false]

    if false
        dem_ids = [:best]
        binning_methods = ["meanmadnorm3"]
        curvature_corrects = [true]
        amplitude_corrects = [true]
    end

    if mask == :land
         dem_ids = [:best]
        binning_methods = ["meanmadnorm3"]
        curvature_corrects = [true]
        amplitude_corrects = [false]
    end
        
    #plotly()
end


# ~5 min for filling and plotting per iteration, 1.5 hours for all permutations and combinations
for binning_method = binning_methods
    for dem_id in dem_ids 
       for curvature_correct in curvature_corrects
            for amplitude_correct = amplitude_corrects

                if false
                    binning_method = "meanmadnorm3"
                    dem_id = :best 
                    curvature_correct = true
                    amplitude_correct = false
                end

                # bin method
                if curvature_correct
                    runid = "$(mask)_dh_$(dem_id)_cc_$(binning_method)_$(project_id)"
                else
                    runid = "$(mask)_dh_$(dem_id)_$(binning_method)_$(project_id)"
                end

                binned_file = joinpath(binned_folder, "$(runid).jld2");

                if amplitude_correct
                    binned_filled_file = joinpath(binned_folder, "$(runid)_filled_ac.jld2")
                    out_id = replace(runid, "dh" => "dv")
                    figure_suffix = "$(out_id)_ac"
                else
                    binned_filled_file = joinpath(binned_folder, "$(runid)_filled.jld2")
                    figure_suffix = replace(runid, "dh" => "dv")
                end

                #fill_dh = true; showplots = true; (using Plots)
                if fill_dh || !isfile(binned_filled_file)
                    dh1 = load(binned_file, "dh_hyps");
                    nobs1 = load(binned_file, "nobs_hyps");

                    # align geotile dataframe with DimArrays
                    gt = collect(dims(dh1[first(keys(dh1))], :geotile))
                    gt_ind = [findfirst(geotiles.id .== g) for g in gt]
                    geotiles = geotiles[gt_ind,:]

                    # create a data frame to store model parameters
                    # initialize dimensional arrays
                    params = Dict();
                    for  mission in keys(dh1)
                        n = length(dims(dh1[mission], :geotile))
                        push!(params, mission => DataFrame(geotile = val(dims(dh1[mission], :geotile)), nobs_raw = zeros(n), nbins_raw = zeros(n), nobs_final = zeros(n), nbins_filt1 = zeros(n), param_m1 = [fill(NaN, size(p1)) for i in 1:n], h0 = fill(NaN,n), t0 = fill(NaN,n), dh0 = fill(NaN,n), bin_std = fill(NaN,n), bin_anom_std =  fill(NaN,n)));
                    end

                    showplots && (i = 652);
                    showplots && plot_height_time(dh1; geotile = geotiles[i,:], fig_suffix = "raw", fig_folder, figure_suffix, mask)

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

                        fname = joinpath(fig_folder, "variogram_range_ratio_$(figure_suffix).png")
                        save(fname, f)
                        display(f)
                    end

                    # filter and fit model to geotiles [4 min for for all glacierized geotiles and 4 missions]        
                    showplots && plot_height_time(dh1; geotile = geotiles[i,:], fig_suffix = "filled", fig_folder, figure_suffix, mask)
                    
                    # apply seasonal amplitude normalization
                    if amplitude_correct
                        hyps_amplitude_normalize!(dh1, params; mission_reference)
                    end
                    
                    showplots && plot_height_time(dh1; geotile = geotiles[i,:], fig_suffix = "filled_ampnorm", fig_folder, figure_suffix, mask)

                    hyps_fill_empty!(dh1, params, geotiles; mask)

                    hyps_fill_updown!(dh1, geotiles; mask)

                    showplots && plot_height_time(dh1; geotile = geotiles[i,:], fig_suffix = "filled_updown", fig_folder, figure_suffix, mask)

                    # save filled geotiles
                    save(binned_filled_file, Dict("dh_hyps" => dh1, "nobs_hyps" => nobs1, "model_param" =>   params ));
                end
            end
        end
    end
end


# find offsets to low vegitated land surfaces

mask_land = :land
binning_method = "meanmadnorm3"
dem_id = :best 
curvature_correct = true
amplitude_correct = false

# bin method
if curvature_correct
    runid = "$(mask_land)_dh_$(dem_id)_cc_$(binning_method)_$(project_id)"
else
    runid = "$(mask_land)_dh_$(dem_id)_$(binning_method)_$(project_id)"
end

landfit_param_file = joinpath(binned_folder, "$(runid)_fit_param.jld2");
if !isfile(landfit_param_file)

    binned_file = joinpath(binned_folder, "$(runid).jld2");

    if amplitude_correct
        binned_filled_file = joinpath(binned_folder, "$(runid)_filled_ac.jld2")
        figure_suffix = "$(runid)_ac";
    else
        binned_filled_file = joinpath(binned_folder, "$(runid)_filled.jld2")
        figure_suffix = runid;
    end

    model_param = load(binned_filled_file, "model_param")
    dh1 = load(binned_filled_file, "dh_hyps");
    nobs1 = load(binned_filled_file, "nobs_hyps");

    drgi = dims(dv_reg, :rgi)
    ddate = dims(dv_reg, :date)
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

        fname = joinpath(fig_folder, "$(figure_suffix)_$(rgi).png")
        save(fname, f)
        showplots && display(f)
    end

    save(landfit_param_file, Dict("landfit" => fit_param));
else
    landfit = load(landfit_param_file, "landfit");
end


begin
    # firn correction [same for all iterations so simply calculate once]
    binning_method = "meanmadnorm3"
    dem_id = :best
    curvature_correct = true
    amplitude_correct = true

    # bin method
    if curvature_correct
        runid = "$(mask)_dh_$(dem_id)_cc_$(binning_method)_$(project_id)"
    else
        runid = "$(mask)_dh_$(dem_id)_$(binning_method)_$(project_id)"
    end

    binned_file = joinpath(binned_folder, "$(runid).jld2")

    if amplitude_correct
        binned_filled_file = joinpath(binned_folder, "$(runid)_filled_ac.jld2")
        out_id = replace(runid, "dh" => "dv")
        figure_suffix = "$(out_id)_ac"
    else
        binned_filled_file = joinpath(binned_folder, "$(runid)_filled.jld2")
        figure_suffix = replace(runid, "dh" => "dv")
    end

    dh1 = load(binned_filled_file, "dh_hyps")

    binning_method_gemb = "mean"
    runid_gemb = "$(mask)_gemb_$(binning_method_gemb)_$(project_id)"
    gemb_file = joinpath(binned_folder, "$(runid_gemb).jld2")

    fac = load(gemb_file, "fac_hyps")
    nobs_gemb = load(gemb_file, "nobs_hyps")
    smb = load(gemb_file, "smb_hyps")

    # align geotile dataframe with DimArrays
    gt = collect(dims(dh1[first(keys(dh1))], :geotile))
    gt_ind = [findfirst(geotiles.id .== g) for g in gt]
    geotiles = geotiles[gt_ind, :]

    facv_reg0, smbv_reg0 = hyps_fac_correction(fac, smb, nobs_gemb, dh1, geotiles, reg)

    drgi = Dim{:rgi}(reg)
    area_reg = DimArray(fill(0.0, length(reg)), drgi)
    for rgi in reg
        rgi_ind = geotiles[:, rgi] .> 0
        area_reg[At(rgi)] = sum(sum(geotiles[rgi_ind, "$(mask)_area_km2"]))
    end

    # initialize
    df = DataFrame()
    dm_reg_masked_all = []
    dv_reg_masked_all = []
    nobs_reg_masked_all = []
end;

# calculate regional volume and mass changes
# this takes about 4 min for all permutations and combinations
@time for binning_method = binning_methods
    for dem_id in dem_ids 
        for curvature_correct in curvature_corrects
            for amplitude_correct = amplitude_corrects

                if false
                    binning_method = "meanmadnorm3"
                    dem_id = :best 
                    curvature_correct = true
                    amplitude_correct = true
                end

                # bin method
                if curvature_correct
                    runid = "$(mask)_dh_$(dem_id)_cc_$(binning_method)_$(project_id)"
                else
                    runid = "$(mask)_dh_$(dem_id)_$(binning_method)_$(project_id)"
                end

                binned_file = joinpath(binned_folder, "$(runid).jld2");

                if amplitude_correct
                    binned_filled_file = joinpath(binned_folder, "$(runid)_filled_ac.jld2")
                    out_id = replace(runid, "dh" => "dv")
                    figure_suffix = "$(out_id)_ac"
                else
                    binned_filled_file = joinpath(binned_folder, "$(runid)_filled.jld2")
                    figure_suffix = replace(runid, "dh" => "dv")
                end

                model_param = load(binned_filled_file, "model_param")
                dh1 = load(binned_filled_file, "dh_hyps");
                nobs1 = load(binned_filled_file, "nobs_hyps");

                # apply offset to dh
                apply_landoffset = true;
                if apply_landoffset
                    dh1 = landoffset!(dh1, geotiles, landfit)
                end

                # calculate geotile volume change
                dv, nobs2 = hyps_volume_change(dh1, nobs1, geotiles; mask)

                dv_reg = hyps_geotile_aggrigate(dv, geotiles, reg; fun = sum)
                nobs_reg = hyps_geotile_aggrigate(nobs2, geotiles, reg; fun = sum)

                # allign to mission_reference
                #if regional_offset_apply
                #    dv_offset, dv_reg = hyps_offset2reference!(dv_reg, mission_reference)
                #end

                ## scale fac
                if fac_scale_apply 
                    fac_scale, facv_reg, smbv_reg  = hyps_scale_fac!(copy(facv_reg0), copy(smbv_reg0), dv_reg)
                end

                # calculate regional mass change 
                dm_reg = hyps_mass_change(dv_reg, facv_reg)

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
fName = "dvdm_reg_$(project_id).jld2"
outfile = joinpath(binned_folder, fName);
save(outfile, Dict("dm" => Dict(dm_reg_masked_all...), "dv" => Dict(dv_reg_masked_all...), "nobs" => Dict(nobs_reg_masked_all...), "df" => df, "area" => area_reg));