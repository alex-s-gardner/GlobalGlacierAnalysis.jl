begin
    using Altim
    using FileIO
    using JLD2
    using DimensionalData
    using ColorSchemes
    using CairoMakie
    using NCDatasets
    using MAT
    using ImageShow
    using FastRunningMedian
    include("utilities_hyps.jl")

    showplots = true
    project_id = :v01;
    geotile_width = 2;

    binned_folder = analysis_paths(; geotile_width).binned
    fig_folder = joinpath(binned_folder, "figures")
end

paramater_sets = 1:5
surface_masks = [:glacier, :glacier_b1km, :land, :glacier_b10km]


#for     paramater_set in paramater_sets
    paramater_set = first(paramater_sets)
    #for surface_mask in surface_masks
    surface_mask = first(surface_masks)
        begin

            # path to file
            dvdm_reg_file = dvdm_reg_filepath(binned_folder, surface_mask, project_id, paramater_set)

            dv = load(dvdm_reg_file, "dv");
            dm = load(dvdm_reg_file, "dm");
            nobs0 = load(dvdm_reg_file, "nobs");
            area = load(dvdm_reg_file, "area");
            df = load(dvdm_reg_file, "df");

            binning_method = "meanmadnorm3"
            dem_id = :best
            curvature_correct = true
            amplitude_correct = true

            # path to files
            binned_filled_file, figure_suffix = binned_filled_filepath(binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)

            dm_reg_masked = dm[figure_suffix]
            dv_reg_masked = dv[figure_suffix]
            nobs_reg_masked = nobs0[figure_suffix]

            t = collect(dims(dm_reg_masked, :date))
            t = Altim.decimalyear.(t)

            colors = ColorSchemes.Set1_4.colors[1:length(dims(dm_reg_masked, :mission))];
            fontsize = 18
            date_intercept = 2012.;
        end

        if true
            dvX = copy(dv_reg_masked[1, 1, :])
            dvX .= 0.
            dmX = copy(dvX)
            nobsX = copy(nobs_reg_masked[:, 1, :])

            for rgi in lookup(dm_reg_masked, :rgi)
                title = "Randolph Glacier Inventory: Region $(rgi[4:end])"

                f = plot_dvdm(dv_reg_masked[:, At(rgi), :], dm_reg_masked[:, At(rgi), :], nobs_reg_masked[:, At(rgi), :]; title, fontsize, colors, date_intercept, area = area[At(rgi)])

                fname = joinpath(fig_folder, "$(figure_suffix)_$(rgi).png")
                save(fname, f)
                showplots && display(f)

                valid = .!isnan.(dm_reg_masked[:, At(rgi), :])
                dvX .+= vec(sum(dv_reg_masked[:, At(rgi), :] .* valid.data, dims=1) ./ sum(valid.data, dims=1))
                dmX .+= vec(sum(dm_reg_masked[:, At(rgi), :] .* valid.data, dims=1) ./ sum(valid.data, dims=1))
                nobsX .+= nobs_reg_masked[:, At(rgi), :];
            end

            title = "Global"

            f = plot_dvdm(dvX, dmX, nobsX; title, fontsize, colors, area=sum(area.data))

            fname = joinpath(fig_folder, "$(figure_suffix)_global.png")
            save(fname, f)
            showplots && display(f)

            for rgi in vcat(collect(lookup(dm_reg_masked, :rgi)), "global")

                if rgi == "global"
                    title = "Global"
                else
                    title = "Randolph Glacier Inventory: Region $(rgi[4:end])"
                end

                dm = df[(df.rgi.==rgi).&df.amplitude_correct, :]
                f = plot_dm_permutations(t, dm, title; colors=colorschemes[:tab20c])

                fname = joinpath(fig_folder, "$(figure_suffix)_permutations_$(rgi).png")
                save(fname, f)
                showplots && display(f)
            end
        end


        if surface_mask != :glacier
            continue
        end

        # plot using a density of 850 kg/m
        #if true
            δ = .850
            dm850X = copy(dv_reg_masked[1, 1, :])
            dm850X .= 0.0
            dmFacX = copy(dm850X)
            dvX = copy(dm850X)

            nobsX = copy(nobs_reg_masked[:, 1, :])

            #for rgi in lookup(dm_reg_masked, :rgi)
                rgi = first(lookup(dm_reg_masked, :rgi))

                title = "Randolph Glacier Inventory: Region $(rgi[4:end])"
                dv = dv_reg_masked[:, At(rgi), :]
                dmFac = dm_reg_masked[:, At(rgi), :]
                dm850 = dv .* δ

                dv_mean, dv_fit = region_fit(dv, date_intercept)
                dmFac_mean, dm_fit = region_fit(dmFac, date_intercept)
                δ_effective = round(Int, dm_fit.param[2]./dv_fit.param[2] * 1000)

                # deteremine monthly fac equivelent
                dFac_mean = dmFac_mean[1:end-1] .- dmFac_mean[2:end]
                dv_mean = dv_mean[1:end-1] .- dv_mean[2:end]

                f = dFac_mean./dv_mean * 1000
                
                #δ_effective_monthly = running_median(f, 3)      
  
                #hist(dFac_mean[:], bins=0.:100.:2000.)

                f = plot_dvdm(dm850, dmFac, nobs_reg_masked[:, At(rgi), :]; title, fontsize, colors, date_intercept, area=area[At(rgi)], dmdm_flag = true, δ_effective)

                fname = joinpath(fig_folder, "$(figure_suffix)_d850_$(rgi).png")
                save(fname, f)
                showplots && display(f)

                valid = .!isnan.(dm_reg_masked[:, At(rgi), :])
                dvX .+= vec(sum(dv .* valid.data, dims=1) ./ sum(valid.data, dims=1))
                dm850X .+= vec(sum(dm850 .* valid.data, dims=1) ./ sum(valid.data, dims=1))
                dmFacX .+= vec(sum(dmFac .* valid.data, dims=1) ./ sum(valid.data, dims=1))
                nobsX .+= nobs_reg_masked[:, At(rgi), :]
            end

            title = "Global"

            v_mean, dv_fit = region_fit(dvX, date_intercept)
            dmFac_mean, dm_fit = region_fit(dmFacX, date_intercept)
            δ_effective = round(Int, dm_fit.param[2] ./ dv_fit.param[2] * 1000)

            f = plot_dvdm(dm850X, dmFacX, nobsX; title, fontsize, colors, area=sum(area.data), dmdm_flag=true, δ_effective)

            fname = joinpath(fig_folder, "$(figure_suffix)_d850_global.png")
            save(fname, f)
            showplots && display(f)
        end


        begin
            # plot with GRACE data 
            grace0 = read_grace_rgi()

            df0 = df[(df.binning_method .== binning_method) .& (df.dem_id .== dem_id) .& (df.curvature_correct .== curvature_correct) .& (df.amplitude_correct .== amplitude_correct), :]
            xlim = (2000, 2024)

            for rgi in vcat(collect(dims(dm_reg_masked, :rgi)))

                #rgi = vcat(collect(dims(dm_reg_masked, :rgi)), "global")[2]
                
                if rgi == "rgi12"
                    continue
                end

                dateN0  = vec(Altim.datenum2date.(grace0[rgi]["datenum"]))
                dateN0  = Altim.decimalyear.(dateN0 )

                # create gap between missions for plotting
                dM0 = vec(grace0[rgi]["dM_gt"]);
                sigma0 = vec(grace0[rgi]["dM_sigma_gt"]);
                
                ind = findfirst(dateN0 .> 2018.);
                
                if ~isempty(ind)
                    dateN0 = vcat(dateN0[1:ind-1], 2018., 2018., dateN0[ind:end]);
                    dM0 = vcat(dM0[1:ind-1], NaN, NaN, dM0[ind:end]);
                    sigma0 = vcat(sigma0[1:ind-1], NaN, NaN, sigma0[ind:end]);
                end
                
                grace = Dict(
                    "date" => dateN0,
                    "dm_gt" => dM0,
                    "2sigma_gt" =>  sigma0.*2, 
                    "label" => "gravimetry"
                )

                altim = Dict(
                    "date" => t,
                    "dm_gt" => df0[df0.rgi .== rgi, :dm_gt][1],
                    "2sigma_gt" =>  sqrt.((df0[df0.rgi .== rgi, :err_2sigma_gt][1]^2) .+ (0.3.*(df0[df0.rgi .== rgi, :dm_gt][1] .- df0[df0.rgi .== rgi, :dv_km3][1])).^2),
                    "label" => "altimetry"
                ) 

                title = "Randolph Glacier Inventory: Region $(rgi[4:end])"

                area0 = area[At(rgi)]
                f = plot_dm([grace, altim]; title, fontsize=18, date_intercept=2012, area=sum(area[At(rgi)]), xlim, date_alignall=(2003, 2004), colors=colorschemes[:Set1_9], scale_color=false, stats_in_label=true,)

                fname = joinpath(fig_folder, "$(figure_suffix)_altim_grace_$(rgi).png")
                save(fname, f)
                showplots && display(f)
            end

            # plot global
            rgi = unique(df.rgi)
            rgi = rgi[(rgi.!="global") .& (rgi.!="rgi5") .& (rgi.!="rgi19") .& (rgi.!="rgi12")]

            dateN0  = vec(Altim.datenum2date.(grace0[rgi[1]]["datenum"]))
            dateN0  = Altim.decimalyear.(dateN0 )

            # create gap between missions for plotting
            dM0 = [vec(grace0[r]["dM_gt"]) for r in rgi];
            sigma0 = [vec(grace0[r]["dM_sigma_gt"]) for r in rgi];

            dM0 = vec(sum(hcat(dM0...), dims=2))
            sigma0 = vec(sqrt.(sum(hcat(sigma0...) .^ 2, dims=2)))
            ind = findfirst(dateN0 .> 2018.);

            if ~isempty(ind)
                dateN0 = vcat(dateN0[1:ind-1], 2018., 2018., dateN0[ind:end]);
                dM0 = vcat(dM0[1:ind-1], NaN, NaN, dM0[ind:end]);
                sigma0 = vcat(sigma0[1:ind-1], NaN, NaN, sigma0[ind:end]);
            end

            grace = Dict(
                "date" => dateN0,
                "dm_gt" => dM0,
                "2sigma_gt" =>  sigma0.*2,
                "label" => "gravimetry"
            )

            ia, ib = Altim.intersectindices(df0.rgi, rgi; bool=false)

            altim = Dict(
                "date" => t,
                "dm_gt" => vec(sum(reduce(hcat, df0[ia, :dm_gt]), dims=2)),
                "2sigma_gt" => sqrt.((sum(df0[ia, :err_2sigma_gt] .^ 2)) .+ (0.3 .* (vec(sum(reduce(hcat, df0[ia, :dm_gt]), dims=2)) .- vec(sum(reduce(hcat, df0[ia, :dv_km3]), dims=2)))) .^ 2),
                "label" => "altimetry"
            ) 

            title = "Global Excluding Greenland and Antarctica"

            area0 = sum(area[At(rgi)])
            f = plot_dm([grace, altim]; title, fontsize=18, date_intercept=2012, xlim, date_alignall=(2003, 2004), colors=colorschemes[:Set1_9], scale_color=false, stats_in_label=true)

            fname = joinpath(fig_folder, "$(figure_suffix)_global_altim_grace.png")
            save(fname, f)
            showplots && display(f)

            # which runs provide the best agreement between icesat and hugonnet?
            # 5 runs have < 50Gt disagreement, one is madnorm3, true, true, true
            valid = .!isnan.(df.icesat_hugonnet_offset_Gt) .& (df.rgi .== "global")
            hist(df.icesat_hugonnet_offset_Gt[valid], bins=0:1:200)
            idx = findall((df.icesat_hugonnet_offset_Gt .< 30) .& (df.rgi .== "global"))


            # compare rate over each decade to compare with Hugonnet
            didx1 = (altim["date"] .>= 2000) .& (altim["date"] .< 2010) .& .!isnan.(altim["dm_gt"])
            didx2 = (altim["date"] .>= 2010) .& (altim["date"] .< 2020) .& .!isnan.(altim["dm_gt"])

            modelfoo(t, p) = p[1] .+ p[2] .* t .+  p[3] .* sin.(2 .* pi .* (t .+ p[4]))
            pfoo = [0., 0., 0., 0.]

            dm_fit1 = curve_fit(modelfoo, (altim["date"][didx1] .- mean(altim["date"][didx1])), altim["dm_gt"][didx1], pfoo)
            dm_fit2 = curve_fit(modelfoo, (altim["date"][didx2] .- mean(altim["date"][didx2])), altim["dm_gt"][didx2], pfoo)

            println("Global Glacier Mass Change Rate Excluding Greenland and Antarctica")
            println("    2000-2010 = $(round(Int, dm_fit1.param[2])) Gt/yr")
            println("    2010-2020 = $(round(Int, dm_fit2.param[2])) Gt/yr")

        end

        # plot with zemp
        begin
            zemp19 = read_zemp2019()

            xlim = (2000, 2024)
            for rgi in vcat(collect(dims(dm_reg_masked, :rgi)))

                if rgi == "rgi12"
                    continue
                end

                date = t;
                didx = (date .>= xlim[1]) .& (date .<= xlim[2])
                date = date[didx]

                altim = Dict(
                    "date" => date,
                    "dm_gt" => df0[df0.rgi.==rgi, :dm_gt][1][didx],
                    "2sigma_gt" => sqrt.((df0[df0.rgi.==rgi, :err_2sigma_gt][1]^2) .+ (0.3 .* (df0[df0.rgi.==rgi, :dm_gt][1][didx] .- df0[df0.rgi.==rgi, :dv_km3][1][didx])) .^ 2),
                    "label" => "altimetry"
                )

                date = vec(Altim.decimalyear.(collect(dims(zemp19.dm_gt, :date))))
                didx = (date .>= xlim[1]) .& (date .<= xlim[2])
                date = date[didx]

                err0 = vec(zemp19.err_gt[At(rgi), didx] .^ 2)
                err0[isnan.(err0)] .= 0
                err0 = sqrt.(cumsum(err0))

                zemp = Dict(
                    "date" => date,
                    "dm_gt" => vec(zemp19.dm_gt[At(rgi), didx]),
                    "2sigma_gt" => err0,
                    "label" => "in situ"
                )

                title = "Randolph Glacier Inventory: Region $(rgi[4:end])"

                f = plot_dm([zemp, altim]; title, fontsize=18, date_intercept=2012, xlim, date_alignall=(2003, 2004), colors=colorschemes[:Set1_9], scale_color=false, stats_in_label=true,)

                fname = joinpath(fig_folder, "$(rgi)_altim_zemp.png")
                save(fname, f)
                showplots && display(f)
            end

            title = "Global Excluding Greenland and Antarctica"
            rgi = unique(df.rgi)
            rgi = rgi[(rgi.!="global").&(rgi.!="rgi5").&(rgi.!="rgi19").&(rgi.!="rgi12")]

            xlim = (2000, 2024)
            ia, ib = Altim.intersectindices(df0.rgi, rgi; bool=false)

            altim = Dict(
                "date" => t,
                "dm_gt" => vec(sum(reduce(hcat, df0[ia, :dm_gt]), dims=2)),
                "2sigma_gt" => sqrt.((sum(df0[ia, :err_2sigma_gt] .^ 2)) .+ (0.3 .* (vec(sum(reduce(hcat, df0[ia, :dm_gt]), dims=2)) .- vec(sum(reduce(hcat, df0[ia, :dv_km3]), dims=2)))) .^ 2),
                "label" => "altimetry"
            ) 

            date = vec(Altim.decimalyear.(collect(dims(zemp19.dm_gt, :date))));
            didx = (date .>= xlim[1]) .& (date .<= xlim[2])
            date = date[didx]

            err0 = vec(sqrt.(sum(zemp19.err_gt[At(rgi), didx] .^ 2, dims=:rgi)))
            err0[isnan.(err0)] .= 0
            err0 = sqrt.(cumsum(err0.^2))

            zemp = Dict(
                "date" => date,
                "dm_gt" => vec(sum(zemp19.dm_gt[At(rgi), didx], dims=1)),
                "2sigma_gt" => err0,
                "label" => "in situ"
            )

            f = plot_dm([grace, altim, zemp]; title, fontsize=18, date_intercept=2012, area=sum(area[At(rgi)]), xlim, date_alignall=(2003, 2004), colors=colorschemes[:Set1_9], scale_color=false, stats_in_label = true);
            fname = joinpath(fig_folder, "$(figure_suffix)_global_altim_grace_zemp.png")
            save(fname, f)
            showplots && display(f)
        end

        begin
            title = "Global"

            rgi = unique(df.rgi)
            rgi = rgi[(rgi.!="global")]

            xlim = (2000, 2024)
            ia, ib = Altim.intersectindices(df0.rgi, rgi; bool=false)

            altim = Dict(
                "date" => t,
                "dm_gt" => vec(sum(reduce(hcat, df0[ia, :dm_gt]), dims=2)),
                "2sigma_gt" => sqrt.((sum(df0[ia, :err_2sigma_gt] .^ 2)) .+ (0.3 .* (vec(sum(reduce(hcat, df0[ia, :dm_gt]), dims=2)) .- vec(sum(reduce(hcat, df0[ia, :dv_km3]), dims=2)))) .^ 2),
                "label" => "altimetry"
            )

            date = vec(Altim.decimalyear.(collect(dims(zemp19.dm_gt, :date))));
            didx = (date .>= xlim[1]) .& (date .<= xlim[2])
            date = date[didx]

            err0 = vec(sqrt.(sum(zemp19.err_gt[At(rgi), didx] .^ 2, dims=:rgi)))
            err0[isnan.(err0)] .= 0
            err0 = sqrt.(cumsum(err0 .^ 2))

            zemp = Dict(
                "date" => date,
                "dm_gt" => vec(sum(zemp19.dm_gt[At(rgi), didx], dims=1)),
                "2sigma_gt" => err0,
                "label" => "in situ"
            )

            f = plot_dm([altim, zemp]; title, fontsize=18, date_intercept=2012, area=sum(area[At(rgi)]), xlim, date_alignall=(2000, 2002), colors=colorschemes[:Set1_9], scale_color=false, stats_in_label=true);
            fname = joinpath(fig_folder, "$(figure_suffix)_global_altim_zemp.png")
            save(fname, f)
            showplots && display(f)
        end

        begin
            title = "Global: Marzeion 2012"
            marzeion12 = read_marzeion2012()
            scenarios = dims(marzeion12.dm_gt, :scenario)

            m12 = []
            xlim = (2000, 2100)
            for scenario in scenarios
                has_data = vec(any(.!isnan.(marzeion12.dm_gt[:, :, :, At(scenario)]), dims=(:rgi, :date )))

                date = vec(Altim.decimalyear.(collect(dims(marzeion12.dm_gt, :date))))
                didx = (date .>= xlim[1]) .& (date .<= xlim[2])
                date = date[didx]

                dm_gt = vec(median(sum(marzeion12.dm_gt[:, didx, has_data, At(scenario)], dims=:rgi), dims=:climate_model))

                # align all to start (zero spread in xlim[1])
                sigma2_gt = sum(marzeion12.dm_gt[:, didx, has_data, At(scenario)], dims=:rgi)
                for cm in dims(sigma2_gt, :climate_model)
                    Δdm = dm_gt[1] - sigma2_gt[:,:, At(cm)][1]
                    sigma2_gt[:,:, At(cm)] .+= Δdm
                end
                sigma2_gt = vec(std(sigma2_gt, dims=:climate_model))

                m12 = push!(m12, 
                    Dict(
                        "date" => date,
                        "dm_gt" => dm_gt,
                        "2sigma_gt" => sigma2_gt,
                        "label" => scenario
                        ))
            end
            f = plot_dm(vcat(m12, [altim]); title, fontsize=18, date_intercept=2012, xlim, stats_in_label=true);
            fname = joinpath(fig_folder, "$(figure_suffix)_global_altim_marzeion2012.png")
            save(fname, f)
            showplots && display(f)
        end

        begin
            title = "Global: Marzeion 2020"
            marzeion20 = read_marzeion2020()
            scenarios = dims(marzeion20.dm_gt, :scenario)
            gmodels = dims(marzeion20.dm_gt, :glacier_model)

            m20 =[]
            xlim = (2001, 2100)
            for scenario in scenarios
                
                date = Altim.decimalyear.(collect(dims(marzeion20.dm_gt, :date)));
                didx = (date .>= xlim[1]) .& (date .<= xlim[2])
                date = date[didx]

                for gmodel in gmodels
                    has_data_clim = vec(any(.!isnan.(marzeion20.dm_gt[:, didx, :, At(gmodel), At(scenario)]), dims=(:rgi, :date)))
                
                    if !(any(has_data_clim))
                        continue
                    end
                
                    dm_gt = vec(median(sum(marzeion20.dm_gt[:, didx, has_data_clim, At(gmodel), At(scenario)], dims=:rgi), dims=(:climate_model)))
                    if all(isnan.(dm_gt))
                        continue
                    end

                    # align all to start (zero spread in xlim[1])
                    sigma2_gt = sum(marzeion20.dm_gt[:, didx, has_data_clim, At(gmodel), At(scenario)], dims=:rgi)
                    for cm in dims(sigma2_gt, :climate_model)
                        Δdm = dm_gt[1] - sigma2_gt[:, :, At(cm)][1]
                        sigma2_gt[:, :, At(cm)] .+= Δdm
                    end
                    sigma2_gt = vec(std(sigma2_gt, dims=:climate_model))
                    
                    m20 = push!(m20,
                        Dict(
                            "date" => date,
                            "dm_gt" => dm_gt,
                            "2sigma_gt" => sigma2_gt,
                            "label" => "$scenario: $gmodel"
                        ))
                end
            end

            f = plot_dm(vcat(m20, [altim]); title, fontsize=18, date_intercept=2012, xlim, date_alignall=(2005, 2010), stats_in_label=true);
            fname = joinpath(fig_folder, "$(figure_suffix)_global_altim_marzeion2020.png")
            save(fname, f)
            showplots && display(f)
        end

        begin
            title = "Global: IPCC AR6"
            IPCCAR6 = read_ipccar6()
            scenarios = dims(IPCCAR6.dm_gt, :scenario)


            AR6 = []
            xlim = (2001, 2100)
            for scenario in scenarios  
                date = Altim.decimalyear.(collect(dims(IPCCAR6.dm_gt, :date)))
                didx = (date .>= xlim[1]) .& (date .<= xlim[2])
                date = date[didx]


                dm_gt = vec(sum(IPCCAR6.dm_gt[:, didx, At(scenario)], dims=:rgi))
                sigma2_gt =sqrt.(sum(IPCCAR6.err_gt[:, didx, At(scenario)].^2, dims=:rgi))
                
                AR6 = push!(AR6,
                    Dict(
                        "date" => date,
                        "dm_gt" => dm_gt,
                        "2sigma_gt" => sigma2_gt,
                        "label" => "$scenario"
                    ))
            end

            f = plot_dm(vcat(AR6, [altim]); title, fontsize=18, date_intercept=2020, xlim, date_alignall=(2014.5, 2015.5), stats_in_label=true)
            fname = joinpath(fig_folder, "$(figure_suffix)_global_altim_ipccar6_fig9_21.png")
            save(fname, f)
            showplots && display(f)
        end


        begin
            start2007 = true;
            title = "Global: IPCC AR6"
            IPCCAR6 = read_ipccar6(; start2007)
            scenarios = dims(IPCCAR6.dm_gt, :scenario)


            AR6 = []
            xlim = (2001, 2100)
            for scenario in scenarios  
                date = Altim.decimalyear.(collect(dims(IPCCAR6.dm_gt, :date)))
                didx = (date .>= xlim[1]) .& (date .<= xlim[2])
                date = date[didx]


                dm_gt = vec(sum(IPCCAR6.dm_gt[:, didx, At(scenario)], dims=:rgi))
                sigma2_gt =sqrt.(sum(IPCCAR6.err_gt[:, didx, At(scenario)].^2, dims=:rgi))
                
                AR6 = push!(AR6,
                    Dict(
                        "date" => date,
                        "dm_gt" => dm_gt,
                        "2sigma_gt" => sigma2_gt,
                        "label" => "$scenario"
                    ))
            end

            f = plot_dm(vcat(AR6, [altim]); title, fontsize=18, date_intercept=2020, xlim, date_alignall=(2006.5, 2007.5), stats_in_label=true)
            fname = joinpath(fig_folder, "$(figure_suffix)_global_altim_ipccar6_fig9_21_start2007.png")
            save(fname, f)
            showplots && display(f)
        end
    end
end




begin

project_id = :v01;
geotile_width = 2;


binned_folder = analysis_paths(; geotile_width).binned
fig_folder = joinpath(binned_folder, "figures")

binning_method = "meanmadnorm3"
dem_id = :best
curvature_correct = true
amplitude_correct = true

end
# path to files
rgi = "global"

rgi = "rgi2"
surface_mask = :glacier

dvdm_reg_file = dvdm_reg_filepath(binned_folder, surface_mask, project_id, 1)
dv = load(dvdm_reg_file, "dv");
reg = collect(dims(dv[first(keys(dv))], :rgi))
reg = push!(reg, "global")

for rgi in reg
    for paramater_set in 1:5
        binned_filled_file, figure_suffix = binned_filled_filepath(binned_folder, surface_mask, dem_id, binning_method, project_id, curvature_correct, amplitude_correct, paramater_set)
        fname = joinpath(fig_folder, "$(figure_suffix)_altim_grace_$(rgi).png")
        #fname = joinpath(fig_folder, "$(figure_suffix)_$(rgi).png")
        if isfile(fname)
        img = load(fname)
        display(img)
        end
end
end