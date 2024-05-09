# add amplitude scaling and aplitude magnitude to image 
begin
    using Arrow
    using Altim
    using DataFrames
    using Extents
    using DimensionalData
    using DimensionalData.LookupArrays
    using Statistics
    using Dates
    using Plots
    using LsqFit
    using DataInterpolations
    using JLD2
    using FileIO
    using RollingFunctions
    using NearestNeighbors
    using Interpolations
    using ScatteredInterpolation
    using DataInterpolations
    using CairoMakie

    project_id = :v01;
    geotile_width = 2;
    fill_dh = false
    showplots = false;

    # Define model that will be fit to all data binned by hypsometry
    #model::Function = model(t, h; t_intercept=2010) = hcat(ones(size(t)), (t .- t_intercept), h, h.^2, cos.(2 * pi * t), sin.(2 * pi * t))
    model1::Function = model1(x, p) = 
        p[1] .+ 
        p[2] .* x[:, 1] .+ 
        p[3] .* x[:, 1] .^ 2 .+ 
        p[4] .* x[:, 2] .+ 
        p[5] .* x[:, 2] .^ 2 .+ 
        sin.(2 .* pi .* (x[:, 1] .+ p[6])) .* (p[7] .+ p[8] .* x[:, 2] .* p[9] .* x[:, 2].^2)
        
        p1 = zeros(9);
        lb1 = [-10.0, -3.0, -2.0, -0.05, -0.0001, -1.0, -7.0, -0.05, -.001];
        ub1 = [+10.0, +3.0, +2.0, +0.05, +0.0001, +1.0, +7.0, +0.05, +.001];

    # seasonal only 
    model1_seasonal::Function = model1_seasonal(x,p) = sin.(2 .* pi .* (x[:, 1] .+ p[6])) .* (p[7] .+ p[8] .* x[:, 2] .* p[9] .* x[:, 2] .^ 2)

    # including quadratic for seasonal does not improve std(anom)
    #(p[6] .* cos.(2 .* pi .* x[:, 1]) .+  p[7].* sin.(2 .* pi .* x[:, 1])) .* (1 .+ p[8] .* x[:, 2] .+ p[9] .* x[:, 2] .^ 2)
    #p1 = zeros(9);
    #lb1 = [-10., -3., -2., -.05, -0.0001, -10., -10., -0.01, -0.0001];
    #ub1 = [+10., +3., +2., +.05, +0.0001, +10., +10., +0.01, +0.0001];

    # model fit across all geotiles for a region for a given year
    model2::Function = model2(h, p) = p[1] .+ p[2] .* h .+  p[3] .* h.^2;
    p2 = zeros(3);
    lb2 = [-30., -.1, -.01];
    ub2 = [+30., +.1, .01];

    model3::Function = model3(t, p) = p[1] .+ p[2] .* t .+ p[3] .* t .^ 2 .+ p[4] .* sin.(2 .* pi .* (t .+ p[5]))
    p3 = zeros(5)

    binned_folder = analysis_paths(; geotile_width).binned
    fig_folder = joinpath(binned_folder, "figures")
    !isdir(fig_folder) && mkdir( fig_folder)


    mask_project  = :glacier;
    gt_file = joinpath(binned_folder, "geotile_$(mask_project )_hyps.arrow");
    geotiles = DataFrame(Arrow.Table(gt_file));
    geotiles.extent = Extent.(getindex.(geotiles.extent, 1));
    geotiles = copy(geotiles);

    # fill dh by region
    reg = names(geotiles);
    reg = reg[contains.(reg, "rgi")];
    rgi_id = replace.(reg,Ref("rgi" => ""));
    rgi_id = parse.(Int64, rgi_id);

    # make regions mutually exclusive by assigning geotiles to region of largest overlap
    for geotile in eachrow(geotiles[!,reg])
        maxind = findfirst(maximum(geotile) .== collect(geotile))
        for r in eachindex(reg)
            if r == maxind
                continue
            else
                geotile[r] = 0
            end
        end
    end

    function plot_dvdm(dv, dm, nobs; title=nothing, fontsize=18, colors=palette(:Set1_4, length(keys(dh1))))

        decdate_intercept = 2010;

        f = Figure(backgroundcolor=RGB(0.98, 0.98, 0.98), size=(1000, 700), fontsize=fontsize)

        ga = f[1:5, 1] = GridLayout()
        gb = f[6, 1] = GridLayout()

        if !isnothing(title)
            Label(ga[1, 1, Top()],
                title, valign=:bottom,
                font=:bold,
                padding=(0, 0, 5, 0)
            )
        end

        decdate = Altim.decimalyear.(dims(dv, :date).val)
        kwargs = (; xminorgridvisible=true, xminorticks=IntervalsBetween(5))
        axmain = Axis(ga[1, 1]; ylabel="mass/volume change [Gt/km³]", kwargs...)

        hidexdecorations!(axmain; grid=false, minorgrid=false)

        axbottom = Axis(gb[1, 1]; yticks=(1:4, dims(nobs, :mission).val), yminorgridvisible=true, yminorticks=IntervalsBetween(2), yminorgridwidth=3, kwargs...)
        hideydecorations!(axbottom; grid=true, minorgrid=false, ticklabels=false)

        linkxaxes!(axmain, axbottom)

        # fit model to average across all missions
        CairoMakie.xlims!(axbottom, 2000, 2024)
        if hasdim(dm, :mission)
            valid = .!isnan.(dm)
            dmm = sum(dm.data .* valid.data, dims=1) ./ sum(valid.data, dims=1)

            valid = vec(.!isnan.(dmm))
            dm_fit = curve_fit(model3, decdate[valid] .- decdate_intercept, dmm[valid], p3)

            for (i, mission) in enumerate(dims(dm,:mission))
                dv0 = dv[At(mission),:]
                dm0 = dm[At(mission),:]

                CairoMakie.lines!(axmain, decdate, vec(dv0), color=(colors[i], 0.2), linewidth=2)
                CairoMakie.lines!(axmain, decdate, vec(dm0), label="$(mission)", color=(colors[i], 1), linewidth=2)
            end
        else
            valid = .!isnan.(dm)
            dm_fit = curve_fit(model3, decdate[valid] .- decdate_intercept, dm[valid], p3)

            CairoMakie.lines!(axmain, decdate, vec(dv), color=(:black, 0.2), linewidth=2)
            CairoMakie.lines!(axmain, decdate, vec(dm), color=(:black, 1), linewidth=2)
        end


        text!(
            axmain, 0.70, 0.95,
            text=
"trend:             $(round(dm_fit.param[2], digits = 1)) Gt yr⁻¹
acceleration: $(round(dm_fit.param[3], digits = 1)) Gt yr⁻²
amplitude:     $(round(dm_fit.param[4], digits = 1)) Gt",
            font=:bold,
            align=(:left, :top),
            space=:relative,
            fontsize=fontsize
        )

        var0 = log.(nobs.data');
        var0[isinf.(var0)] .= NaN

        hm = CairoMakie.heatmap!(axbottom, decdate, 1:4, var0, colormap=Reverse(:magma))

        rowgap!(ga, 1)
        rowgap!(gb, 1)

        CairoMakie.Colorbar(f[6, 2], hm; label="log(count)"); #, vertical=false, bbox=axbottom.scene.px_area)

        return f
    end
end


#for dem_id in [:best, :cop30_v2]
dem_id = :cop30_v2
     #for binning_method = ["median", "meanmadnorm3", "meanmadnorm5"];
     binning_method = "meanmadnorm3"
         #for curvature_correct in [true, false];
            curvature_correct = true

            #for amplitude_correct = [true, false];
                amplitude_correct = true
                
                
                regional_offsets = true;

                # run parameters
                min_n = 5; # icesat is highly sensitive to this number

                # bin method
                if curvature_correct
                    runid = "glacier_dh_$(dem_id)_cc_$(binning_method)_$(project_id)"
                else
                    runid = "glacier_dh_$(dem_id)_$(binning_method)_$(project_id)"
                end

                param = (
                    # bin method
                    
                    runid = runid,

                    # filter parameters 

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

                    smooth_h2t_length_scale = 2000, # 300 m = 1 year in distance for anomaly
                    model1_madnorm_max = 5, # this is a sigma-equivelent threshold

                    model_fit_weight = 5, 
                    min_region_coverage = 0.5,
                )
              
                binned_file = joinpath(binned_folder, "$(param.runid).jld2");

                if amplitude_correct
                    binned_filled_file = joinpath(binned_folder, "$(param.runid)_filled_ac.jld2")
                    out_id = replace(param.runid, "dh" => "dv")
                    figure_suffix = "$(out_id)_ac"
                else
                    binned_filled_file = joinpath(binned_folder, "$(param.runid)_filled.jld2")
                    figure_suffix = replace(param.runid, "dh" => "dv")
                end

                dh1 = load(binned_file, "dh_hyps");
                nobs1 = load(binned_file, "nobs_hyps");

                #dhX = load(binned_file, "dh_hyps");

                # align geotile dataframe
                gt = collect(dims(dh1[first(keys(dh1))], :geotile))
                gt_ind = [findfirst(geotiles.id .== g) for g in gt]
                geotiles = geotiles[gt_ind,:]

                # create a data frame to store model parameters
                # initialize dimensional arrays
                params = Dict();
                for  mission in keys(dh1)
                    n = length(dims(dh1[mission], :geotile))
                    push!(params, mission => DataFrame(geotile = val(dims(dh1[mission], :geotile)), nobs_raw = zeros(n), nbins_raw = zeros(n), nobs_final = zeros(n), nbins_filt1 = zeros(n), param_m1 = [fill(NaN, size(p1)) for i in 1:n], h0 = fill(NaN,n), t0 = fill(NaN,n), bin_std = fill(NaN,n), bin_anom_std =  fill(NaN,n)));
                end

                #  <><><><><><><><><><><><><><><><><><> FOR TESTING  <><><><><><><><><><><><><><><><><><>
                # products = products[keys(products)[2:2]]
                #  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

                # filter and fill on individual geotiles [4 min for for all glacierized geotiles and 4 missions]
                # if fill_dh || !isfile(binned_filled_file)
                
                    t = Altim.decimalyear.(dims(dh1[first(keys(dh1))], :date))
                    t = repeat(t, 1, length(dims(dh1[first(keys(dh1))], :height)))

                    h = val(dims(dh1[first(keys(dh1))], :height))'
                    h = repeat(h, length(dims(dh1[first(keys(dh1))], :date)), 1)

                    # for mission in keys(dh1)

                        # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><>
                        mission = "icesat"
                        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
                        # Threads.@threads for geotile in dims(dh1[mission], :geotile)

                            # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><>
                            geotile =   "lat[+74+76]lon[+060+062]"
                            # geotile = first(dims(dh1[mission], :geotile))
                            
                            # geotile = geotiles[findfirst((geotiles.rgi1 .> 0.) .& (geotiles.glacier_frac .> 0.3)),:]
                            # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

                            # println("$mission: $geotile")

                            k = findfirst(params[mission].geotile .== geotile)
                            df = @view params[mission][k,:]
                    
                            dh0 = @view dh1[mission][At(geotile),:,:];
                            nobs0 = @view nobs1[mission][At(geotile),:,:];
                            df.nobs_raw = sum(nobs0)
                            df.nbins_raw = sum(nobs0.>0)

                            ###################################### FILTER 1 ################################
                            valid = .!isnan.(dh0) .& (nobs0 .> param.bincount_min[mission]) .& (abs.(dh0) .< 200)
                            ################################################################################
                            
                            dh0[.!valid] .= NaN
                            nobs0[.!valid] .= 0
                            
                            # if there are not enough points to fit a model the set all to NaNs
                            va = sum(valid)
                            if va <= (length(p1) + 2)
                                dh0[:] .= NaN
                                nobs0[:] .= 0
                                continue
                            end

                            # determine valid range of data
                            (rrange, crange) = Altim.validrange(valid);
                            dh0 = dh0[rrange,crange];
                            valid = valid[rrange,crange];
                            nobs0 = nobs0[rrange,crange];
                            t0 = t[rrange,crange];
                            h0 = h[rrange,crange];

                            showplots && Plots.heatmap(dh0)

                            # center predictors and observations
                            t0_mean = round(mean(t0)); # remove an exact integer to keep phase from moving
                            df.t0 = t0_mean;
                            t0 = t0 .- t0_mean; 
                            h0_mean = mean(h0);
                            df.h0 = h0_mean;
                            h0 = h0 .- h0_mean;
                            dh0_median = median(dh0[valid]) + 0.0000001; # add a small offset to prevent numerical instability
                            dh0 = dh0 .- dh0_median;

                            df.bin_std = std(dh0[valid])

                            # fit global model 
                            fit1=[]
                            try
                                fit1 = curve_fit(model1, hcat(t0[valid], h0[valid]), dh0[valid], nobs0[valid], p1; lower=lb1, upper=ub1)
                            catch 
                                fit1 = curve_fit(model1, hcat(t0[valid], h0[valid]), dh0[valid], nobs0[valid], p1)
                            end

                            dh0_mdl = model1(hcat(t0[valid], h0[valid]), fit1.param);
                            dh0_anom = dh0[valid] .- dh0_mdl;

                            ###################################### FILTER 2 ####################################
                            # filter model1_madnorm_max sigma outliers
                            valid[valid] = Altim.madnorm(dh0_anom) .<= param.model1_madnorm_max
                            vb = sum(valid)
                            df.nbins_filt1 = vb;

                            if vb <= (length(p1) + 2)
                                dh0[:] .= NaN
                                nobs0[:] .= 0
                                continue
                            end

                            if vb < va
                                dh0[.!valid] .= NaN
                                nobs0[.!valid] .= 0
                                fit1 = curve_fit(model1, hcat(t0[valid], h0[valid]), dh0[valid], nobs0[valid], p1; lower=lb1, upper=ub1)

                                dh0_mdl = model1(hcat(t0[valid], h0[valid]), fit1.param);
                                dh0_anom = dh0[valid] .- dh0_mdl;
                            end
                            ####################################################################################
                            df.bin_anom_std = std(dh0_anom)

                            # take the median of the x closest neighbors 
                            if sum(valid) < param.smooth_n[mission]
                                anom_smooth = zeros(length(dh0))
                            else
                                # scale height distance relative to time (i.e. length-scale)
                                pts = hcat(t0[valid], h0[valid]/param.smooth_h2t_length_scale)'
                                kdtree = KDTree(pts)
                                (idxs, _) = knn(kdtree, pts, param.smooth_n[mission])
                                anom0 = map(ind -> median(dh0_anom[ind]), idxs)

                                # extrema(anom0)
                                # interpolate anomalies using weighted distance (Shepard(2))
                                itp = ScatteredInterpolation.interpolate(Shepard(2), pts, anom0);
                                pts = hcat(t0[:], h0[:]/param.smooth_h2t_length_scale)'
                                anom_smooth = vec(evaluate(itp, pts))
                            end

                            # fill out valid range (no extraploation beyond (rrange,crange) of geotile with the model only
                            dh1[mission][At(geotile),rrange,crange] = model1(hcat(t0[:], h0[:]), fit1.param) .+ dh0_median .+ anom_smooth
                            nobs1[mission][At(geotile), rrange, crange] = nobs0;

                            showplots && Plots.heatmap(dh1[mission][At(geotile),rrange,crange])

                            # println("granule interp: $(mission) - $(geotile)")

                            # add final parameters to DataFrame
                            df.nobs_final = sum(nobs0);
                            df.param_m1 = fit1.param;
                        end
                    end

                    # apply seasonal amplitude normalization
                    mission_ref = "icesat2"
                    if amplitude_correct
                        for mission in keys(dh1)
                            
                            if mission == mission_ref
                                continue
                            end
                            
                            for geotile in dims(dh1[mission], :geotile)
                            #geotile = first(dims(dh1[mission], :geotile))
                                k = findfirst(params[mission].geotile .== geotile)
                                df0 = params[mission][k,:]
                                dfr = params[mission_ref][k, :]

                                if any(isnan.(df0.param_m1)) || any(isnan.(dfr.param_m1))
                                    continue
                                end
                        
                                dh0 = dh1[mission][At(geotile),:,:];
                                valid = .!isnan.(dh0)
                                
                                (rrange, crange) = Altim.validrange(valid);
                                dh0 = dh0[rrange,crange];
                                valid = valid[rrange,crange];
                                t0 = t[rrange,crange] .- df0.t0;
                                h0 = h[rrange,crange] .- df0.h0;

                                p0 = df0.param_m1
                                pr = dfr.param_m1
                                Δp = pr .- p0
                                Δp[1:6] = p0[1:6];
                                
                                dh1[mission][At(geotile), rrange, crange] = dh0[:] .+ model1_seasonal(hcat(t0[:], h0[:]), Δp)
                            end
                        end
                    end
                    # 

                    # fill missing values with dh hysometric mean for full region within valid range of data.
                    # extrapolate first and last values outside of valid elevation range for each year.
                    # 19s for all 4 missions for all glacierized geotiles
                    for mission in keys(dh1)

                        # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><
                        # mission = "icesat"
                        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
                        for rgi in reg
                            # <><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><
                            # rgi = "rgi2"
                            # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

                            rgi_ind = geotiles[:, rgi] .> 0;
                            geotile_ids = geotiles[rgi_ind, :].id;

                            # find date and hight limits
                            dh0 = dh1[mission][At(geotile_ids), :, :]
                            valid = .!isnan.(dh0)

                            if .!any(valid)
                                continue
                            end

                            (rrange, crange, zrange) = Altim.validrange(valid)
                            dh0 = dh0[:,crange, zrange]
                            nobs0 = nobs1[mission][At(geotile_ids), crange, zrange]
                            valid0 = valid[:,crange, zrange]

                            if showplots 
                                p = plot()
                                for i in 1:length(dh0[:,1,1])
                                    if any(valid0[i,:,:])
                                        p = plot(dh0[i,:,:])
                                        display(p)
                                    end
                                end
                            end
                            
                            area = reduce(hcat, geotiles[rgi_ind, :glacier_area_km2])[zrange,:]
                            area_total = sum(area, dims=1)
                        
                            h0 = val(dims(dh0, :height))'
                            h0 = repeat(h0, length(dims(dh0, :geotile)), 1)

                            # loop through each date and fill empty geotiles with second order polynomia fit to all data
                            for i in eachindex(dims(dh0, :date))
                                # ----------------------------------- FOR TESTING ------------------------------
                                # i = 71
                                #println("$(rgi): i = $(i)")
                                # ------------------------------------------------------------------------------

                                dh00 = @view dh0[:,i,:]
                                valid00 = valid0[:,i,:]
                                nobs00 = @view nobs0[:,i,:];
                                
                                if sum(area_total[vec(any(valid00, dims=2))]) ./ sum(area_total) .< param.min_region_coverage
                                    # println(i)
                                    dh00[:] .= NaN;
                                    nobs00[:] .= 0;
                                else
                                    if (sum(valid00) > (length(p2)+1))  &&  !all(valid00)
                                        dh_med = median(dh00[valid00])

                                        # fit 2nd order plynomial to all data within region
                                        # give a weighting of 5 to modeled values, this is needed to keep the fit in check
                                        w = nobs00 .+ param.model_fit_weight

                                        fit2 = curve_fit(model2, h0[valid00], (dh00[valid00] .- dh_med), w[valid00], p2)
                                        

                                        showplots && plot(h0', model2(h0, fit2.param)', legend= false; )
                                        showplots && plot(h0[valid00], dh00[valid00], legend= false; seriestype=:scatter)

                                        # replace missing data
                                        dh00[.!valid00] = model2(h0[.!valid00] .- dh_med, fit2.param) .+ dh_med
                                    end
                                end
                            end
                            dh1[mission][At(geotile_ids), crange, zrange] = dh0
                            nobs1[mission][At(geotile_ids), crange, zrange] = nobs0
                        end
                        
                        # extrapolate first and last values for unmeasured elevation bins in a given year
                        Threads.@threads for geotile in dims(dh1[mission], :geotile)
                            dh0 = dh1[mission][At(geotile),:,:];
                            valid = .!isnan.(dh0.data)
                            if .!any(valid)
                                continue
                            end

                            (crange, _) = Altim.validrange(valid)

                            for i in crange
                                if any(valid[i,:])
                                    f = findfirst(vec(valid[i, :]))
                                    l = findlast(vec(valid[i, :]))
                                    dh0[i, 1:f] .= dh0[i, f]
                                    dh0[i, l:end] .= dh0[i, l]
                                end
                            end
                            dh1[mission][At(geotile),:,:] = dh0   
                        end
                    end
                    save(binned_filled_file, Dict("dh_hyps" => dh1, "nobs_hyps" => nobs1, "model_param" =>   params ));
                else
                    model_param = load(binned_filled_file, "model_param")
                    dh1 = load(binned_filled_file, "dh_hyps");
                    nobs1 = load(binned_filled_file, "nobs_hyps");
                end

                # calculate volume change for each geotile and regional offsets for hugonnet dataset
                # 3s for all 4 missions for all glacierized geotiles
                begin 
                    (ngeotile, ndate, nheight) = size(dh1[first(keys(dh1))])

                    missions = collect(keys(dh1))
                    nmissions = length(missions)
                    dmission = Dim{:mission}(missions)
                    drgi = Dim{:rgi}(reg)

                    dv = DimArray(fill(NaN, nmissions, ngeotile, ndate), (dmission, dims(dh1[first(missions)], :geotile), dims(dh1[first(missions)], :date)))
                    dv_reg = DimArray(fill(NaN, nmissions, length(reg), ndate), (dmission, drgi, dims(dh1[first(missions)], :date)))
                    dv_offset = DimArray(fill(0.0, nmissions, length(reg)), (dmission, drgi,))

                    glacier_area = reduce(hcat, geotiles.glacier_area_km2[:])'
                    glacier_area = permutedims(repeat(glacier_area, 1, 1, ndate), (1, 3, 2))

                    for mission in keys(dh1)
                        dv[At(mission),:,:] = dropdims(sum(dh1[mission] .* glacier_area, dims=3), dims=3) / 1000
                    end

                    for mission in dims(dv_offset, :mission)
                        for rgi in dims(dv_offset, :rgi)
                            rgi_ind = geotiles[:, rgi] .> 0
                            geotile_ids = geotiles[rgi_ind, :].id
                            dv_reg[At(mission),At(rgi),:] = vec(sum(dv[At(mission), At(geotile_ids), :, 1], dims=1))
                        end
                    end

                    # align all to icesat_2
                    base_mission = "icesat2"

                    for mission in dims(dv_offset, :mission)
                        if mission ==  base_mission
                            continue
                        end
                        
                        for rgi in dims(dv_offset, :rgi)

                            dv1 = dv_reg[At(base_mission),At(rgi),:]
                            valid1 = .!isnan.(dv1)

                            dv2 = dv_reg[At(mission),At(rgi),:]
                            valid2 = .!isnan.(dv2)

                            overlap = valid1 .& valid2

                            if any(overlap)
                                dv_offset[At(mission),At(rgi),:] = median(dv1[overlap]) - median(dv2[overlap])
                            end
                        end
                    end
                end

                # firn correction
                begin
                    binning_method_gemb = "mean"; 
                    runid_gemb = "glacier_gemb_$(binning_method_gemb)_$(project_id)"
                    gemb_file = joinpath(binned_folder, "$(runid_gemb).jld2");

                    fac = load(gemb_file, "fac_hyps");
                    nobs_gemb = load(gemb_file, "nobs_hyps");
                    smb = load(gemb_file, "smb_hyps");

                    nrgi = length(reg);
                    ndate = length(dims(dh1[first(keys(dh1))], :date));
                    nheight = length(dims(dh1[first(keys(dh1))], :height));
                    h0 = collect(dims(dh1[first(keys(dh1))], :height));
                    facv_reg = DimArray(fill(0.0, nrgi, ndate, nheight), (rgi=reg, date=collect(dims(dh1[first(keys(dh1))], :date)), height=h0));
                    smbv_reg = DimArray(fill(0.0, nrgi, ndate, nheight), (rgi=reg, date=collect(dims(dh1[first(keys(dh1))], :date)), height=h0));

                    # map gemb date into dh date
                    d1 = Altim.decimalyear.(collect(dims(fac, :date)))
                    kdtree = KDTree(d1')
                    d2 = Altim.decimalyear.(collect(dims(dh1[first(keys(dh1))], :date)))
                    idxs, dists = nn(kdtree, d2')

                    for rgi in reg
                        #  <><><><><><><><><><><><><><><><><><> FOR TESTING  <><><><><><><><><><><><><><><><><><>
                        # rgi = "rgi10"
                        #  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

                        rgi_ind = geotiles[:, rgi] .> 0
                        geotile = geotiles[rgi_ind, :]

                        # remove climatolotgy (assume some rate of vertical flow (divergence))

                        fac1 = fac[At(geotile.id),idxs,:]
                        nobs_gemb1 = nobs_gemb[At(geotile.id),idxs,:]
                        smb1 = smb[At(geotile.id), idxs, :]
                        valid1 = .!isnan.(fac1)

                        if !any(valid1)
                            continue
                        end
                        area_rgi = vec(sum(reduce(hcat, geotile.glacier_area_km2), dims=2))

                        # average accross all geotiles
                        fac_foo = fac1 .* nobs_gemb1;
                        smb_foo = smb1 .* nobs_gemb1
                        fac_foo[.!valid1] .= 0;
                        smb_foo[.!valid1] .= 0;

                        nobs_gemb_foo = valid1 .* nobs_gemb1
                        fac_rgi = dropdims(sum(fac_foo, dims=1) ./ sum(nobs_gemb_foo, dims=1), dims=1)
                        smb_rgi = dropdims(sum(smb_foo, dims=1) ./ sum(nobs_gemb_foo, dims=1), dims=1)
                        
                        valid = .!isnan.(fac_rgi)

                        # remove first values for each elevation band to get fac anomaly
                        for i = 1:size(fac_rgi,2)
                            if any(valid[:, i])
                                f = findfirst(vec(valid[:, i]))
                                fac_rgi[:, i] .-= fac_rgi[f, i]
                            end
                        end

                        # take cumulitive sum of smb
                        for i = 1:size(smb_rgi, 2)
                            if any(valid[:, i])
                                smb_rgi[:, i] = cumsum(smb_rgi[:, i])
                            end
                        end

                        # set unmodeled elevations == to first and last modeled levelations
                        valid = .!isnan.(fac_rgi);
                        (crange,_) = Altim.validrange(valid)

                        for i in crange
                            if any(valid[i,:])
                                f = findfirst(vec(valid[i,:]))
                                l = findlast(vec(valid[i, :]))

                                fac_rgi[i, 1:f] .= fac_rgi[i,f]
                                fac_rgi[i, l:end] .= fac_rgi[i,f]

                                smb_rgi[i, 1:f] .= smb_rgi[i,f]
                                smb_rgi[i, l:end] .= smb_rgi[i,f]

                            end
                        end

                        # convert from m to km.^3
                        for r in eachrow(fac_rgi)
                            r[:] = r .* area_rgi / 1000
                        end

                        for r in eachrow(smb_rgi)
                            r[:] = r .* area_rgi ./ (1000^2)
                        end

                        # set unmodeled times == to first and last modeled levelations
                        valid = .!isnan.(fac_rgi)
                        (_, rrange) = Altim.validrange(valid)
                        for i in rrange
                            if any(valid[:,i])
                                f = findfirst(vec(valid[:, i]))
                                l = findlast(vec(valid[:, i]))

                                fac_rgi[1:f,i] .= fac_rgi[f,i]
                                fac_rgi[l:end,i] .= fac_rgi[l,i]

                                smb_rgi[1:f,i] .= smb_rgi[f,i]
                                smb_rgi[l:end,i] .= smb_rgi[l,i]

                            end
                        end

                        # fill gap if it exists
                        for r in eachrow(fac_rgi)
                            if any(isnan.(r))
                                valid = .!isnan.(r)
                                itp = DataInterpolations.LinearInterpolation(r[valid].data, h0[valid])
                                r[.!valid] = itp(h0[.!valid])
                            end
                        end

                        for r in eachrow(smb_rgi)
                            if any(isnan.(r))
                                valid = .!isnan.(r)
                                itp = DataInterpolations.LinearInterpolation(r[valid].data, h0[valid])
                                r[.!valid] = itp(h0[.!valid])
                            end
                        end

                        # interpolate to date of dh data 
                        facv_reg[At(rgi), :, :] = fac_rgi;
                        smbv_reg[At(rgi), :, :] = smb_rgi;
                    end
                end

                ## scale fac
                begin
                    mission_ref = "icesat2"
                    fac_scale= DimArray(fill(NaN, length(reg)), (rgi = reg,))

                    for rgi in dims(dv_offset, :rgi)
                        fac1 = sum(facv_reg[At(rgi), :, :], dims=2)
                        smb1 = sum(smbv_reg[At(rgi), :, :], dims=2)
                        geb_dv1 = (smb1 ./ 0.91) .+ fac1
                        dv1 = copy(dv_reg[At(mission_ref),At(rgi), :])

                        # last valid date of gemb that was not interpolated
                        lastvalid = dims(fac,:date)[findlast(any(.!isnan.(fac), dims = (1,3)))[2]]

                        # findoverlap
                        valid1 = .!isnan.(fac1)
                        valid2 = .!isnan.(dv1)
                        overlap = vec((valid1 .& valid2) .& (dims(dv1, :date) .<= lastvalid))
                        t = decimalyear.(dims(dv1, :date))
                        # find amplitude
                        dv_fit = curve_fit(model3, t[overlap], dv1[overlap], p3)

                        dgemb_fit = curve_fit(model3, t[overlap], geb_dv1[overlap], p3)

                        sf = dv_fit.param[4] / dgemb_fit.param[4];
                        if isinf(sf)
                            sf = 0
                        end
                        fac_scale[At(rgi)] = sf
                    end
                end

                # calculate regional mass change change number of observations
                begin
                    dh_keys = keys(dh1)

                    dm_reg = copy(dv_reg)
                    nobs_reg = copy(dv_reg)

                    for rgi in dims(dm_reg, :rgi)

                        rgi_ind = geotiles[:, rgi] .> 0
                        geotile_ids = geotiles[rgi_ind, :].id

                        dfac = sum(facv_reg[At(rgi), :, :], dims=2) .* fac_scale[At(rgi)]

                        for mission in dims(dm_reg, :mission)
                            
                            dv0 = @view dv_reg[At(mission),At(rgi), :]
                                if regional_offsets
                                    dv0 .+= dv_offset[At(mission),At(rgi), :]
                                end
                            dm_reg[At(mission),At(rgi), :]  = (vec(dv0) .- vec(dfac)) * 0.910;
                            nobs_reg[At(mission),At(rgi), :]  = vec(sum(nobs1[mission][At(geotile_ids),:,:], dims=(1,3)))
                        end
                    end
                end

                begin
                
                    colors = palette(:Set1_4, length(keys(dh1)))
                    fontsize = 18

                    dvX = copy(dv_reg[1, 1, :])
                    dvX .= 0.
                    dmX = copy(dvX)
                    nobsX = copy(nobs_reg[:, 1, :])

                    for rgi in reg
                        title = "Randolph Glacier Inventory: Region $(rgi[4:end])"
                        dv = dv_reg[:, At(rgi), :, :]
                        dm = dm_reg[:, At(rgi), :, :]
                        nobs1 = nobs_reg[:, At(rgi), :, :]

                        f = plot_dvdm(dv, dm, nobs1; title, fontsize, colors)

                        fname = joinpath(fig_folder, "$(rgi)_$(figure_suffix).png")
                        save(fname, f)
                        display(f)

                        valid = .!isnan.(dm)
                        dvX .+= vec(sum(dv.data .* valid.data, dims=1) ./ sum(valid.data, dims=1))
                        dmX .+= vec(sum(dm.data .* valid.data, dims=1) ./ sum(valid.data, dims=1))
                        nobsX .+= nobs1;

                        # -------------------------------------- FOR TESTING -------------------------------
                        #rgi = reg[11]
                        # ----------------------------------------------------------------------------------
                    end

                    title = "Global"
                    f = plot_dvdm(dvX, dmX, nobsX; title, fontsize, colors)

                    fname = joinpath(fig_folder, "global_$(figure_suffix).png")
                    save(fname, f)
                    display(f)
                end
            end
        end
    end
end
