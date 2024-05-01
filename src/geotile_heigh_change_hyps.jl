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

    project_id = :v01;
    geotile_width = 2;
    
    #binning_method = "median"; # "meanmadnorm3"
    binning_method = "meanmadnorm3";

    # run parameters
    min_n = 5; # icesat is highly sensitive to this number

    dem_id = :best # [:rema_v2_10m]#[:cop30_v2]; #, :nasadem_v1, :rema_v2_10m, :arcticdem_v4_10m, :best]
    curvature_correct = true;
    regional_offsets = false;

    # bin method
    if  curvature_correct
        runid = "glacier_$(dem_id)_dh_cc_$(binning_method)_$(project_id)"
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

        min_region_coverage = 0.5,
    )

    binned_folder = analysis_paths(; geotile_width).binned
    binned_file = joinpath(binned_folder, "$(param.runid).jld2");

    mask = :glacier;
    gt_file = joinpath(binned_folder, "geotile_$(mask)_hyps.arrow");
    geotiles = DataFrame(Arrow.Table(gt_file));
    geotiles.extent = Extent.(getindex.(geotiles.extent, 1));

    # fill dh by region
    reg = names(geotiles);
    reg = reg[contains.(reg, "rgi")];
    reg = reg[1:1]
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

    showplots = false;

    dh1 = load(binned_file, "dh_hyps");
    nobs1 = load(binned_file, "nobs_hyps");

    gt = collect(dims(dh1[first(keys(dh1))], :geotile))
    gt_ind = [findfirst(geotiles.id .== g) for g in gt]
    geotiles = geotiles[gt_ind,:]

    # create a data frame to store model parameters
    # initialize dimensional arrays
    params = Dict();
    for  mission in keys(dh1)
        n = length(dims(dh1[mission], :geotile))
        push!(params, mission => DataFrame(geotile = val(dims(dh1[mission], :geotile)), nobs_raw = zeros(n), nbins_raw = zeros(n), nobs_final = zeros(n), nbins_filt1 = zeros(n), param_m1 = [fill(NaN, size(p1)) for i in 1:n], bin_std = fill(NaN,n), bin_anom_std =  fill(NaN,n)));
    end

    #  <><><><><><><><><><><><><><><><><><> FOR TESTING  <><><><><><><><><><><><><><><><><><>
    # products = products[keys(products)[2:2]]
    #  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
end

# filter and fill on individual geotiles [4 min for for all glacierized geotiles and 4 missions]
begin
    @time for mission in keys(dh1)

        # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><>
        # mission = "icesat"
        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

        t = Altim.decimalyear.(dims(dh1[mission], :date))
        t = repeat(t, 1, length(dims(dh1[mission], :height)))
        
        h = val(dims(dh1[mission], :height))'
        h = repeat(h, length(dims(dh1[mission], :date)),1)

        Threads.@threads for geotile in dims(dh1[mission], :geotile)
            # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><>
            # geotile = first(dims(dh1[mission], :geotile))
            # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            k = findfirst(params[mission].geotile .== geotile)
            df = @view params[mission][k,:]
    
            dh0 = dh1[mission][At(geotile),:,:];
            nobs0 = nobs1[mission][At(geotile),:,:];
            df.nobs_raw = sum(nobs0)
            df.nbins_raw = sum(nobs0.>0)

            ###################################### FILTER 1 ####################################
            valid = .!isnan.(dh0) .& (nobs0 .> param.bincount_min[mission]);
            dh0[.!valid] .= NaN
            nobs0[.!valid] .= 0
            ####################################################################################
            
            if sum(valid) < (length(p1)+1)
                continue
            end

            # determine valid range of data
            (rrange, crange) = Altim.validrange(valid);
            dh0 = dh1[mission][At(geotile),rrange,crange];
            nobs0 = nobs1[mission][At(geotile),rrange,crange];
            t0 = t[rrange,crange];
            h0 = h[rrange,crange];

            # center predictors and observations
            t0_mean = round(mean(t0)); # remove an exact integer to keep phase from moving
            t0 = t0 .- t0_mean; 
            h0_mean = mean(h0);
            h0 = h0 .- h0_mean;

            valid = .!isnan.(dh0);

            if sum(valid) < (length(p1)+1)
                continue
            end

            dh0_median = median(dh0[valid]);
            dh0 = dh0 .- dh0_median;

            df.bin_std = std(dh0[valid])

            # fit global model 
            va = sum(valid)
            fit1 = curve_fit(model1, hcat(t0[valid], h0[valid]), dh0[valid], nobs0[valid], p1; lower=lb1, upper=ub1)
            dh0_mdl = model1(hcat(t0[valid], h0[valid]), fit1.param);
            dh0_anom = dh0[valid] .- dh0_mdl;

            ###################################### FILTER 2 ####################################
            # filter model1_madnorm_max sigma outliers
            valid[valid] = Altim.madnorm(dh0_anom) .<= param.model1_madnorm_max
            vb = sum(valid)
            df.nbins_filt1 = vb;

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
                anom_smooth = zeros(size(dh0))
            else
                # scale height distance relative to time (i.e. length-scale)
                pts = hcat(t0[valid], h0[valid]/param.smooth_h2t_length_scale)'
                kdtree = KDTree(pts)
                (idxs, _) = knn(kdtree, pts, param.smooth_n[mission])
                anom0 = map(ind -> median(dh0_anom[ind]), idxs)

                # extrema(anom0)
                # interpolate anomalies using weighted distance
                itp = ScatteredInterpolation.interpolate(Shepard(2), pts, anom0);
                pts = hcat(t0[:], h0[:]/param.smooth_h2t_length_scale)'
                anom_smooth = vec(evaluate(itp, pts))
            end

            # fill out valid range (no extraploation beyond (rrange,crange) of geotile with the model only
            dh1[mission][At(geotile),rrange,crange] = model1(hcat(t0[:], h0[:]), fit1.param) .+ dh0_median .+ anom_smooth
            nobs1[mission][At(geotile), rrange, crange] = nobs0;

            # println("granule interp: $(mission) - $(geotile)")

            # add final parameters to DataFrame
            df.nobs_final = sum(nobs0);
            df.param_m1 = fit1.param;
        end
    end

    println("#############################################################################")
    println(param)
    println("Δstd = std(original) minus std(anomaly) [larger = better]")
    for mission in keys(dh1)
        valid = .!isnan.(params[mission].bin_std)
        Δ = mean(params[mission].bin_std[valid]) .- mean(params[mission].bin_anom_std[valid])
        Δ = round(Δ, digits=3)
        println("$(mission): Δstd = $(Δ)")
    end
    println("#############################################################################")
end

# fill missing values with dh hysometric mean for full region within valid range of data.
# extrapolate first and last values outside of valid elevation range for each year.
# 19s for all 4 missions for all glacierized geotiles
for mission in keys(dh1)

    # <><><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><><>
    # mission = "icesat"
    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
    for rgi in reg
        # <><><><><><><><><><><><><><><><><> FOR TESTING <><><><><><><><><><><><><><><><><>
        # rgi = first(reg)
        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
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
        valid = valid[:,crange, zrange]

        area = reduce(hcat, geotiles[rgi_ind, :glacier_area_km2])[zrange,:]
        area_total = sum(area, dims=1)
        area = repeat(area, 1, 1, length(crange))
        area = permutedims(area,(2,3,1))

        dh_specific = copy(dh0)
        dh_specific[.!valid] .= 0;
        area0 = copy(area)
        area0[.!valid] .= 0;

        # this is area0 averaged elevation change (NOTE: not a rate)
        dh_specific = dropdims(sum(dh_specific .* area0, dims=3) ./ sum(area0, dims=3), dims = 3)
        baddata = isnan.(dh_specific)

        # check how much of total regional is sampled
        showplots && heatmap(collect(dims(dh0, :date)), collect(dims(dh0, :geotile)), dh_specific)

        if all(baddata)
            println("all bad data: $(rgi)")
            continue
        end

        # give modeled values a wight of 1 by adding one to nobs
        h0 = val(dims(dh0, :height))'
        h0 = repeat(h0, length(dims(dh0, :geotile)), 1)

        # loop through each date and fill empty geotiles with second order polynomia fit to all data
        for i in eachindex(dims(dh0, :date))
            # ----------------------------------- FOR TESTING ------------------------------
            # i = 1
            #println("$(rgi): i = $(i)")
            # ------------------------------------------------------------------------------

            # give a wighting of 1 to model
            dh00 = @view dh0[:,i,:]
            valid00 = .!isnan.(dh00)

            if sum(area_total[vec(any(valid00, dims=2))]) ./ sum(area_total) .< param.min_region_coverage
                dh0[:,i,:] .= NaN;
            else
                if (sum(valid00) > (length(p2)+1))  &&  !all(valid00)
                    dh_med = median(dh00[valid00])

                    # fit 2nd order plynomial to all data within region
                    fit2 = curve_fit(model2, h0[valid00], dh00[valid00] .- dh_med, p2)
                    
                    #plot(h0[valid00], dh_fill00[valid00]; seriestype=:scatter)

                    # replace missing data
                    dh00[.!valid00] = model2(h0[.!valid00], fit2.param) .+ dh_med
                end
            end
        end
        dh1[mission][At(geotile_ids), crange, zrange] = dh0
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

# calculate volume change for each geotile and regional offsets for hugonnet dataset
# 3s for all 4 missions for all glacierized geotiles
begin
    dv = Dict();
    (ngeotile, ndate, nheight) = size(dh1[first(keys(dh1))]);
    for mission in keys(dh1)
        push!(dv, String(mission) => DimArray(fill(NaN, ngeotile, ndate), (dims(dh1[mission], :geotile), dims(dh1[mission], :date))))
    end

    glacier_area = reduce(hcat, geotiles.glacier_area_km2[:])';
    glacier_area = permutedims(repeat(glacier_area,1,1,ndate), (1,3,2));

    for mission in keys(dh1)
        dv[mission] = dropdims(sum(dh1[mission] .* glacier_area, dims=3), dims=3)/1000
    end

    # align hugonnet to icesat
    if regional_offsets
        region_offsets = DataFrame(rgi = reg, icesat_minus_hugonnet = zeros(size(reg)))
        for (i, rgi) in enumerate(reg)
            rgi_ind = geotiles[:, rgi] .> 0
            geotile_ids = geotiles[rgi_ind, :].id

            dv1 = sum(dv["icesat"][At(geotile_ids),:,1], dims = 1)
            valid1 = .!isnan.(dv1)

            dv2 = sum(dv["hugonnet"][At(geotile_ids),:,1], dims = 1)
            valid2 = .!isnan.(dv2)

            overlap = valid1 .& valid2;

            region_offsets[i, :icesat_minus_hugonnet] = median(dv1[overlap]) - median(dv2[overlap])
        end
    end
end


begin
    fig_folder = joinpath(binned_folder, "figures")
    if !isdir(fig_folder)
        mkpath(fig_folder)
    end

    out_id = replace(param.runid, "dh" => "dv")

    clr = palette(:Set1_4, length(keys(dh1)));
    for rgi in reg
        # -------------------------------------- FOR TESTING -------------------------------
        # rgi = first(reg)
        # ----------------------------------------------------------------------------------

        rgi_ind = geotiles[:, rgi] .> 0
        geotile_ids = geotiles[rgi_ind, :].id
        p = plot(ylabel="volume change [km³]", xlims=(DateTime(2000), DateTime(2024)), title="$rgi")
    
        for (i, mission) in enumerate(keys(dh1))
            dv0 = sum(dv[mission][At(geotile_ids),:,1], dims=1)

            if regional_offsets && mission == "hugonnet"
                dv0 .+= region_offsets[region_offsets.rgi .== rgi, :icesat_minus_hugonnet]
            end
            plot!(val(dims(dv0, :date)), dv0.data', label="$(mission)", color = clr[i])
        end

        fname = joinpath(fig_folder, "$(rgi)_$(out_id)")
        savefig(p, fname)
        display(p)
    end
end

using CarioMakie

f = Figure(backgroundcolor=RGBf(0.98, 0.98, 0.98),
    size=(1000, 700))


Label(ga[1, 1:2, Top()], "Stimulus ratings", valign = :bottom,
    font = :bold,
    padding = (0, 0, 5, 0))

