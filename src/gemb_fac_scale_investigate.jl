
begin
    using Altim
    using FileIO
    using DimensionalData
    using Statistics
    using LsqFit
    using Plots
    using Arrow
    using DataFrames
    using Extents

    project_id = :v01;
    geotile_width = 2;
    binned_folder = analysis_paths(; geotile_width).binned
    binning_method_gemb = "mean"; 
    runid_gemb = "glacier_gemb_$(binning_method_gemb)_$(project_id)"
    gemb_file = joinpath(binned_folder, "$(runid_gemb).jld2");
    gemb = load(gemb_file);

    mask_project  = :glacier
    gt_file = joinpath(binned_folder, "geotile_$(mask_project )_hyps.arrow")
    geotiles = DataFrame(Arrow.Table(gt_file))

    valid_geotiles = vec(any(.!isnan.(first(gemb)[2]), dims=(2, 3)))
    for k in keys(gemb)
        if k == "t_air_hyps"
            gemb[k] .-= 273.15
        elseif (k == "nobs_hyps") || (k == "fac_hyps") || (k == "fac_to_depth_hyps")
            
        else
            gemb[k] ./= 1000
            gemb[k] = cumsum(gemb[k], dims=2)
        end

        gemb[k] =  gemb[k][valid_geotiles,:,:]
    end

    # add smb and fac to get height change
    push!(gemb, "height_hyps" => gemb["smb_hyps"] .+ gemb["fac_hyps"])

    gt = collect(dims(gemb[first(keys(gemb))], :geotile))
    gt_ind = [findfirst(geotiles.id .== g) for g in gt]
    geotiles = geotiles[gt_ind,:]

    model3::Function = model3(t, p) = p[1] .+ p[2] .* t .+ p[3] .* t .^ 2 .+ p[4].*cos.(2 * pi * t) .+ p[5].*sin.(2 * pi * t)
    p3 = zeros(5)
    t = Altim.decimalyear.(dims(first(gemb)[2], :date))

    # center time with integer value
    t = t .- round(mean(t))

    #
    dparameter = Dim{:parameter}(1:length(p3))
        ngeotile = length(dims(first(gemb)[2], :geotile))
    nheight = length(dims(first(gemb)[2], :height))
    nparameter = length(dparameter)

    param = Dict()
    for k in keys(gemb)
        param[k] = DimArray(fill(NaN, ngeotile, nheight, nparameter), (dims(first(gemb)[2], :geotile), dims(first(gemb)[2], :height), dparameter));
    end;

    # 9s
    Threads.@threads for geotile in dims(first(gemb)[2], :geotile)
        #geotile = dims(fac, :geotile)[565]
        for h in dims(first(gemb)[2], :height)
        #h = dims(fac, :height)[15]
            valid = .!isnan.(gemb[first(keys(gemb))][At(geotile), :, At(h)])
            if any(valid)
                for k in keys(gemb)
                    foo = gemb[k][At(geotile), :, At(h)][valid]
                    param[k][At(geotile), At(h), :] = curve_fit(model3, t[valid], foo, p3).param
                end
            end
        end
    end

    phase(model_hyps) = atan.(model_hyps[:, :, 4], model_hyps[:, :, 5]) / (2 * pi)
    amplitude(model_hyps) = sqrt.(model_hyps[:, :, 4] .^ 2 .+ model_hyps[:, :, 5] .^ 2)

    amp_phase = Dict()
    dvar = Dim{:parameter}(["amplitude", "phase"])
    for k in keys(gemb)
        amp_phase[k] = DimArray(fill(NaN, ngeotile, nheight, length(dvar)), (dims(first(gemb)[2], :geotile), dims(first(gemb)[2], :height), dvar))
    end;

    for k in keys(gemb)
        amp_phase[k][:, :, At("amplitude")] = amplitude(param[k])
        amp_phase[k][:, :, At("phase")] = phase(param[k])
    end
end

# Question: How much error does scaling FAC introduce
valid = .!isnan.(amp_phase[first(keys(amp_phase))])
# refreeze and fac amplitude are stongly related
x = param["refreeze_hyps"][:, :, 1][valid[:, :, 1]]
y = amp_phase["fac_hyps"][:, :, At("amplitude")][valid[:, :, 1]]
plot(x, y; legend=false, seriestype=:scatter)

reg_ind = geotiles.rgi1 .> 0;
x = param["refreeze_hyps"][reg_ind, :, 1][valid[reg_ind, :, 1]]
#x = amp_phase["acc_hyps"][reg_ind, :, At("amplitude")][valid[reg_ind, :, 1]]
y = amp_phase["fac_hyps"][reg_ind, :, At("amplitude")][valid[reg_ind, :, 1]]
plot(x, y; legend=false, seriestype=:scatter, xlims = (0,5))

model1(x,p) = x[:,1] .* p[1] .+ x[:,2] .* p[2] .+ x[:,3] .* p[3] .+ x[:,4] .* p[4]

x = hcat(
    param["t_air_hyps"][reg_ind, :, 1][valid[reg_ind, :, 1]],
    amp_phase["t_air_hyps"][reg_ind, :, At("amplitude")][valid[reg_ind, :, 1]],
    param["acc_hyps"][reg_ind, :, 1][valid[reg_ind, :, 1]],
    amp_phase["acc_hyps"][reg_ind, :, At("amplitude")][valid[reg_ind, :, 1]]
)

y = amp_phase["height_hyps"][reg_ind, :, 1][valid[reg_ind, :, 1]]

p1 = zeros(4)
fit1 = curve_fit(model1, x, y, p1)
plot(model1(x, fit1.param), y; seriestype=:scatter)