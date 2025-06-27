
import GlobalGlacierAnalysis as GGA
using FileIO
using DimensionalData
using DataFrames
using Statistics

fill_param = 2
amplitude_correct = true
plots_show = false
update_geotile = true
mission_reference_for_amplitude_normalization = "icesat2"

binned_file = "/mnt/bylot-r3/data/binned/2deg/glacier_dh_cop30_v2_nmad3_v01.jld2"
binned_filled_file, figure_suffix = GGA.binned2filled_filepath(; binned_file, amplitude_correct, fill_param)
param = GGA.binned_filled_fileparts(binned_filled_file)
geotiles = GGA.geotiles_mask_hyps(param.surface_mask, 2)
single_geotile_test = GGA.geotiles_golden_test[1]

 # check if reference mission is present, if not, add
if update_geotile && (!any(keys(dh1) .== mission_reference_for_amplitude_normalization)) && isfile(binned_filled_file)
    params_fill_reference = FileIO.load(binned_filled_file, "model_param")[mission_reference_for_amplitude_normalization]
else
    params_fill_reference = nothing
end

dh1 = FileIO.load(binned_file, "dh_hyps");
nobs1 = FileIO.load(binned_file, "nobs_hyps");
param1 = param;


# align geotile dataframe with DimArrays
surface_mask_area_km2 = DimArray((reduce(hcat, geotiles[:, "$(param.surface_mask)_area_km2"])'), (Dim{:geotile}(geotiles.id), dims(dh1[first(keys(dh1))], :height)))
geotile_extent = DimArray(geotiles[:, "extent"], Dim{:geotile}(geotiles.id))

# filter geotiles to those that need some replacing with model
keep = falses(nrow(geotiles))

for rgi in regions2replace_with_model
    keep = keep .| (geotiles[:, rgi] .> 0)
end
geotiles2replace = geotiles.id[keep]

if .!isnothing(single_geotile_test)
    for k in keys(dh1)
        ind = findfirst(dims(dh1[k], :geotile) .== single_geotile_test)
        dh1[k] = dh1[k][ind:ind, :, :]
        nobs1[k] = nobs1[k][ind:ind, :, :]
    end

    # remove empty missions
    for k in keys(dh1)
        if all(isnan.(dh1[k]))
            delete!(dh1, k)
            delete!(nobs1, k)
        end
    end
end

if update_geotile && isfile(binned_filled_file)
    old_keys = setdiff(keys(dh1), update_geotile_missions)
    for k in old_keys
        delete!(dh1, k)
    end
end


# start of function
function geotile_binned_fill!(dh1, nobs1, geotile_extent, surface_mask_area_km2;
    fill_param = 2,
    amplitude_correct = true,
    params_fill_reference = nothing,
    remove_land_surface_trend = GGA.mission_land_trend(),
    missions2align2 = ["icesat2", "icesat"],
    geotiles2replace = nothing,
    mission_replace_with_model = "hugonnet",
    process_time_show = false,
    fig_folder = nothing,
    plots_show = false,
    plots_save = false,
    geotiles2plot = GGA.geotiles_golden_test[1:1],
    plot_save_format = ".png"
)

    param_filling = GGA.binned_filling_parameters[fill_param]
    process_time_show ? t6 = time() : nothing

    # create a data frame to store model parameters
    params_fill = Dict()
    for mission in keys(dh1)
        n = length(dims(dh1[mission], :geotile))

        params_fill[mission] = DataFrame(
            geotile=val(dims(dh1[mission], :geotile)), 
            nobs_raw=zeros(n), nbins_raw=zeros(n), 
            nobs_final=zeros(n), nbins_filt1=zeros(n), 
            param_m1=[fill(NaN, size(GGA.p1)) for i in 1:n], 
            h0=fill(NaN, n), 
            t0=fill(NaN, n), 
            dh0=fill(NaN, n), 
            bin_std=fill(NaN, n), 
            bin_anom_std=fill(NaN, n),
            )

        for mission_ref in missions2align2
            params_fill[mission][!, "offset"] = zeros(n)
            params_fill[mission][!, "offset_$mission_ref"] = fill(NaN, n)
            params_fill[mission][!, "offset_nmad_$mission_ref"] = fill(NaN, n)
            params_fill[mission][!, "offset_nobs_$mission_ref"] = zeros(Int64,n)
        end
    end

    process_time_show ? t7 = time() : nothing
    process_time_show && printstyled("params_fill initialized: $(round(Int16, t7 - t6))s\n", color=:green)

    # plot raw binned height anomalies
    if plots_show || plots_save
        colorbar_label = "binned height anomalies [m]"
        plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")

        plot_elevation_time_multimission_geotiles(
            dh1,
            geotiles2plot;
            area_km2=surface_mask_area_km2,
            colorrange=(-20, 20),
            colorbar_label,
            hypsometry=true,
            area_averaged=true,
            plots_show,
            plots_save,
            plot_save_path_prefix,
            plot_save_format,
        )
    end

    # interpolate height anomalies
    begin
        GGA.hyps_model_fill!(dh1, nobs1, params_fill; bincount_min=param_filling.bincount_min,
            model1_nmad_max=param_filling.model1_nmad_max, smooth_n=param_filling.smooth_n,
            smooth_h2t_length_scale=param_filling.smooth_h2t_length_scale, show_times = false)

        if plots_show || plots_save
            colorbar_label = "interpolated height anomalies [m]"
            plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")

            plot_elevation_time_multimission_geotiles(
                dh1,
                geotiles2plot;
                area_km2=surface_mask_area_km2,
                colorrange=(-20, 20),
                colorbar_label,
                hypsometry=true,
                area_averaged=true,
                plots_show,
                plots_save,
                plot_save_path_prefix,
                plot_save_format,
            )
        end


        process_time_show ? t8 = time() : nothing
        process_time_show && printstyled("model_fill: $(round(Int16, t8 - t7))s\n", color=:green)
    end

    # apply seasonal amplitude normalization
    if amplitude_correct

        # check if reference mission is present, if not, add
        if isnothing(params_fill_reference)
            params_fill_reference = params_fill[mission_reference_for_amplitude_normalization]
        end

        for mission in setdiff(keys(dh1), [mission_reference_for_amplitude_normalization])
            GGA.hyps_amplitude_normalize!(dh1[mission], params_fill[mission], params_fill_reference)
        end


        if plots_show || plots_save
            colorbar_label = "normalized height anomalies [m]"
            plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")
            plot_elevation_time_multimission_geotiles(
                dh1,
                geotiles2plot;
                area_km2=surface_mask_area_km2,
                colorrange=(-20, 20),
                colorbar_label,
                hypsometry=true,
                area_averaged=true,
                plots_show,
                plots_save,
                plot_save_path_prefix,
                plot_save_format,
            )

        end
    end

    process_time_show ? t9 = time() : nothing
    process_time_show && printstyled("amplitude_normalize: $(round(Int16, t9 - t8))s\n", color=:green)

    # fill geotiles
    begin
        # hyps_fill_empty! can add mission data to geotiles that would 
        # otherwise be empty. an example of this is lat[+60+62]lon[-142-140] 
        # which has not GEDI data but GEDI data is added after hyps_fill_empty! 
        # becuase at least on of its 5 closest neighbors have GEDI data
        dh1 = GGA.hyps_fill_empty!(dh1, params_fill, geotile_extent, surface_mask_area_km2)

        process_time_show ? t10 = time() : nothing
        process_time_show && printstyled("hyps_fill_empty: $(round(Int16, t10 - t9))s\n", color=:green)

        dh1 = GGA.hyps_fill_updown!(dh1, surface_mask_area_km2)

        process_time_show ? t11 = time() : nothing
        process_time_show && printstyled("hyps_fill_updown: $(round(Int16, t11 - t10))s\n", color=:green)

        if plots_show || plots_save
            colorbar_label = "extrapolated height anomalies [m]"
            plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")
            plot_elevation_time_multimission_geotiles(
                dh1,
                geotiles2plot;
                area_km2=surface_mask_area_km2,
                colorrange=(-20, 20),
                colorbar_label,
                hypsometry=true,
                area_averaged=true,
                plots_show,
                plots_save,
                plot_save_path_prefix,
                plot_save_format,
            )
        end
    end

    # correct for any erronious trends found over land
    begin
        dgeotile = dims(dh1[first(keys(dh1))], :geotile)
        dheight = dims(dh1[first(keys(dh1))], :height)
        
        for mission in keys(dh1)
            if !isnothing(remove_land_surface_trend) && (remove_land_surface_trend[At(mission)] != 0)
                ddate = dims(dh1[mission], :date)
                decyear = GGA.decimalyear.(ddate)

                # center date arround mission center date
                _, date_range, _ = GGA.validrange(.!isnan.(dh1[mission]))
                mid_date = mean(decyear[date_range])
                delta = (decyear .- mid_date) .* remove_land_surface_trend[At(mission)]

                for geotile in dgeotile
                    for height in dheight
                        dh1[mission][geotile=At(geotile), height=At(height)] .-= delta
                    end
                end
            end
        end

        if plots_show || plots_save
            colorbar_label = "land surface trend corrected height anomalies [m]"
            plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")

            plot_elevation_time_multimission_geotiles(
                dh1,
                geotiles2plot;
                area_km2=surface_mask_area_km2,
                colorrange=(-20, 20),
                colorbar_label,
                hypsometry=true,
                area_averaged=true,
                plots_show,
                plots_save,
                plot_save_path_prefix,
                plot_save_format,
            )
        end
    end

    # align height nomalies to reference missions
    begin
        
        dh1, params_fill = GGA.hyps_align_dh!(dh1, nobs1, params_fill, area_km2; missions2align2)

        if plots_show || plots_save
            colorbar_label = "adjusted height anomalies [m]"
            plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")

            plot_elevation_time_multimission_geotiles(
                dh1,
                geotiles2plot;
                area_km2=surface_mask_area_km2,
                colorrange=(-20, 20),
                colorbar_label,
                hypsometry=true,
                area_averaged=true,
                plots_show,
                plots_save,
                plot_save_path_prefix,
                plot_save_format,
            )
        end
    end

    # fill geotiles with model
    begin
        dh1, nobs1 = GGA.replace_with_model!(dh1, nobs1, geotiles2replace; mission_replace=mission_replace_with_model, missions2align2)

        if plots_show || plots_save
            colorbar_label = "model-filled height anomalies [m]"
            plot_save_path_prefix = joinpath(fig_folder, "$(figure_suffix)_$(replace(colorbar_label, " " => "_"))")

            plot_elevation_time_multimission_geotiles(
                dh1,
                geotiles2plot;
                area_km2=surface_mask_area_km2,
                colorrange=(-20, 20),
                colorbar_label,
                hypsometry=true,
                area_averaged=true,
                plots_show,
                plots_save,
                plot_save_path_prefix,
                plot_save_format,
            )
        end
    end

    return (dh1, nobs1, params_fill)
end


function geotile_binned_fill!(dh1, nobs1, geotile_extent, surface_mask_area_km2;
    fill_param = 2,
    amplitude_correct = true,
    mission_reference_for_amplitude_normalization = "icesat2",
    params_fill_reference = nothing,
    remove_land_surface_trend = GGA.mission_land_trend(),
    missions2align2 = ["icesat2", "icesat"],
    geotiles2replace = nothing,
    mission_replace_with_model = "hugonnet",
    process_time_show = false,
    fig_folder = nothing,
    plots_show = false,
    plots_save = false,
    geotiles2plot = GGA.geotiles_golden_test[1:1],
    plot_save_format = ".png"
)


# update dh1 and nobs1 with new missions
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
if isnothing(single_geotile_test)
    save(binned_filled_file, Dict("dh_hyps" => dh1, "nobs_hyps" => nobs1, "model_param" => params_fill))

    process_time_show ? t12 = time() : nothing
    process_time_show && printstyled("saving binned_filled_file: $(round(Int16, t12 - t11))s\n", color=:green)
end