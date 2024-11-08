using Altim
using FileIO
using DimensionalData
using Plots
using Statistics
using LsqFit
#function geotile_align_replace(;
    mission_ref1="icesat2"
    mission_ref2="icesat"
    refperiod=(2018, 2024)
    min_trend_count=5
    remove_land_surface_trend=Altim.mission_land_trend()
    project_id=:v01
    surface_masks=[:glacier, :glacier_rgi7, :glacier_b1km, :glacier_b10km]
    binned_folders=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg")
    dem_ids=[:best, :cop30_v2]
    binning_methods=["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"]
    curvature_corrects=[true, false]
    paramater_sets=[1, 2]
    amplitude_corrects=[true, false]
    geotile_width = 2
    regions2replace_with_model = ["rgi19"]
    mission_replace_with_model = "hugonnet"
    showplots=false
    force_remake=false
#)

    # compute regional volume change, firn correction and mass change
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


    geotile = "lat[+28+30]lon[+082+084]"

    #for param in params
        param = first(params)

        # paths to files
        binned_filled_file, figure_suffix = Altim.binned_filled_filepath(; param...)
        binned_aligned_file = replace(binned_filled_file, ".jld2" => "_aligned.jld2")

        #if isfile(binned_filled_file) && ((!isfile(binned_aligned_file)) || force_remake)
            t1 = time()

            dh = FileIO.load(binned_filled_file, "dh_hyps")
            nobs_hyps = FileIO.load(binned_filled_file, "nobs_hyps")
            model_param = FileIO.load(binned_filled_file, "model_param")

            heatmap(dh["hugonnet"][At(geotile),:,:])
            heatmap!(dh["icesat2"][At(geotile),:,:])
            heatmap!(dh["icesat"][At(geotile),:,:])
            heatmap!(dh["gedi"][At(geotile),:,:])
            

        #function geotile_align!(

            mission_ref1="icesat2"
            mission_ref2="icesat"
            refperiod=(2019.5, 2023.5)
            min_trend_count=5
            remove_land_surface_trend=Altim.mission_land_trend()
            showplots=false
        #)

            missions = collect(keys(dh))
            ddata = dims(dh[missions[1]], :date)
            dheight = dims(dh[missions[1]], :height)
            dgeotile = dims(dh[missions[1]], :geotile)

            decyear = Altim.decimalyear.(ddata)
            Δdecyear = decyear .- mean(decyear)
            index_refperiod = (decyear .>= min(refperiod...)) .& (decyear .< max(refperiod...))

            offset = Dict()
            trend = Dict()
            for mission in missions
                offset[mission] = zeros(dgeotile, dheight)
                trend[mission] = zeros(dgeotile, dheight)
            end

            if showplots
                geotile0 = "lat[+28+30]lon[+082+084]"
                height0 = 4050
                p = Plots.plot(dh["hugonnet"][At(geotile0), :, At(height0)]);
                Plots.plot!(dh["icesat"][At(geotile0), :, At(height0)])
                Plots.plot!(dh["icesat2"][At(geotile0), :, At(height0)])
            end


        



            #Threads.@threads for geotile in dgeotile
                #geotile = dgeotile[700]

                (ridx, cidx) = Altim.validrange(.!isnan.(dh[mission_ref1][At(geotile), index_refperiod, :]))
                offset_ref = median(dh[mission_ref1][At(geotile), index_refperiod, cidx], dims = :date)
                fit = curve_fit(Altim.offset_trend, collect(dheight[cidx]), offset_ref[:], Altim.offset_trend_p)
                offset_ref = DimArray(Altim.offset_trend(dheight, fit.param), dheight)

                for height in dheight[cidx]
                    #height = 4050;

                    ref1 = dh[mission_ref1][At(geotile), :, At(height)]
                    ref2 = dh[mission_ref2][At(geotile), :, At(height)]
                    valid_ref1 = .!isnan.(ref1)
                    valid_ref2 = .!isnan.(ref2)

                    ref = copy(ref2)
                    ref[valid_ref1] = ref1[valid_ref1]
                    valid_ref = .!isnan.(ref)

                    # set reference time period to zero (anomalies are calculated)
                    ref .-= offset_ref[At(height)]
                    dh[mission_ref1][At(geotile), :, At(height)] .-= offset_ref[At(height)]
                    dh[mission_ref2][At(geotile), :, At(height)] .-= offset_ref[At(height)]

                    # align overlap in mean and trend
                    for mission in setdiff(missions, hcat(mission_ref1, mission_ref2))
                        #mission = "hugonnet"

                        delta = dh[mission][At(geotile), :, At(height)] .- ref
                        valid_delta = .!isnan.(delta)

                        if !isnothing(remove_land_surface_trend)
                            trend[mission][At(geotile), At(height)] += remove_land_surface_trend[At(mission)]
                            delta .-= (trend[mission][At(geotile), At(height)] .* Δdecyear)
                        end

                        if .!any(valid_delta)
                            continue
                        end

                        offset[mission][At(geotile), At(height)] += median(delta[valid_delta])

                        if (sum(valid_ref1 .& valid_delta) >= min_trend_count) && (sum(valid_ref2 .& valid_delta) >= min_trend_count)
                            # correct offset and slope

                            ## we could weight by observation here, but for now just do simply fit ##
                            fit = curve_fit(Altim.offset_trend, Δdecyear[valid_delta], (delta[valid_delta] .- offset[mission][At(geotile), At(height)]), Altim.offset_trend_p)

                            #offset_trend = delta[valid_delta] \ hcat(ones(sum(valid_delta)))
                            offset[mission][At(geotile), At(height)] += fit.param[1]
                            trend[mission][At(geotile), At(height)] += fit.param[2]
                        end

                        dh[mission][At(geotile), :, At(height)] .-= offset[mission][At(geotile), At(height)]
                        dh[mission][At(geotile), :, At(height)] -= (trend[mission][At(geotile), At(height)] .* Δdecyear)
                    end
                end
            end

            if showplots
                p = Plots.plot(dh["hugonnet"][At(geotile0), :, At(height0)]);
                Plots.plot!(dh["icesat"][At(geotile0), :, At(height0)])
                Plots.plot!(dh["icesat2"][At(geotile0), :, At(height0)])
                display(p)
            end

            return dh
        end






            dh = Altim.geotile_align!(
                dh;
                mission_ref1,
                mission_ref2,
                min_trend_count,
                remove_land_surface_trend=Altim.mission_land_trend(),
                showplots,
            )

            heatmap(dh["hugonnet"][At(geotile),:,:])
            heatmap!(dh["icesat2"][At(geotile),:,:])
            heatmap!(dh["icesat"][At(geotile),:,:])
            heatmap!(dh["gedi"][At(geotile),:,:])

            # filter geotiles to those that need some replacing with model 
            geotiles = Altim.geotiles_w_mask(geotile_width)
            keep = falses(nrow(geotiles))

            for rgi in regions2replace_with_model
                keep = keep .| (geotiles[:,rgi].> 0)
            end
            geotiles2replace = geotiles[keep, :]

            dh = replace_with_model!(dh, geotiles2replace.id; mission_replace=mission_replace_with_model, mission_ref1, mission_ref2)
           
            save(binned_aligned_file, Dict("dh_hyps" =>
                    dh, "nobs_hyps" => nobs_hyps, "model_param" => model_param, "offset" => offset, "trend" => trend))
            
            println("$binned_aligned_file aligned and replaced: $(round(Int,time() -t1))s")
        end
    end
end