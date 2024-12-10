# This script processes GEMB (Glacier Energy Mass Balance) model data into geotiles and regions.
# Key processing steps:
#
# 1. Initialization (~1 min):
#    - Sets up project parameters (dates, heights, precipitation scales)
#    - Loads required packages and geotile definitions
#    - Configures input/output paths
#
# 2. GEMB Data Merging (~3 min):
#    - Combines multiple GEMB .mat files into single dataset
#    - Handles precipitation scaling and elevation delta variations
#    - Processes coordinates and saves merged data
#
# 3. Geotile Processing (~3.5 min):
#    - Bins GEMB data into time-height geotiles
#    - Buffers geotile extents to ensure coverage
#    - Tracks observation counts
#    - Handles multiple precipitation scales
#
# 4. Gap Filling and Height Classes (~9 hrs):
#    - Fills missing values across elevations using LOESS smoothing
#    - Creates new elevation classes through height shifts
#    - Applies physical constraints on variables
#    - Calculates volume changes using glacier areas
#    - Saves filled data for each geotile
#
# 5. Regional Aggregation (~11 sec):
#    - Aggregates geotile data into RGI regions
#    - Calculates regional glacier areas
#    - Sums volume changes by region
#    - Saves regional metrics
#
# Total runtime is approximately 9.5 hours, with most time spent on gap filling
# and creating height classes in step 4.


@time begin
    # set force_remake == true to redo all steps from scratch
    force_remake = true;

    # This section initializes the GEMB (Glacier Energy Mass Balance) model:
    #
    # Key steps:
    # 1. Loads required packages for data processing, statistics, and file I/O
    # 2. Sets model parameters:
    #    - Project ID and version
    #    - Geotile width and buffer distances
    #    - Date ranges from 1940 onwards
    #    - Height bins and elevation deltas
    #    - Precipitation scaling factors
    # 3. Configures input/output file paths
    # 4. Defines data binning functions (mean, median, etc)
    # 5. Loads and processes geotile definitions:
    #    - Filters to glacier-containing tiles
    #    - Makes RGI regions mutually exclusive
    #    - Standardizes area column names
    # 6. Sets minimum required GEMB coverage threshold
    begin
        using Altim 
        using Arrow
        using DataFrames
        using Extents
        using Rasters
        using GeoInterface
        using Shapefile
        using DimensionalData
        using Statistics
        using Dates
        using Plots
        using JLD2
        using FileIO
        using MAT
        using GeoTiles
        using DimensionalData
        using DataInterpolations
        using ProgressMeter
        using Loess

        # run parameters
        project_id = :v01;
        geotile_width = 2;
        binning_method = "mean";
        surface_mask = :glacier
        geotile_buffer = 50000 # distance in meters outside of geotiles to look for data

        # define date and hight binning ranges 
        date_range, date_center = Altim.project_date_bins()

        # expand daterange to 1940 by make sure to match exisiting project ranges 
        foo = collect(date_range[1]:-Day(30):DateTime(1940))

        date_range = foo[end]:Day(30):date_range[end]
        date_center = date_range[1:end-1] .+ Day(Day(30) / 2)

        height_range, height_center = Altim.project_height_bins()

        #Δheight simulates changing model elevation to increase / decrease melt, this is done in the regional Δvolume calculation
        height_bin_interval = height_center[2] - height_center[1]
        dΔheight = Dim{:Δheight}(-3000:height_bin_interval:3000) 

        # exclude derived variables of smb and runoff (these are calculated later to ensure mass conservation after interpolation)
        vars = ["latitude", "longitude", "date", "fac", "acc", "refreeze", "height", "melt", "rain", "ec"]

        if  false
            gemb_folder = ["/home/schlegel/Share/GEMBv1/"];
            file_uniqueid = "rv1_0_19500101_20231231"
            elevation_delta = [0]
            precipitation_scale = [1]
            filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/$file_uniqueid.jld2"

        elseif false
            gemb_folder = "/home/schlegel/Share/GEMBv1/Alaska_sample/v1/"
            file_uniqueid = "1979to2023_820_40_racmo_grid_lwt"
            elevation_delta = [-1000, -750 ,-500, -250, 0, 250, 500, 750, 1000]
            precipitation_scale = [0.5, 1, 1.5, 2, 5, 10]
            filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/$file_uniqueid.jld2"
        else
            gemb_folder = ["/home/schlegel/Share/GEMBv1/NH_sample/", "/home/schlegel/Share/GEMBv1/SH_sample/"]
            file_uniqueid = nothing
            elevation_delta = [-200, 0, 200, 500, 1000]
            precipitation_scale = [.75 1 1.25 1.5 2]
            filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0.jld2"

        end

        filename_gemb_geotile = replace(filename_gemb_combined, ".jld2" => "_geotile.jld2")
        filename_gemb_geotile_filled_dv = replace(filename_gemb_geotile, ".jld2" => "_filled_dv.jld2")
        filename_gemb_geotile_filled_dv_reg = replace(filename_gemb_geotile, ".jld2" => "_filled_dv_reg.jld2")

        # funtion used for binning data
        if binning_method == "meanmadnorm3"
            binfun::Function = binfun(x) = mean(x[Altim.madnorm(x).<3])
        elseif binning_method == "median"
            binfun::Function = binfun(x) = median(x)
        elseif binning_method == "mean"
            binfun::Function = binfun(x) = mean(x)
        else
            error("unrecognized binning method")
        end

        showplots = false;

        # load geotile definitions with corresponding hypsometry
        geotiles = Altim.geotiles_mask_hyps(surface_mask, geotile_width)

        # filter geotiles
        geotiles = geotiles[geotiles.glacier_frac.>0.0, :]

        # make geotile rgi regions mutually exexclusive 
        geotiles, reg = Altim.geotiles_mutually_exclusive_rgi!(geotiles)

        if !any(Base.contains.("area_km2", names(geotiles)))
            original_area_name = string(surface_mask)*"_area_km2"
            generic_area_name = "area_km2"
            rename!(geotiles, original_area_name => generic_area_name)
        end

        # min gemb coverage =  0.5
        min_gemb_coverage =  0.75
    end;


    # Merges multiple GEMB data files into a single combined file.
    #
    # This section of code:
    # 1. Reads all .mat GEMB files from specified folders
    # 2. Extracts precipitation scale and elevation delta indices from filenames
    # 3. Combines data from all files, handling both single and multiple class cases
    # 4. Processes coordinates to ensure longitudes are in -180 to 180 range
    # 5. Merges all variables into a single dictionary
    # 6. Saves the combined data to a JLD2 file
    #
    # The combined file contains:
    # - Dates (as DateTime)
    # - Latitude/longitude coordinates
    # - Heights
    # - Precipitation scale values
    # - Elevation delta values
    # - All other GEMB variables specified in vars
    #
    # Takes ~3 minutes to run for all data.
    if .!(isfile(filename_gemb_combined)) || force_remake

        println("merging all gemb data into a single file, this will take ~3 min")

        gemb_files = allfiles.(gemb_folder; subfolders=false, fn_endswith=".mat", fn_contains=file_uniqueid)
        gemb_files = vcat(gemb_files...)

        expected_num_files = length(elevation_delta) * length(precipitation_scale) * length(gemb_folder)
        #if length(gemb_files) !== expected_num_files
        #    error("the number of files does not match the number of classes")
        #end

        datebin_edges = Altim.decimalyear.(date_range)
        # this takes 4 m
        gemb = []
        
        Threads.@threads for gemb_file in gemb_files
        #gemb_file = first(gemb_files)

            gemb0 = Altim.gemb_read2(gemb_file; vars, datebin_edges)

            # no classes
            if (length(elevation_delta) == 1) .&  (length(precipitation_scale) == 1)
                precip_scale_ind = 1
                elev_delta_ind = 1
            else
                ind1 = findlast("_p", gemb_file)
                ind2 = findlast("_t", gemb_file)
                ind3 = findlast(".mat", gemb_file)
            
                precip_scale_ind = parse(Int8, gemb_file[ind1[end]+1 : ind2[1]-1])
                elev_delta_ind = parse(Int8, gemb_file[ind2[end]+1 : ind3[1]-1])
            end

            gemb0["elevation_delta"] = elevation_delta[elev_delta_ind]
            gemb0["precipitation_scale"] = precipitation_scale[precip_scale_ind]

            push!(gemb, gemb0)
        end

        if (length(elevation_delta) == 1) .&  (length(precipitation_scale) == 1)
            gemb0 = Dict()
            for k in keys(gemb[1])
                if k !== "date"
                    gemb0[k] = reduce(vcat,[g[k] for g in gemb])
                else
                    gemb0[k] = gemb[1][k] 
                end
            end
            gemb = [gemb0]
        end

        # merge all simulations into a single file
        dates = gemb[1]["date"]
        latitude = reduce(vcat, getindex.(gemb, Ref("latitude")))
        longitude = reduce(vcat, getindex.(gemb, Ref("longitude")))
        height0 = reduce(vcat, getindex.(gemb, Ref("height")))

        # make sure longitude is in -180 to 180 range
        index = longitude .< -180
        longitude[index] = longitude[index] .+ 360
        index = longitude .> 180
        longitude[index] = longitude[index] .- 360

        precipitation_scale0 = Float64[]
        elevation_delta0 = Float64[]
        for g in gemb
            precipitation_scale0 = vcat(precipitation_scale0, fill(g["precipitation_scale"], length(g["latitude"])))
            elevation_delta0 = vcat(elevation_delta0, fill(g["elevation_delta"], length(g["latitude"])))
        end

        gemb0 = Dict(
            "date" => Altim.decimalyear2datetime.(vec(dates)),
            "latitude" => vec(latitude), 
            "longitude" => vec(longitude),
            "height" => vec(height0), 
            "precipitation_scale" =>  vec(precipitation_scale0), 
            "elevation_delta" =>  vec(elevation_delta0)
            )

        foo = Dict()
        for k in setdiff(vars, ["date", "latitude", "longitude", "height"])
        #k = first(setdiff(vars, ["date", "latitude", "longitude", "height"]))
            foo[k] = reduce(vcat, getindex.(gemb, Ref(k)))
        end
        gemb = merge(foo, gemb0)

        # save output
        save(filename_gemb_combined, gemb);
    end


    # This section combines GEMB model runs into time-height geotiles that match the dh processing:
    # - Loads combined GEMB data from file
    # - Processes elevation_delta runs to fill time-height space
    # - Handles precipitation_scale runs as unique model runs
    # - Creates dimensions for height, date, geotile ID and precipitation scale
    # - For each geotile:
    #   - Buffers the geotile extent until sufficient data coverage
    #   - Bins data into height classes for each precipitation scale
    #   - Calculates means when multiple points exist
    #   - Tracks number of observations
    # - Saves processed data to JLD2 file
    # Takes ~3.5 minutes to run
    if .!(isfile(filename_gemb_geotile)) || force_remake
        println("putting gemb data into geotiles, this will take ~3.5 min")
        # combine into geotiles
        gemb = load(filename_gemb_combined);

        # add elevation_delta to height to get height_effective
        gemb["height_effective"] = gemb["height"] .+ gemb["elevation_delta"]
        
        unique_precipitation_scale = unique(gemb["precipitation_scale"])

        vars = setdiff(collect(keys(gemb)), ("precipitation_scale", "elevation_delta", "latitude", "longitude", "height", "date", "height_effective"))

        # create or extract dimensions
        dheight = Dim{:height}(height_center)
        ddate = Dim{:date}(date_center)
        dgeotile = Dim{:geotile}(geotiles.id)
        dpscale = Dim{:pscale}(unique_precipitation_scale)

        gemb0 = Dict();
        for k in vars
            gemb0[k] = fill(NaN, (dgeotile, ddate, dheight, dpscale))
        end
        gemb0["nobs"] = fill(0, (dgeotile, ddate, dheight, dpscale))

        Threads.@threads for geotile in eachrow(geotiles)
        #geotile = geotiles[500, :]

            # itteratively buffer by geotile_buffer meters untill 
            ext = geotile.extent
            latitude_distance, longitude_distance = Altim.meters2lonlat_distance(geotile_buffer, mean(ext.Y))
            has_data = falses(size(height_center))
            in_geotile = falses(size(gemb["longitude"]))

            total_area = sum(geotile.area_km2);
            if total_area == 0
                for k in vars
                    gemb0[k][At(geotile.id), :, :, :] .= 0.
                end

                continue
            end

            while (sum(geotile.area_km2[has_data]) / total_area) < min_gemb_coverage
                in_geotile = vec([Altim.within(ext, x, y) for (x, y) in zip(gemb["longitude"], gemb["latitude"])])
                h = gemb["height_effective"][in_geotile]

                for i in 1:(length(height_range)-1)
                    has_data[i] = any((h .>= height_range[i]) .& (h .< height_range[i+1]))
                end

                # buffer for next iteration
                ext = Extents.buffer(ext, (X=longitude_distance, Y=latitude_distance))
            end
                    
            for i in 1:(length(height_range)-1)
            #i = 50

                index_height = (gemb["height_effective"] .>= height_range[i]) .& (gemb["height_effective"] .< height_range[i+1])
                height0 = height_center[i]

                if !any(index_height)
                    continue
                end

                #loop for each precipitation scale
                for pscale in dpscale
                #pscale = dpscale[2]

                    index_precipitation = gemb["precipitation_scale"] .== pscale
                    index = index_precipitation .& in_geotile .& index_height
                    sum_index = sum(index)

                    if sum_index == 0
                        continue
                    elseif sum_index == 1
                        for k in vars
                            gemb0[k][At(geotile.id), :, At(height0), At(pscale)] = gemb[k][index,:];
                        end
                    else
                        for k in vars
                            gemb0[k][At(geotile.id), :, At(height0), At(pscale)] = mean(gemb[k][index,:], dims=1);
                        end
                    end

                    valid_date = .!isnan.(gemb0[vars[1]][At(geotile.id), :, At(height0), At(pscale)])
                    gemb0["nobs"][At(geotile.id), valid_date, At(height0), At(pscale)] .= sum(index);
                end
            end
        end
        save(filename_gemb_geotile, gemb0)
    end

    # This section fills missing GEMB variables across elevations and creates new elevation classes.
    # Key steps:
    #
    # 1. Loads GEMB geotile data and extracts dimensions
    # 2. For each variable (except nobs, smb, runoff):
    #    - Smooths and extrapolates values across heights using LOESS
    #    - Calculates rates and handles missing values
    #    - Applies physical constraints (e.g. positive rates for acc, refreeze, melt, rain)
    #    - Uses linear extrapolation to fill gaps
    #
    # 3. Creates volume change results dictionary
    #
    # 4. For each geotile:
    #    - Shifts elevations to create new height classes
    #    - Applies physical constraints
    #    - Converts rates back to cumulative variables
    #    - Saves filled data to disk
    #    - Calculates volume change time series by multiplying by glacier area
    #
    # The process takes approximately 9 hours to complete and handles:
    # - Smoothing and gap filling across elevations
    # - Creation of new elevation classes
    # - Volume change calculations
    # - Physical constraints on variables
    #
    # NOTE: Linear extrapolation over long distances is not ideal and requires
    # running GEMB for a wider range of elevation classes to improve results.
    if .!(isfile(filename_gemb_geotile_filled_dv)) || force_remake

        printstyled("!!!!! FILLING GEMB GEOTILES, THIS WILL TAKE ~9 HRS !!!!!\n", color=:red)
        
        # fill filename_gemb_geotile
        gemb = load(filename_gemb_geotile)
        dgeotile = dims(gemb[first(keys(gemb))], :geotile)
        ddate = dims(gemb[first(keys(gemb))], :date)
        dheight = dims(gemb[first(keys(gemb))], :height)
        dpscale = dims(gemb[first(keys(gemb))], :pscale)
        x = collect(dheight);

        # NOTE: smb and runoff are excluded from the filling and extrapolation as they are derivatives 
        # of melt, refreeze and acc. To maitain mass conservation "smb" and "runoff" are 
        # reonstucted from native variables after the filling and extrapolation is complete.
        vars = setdiff(collect(keys(gemb)), ["nobs", "smb", "runoff"]);
        
        # test subset of geotiles
        ### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ###
        if false
            gtids = ["lat[-14-12]lon[-076-074]", "lat[-10-08]lon[-078-076]", "lat[+28+30]lon[+082+084]", "lat[+58+60]lon[-136-134]", "lat[+68+70]lon[-070-068]", "lat[+62+64]lon[-044-042]", "lat[+32+34]lon[+076+078]", "lat[+28+30]lon[+092+094]", "lat[-46-44]lon[+168+170]", "lat[-70-68]lon[-068-066]", "lat[-52-50]lon[-074-072]"]

            test_ind = [findfirst(isequal(gtid), dgeotile.val) for gtid in gtids]
            dgeotile = dgeotile[test_ind]
        end
        ### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ###

        # This take 36 minutes
        for k in vars
        #k = "acc"
            @showprogress desc="Smoothing and extrapolating $k for all heights... "  Threads.@threads for geotile in dgeotile   
            #geotile = "lat[+34+36]lon[+076+078]"

                for pscale in dpscale
                #pscale = 1

                    # there is an issue with the interpolating the cumulitive variables... interpolation needs to be done on rates to ensure they are physical
                    var0 = gemb[k][At(geotile), :, :, At(pscale)]

                    # sanity check
                    #CairoMakie.heatmap(var0)

                    # if all values are 0 then skip
                    if all(var0 .== 0)
                        continue
                    end

                    # calculate rates
                    ind = findfirst(.!vec(all(isnan.(var0), dims=:height))) - 1
                    var0[ind, :] .= 0
                    var0_rate = var0[2:end,:] .- collect(var0[1:end-1,:])

                    # sanity check
                    #CairoMakie.heatmap(var0_rate)

                    for date in dims(var0_rate, :date)
                    #date = ddate[925]

                        #println((; k, geotile, date, pscale))
                        y = var0_rate[At(date), :]
                        
                        # sanity check
                        #p = lines(y)
                
                        valid = .!isnan.(collect(y))
                    
                        # only positive rates are valid for these variables
                        if k in ("acc", "refreeze", "melt", "rain")
                            valid = valid .& (y .>= 0)
                        end

                        if all(valid) || !any(valid)
                            continue
                        end

                        # smooth variables with height and interpolate gaps
                        validgap = Altim.validgaps(valid)

                        # NOTE: interplation increases per variable time from 1 min to 6 min 
                        if sum(valid) > 3
                            model = loess(x[valid], y[valid], span=0.3)
                            y[valid .| validgap]  = predict(model, x[valid .| validgap])
                        end

                        fill_index = isnan.(collect(y))
                        if (k == "ec") || (k == "fac")
                            valid_interp_index = .!fill_index
                        else
                            # do not include zeros in linear extrapolation as a funciton of elevation
                            valid_interp_index = .!fill_index .& (y .!= 0)
                        end

                        # linearly extrapolate the data
                        # simple linear extrapolation is crude but robust... more precipitation and melt classes are needed in future 
                        if any(fill_index) &&  sum(valid_interp_index) > 2
                            x0 = (x .- mean(x[valid_interp_index]))/1000
                            # apply a functional fit to extrapolate data [this does not work well for melt]
                            M = hcat(ones(sum(valid_interp_index)), x0[valid_interp_index]);
                            param1 = M\y[valid_interp_index]

                            y[fill_index] = max.(param1[1] .+ x0[fill_index].*param1[2], 0)
                        else
                            y[fill_index] .= 0
                        end

                        gemb[k][At(geotile), At(date), :, At(pscale)] = y;
                    end
                end
            end
        end

        # set some physcial constraints
        gemb = Altim.gemb_rate_physical_constraints!(gemb)

        # create a DD for volume change results
        gemb_dv = Dict()
        for k in keys(gemb)
            gemb_dv[k] = fill(NaN, (dgeotile, ddate, dpscale, dΔheight))
        end

        # add height classes and save to disk for each geotile [THIS TAKES 8 HRS]
        @showprogress desc="Adding height classes [expect 9 hrs]..." for geotile in dgeotile

            filename_gemb_geotile_filled = replace(filename_gemb_geotile, ".jld2" => "_$(geotile)_filled.jld2")

            gemb_gt = Dict()
            for k in vars
                gemb_gt[k] = fill(NaN, (ddate, dheight, dpscale, dΔheight))
            end

            Threads.@threads for k in vars
            #k = "runoff"
                for Δheight in dΔheight
                #Δheight = dΔheight[7]
            
                    for pscale in dpscale
                    #pscale = dpscale[5]
        
                        if k == "nobs"
                            gemb_gt[k][:,:,At(pscale),At(Δheight)] = gemb[k][At(geotile), :, :, At(pscale)]
                        else
                            var0 = copy(gemb[k][At(geotile), :, :, At(pscale)])
            
                            # shift elevation
                            height_shift = round(Int16(Δheight / height_bin_interval))

                            if (k == "ec") || (k == "fac")
                                gemb_gt[k][:,:,At(pscale),At(Δheight)] = Altim._matrix_shift_ud!(var0, height_shift; exclude_zeros_in_extrapolation = false)
                            else
                                gemb_gt[k][:,:,At(pscale),At(Δheight)] = Altim._matrix_shift_ud!(var0, height_shift; exclude_zeros_in_extrapolation = true)
                            end
                        end
                    end
                end
            end

            gemb_gt = Altim.gemb_rate_physical_constraints!(gemb_gt)
            
            # transform back into cumulative variables
            for k in setdiff(collect(keys(gemb_gt)), ["nobs"])
                
                for height in dheight
                #height = dheight[20]

                    for pscale in dpscale
                    #pscale = first(dpscale)
                        
                        for Δheight in dΔheight
                            y = gemb_gt[k][:,At(height),At(pscale),At(Δheight)]
                            valid = .!isnan.(y)
                            if any(valid)
                                gemb_gt[k][valid, At(height), At(pscale), At(Δheight)] = cumsum(y[valid], dims = :date);
                            end
                        end
                    end
                end
            end     

            # save to disk
            # NOTE: files saved to disk are only used if the hypsometric method is used for 
            # geotile2glacier downscaling].. as of now the defualt method is by AREA so the 
            # files are not used.
            save(filename_gemb_geotile_filled, gemb_gt)

            # create volume change time series for each geotile
            df_index = findfirst(isequal(geotile), geotiles.id)
            area = ones(length(ddate), 1) * hcat(geotiles[df_index, :area_km2])' ./ 1000 # divide by 1000 so that result in km3 (NOTE: vars are in mwe)
            for k in collect(keys(gemb_gt))
                for pscale in dpscale
                #pscale = first(dpscale)            
                    for Δheight in dΔheight
                        if k == "nobs"
                            gemb_dv[k][At(geotile), :, At(pscale), At(Δheight)] =  sum(gemb_gt[k][:,:,At(pscale),At(Δheight)], dims = :height);
                        else
                            gemb_dv[k][At(geotile), :, At(pscale), At(Δheight)] =  sum(gemb_gt[k][:,:,At(pscale),At(Δheight)] .* area , dims = :height);
                        end
                    end
                end
            end
        end
        save(filename_gemb_geotile_filled_dv, gemb_dv)
    end


    # Aggregates GEMB data into regional volume metrics
    #
    # This section:
    # 1. Loads the geotile-level volume change data
    # 2. Creates dimensions for RGI regions, dates, precipitation scales and height deltas
    # 3. Aligns geotile dataframe indices with DimArrays
    # 4. For each RGI region:
    #    - Calculates total glacier area in km2
    #    - Sums volume changes across all geotiles in the region
    # 5. Saves regional aggregated data to file
    #
    # Takes ~11 seconds to run
    if .!(isfile(filename_gemb_geotile_filled_dv_reg)) || force_remake
        println("calculate volume anomaly for regions, this will take ~11 sec")

        gemb = load(filename_gemb_geotile_filled_dv)
        drgi = Dim{:rgi}(reg)
        dpscale = dims(gemb[first(keys(gemb))], :pscale)
        dΔheight = dims(gemb[first(keys(gemb))], :Δheight)
        ddate = dims(gemb[first(keys(gemb))], :date)
        vars = collect(keys(gemb))

        # align geotile dataframe with DimArrays [this is overly cautious]
        gt = collect(dims(gemb[first(keys(gemb))], :geotile))
        gt_ind = [findfirst(geotiles.id .== g0) for g0 in gt]
        geotiles = geotiles[gt_ind, :]

        # create or extract dimensions
        gemb0 = Dict();
        for k in vars
            gemb0[k] = fill(NaN, (drgi, ddate, dpscale, dΔheight))
        end

        area_km2 = fill(0., (drgi))

        for rgi in drgi
            #rgi = first(drgi)
            index = geotiles[:, rgi] .> 0.0

            area_km2[At(rgi)] = sum(reduce(hcat, geotiles.area_km2[index]))

            for k in vars
            #k = "runoff"
                gemb0[k][At(rgi), :, :, :] = sum(gemb[k][index, :, :, :], dims = :geotile)
            end
        end

        gemb0["area_km2"] = area_km2;

        save(filename_gemb_geotile_filled_dv_reg, gemb0)
    end
end