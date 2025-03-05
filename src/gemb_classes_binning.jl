## NOTE: raw REMB output is in unit of m ice equivalent assuing an ice density of 910 kg/m3

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
#    - Calculates volume changes using glacier hysometry
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
#@time begin
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
        force_remake = true; # set force_remake == true to redo all steps from scratch 
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
        dpscale = Dim{:pscale}(sort(unique_precipitation_scale))

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

        printstyled("!!!!! FILLING GEMB GEOTILES ADD ADDING ΔHEIGHT CLASSES, THIS WILL TAKE ~4 HRS using 48 threads !!!!!\n", color=:red)
        
        # fill filename_gemb_geotile
        gemb = load(filename_gemb_geotile)
        k = first(keys(gemb))
        dgeotile = dims(gemb[k], :geotile)
        ddate = dims(gemb[k], :date)
        dheight = dims(gemb[k], :height)
        dpscale = dims(gemb[k], :pscale)
        x = collect(dheight);
        npts_linear_extrapolation = 7;

        # NOTE: smb and runoff are excluded from the filling and extrapolation as they are derivatives 
        # of melt, refreeze and acc. To maitain mass conservation "smb" and "runoff" are 
        # reonstucted from native variables after the filling and extrapolation is complete.
        vars = setdiff(collect(keys(gemb)), ["nobs", "smb", "runoff"]);
        
        # test subset of geotiles
        ### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ###
        if false
            gtids = ["lat[-14-12]lon[-076-074]", "lat[-10-08]lon[-078-076]", "lat[+28+30]lon[+082+084]", "lat[+58+60]lon[-136-134]", "lat[+68+70]lon[-070-068]", "lat[+62+64]lon[-044-042]", "lat[+32+34]lon[+076+078]", "lat[+28+30]lon[+092+094]", "lat[-46-44]lon[+168+170]", "lat[-70-68]lon[-068-066]", "lat[-52-50]lon[-074-072]"]
            
            #gtids = ["lat[-44-42]lon[+170+172]"]
            
            test_ind = [findfirst(isequal(gtid), dgeotile.val) for gtid in gtids]
            dgeotile = dgeotile[test_ind]
        end
        ### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ###

        # This take 36 minutes        
        @showprogress desc="Smoothing and extrapolating across all heights... "  Threads.@threads for geotile in dgeotile   
        #geotile = first(dgeotile)

            for k in vars
            #k = "fac"
                for pscale in dpscale
                #pscale = 1

                    # there is an issue with the interpolating the cumulitive variables... 
                    # interpolation needs to be done on rates to ensure they are physical
                    var0 = gemb[k][At(geotile), :, :, At(pscale)]

                    dheight = dims(var0, :height)

                    # for LSQ fitting
                    x0 = (1:length(dheight))./length(dheight) .- 0.5
                    M0 =  M = hcat(ones(length(x0)), x0);
                    npts = npts_linear_extrapolation; # just so I can use as shorthand

                    # sanity check
                    #CairoMakie.heatmap(var0)

                    # if all values are 0 then skip
                    if all(var0 .== 0)
                        continue
                    end

                    # calculate rates
                    if k == "fac"
                        # do not compute rates for interpolating FAC
                        var0_rate = var0[2:end,:];
                    else
                        ind = findfirst(.!vec(all(isnan.(var0), dims=:height))) - 1
                        var0[ind, :] .= 0
                        var0_rate = var0[2:end,:] .- collect(var0[1:end-1,:])
                    end

                    # sanity check
                    #CairoMakie.heatmap(var0_rate)

                    for date in dims(var0_rate, :date)
                    #date = ddate[920]

                        #println((; k, geotile, date, pscale))
                        y = var0_rate[At(date), :]
                        
                        # sanity check
                        #p = lines(y)
                
                        valid = .!isnan.(collect(y))
                    
                        if all(valid) || !any(valid)
                            continue
                        end
                    
                        
                        if k in ("acc", "refreeze", "melt", "rain", "fac")
                            # only positive rates are valid for these variables
                            # remember, "fac" is not a rate
                            valid = valid .& (y .>= 0)
                        end

                        if all(valid) || !any(valid)
                            y[:] .= 0
                        else
                            # smooth variables with height and interpolate gaps
                            validgap = Altim.validgaps(valid)

                            # NOTE: interplation increases per variable time from 1 min to 6 min 
                            if sum(valid) > 4
                                model = loess(x[valid], y[valid], span=0.3)
                                y[valid .| validgap]  = predict(model, x[valid .| validgap])
                            end

                            # sanity check
                            #lines!(y)

                            fill_index = isnan.(collect(y))
                            if (k == "ec")
                                valid_interp_index = .!fill_index
                            else
                                # do not include zeros in linear extrapolation as a funciton of elevation
                                valid_interp_index = .!fill_index .& (y .!= 0)
                            end

                            # linearly extrapolate the data
                            # simple linear extrapolation is crude but robust... 
                            # more precipitation and melt classes are needed in future 
                            if any(fill_index) &&  (sum(valid_interp_index) > 3)
                                if sum(valid_interp_index) < npts+2
                                    # apply a functional fit to extrapolate data [this does not work well for melt]
                                    M = @view M0[valid_interp_index,:]
                                    param1 = M\y[valid_interp_index]

                                    y[fill_index] = param1[1] .+ x0[fill_index].*param1[2]
              
                                else
                                    # split into upper and lower extrapolations
                                    fill_index_lower = copy(fill_index)
                                    fill_index_lower[findfirst(.!fill_index):end] .= false
                                    
                                    fill_index_upper = copy(fill_index)
                                    fill_index_upper[1:findfirst(.!fill_index_upper)] .= false

                                    valid_interp_index_lower = findall(valid_interp_index)[1:npts]
                                    valid_interp_index_upper = findall(valid_interp_index)[end-npts+1:end]

                                    if any(fill_index_lower)
                                        M = @view M0[valid_interp_index_lower,:]
                                        param1 = M\y[valid_interp_index_lower]

                                        y[fill_index_lower] = param1[1] .+ x0[fill_index_lower].*param1[2]
                                    end

                                    if any(fill_index_upper)
                                        M = @view M0[valid_interp_index_upper,:]
                                        param1 = M\y[valid_interp_index_upper]
                                        y[fill_index_upper] = param1[1] .+ x0[fill_index_upper].*param1[2]
                                    end
                                end

                                if k in ("acc", "refreeze", "melt", "rain", "fac")
                                    y[y.<0] .= 0
                                end

                            else
                                y[fill_index] .= 0
                            end

                            # sanity check
                            #lines!(y)
                            #p
                        end

                        gemb[k][At(geotile), At(date), :, At(pscale)] = y;
                    end
                end
            end
            #sanity check
            #for k in vars
            #    p = CairoMakie.heatmap(gemb[k][At(geotile), :, :, At(pscale)])
            #    display(p)
            #end      
        end

        # set some physcial constraints [I think this is redundent now]
        gemb = Altim.gemb_rate_physical_constraints!(gemb)
        gemb["fac"][gemb["fac"] .< 0] .= 0 #fac is handled seperately as it is not a rate

        # create a DD for volume change results
        gemb_dv = Dict()
        for k in keys(gemb)
            gemb_dv[k] = fill(NaN, (dgeotile, ddate, dpscale, dΔheight))
        end

        # add height classes and save to disk for each geotile [THIS TAKES 8 HRS]
        @showprogress desc="Adding Δheight classes [expect 4 on 48 threads]..." Threads.@threads for geotile in dgeotile
        #geotile = first(dgeotile)
            
            filename_gemb_geotile_filled = replace(filename_gemb_geotile, ".jld2" => "_$(geotile)_filled.jld2")

            if !isfile(filename_gemb_geotile_filled) || force_remake

                gemb0 = Dict()
                for k in vars
                    gemb0[k] = deepcopy(gemb[k][At(geotile), :, :, :])
                end

                gemb_gt = Dict()
                for k in vars
                    gemb_gt[k] = fill(NaN, (ddate, dheight, dpscale, dΔheight))
                end

                for k in vars
                #k = "melt"

                    if (k == "ec")
                        exclude_zeros_in_extrapolation = false
                    else
                        exclude_zeros_in_extrapolation = true
                    end

                    for Δheight in dΔheight
                    #Δheight = +2000
                
                        for pscale in dpscale
                        #pscale = 1
                            if k == "nobs"
                                gemb_gt[k][:,:,At(pscale),At(Δheight)] = gemb0[k][:, :, At(pscale)]
                            else              
                                # shift elevation
                                height_shift = round(Int16(Δheight / height_bin_interval))
   
                                gemb_gt[k][:,:,At(pscale),At(Δheight)] = Altim._matrix_shift_ud!(deepcopy(gemb0[k][:, :, At(pscale)]), height_shift; exclude_zeros_in_extrapolation, npts_linear_extrapolation)
                            end
                        end
                    end
                end

                gemb_gt = Altim.gemb_rate_physical_constraints!(gemb_gt)
                gemb["fac"][gemb["fac"] .< 0] .= 0 #fac is handled seperately as it is not a rate
                
                #sanity check
                #for k in vars
                #    p = CairoMakie.heatmap(gemb_gt[k][:, :, At(1), At(+200)]); display(p)
                #    p = CairoMakie.heatmap(gemb_gt[k][:, :, At(1), At(-200)]); display(p)
                #end  

                # transform back into cumulative variables
                for k in setdiff(collect(keys(gemb_gt)), ["nobs", "fac"])
                    sindex = findfirst(.!isnan.(gemb_gt[k][:, 1, 1, 1]))
                    gemb_gt[k][sindex:end, :, :, :] = cumsum(gemb_gt[k][sindex:end, :, :, :], dims=:date)
                end

                # save to disk
                # NOTE: files saved to disk are used calculating ice discahge for unmeasured glaciers 
                # AND for geotile2glacier downscaling if the hypsometric method is used
                save(filename_gemb_geotile_filled, gemb_gt)
            else
                println("direct read from $filename_gemb_geotile_filled")
                gemb_gt = load(filename_gemb_geotile_filled)
            end

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

#end

    # Try creating new `pscale` classes
   begin
        printstyled("!!!!! ADDING PSCALE CLASSES, THIS WILL TAKE ~4 HRS using 48 threads !!!!!\n", color=:red)
            
        dpscale_new = Dim{:pscale}(0.25:0.25:4.0)

        filename_gemb_geotile_filled = replace(filename_gemb_geotile, ".jld2" => "_$(geotiles.id[1])_filled.jld2")
        gembX = load(filename_gemb_geotile_filled)
        
        k = first(keys(gembX))
        dpscale = dims(gembX[k], :pscale)
        dΔheight = dims(gembX[k], :Δheight)
        ddate = dims(gembX[k], :date)
        dgeotile = Dim{:geotile}(geotiles.id)

        # create a DD for volume change results
        gemb_dv = Dict()
        for k in keys(gembX)
            gemb_dv[k] = fill(NaN, (dgeotile, ddate, dpscale_new, dΔheight))
        end 

        filename_gemb_geotile_filled_extra_dv = replace(filename_gemb_geotile, ".jld2" => "_filled_extra_dv.jld2")
        
        #index = Altim.within.(Ref(Altim.Extent(X=(-154,-131), Y=(57, 64))), mean.(getindex.(geotiles.extent,1)), mean.(getindex.(geotiles.extent,2)))
        #geotiles = geotiles[index, :]
        

        # Threads.@threads is being used lower down... 
        @showprogress desc = "Adding pscale classes [expect 4 on 48 threads]..."  for geotile in geotiles.id

            println("\n$geotile")

            filename_gemb_geotile_filled = replace(filename_gemb_geotile, ".jld2" => "_$(geotile)_filled.jld2")
            filename_gemb_geotile_filled_extra = replace(filename_gemb_geotile, ".jld2" => "_$(geotile)_filled_extra.jld2")
            
            if !isfile(filename_gemb_geotile_filled_extra)
                gemb = load(filename_gemb_geotile_filled)
                vars = setdiff(collect(keys(gemb)), ["nobs", "smb", "runoff"]);
        
                gemb_new = Dict()
                for k in vars
                    if k in ("acc", "refreeze", "melt", "rain", "fac")
                        allow_negative = false
                    else
                        allow_negative = true
                    end

                    # all variables are converted to rates except for fac
                    if k != "fac"
                        sind = findfirst(.!isnan.(gemb[k][:, 1, 1, 1]))-1
                        if sind > 2
                            gemb[k][findfirst(.!isnan.(gemb[k][:, 1, 1, 1]))-1, :, :, :] .= 0
                        end

                        gemb[k][2:end, :, :, :] = gemb[k][2:end, :, :, :] .- collect(gemb[k][1:end-1, :, :, :])

                        if sind > 2
                            gemb[k][findfirst(.!isnan.(gemb[k][:, 1, 1, 1]))-1, :, :, :] .= NaN
                        end
                    end

                    gemb_new[k] = Altim.add_pscale_classes(deepcopy(gemb[k][:, :, :, :]), dpscale_new; allow_negative=false)
                end

                gemb_new = Altim.gemb_rate_physical_constraints!(gemb_new)
                gemb_new["fac"][gemb_new["fac"] .< 0] .= 0 #fac is handled seperately as it is not a rate

                # transform back into cumulative variables
                for k in setdiff(collect(keys(gemb_new)), ["nobs", "fac"])
                    # restore to cumulative variable if it was converted to a rate
                    sindex = findfirst(.!isnan.(gemb_new[k][:, 1, 1, 1]))
                    gemb_new[k][sindex:end, :, :, :] = cumsum(gemb_new[k][sindex:end, :, :, :], dims=:date)
                end

                save(filename_gemb_geotile_filled_extra, gemb_new)
            else
                gemb_new = load(filename_gemb_geotile_filled_extra)
            end

            # create volume change time series for each geotile
            df_index = findfirst(isequal(geotile), geotiles.id)
            area = ones(length(ddate), 1) * hcat(geotiles[df_index, :area_km2])' ./ 1000 # divide by 1000 so that result in km3 (NOTE: vars are in mwe)
            for k in collect(keys(gemb_new))
                for pscale in dpscale_new
                #pscale = first(dpscale)            
                    for Δheight in dΔheight
                        if k == "nobs"
                            gemb_dv[k][At(geotile), :, At(pscale), At(Δheight)] =  sum(gemb_new[k][:,:,At(pscale),At(Δheight)], dims = :height);
                        else
                            gemb_dv[k][At(geotile), :, At(pscale), At(Δheight)] =  sum(gemb_new[k][:,:,At(pscale),At(Δheight)] .* area , dims = :height);
                        end
                    end
                end
            end
        end

        save(filename_gemb_geotile_filled_extra_dv, gemb_dv)
    end
end