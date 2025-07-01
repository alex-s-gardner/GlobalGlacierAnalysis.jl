#    gemb_classes_binning.jl
#
#Process GEMB (Glacier Energy and Mass Balance) model output for global glacier analysis.

# This script:
# 1. Loads and combines raw GEMB output from multiple simulations
# 2. Organizes data into geotiles with consistent spatial and temporal dimensions
# 3. Fills data gaps through interpolation and extrapolation across elevation bands
# 4. Extends model results to cover additional precipitation scaling factors
# 5. Extrapolates data temporally using climatological means

# The processed data is saved at multiple stages to enable efficient reuse and analysis.
# GEMB output is in units of meters ice equivalent (m i.e.) assuming an ice density of 910 kg/m³.

#begin
begin
    import GlobalGlacierAnalysis as GGA
    using Arrow
    using DataFrames
    using Extents
    using Rasters
    using GeoInterface
    using Shapefile
    using DimensionalData
    using Statistics
    using Dates
    using JLD2
    using FileIO
    using MAT
    using GeoTiles
    using DimensionalData
    using DataInterpolations
    using ProgressMeter
    using Loess
    using NaNStatistics

    # run parameters 
    force_remake_before = DateTime(2025, 6, 23, 9, 0, 0)
    project_id = :v01;
    geotile_width = 2;
    binning_method = "mean";
    surface_mask = :glacier
    geotile_buffer = 50000 # distance in meters outside of geotiles to look for data
    showplots = false;
    gemb_run_id=4;

     # exclude derived variables of smb and runoff (these are calculated later to ensure mass conservation after interpolation)
    vars = ["latitude", "longitude", "date", "fac", "acc", "refreeze", "height", "melt", "rain", "ec"]

    # min gemb coverage 
    min_gemb_coverage =  0.75

    @warn "IT TAKES SEVERAL DAYS, POSSIBLY A WEEK, TO PREPROCESS ALL GEMB DATA"
end

begin
    gembinfo = GGA.gemb_info(; gemb_run_id)

    # define date and hight binning ranges 
    date_range, date_center = GGA.project_date_bins()

    # expand daterange to 1940 by make sure to match exisiting project ranges 
    foo = collect(date_range[1]:-Day(30):DateTime(1940))

    date_range = foo[end]:Day(30):date_range[end]
    date_center = date_range[1:end-1] .+ Day(Day(30) / 2)

    height_range, height_center = GGA.project_height_bins()

    #Δheight simulates changing model elevation to increase / decrease melt, this is done in the regional Δvolume calculation
    height_bin_interval = height_center[2] - height_center[1]
    dΔheight = Dim{:Δheight}(-2000:height_bin_interval:2000) 


    filename_gemb_geotile = replace(gembinfo.filename_gemb_combined, ".jld2" => "_geotile.jld2")
    filename_gemb_geotile_filled_dv = replace(filename_gemb_geotile, ".jld2" => "_filled_dv.jld2")
    filename_gemb_geotile_filled_extra_dv = replace(filename_gemb_geotile, ".jld2" => "_filled_extra_dv.jld2")
    filename_gemb_geotile_filled_extra_extrap_dv = replace(filename_gemb_geotile, ".jld2" => "_filled_extra_extrap_dv.jld2")

    binfun = GGA.binningfun_define(binning_method)::Function
    
    # load geotile definitions with corresponding hypsometry
    geotiles = GGA.geotiles_mask_hyps(surface_mask, geotile_width)

    # filter geotiles
    geotiles = geotiles[geotiles.glacier_frac.>0.0, :]

    # make geotile rgi regions mutually exexclusive 
    geotiles, reg = GGA.geotiles_mutually_exclusive_rgi!(geotiles)

    if !any(Base.contains.("area_km2", names(geotiles)))
        original_area_name = string(surface_mask)*"_area_km2"
        generic_area_name = "area_km2"
        rename!(geotiles, original_area_name => generic_area_name)
    end
end;

# Merge GEMB model data from multiple files into a single combined file.
# 
# This function processes GEMB (Glacier Energy and Mass Balance) model data by:
# 1. Collecting all relevant .mat files from specified folders
# 2. Reading and combining data across different elevation delta and precipitation scale classes
# 3. Standardizing longitude values to -180 to 180 range
# 4. Merging all simulations into a single dictionary with consistent structure
# 5. Saving the combined data to a JLD2 file
# 
# The process takes approximately 3 minutes to complete for the full dataset.
# Takes ~3 minutes to run for all data.

if isfile(gembinfo.filename_gemb_combined) && ((Dates.unix2datetime(mtime(gembinfo.filename_gemb_combined)) > force_remake_before) || isnothing(force_remake_before))
    printstyled("\n    -> Skipping $(gembinfo.filename_gemb_combined) as it was created after the force_remake_before date: $force_remake_before \n"; color=:light_green)
else
    println("merging all gemb data into a single file, this will take ~3 min")

    gemb_files = GGA.allfiles.(gembinfo.gemb_folder; subfolders=false, fn_endswith=".mat", fn_contains=gembinfo.file_uniqueid)
    gemb_files = vcat(gemb_files...)

    expected_num_files = length(gembinfo.elevation_delta) * length(gembinfo.precipitation_scale) * length(gembinfo.gemb_folder)
    #if length(gemb_files) !== expected_num_files
    #    error("the number of files does not match the number of classes")
    #end

    datebin_edges = GGA.decimalyear.(date_range)
    # this takes 4 m
    gemb = []
    
    Threads.@threads for gemb_file in gemb_files
    #gemb_file = first(gemb_files)

        gemb0 = GGA.gemb_read2(gemb_file; vars, datebin_edges)

        # no classes
        if (length(gembinfo.elevation_delta) == 1) .&  (length(gembinfo.precipitation_scale) == 1)
            precip_scale_ind = 1
            elev_delta_ind = 1
        else
            ind1 = findlast("_p", gemb_file)
            ind2 = findlast("_t", gemb_file)
            ind3 = findlast(".mat", gemb_file)
        
            precip_scale_ind = parse(Int8, gemb_file[ind1[end]+1 : ind2[1]-1])
            elev_delta_ind = parse(Int8, gemb_file[ind2[end]+1 : ind3[1]-1])
        end

        gemb0["elevation_delta"] = gembinfo.elevation_delta[elev_delta_ind]
        gemb0["precipitation_scale"] = gembinfo.precipitation_scale[precip_scale_ind]

        push!(gemb, gemb0)
    end

    if (length(gembinfo.elevation_delta) == 1) .&  (length(gembinfo.precipitation_scale) == 1)
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
        "date" => GGA.decimalyear2datetime.(vec(dates)),
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
    save(gembinfo.filename_gemb_combined, gemb);
end

# Process GEMB data into geotiles with elevation binning
#
# Takes raw GEMB model output and organizes it into geotiles with consistent spatial,
# temporal and elevation dimensions. For each geotile, data is binned by elevation class
# and precipitation scaling factor. The function continues to buffer around each geotile
# until sufficient data coverage is achieved.
#
# Processing time is approximately 3.5 minutes with multi-threading.
# Takes ~3.5 minutes to run
if isfile(filename_gemb_geotile) && ((Dates.unix2datetime(mtime(filename_gemb_geotile)) > force_remake_before) || isnothing(force_remake_before))
    printstyled("    -> Skipping $(filename_gemb_geotile) as it was created after the force_remake_before date: $force_remake_before \n"; color=:light_green)
else
    println("putting gemb data into geotiles, this will take ~3.5 min")
    # combine into geotiles
    gemb = load(gembinfo.filename_gemb_combined);

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
        latitude_distance, longitude_distance = GGA.meters2lonlat_distance(geotile_buffer, mean(ext.Y))
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
            in_geotile = vec([GGA.within(ext, x, y) for (x, y) in zip(gemb["longitude"], gemb["latitude"])])
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



# Fill gaps in GEMB geotile data and add elevation change (Δheight) classes.
#
# This function processes GEMB geotile data by:
# 1. Filling data gaps across elevation bands using interpolation
# 2. Extrapolating data to cover additional elevation ranges
# 3. Creating new datasets with various elevation change scenarios (Δheight)
# 4. Calculating volume changes for each geotile, precipitation scale, and Δheight
#
# The process involves smoothing and extrapolating data across elevation bands,
# applying physical constraints to ensure realistic values, and converting rates
# back to cumulative variables. Results are saved at multiple stages.
#
# Note: Linear extrapolation over long distances is not ideal and requires
# running GEMB for a wider range of elevation classes to improve results.

if isfile(filename_gemb_geotile_filled_dv) && ((Dates.unix2datetime(mtime(filename_gemb_geotile_filled_dv)) > force_remake_before) || isnothing(force_remake_before))
    printstyled("    -> Skipping $(filename_gemb_geotile_filled_dv) as it was created after the force_remake_before date: $force_remake_before \n"; color=:light_green)
else

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
                        validgap = GGA.validgaps(valid)

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
    gemb = GGA.gemb_rate_physical_constraints!(gemb)
    gemb["fac"][gemb["fac"] .< 0] .= 0 #fac is handled seperately as it is not a rate

    # create a DD for volume change results
    gemb_dv = Dict()
    for k in keys(gemb)
        gemb_dv[k] = fill(NaN, (dgeotile, ddate, dpscale, dΔheight))
    end

    # add height classes and save to disk for each geotile [THIS TAKES 8 HRS]
    # TODO: I should look into moving the @threads lower down to minimize memory usage
    @showprogress desc="Adding Δheight classes [expect 4 hrs on 128 threads & consumes > 900GB of memory]..." Threads.@threads for geotile in dgeotile
    #geotile = first(dgeotile)
        
        filename_gemb_geotile_filled = replace(filename_gemb_geotile, ".jld2" => "_$(geotile)_filled.jld2")

        if !isfile(filename_gemb_geotile_filled) || (!isnothing(force_remake_before) && Dates.unix2datetime(mtime(filename_gemb_geotile_filled)) < force_remake_before)

            printstyled("    -> Creating $(filename_gemb_geotile_filled) \n"; color=:light_black)

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

                            gemb_gt[k][:,:,At(pscale),At(Δheight)] = GGA._matrix_shift_ud!(deepcopy(gemb0[k][:, :, At(pscale)]), height_shift; exclude_zeros_in_extrapolation, npts_linear_extrapolation)
                        end
                    end
                end
            end

            gemb_gt = GGA.gemb_rate_physical_constraints!(gemb_gt)
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
            printstyled("    -> Loading from disk $(filename_gemb_geotile_filled_extra) as it was created after the force_remake_before date: $force_remake_before \n"; color=:light_green)
            gemb_gt = load(filename_gemb_geotile_filled)
        end
    
        # create volume change time series for each geotile
        df_index = findfirst(isequal(geotile), geotiles.id)
        area = ones(length(ddate), 1) * hcat(geotiles[df_index, :area_km2])' ./ 1000 # divide by 1000 so that result in km3 (NOTE: vars are in mie)
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



# rm("/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2024_820_40_racmo_grid_lwt_e97_0_geotile_lat[+76+78]lon[-036-034]_filled_extra.jld2")


# Add intermediate precipitation scaling classes to GEMB geotile data.
#
# This function:
# 1. Creates a finer resolution precipitation scaling grid (0.25 to 4.0 in 0.25 increments)
# 2. Interpolates existing GEMB data to these new precipitation scaling values
# 3. Applies physical constraints to ensure realistic values
# 4. Converts rates back to cumulative variables
# 5. Calculates volume changes for each geotile, precipitation scale, and elevation change
#
# Processing time is approximately 4 hours using 48 threads.
# Takes approximately 4 hours to run using 48 threads.

if isfile(filename_gemb_geotile_filled_extra_dv) && ((Dates.unix2datetime(mtime(filename_gemb_geotile_filled_extra_dv)) > force_remake_before) || isnothing(force_remake_before))
    printstyled("    -> Skipping $(filename_gemb_geotile_filled_extra_dv) as it was created after the force_remake_before date: $force_remake_before \n"; color=:light_green)
else
    printstyled("!!!!!Adding pscale classes [expect > 1 day on 128 threads, consumes > 1000GB of memory] !!!!!\n", color=:red)
        
    dpscale_new = Dim{:pscale}(0.25:0.25:4.0)

    filename_gemb_geotile_filled = replace(filename_gemb_geotile, ".jld2" => "_$(geotiles.id[1])_filled.jld2")
    gembX = load(filename_gemb_geotile_filled)
    
    k = first(keys(gembX))
    dpscale = dims(gembX[k], :pscale)
    dΔheight = dims(gembX[k], :Δheight)
    ddate = dims(gembX[k], :date)
    dgeotile = Dim{:geotile}(geotiles.id)
    dheight = dims(gembX[k], :height)
    
    # create a DD for volume change results
    gemb_dv = Dict{String, DimArray{Float64, 4}}()
    for k in keys(gembX)
        gemb_dv[k] = fill(NaN, (dgeotile, ddate, dpscale_new, dΔheight))
    end 

    # @threads is used lower down to minimize memory usage [was killing Devon]
    @showprogress desc = "Adding pscale classes [expect > 1 day on 128 threads, consumes > 100GB of memory]..."  for geotile in geotiles.id
    #geotile = first(geotiles.id)
    
        filename_gemb_geotile_filled = replace(filename_gemb_geotile, ".jld2" => "_$(geotile)_filled.jld2")
        filename_gemb_geotile_filled_extra = replace(filename_gemb_geotile, ".jld2" => "_$(geotile)_filled_extra.jld2")
        
        if !isfile(filename_gemb_geotile_filled_extra) ||  (!isnothing(force_remake_before) && Dates.unix2datetime(mtime(filename_gemb_geotile_filled_extra)) < force_remake_before)

            printstyled("\n    -> Creating $(filename_gemb_geotile_filled_extra) \n"; color=:light_black)

            gemb = load(filename_gemb_geotile_filled)
            vars = setdiff(collect(keys(gemb)), ["nobs", "smb", "runoff"]);
    
            gemb_new = Dict()

            # would be better to use @threads here ?
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

                gemb_new[k] = GGA.add_pscale_classes(deepcopy(gemb[k][:, :, :, :]), dpscale_new; allow_negative=false)
            end

            gemb_new = GGA.gemb_rate_physical_constraints!(gemb_new)
            gemb_new["fac"][gemb_new["fac"] .< 0] .= 0 #fac is handled seperately as it is not a rate

            # transform back into cumulative variables
            for k in setdiff(collect(keys(gemb_new)), ["nobs", "fac"])
                # restore to cumulative variable if it was converted to a rate
                sindex = findfirst(.!isnan.(gemb_new[k][:, 1, 1, 1]))
                gemb_new[k][sindex:end, :, :, :] = cumsum(gemb_new[k][sindex:end, :, :, :], dims=:date)
            end

            save(filename_gemb_geotile_filled_extra, gemb_new)
        else
            printstyled("\n    -> Loading from disk $(filename_gemb_geotile_filled_extra) as it was created after the force_remake_before date: $force_remake_before \n"; color=:light_green)
            gemb_new = load(filename_gemb_geotile_filled_extra)
        end

        # create volume change time series for each geotile
        df_index = findfirst(isequal(geotile), geotiles.id)
        area = DimArray(vec(hcat(geotiles[df_index, :area_km2])' ./ 1000), dheight) # divide by 1000 so that result in km3 (NOTE: vars are in mie)
        Threads.@threads for k in keys(gemb_new)
            if k == "nobs"
                gemb_dv[k][geotile = At(geotile)] =  dropdims(sum(gemb_new[k], dims = :height), dims = :height);
            else
                var0 = (@d gemb_new[k] .* area)
                gemb_dv[k][geotile = At(geotile)] = dropdims(sum(var0, dims = :height), dims = :height);
            end
        end
    end

    save(filename_gemb_geotile_filled_extra_dv, gemb_dv)
end




# Extrapolate GEMB data beyond the last valid date using monthly climatology.
#
# This function:
# 1. Loads volume change time series data for each geotile
# 2. Calculates monthly differences between time steps
# 3. Computes the mean monthly climatology from available data
# 4. Extends the time series beyond the last valid date using this climatology
# 5. Ensures continuity by adding the last valid value to the extrapolated differences
# 6. Saves the extended dataset to a new file
#
# Processing time is approximately 10 minutes.
# Takes approximately 10 min to run.

if isfile(filename_gemb_geotile_filled_extra_extrap_dv) && (isnothing(force_remake_before) || Dates.unix2datetime(mtime(filename_gemb_geotile_filled_extra_extrap_dv)) > force_remake_before)
    printstyled("\n    -> Skipping $(filename_gemb_geotile_filled_extra_extrap_dv) as it was created after the force_remake_before date: $force_remake_before \n"; color=:light_green)
else
    printstyled("\n    -> Creating $(filename_gemb_geotile_filled_extra_extrap_dv) \n"; color=:light_black)
    # repeat mean climatology for later years
    gemb_dv = load(filename_gemb_geotile_filled_extra_dv)

    for k in keys(gemb_dv)
        #k = first(keys(gemb_dv))
        v = gemb_dv[k]

        # take the difference between between time steps
        v_diff = diff(v; dims=:date)

        (dgeotile, ddate, dpscale, dΔheight) = DimensionalData.dims(v)
    
        # `diff` shifts time to the left ... instead we want to shift time to the right
        v_diff = DimArray(parent(v_diff), (dgeotile, ddate[2:end], dpscale, dΔheight))

        # calculate climatology for each month
        v_diff_monthly = groupby(v_diff, :date => month)
        v_diff_mean = cat([nanmean(M; dims=dimnum(M, :date)) for M in v_diff_monthly]...; dims = 2)
        dmonth = Dim{:month}(val(dims(v_diff_monthly, :date)))
        v_diff_mean = DimArray(v_diff_mean, (dgeotile, dmonth, dpscale, dΔheight))

        gemb_last_valid_date = ddate[findlast(.!isnan.(v_diff[1,:,1,1]))];
        #fill_index = isnan.(v_diff)
        #fill_index[:,findall(ddate .< gemb_last_valid_date),:,:] .= false

        v_fill = v_diff_mean[:, At(month.(ddate)), :, :]
        v_fill = DimArray(parent(v_fill), (dgeotile, ddate, dpscale, dΔheight))

        #v_diff[fill_index] = v_fill[fill_index]

        v_diff = cumsum(v_fill[:,DimensionalData.Between(gemb_last_valid_date,last(ddate)),:,:]; dims =:date)

        for d in dims(v_diff, :date)
            v_diff[:, At(d),:,:] .+= parent(v[:,At(gemb_last_valid_date),:,:])
        end
        
        v[:,DimensionalData.Between(gemb_last_valid_date,last(ddate)),:,:] = v_diff
            
        # sanity check
        # plot(v[100,:,1,1])
        #plot(gemb_dv[k][100,:,1,1])
    end
    
    save(filename_gemb_geotile_filled_extra_extrap_dv, gemb_dv)
end
end