# load packages
begin
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

    # run parameters
    force_remake = true;
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
    dΔheight = Dim{:Δheight}(-2000:height_bin_interval:2000) 

    vars = ["latitude", "longitude", "date", "smb", "fac", "acc", "runoff", "refreeze", "height", "melt", "rain", "ec"]


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
    filename_gemb_geotile_filled = replace(filename_gemb_geotile, ".jld2" => "_filled.jld2")
    filename_gemb_geotile_filled_dv = replace(filename_gemb_geotile, ".jld2" => "_filled_dv.jld2")
    filename_gemb_geotile_filled_dv_reg = replace(filename_gemb_geotile, ".jld2" => "_filled_d_reg.jld2")

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

# 2.7 min for all data
# merging all gemb data into a single file
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
    for var0 in setdiff(vars, ["date", "latitude", "longitude", "height"])
    #var0 = first(setdiff(vars, ["date", "latitude", "longitude", "height"]))
        foo[var0] = reduce(vcat, getindex.(gemb, Ref(var0)))
    end
    gemb = merge(foo, gemb0)

    # save output
    save(filename_gemb_combined, gemb);
end


# combine all interations into time-height geotiles that match the dh processing
# elevation_delta runs are used to fill in time-height space 
# precipitation_scale runs are treated as unique model runs

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
    for var0 in vars
        gemb0[var0] = fill(NaN, (dgeotile, ddate, dheight, dpscale))
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
            for var0 in vars
                gemb0[var0][At(geotile.id), :, :, :] .= 0.
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
                    for var0 in vars
                        gemb0[var0][At(geotile.id), :, At(height0), At(pscale)] = gemb[var0][index,:];
                    end
                else
                    for var0 in vars
                        gemb0[var0][At(geotile.id), :, At(height0), At(pscale)] = mean(gemb[var0][index,:], dims=1);
                    end
                end

                valid_date = .!isnan.(gemb0[vars[1]][At(geotile.id), :, At(height0), At(pscale)])
                gemb0["nobs"][At(geotile.id), valid_date, At(height0), At(pscale)] .= sum(index);
            end
        end
    end
    save(filename_gemb_geotile, gemb0)
end

# filling gemb geotiles
if .!(isfile(filename_gemb_geotile_filled)) || force_remake
    println("filling gemb geotiles, this will take ~2 min")
    # fill filename_gemb_geotile
    gemb = load(filename_gemb_geotile)
    dgeotile = dims(gemb[first(keys(gemb))], :geotile)
    ddate = dims(gemb[first(keys(gemb))], :date)
    dheight = dims(gemb[first(keys(gemb))], :height)
    dpscale = dims(gemb[first(keys(gemb))], :pscale)
    x = collect(dheight);
    vars = setdiff(collect(keys(gemb)), ["nobs","smb"]);
    
    Threads.@threads for k in collect(vars)
    #k = first(keys(gemb))

        for geotile in dgeotile   
        #geotile = first(dgeotile)
         
            for date in ddate
            #date = first(ddate)

                for pscale in dpscale
                    #pscale = first(dpscale)
                    #println((; k, geotile, date, pscale))
                    y = gemb[k][At(geotile), At(date), :, At(pscale)]
                    notnan = .!isnan.(collect(y))

                    if all(notnan) || !any(notnan)
                        continue
                    end

                    interp = LinearInterpolation(collect(y[notnan]), x[notnan]; extrapolate = true)
                    y[.!notnan] = interp(x[.!notnan])

                    # set some physcial constraints
                    if any(k .== ["runoff", "rain", "acc", "melt", "refreeze"])
                        y[y.<0] .= 0
                    end

                    gemb[k][At(geotile), At(date), :, At(pscale)] = y;
                end
            end
        end
    end

    # SMB must equal acc - runoff - ec
    k = "smb"
    
    for geotile in dgeotile   
    #geotile = first(dgeotile)
        
        for date in ddate
        #date = first(ddate)

            for pscale in dpscale
                
                acc = gemb["acc"][At(geotile), At(date), :, At(pscale)]
                runoff = gemb["runoff"][At(geotile), At(date), :, At(pscale)]
            
                gemb[k][At(geotile), At(date), :, At(pscale)] = acc .- runoff;
            end
        end
    end

    save(filename_gemb_geotile_filled, gemb)
end


# multiply glacier area to get volume change.
if .!(isfile(filename_gemb_geotile_filled_dv)) || force_remake
    println("calculate volume anomaly for gemb geotiles, this will take ~0.5 min")

    gemb = load(filename_gemb_geotile_filled)
    
    ddate = dims(gemb[first(keys(gemb))], :date)
    dpscale = dims(gemb[first(keys(gemb))], :pscale)
    dgeotile = dims(gemb[first(keys(gemb))], :geotile)
    vars = collect(keys(gemb))
    
    # align geotile dataframe with DimArrays [this is overly cautious]
    gt = collect(dims(gemb[first(keys(gemb))], :geotile))
    gt_ind = [findfirst(geotiles.id .== g0) for g0 in gt]
    geotiles = geotiles[gt_ind, :]

    gemb_dv = Dict()
    for var0 in vars
        gemb_dv[var0] = fill(NaN, (dgeotile, ddate, dpscale, dΔheight))
    end

    # this takes only 0.5 min using @threads
    Threads.@threads for geotile in eachrow(geotiles)
    #geotile = geotiles[604,:]
        
        area = ones(length(ddate), 1) * hcat(geotile.area_km2)';
        for Δheight in dΔheight
        #Δheight = dΔheight[7]

            for var0 in vars
            #var0 = "runoff"
                
                for pscale in dpscale
                #pscale = dpscale[5]
     
                    if var0 == "nobs"
                        # do not multiply by area
                        v0 = gemb[var0][At(geotile.id), :, :, At(pscale)]
                        gemb_dv[var0][At(geotile.id), :, At(pscale), At(Δheight)] =  sum(v0, dims = :height);
                    else

                        v0 = gemb[var0][At(geotile.id), :, :, At(pscale)] ./ 1000
                        
                        # shift elevation
                        dh = round(Int16(Δheight / height_bin_interval))
                        if dh > 0
                            v0[:,1:end-dh] = v0[:,1+abs(dh):end]
                        else
                            v0[:, -dh+1:end] = v0[:,1:end+dh]
                        end

                        gemb_dv[var0][At(geotile.id), :, At(pscale), At(Δheight)] =  sum(v0 .* area, dims = :height);
                    end
                end
            end
        end
    end

    save(filename_gemb_geotile_filled_dv, gemb_dv)
end


# now create regional volume metrics (sum across elevation bands)
# multiply glacier area to get volume change.
if .!(isfile(filename_gemb_geotile_filled_dv_reg)) || force_remake
    println("calculate volume anomaly for regions, this will take ~11 sec")

    gemb = load(filename_gemb_geotile_filled_dv)
    drgi = Dim{:rgi}(reg)
    dpscale = dims(gemb[first(keys(gemb))], :pscale)
    dΔheight = dims(gemb[first(keys(gemb))], :Δheight)
    ddate = dims(gemb[first(keys(gemb))], :date)
    vars = keys(gemb)

    # align geotile dataframe with DimArrays [this is overly cautious]
    gt = collect(dims(gemb[first(keys(gemb))], :geotile))
    gt_ind = [findfirst(geotiles.id .== g0) for g0 in gt]
    geotiles = geotiles[gt_ind, :]

    # create or extract dimensions
    gemb0 = Dict();
    for var0 in vars
        gemb0[var0] = fill(NaN, (drgi, ddate, dpscale, dΔheight))
    end

    area_km2 = fill(0., (drgi))


    for rgi in drgi
        #rgi = first(drgi)
        index = geotiles[:, rgi] .> 0.0

        area_km2[At(rgi)] = sum(reduce(hcat, geotiles.area_km2[index]))

        for var0 in vars
        #var0 = "runoff"
            gemb0[var0][At(rgi), :, :, :] = sum(gemb[var0][index, :, :, :], dims = :geotile)
        end
    end

    gemb0["area_km2"] = area_km2;

    save(filename_gemb_geotile_filled_dv_reg, gemb0)
end
end