# load packages
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

    # run parameters
    force_remake = true;
    project_id = :v01;
    geotile_width = 2;
    binning_method = "mean"; # "meanmadnorm3"
    surface_mask = :glacier

    # define date and hight binning ranges 
    Δd = 30
    date_range = DateTime(1990):Day(Δd):DateTime(2026, 1, 1)
    date_center = date_range[1:end-1] .+ Day(Δd / 2)

    if  false
        gemb_folder = ["/home/schlegel/Share/GEMBv1/"];
        file_uniqueid = "rv1_0_19500101_20231231"
        elevation_delta = [0]
        precipitation_scale = [1]
        gemb_file_binned = "/mnt/bylot-r3/data/gemb/raw/$file_uniqueid.jld2"

    elseif false
        gemb_folder = "/home/schlegel/Share/GEMBv1/Alaska_sample/v1/"
        file_uniqueid = "1979to2023_820_40_racmo_grid_lwt"
        elevation_delta = [-1000 -750 -500 -250 0 250 500 750 1000]
        precipitation_scale = [0.5 1 1.5 2 5 10]
        gemb_file_binned = "/mnt/bylot-r3/data/gemb/raw/$file_uniqueid.jld2"

    else
        gemb_folder = ["/home/schlegel/Share/GEMBv1/NH_sample/", "/home/schlegel/Share/GEMBv1/SH_sample/"]
        file_uniqueid = nothing
        elevation_delta = [-200, 0, 200, 500, 1000]
        precipitation_scale = [.75 1 1.25 1.5 2]
        gemb_file_binned = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0.jld2"

    end

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

    # DO NOT MODIFY: these need to match elevations in hypsometry
    Δh = 100;
    height_range = 0:100:10000;
    height_center = height_range[1:end-1] .+ Δh/2;

    # 3.2 hours for all glaciers, all missions/datasets on 96 threads
    # 45 min for icesat and icesat-2 missions only
    # 2.3 hours for hugonnet only
    showplots = false;
end;

# 30s for all data
@time begin
    # load geotile definitions with corresponding hypsometry
    geotiles = Altim.geotiles_mask_hyps(surface_mask, geotile_width)

    # filter geotiles
    geotiles = geotiles[geotiles.glacier_frac.>0.0, :]

    if .!(isfile(gemb_file_binned)) || force_remake   
        gemb_files = allfiles.(gemb_folder; subfolders=false, fn_endswith=".mat", fn_contains=file_uniqueid)
        gemb_files = vcat(gemb_files...)

        expected_num_files = length(elevation_delta) * length(precipitation_scale) * length(gemb_folder)
        #if length(gemb_files) !== expected_num_files
        #    error("the number of files does not match the number of classes")
        #end

        vars = ["latitude","longitude", "datetime", "smb", "fac", "height_ref"]
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
                if k !== "datetime"
                    gemb0[k] = reduce(vcat,[g[k] for g in gemb])
                else
                    gemb0[k] = gemb[1][k] 
                end
            end
            gemb = [gemb0]
        end

        # merge all simulations into a single directory
        datetime = gemb[1]["datetime"]
        latitude = reduce(vcat, getindex.(gemb, Ref("latitude")))
        longitude = reduce(vcat, getindex.(gemb, Ref("longitude")))

        # make sure longitude is in -180 to 180 range
        index = longitude .< -180
        longitude[index] = longitude[index] .+ 360
        index = longitude .> 180
        longitude[index] = longitude[index] .- 360

        smb = reduce(vcat, getindex.(gemb, Ref("smb")))
        fac = reduce(vcat, getindex.(gemb, Ref("fac")))

        precipitation_scale = []
        elevation_delta = []
        for g in gemb
            precipitation_scale = reduce(vcat, fill(g["precipitation_scale"], size(g["latitude"])))
            elevation_delta = reduce(vcat, fill(g["elevation_delta"], size(g["latitude"])))
        end

        gemb = Dict(
            "datetime" => datetime,
            "latitude" => latitude, 
            "longitude" => longitude, 
            "smb" => smb, 
            "fac" => fac, 
            "precipitation_scale" =>  precipitation_scale, 
            "elevation_delta" =>  elevation_delta
            )

        # save output
        save(gemb_file_binned, Dict("gemb" => gemb));
    end
end