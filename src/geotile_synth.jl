begin
using Altim
using FileIO
using DimensionalData
using Statistics
using Dates
using LsqFit
using Distances
using Plots
using DataFrames
using Rasters
using GeoDataFrames
using GeometryOps
using BinStatistics
using GeoFormatTypes
using StatsBase
using DataInterpolations
#using ProfileView

#include("mapzonal.jl")

paths = Altim.pathlocal
geotile_width = 2;

# to include in uncertainty
project_id = ["v01"]
surface_mask=["glacier", "glacier_rgi7", "glacier_b1km"]
dem_id=["best" "cop30_v2"]
curvature_correct=[false true]
amplitude_correct=[true]
binning_method=["median" "meanmadnorm5" "meanmadnorm3" "meanmadnorm10"] # ["median" "meanmadnorm10" "meanmadnorm5" "meanmadnorm3"]
paramater_set=[1, 2, 3, 4]
binned_folder=("/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg")

filename_gemb_combined = "/mnt/bylot-r3/data/gemb/raw/FAC_forcing_glaciers_1979to2023_820_40_racmo_grid_lwt_e97_0.jld2"

# for latitudes below this set discharge2smb
discharge2smb_max_latitude = -60;
discharge2smb_equilibrium_period = (Date(1979), Date(2000))

param_nt = (;project_id, surface_mask, dem_id, curvature_correct, amplitude_correct, binning_method, paramater_set, binned_folder)
params = Altim.ntpermutations(param_nt)

# only include files that exist
path2runs = String[]
for param in params
    binned_aligned_file = Altim.binned_aligned_filepath(; param...)
    if isfile(binned_aligned_file)
        push!(path2runs, binned_aligned_file)
    end
end

# calculate individual glacier hypsometry
glacier_geom = Dict();
geomfile_rgi6 = joinpath(paths.data_dir, "GlacierOutlines/rgi60/rgi60_Global.gpkg")
geomfile_rgi7 = joinpath(paths.data_dir, "GlacierOutlines/RGI2000-v7.0-G-global-fix/rgi70_Global.gpkg")

for geomfile_rgi = [geomfile_rgi6, geomfile_rgi7]
    glacier_geotile_hyps_fn = replace(geomfile_rgi, ".gpkg" => "geotile_hyps.jld2")

    if Base.contains(geomfile_rgi, "rgi6")
        persistent_attribute = :RGIId
    else
        persistent_attribute = :rgi_id
    end

    if !isfile(glacier_geotile_hyps_fn)
        glacier_geom = GeoDataFrames.read(geomfile_rgi)
        h = Raster(paths.cop30_v2, lazy=true)
        height_range, height_center = Altim.project_height_bins()
        geotiles = Altim.geotiles_w_mask(geotile_width)
        geotiles = geotiles[(geotiles.glacier_frac.>0.0), :]

        # this takes 37 min for all tiles
        @time geotile_zonal_area_hyps(h, height_range, glacier_geom, geotiles.id; persistent_attribute)
        FileIO.save(glacier_geotile_hyps_fn, Dict("glaciers" => glacier_geom))
    end
end

# measured dischage and dischage for unmeasured areas by setting dischage == unscaled smb
globaldischage_fn =joinpath(paths.data_dir, "GlacierOutlines/GlacierDischarge/global_glacier_dischage.jld2")
if .!isfile(globaldischage_fn)
    glacier_geotile_hyps_fn = joinpath(paths.data_dir, "GlacierOutlines/rgi60/rgi60_Global.gpkg")
    glacier_geotile_hyps_fn = replace(glacier_geotile_hyps_fn, ".gpkg" => "geotile_hyps.jld2")
    glaciers = load(glacier_geotile_hyps_fn, "glaciers")

    filename_gemb_geotile_filled = replace(filename_gemb_combined, ".jld2" => "_geotile_filled.jld2")

    gemb = Dict()
    gemb["smb"] = load(filename_gemb_geotile_filled, "smb")

    geotiles = Altim.geotiles_w_mask(geotile_width)
    geotiles = geotiles[(geotiles.glacier_frac.>0.0), :]
    geotiles[!, :pscale] .= 1.0 # no prcipitaiton scaling
    geotiles[!, :Δheight] .= 0.0

    @time glaciers = Altim.geotile2glacier!(glaciers, gemb, geotiles)
    discharge0 = Altim.discharge2smb(glaciers[sum.(glaciers.area_km2).>0, :]; discharge2smb_max_latitude, discharge2smb_equilibrium_period)

    # append to measured dischage
    discharge = Altim.glacier_discharge(; datadir=Altim.pathlocal[:data_dir])
    discharge = vcat(discharge, discharge0)
    discharge[discharge.discharge_gtyr.<0, :discharge_gtyr] .= 0;

    save(globaldischage_fn, Dict("discharge"=>discharge))
else
    discharge = load(globaldischage_fn, "discharge")
end


function geotile_synthesis_error(;
    path2runs,
    outfile = "/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2",
    minimum_error = 0.05,
    force_remake = false
    )

    #outfile = "/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2"

    if !isfile(outfile) || force_remake

        #load example file to get dimension and mission info
        dh = FileIO.load(path2runs[1], "dh_hyps");

        missions = collect(keys(dh))
        ddate = dims(dh[missions[1]], :date)
        dheight = dims(dh[missions[1]], :height)
        dgeotile = dims(dh[missions[1]], :geotile)
        dfile = Dim{:file}(path2runs)

        dh_all = Dict()
        for mission in missions
            dh_all[mission] = fill(NaN, (dfile, dgeotile, ddate, dheight))
        end

        # this takes about 8 min for length(files) = 96
        Threads.@threads for filepath in path2runs
            dh = FileIO.load(filepath, "dh_hyps");
            for mission in missions
                dh_all[mission][At(filepath),:,:,:] = dh[mission]
            end
        end

        # calculate standard deviation across runs as an error metric
        dh_all_std = Dict()
        for mission in missions
            dh_all_std[mission] = fill(NaN, (dgeotile, ddate, dheight))
        end

        for mission in missions
        #mission = "gedi"
        
            Threads.@threads for geotile in dgeotile
            #geotile = dgeotile[510]
                
                # !!! I THINK THIS IS FIXED NOW !!!!
                # TODO: There is an issue where some missions (ICESat2) has different hight range than the other missions... this cuases issues with non-rectangular error matrix
                # See this example
                # heatmap(dh_err["hugonnet"][At("lat[-28-26]lon[-070-068]"),:,:])
                # heatmap!(dh_err["icesat2"][At("lat[-28-26]lon[-070-068]"),:,:])

                valid1 = dropdims(any(.!isnan.(dh_all[mission][:, At(geotile), :, :]), dims = :file), dims = :file)
                
                if !(any(valid1))
                    continue
                end

                vdate, vheight  = Altim.validrange(valid1)
                
                for date in ddate[vdate]

                    for height in dheight[vheight]
                        var0 = dh_all[mission][:, At(geotile), At(date), At(height)];
                        valid2 = .!isnan.(var0);

                        if any(valid2)
                            dh_all_std[mission][At(geotile), At(date), At(height)] = max(std(var0[valid2]),minimum_error)
                        end
                    end
                end
            end
        end

        # save the error so that it can be used in the synthesis of individual model runs
        save(outfile, Dict("dh_hyps_error" => dh_all_std, "files_included" => path2runs)); 
    end
    return outfile
end

## synthesize missions weighting by synthesis_error
# this take about 2 hours for 252 run files
path2geotile_synthesis_error = geotile_synthesis_error(; path2runs, outfile="/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2", force_remake=false)

function geotile_synthesize_runs(;
    path2runs,
    path2geotile_synthesis_error,
    missions2include = ["hugonnet", "gedi", "icesat", "icesat2"],
    showplots = false,
    force_remake = false
    )

    dh_err = load(path2geotile_synthesis_error, "dh_hyps_error")
    missions = missions2include
    dgeotile = dims(dh_err[missions2include[1]], :geotile)

    # convert error to weights
    w = copy(dh_err)
    for mission in missions
        #w0 = w[mission];
        w[mission] = 1 ./ (dh_err[mission].^2)

        # mask GEDI outside of observational latitude limits
        if mission == "gedi"
            lat_min = getindex.(Altim.geotile_extent.(dgeotile), :min_y)
            lat_max = getindex.(Altim.geotile_extent.(dgeotile), :max_y)
            exclude_gedi = (lat_max .> 51.6) .| (lat_min .< -51.6)
            w[mission][exclude_gedi,:,:] .= 0
        end

        w[mission][isnan.(w[mission])] .= 0;
    end
  
    Threads.@threads for binned_aligned_file in path2runs
    #binned_aligned_file = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_b1km_dh_best_cc_meanmadnorm3_v01_filled_ac_p1_aligned.jld2"    

        if showplots
            file_parts = splitpath(binned_aligned_file)
            binned_folder = joinpath(file_parts[1:findfirst(file_parts.=="2deg")])

            fig_folder = joinpath(binned_folder, "figures")
            figure_suffix = replace(file_parts[end], ".jld2" => "")
            binned_synthesized_file = replace(binned_aligned_file, "aligned.jld2" => "synthesized.jld2")
            params = Altim.binned_filled_fileparts(binned_aligned_file)
            geotiles = Altim.geotiles_mask_hyps(params.surface_mask, geotile_width)
        end
        
        binned_synthesized_file = replace(binned_aligned_file, "aligned.jld2" => "synthesized.jld2")

        if !(isfile(binned_synthesized_file)) || force_remake
            t1 = time()
            dh = FileIO.load(binned_aligned_file, "dh_hyps");

            #heatmap(Altim.dh2dv(dh["hugonnet"], geotiles["glacier_b1km"])[510:511,:])
            #heatmap!(Altim.dh2dv(dh["icesat2"], geotiles["glacier_b1km"])[510:511,:])
            #heatmap!(Altim.dh2dv(dh["gedi"], geotiles["glacier_b1km"])[510:511,:])

            if showplots
              p = Altim.plot_height_time(dh; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="final", fig_folder, figure_suffix, mask=params.surface_mask, showplots)
            end
            
            dh_synth = fill(0., dims(dh[missions[1]]))
            w_synth = copy(dh_synth)
            
            for mission in missions
                w0 = copy(w[mission])

                var0 = dh[mission];
                nanidx = isnan.(var0);
                var0[nanidx] .= 0
                w0[nanidx] .= 0

                # weighting needs to be non-zero for valid 
                w0[.!nanidx .& (w0 .== 0)] .= 0.001

                dh_synth += (var0 .* w0)

                w_synth += w0
            end

            dh_synth = dh_synth ./ w_synth            
            dh_synth_err = sqrt.(1 ./ w_synth)
            dh_synth_err[isnan.(dh_synth)] .= NaN

            #if plot_dh_as_function_of_time_and_elevation
            if showplots
                # load geotiles
                p = Altim.plot_height_time(dh_synth; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="raw", fig_folder, figure_suffix, mask=params.surface_mask, mission = "synthesis", showplots)              
            end

            save(binned_synthesized_file, Dict("dh_hyps" => dh_synth, "dh_hyps_err" => dh_synth_err))
            println("$binned_aligned_file synthesized: $(round(Int,time() -t1))s")
        end
    end
end

# this takes about 2 hours for 252 run files
geotile_synthesize_runs(;
    path2runs,
    path2geotile_synthesis_error,
    missions2include=["hugonnet" "gedi" "icesat" "icesat2"],
    force_remake=false
)

## find optimal fit to GEMB data

#### add fac correction here
path2runs = replace.(path2runs, "aligned.jld2" => "synthesized.jld2")

geotile_width = 2;
force_remake = false

dh = FileIO.load(path2runs[1], "dh_hyps")

# load and align geotile dataframe with DimArrays
dgeotile = dims(dh, :geotile)

run_parameters_all = Altim.binned_filled_fileparts.(path2runs)
surface_masks = unique(getindex.(run_parameters_all, :surface_mask))

geotiles = Dict()
for surface_mask in surface_masks
    geotiles0 = Altim.geotiles_mask_hyps(surface_mask, geotile_width)
    geotiles0, reg = Altim.geotiles_mutually_exclusive_rgi!(copy(geotiles0))
    gt_ind = [findfirst(geotiles0.id .== g0) for g0 in collect(dgeotile)]
    rename!(geotiles0,"$(surface_mask)_area_km2" => "area_km2")
    geotiles[surface_mask] = geotiles0[gt_ind, :]
end

# data was not filled properly ... so apply here until geotile_synthesis_error is rerun
# this is a helper function
function fill_nans_to_make_valid_rectangle_old(dh)
    dgeotile = dims(dh, :geotile)
    for geotile in dgeotile
        #geotile = first(dgeotile)
        #geotile = "lat[-28-26]lon[-070-068]"
        dh0 = @view(dh[At(geotile),:,:])
        if all(isnan.(dh0))
            continue
        end
        rrange, crange = Altim.validrange(.!isnan.(dh0))
        if any(isnan.(dh0[rrange, crange]))
            dh1 = @view dh0[rrange, crange]
            valid = .!isnan.(dh1)
            for i in 1:length(dims(dh1, 1))
                if any(valid[i, :])
                    f = findfirst(vec(valid[i, :]))
                    l = findlast(vec(valid[i, :]))
                    dh1[i, 1:f] .= dh1[i, f]
                    dh1[i, l:end] .= dh1[i, l]
                end
            end
        end
    end
    return dh
end

# load gemb data
filename_gemb_geotile_filled_dv = replace(filename_gemb_combined, ".jld2" => "_geotile_filled_dv.jld2")
dv_gemb = Dict()
smb, fac = load(filename_gemb_geotile_filled_dv, ("smb", "fac"))

for binned_synthesized_file in path2runs
#binned_synthesized_file = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_b1km_dh_best_cc_meanmadnorm3_v01_filled_ac_p1_synthesized.jld2"

#binned_synthesized_file = path2runs[1]
#println(binned_synthesized_file)

    synthesized_gemb_fit = replace(binned_synthesized_file,  ".jld2" => "_gemb_fit.jld2")

    if !(isfile(synthesized_gemb_fit)) || force_remake
        t1 = time()

        run_parameters = Altim.binned_filled_fileparts(synthesized_gemb_fit)
        surface_mask = run_parameters.surface_mask;
        dh = FileIO.load(binned_synthesized_file, "dh_hyps")
 

        # for the Antarctic we only want to fit to ICESat/-2
        #aligned = replace(binned_synthesized, "synthesized.jld2" => "aligned.jld2")
        #dh0 = FileIO.load(binned_synthesized_file, "dh_hyps")

        dh = fill_nans_to_make_valid_rectangle_old(dh)

        # convert elevation change to volume change
        dv_altim = Altim.dh2dv(dh, geotiles[surface_mask]);

        #=
        # ERROR IS IN THE SYNTHESIS !!!!
        heatmap(Altim.dh2dv(dh, geotiles[surface_mask])[510:511,:])


        heatmap(dv_altim[510:511,:])

        binned_aligned_file = replace(binned_synthesized_file,  "synthesized" => "aligned")
        dhX = load(binned_aligned_file, "dh_hyps")
        heatmap(Altim.dh2dv(dhX["icesat"], geotiles[surface_mask])[510:511,:])
        heatmap(Altim.dh2dv(dhX["icesat2"], geotiles[surface_mask])[510:511,:])
        heatmap(Altim.dh2dv(dhX["hugonnet"], geotiles[surface_mask])[510:511,:])
        heatmap(Altim.dh2dv(dhX["gedi"], geotiles[surface_mask])[510:511,:])

        binned_filled_file = replace(binned_aligned_file,  "_aligned"=>"")
        dhX = load(binned_filled_file, "dh_hyps")
        heatmap(Altim.dh2dv(dhX["icesat"], geotiles[surface_mask]))
        heatmap(Altim.dh2dv(dhX["icesat2"], geotiles[surface_mask]))
        heatmap(Altim.dh2dv(dhX["hugonnet"], geotiles[surface_mask]))
        heatmap(Altim.dh2dv(dhX["gedi"], geotiles[surface_mask]))
=#
        if any(all(isnan.(dv_altim), dims = :date))
            dgeotile = dims(dv_altim, :geotile)
            all_nan_geotiles = dgeotile[findall(dropdims(all(isnan.(dv_altim), dims = :date), dims = :date))]
            printstyled("Geotiles that only contain NaNs:\n"; color = :red)
            println("$(collect(all_nan_geotiles)))")
            printstyled("Possible sources of error include:\n"; color = :red)
            printstyled("  [1] synthesis_error_file was run on a subset of files that do not include the full error range of the data, to fix this you need to delete the exisiting error file and any synthesis files that were create :\n"; color = :red)
            error("geotile volume change contains all NaNs")
        end

        # find optimal fit to gemb data
        # df = Altim.gemb_bestfit(dv_altim, smb, fac, discharge, geotiles[surface_mask])

        # there are issues with calibrating the smb model to individual geotiles as dischage 
        # is assigned at a point but can cross multiple geotiles. To tackle this issue we 
        # calibrate the model at a regional scale.

        # TODO: If we could identify groupings of geotiles that enclosed groupings of 
        # glaicers then smb calibration could be done on a more granular scale but this 
        # would be hard to do and might look very similar to what we have now. Another 
        # option is to calculate flux across geotiles but this would be a v2 project. 
        try
            df = Altim.gemb_bestfit_regional(dv_altim, smb, fac, discharge, geotiles[surface_mask])

            df[!,:area_km2] = sum.(geotiles[surface_mask].area_km2)
            df.mad = df.mad ./ df.area_km2
            rename!(df, "mad"=>"mad_m")

            save(synthesized_gemb_fit, Dict("gemb_fit" => df))
            println("$binned_synthesized_file optimal GEMB fit found: $(round(Int,time() -t1))s")
        catch
              printstyled("$binned_synthesized_file error in optimal GEMB fit \n"; color = :red)
        end

    end
end


function geotile_zonal_area_hyps(ras, ras_range, zone_geom, geotile_ids; persistent_attribute = :RGIId)
    
    df = DataFrame()

    for geotile_id0 in geotile_ids
    #geotile = first(geotile_ids)
        t1 = time()
        bounding_polygon = Altim.extent2rectangle(Altim.GeoTiles.extent(geotile_id0))
        index = GeometryOps.intersects.(zone_geom[!, :geometry], Ref(bounding_polygon))
        if !any(index)
            println("skipping $(geotile_id0): no intersecting polygons")
            continue
        end

        zone_geom0 = DataFrame(zone_geom[index, :])
        zone_geom0[!, :geotile] .= geotile_id0

        # do a double crop as cropping the polygons themselves is just way to complicated right now
        ras0 = Rasters.crop(ras, to = zone_geom0)
        ras0 = Rasters.crop(ras0, to = bounding_polygon)

        ras0 = ras0[:,:]# having issues going from lazy to inmem ... not logical but this works

        rs = RasterStack(ras0, Rasters.cellarea(ras0))

        area_m2 = Altim.mapzonal(Base.Fix2(geotile_zonal_area_hyps, ras_range), identity, rs; of=zone_geom0[!, :geometry])
        zone_geom0[!, :area_km2] = area_m2 * (1E-3)^2

        println("$(geotile_id0) zonal area hyps done: $(round(Int,time() -t1))s")
        
        append!(df, zone_geom0[:, [:geotile, persistent_attribute, :area_km2]]; promote=true)
    end
    return df
end

# for a single geotile
# Define binning functions
function geotile_zonal_area_hyps(x, y, bins::StepRange)
    n = length(bins) - 1
    bin_index = @. ceil(Int64, ($n) * (x - $first(bins)) / ($last(bins) - $first(bins)))

    binned = zeros(eltype(y), n)
    for (idx, bin_idx) in enumerate(bin_index)
        if bin_idx > 0 && bin_idx < (n + 1)
            binned[bin_idx] += y[idx]
        end
    end

    return binned
end

function geotile_zonal_area_hyps(nt, x_bin_edges)
    binned = geotile_zonal_area_hyps(getindex.(nt, 1), getindex.(nt, 2), x_bin_edges)
    return binned
end

#=
function geotile_zonal_area_hyps(ras, ras_range, zone_geom, geotile_id::String)

    bounding_polygon = Altim.extent2rectangle(Altim.GeoTiles.extent(geotile_id))
    ras0 = Rasters.crop(ras, to = bounding_polygon)

    index = GeometryOps.intersects.(zone_geom[:, :geom], Ref(bounding_polygon))

    if !any(index)
        return nothing
    end

    zone_geom0 = DataFrame(zone_geom[index, :])
    zone_geom0[!, :local_id] = Int16.(1:nrow(zone_geom0))
    zone_geom0[!, :geotile] .= geotile_id

     # find local UTM/PS zone 
    epsg = Altim.utm_epsg(GeometryOps.centroid(bounding_polygon)...)
    target_crs = GeoFormatTypes.EPSG(epsg)
    source_crs = crs(ras0)

    # reproject elevation data
    utm_res = 50;
    ras0 = Rasters.resample(ras0; res = utm_res, crs = target_crs, method = :near)

    # reproject zonal_geom0
    rename!(zone_geom0, "geom" => "geometry")
    zone_geom0[!,:geom] = GeometryOps.reproject(zone_geom0[:,:geometry]; source_crs, target_crs)
    
    # bin data 
    bin_index = ceil.(Int16, (length(ras_range)-1) .* (ras0 .- first(ras_range)) ./ (last(ras_range) - first(ras_range)))
    
    glacier_id = rasterize!(last, zeros(Int16,dims(ras0)), zone_geom0; fill = :local_id)

    # index into raterized geometries that call within bin range
    index = (glacier_id .> 0) .& (bin_index .> 0) .& (bin_index .<= length(ras_range));

    df = DataFrame(id=glacier_id[index], bin=bin_index[index])
    gdf = DataFrames.groupby(df, [:id,:bin])

    count = zeros(nrow(zone_geom0), length(ras_range)-1)
    for k =  eachindex(gdf)
        ind = CartesianIndex(collect(k)...)
        count[ind] = nrow(gdf[k])
    end

    area_km2 =  count .* (utm_res/1000)^2
    zone_geom0[!, :area_km2] = eachrow(area_km2)
    return zone_geom0
end

#@time begin
# compute glacier wide statisitics
geomfile_rgi6 = joinpath(paths.data_dir, "GlacierOutlines/rgi60/rgi60_Global.gpkg")
    glacier_geotile_hyps_fn = replace(geomfile_rgi6, ".gpkg" => "geotile_hyps.jld2")
glacier_geom = GeoDataFrames.read(geomfile_rgi6)

h = Raster(paths.cop30_v2, lazy=true)
height_range, height_center = Altim.project_height_bins()
geotiles = Altim.geotiles_w_mask(geotile_width)
geotiles = geotiles[(geotiles.glacier_frac.>0.0), :]

if !isfile(glacier_geotile_hyps_fn)
    df = geotile_zonal_area_hyps(h, height_range, glacier_geom, geotiles.id; persistent_attribute=:RGIId)

    glacier_geom = innerjoin(glacier_geom, df, on = :RGIId)

    # this can be removed next interation
    md = DataFrames.metadata(glacier_geomX)
    for k in keys(md)
        DataFrames.metadata!(glacier_geom, k, md[k])
    end
    #
    glacier_geom[!, :geom] = GeometryOps.tuples(glacier_geom.geom)

    # save as a jld2 as there are issues with saving and loading a vector (area_km2) to a geopkg right now
    FileIO.save(glacier_geotile_hyps_fn, Dict("glacier" => glacier_geom))
else
    # this is slow... but hey, it works
    glacier = FileIO.load(glacier_geotile_hyps_fn, "glacier")
end
=#

# this is slow... but hey, it works
geomfile_rgi = geomfile_rgi6
glacier_geotile_hyps_fn = replace(geomfile_rgi, ".gpkg" => "geotile_hyps.jld2")
glaciers0 = FileIO.load(glacier_geotile_hyps_fn, "glaciers")

# once I have the hypsometry of each glacier for all geotiles I think simply loop over geotiles and mutiply hysometry times variable and sum over elevation range
filename_gemb_geotile_filled = replace(filename_gemb_combined, ".jld2" => "_geotile_filled.jld2")
gemb = load(filename_gemb_geotile_filled)

#filename_gemb_geotile_filled = replace(filename_gemb_combined, ".jld2" => "_geotile.jld2")
#gemb = load(filename_gemb_geotile_filled)


# loop for each model run
if false
    # Downscale geotile-level data to individual glaciers by:
    # 1. Loading the GEMB fit parameters and height change data for each geotile
    # 2. Converting geotile height changes to per-glacier values
    # 3. Scaling GEMB variables (SMB, FAC, runoff) using fit parameters and converting to per-glacier values
    # Takes ~8 min per iteration
    println("start downscaling from geotile to glacier, each iteration takes about 8 min so be patient")
    for binned_synthesized_file in path2runs
    #binned_synthesized_file = path2runs[1]

        synthesized_gemb_fit = replace(binned_synthesized_file, ".jld2" => "_gemb_fit.jld2")
        perglacier_synthesized_file = replace(binned_synthesized_file, ".jld2" => "_perglacier.jld2")
        
        if !isfile(perglacier_synthesized_file)
            t1 = time()
            gemb_fit = load(synthesized_gemb_fit, "gemb_fit")
            dh = load(binned_synthesized_file, "dh_hyps")
            glaciers = copy(glaciers0);

            glaciers = Altim.geotile2glacier!(glaciers, dh; varname = "dh")
            for varname = ["smb", "fac", "runoff"]
                glaciers = Altim.geotile2glacier!(glaciers, Altim.gembscale(gemb[varname], gemb_fit); varname)
            end
            
            FileIO.save(perglacier_synthesized_file, Dict("glaciers" => glaciers))

            fn = splitpath(perglacier_synthesized_file)
            println("downscaling from geotile to glacier complete $(round(Int,time() -t1))s: $(fn[end])")
        end
    end
end

end


binned_synthesized_file  = "/mnt/bylot-r3/data/binned/2deg/glacier_dh_best_meanmadnorm3_v01_filled_ac_p1_synthesized.jld2"
synthesized_gemb_fit = replace(binned_synthesized_file, ".jld2" => "_gemb_fit.jld2")
perglacier_synthesized_file = replace(binned_synthesized_file, ".jld2" => "_perglacier.jld2")
glaciers = FileIO.load(perglacier_synthesized_file, "glaciers")

# remove glaciers without area
zero_area_index  = sum.(glaciers.area_km2) .== 0
glaciers = glaciers[.!zero_area_index,:]

# remove glaciers where the elevation profile fell outside of valid range (only 12 glaciers)
glaciers = glaciers[.!map(x -> all(isnan.(x)), glaciers.dh),:]

# read in regions
rgi_regions = GeoDataFrames.read(paths.rgi6_regions_shp)

# read in GRACE data
grace = Altim.read_grace_rgi(; datadir=setpaths()[:grace_rgi])

# select a region
begin
# TODO: rgi 13/14/15 gemb and obs do not match.

if false
    geotile = "lat[+28+30]lon[+082+084]"# "lat[+56+58]lon[-134-132]"
    extent = geotiles["glacier"][findfirst(geotiles["glacier"].id .== geotile), :extent]
    region_geom = Altim.extent2rectangle(extent)
    in_region = GeometryOps.within.(tuple.(glaciers.CenLon, glaciers.CenLat), Ref(region_geom))
    in_region_dischage = GeometryOps.within.(tuple.(discharge.longitude, discharge.latitude), Ref(region_geom))
    
    title = geotile
else
    #rgi = 13
    rgi = 19
    in_region = falses(nrow(glaciers))
    in_region_dischage = falses(nrow(discharge))
    for region_geom = rgi_regions[rgi_regions.RGI_CODE.==rgi,:geometry]
        in_region .|= GeometryOps.within.(tuple.(glaciers.CenLon, glaciers.CenLat), Ref(region_geom))
        in_region_dischage .|= GeometryOps.within.(tuple.(discharge.longitude, discharge.latitude), Ref(region_geom))
    end
    title = Altim.rgi2label["rgi$rgi"]
end

# export for plotting
if false
    var_name = "runoff"
    v0 = copy(glaciers[in_region, :])
    start_date = DateTime(2000,04,1)
    end_date= DateTime(2023,04,1)

    for var_name in ["runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze", "dh"]
        sindex = findfirst(dims(v0[1, var_name], :date) .> start_date)
        eindex = findlast(dims(v0[1, var_name], :date) .< end_date)
        v1 = getindex.(v0[:, var_name],  eindex ) .- getindex.(v0[:, var_name],  sindex)
        v1 = v1.*sum.(glaciers[in_region, :area_km2])/1000
        v0[!, var_name] = v1
    end

    rename!(v0, "geom" => "geometry")
    #v0[!, :geometry] = tuple.(v0.CenLon,v0.CenLat)
    GeoDataFrames.write("glacier_2000_2023_anomaly.gpkg", v0[:, ["geometry", "RGIId", "runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze", "dh"]])
end

# calculate retional volume sums all sampled to the dh data
reg = Dict()
ddate = dims(glaciers[1, :dh], :date)
decyear = Altim.decimalyear.(ddate)
Δdecyear = decyear .- decyear[1]

for var_name in ["dh", "runoff", "fac", "smb"]; #, "rain", "acc", "melt", "ec", "refreeze", "dh"]
#var_name = "dh"
    
    v0 = glaciers[in_region , var_name]

    decyear0 = Altim.decimalyear.(dims(v0[1], :date))
    area_km2 = sum.(glaciers[in_region, :area_km2])
    v0 = v0 ./ 1000 .* area_km2
    v0 = vec(sum(reduce(hcat, v0), dims=2))
   
    notnan = .!isnan.(v0)
    decyear0 = decyear0[notnan]
    v0 = v0[notnan]
    index = (decyear .>= minimum(decyear0)) .& (decyear .<= maximum(decyear0))

    interp = DataInterpolations.LinearInterpolation(v0, decyear0)
    out = fill(NaN, ddate)
    out[index] = interp(decyear[index])

    reg[var_name] = out
end

reg["dischage"] = sum(discharge[in_region_dischage, :discharge_gtyr])

# add grace data
begin
    var_name = "grace"

    grace_date = Dim{:date}(Altim.datenum2date.(vec(grace["rgi$(rgi)"]["dM_gt_mdl_fill_date"])))
    grace_data = DimArray(vec(grace["rgi$(rgi)"]["dM_gt_mdl_fill"]), grace_date)

    grace_notnan = .!isnan.(grace_data)
    grace_decyear = Altim.decimalyear.(grace_date[grace_notnan])
    index = (decyear .>= minimum(grace_decyear)) .& (decyear .<= maximum(grace_decyear))

    grace_interp = fill(NaN, ddate)
    interp = DataInterpolations.LinearInterpolation(grace_data[grace_notnan], Altim.decimalyear.(grace_date[grace_notnan]))
    grace_interp[index] = interp(decyear[index])
    reg[var_name] = grace_interp
end



begin
    volume2mass = Altim.δice / 1000

    if false
        A = reg["dh"] .- reg["fac"]
        A_label = "this study"
        B = reg["grace"]
        B_label = "GRACE/-FO"
        y_label = "mass anomaly [Gt]"
    else
        A = reg["dh"] 
        A_label = "observed"
        B = reg["smb"] ./ volume2mass .+ reg["fac"] .- (reg["dischage"] .* Δdecyear) ./ volume2mass
        B_label = "modeled"
        y_label = "volume anomaly [km³]"
    end

    # align B to A
    index = .!isnan.(A) .& .!isnan.(B)
    B = B .- mean(B[index].-A[index])

    p = plot(A; title, label=A_label)
    plot!(B;title, label=B_label)
    ylabel!(p,y_label)
    xlims!((DateTime(2000, 1, 1), DateTime(2024, 1, 1)))
    display(p)
end

end