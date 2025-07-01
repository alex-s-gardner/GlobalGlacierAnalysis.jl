"""
    geotile_synthesize(path2runs; error_file, mission_error, update_missions, force_remake_before, missions2include)

Synthesize elevation change data across multiple satellite missions for improved spatial coverage and temporal consistency.

# Arguments
- `path2runs`: Vector of file paths to binned altimetry data files
- `error_file`: Path to save geotile synthesis error assessment (default: "/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2")
- `mission_error`: Mission-specific error estimates for data quality assessment
- `update_missions`: Force update of mission-specific geotile data (default: nothing)
- `force_remake_before`: Date before which to force recreation of existing files (default: nothing)
- `missions2include`: Array of mission names to include in synthesis (default: ["hugonnet", "gedi", "icesat", "icesat2"])

# Description
This function performs a two-step synthesis process:
1. Calculates geotile synthesis error to assess data quality and consistency across missions
2. Combines elevation change data from different satellite missions into a single dataset

The synthesis improves spatial coverage and temporal consistency by leveraging multiple altimetry datasets.
"""
function geotile_synthesize(path2runs;
    error_file="/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2",
    mission_error = fill(0.0, Dim{:mission}(["hugonnet", "gedi", "icesat", "icesat2"])),
    update_missions=nothing,
    single_geotile_test=nothing,
    missions2include=["hugonnet" "gedi" "icesat" "icesat2"],
    force_remake_before=nothing,
)
    # Calculate geotile synthesis error to assess data quality and consistency across missions
    dh_hyps_error, _ = geotile_synthesis_error(;
        path2runs,
        outfile=error_file,
        mission_error,
        update_missions,
        force_remake_before,
        single_geotile_test,
    )


    # Synthesize geotile runs across multiple missions
    # This combines elevation change data from different satellite missions into a single dataset
    # for improved spatial coverage and temporal consistency
    geotile_synthesize_runs(;
        path2runs,
        dh_err=dh_hyps_error,
        missions2include,
        force_remake_before = Dates.unix2datetime(mtime(error_file)), # synthesis files need to be update if synthesis_error was updated
    )
end

"""
    error_over_land(binned_file="/mnt/bylot-r3/data/binned_unfiltered/2deg/land_dh_best_cc_nmad5_v01.jld2")

Calculate mission-specific error estimates over land surfaces using median absolute deviation (MAD).

# Arguments
- `binned_file`: Path to binned elevation change data file containing land surface measurements

# Returns
- `land_error`: DimArray containing MAD-based error estimates for each satellite mission over land surfaces

# Description
Loads binned elevation change data and computes mission-specific error estimates using the median 
absolute deviation (MAD) statistic. This provides robust error estimates for each satellite mission
when measuring elevation changes over stable land surfaces, which can be used to assess data quality
and inform synthesis weighting schemes.
"""
function binned_mad_mission(binned_file="/mnt/bylot-r3/data/binned_unfiltered/2deg/land_dh_best_cc_nmad5_v01.jld2")

    f = load(binned_file, "dh_hyps")

    dmission = Dim{:mission}(collect(keys(f)))
    land_error = fill(NaN, dmission)

    for mission in dmission
        valid = .!isnan.(f[mission])
        land_error[At(mission)] = mad(f[mission][valid])
    end

    return land_error
end

"""
    geotile_synthesis_error(;
        path2runs,
        outfile="/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2",
        mission_error,
        update_missions=nothing,
        single_geotile_test=nothing,
        force_remake_before=nothing
    )

Calculate standard deviation across model runs as an error metric for each satellite mission.

# Arguments
- `path2runs`: Vector of paths to input data files containing binned elevation change data
- `outfile`: Path for saving computed error metrics and file tracking information
- `mission_error`: Floor error values for each mission to ensure minimum error estimates
- `update_missions`: Optional list of specific missions to update instead of processing all missions
- `single_geotile_test`: Single geotile identifier for testing purposes (output not saved to file)
- `force_remake_before`: DateTime threshold for forcing recalculation of existing results

# Returns
- `dh_hyps_error`: Dictionary containing error estimates for each mission across geotiles, dates, and heights
- `files_included`: Vector of file paths that were processed

# Description
Computes error metrics by calculating standard deviation across multiple model runs for each satellite mission.
The function handles incremental updates, file caching, and single geotile testing scenarios.
"""
function geotile_synthesis_error(;
    path2runs,
    outfile="/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2",
    mission_error = fill(0.0, Dim{:mission}(["hugonnet", "gedi", "icesat", "icesat2"])),
    update_missions=nothing,
    single_geotile_test=nothing, #geotiles_golden_test[1], #geotiles_golden_test[2],
    geotile_width = 2,
    force_remake_before=nothing,
)

    if .!isnothing(single_geotile_test)
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)], OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!"
    end

    if isfile(outfile) && isnothing(force_remake_before) && isnothing(update_missions) && isnothing(single_geotile_test)

        if !issetequal(path2runs, load(outfile, "files_included"))
            @warn "Updating all missions in geotile_synthesis_error as input file list does not match exisiting"
        else
            printstyled("    -> Loading $(outfile) as it already exists \n"; color=:light_green)
            (dh_hyps_error, files_included) = load(outfile,("dh_hyps_error", "files_included" ))
            return (dh_hyps_error, files_included)
        end

    elseif isfile(outfile) && !isnothing(force_remake_before) && (Dates.unix2datetime(mtime(outfile)) >= force_remake_before) && isnothing(single_geotile_test)
         if !issetequal(path2runs, load(outfile, "files_included"))
            @warn "Updating all missions in geotile_synthesis_error as input file list does not match exisiting"
         else
            printstyled("    -> Loading $(outfile) as it was created after the force_remake_before date: $force_remake_before \n"; color=:light_green)
            (dh_hyps_error, files_included) = load(outfile,("dh_hyps_error", "files_included" ))
            return (dh_hyps_error, files_included)
         end
    elseif isfile(outfile) && !isnothing(force_remake_before) && (Dates.unix2datetime(mtime(outfile)) < force_remake_before) && isnothing(single_geotile_test)
        printstyled("    -> Recomputing $(outfile) as it was created before the force_remake_before date: $force_remake_before \n"; color=:light_green)
    elseif isnothing(single_geotile_test)
        printstyled("    -> Computing $(outfile) \n"; color=:light_green)
    end

    # load example file to get dimension and mission info
    dh = FileIO.load(path2runs[1], "dh_hyps")
    missions = collect(keys(dh))
    ddate = dims(dh[missions[1]], :date)
    dheight = dims(dh[missions[1]], :height)
    dgeotile = dims(dh[missions[1]], :geotile)
    dfile = Dim{:file}(path2runs)

    missions2update = missions

    if isfile(outfile) && !isnothing(update_missions) && issetequal(path2runs, load(outfile, "files_included")) && isnothing(single_geotile_test)
        missions2update = update_missions
        printstyled("\n    -> Updating selected missions in geotile_synthesis_error: $(update_missions)")
    elseif isfile(outfile) && !isnothing(update_missions) && !issetequal(path2runs, load(outfile, "files_included"))
        @warn "Updating all missions in geotile_synthesis_error as input file list does not match exisiting: $(outfile)"
    end

    if .!isnothing(single_geotile_test)
        dgeotile = Dim{:geotile}([single_geotile_test])
    end

    # initialize output
    dh_all_std = Dict()
    if !isnothing(update_missions)
        dh_all_std = load(outfile, "dh_hyps_error")
    else
        dh_all_std = Dict()
        for mission in missions
            dh_all_std[mission] = fill(NaN, (dgeotile, ddate, dheight)) 
        end
    end

    # loop over missions to reduce memory usage
    for mission in missions2update

        printstyled("    -> loading all $(mission) data... takes ~1 min \n"; color=:light_grey)
        
        dh_all = fill(NaN, (dfile, dgeotile, ddate, dheight))

        # this takes about 8 min for length(files) = 96 # this blows up memory usage
        Threads.@threads for filepath in path2runs
            dh = FileIO.load(filepath, "dh_hyps")[mission][geotile=At(collect(dgeotile))]
            dh_all[file = At(filepath)] = dh
        end

        @showprogress dt = 1 desc = "Calculating standard deviation (error) across runs for $(mission) ..." Threads.@threads for geotile in dgeotile

            valid1 = dropdims(any(.!isnan.(dh_all[geotile=At(geotile)]), dims=:file), dims=:file)

            if !(any(valid1))
                continue
            end

            vdate, vheight = validrange(valid1)

            # for some geotiles and some missions there are no valid data so this is needed over using dims
            for date in ddate[vdate]

                for height in dheight[vheight]
                    var0 = dh_all[:, At(geotile), At(date), At(height)]
                    valid2 = .!isnan.(var0)

                    if sum(valid2) > 1 # std returns NaN if there is only one valid data point
                        std0 = std(var0[valid2]) + mission_error[At(mission)]
                        dh_all_std[mission][At(geotile), At(date), At(height)] = std0 
                    end
                end
            end
        end
    end

    dh_hyps_error = dh_all_std
    files_included = path2runs

    @warn "testing"
    return dh_hyps_error, files_included

    if isnothing(single_geotile_test)
        # save the error so that it can be used in the synthesis of individual model runs
        save(outfile, Dict("dh_hyps_error" => dh_hyps_error, "files_included" => files_included))
    else
        params = binned_filled_fileparts.(files_included)
        surface_masks = unique(getindex.(params, :surface_mask))
        surface_mask = surface_masks[1]

        area_km2 = _geotile_area_km2(surface_mask, geotile_width)
       
        outfile_parts = splitpath(outfile)
        plot_save_path_prefix = joinpath(pathlocal[:figures], replace(outfile_parts[end], ".jld2" => ""))

        f = plot_elevation_time_multimission_geotiles(
            dh_hyps_error,
            geotiles2plot=[single_geotile_test];
            area_km2,
            colorrange=(-10, 10),
            label=nothing,
            colorbar_label="height anomaly error",
            hypsometry=true,
            area_averaged=true,
            plots_show=true,
            plots_save=true,
            plot_save_path_prefix,
            plot_save_format=".png"
        )        
    end

    return (dh_hyps_error, files_included)    
end

"""
    geotile_synthesize_runs(;
        path2runs,
        dh_err,
        missions2include=["hugonnet", "gedi", "icesat", "icesat2"],
        geotile2plot=nothing,
        showplots=false,
        single_geotile_test=nothing,
        force_remake_before=nothing
    )

Synthesize elevation change data from multiple satellite missions into a single dataset.

# Arguments
- `path2runs`: Paths to input data files
- `dh_err`: Mission error estimates for each mission
- `missions2include`: Mission names to include in synthesis
- `geotile2plot`: Specific geotile to plot (optional)
- `showplots`: Whether to generate diagnostic plots
- `single_geotile_test`: Single geotile for testing (optional)
- `force_remake_before`: Recompute if output file exists but was created before this date

# Returns
Path to synthesized output file with weighted average elevation change data and uncertainty estimates
"""
function geotile_synthesize_runs(;
    path2runs,
    dh_err,
    missions2include=["hugonnet", "gedi", "icesat", "icesat2"],
    geotiles2plot = nothing,
    single_geotile_test=nothing, #geotiles_golden_test[1], #geotiles_golden_test[2],
    dh_override = nothing, # used to override the dh_hyps data with a different set of data
    geotile_width = 2,
    force_remake_before=nothing,
)

    if !isnothing(single_geotile_test)
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)], OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!\n"
        geotiles2plot = [single_geotile_test]
    end

    if !isnothing(dh_override)
        @warn "!!!!!!!!!!!!!! USING dh_override, NOT READING FROM FILE, OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!\n"
        path2runs = ["no filenames needed"]
    end

    missions = missions2include
    dgeotile = dims(dh_err[missions2include[1]], :geotile)

    if !isnothing(single_geotile_test)
        dgeotile = Dim{:geotile}([single_geotile_test])
        dheight = dims(dh_err[missions2include[1]], :height)
        ddate = dims(dh_err[missions2include[1]], :date)
        for mission in missions
            foo = zeros(dgeotile, ddate, dheight)
            foo[1,:,:] =  dh_err[mission][geotile=At(single_geotile_test)]
            dh_err[mission] = foo
        end
    end

    # convert error to weights
    w = copy(dh_err)

    for mission in missions
        #w0 = w[mission];
        w[mission] = 1 ./ (dh_err[mission] .^ 2)

        ## mask GEDI outside of observational latitude limits
        if mission == "gedi"
            lat_min = minimum.(getindex.(geotile_extent.(dgeotile), :Y))
            lat_max = maximum.(getindex.(geotile_extent.(dgeotile), :Y))
            exclude_gedi = (lat_max .> 51.6) .| (lat_min .< -51.6)
            w[mission][exclude_gedi, :, :] .= 0
        end

        w[mission][isnan.(w[mission])] .= 0
    end

    # can not use Threads.@threads here as it will not return values when dh_override != nothing
    @showprogress dt = 1 desc = "Synthesizing runs ..." for binned_aligned_file in path2runs
        #binned_aligned_file = path2runs[1]   

        binned_synthesized_file = replace(binned_aligned_file, "aligned.jld2" => "synthesized.jld2")

        if isfile(binned_synthesized_file) && (isnothing(force_remake_before) || Dates.unix2datetime(mtime(binned_synthesized_file)) > force_remake_before) && isnothing(single_geotile_test) && isnothing(dh_override)
            printstyled("    -> Skipping $(binned_synthesized_file) because it was created after the latest synthesis_error: $force_remake_before\n"; color=:light_green)
            continue
        else
            t1 = time()
            

            if !isnothing(dh_override)
                dh = dh_override

                w00 = deepcopy(w)
                dgeotile = dims(dh[missions[1]], :geotile)

                for mission in missions
                    w00[mission] = w00[mission][geotile=At(collect(dgeotile))]
                end
            else
                dh = FileIO.load(binned_aligned_file, "dh_hyps")
                w00 = deepcopy(w)
            end

            if !isnothing(single_geotile_test)
                for mission in missions
                    foo = zeros(dgeotile, ddate, dheight)
                    foo[1,:,:] = dh[mission][geotile=At(single_geotile_test)]
                    dh[mission] = foo
                end
            end

            dh_synth = fill(0.0, dims(dh[missions[1]]))
            w_synth = copy(dh_synth)

            for mission in missions
                w0 = w00[mission]
                var0 = dh[mission]

                nanidx = isnan.(var0)
                var0[nanidx] .= 0
                w0[nanidx] .= 0

                # weighting needs to be non-zero for valid 
                w0[.!nanidx.&(w0.==0)] .= 0.001

                dh_synth += (var0 .* w0)

                w_synth += w0
            end

            dh_synth = dh_synth ./ w_synth
            dh_synth_err = sqrt.(1 ./ w_synth)
            dh_synth_err[isnan.(dh_synth)] .= NaN

            #if plot_dh_as_function_of_time_and_elevation
            if !isnothing(geotiles2plot)

                println(binned_aligned_file)
                param = binned_filled_fileparts(binned_aligned_file)
                surface_mask =param[:surface_mask]

                area_km2 = _geotile_area_km2(surface_mask, geotile_width)

                for mission in missions
                    dh[mission][dh[mission].==0] .= NaN
                end

                # combine icesat and icesat2 for a 4 panel figure
                icesat12 = "ICESat & ICESat 2"
                dh[icesat12] = dh["icesat2"]
                notvalid = isnan.(dh[icesat12])
                dh[icesat12][notvalid] = dh["icesat"][notvalid]

                delete!(dh, "icesat")
                delete!(dh, "icesat2")

                dh["Synthesis"] = dh_synth

                fig_folder = pathlocal[:figures]
                
                if occursin("binned_unfiltered", param.binned_folder)
                    fig_folder = joinpath(fig_folder, "binned_unfiltered")
                else
                    fig_folder = joinpath(fig_folder, "binned")
                end
                
                binned_synthesized_file_parts = splitpath(binned_synthesized_file)
                plot_save_path_prefix = joinpath(fig_folder, replace(binned_synthesized_file_parts[end], ".jld2" => ""))

                plot_elevation_time_multimission_geotiles(
                    dh;
                    geotiles2plot,
                    area_km2,
                    colorrange=(-20, 20),
                    colorbar_label="height anomaly",
                    hypsometry=true,
                    area_averaged=true,
                    plots_show=true,
                    plots_save=true,
                    plot_save_path_prefix,
                    plot_save_format=".png",
                    mission_order=plot_order["synthesis"],
                )
            end

            if isnothing(single_geotile_test) && isnothing(dh_override)
                save(binned_synthesized_file, Dict("dh_hyps" => dh_synth, "dh_hyps_err" => dh_synth_err))
                println("$binned_aligned_file synthesized: $(round(Int,time() -t1))s")
            else
                return dh_synth, dh_synth_err
            end
        end
    end
end

"""
    geotile_zonal_area_hyps(ras, ras_range, zone_geom, geotile_ids; persistent_attribute=:RGIId)

Calculate hypsometric area distribution for each geotile and glacier polygon intersection.

# Arguments
- `ras`: Input raster (typically elevation data)
- `ras_range`: Range of elevation bins for hypsometric binning
- `zone_geom`: DataFrame containing glacier polygons with :geometry column
- `geotile_ids`: List of geotile IDs to process
- `persistent_attribute`: Column name in zone_geom to preserve in output (default: :RGIId)

# Returns
- DataFrame with columns:
  - geotile: Geotile ID
  - persistent_attribute: Preserved identifier from input (e.g. RGIId)
  - area_km2: Area in km² for each geotile-glacier intersection

# Notes
- Crops raster data to intersection of geotiles and glacier polygons
- Calculates areas using cell-by-cell accumulation into elevation bins
- Returns areas in km² after converting from m²
"""
function geotile_zonal_area_hyps(ras, ras_range, zone_geom, geotile_ids; persistent_attribute=:RGIId)

    df = DataFrame()

    @showprogress dt = 1 desc = "Calculating hypsometry for each geotile and glacier ..." for geotile_id0 in geotile_ids
        #geotile = first(geotile_ids)
        t1 = time()
        bounding_polygon = extent2rectangle(GeoTiles.extent(geotile_id0))
        
        geom_name = GI.geometrycolumns(zone_geom)[1]

        index = GO.intersects.(zone_geom[!, geom_name], Ref(bounding_polygon))
        if !any(index)
            println("skipping $(geotile_id0): no intersecting polygons")
            continue
        end

        zone_geom0 = DataFrame(zone_geom[index, :])
        zone_geom0[!, :geometry] .= geotile_id0

        # do a double crop as cropping the polygons themselves is just way to complicated right now
        ras0 = Rasters.crop(ras, to=zone_geom0)
        ras0 = Rasters.crop(ras0, to=bounding_polygon)

        ras0 = ras0[:, :]# having issues going from lazy to inmem ... not logical but this works

        rs = RasterStack(ras0, Rasters.cellarea(ras0))

        geo_column, = GeoInterface.geometrycolumns(zone_geom0)
        area_m2 = mapzonal(Base.Fix2(geotile_zonal_area_hyps, ras_range), identity, rs; of=zone_geom0[!, geo_column])
        zone_geom0[!, :area_km2] = area_m2 * (1E-3)^2
        zone_geom0[!, :geotile] .= geotile_id0

        append!(df, zone_geom0[:, [:geotile, persistent_attribute, :area_km2]]; promote=true)
    end
    return df
end
"""
    geotile_zonal_area_hyps(x, y, bins::StepRange)

Accumulate values into elevation bins.

# Arguments
- `x`: Elevation values
- `y`: Values to accumulate (typically areas)
- `bins`: Elevation bin edges

# Returns
Vector of accumulated values per elevation bin
"""
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

"""
    geotile_zonal_area_hyps(nt, x_bin_edges)

Accumulate values into elevation bins from a NamedTuple input.

# Arguments
- `nt`: NamedTuple with elevation values in first field and values to accumulate in second field
- `x_bin_edges`: Elevation bin edges as a StepRange

# Returns
Vector of accumulated values per elevation bin
"""
function geotile_zonal_area_hyps(nt, x_bin_edges)
    binned = geotile_zonal_area_hyps(getindex.(nt, 1), getindex.(nt, 2), x_bin_edges)
    return binned
end
"""
    geotile_grouping!(geotiles0, glaciers, min_area_km2; geotile_groups_manual=nothing)

Group geotiles that share large glaciers and apply manual grouping overrides.

# Arguments
- `geotiles0`: DataFrame with geotile information
- `glaciers`: DataFrame with glacier geometries and areas
- `min_area_km2`: Minimum glacier area threshold for grouping
- `geotile_groups_manual`: Optional manual grouping specifications

# Returns
Modified DataFrame with group assignments and rectangular geometries
"""
function geotile_grouping!(geotiles0, glaciers, min_area_km2; geotile_groups_manual=nothing)
    # Initialize columns to store discharge indices and geotile intersections
    geotiles0[!, :glacier_overap_ind] .= [Int64[]]
    geotiles0[!, :geotile_intersect] .= [Int64[]]

    glaciers_large = glaciers[glaciers.Area.>min_area_km2, :]
    geometry_column = first(GI.geometrycolumns(glaciers_large))

    # Build spatial index tree for efficient intersection testing
    tree = SortTileRecursiveTree.STRtree(glaciers_large[:, geometry_column]; nodecapacity=3)

    # For each geotile, find intersecting glaciers
    @showprogress dt = 1 desc = "Match glaciers_large with geotiles ..." Threads.@threads for gt in eachrow(geotiles0)
        # Query tree for potential intersecting glaciers
        potential_idxs = SortTileRecursiveTree.query(tree, gt.extent)

        # Find glaciers that actually intersect the geotile
        intersecting = findall(Base.Fix1(GO.intersects, extent2rectangle(gt.extent)),
            view(glaciers_large[:, geometry_column], potential_idxs))

        gt.glacier_overap_ind = potential_idxs[intersecting]
    end

    # Group geotiles that share glaciers
    @showprogress dt = 1 desc = "Group geotiles by glaciers_large ..." for (i, gt) in enumerate(eachrow(geotiles0))
        if .!isempty(gt.glacier_overap_ind)
            # Find other geotiles that share glaciers with this one
            intersecting_geotiles = .!isdisjoint.(Ref(gt.glacier_overap_ind), geotiles0.glacier_overap_ind)
            intersecting_geotiles[i] = false

            if any(intersecting_geotiles)
                gt.geotile_intersect = findall(intersecting_geotiles)
            end
        end
    end

    # Identify connected groups of geotiles
    connectivity = geotiles0.geotile_intersect
    geotiles0[!, :group] = connected_groups(geotiles0.geotile_intersect)


    # Assign new group numbers to manual overrides
    if !isnothing(geotile_groups_manual)
        group0 = maximum(geotiles0.group)
        for grp in geotile_groups_manual
            for g in grp
                geotiles0[findfirst(geotiles0.id .== g), :group] = group0
            end
            group0 += 1
        end
    end

    # Convert geotile extents to rectangles 
    geometry_column = first(GI.geometrycolumns(geotiles0))
    geotiles0[!, geometry_column] = extent2rectangle.(geotiles0.extent)

    return geotiles0[:, Not(:extent, :glacier_overap_ind, :geotile_intersect)]
end

"""
    geotile_grouping(; surface_mask="glacier", min_area_km2=100, geotile_width=2, force_remake_before=nothing)

Group geotiles based on glacier connectivity and save results to file.

# Arguments
- `surface_mask`: Surface mask type to use for filtering (default: "glacier")
- `min_area_km2`: Minimum glacier area threshold for grouping (default: 100)
- `geotile_width`: Width of geotiles in degrees (default: 2)
- `force_remake_before`: Date before which to force recalculation (default: nothing)

# Returns
- `geotiles_out`: GeoDataFrame containing geotile groups with columns :geometry, :id, :group

# Description
This function identifies groups of geotiles that are connected by large overlapping glaciers.
It either loads existing grouping results from file or calculates new groupings by:
1. Filtering geotiles to only include those with glacier coverage
2. Loading glacier geometry data
3. Calling geotile_grouping! to perform the actual grouping logic
4. Saving results to an Arrow file for future use

The function uses caching to avoid recalculating groupings unless the output file
doesn't exist or is older than the force_remake_before date.
"""
function geotile_grouping(; surface_mask="glacier", min_area_km2=100, geotile_width=2, force_remake_before=nothing)

    geotile_groups_fn = joinpath(pathlocal[:data_dir], "project_data", "geotile_groups_$(surface_mask)_$(min_area_km2).arrow")
    if !isfile(geotile_groups_fn) || (!isnothing(force_remake_before) && Dates.unix2datetime(mtime(geotile_groups_fn)) < force_remake_before)
        geotiles0 = geotiles_w_mask(geotile_width)
        geotiles0 = geotiles0[geotiles0[:, "$(surface_mask)_frac"].>0.0, :]

        geometry0 = GeoDataFrames.read(pathlocal[Symbol("$(surface_mask)_individual")])

        # identify groups of geotiles that are connected by large overlapping glaciers
        geotiles0 = geotile_grouping!(geotiles0, geometry0, min_area_km2; geotile_groups_manual=geotile_groups_forced())

        geotiles_out = geotiles0[:, [:geometry, :id, :group]]

        # Write geotile groups to file, excluding some columns
        GeoDataFrames.write(geotile_groups_fn, geotiles_out; crs=GFT.EPSG(4326))
    else
        geotiles_out = GeoDataFrames.read(geotile_groups_fn)
    end

    return geotiles_out

end


"""
    discharge2smb(glaciers; discharge2smb_max_latitude=-60, discharge2smb_equilibrium_period=(Date(1979), Date(2000)))

Calculate discharge from surface mass balance (SMB) trends for glaciers below a specified latitude.

# Arguments
- `glaciers`: DataFrame containing glacier data with columns :CenLat, :CenLon, :smb, and :area_km2
- `discharge2smb_max_latitude`: Maximum latitude threshold for calculating discharge (default: -60)
- `discharge2smb_equilibrium_period`: Time period for equilibrium calculation (default: 1979-2000)

# Returns
- `discharge0`: DataFrame with columns :latitude, :longitude, :discharge_gtyr, :discharge_err_gtyr, :frontal_ablation_gtyr

# Description
This function estimates glacier discharge by calculating SMB trends during an equilibrium period.
For glaciers below the specified latitude threshold, it fits a linear trend to the SMB data
and converts this to discharge in gigatons per year based on glacier area.
"""
function discharge2smb(glaciers; discharge2smb_max_latitude=-60, discharge2smb_equilibrium_period=(Date(1979), Date(2000)))
    index_discharge2smb = glaciers.CenLat .< discharge2smb_max_latitude

    discharge0 = DataFrame(latitude=glaciers[index_discharge2smb, :CenLat], longitude=glaciers[index_discharge2smb, :CenLon], discharge_gtyr=NaN, discharge_err_gtyr=NaN, frontal_ablation_gtyr=NaN)

    ddate = dims(glaciers[1, :smb], :date)
    decyear = decimalyear.(ddate)
    Δdecyear = decyear .- decyear[1]
    index_date = (ddate .>= discharge2smb_equilibrium_period[1]) .& (ddate .<= discharge2smb_equilibrium_period[2])

    for (i, glacier) in enumerate(eachrow(glaciers[index_discharge2smb, :]))
        y = glacier.smb[index_date]
        fit = curve_fit(offset_trend, Δdecyear[index_date], y .- mean(y), offset_trend_p)
        discharge0[i, :discharge_gtyr] = fit.param[2] .* sum(glacier.area_km2) / 1000
    end
    return discharge0
end

"""
    dh2dv(dh, area_km2) -> dv

Convert elevation change (dh) to volume change (dv) using area information from geotiles.

# Arguments
- `dh`: DimArray of elevation changes with dimensions (geotile, date, height)
- `geotiles`: DataFrame containing geotile information with columns :id and :area_km2

# Returns
- `dv`: DimArray of volume changes with dimensions (geotile, date)

# Description
This function calculates volume change by multiplying elevation change by the corresponding 
area for each geotile and height bin, then summing across all height bins. The result is 
converted from km³ to Gt by dividing by 1000.
"""
function dh2dv(dh, area_km2)

    dv = fill(NaN, dims(dh, :date))
    if all(isnan.(dh))
        return dv
    else
        (index_date, index_height) = validrange(.!isnan.(dh))

        foo = @d (dh[index_date, index_height] .* (area_km2[index_height] ./ 1000))
        dv[index_date] = vec(sum(foo; dims=:height))
        return dv
    end
end

function dh2dv_geotile(dh, area_km2)

    dgeotile = dims(dh, :geotile)
    ddate = dims(dh, :date)
    dv = fill(NaN, dgeotile, ddate)
    for geotile in dgeotile
        dv[geotile = At(geotile)] = dh2dv(dh[geotile = At(geotile)], area_km2[geotile = At(geotile)])
    end
    return dv
end


"""
    geotile2glacier!(glaciers, var0; varname)

Aggregates geotile data to glacier-level time series using area-weighted averaging.

# Arguments
- `glaciers`: DataFrame with glacier metadata including geotile assignments and area_km2
- `var0`: DimArray of geotile data with dimensions (geotile, date, height)
- `varname`: Column name for storing the resulting glacier time series

# Returns
Modified `glaciers` DataFrame with new `varname` column containing glacier-level time series
"""
function geotile2glacier!(glaciers, var0; varname)
    # Pre-allocate output arrays
    ddate = dims(var0, :date)
    glaciers[!, varname] = [fill(NaN, ddate) for _ in 1:nrow(glaciers)]

    # Group glaciers by geotile for faster lookup
    glaciers_by_geotile = DataFrames.groupby(glaciers, :geotile)
    geotiles = getindex.(keys(glaciers_by_geotile), "geotile")

    # Process each geotile in parallel
    @showprogress dt = 1 desc = "Geotile \"variable\" to glacier..." Threads.@threads for geotile in geotiles
        # Get glaciers in this geotile
        geotile_glaciers = glaciers_by_geotile[(geotile,)]

        # Get valid data range for this geotile
        v0 = var0[At(geotile), :, :]
        valid_rows, _ = validrange(.!isnan.(v0))

        # Skip if no valid data
        isempty(valid_rows) && continue

        # Pre-allocate temporary array for this geotile's calculations
        out = Vector{Float64}(undef, length(valid_rows))

        # Process each glacier in this geotile
        for glacier in eachrow(geotile_glaciers)
            # Get valid area indices
            area_index = glacier.area_km2 .> 0
            any(area_index) || continue

            # Get valid column range and total area
            valid_cols, = validrange(area_index)
            total_area = sum(view(glacier.area_km2, valid_cols))

            # Calculate area weights once
            area_weights = view(glacier.area_km2, valid_cols) ./ total_area

            # Calculate weighted averages for each row
            @views for (i, row) in enumerate(eachrow(v0[valid_rows, valid_cols]))
                out[i] = sum(row .* area_weights)
            end

            # Assign results to glacier
            glacier[varname][valid_rows] = out
        end
    end

    return glaciers
end


"""
    global_discharge_filled(;
        surface_mask="glacier",
        discharge_global_fn=paths[:discharge_global],
        downscale_to_glacier_method="area",
        gemb_run_id=4,
        discharge2smb_max_latitude=-60,
        discharge2smb_equilibrium_period=(Date(1979), Date(2000)),
        pscale=1,
        Δheight=0,
        geotile_width=2,
        force_remake_before=nothing
    )

Calculate global glacier discharge by combining measured discharge data with estimated discharge 
from surface mass balance for unmeasured glaciers.

This function processes glacier discharge data by:
1. Loading existing discharge data if available and recent enough
2. Filtering geotiles to only those containing glaciers
3. Loading surface mass balance (SMB) data from GEMB model runs
4. Downscaling SMB data to individual glaciers using either hypsometry or area-based methods
5. Estimating discharge for unmeasured glaciers using SMB-to-discharge relationships
6. Combining measured and estimated discharge data
7. Saving the complete global discharge dataset

# Arguments
- `surface_mask`: Surface mask type to filter geotiles (default: "glacier")
- `discharge_global_fn`: Path to save/load global discharge data
- `downscale_to_glacier_method`: Method for downscaling geotile data to glaciers ("hypsometry" or "area")
- `gemb_run_id`: GEMB model run identifier
- `discharge2smb_max_latitude`: Maximum latitude for SMB-to-discharge conversion (default: -60°)
- `discharge2smb_equilibrium_period`: Date range for equilibrium period calculations
- `pscale`: Precipitation scaling factor (default: 1.0)
- `Δheight`: Height adjustment in meters (default: 0)
- `geotile_width`: Width of geotiles in degrees (default: 2)
- `force_remake_before`: Optional DateTime to force regeneration of files created before this date

# Returns
- DataFrame containing global glacier discharge data with columns for glacier identifiers and discharge values

# Notes
- For glaciers without direct discharge measurements, discharge is estimated using surface mass balance data
- Negative discharge values are set to zero
- The function supports both hypsometry-based and area-based downscaling methods
"""
function global_discharge_filled(;
    surface_mask="glacier",
    discharge_global_fn=paths[:discharge_global],
    downscale_to_glacier_method="area",
    gemb_run_id = 4,
    discharge2smb_max_latitude=-60,
    discharge2smb_equilibrium_period=(Date(1979), Date(2000)),
    pscale=1,
    Δheight=0,
    geotile_width=2,
    force_remake_before=nothing
)

    if isfile(discharge_global_fn) && isnothing(force_remake_before)
        discharge = FileIO.load(discharge_global_fn, "discharge")
    elseif isfile(discharge_global_fn) && Dates.unix2datetime(mtime(discharge_global_fn)) > force_remake_before
        discharge = FileIO.load(discharge_global_fn, "discharge")
    else
        # Get geotiles containing glaciers and initialize parameters:
        # - Filter to only tiles with glaciers
        # - Set precipitation scaling to 1.0 (no scaling)
        # - Set height adjustment to 0

        # Load geotiles with glacier mask and filter to only include tiles with glacier coverage
        geotiles = geotiles_w_mask(geotile_width)
        geotiles = geotiles[(geotiles[!, "$(surface_mask)_frac"].>0.0), :]

         gembinfo = gemb_info(; gemb_run_id)
        filename_gemb_geotile_filled_dv = replace(gembinfo.filename_gemb_combined, ".jld2" => "_geotile_filled_extra_extrap_dv.jld2")

        # load example geotile to get dimensions
        if downscale_to_glacier_method == "hypsometry"

            surface_mask_hypsometry = geotile_hypsometry(geotiles, surface_mask; dem_id=:cop30_v2, force_remake_before=nothing)

            filename_gemb_geotile_filled = replace(filename_gemb_geotile, ".jld2" => "_$(surface_mask_hypsometry.geotile[1])_filled.jld2")
            smb = load(filename_gemb_geotile_filled, "smb")

            dgeotile = Dim{:geotile}(geotiles.id)
            ddate = dims(smb, :date)
            dheight = dims(smb, :height)

            var1 = fill(NaN, dgeotile, ddate, dheight)

            @time for varname = ["smb"]
                @showprogress desc = "Populate $(varname)..." Threads.@threads for geotile in geotiles.id
                    filename_gemb_geotile_filled = replace(filename_gemb_geotile, ".jld2" => "_$(geotile)_filled.jld2")
                    var0 = load(filename_gemb_geotile_filled, varname)
                    var1[At(geotile), :, :] = var0[:, :, At(pscale), At(Δheight)]
                end

                # Convert geotile-level data to per-glacier values
                surface_mask_hypsometry = geotile2glacier!(surface_mask_hypsometry, var1; varname)
            end
        elseif downscale_to_glacier_method == "area"

            smb = load(filename_gemb_geotile_filled_dv, "smb")
            smb = smb[:, :, At(pscale), At(Δheight)]

            ddate = dims(smb, :date)
            surface_mask_hypsometry[!, "smb"] = [fill(NaN, ddate) for _ in 1:nrow(surface_mask_hypsometry)]

            for geotile in unique(surface_mask_hypsometry.geotile)
                #geotile = "lat[+28+30]lon[+082+084]"
                gindex = findall(geotile .== surface_mask_hypsometry.geotile)
                gt_dv = smb[At(geotile), :]
                garea = sum.(surface_mask_hypsometry[gindex, :area_km2])
                gweighting = garea ./ sum(garea)
                for (i, ig) in enumerate(gindex)
                    surface_mask_hypsometry[ig, "smb"][:] = gt_dv * gweighting[i] ./ garea[i] * 1000 # convert from km3 to mwe
                end
            end
        else
            error("downscale_to_glacier_method must be either \"hypsometry\" or \"area\"")
        end

        # For unmeasured glaciers, estimate discharge using surface mass balance
        # Only process glaciers with non-zero area
        discharge0 = discharge2smb(
            surface_mask_hypsometry[sum.(surface_mask_hypsometry.area_km2).>0, :];
            discharge2smb_max_latitude,
            discharge2smb_equilibrium_period
        )

        # Combine estimated discharge with measured discharge data
        discharge = glacier_discharge()
        discharge = vcat(discharge, discharge0)

        # Set any negative discharge values to zero
        discharge[discharge.discharge_gtyr.<0, :discharge_gtyr] .= 0

        # Save the combined measured and estimated discharge data
        save(globaldischarge_fn, Dict("discharge" => discharge))
    end
    return discharge
end



"""
    geotile_synthesis_gembfit_dv(; path2runs, discharge, gemb_run_id, geotile_width, force_remake_before)

Synthesize geotile-level data with GEMB fit parameters and compute calibrated timeseries.

This function processes multiple runs to create calibrated geotile-level datasets by:
1. Loading and aligning geotiles for each surface mask
2. Calculating glacier hypsometry and identifying geotiles with glaciers
3. Applying GEMB fit parameters (pscale, Δheight) to calibrate variables
4. Computing derived variables (discharge, volume change, mass change)
5. Averaging certain variables by geotile groups
6. Adding regional glacier inventory (RGI) identifiers

# Arguments
- `path2runs`: Vector of file paths to binned synthesized data files
- `discharge`: DataFrame containing glacier discharge measurements
- `gemb_run_id`: String identifier for the GEMB run
- `geotile_width`: Width of geotiles in degrees
- `force_remake_before`: Optional DateTime to force regeneration of files created before this date

# Returns
- None (saves processed data to files with "_gembfit_dv.jld2" suffix)

# Variables Processed
- Input variables: dv_altim, runoff, fac, smb, rain, acc, melt, ec, refreeze
- Derived variables: discharge, dv (volume change), dm (mass change), dm_altim (altimetry mass change)
"""
function geotile_synthesis_gembfit_dv(path2runs, discharge; gemb_run_id, geotile_width, force_remake_before)

    gembinfo = gemb_info(; gemb_run_id)
    filename_gemb_geotile_filled_dv = replace(gembinfo.filename_gemb_combined, ".jld2" => "_geotile_filled_extra_extrap_dv.jld2")

    # Load elevation change data to get geotile dimensions
    dgeotile = dims(FileIO.load(path2runs[1], "dh_hyps"), :geotile)

    # Initialize dictionaries to store geotiles and glaciers for each surface mask
    geotiles = Dict()
    has_glacier = Dict()
    area_km2 = Dict()

    params = binned_filled_fileparts.(path2runs)
    surface_masks = unique(getindex.(params, :surface_mask))

    # Process each surface mask to load aligned geotiles and calculate hypsometry
    for surface_mask in surface_masks
        # Load and align geotiles with the synthesized data dimensions
        geotiles[surface_mask] = _geotile_load_align(;
            surface_mask="glacier",
            geotile_width,
            geotile_order=collect(dgeotile)
        )

        area_km2[surface_mask] = _geotile_area_km2(surface_mask, geotile_width)

        # Calculate glacier hypsometry for the aligned geotiles
        glaciers0 = geotile_hypsometry(
            geotiles[surface_mask],
            surface_mask;
            dem_id=:cop30_v2,
            force_remake_before=nothing
        )

        has_glacier[surface_mask] = [id in unique(glaciers0.geotile) for id in geotiles[surface_mask].id]
    end

    # Load synthesized data and GEMB fit parameters
    volume2mass = δice / 1000

    # removed Threads.@threads due to memory issues
    @showprogress desc = "Computing calibrated geotile level data timeseries for all runs..." for binned_synthesized_file in path2runs

        binned_synthesized_dv_file = replace(binned_synthesized_file, ".jld2" => "_gembfit_dv.jld2")

        run_parameters = binned_filled_fileparts(binned_synthesized_file)
        surface_mask = run_parameters.surface_mask

        if isfile(binned_synthesized_dv_file) || (isnothing(force_remake_before) || Dates.unix2datetime(mtime(binned_synthesized_dv_file)) > force_remake_before)
            printstyled("    -> Skipping $(binned_synthesized_dv_file) because it was created after force_remake_before:$force_remake_before\n"; color=:light_green)
            continue
        else

            synthesized_gemb_fit = replace(binned_synthesized_file, ".jld2" => "_gemb_fit.arrow")

            gemb_fit = GeoDataFrames.read(synthesized_gemb_fit)

            geotiles0 = copy(geotiles[surface_mask])
            geotile_groups = GeoDataFrames.read(geotile_groups_fn)

            # Filter to only include geotiles containing glaciers
            geotiles0 = geotiles0[has_glacier[surface_mask], :]

            # Define variables to process
            varnames = ["dv_altim", "runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze"]

            # Add group assignments from geotile_groups
            geotiles0[!, :group] .= 0
            for geotile in eachrow(geotiles0)
                geotile.group = geotile_groups[findfirst(isequal(geotile.id), geotile_groups.id), :group]
            end

            # Add GEMB fit parameters and mass conversion factors
            geotiles0[!, :pscale] .= 0.0
            geotiles0[!, :Δheight] .= 0.0
            geotiles0[!, :mie2cubickm] .= 0.0

            # load to get date vector
            var0 = load(filename_gemb_geotile_filled_dv, "smb")
            ddate = dims(var0, :date)
            Δdecyear = decimalyear.(ddate) .- decimalyear.(ddate[1])
            geotiles0[!, :discharge] .= [fill(NaN, length(ddate)) for _ in 1:nrow(geotiles0)]

            for geotile in eachrow(geotiles0)
                fit_index = findfirst(isequal(geotile.id), gemb_fit.id)
                geotile.pscale = gemb_fit[fit_index, :pscale]
                geotile.Δheight = gemb_fit[fit_index, :Δheight]
                area_km20 = sum(geotile.area_km2)
                geotile.mie2cubickm = area_km20 / 1000  # Convert meters ice equivalent to cubic kilometers

                # include discharge average over total glacier area [i.e. save in units of mie]
                index = within.(Ref(geotile.extent), discharge.longitude, discharge.latitude)
                geotile.discharge = sum(discharge.discharge_gtyr[index]) / volume2mass * Δdecyear # Gt to units of mie
            end

            # add discharge date metadata
            colmetadata!(geotiles0, "discharge", "date", collect(ddate), style=:note)
            colmetadata!(geotiles0, "discharge", "units", "km3 [i.e.]", style=:note)

            # Process each variable

            for varname in varnames

                # Special handling for height change data
                if varname == "dv_altim"
                    dh = load(binned_synthesized_file, "dh_hyps")
                    var0 = dh2dv_geotile(dh, area_km2[surface_mask])
                else
                    var0 = load(filename_gemb_geotile_filled_dv, varname)
                end

                # Initialize time series arrays
                ddate = dims(var0, :date)
                geotiles0[!, varname] = [fill(NaN, length(ddate)) for _ in 1:nrow(geotiles0)]

                # Add date metadata for time series
                colmetadata!(geotiles0, varname, "date", collect(ddate), style=:note)
                colmetadata!(geotiles0, varname, "units", "km3 [i.e.]", style=:note)

                # Apply GEMB scaling and convert units for each geotile
                for geotile in eachrow(geotiles0)
                    geotile[varname] = var0[At(geotile.id), :, At(geotile.pscale), At(geotile.Δheight)]
                end
            end

            # modeled dv
            geotiles0[!, :dv] = (geotiles0[:, :smb] .- geotiles0[:, :discharge] .+ geotiles0[:, :fac])
            colmetadata!(geotiles0, "dv", "date", collect(ddate), style=:note) # to units of km3 assuming an ice density of 910 kg/m3
            colmetadata!(geotiles0, "dv", "units", "km3 [i.e.]", style=:note)

            # modeled dm
            geotiles0[!, :dm] = (geotiles0[:, :smb] .- geotiles0[:, :discharge]) .* volume2mass # to units of Gt
            colmetadata!(geotiles0, "dm", "date", collect(ddate), style=:note)
            colmetadata!(geotiles0, "dm", "units", "Gt", style=:note)

            # altimetry dm
            geotiles0[!, :dm_altim] = copy(geotiles0[!, :dv_altim])
            for r in eachrow(geotiles0)
                model = LinearInterpolation(r.fac, decimalyear.(colmetadata(geotiles0, "fac", "date")))
                fac = model(decimalyear.(colmetadata(geotiles0, "dv_altim", "date")))
                r.dm_altim = (r.dv_altim .- fac) .* volume2mass # in units of Gt
            end

            colmetadata!(geotiles0, "dm_altim", "date", colmetadata(geotiles0, "dv_altim", "date"), style=:note)
            colmetadata!(geotiles0, "dm_altim", "units", "Gt", style=:note)

            varnames_all = vcat(varnames, "discharge", "dv", "dm", "dm_altim")

            # `dichage`, `dm` and `dv` must be average by geotile groups as they are not valid for single geotiles
            vars2average = ["discharge", "dv", "dm", "dv_altim", "dm_altim"]
            gdf = groupby(geotiles0, :group)

            for g in gdf
                #g = gdf[4]
                if nrow(g) > 1
                    for varn in vars2average
                        #varn = "dv_altim"
                        foo = zeros(eltype(g[1, varn]), length(g[1, varn]))
                        for r in eachrow(g)
                            foo .+= r[varn]
                        end
                        foo = foo ./ sum(g[:, :mie2cubickm])

                        for r in eachrow(g)
                            r[varn] = foo .* r[:mie2cubickm]
                        end
                    end
                end
            end

            #### now compute regional estimates
            geotiles0[!, :rgiid] .= Int8(0)
            for i = 1:19
                rgi = "rgi$i"
                geotiles0[geotiles0[:, rgi].>0, :rgiid] .= i
                geotiles0 = geotiles0[:, Not(rgi)]
            end

            # sanity check
            for varname in varnames_all
                if any([all(isnan.(v)) for v in geotiles0[:, varname]])
                    println("$(varname) has all NaNs")
                end
            end

            FileIO.save(binned_synthesized_dv_file, Dict("geotiles" => geotiles0))
        end
    end
end


