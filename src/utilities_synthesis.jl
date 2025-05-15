"""
    geotile_synthesis_error(;
        path2runs,
        outfile="/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2",
        minimum_error=0.05,
        force_remake_before=nothing,
        update_geotile_missions=nothing
    )

Calculate standard deviation across model runs as an error metric for each satellite mission.

# Arguments
- `path2runs`: Vector of paths to input data files
- `outfile`: Path for saving error metrics
- `minimum_error`: Floor value for error estimates
- `force_remake_before`: Date threshold for forcing recalculation
- `update_geotile_missions`: Optional list of specific missions to update

# Returns
Path to the output file containing error estimates
"""
function geotile_synthesis_error(;
    path2runs,
    outfile="/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2",
    minimum_error=0.05,
    force_remake_before=nothing,
    update_geotile_missions = nothing
)

    #outfile = "/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2"

    if !isfile(outfile) || !isnothing(force_remake_before)

        if isfile(outfile) && !isnothing(force_remake_before)
            if Dates.unix2datetime(mtime(outfile)) > force_remake_before
                @warn "Skipping $(outfile) because it was created after $force_remake_before"
                return outfile
            end
        end

        #load example file to get dimension and mission info
        dh = FileIO.load(path2runs[1], "dh_hyps")

        missions = collect(keys(dh))
        ddate = dims(dh[missions[1]], :date)
        dheight = dims(dh[missions[1]], :height)
        dgeotile = dims(dh[missions[1]], :geotile)
        dfile = Dim{:file}(path2runs)

        if isfile(outfile) && !isnothing(update_geotile_missions) && issetequal(path2runs, load(outfile, "files_included"))
            update_selected_missions = true
            missions = update_geotile_missions
            @warn "Updating selected missions in geotile_synthesis_error: $(update_geotile_missions)"
        else
            update_selected_missions = false
            @warn "Updating all missions in geotile_synthesis_error as input file list does not match exisiting: $(outfile)"
        end

        # THIS SHOULD BE LOOPED FOR EACH MISSION AS READING THE FILES IS THE SLOWEST PART [90% of time] as memory is an issue
        dh_all = Dict()
        for mission in missions
            dh_all[mission] = fill(NaN, (dfile, dgeotile, ddate, dheight))
        end

        # this takes about 8 min for length(files) = 96 # this blows up memory usage
        Threads.@threads for filepath in path2runs
            dh = FileIO.load(filepath, "dh_hyps")
            for mission in missions
                dh_all[mission][At(filepath), :, :, :] = dh[mission]
            end
        end

        # calculate standard deviation across runs as an error metric
        dh_all_std = Dict()
        if update_selected_missions
            dh_all_std = load(outfile, "dh_hyps_error")
        else
            dh_all_std = Dict()
            for mission in missions
                dh_all_std[mission] = fill(NaN, (dgeotile, ddate, dheight))
            end
        end
        
        for mission in missions
            #mission = "gedi"

            @showprogress dt = 1 desc = "Calculating standard deviation (error) across runs for $(mission) ..." Threads.@threads for geotile in dgeotile
                #geotile = dgeotile[510]

                # !!! I THINK THIS IS FIXED NOW !!!!
                # TODO: There is an issue where some missions (ICESat2) has different hight range than the other missions... this cuases issues with non-rectangular error matrix
                # See this example
                # heatmap(dh_err["hugonnet"][At("lat[-28-26]lon[-070-068]"),:,:])
                # heatmap!(dh_err["icesat2"][At("lat[-28-26]lon[-070-068]"),:,:])

                valid1 = dropdims(any(.!isnan.(dh_all[mission][:, At(geotile), :, :]), dims=:file), dims=:file)

                if !(any(valid1))
                    continue
                end

                vdate, vheight = Altim.validrange(valid1)

                for date in ddate[vdate]

                    for height in dheight[vheight]
                        var0 = dh_all[mission][:, At(geotile), At(date), At(height)]
                        valid2 = .!isnan.(var0)

                        if any(valid2)
                            dh_all_std[mission][At(geotile), At(date), At(height)] = max(std(var0[valid2]), minimum_error)
                        end
                    end
                end
            end
        end

        # save the error so that it can be used in the synthesis of individual model runs
        save(outfile, Dict("dh_hyps_error" => dh_all_std, "files_included" => path2runs))
    end
    return outfile
end
"""
    geotile_synthesize_runs(;
        path2runs,
        path2geotile_synthesis_error,
        missions2include=["hugonnet", "gedi", "icesat", "icesat2"],
        showplots=false,
        force_remake_before=nothing
    )

Synthesize elevation change data from multiple satellite missions into a single dataset.

# Arguments
- `path2runs`: Paths to input data files
- `path2geotile_synthesis_error`: Path to file with mission error estimates
- `missions2include`: Mission names to include in synthesis
- `showplots`: Whether to generate diagnostic plots
- `force_remake_before`: Recompute if output file exists but was created before this date

# Returns
Path to synthesized output file with weighted average elevation change data and uncertainty estimates
"""
function geotile_synthesize_runs(;
    path2runs,
    path2geotile_synthesis_error,
    missions2include=["hugonnet", "gedi", "icesat", "icesat2"],
    showplots=false,
    force_remake_before=nothing
)

    dh_err = load(path2geotile_synthesis_error, "dh_hyps_error")
    missions = missions2include
    dgeotile = dims(dh_err[missions2include[1]], :geotile)

    # convert error to weights
    w = copy(dh_err)
    for mission in missions
        #w0 = w[mission];
        w[mission] = 1 ./ (dh_err[mission] .^ 2)

        # mask GEDI outside of observational latitude limits
        if mission == "gedi"
            lat_min = getindex.(Altim.geotile_extent.(dgeotile), :min_y)
            lat_max = getindex.(Altim.geotile_extent.(dgeotile), :max_y)
            exclude_gedi = (lat_max .> 51.6) .| (lat_min .< -51.6)
            w[mission][exclude_gedi, :, :] .= 0
        end

        w[mission][isnan.(w[mission])] .= 0
    end

    @showprogress dt = 1 desc = "Synthesizing runs ..." Threads.@threads for binned_aligned_file in path2runs
        #binned_aligned_file = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_b1km_dh_best_cc_meanmadnorm3_v01_filled_ac_p1_aligned.jld2"    

        if showplots
            file_parts = splitpath(binned_aligned_file)
            binned_folder = joinpath(file_parts[1:findfirst(file_parts .== "2deg")])

            fig_folder = joinpath(binned_folder, "figures")
            figure_suffix = replace(file_parts[end], ".jld2" => "")
            binned_synthesized_file = replace(binned_aligned_file, "aligned.jld2" => "synthesized.jld2")
            params = Altim.binned_filled_fileparts(binned_aligned_file)
            geotiles = Altim.geotiles_mask_hyps(params.surface_mask, geotile_width)
        end

        binned_synthesized_file = replace(binned_aligned_file, "aligned.jld2" => "synthesized.jld2")

        if !(isfile(binned_synthesized_file)) || !isnothing(force_remake_before)

             if isfile(binned_synthesized_file) && !isnothing(force_remake_before)
                if Dates.unix2datetime(mtime(binned_synthesized_file)) > force_remake_before
                    @warn "Skipping $(binned_synthesized_file) because it was created after $force_remake_before"
                    continue
                end
            end

            t1 = time()
            dh = FileIO.load(binned_aligned_file, "dh_hyps")

            #heatmap(Altim.dh2dv(dh["hugonnet"], geotiles["glacier_b1km"])[510:511,:])
            #heatmap!(Altim.dh2dv(dh["icesat2"], geotiles["glacier_b1km"])[510:511,:])
            #heatmap!(Altim.dh2dv(dh["gedi"], geotiles["glacier_b1km"])[510:511,:])

            if showplots
                p = Altim.plot_height_time(dh; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="final", fig_folder, figure_suffix, mask=params.surface_mask, showplots)
            end

            dh_synth = fill(0.0, dims(dh[missions[1]]))
            w_synth = copy(dh_synth)

            for mission in missions
                w0 = copy(w[mission])

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
            if showplots
                # load geotiles
                p = Altim.plot_height_time(dh_synth; geotile=geotiles[dh_time_elevation_idx, :], fig_suffix="raw", fig_folder, figure_suffix, mask=params.surface_mask, mission="synthesis", showplots)
            end

            save(binned_synthesized_file, Dict("dh_hyps" => dh_synth, "dh_hyps_err" => dh_synth_err))
            println("$binned_aligned_file synthesized: $(round(Int,time() -t1))s")
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

    for geotile_id0 in geotile_ids
        #geotile = first(geotile_ids)
        t1 = time()
        bounding_polygon = Altim.extent2rectangle(Altim.GeoTiles.extent(geotile_id0))
        
        geom_name = GI.geometrycolumns(zone_geom)[1]

        index = GeometryOps.intersects.(zone_geom[!, geom_name], Ref(bounding_polygon))
        if !any(index)
            println("skipping $(geotile_id0): no intersecting polygons")
            continue
        end

        zone_geom0 = DataFrame(zone_geom[index, :])
        zone_geom0[!, :geotile] .= geotile_id0

        # do a double crop as cropping the polygons themselves is just way to complicated right now
        ras0 = Rasters.crop(ras, to=zone_geom0)
        ras0 = Rasters.crop(ras0, to=bounding_polygon)

        ras0 = ras0[:, :]# having issues going from lazy to inmem ... not logical but this works

        rs = RasterStack(ras0, Rasters.cellarea(ras0))

        area_m2 = Altim.mapzonal(Base.Fix2(geotile_zonal_area_hyps, ras_range), identity, rs; of=zone_geom0[!, :geometry])
        zone_geom0[!, :area_km2] = area_m2 * (1E-3)^2

        println("$(geotile_id0) zonal area hyps done: $(round(Int,time() -t1))s")

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
    tree = STRtree(glaciers_large[:, geometry_column]; nodecapacity=3)

    # For each geotile, find intersecting glaciers
    @showprogress dt = 1 desc = "Match glaciers_large with geotiles ..." Threads.@threads for gt in eachrow(geotiles0)
        # Query tree for potential intersecting glaciers
        potential_idxs = SortTileRecursiveTree.query(tree, gt.extent)

        # Find glaciers that actually intersect the geotile
        intersecting = findall(Base.Fix1(GO.intersects, Altim.extent2rectangle(gt.extent)),
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
    geotiles0[!, :group] = Altim.connected_groups(geotiles0.geotile_intersect)


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
    geotiles0[!, geometry_column] = Altim.extent2rectangle.(geotiles0.extent)


    return geotiles0[:, Not(:extent, :area_km2, :glacier_overap_ind, :geotile_intersect)]
end