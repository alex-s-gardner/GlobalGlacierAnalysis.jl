"""
    geotile_synthesize(path2runs; error_file, mission_error, missions2update, force_remake_before, missions2include)

Synthesize elevation change data across multiple satellite missions for improved spatial coverage and temporal consistency.

# Arguments
- `path2runs`: Vector of file paths to binned altimetry data files
- `error_file`: Path to save geotile synthesis error assessment (default: "/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2")
- `mission_error`: Mission-specific error estimates for data quality assessment
- `missions2update`: Force update of mission-specific geotile data (default: nothing)
- `force_remake_before`: Date before which to force recreation of existing files (default: nothing)
- `missions2include`: Array of mission names to include in synthesis (default: ["hugonnet", "gedi", "icesat", "icesat2"])

# Description
This function performs a two-step synthesis process:
1. Calculates geotile synthesis error to assess data quality and consistency across missions
2. Combines elevation change data from different satellite missions into a single dataset

The synthesis improves spatial coverage and temporal consistency by leveraging multiple altimetry datasets.

# Returns
- Nothing; writes synthesized files to disk (or returns dh_synth, dh_synth_err when single_geotile_test is set).

# Examples
```julia
julia> geotile_synthesize(path2runs; missions2include=["hugonnet", "gedi", "icesat", "icesat2"])
julia> # Writes synthesis files next to each run file
```
"""
function geotile_synthesize(path2runs;
    error_file="/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2",
    mission_error = fill(0.0, Dim{:mission}(["hugonnet", "gedi", "icesat", "icesat2"])),
    missions2update=nothing,
    single_geotile_test=nothing,
    missions2include=nothing,
    force_remake_before=nothing,
)

    # Calculate geotile synthesis error to assess data quality and consistency across missions
    dh_hyps_error, _ = geotile_synthesis_error(;
        path2runs,
        outfile=error_file,
        mission_error,
        missions2update,
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

# Examples
```julia
julia> land_error = binned_mad_mission(binned_file)
julia> gedi_error = land_error[At("gedi")]
```
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
        missions2update=nothing,
        single_geotile_test=nothing,
        force_remake_before=nothing
    )

Calculate standard deviation across model runs as an error metric for each satellite mission.

# Arguments
- `path2runs`: Vector of paths to input data files containing binned elevation change data
- `outfile`: Path for saving computed error metrics and file tracking information
- `mission_error`: Floor error values for each mission to ensure minimum error estimates
- `missions2update`: Optional list of specific missions to update instead of processing all missions
- `single_geotile_test`: Single geotile identifier for testing purposes (output not saved to file)
- `force_remake_before`: DateTime threshold for forcing recalculation of existing results

# Returns
- `dh_hyps_error`: Dictionary containing error estimates for each mission across geotiles, dates, and heights
- `files_included`: Vector of file paths that were processed

# Description
Computes error metrics by calculating standard deviation across multiple model runs for each satellite mission.
The function handles incremental updates, file caching, and single geotile testing scenarios.

# Examples
```julia
julia> (dh_hyps_error, files_included) = geotile_synthesis_error(; path2runs, outfile=error_file)
julia> # Saves to outfile unless single_geotile_test is set
```
"""
function geotile_synthesis_error(;
    path2runs,
    outfile="/mnt/bylot-r3/data/binned/2deg/geotile_synthesis_error.jld2",
    mission_error = fill(0.0, Dim{:mission}(["hugonnet", "gedi", "icesat", "icesat2"])),
    missions2update=nothing,
    single_geotile_test=nothing, #geotiles_golden_test[1], #geotiles_golden_test[2],
    geotile_width = 2,
    force_remake_before=nothing,
)

    if .!isnothing(single_geotile_test)
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)], OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!"
    end

    if isfile(outfile) && isnothing(force_remake_before) && isnothing(missions2update) && isnothing(single_geotile_test)

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

    

    if isfile(outfile) && !isnothing(missions2update) && issetequal(path2runs, load(outfile, "files_included")) && isnothing(single_geotile_test)
        missions2update = missions2update
        printstyled("\n    -> Updating selected missions in geotile_synthesis_error: $(missions2update)\n")
    elseif isfile(outfile) && !isnothing(missions2update) && !issetequal(path2runs, load(outfile, "files_included"))
        @warn "Updating all missions in geotile_synthesis_error as input file list does not match exisiting: $(outfile)"
    end

    if .!isnothing(single_geotile_test)
        dgeotile = Dim{:geotile}([single_geotile_test])
    end

    # load example file to get dimension and mission info
    dh = FileIO.load(path2runs[1], "dh_hyps")
    missions = collect(keys(dh))
    ddate = dims(dh[missions[1]], :date)
    dheight = dims(dh[missions[1]], :height)
    dgeotile = dims(dh[missions[1]], :geotile)
    dfile = Dim{:file}(path2runs)

    if isnothing(missions2update)
        missions2update = missions
    end

    # initialize output
    dh_all_std = Dict()
    if !isnothing(missions2update)
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

    if isnothing(single_geotile_test)
        # save the error so that it can be used in the synthesis of individual model runs
        save(outfile, Dict("dh_hyps_error" => dh_hyps_error, "files_included" => files_included))
    else
        params = binned_filled_fileparts.(files_included)
        surface_masks = unique(getindex.(params, :surface_mask))
        surface_mask = surface_masks[1]

        area_km2 = _geotile_area_km2(;surface_mask, geotile_width)
       
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
        plots_show=false,
        plots_save=false,
        single_geotile_test=nothing,
        force_remake_before=nothing
    )

Synthesize elevation change data from multiple satellite missions into a single dataset.

# Arguments
- `path2runs`: Paths to input data files
- `dh_err`: Mission error estimates for each mission
- `missions2include`: Mission names to include in synthesis
- `geotile2plot`: Specific geotile to plot (optional)
- `plots_show`: Whether to generate diagnostic plots
- `plots_save`: Whether to save diagnostic plots
- `single_geotile_test`: Single geotile for testing (optional)
- `force_remake_before`: Recompute if output file exists but was created before this date

# Returns
- Path to synthesized output file with weighted average elevation change data and uncertainty estimates (or dh_synth, dh_synth_err when single_geotile_test is set).

# Examples
```julia
julia> geotile_synthesize_runs(; path2runs, dh_err=dh_hyps_error, missions2include=["hugonnet", "gedi", "icesat", "icesat2"])
julia> # Writes synthesis files for each run
```
"""
function geotile_synthesize_runs(;
    path2runs,
    dh_err,
    missions2include=["hugonnet", "gedi", "icesat", "icesat2"],
    geotiles2plot = nothing,
    single_geotile_test=nothing, #geotiles_golden_test[1], #geotiles_golden_test[2],
    dh_override = nothing, # used to override the dh_hyps data with a different set of data
    geotile_width = 2,
    plots_show=false,
    plots_save=false,
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

    if isnothing(missions2include)
        missions2include = collect(keys(dh_err))
    end
    mission_specs = project_products(; project_id=:v01)
    
    
    if isnothing(single_geotile_test)
        dgeotile = dims(dh_err[first(missions2include)], :geotile)
    else
        dgeotile = Dim{:geotile}([single_geotile_test])
        dheight = dims(dh_err[first(missions2include)], :height)
        ddate = dims(dh_err[first(missions2include)], :date)
        for mission in missions2include
            foo = zeros(dgeotile, ddate, dheight)
            foo[1,:,:] =  dh_err[mission][geotile=At(single_geotile_test)]
            dh_err[mission] = foo
        end
    end

    geotile_rectangles = extent2rectangle.(geotile_extent.(collect(dgeotile)))
    
    # convert error to weights
    w = copy(dh_err)

    for mission in missions2include
        #w0 = w[mission];
        w[mission] = 1 ./ (dh_err[mission] .^ 2)

        ## mask Goutside of observational limits

        mission_extent = (X=mission_specs[Symbol(mission)].longitude_limits, Y=mission_specs[Symbol(mission)].latitude_limits)
        mission_rectangle = extent2rectangle(mission_extent)

        exclude_mission = .!GO.intersects.(geotile_rectangles, Ref(mission_rectangle))
        w[mission][exclude_mission, :, :] .= 0
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
                dgeotile = dims(dh[first(missions2include)], :geotile)

                for mission in missions2include
                    w00[mission] = w00[mission][geotile=At(collect(dgeotile))]
                end
            else
                dh = FileIO.load(binned_aligned_file, "dh_hyps")
                w00 = deepcopy(w)
            end

            if !isnothing(single_geotile_test)
                for mission in missions2include
                    foo = zeros(dgeotile, ddate, dheight)
                    foo[1,:,:] = dh[mission][geotile=At(single_geotile_test)]
                    dh[mission] = foo
                end
            end

            dh_synth = fill(0.0, dims(dh[first(missions2include)]))
            w_synth = copy(dh_synth)

            for mission in missions2include
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

                area_km2 = _geotile_area_km2(;surface_mask, geotile_width)

                for mission in missions2include
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
                    plots_show,
                    plots_save,
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

# Examples
```julia
julia> df = geotile_zonal_area_hyps(ras, 0:100:5000, zone_geom, geotile_ids; persistent_attribute=:RGIId)
julia> # DataFrame with geotile, RGIId, area_km2 per intersection
```

# Notes
- Crops raster data to intersection of geotiles and glacier polygons
- Calculates areas using cell-by-cell accumulation into elevation bins
- Returns areas in km² after converting from m²
"""
function geotile_zonal_area_hyps(ras, ras_range, zone_geom, geotile_ids; persistent_attribute=:RGIId)

    df = DataFrame()

    @showprogress dt = 1 desc = "Calculating hypsometry for each geotile and glacier ..." for geotile_id0 in geotile_ids
        #geotile = first(geotile_ids)
        bounding_polygon = extent2rectangle(GeoTiles.extent(geotile_id0))
        
        geom_name = GI.geometrycolumns(zone_geom)[1]

        index = GO.intersects.(zone_geom[!, geom_name], Ref(bounding_polygon))
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

        geo_column, = GeoInterface.geometrycolumns(zone_geom0)
        area_m2 = mapzonal(Base.Fix2(geotile_zonal_area_hyps, ras_range), identity, rs; of=zone_geom0[!, geo_column])
        zone_geom0[!, :area_km2] = area_m2 * (1E-3)^2
        zone_geom0[!, :geotile] .= geotile_id0

        append!(df, zone_geom0[:, [:geotile, persistent_attribute, :CenLon, :CenLat, :area_km2]]; promote=true)
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
- Vector of accumulated values per elevation bin

# Examples
```julia
julia> binned = geotile_zonal_area_hyps(elevations, areas, 0:100:5000)
```
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
- Vector of accumulated values per elevation bin

# Examples
```julia
julia> binned = geotile_zonal_area_hyps((elevations, areas), 0:100:5000)
```
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
- Modified DataFrame with group assignments and rectangular geometries

# Examples
```julia
julia> geotile_grouping!(geotiles0, glaciers, 100; geotile_groups_manual=geotile_groups_forced())
julia> # geotiles0 now has :group column
```
"""
function geotile_grouping!(geotiles0, glaciers, min_area_km2; geotile_groups_manual=nothing)
    # Initialize columns to store discharge indices and geotile intersections
    geotiles0[!, :glacier_overap_ind] .= [Int64[]]
    geotiles0[!, :geotile_intersect] .= [Int64[]]

    if "Area" in names(glaciers)
        glaciers_large = glaciers[glaciers.Area.>min_area_km2, :]
    else
        glaciers_large = glaciers[glaciers.area_km2.>min_area_km2, :]
    end
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

# Examples
```julia
julia> geotiles_out = geotile_grouping(; surface_mask="glacier", min_area_km2=100)
julia> # GeoDataFrame with :geometry, :id, :group
```
"""
function geotile_grouping(; surface_mask="glacier", min_area_km2=100, geotile_width=2, force_remake_before=nothing)

    geotile_groups_fn = joinpath(pathlocal[:data_dir], "project_data", "geotile_groups_$(surface_mask)_$(min_area_km2).arrow")
    if !isfile(geotile_groups_fn) || (!isnothing(force_remake_before) && Dates.unix2datetime(mtime(geotile_groups_fn)) < force_remake_before)
        geotiles0 = geotiles_mask_hyps(surface_mask, geotile_width)
        geotiles0[!, :area_km2] = sum.(geotiles0[:, "$(surface_mask)_area_km2"])

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

# Examples
```julia
julia> discharge0 = discharge2smb(glaciers; discharge2smb_max_latitude=-60)
julia> # DataFrame with latitude, longitude, discharge_gtyr
```
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

# Examples
```julia
julia> dv = dh2dv(dh, area_km2)
julia> # DimArray (geotile, date) in Gt
```
"""
function dh2dv(dh, area_km2)

    dv = fill(NaN, dims(dh, :date))
    index_height = area_km2 .> 0
    if !any(index_height)
        dv .= 0
        return dv
    elseif all(isnan.(dh))
        warn("dh2dv: all dh is NaN")
        return dv
    else
        (index_date, _) = validrange(.!isnan.(dh))
        foo = @d (dh[index_date, index_height] .* (area_km2[index_height] ./ 1000))
        dv[index_date] = dropdims(sum(foo; dims=:height), dims=:height)
        return dv
    end
end

"""
    dh2dv_geotile(dh, area_km2)

Convert elevation change to volume change for each geotile.

# Arguments
- `dh`: DimArray of elevation changes with dimensions (geotile, date, height)
- `area_km2`: DimArray of areas with dimensions (geotile, height)

# Returns
- DimArray of volume changes with dimensions (geotile, date) in km³ (ice equivalent).

# Examples
```julia
julia> dv = dh2dv_geotile(dh, area_km2)
julia> dv_geotile_1 = dv[At(geotile_id), :]
```
"""
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
- Modified `glaciers` DataFrame with new `varname` column containing glacier-level time series

# Examples
```julia
julia> geotile2glacier!(glaciers, var0; varname=:runoff)
julia> # glaciers now has glaciers.runoff filled
```
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
        gemb=DimStack(),
        discharge2smb_max_latitude=-60,
        discharge2smb_equilibrium_period=(Date(1979), Date(2000)),
        pscale=1,
        mscale=1,
        geotile_width=2,
        force_remake_before=nothing,
        force_remake_before_hypsometry=nothing
    )

Compute and return a global glacier discharge DataFrame by merging observed and modeled values.

This function assembles a best-estimate, comprehensive global discharge dataset by combining:
1. Measured glacier discharge (from direct observations where available)
2. Estimated discharge for all other glaciers based on GEMB-modeled surface mass balance (SMB)

The routine performs the following steps:
- Loads a saved discharge file if it exists and is up to date;
- Otherwise, loads geotiles and identifies those with glacier coverage,
- Loads SMB results from a GEMB run, applying requested precipitation and temperature scaling,
- Aggregates GEMB SMB timeseries for each glacier, area-weighted by the glacier's share of its geotile,
- Converts SMB to discharge for unmeasured glaciers using a specified empirical relationship,
- Concatenates the observed and estimated discharge values,
- Ensures no negative values (clips at zero) and that Antarctic total discharge is within a plausible range,
- Saves the result for future fast loading.

# Arguments
- `surface_mask`: (String) Geotile/glacier mask type (e.g. "glacier"). Default: "glacier".
- `discharge_global_fn`: (String) Path to read/write the combined discharge file. Default: `paths[:discharge_global]`
- `gemb`: The GEMB output dictionary (e.g., as constructed with `gemb_ensemble_dv`). Default: empty `DimStack()`.
- `discharge2smb_max_latitude`: (Float) Only glaciers south of this latitude have discharge estimated via SMB-to-discharge when not observed. Default: -60.
- `discharge2smb_equilibrium_period`: (Tuple{Date,Date}) Date range for SMB equilibrium (long-term mean). Default: (1979, 2000).
- `pscale`: (Real) Precipitation scaling factor used for GEMB data extraction. Default: 1.
- `mscale`: (Real) Elevation offset (m) to sample GEMB data at when extracting timeseries. Default: 1.
- `geotile_width`: (Integer) Geotile side length in degrees. Default: 2.
- `force_remake_before`: (Union{Nothing,DateTime}) Force regeneration if file is older than this. Default: nothing.
- `force_remake_before_hypsometry`: (Union{Nothing,DateTime}) Like above, but only for internal DEM/hypsometry products.

# Returns
- DataFrame: Single table with one row per glacier (or region), with columns such as `id`, `discharge_gtyr`, geographic location, and other metadata.

# Details
- For glaciers with direct discharge measurements, those values are used.
- For unmeasured glaciers, discharge is estimated by downscaling GEMB-modeled SMB based on glacier area and converting SMB to discharge with a region-dependent empirical relation.
- Results are clipped so that discharge cannot be negative.
- Antarctic discharge sanity check: final total must be 80–110 Gt/yr, otherwise an error is raised.

# Examples
```julia
julia> discharge = global_discharge_filled(; surface_mask="glacier", geotile_width=2)
julia> # DimStack or DataFrame with discharge per geotile
```

# Notes
- Discharge is in units of Gt/yr.
- The function can be forced to remake its output if `force_remake_before` is newer than the last save time.
- Downscaling to individual glaciers uses area-weighting.
"""
function global_discharge_filled(;
    surface_mask="glacier",
    discharge_global_fn=paths[:discharge_global],
    gemb = DimStack(),
    discharge2smb_max_latitude=-60,
    discharge2smb_equilibrium_period=(Date(1979), Date(2000)),
    pscale=1,
    mscale=1,
    geotile_width=2,
    force_remake_before=nothing,
    force_remake_before_hypsometry=nothing
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
        geotiles = _geotile_load_align(;
            surface_mask,
            geotile_width,
            only_geotiles_w_area_gt_0=true
        )

        smb = gemb[:smb][pscale=At(pscale), mscale=At(mscale)]

        surface_mask_hypsometry = geotile_hypsometry(geotiles, surface_mask; dem_id=:cop30_v2, force_remake_before=force_remake_before_hypsometry)

        ddate = dims(smb, :date)
        surface_mask_hypsometry[!, "smb"] = [fill(0.0, ddate) for _ in 1:nrow(surface_mask_hypsometry)]

        Threads.@threads for geotile in geotiles.id
            gindex = findall(geotile .== surface_mask_hypsometry.geotile)
            
            if isempty(gindex)
                continue
            end

            gt_dv = smb[geotile=At(geotile)]

            if all(isnan.(gt_dv))
                continue
            end

            garea = sum.(surface_mask_hypsometry[gindex, :area_km2])
            if sum(garea) == 0
                continue
            end

            gweighting = garea ./ sum(garea)
            for (i, ig) in enumerate(gindex)
                if garea[i] == 0
                    continue
                end
                surface_mask_hypsometry[ig, "smb"][:] = gt_dv * gweighting[i] ./ garea[i] .* 1000 # convert from km3 to mwe
            end
        end

        # For unmeasured glaciers, estimate discharge using surface mass balance
        # Only process glaciers with non-zero area
        discharge0 = discharge2smb(
            surface_mask_hypsometry;
            discharge2smb_max_latitude,
            discharge2smb_equilibrium_period
        )

        # Combine estimated discharge with measured discharge data
        discharge = glacier_discharge()
        discharge = vcat(discharge, discharge0)

        # Set any negative discharge values to zero
        discharge[discharge.discharge_gtyr.<0, :discharge_gtyr] .= 0

        # check that the Antarctic discharge is between 80 and 110 Gt/yr
        antarctic_discharge = sum(discharge[discharge.latitude.<-60, :discharge_gtyr])
        if antarctic_discharge < 80 || antarctic_discharge > 110
            error("Antarctic discharge =  $antarctic_discharge Gt/yr, should be between 80 and 110 Gt/yr")
        end

        # Save the combined measured and estimated discharge data
        save(discharge_global_fn, Dict("discharge" => discharge))
    end


    discharge = discharge2geotile(discharge, dims(gemb, :geotile); mass2volume=true)
    
    return discharge
end



"""
    geotile_synthesis_gembfit_dv(path2runs, discharge, gemb; geotile_width, force_remake_before, geotile_grouping_min_feature_area_km2=100)

Synthesize geotile-level data with GEMB fit parameters and compute calibrated timeseries.

This function processes multiple runs to create calibrated geotile-level datasets by:
1. Loading and aligning geotiles for each surface mask.
2. Calculating glacier hypsometry and identifying geotiles with glaciers.
3. Applying GEMB fit parameters (pscale, mscale) to calibrate variables.
4. Computing derived variables (discharge, volume change, mass change).
5. Averaging certain variables by geotile groups.
6. Adding regional glacier inventory (RGI) identifiers.

# Arguments
- `path2runs::Vector{String}`: Vector of file paths to binned synthesized data files.
- `discharge::DataFrame`: DataFrame containing glacier discharge measurements.
- `gemb`: GEMB model output (dictionary of variables, e.g., from `gemb_ensemble_dv`).
- `geotile_width::Int`: Width of geotiles in degrees.
- `force_remake_before::Union{Nothing, DateTime}`: Optional DateTime to force regeneration of files created before this date.
- `geotile_grouping_min_feature_area_km2::Int=100`: Minimum feature area (km²) for geotile grouping.

# Returns
- Nothing. Saves processed data to files with "_gembfit_dv.jld2" suffix.

# Variables Processed
- Input variables: dv_altim, runoff, fac, smb, rain, acc, melt, ec, refreeze
- Derived variables: discharge, dv (volume change), dm (mass change), dm_altim (altimetry mass change)

# Examples
```julia
julia> geotile_synthesis_gembfit_dv(path2runs, discharge, gemb; geotile_width=2, force_remake_before=nothing)
julia> # Writes _gembfit_dv.jld2 files for each run
```
"""
function geotile_synthesis_gembfit_dv(path2runs, discharge, gemb; geotile_width, force_remake_before, geotile_grouping_min_feature_area_km2 = 100)
   
    # Load elevation change data to get geotile dimensions

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
            surface_mask,
            geotile_width,
            only_geotiles_w_area_gt_0=true
        )

        area_km2[surface_mask] = _geotile_area_km2(; surface_mask, geotile_width)

        has_glacier[surface_mask] = dropdims(sum(area_km2[surface_mask], dims=:height), dims=:height) .> 0
    end

    # Load synthesized data and GEMB fit parameters
    @showprogress desc = "Computing calibrated geotile level data timeseries for all runs..." Threads.@threads for binned_synthesized_file in path2runs

        binned_synthesized_dv_file = replace(binned_synthesized_file, ".jld2" => "_gembfit_dv.jld2")

        if (isfile(binned_synthesized_dv_file) && (isnothing(force_remake_before)) || ((Dates.unix2datetime(mtime(binned_synthesized_dv_file)) > force_remake_before)))
            printstyled("    -> Skipping $(binned_synthesized_dv_file) because it was created after force_remake_before: $force_remake_before\n"; color=:light_green)
            continue
        else
            run_parameters = binned_filled_fileparts(binned_synthesized_file)
            surface_mask = run_parameters.surface_mask
            gemb_fit = GeoDataFrames.read(replace(binned_synthesized_file, ".jld2" => "_gembfit.arrow"))

            geotiles0 = individual_geotile_synthesis_gembfit_dv(binned_synthesized_file, gemb, gemb_fit, discharge, geotiles[surface_mask], area_km2[surface_mask])

            FileIO.save(binned_synthesized_dv_file, Dict("geotiles" => geotiles0))
        end
    end
end

"""
    ensemble_area_average_height_anomalies(path2runs_synthesized_all_ensembles; geotiles2extract=nothing)

Compute area-averaged height anomaly ensemble time series for selected geotiles over all runs.

# Arguments
- `path2runs_synthesized_all_ensembles`: Vector of file paths to synthesized run output files (e.g., `.jld2`) for each ensemble member/run.
- `geotiles2extract`: Vector of geotile identifiers to extract and compute area-averaged anomalies for.
    If `nothing`, uses all geotiles found in the first run file.

# Returns
- `dh_area_averaged`: A DimArray of area-averaged height anomalies for each run and geotile, with dimensions (`:run`, `:geotile`, `:date`).

# Description
For each run (i.e., each path in `path2runs_synthesized_all_ensembles`), this function:
- Loads the height anomaly array (`dh_hyps`) for that run.
- For each geotile, computes the area-weighted average height anomaly across heights.
- Assembles the results into a 3D array whose dimensions are (run, geotile, date).

Area weighting is based on surface-specific geotile areas loaded per unique surface mask.

# Examples
```julia
julia> dh_area_averaged = ensemble_area_average_height_anomalies(path2runs; geotiles2extract=["lat[+60+62]lon[-142-140]"])
julia> # DimArray (run, geotile, date)
```
"""
function ensemble_area_average_height_anomalies(path2runs_synthesized_all_ensembles; geotiles2extract = nothing)
    # read in example file
    dh0 = load(path2runs_synthesized_all_ensembles[1], "dh_hyps")
    (dgeotile, ddate, dheight) = dims(dh0)

    drun = Dim{:run}(path2runs_synthesized_all_ensembles)

    params = binned_filled_fileparts.(path2runs_synthesized_all_ensembles);
    surface_masks_unique = unique(getindex.(params, :surface_mask))

    area_km2 = Dict()
    for surface_mask in surface_masks_unique
        area_km2[surface_mask] = _geotile_area_km2(; surface_mask, params[1].geotile_width)
    end

    if !isnothing(geotiles2extract)
        dgeotile = Dim{:geotile}(geotiles2extract)
    end

    dh_area_averaged = fill(NaN, drun, dgeotile, ddate)
    # this takes about 2 minutes to run
    for (param, path2file) in zip(params, path2runs_synthesized_all_ensembles)
        dh0 = load(path2file, "dh_hyps")

        for geotile in dgeotile
            dh_area_averaged[run=At(path2file), geotile=At(geotile)] = dh_area_average(dh0[geotile=At(geotile)], area_km2[param.surface_mask][geotile=At(geotile)])
        end
    end
    return dh_area_averaged
end

"""
    discharge2geotile(discharge, geotiles)

Aggregate discharge variables from a DataFrame to geotile regions.

# Arguments
- `discharge`: DataFrame containing discharge data with columns for longitude, latitude, and one or more discharge variables.
- `geotiles`: DataFrame containing geotile information, including :id and :extent columns.

# Returns
- `discharge0`: DimArray with dimensions (:varname, :geotile), where each entry is the sum of the discharge variable within the geotile extent.

# Description
For each geotile, this function sums the values of each discharge variable for all points within the geotile's extent.
The result is a DimArray indexed by variable name and geotile id.

# Examples
```julia
julia> discharge0 = discharge2geotile(discharge, geotiles)
julia> runoff_by_geotile = discharge0[At("discharge_gtyr"), :]
```
"""
function discharge2geotile(discharge, geotiles)
    dgeotile = Dim{:geotile}(geotiles.id)
    dvarname = Dim{:varname}(setdiff(names(discharge), ["longitude", "latitude", "extent"]))

    discharge0 = zeros(dvarname, dgeotile)

    for geotile_row in eachrow(geotiles)
        index = within.(Ref(geotile_row.extent), discharge.longitude, discharge.latitude)
        if any(index)
            for varname in dvarname
                discharge0[varname=At(varname), geotile=At(geotile_row.id)] = sum(discharge[index, varname])
            end
        end
    end

    return discharge0
end

"""
    discharge2geotile(discharge, dgeotile::Dim; mass2volume=false)

Aggregate discharge from a DataFrame to geotiles using a dimension of geotile IDs.

# Arguments
- `discharge`: DataFrame with longitude, latitude, and discharge columns (e.g. discharge_gtyr)
- `dgeotile`: Dim{:geotile} of geotile identifiers
- `mass2volume`: If true, convert mass (Gt/yr) to volume (km³/yr) using ice density (default: false)

# Returns
- DimStack with one layer per discharge variable, dimensioned by geotile (units Gt/yr or km³/yr).

# Examples
```julia
julia> dgeotile = Dim{:geotile}(geotiles.id)
julia> discharge0 = discharge2geotile(discharge, dgeotile; mass2volume=true)
```
"""
function discharge2geotile(discharge, dgeotile::Dim; mass2volume = false)

    if mass2volume
        mass2volume = 1000 / δice
        units = "km3/yr"
    else
        mass2volume = 1.0
        units = "Gt/yr"
    end

    dvarname = setdiff(names(discharge), ["longitude", "latitude", "extent"])
    outvarname = Symbol.(replace.(dvarname, "_gtyr" => ""))
    

    discharge0 = DimStack([zeros(dgeotile; name) for name in outvarname]; metadata=Dict("units" => units))

    for geotile in dgeotile
        geotile_extent0 = geotile_extent(geotile)
        index = within.(Ref(geotile_extent0), discharge.longitude, discharge.latitude)
        if any(index)
            for i in eachindex(dvarname)
                discharge0[outvarname[i]][geotile=At(geotile)] = sum(discharge[index, dvarname[i]]) * mass2volume
            end
        end
    end   
    return discharge0
end

"""
    individual_geotile_synthesis_gembfit_dv(binned_synthesized_file, gemb, gemb_fit, discharge, geotiles0, area_km2)

Compute calibrated geotile-level time series for a single run using GEMB fit parameters.

Combines altimetry-derived volume change with GEMB-scaled SMB, discharge, and derived variables
for each geotile, applying per-geotile pscale and mscale from gemb_fit.

# Arguments
- `binned_synthesized_file`: Path to binned synthesized altimetry file
- `gemb`: GEMB model output (e.g. DimStack from gemb_ensemble_dv)
- `gemb_fit`: DataFrame with per-geotile fit parameters (pscale, mscale, group)
- `discharge`: Discharge aggregated to geotiles (e.g. from discharge2geotile)
- `geotiles0`: DataFrame of geotiles to process
- `area_km2`: DimArray of area per geotile and height

# Returns
- DataFrame of geotiles with time series columns (dv_altim, runoff, fac, smb, etc.) and metadata.

# Examples
```julia
julia> geotiles0 = individual_geotile_synthesis_gembfit_dv(binned_synthesized_file, gemb, gemb_fit, discharge, geotiles0, area_km2)
```
"""
function individual_geotile_synthesis_gembfit_dv(binned_synthesized_file, gemb, gemb_fit, discharge, geotiles0, area_km2)
    # Filter to only include geotiles containing glaciers
    volume2mass = δice / 1000
    geotiles0 = deepcopy(geotiles0) # geotiles0 is modified within the function

    # Add GEMB fit parameters and mass conversion factors
    geotiles0[!, :pscale] .= 1.0
    geotiles0[!, :mscale] .= 0.0
    geotiles0[!, :mie2cubickm] .= 0.0
    geotiles0[!, :group] .= 0

    dpscale = dims(gemb, :pscale)
    dmscale = dims(gemb, :mscale)

    dv_gemb0 = gemb_dv_sample(dpscale[1], dmscale[1], gemb[:smb][geotile=At(geotiles0.id[1])])
    ddate_gemb = dims(dv_gemb0, :date)
    Δdecyear = decimalyear.(ddate_gemb) .- decimalyear.(ddate_gemb)[1]
    geotiles0[!, :discharge] .= [fill(0.0, length(ddate_gemb)) for _ in 1:nrow(geotiles0)]

    for geotile in eachrow(geotiles0)
        fit_index = findfirst(isequal(geotile.id), gemb_fit.id)
        if isnothing(fit_index)
            continue
        end
        
        geotile.pscale = gemb_fit[fit_index, :pscale]
        geotile.mscale = gemb_fit[fit_index, :mscale]
        geotile.group = gemb_fit[fit_index, :group]
        area_km20 = sum(area_km2[geotile = At(geotile.id)])
        geotile.mie2cubickm = area_km20 / 1000  # Convert meters ice equivalent to cubic kilometers

        # include discharge average over total glacier area [i.e. save in units of mie]
        geotile.discharge = discharge[:discharge][geotile = At(geotile.id)] / volume2mass * Δdecyear # Gt to units of mie
    end

    # add discharge date metadata
    colmetadata!(geotiles0, "discharge", "date", collect(val(ddate_gemb)), style=:note)
    colmetadata!(geotiles0, "discharge", "units", "km3 [i.e.]", style=:note)

    # Process each variable

    ## Populate with altimetry data
    # Special handling for height change data
   
    dh = FileIO.load(binned_synthesized_file, "dh_hyps")
    dv_altim = dh2dv_geotile(dh, area_km2)

    # Initialize time series arrays
    ddate_altim = dims(dv_altim, :date)
    geotiles0[!, :dv_altim] = [fill(0.0, length(ddate_altim)) for _ in 1:nrow(geotiles0)]

    # Add date metadata for time series
    colmetadata!(geotiles0, :dv_altim, "date", collect(val(ddate_altim)), style=:note)
    colmetadata!(geotiles0, :dv_altim, "units", "km3 [i.e.]", style=:note)

    # Apply GEMB scaling and convert units for each geotile
    dgeotile = dims(dv_altim, :geotile)

    for geotile in eachrow(geotiles0)
        if in(geotile.id, dgeotile)
            geotile[:dv_altim] = dv_altim[geotile=At(geotile.id)]
        end
    end

    ## Polulate with GEMB data
    # Initialize time series arrays
    for k in keys(gemb)

        geotiles0[!, Symbol(k)] = [zeros(length(ddate_gemb)) for _ in 1:nrow(geotiles0)]

        # Add date metadata for time series
        colmetadata!(geotiles0, Symbol(k), "date", collect(val(ddate_gemb)), style=:note)
        colmetadata!(geotiles0, Symbol(k), "units", "km3 [i.e.]", style=:note)
    end

    for geotile in eachrow(geotiles0)
        # Apply GEMB scaling and convert units for each geotile
        if in(geotile.id, dgeotile)
            for k in keys(gemb)
                geotile[Symbol(k)][:] = gemb_dv_sample(geotile.pscale, geotile.mscale, gemb[k][geotile=At(geotile.id)]).data
            end
        end
    end

    # modeled dm
    geotiles0[!, :dm] = (geotiles0[:, :dv] .- geotiles0[:, :fac]) .* volume2mass # to units of Gt
    colmetadata!(geotiles0, "dm", "date", collect(val(ddate_gemb)), style=:note)
    colmetadata!(geotiles0, "dm", "units", "Gt", style=:note)

    # altimetry dm
    geotiles0[!, :dm_altim] = copy(geotiles0[!, :dv_altim])

    decyear_fac = decimalyear.(colmetadata(geotiles0, "fac", "date"))
    decyear_altim = decimalyear.(colmetadata(geotiles0, "dv_altim", "date"))
    for r in eachrow(geotiles0)
        itp = DataInterpolations.LinearInterpolation(r.fac, decyear_fac)
        fac = itp(decyear_altim)
        r.dm_altim = (r.dv_altim .- fac) .* volume2mass # in units of Gt
    end

    colmetadata!(geotiles0, "dm_altim", "date", colmetadata(geotiles0, "dv_altim", "date"), style=:note)
    colmetadata!(geotiles0, "dm_altim", "units", "Gt", style=:note)

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


    varnames_all = vcat("dv_altim", String.(keys(gemb))..., "discharge", "dm", "dm_altim")
    # sanity check
    for varname in varnames_all
        if any([all(isnan.(v)) for v in geotiles0[:, varname]])
            println("$(varname) has all NaNs")
        end
    end

    return geotiles0
end
        