#WARNING: NOT TESTED - direct AI translation from PyGEM


"""
# Translated from the Python Glacier Evolution Model (PyGEM)

copyright © 2018 David Rounce <drounce@cmu.edu>

Mass redistribution curves for glacier evolution modeling
"""
# module MassRedistribution

using OrderedCollections
using Dates
using Statistics
using LinearAlgebra

import DimensionalData as DD
import GeoInterface as GI

# Constants
const δice = 910  # density of glacier ice [kg m-3]

"""
    MassRedistributionCurveModel

Glacier geometry updated using mass redistribution curves (delta-h method).

This uses mass redistribution curves from Huss et al. (2010) to update glacier geometry.
"""
Base.@kwdef mutable struct MassRedistributionCurveModel
    flowlines::Vector{Any}
    mb_model::Any = nothing
    y0::Float64 = 0.0
    yr::Float64 = 0.0
    glen_a::Union{Float64,Nothing} = nothing
    fs::Float64 = 0.0
    is_tidewater::Bool = false
    water_level::Float64 = 0.0
    option_areaconstant::Bool = false
    spinupyears::Int = 0
    constantarea_years::Int = 0
    glac_idx_initial::Vector{Vector{Int}} = []
    check_for_boundaries::Bool = true
    calving_k::Float64 = 0.0
    calving_m3_since_y0::Float64 = 0.0
    calving_rate_myr::Float64 = 0.0
    fls::Vector{Any} = []
    mb_step::String = "annual"
    mb_elev_feedback::String = "annual"
end

"""
    run_until(model, y1; run_single_year=false)

Runs the model from the current year up to a given year date y1.
"""
function run_until(model::MassRedistributionCurveModel, y1; run_single_year=false)
    if run_single_year
        updategeometry(model, y1)
    else
        years = floor(Int, model.yr):floor(Int, y1-1)
        for year in years
            updategeometry(model, year)
        end
    end
    
    # Check for domain bounds
    if model.check_for_boundaries
        if model.fls[end].thick[end] > 10
            error("Glacier exceeds domain boundaries, at year: $(model.yr)")
        end
    end
    
    # Check for NaNs
    for fl in model.fls
        if any(.!isfinite.(fl.thick))
            error("NaN in numerical solution.")
        end
    end
    
    model.yr = y1
    return nothing
end

"""
    run_until_and_store(model, y1; run_path=nothing, diag_path=nothing, store_monthly_step=nothing)

Runs the model and returns intermediate steps in data structures.
"""
function run_until_and_store(model::MassRedistributionCurveModel, y1; 
                             run_path=nothing, diag_path=nothing, store_monthly_step=nothing)
    if floor(Int, y1) != y1
        error("run_until_and_store only accepts integer year dates.")
    end

    if !model.mb_model.hemisphere
        error("run_until_and_store needs a mass-balance model with an unambiguous hemisphere.")
    end
    
    # Time setup
    yearly_time = floor(Int, model.yr):floor(Int, y1)

    if store_monthly_step === nothing
        store_monthly_step = model.mb_step == "monthly"
    end

    if store_monthly_step
        monthly_time = monthly_timeseries(model.yr, y1)
    else
        monthly_time = floor(Int, model.yr):floor(Int, y1)
    end

    sm = 1  # Hydro month placeholder - would be replaced with actual value

    # Initialize outputs
    ny = length(yearly_time)
    nm = length(monthly_time)
    sects = [fill(NaN, (ny, fl.nx)) for fl in model.fls]
    widths = [fill(NaN, (ny, fl.nx)) for fl in model.fls]
    bucket = [fill(NaN, ny) for _ in model.fls]
    
    # Create dimensional dataset
    run_ds = []
    diag_ds = Dict{String,Any}()
    
    # Dataset attributes
    diag_ds["description"] = "Altim.jl model output"
    diag_ds["version"] = "0.1.0"
    diag_ds["calendar"] = "365-day no leap"
    diag_ds["creation_date"] = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    diag_ds["hemisphere"] = model.mb_model.hemisphere
    diag_ds["water_level"] = model.water_level

    # Data variables
    diag_ds["volume_m3"] = fill(NaN, nm)
    diag_ds["area_m2"] = fill(NaN, nm)
    diag_ds["length_m"] = fill(NaN, nm)
    diag_ds["ela_m"] = fill(NaN, nm)
    diag_ds["volume_bsl_m3"] = fill(NaN, nm)
    diag_ds["volume_bwl_m3"] = fill(NaN, nm)
    
    if model.is_tidewater
        diag_ds["calving_m3"] = fill(NaN, nm)
        diag_ds["calving_rate_myr"] = fill(NaN, nm)
    end

    # Run model with recording
    j = 1
    for (i, yr) in enumerate(yearly_time[1:end-1])
        # Record initial parameters
        if i == 1
            diag_ds["volume_m3"][i] = volume_m3(model)
            diag_ds["area_m2"][i] = area_m2(model)
            diag_ds["length_m"][i] = length_m(model)
            
            if model.is_tidewater
                diag_ds["volume_bsl_m3"][i] = volume_bsl_m3(model)
                diag_ds["volume_bwl_m3"][i] = volume_bwl_m3(model)
            end
        end
        
        # Run model for this year
        run_until(model, yr, run_single_year=true)
        
        # Store section data
        for (s, w, b, fl) in zip(sects, widths, bucket, model.fls)
            s[j, :] = fl.section
            w[j, :] = fl.widths_m
            if model.is_tidewater
                try
                    b[j] = fl.calving_bucket_m3
                catch
                    # Ignore if attribute doesn't exist
                end
            end
        end
        j += 1
        
        # Update diagnostics
        diag_ds["volume_m3"][i+1] = volume_m3(model)
        diag_ds["area_m2"][i+1] = area_m2(model)
        diag_ds["length_m"][i+1] = length_m(model)

        if model.is_tidewater
            diag_ds["calving_m3"][i+1] = model.calving_m3_since_y0
            diag_ds["calving_rate_myr"][i+1] = model.calving_rate_myr
            diag_ds["volume_bsl_m3"][i+1] = volume_bsl_m3(model)
            diag_ds["volume_bwl_m3"][i+1] = volume_bwl_m3(model)
        end
    end

    # Prepare output datasets
    for (i, (s, w, b)) in enumerate(zip(sects, widths, bucket))
        ds = Dict{String,Any}()
        ds["description"] = "Altim.jl model output"
        ds["version"] = "0.1.0"
        ds["calendar"] = "365-day no leap"
        ds["creation_date"] = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
        ds["ts_section"] = s
        ds["ts_width_m"] = w
        if model.is_tidewater
            ds["ts_calving_bucket_m3"] = b
        end
        push!(run_ds, ds)
    end

    # Write outputs if paths provided
    if run_path !== nothing || diag_path !== nothing
        # Implementation would depend on your file format preferences
        # For example, using JLD2 or NetCDF.jl
    end

    return run_ds, diag_ds
end

"""
    updategeometry(model, year; debug=false)

Update glacier geometry for a given year.
"""
function updategeometry(model::MassRedistributionCurveModel, year; debug=false)
    if debug
        println("year: $year")
    end
    
    # Loop over flowlines
    for fl_id in eachindex(model.fls)
        fl = model.fls[fl_id]
        
        # Flowline state
        heights = copy(fl.surface_h)
        section_t0 = copy(fl.section)
        thick_t0 = copy(fl.thick)
        width_t0 = copy(fl.widths_m)
        
        # CONSTANT AREAS
        # Mass redistribution ignored for calibration and spinup years (glacier properties constant)
        if model.option_areaconstant || year < model.spinupyears || year < model.constantarea_years
            # Run mass balance
            glac_bin_massbalclim_annual = get_annual_mb(model.mb_model, heights, fls=model.fls, fl_id=fl_id, 
                                                       year=year, debug=false)
        # MASS REDISTRIBUTION
        else
            # FRONTAL ABLATION
            if model.is_tidewater
                # Frontal ablation (m³ ice)
                fa_m3 = _get_annual_frontalablation(model, heights, fls=model.fls, fl_id=fl_id, 
                                                   year=year, debug=false)
                if debug
                    println("fa_m3_init: $fa_m3")
                    vol_init = sum(model.fls[fl_id].section .* fl.dx_meter)
                    println("  volume init: $(round(vol_init))")
                    println("  volume final: $(round(vol_init-fa_m3))")
                end
                
                # First, remove volume lost to frontal ablation
                glac_idx_bsl = findall((thick_t0 .> 0) .& (fl.bed_h .< model.water_level))
                while fa_m3 > 0 && !isempty(glac_idx_bsl)
                    if debug
                        println("fa_m3_remaining: $fa_m3")
                    end
                    
                    last_idx = glac_idx_bsl[end]
                    
                    if debug
                        println("before: $(round(model.fls[fl_id].section[last_idx], digits=0)) " * 
                                "$(round(model.fls[fl_id].thick[last_idx], digits=0)) " * 
                                "$(round(heights[last_idx], digits=0))")
                    end
                    
                    vol_last = section_t0[last_idx] * fl.dx_meter
                    
                    # If frontal ablation more than bin volume, remove entire bin
                    if fa_m3 > vol_last
                        # Record frontal ablation (m³ w.e.) in mass balance model for output
                        model.mb_model.glac_bin_frontalablation[last_idx, Int(12*(year+1)-1)] = 
                            vol_last * 910 / 1000  # ice density / water density
                        # Update ice thickness and section area
                        section_t0[last_idx] = 0
                        model.fls[fl_id].section = section_t0
                        # Remove volume from frontal ablation "bucket"
                        fa_m3 -= vol_last
                        
                    # Otherwise, remove ice from the section
                    else
                        # Update section to remove frontal ablation
                        section_t0[last_idx] = section_t0[last_idx] - fa_m3 / fl.dx_meter
                        model.fls[fl_id].section = section_t0
                        # Record frontal ablation(m³ w.e.)
                        model.mb_model.glac_bin_frontalablation[last_idx, Int(12*(year+1)-1)] = 
                            fa_m3 * 910 / 1000  # ice density / water density
                        # Frontal ablation bucket now empty
                        fa_m3 = 0
                    end
                    
                    # Update flowline
                    heights = copy(model.fls[fl_id].surface_h)
                    section_t0 = copy(model.fls[fl_id].section)
                    thick_t0 = copy(model.fls[fl_id].thick)
                    width_t0 = copy(model.fls[fl_id].widths_m)
                    
                    if debug
                        println("after: $(round(model.fls[fl_id].section[last_idx], digits=0)) " * 
                                "$(round(model.fls[fl_id].thick[last_idx], digits=0)) " * 
                                "$(round(heights[last_idx], digits=0))")
                        println("  vol final: $(sum(model.fls[fl_id].section .* fl.dx_meter))")
                    end
                    
                    glac_idx_bsl = findall((thick_t0 .> 0) .& (fl.bed_h .< model.water_level))
                end
                
                # Redistribute mass if glacier was not fully removed by frontal ablation
                if !isempty(findall(section_t0 .> 0))
                    # Mass redistribution according to Huss empirical curves
                    # Annual glacier mass balance [m ice s⁻¹]
                    glac_bin_massbalclim_annual = get_annual_mb(model.mb_model, heights, fls=model.fls, fl_id=fl_id, 
                                                               year=year, debug=false)
                    sec_in_year = sum(model.mb_model.dates_table[12*year:12*(year+1)-1, :daysinmonth]) * 24 * 3600
                    
                    _massredistributionHuss(model, section_t0, thick_t0, width_t0, glac_bin_massbalclim_annual, 
                                           model.glac_idx_initial[fl_id], heights, sec_in_year=sec_in_year)
                end
            # Non-tidewater glacier mass redistribution
            else
                # Annual glacier mass balance [m ice s⁻¹]
                glac_bin_massbalclim_annual = get_annual_mb(model.mb_model, heights, fls=model.fls, fl_id=fl_id, 
                                                           year=year, debug=false)
                sec_in_year = sum(model.mb_model.dates_table[12*year:12*(year+1)-1, :daysinmonth]) * 24 * 3600
                
                _massredistributionHuss(model, section_t0, thick_t0, width_t0, glac_bin_massbalclim_annual, 
                                       model.glac_idx_initial[fl_id], heights, sec_in_year=sec_in_year)
            end
        end
        
        # Record glacier properties (volume [m³], area [m²], thickness [m], width [km])
        # Record the next year's properties as well
        # 'year + 1' used so the glacier properties are consistent with mass balance computations
        year_int = Int(year)  # required to ensure proper indexing with run_until_and_store
        glacier_area = fl.widths_m .* fl.dx_meter
        glacier_area[fl.thick .== 0] .= 0
        model.mb_model.glac_bin_area_annual[:, year_int+1] = glacier_area
        model.mb_model.glac_bin_icethickness_annual[:, year_int+1] = fl.thick
        model.mb_model.glac_bin_width_annual[:, year_int+1] = fl.widths_m
        model.mb_model.glac_wide_area_annual[year_int+1] = sum(glacier_area)
        model.mb_model.glac_wide_volume_annual[year_int+1] = sum(fl.section .* fl.dx_meter)
    end
end

"""
    _get_annual_frontalablation(model, heights; year=nothing, fls=nothing, fl_id=nothing, calving_k=nothing, debug=false)

Calculate frontal ablation for a given year.

Returns frontal ablation (m³ ice).
"""
function _get_annual_frontalablation(model::MassRedistributionCurveModel, heights; 
                                    year=nothing, fls=nothing, fl_id=nothing, calving_k=nothing, debug=false)
    # Flowlines and various attributes
    fl = fls[fl_id]
    @assert heights ≈ fl.surface_h
    glacier_area_t0 = fl.widths_m .* fl.dx_meter
    fl_widths_m = fl.widths_m
    fl_section = fl.section
    
    # Ice thickness (average)
    if fl_section !== nothing && fl_widths_m !== nothing
        icethickness_t0 = zeros(size(fl_section))
        icethickness_t0[fl_widths_m .> 0] = fl_section[fl_widths_m .> 0] ./ fl_widths_m[fl_widths_m .> 0]
    else
        icethickness_t0 = nothing
    end

    # Quality control: ensure you only have glacier area where there is ice
    if icethickness_t0 !== nothing
        glacier_area_t0[icethickness_t0 .== 0] .= 0
    end

    # ----- FRONTAL ABLATION -----
    # -- using OGGM's parameterization -----
    # We do calving only if the last glacier bed pixel is below water
    # (this is to avoid calving elsewhere than at the front)
    q_calving = 0.0
    if sum(glacier_area_t0) > 0
        try
            last_above_wl = findlast((fl.surface_h .> model.water_level) .& (fl.thick .> 0))
        catch
            last_above_wl = findlast((fl.bed_h .<= model.water_level) .& (fl.thick .> 0))
        end
        
        if last_above_wl !== nothing
            if fl.bed_h[last_above_wl] < model.water_level
                # Volume [m³] and bed elevation [masl] of each bin
                if debug
                    println("\nyear: $year \n  sea level: $(model.water_level) bed elev: $(round(fl.bed_h[last_above_wl], digits=2))")
                    println("  estimate frontal ablation")
                    println(" min elevation: $(fl.surface_h[last_above_wl])")
                end
                
                # --- The rest is for calving only ---
                model.calving_rate_myr = 0.0
        
                # OK, we're really calving
                section = fl.section
        
                # Calving law
                h = fl.thick[last_above_wl]
                d = h - (fl.surface_h[last_above_wl] - model.water_level)
                k = model.calving_k
                q_calving = k * d * h * fl.widths_m[last_above_wl]
                
                # Max frontal ablation is removing all bins with bed below water level
                glac_idx_bsl = findall((fl.thick .> 0) .& (fl.bed_h .< model.water_level))
                q_calving_max = sum(section[glac_idx_bsl]) * fl.dx_meter
                
                if q_calving > q_calving_max + 1e-6  # tolerance
                    q_calving = q_calving_max
                end
        
                # Add to the bucket and the diagnostics
                model.calving_m3_since_y0 += q_calving
                model.calving_rate_myr = q_calving / section[last_above_wl]
            end
        end
    end

    return q_calving
end

"""
    _massredistributionHuss(model, section_t0, thick_t0, width_t0, glac_bin_massbalclim_annual, 
                            glac_idx_initial, heights; debug=false, hindcast=0, sec_in_year=365*24*3600)

Mass redistribution according to empirical equations from Huss and Hock (2015) accounting for retreat/advance.
"""
function _massredistributionHuss(model::MassRedistributionCurveModel, section_t0, thick_t0, width_t0, 
                                glac_bin_massbalclim_annual, glac_idx_initial, heights; 
                                debug=false, hindcast=0, sec_in_year=365*24*3600)
    # Glacier area [m²]
    glacier_area_t0 = width_t0 .* model.fls[1].dx_meter
    glacier_area_t0[thick_t0 .== 0] .= 0
    
    # Annual glacier-wide volume change [m³]
    # units: [m ice / s] * [s] * [m²] = m³ ice
    glacier_volumechange = sum(glac_bin_massbalclim_annual .* sec_in_year .* glacier_area_t0)
    
    # For hindcast simulations, volume change is the opposite
    if hindcast == 1
        glacier_volumechange = -1 * glacier_volumechange
    end
    
    if debug
        println("\nDebugging Mass Redistribution Huss function\n")
        println("glacier volume change: $glacier_volumechange")
    end
    
    # If volume loss is more than the glacier volume, melt everything and stop here
    glacier_volume_total = sum(model.fls[1].section .* model.fls[1].dx_meter)
    if (glacier_volume_total + glacier_volumechange) < 0
        # Set all to zero and return
        model.fls[1].section .*= 0
        return
    end
    
    # Otherwise, redistribute mass loss/gains across the glacier
    # Determine where glacier exists            
    glac_idx_t0 = findall(model.fls[1].thick .> 0)
    
    # Compute ice thickness [m ice], glacier area [m²], ice thickness change [m ice] after redistribution
    icethickness_change, glacier_volumechange_remaining = _massredistributioncurveHuss(
        model, section_t0, thick_t0, width_t0, glac_idx_t0,
        glacier_volumechange, glac_bin_massbalclim_annual,
        heights, debug=false)
    
    if debug
        println("\nmax icethickness change: $(round(maximum(icethickness_change), digits=3))")
        println("min icethickness change: $(round(minimum(icethickness_change), digits=3))")
        println("volume remaining: $glacier_volumechange_remaining")
        nloop = 0
    end

    # Glacier retreat
    # If glacier retreats (ice thickness == 0), volume change needs to be redistributed over glacier again
    while glacier_volumechange_remaining < 0
        if debug
            println("\n\nGlacier retreating (loop $nloop):")
        end
        
        section_t0_retreated = copy(model.fls[1].section)
        thick_t0_retreated = copy(model.fls[1].thick)
        width_t0_retreated = copy(model.fls[1].widths_m)
        glacier_volumechange_remaining_retreated = copy(glacier_volumechange_remaining)
        glac_idx_t0_retreated = findall(thick_t0_retreated .> 0)
        glacier_area_t0_retreated = width_t0_retreated .* model.fls[1].dx_meter
        glacier_area_t0_retreated[thick_t0 .== 0] .= 0
        
        # Set climatic mass balance for the case when there are less than 3 bins  
        # Distribute the remaining glacier volume change over the entire glacier (remaining bins)
        massbalclim_retreat = zeros(size(thick_t0_retreated))
        massbalclim_retreat[glac_idx_t0_retreated] .= glacier_volumechange_remaining / 
                                                     sum(glacier_area_t0_retreated) / sec_in_year
        
        # Mass redistribution 
        # Apply mass redistribution using Huss' empirical geometry change equations
        icethickness_change, glacier_volumechange_remaining = _massredistributioncurveHuss(
            model, section_t0_retreated, thick_t0_retreated, width_t0_retreated, glac_idx_t0_retreated, 
            glacier_volumechange_remaining_retreated, massbalclim_retreat, heights, debug=false)
        
        # Avoid rounding errors that get loop stuck
        if abs(glacier_volumechange_remaining) < 1
            glacier_volumechange_remaining = 0
        end
        
        if debug
            println("ice thickness change: $icethickness_change")
            println("\nmax icethickness change: $(round(maximum(icethickness_change), digits=3))")
            println("min icethickness change: $(round(minimum(icethickness_change), digits=3))")
            println("volume remaining: $glacier_volumechange_remaining")
            nloop += 1
        end
    end

    # Glacier advances 
    # Based on ice thickness change exceeding threshold
    icethickness_advancethreshold = 0.1  # Placeholder threshold
    
    while any(icethickness_change .> icethickness_advancethreshold)
        if debug
            println("advancing glacier")
        end
        
        # Record glacier area and ice thickness before advance corrections applied
        section_t0_raw = copy(model.fls[1].section)
        thick_t0_raw = copy(model.fls[1].thick)
        width_t0_raw = copy(model.fls[1].widths_m)
        glacier_area_t0_raw = width_t0_raw .* model.fls[1].dx_meter
        
        if debug
            println("\n\nthickness t0: $thick_t0_raw")
            println("glacier area t0: $glacier_area_t0_raw")
            println("width_t0_raw: $width_t0_raw\n\n")
        end
        
        # Index bins that are advancing
        icethickness_change[icethickness_change .<= icethickness_advancethreshold] .= 0
        glac_idx_advance = findall(icethickness_change .> 0)
        
        # Update ice thickness based on maximum advance threshold [m ice]
        model.fls[1].thick[glac_idx_advance] = model.fls[1].thick[glac_idx_advance] .- 
            (icethickness_change[glac_idx_advance] .- icethickness_advancethreshold)
            
        glacier_area_t1 = copy(model.fls[1].widths_m) .* model.fls[1].dx_meter
        
        # Advance volume [m³]
        advance_volume = sum(glacier_area_t0_raw[glac_idx_advance] .* thick_t0_raw[glac_idx_advance]) - 
                         sum(glacier_area_t1[glac_idx_advance] .* model.fls[1].thick[glac_idx_advance])

        if debug
            println("advance volume [m³]: $advance_volume")
        end

        # Handle glacier advance details here...
        # Implementation depends on specific requirements and available data
        # This would include handling terminus position, adding new bins, etc.
        
        # Set icethickness change and advance_volume to break the loop if needed
        if advance_volume <= 0
            icethickness_change[icethickness_change .> 0] .= 0
        end
    end
end

"""
    _massredistributioncurveHuss(model, section_t0, thick_t0, width_t0, glac_idx_t0, glacier_volumechange, 
                                massbalclim_annual, heights; debug=false)

Apply the mass redistribution curves from Huss and Hock (2015).
"""
function _massredistributioncurveHuss(model::MassRedistributionCurveModel, section_t0, thick_t0, width_t0, 
                                     glac_idx_t0, glacier_volumechange, massbalclim_annual, heights; debug=false)
    if debug
        println("\nDebugging mass redistribution curve Huss\n")
    end

    # Apply Huss redistribution if there are at least 3 elevation bands; otherwise, use the mass balance        
    # Glacier area used to select parameters
    glacier_area_t0 = width_t0 .* model.fls[1].dx_meter
    glacier_area_t0[thick_t0 .== 0] .= 0
    
    # Apply mass redistribution curve
    if length(glac_idx_t0) > 3
        # Select the factors for the normalized ice thickness change curve based on glacier area
        glacier_area_sum = sum(glacier_area_t0)
        if glacier_area_sum > 20 * 1e6
            gamma, a, b, c = 6, -0.02, 0.12, 0
        elseif glacier_area_sum > 5 * 1e6
            gamma, a, b, c = 4, -0.05, 0.19, 0.01
        else
            gamma, a, b, c = 2, -0.30, 0.60, 0.09
        end
        
        # Reset variables
        elevrange_norm = zeros(size(glacier_area_t0))
        icethicknesschange_norm = zeros(size(glacier_area_t0))
        
        # Normalized elevation range [-]
        # (max elevation - bin elevation) / (max_elevation - min_elevation)
        if sum(glacier_area_t0) > 0
            max_height = maximum(heights[glac_idx_t0])
            min_height = minimum(heights[glac_idx_t0])
            height_range = max_height - min_height
            
            elevrange_norm[glacier_area_t0 .> 0] = (max_height .- heights[glacier_area_t0 .> 0]) ./ height_range
            
            # Normalized ice thickness change [-]
            icethicknesschange_norm[glacier_area_t0 .> 0] = 
                ((elevrange_norm[glacier_area_t0 .> 0] .+ a).^gamma .+ 
                 b .* (elevrange_norm[glacier_area_t0 .> 0] .+ a) .+ c)
                 
            # Limit the icethicknesschange_norm to between 0 - 1
            icethicknesschange_norm[icethicknesschange_norm .> 1] .= 1
            icethicknesschange_norm[icethicknesschange_norm .< 0] .= 0
            
            # Huss' ice thickness scaling factor, fs_huss [m ice]
            fs_huss = glacier_volumechange / sum(glacier_area_t0 .* icethicknesschange_norm)
            
            if debug
                println("fs_huss: $fs_huss")
            end
            
            # Volume change [m³ ice]
            bin_volumechange = icethicknesschange_norm .* fs_huss .* glacier_area_t0
        else
            bin_volumechange = zeros(size(glacier_area_t0))
        end
    else
        # Compute volume change in each bin based on the climatic mass balance
        bin_volumechange = massbalclim_annual .* glacier_area_t0
    end
    
    if debug
        println("-----\n")
        vol_before = section_t0 .* model.fls[1].dx_meter
    end

    # Update cross sectional area
    # Volume change divided by length (dx); units m²
    section_change = bin_volumechange ./ model.fls[1].dx_meter
    model.fls[1].section = max.(model.fls[1].section .+ section_change, 0)
    
    # Ice thickness change [m ice]
    icethickness_change = model.fls[1].thick .- thick_t0
    
    # Glacier volume
    vol_after = model.fls[1].section .* model.fls[1].dx_meter
    
    if debug
        println("vol_chg_wanted: $(sum(bin_volumechange))")
        println("vol_chg: $(sum(vol_after) - sum(vol_before))")
        println("\n-----")
    end
    
    # Compute the remaining volume change
    bin_volumechange_remaining = bin_volumechange .- (model.fls[1].section .* model.fls[1].dx_meter .- 
                                                     section_t0 .* model.fls[1].dx_meter)
    
    # Remove values below tolerance to avoid rounding errors
    tolerance = 1e-6
    bin_volumechange_remaining[abs.(bin_volumechange_remaining) .< tolerance] .= 0
    
    # Glacier volume change remaining - if less than zero, then needed for retreat
    glacier_volumechange_remaining = sum(bin_volumechange_remaining)
    
    if debug
        println(glacier_volumechange_remaining)
    end

    return icethickness_change, glacier_volumechange_remaining
end

# Helper functions for diagnostics
volume_m3(model::MassRedistributionCurveModel) = sum(sum(fl.section .* fl.dx_meter) for fl in model.fls)
area_m2(model::MassRedistributionCurveModel) = sum(sum(fl.widths_m .* fl.dx_meter) for fl in model.fls)
length_m(model::MassRedistributionCurveModel) = sum(fl.dx_meter * count(fl.thick .> 0) for fl in model.fls)
volume_bsl_m3(model::MassRedistributionCurveModel) = sum(
    sum((fl.section .* fl.dx_meter)[fl.bed_h .< model.water_level]) for fl in model.fls)
volume_bwl_m3(model::MassRedistributionCurveModel) = volume_bsl_m3(model)

# Monthly time series helper
function monthly_timeseries(yr0, yr1)
    start_date = Date(floor(Int, yr0), 1, 1)
    end_date = Date(floor(Int, yr1), 12, 31)
    months = Dates.month(start_date):Dates.month(end_date-Dates.Day(1))+12*(Dates.year(end_date)-Dates.year(start_date))
    return collect(yr0 .+ (months .- 1) ./ 12)
end

export MassRedistributionCurveModel, run_until, run_until_and_store

end # module