# regional_results.jl
#
# This script processes and analyzes glacier mass change data at regional scales.
#
# Main components:
# 1. Data aggregation and uncertainty estimation:
#    - Loads multiple model runs and aggregates geotile data into regional time series
#    - Calculates uncertainties from ensemble spread (95% confidence intervals)
#    - Incorporates GRACE satellite data for comparison
#    - Reformats data for consistent plotting and analysis
#
# 2. Glacier discharge and mass change analysis:
#    - Processes glacier discharge data by region
#    - Maps discharge points to RGI regions using spatial joins
#    - Identifies endorheic (inland-draining) glaciers
#    - Computes mass change rates for each glacier and region
#
# 3. Trend and acceleration analysis:
#    - Fits trends and seasonal cycles to time series
#    - Identifies regions with significant acceleration/deceleration
#    - Calculates mass change under different scenarios (with/without firn air content correction)
#    - Compares altimetry-derived trends with GRACE satellite measurements
#
# 4. Visualization and reporting:
#    - Generates multi-region plots of mass change
#    - Calculates statistics on glacier contribution to sea level rise - Analyzes endorheic basin contributions
#    - Reports regional mass change trends with uncertainties
#
# Key outputs:
# - Regional time series of glacier mass change
# - Trend and acceleration estimates with uncertainties
# - Comparison between different measurement techniques
# - Analysis of endorheic vs. ocean-draining glacier contributions

# NOTE: errors are at the 95% confidence interval from the reference run for this study and 95% confidence interval from GRACE
#function regional_results(path2runs_synthesized, ensemble_reference_file; error_quantile=0.95, error_scaling=1.5)
begin
    error_quantile=0.95
    error_scaling=1.5

    using CSV
    using DataFrames
    using DimensionalData
    using Dates
    using Statistics
    using CairoMakie
    using DataInterpolations
    using Unitful
    import GlobalGlacierAnalysis as GGA

    Unitful.register(GGA.MyUnits)
    using GlobalGlacierAnalysis.MyUnits

    reference_ensemble_file = GGA.reference_ensemble_file;
    ensemble_reference_file = replace(reference_ensemble_file, "_aligned.jld2" => "_synthesized.jld2")
    
    # to include in uncertainty
    paths = GGA.pathlocal
    
    glacier_flux_path = GGA.pathlocal[:glacier_summary]

    discharge_fractional_error = 0.15
    dates2plot = [DateTime(1980, 1, 1),DateTime(2025, 1, 15)] # also dates used for trend fitting
    dates4trend = [DateTime(2000, 3, 1), DateTime(2024, 12, 15)] # use these date for comparing individual model runs to GRACE (only need for rmse anlysis)
    dates_firstdecade = [DateTime(2000, 1, 1), DateTime(2010, 1, 1)] 
    dates_lastdecade = [DateTime(2014, 1, 1), DateTime(2024, 1, 1)] 
    dates_altim_grace_overlap = [DateTime(2002, 5, 15), DateTime(2024, 12, 15)] # use these date for comparing 
    dates_glambie_overlap = [DateTime(2000, 6, 1), DateTime(2023, 1, 1)] # use these date for comparing individual model runs to GRACE (only need for rmse anlysis)
    center_on_dates = DateTime(2002, 5, 1):Month(1):DateTime(2005, 12, 15)

    glacier_summary_file = GGA.pathlocal[:glacier_summary]

    # path2perglacier = replace(ensemble_reference_file, ".jld2" => "_perglacier.jld2")
    path2discharge =  paths[:discharge_global]
    path2rgi_regions = paths[:rgi6_regions_shp]

    path2river_flux = joinpath(paths[:river], "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.arrow")


    path2runs_filled, params = GGA.binned_filled_filepaths(;
        project_id=:v01,
        surface_masks=["glacier", "glacier_rgi7"],
        dem_ids=["best", "cop30_v2"],
        curvature_corrects=[false, true],
        amplitude_corrects=[true],
        binning_methods=["median", "nmad3", "nmad5"],
        fill_params=[1, 2, 3, 4],
        binned_folders=[GGA.analysis_paths(; geotile_width=2).binned, replace(GGA.analysis_paths(; geotile_width=2).binned, "binned" => "binned_unfiltered")],
        include_existing_files_only=true
    )
    path2runs_synthesized = replace.(path2runs_filled, "aligned.jld2" => "synthesized.jld2")
end

begin
    # pre-process data
    #begin #[10 min]
    # load discharge for each RGI [<1s]
    discharge_rgi = GGA.discharge_rgi(path2discharge, path2rgi_regions; fractional_error = discharge_fractional_error);

    # load results for all runs for each RGI [18s]
    runs_rgi = GGA.runs2rgi(path2runs_synthesized);

    # fit trends and acceleration for each model and RGI [33s]
    runs_rgi_fits = GGA.rgi_trends(runs_rgi, discharge_rgi, dates4trend);

    # center data
    runs_rgi = GGA.runs_center!(runs_rgi, center_on_dates);
    GGA.check_for_all_nans(runs_rgi)

    # extract reference run and calculate 95% confidence interval
    regions = GGA.runs_ref_and_err(runs_rgi, ensemble_reference_file; error_quantile, error_scaling);
    GGA.check_for_all_nans(regions)

    region_fits = GGA.region_fit_ref_and_err(runs_rgi_fits, ensemble_reference_file; error_quantile, error_scaling, discharge = discharge_rgi)
    GGA.check_for_all_nans(region_fits)

    # load GRACE data [95% confidence interval]
    grace = GGA.grace_masschange();
    grace[:,:,At(false)] = GGA.regions_center!(grace[:,:,At(false)], center_on_dates)
    regions["dm_grace"] = grace
    GGA.check_for_all_nans(regions)

    # add GlaMBIE results
    regions["dm_glambie"] = GGA.glambie2024()
    regions["dm_glambie"][:,:,At(false)] = GGA.regions_center!(regions["dm_glambie"][:,:,At(false)], center_on_dates)
    GGA.check_for_all_nans(regions)

    # add discarge to region_fits
    using DimensionalData
    (dvarname, drgi, dparameters, derr) = dims(region_fits)
    dvarname_new = Dim{:varname}([collect(dvarname)..., "discharge"])
    region_fits = cat(region_fits, zeros(Dim{:varname}(["discharge"]), drgi, dparameters, derr), dims=dvarname_new)
    region_fits[At("discharge"), :, At("trend"), At(false)] = discharge_rgi[:, At(false)]
    region_fits[At("discharge"), :, At("trend"), At(true)] = discharge_rgi[:, At(true)]
   
    # deteremine endorheic and non-endorheic glacier runoff [3 min]
    # rivers have only been calculated for rgi6 so this will throw an error if using glacier_rgi7
    if !occursin("glacier_rgi7", ensemble_reference_file)
        error("rivers have only been calculated for rgi7 so this will throw an error if using rgi6")
    else
        dm_gt_sinktype = GGA.rgi_endorheic(path2river_flux, glacier_summary_file; dates4trend)
        endorheic_scale_correction = 1 .- dm_gt_sinktype[:,:, At(true)] ./ dm_gt_sinktype[:,:, At(false)]
    end

    # compute trends over grace-altim overlap period
    runs_rgi_fits_grace = GGA.rgi_trends(runs_rgi, discharge_rgi, dates_altim_grace_overlap);
    dvarname2 = Dim{:varname}([collect(dvarname)..., "dm_grace"])
    region_fits_grace = zeros(dvarname2, drgi, dparameters, derr)
    region_fits_grace[At(collect(dvarname)), :, : ,:] = GGA.region_fit_ref_and_err(runs_rgi_fits_grace, ensemble_reference_file; error_quantile, error_scaling, discharge = discharge_rgi)

    # regions 12 and 99 contain NaN so exclude them in trend calculation
    region_fits_grace[At("dm_grace"), At(setdiff(drgi, [12, 99])), : , At(false)] = GGA.rgi_trends(regions["dm_grace"][At(setdiff(drgi, [12, 99])),:, At(false)], dates_altim_grace_overlap);

    # compute trends over the GlaMBIE period
    dvarname2 = Dim{:varname}([collect(dvarname)..., "dm_glambie"])
    region_fits_glambie = zeros(dvarname2, drgi, dparameters, derr)
    runs_rgi_fits_glambie = GGA.rgi_trends(runs_rgi, discharge_rgi, dates_glambie_overlap);
    region_fits_glambie[At(collect(dvarname)), :, :, :] = GGA.region_fit_ref_and_err(runs_rgi_fits_glambie, ensemble_reference_file; error_quantile, error_scaling, discharge = discharge_rgi)
    region_fits_glambie[At("dm_glambie"), :, : , At(false)] = GGA.rgi_trends(regions["dm_glambie"][:,:, At(false)], dates_glambie_overlap);

    runs_rgi_noaltim = filter(p -> p[1] in setdiff(keys(runs_rgi), ["dm_altim", "dv_altim"]), runs_rgi)
    runs_rgi_fits_firstdecade = GGA.rgi_trends(runs_rgi_noaltim, discharge_rgi, dates_firstdecade);
    region_fits_firstdecade = GGA.region_fit_ref_and_err(runs_rgi_fits_firstdecade, ensemble_reference_file; error_quantile, error_scaling, discharge = discharge_rgi);
    runs_rgi_fits_lastdecade = GGA.rgi_trends(runs_rgi_noaltim, discharge_rgi, dates_lastdecade);
    region_fits_lastdecade = GGA.region_fit_ref_and_err(runs_rgi_fits_lastdecade, ensemble_reference_file; error_quantile, error_scaling, discharge = discharge_rgi);
end;

# CREATE OUTPUTS FOR PAPER
begin
    geotiles = GGA._geotile_load_align(; surface_mask="glacier")

    exclude_regions = []#[13, 14, 15]
    params = ["trend", "acceleration", "amplitude", "phase"]
    regional_results = GGA.error_bar_table(region_fits; params = params, rgi_regions = setdiff(drgi, exclude_regions));

    altim_cols = names(regional_results)[occursin.("_altim", names(regional_results))]
    rename!(regional_results, altim_cols .=> replace.(altim_cols, "_altim" => "_synthesis"))
    
    rename!(regional_results, "gsi_trend_[Gt/yr]" => "gsi")

    # add area_km2 to regional_results
    regional_results[!, :area_km2] .= 0.0
    for dfr in eachrow(regional_results)
        if dfr.rgi == 98
            index_rgi = (geotiles[:, "rgi13"] .> 0) .| (geotiles[:, "rgi14"] .> 0) .| (geotiles[:, "rgi15"] .> 0)
        elseif dfr.rgi == 99
            index_rgi = trues(nrow(geotiles))
        else
            index_rgi = geotiles[:, "rgi$(dfr.rgi)"].>0
        end

        dfr.area_km2 = sum(sum(geotiles.area_km2[index_rgi]))
    end

    # dump table to console
    GGA.show_error_bar_table(regional_results)

    # for saving to csv
    regional_results_out = copy(regional_results)
    colnames = names(regional_results_out)
    delete_index = falses(length(colnames))
    delete_cols_containing =["discharge_amplitude", "discharge_phase", "discharge_acceleration", "gsi_phase", "gsi_amplitude", "gsi_acceleration"]
    for delete_col in delete_cols_containing
        delete_index .|= occursin.(Ref(delete_col), colnames)
    end
    regional_results_out = regional_results_out[:, .!delete_index]


    dm_altim_index = findall(occursin.(Ref("dm_altim"), colnames))
    for i in dm_altim_index
        rename!(regional_results_out, colnames[i] => replace(colnames[i], "altim" => "altimetry"))
    end

    output_file = joinpath(paths[:project_dir], "Gardner2025_regional_results.csv")
    CSV.write(output_file, regional_results_out; bom=true)

    endorheic_fraction = DataFrame() 
    endorheic_fraction[!, :rgi] = val(dims(endorheic_scale_correction, :rgi))
    endorheic_fraction[!, :region_name] = GGA.rginum2label.(endorheic_fraction[!, :rgi])
    endorheic_fraction[!, :dm_endorheic_fraction] = GGA.rginum2label.(endorheic_fraction[!, :rgi])
    for varname in dims(endorheic_scale_correction, :varname)
        col_name = "$(varname)_endorheic_fraction"
        endorheic_fraction[!, col_name] = 1 .- endorheic_scale_correction[varname = At(varname)]
    end

    include_regions = .!in.(endorheic_fraction.rgi, Ref(exclude_regions))
    CSV.write(joinpath(paths[:project_dir], "Gardner2025_endorheic_fraction.csv"), endorheic_fraction[include_regions, :]; bom=true)


    # remove regions, trim dates
    regions_out = deepcopy(regions);

    date_range = DateTime(2000, 1, 1).. DateTime(2025, 1, 15)

    # ensure dimensions conform
    var0 = "dm"
    drgi = dims(regions_out[var0], :rgi)
    ddate = dims(regions_out[var0][date=date_range], :date)
    decyear = GGA.decimalyear.(val(ddate))
    derror = dims(regions_out[var0], :error)


    for k in keys(regions_out)
        #k = first(keys(regions_out))

        var1 = regions_out[k][rgi=At(setdiff(drgi, exclude_regions)), date=date_range]
        ddate1 = dims(var1, :date)

        if ddate != ddate1
            var2 = fill(NaN, drgi, ddate, derror)
            decyear2 = GGA.decimalyear.(val(dims(var1, :date)))

            for rgi in drgi

                if !(rgi in val(dims(var1, :rgi)))
                    continue
                end
                
                for error in derror
                    ts = var1[rgi = At(rgi), error = At(error)]
                    ts_interp = DataInterpolations.LinearInterpolation(ts, decyear2;extrapolation = ExtrapolationType.Constant)
                    var2[rgi = At(rgi), error = At(error)] = ts_interp(decyear)
                end
            end
        end

        # add units
        if k == "fac" || occursin("dv", k)
            unit = 1u"km^3"
        else
            unit = 1u"Gt"
        end

        var1 *= unit
        regions_out[k] = var1[rgi = At(setdiff(drgi, exclude_regions))]
    end


    # convert to DimStack
    foo0 = [];
    for k in keys(regions_out)
        
        foo = DimArray(fill(zero(eltype(regions_out[k]))*NaN, drgi, ddate, derror); name=k)
        
        try
            foo[DimSelectors(regions_out[k])] = regions_out[k]
        catch
            @warn "Interpolating $(k) in time to match other regional time series)"
            for rgi in dims(regions_out[k],:rgi)

                for error0 in dims(regions_out[k],:error)

                    # interpolate anomalies using weighted distance (Shepard(2))
                    ts_interp = DataInterpolations.LinearInterpolation(regions_out[k][rgi = At(rgi), error = At(error0)], GGA.decimalyear.(val(dims(regions_out[k],:date))); extrapolation = ExtrapolationType.Linear)

                    foo[rgi = At(rgi), error = At(error0)] = ts_interp(GGA.decimalyear.(val(ddate)))
                end
            end
        end

        push!(foo0, foo)
    end

    foo1 = [];
    for i in eachindex(foo0)

        push!(foo1, dropdims(foo0[i][error = At(false)], dims = (:error,)))
        push!(foo1, DimensionalData.rebuild(dropdims(foo0[i][error = At(true)], dims = (:error,)); name = DimensionalData.name(foo0[i]) * "_error"))
    end

    regions_out_stack = DimStack( foo1... ; metadata = Dict(
        "title" => "glacier mass and volume change time series",
        "version" => "final - " * Dates.format(now(), "yyyy-mm-dd"),))

    filename0 = joinpath(paths[:project_dir], "Gardner2025_regional_timseries.nc")

    GGA.dimstack2netcdf(regions_out_stack, filename0)


    # copy data readme file to project directory
    cp("/home/gardnera/Documents/GitHub/GlobalGlacierAnalysis.jl/src/Gardner2025_DataReadMe.txt", joinpath(paths[:project_dir], "Gardner2025_DataReadMe.txt"); force=true)
end

begin
    exclude_regions = [13, 14, 15]
    
    varnames = ["gsi"]
    params = ["trend"]

    println("----------- from $(dates4trend[1]) to $(dates4trend[2]) -----------")
    regional_results = GGA.error_bar_table(region_fits; varnames, params, rgi_regions=setdiff(drgi, exclude_regions), digits=2)

    GGA.show_error_bar_table(regional_results; cols2display = 1)

end

# dump summary conclusions to console
begin
    drgi = dims(region_fits, :rgi);
    rgis = setdiff(drgi, [99, 13, 14, 15]);

begin
    println("\n-------------------------- ABSTRACT ---------------------------------")
    v = "runoff";
    p = "trend";
    rgi = 99;
    println("glaciers outside the ice sheets contributed: $(round(Int,region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)] + region_fits[varname = At("discharge"), rgi=At(99), parameter= At("trend"), error = At(false)])) ± $(round(Int,hypot(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)], region_fits[varname = At("discharge"), rgi=At(99), parameter= At("trend"), error = At(true)]))) Gt/yr")

    v = "dm";
    p = "trend";
    rgi = 99;
    println("unsustainable contribution (i.e. net loss) of: $(round(Int,region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)])) ± $(round(Int,region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)])) Gt/yr")

    v = "runoff";
    p = "acceleration";
    rgi = 99;
    println("glacier runoff accelerated by: $(round(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)], digits=1)) ± $(round(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)], digits=1)) Gt/yr²")

    v = "acc";
    p = "acceleration";
    rgi = 99;
    println("losses were partially offset by an increase in accumulation: $(round(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)], digits=1)) ± $(round(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)], digits=1)) Gt/yr²")

    v = "dm"
    rgi = 99;
    val0 = region_fits[varname = At(v), rgi = At(rgi), parameter = At("trend"), error = At(false)]*endorheic_scale_correction[At(v), At(rgi)];
    println("we determine that $(round(Int, val0)) Gt/yr of glacier loss contributed $(round((val0 / GGA.ocean_area_km2* 1E6), digits=2)) mm/yr SLR)")

    val0 = region_fits[varname = At(v), rgi = At(rgi), parameter = At("trend"), error = At(false)]*(1. - endorheic_scale_correction[At(v), At(rgi)]);
    println("with $(round(Int, val0)) Gt yr-1 remaining within landlocked endorheic basins")
end

begin
    println("\n----------------------Glacier freshwater flux from satellite data --------------------")
    rgis = setdiff(drgi, [5, 19, 13, 14, 15, 99])
    yvar = "dm"
    xvar = "dm_grace"

    x = region_fits_grace[At(xvar), At(rgis), At("trend"), At(false)]
    y = region_fits_grace[At(yvar), At(rgis), At("trend"), At(false)]
    yscale = endorheic_scale_correction[At("dm"), At(rgis)]
    rmse = sqrt(mean((x .- (y .* yscale)) .^ 2))
    println("rates of loss agree within $(round(rmse, digits=1)) Gt/yr RMSE for the overlapping span")
end

begin
    println("\n------------- Unsustainable glacier freshwater flux --------------------")

    v = "dm"
    p = "trend"
    rgi = 99;
    val0 = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)];
    println("Globally, glaciers lost ice at a net rate of: $(round(Int,val0)) ± $(round(Int,region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)])) Gt/yr")

    v = "dm";
    sf = 1. - endorheic_scale_correction[At(v), At(rgi)];
    val1 = sf*region_fits[varname = At(v), rgi = At(rgi), parameter = At("trend"), error = At(false)]
    println("About $(round(Int, sf*100))% ($(round(Int, val1)) Gt yr-1) of this loss occurred within endorheic basins")
    
    rgi = 98;
    sf = 1. - endorheic_scale_correction[At(v), At(rgi)];
    val2 = sf * region_fits[varname = At(v), rgi = At(rgi), parameter = At("trend"), error = At(false)];
    println("the majority of which ($(round(Int, val2/val1 *100))%; $(round(Int, val2)) Gt yr-1) occurred in High Mountain Asia")
    
    rgi = 12;
    sf = 1. - endorheic_scale_correction[At(v), At(rgi)];
    val3 = sf * region_fits[varname = At(v), rgi = At(rgi), parameter = At("trend"), error = At(false)];
    println("Another $(round(Int, val3/val1 *100))%; ($(round(val3, digits=1)) Gt yr-1) of glacier loss within endorheic basins occurred in North Asia, accounting for $(round(Int,sf*100))% of that region’s loss")


    println("In total, $(round(Int,val0-val1)) Gt yr-1 of the $(round(Int,val0)) Gt yr-1 glacier loss reached the ocean and directly")


    v = "dm";
    p = "acceleration";
    rgi = 99;
    println("Globally, glacier loss accelerated at a rate of $(round(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)], digits=1)) ± $(round(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)], digits=1)) Gt yr-2");

    v = "dm";
    p = "acceleration";
    rgi = 1;
    println("Alaska experiencing the highest accelerated loss at $(round(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)], digits=1)) ± $(round(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)], digits=1)) Gt yr-2");
    
    rgi = 9
    println("followed by the Russian Arctic at $(round(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)], digits=1)) ± $(round(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)], digits=1)) Gt yr-2");

    rgi = 7
    val4 = round(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)], digits=1);
    rgi = 98;
    val5 = round(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)], digits=1);


    println("High Mountain Asia and Svalbard have rates of acceleration of $(round(val5, digits=1)) and $(round(val4, digits=1)) Gt yr-2, respectively")

    rgi = 99;
    v = "dm";
    p = "trend";
    gardner_glambie = region_fits_glambie[varname=At(v), rgi=At(rgi), parameter=At(p), error=At(false)];
    println("Glacier Mass Balance Intercomparison Exercise (GlaMBIE: GlaMBIE Team, 2025) (2000-2024: $(round(Int, region_fits_glambie[varname =At(v), rgi = At(rgi), parameter = At(p), error = At(false)])) ± $(round(Int, region_fits_glambie[varname =At(v), rgi = At(rgi), parameter = At(p), error = At(true)])) Gt yr-1")


    rgi = 5;
    println("largest differences in mean rates between studies are found for Greenland Periphery ($(round(Int, region_fits_glambie[varname =At(v), rgi = At(rgi), parameter = At(p), error = At(false)])) ± $(round(Int, region_fits_glambie[varname =At(v), rgi = At(rgi), parameter = At(p), error = At(true)])) Gt yr-1")

    rgi = 1;
    println("and Alaska ($(round(Int, region_fits_glambie[varname=At(v), rgi=At(rgi), parameter=At(p), error=At(false)]))) ± $(round(Int, region_fits_glambie[varname =At(v), rgi = At(rgi), parameter = At(p), error = At(true)])) Gt yr-1")


    rgi = 99;
    p = "acceleration";
    a1 = region_fits_glambie[varname=At(v), rgi=At(rgi), parameter=At(p), error=At(false)];
    a2 = region_fits_glambie[varname=At(v), rgi=At(rgi), parameter=At(p), error=At(true)];
        println("despite our finding of about  $(round(Int,((a1 - -3.6)/-3.6)*100))% larger rates of accelerated loss ($(round(a1,digits=1)) ± $(round(a2,digits=1)) Gt yr-1")
end

#begin
    println("\n----------------------Glaciers contributions to rivers and oceans--------------------")

    v = "runoff";
    p = "trend";
    rgi = 99;
    
    val0 = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)] + region_fits[varname = At("discharge"), rgi=At(99), parameter= At("trend"), error = At(false)];
    println("Over the study period, glaciers contributed: $(round(Int,val0)) ± $(round(Int,hypot(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)], region_fits[varname = At("discharge"), rgi=At(99), parameter= At("trend"), error = At(true)]))) Gt/yr of freshwater to rivers and oceans")


    v = "discharge";
    p = "trend";
    rgi = 99;
    val1 = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)];
    println("Of this total, $(round(Int,val1/val0*100))% ($(round(Int,val1)) ± $(round(Int,region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)])) Gt/yr) entered the ocean via frontal ablation")


    rgi = 19;
    val1A = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)]
    val1Ae = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)]
    println("Antarctic Islands ($(round(Int,val1A/val1*100))%; $(round(Int,val1A)) ± $(round(Int,val1Ae)) Gt yr⁻¹")

    rgi = [3,4,5,6,7,9];
    val1A = sum(region_fits[varname = At(v), rgi=At([3,4,5,6,7,8,9]), parameter= At(p), error = At(false)], dims=:rgi)[1]
    val1Ae = sum(region_fits[varname = At(v), rgi=At([3,4,5,6,7,8,9]), parameter= At(p), error = At(true)], dims=:rgi)[1]
    println("Arctic Islands including Greenland Periphery ($(round(Int,val1A/val1*100))%; $(round(Int,val1A)) ± $(round(Int,val1Ae)) Gt yr⁻¹")
    
    rgi = 17;
    val1A = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)]
    val1Ae = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)]
    println("Southern Andes ($(round(Int,val1A/val1*100))%; $(round(Int,val1A)) ± $(round(Int,val1Ae)) Gt yr⁻¹")
    
    rgi = 1;
    val1A = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)]
    val1Ae = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)]
    println("and Alaska ($(round(Int,val1A/val1*100))%; $(round(Int,val1A)) ± $(round(Int,val1Ae)) Gt yr⁻¹")


    v = "runoff";
    p = "trend";
    rgi = 99;
    val2 = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)];
    println("The remaining, $(round(Int,val2/val0*100))% ($(round(Int,val2)) ± $(round(Int,region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)])) Gt yr⁻¹) was runoff from surface melt")

    v = "refreeze";
    p = "trend";
    rgi = 99;
    val3 = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)];
    println("that did not refreeze within the firn ($(round(Int,val3)) ± $(round(Int,region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)])) Gt yr⁻¹)")

    v = "runoff"
    sf4 = endorheic_scale_correction[At(v), At(rgi)];
    val4 = val2*sf4;
    println("Of all runoff globally, $(round(Int, sf4*100))% ($(round(Int,val4)) Gt yr⁻¹) reached the ocean, while $(round(Int, (1-sf4)*100))% ($(round(Int,val2*(1-sf4))) Gt yr⁻¹) was retained in endorheic basins")

    rgi = 98;
    v = "runoff";
    p = "trend";
    sf5 = endorheic_scale_correction[At(v), At(rgi)];
    val5 = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)] * (1-sf5);
    println("with HMA accounting for $(round(Int,val5/(val2*(1-sf4))*100))% of that endorheic runoff")


    v = "runoff";
    p = "acceleration";
    rgi = 99;
    val6 = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)];
    println("Glacier runoff increased globally at a rate of $(round(val6, digits=1)) ± $(round(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)], digits=1)) Gt yr⁻² over the past quarter century")


    v = "acc";
    p = "acceleration";
    rgi = 99;
    val7 = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)];
    println(" partially offset by a $(round(val7, digits=1)) ± $(round(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)], digits=1)) Gt yr⁻² increase in accumulation that has occurred ")

    v = "fac";
    p = "trend";
    rgi = 99;
    println("Warming also causes firn densification, which has reduced glacier volume by $(round(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)], digits=1)) ± $(round(region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)], digits=1)) km3 yr-1 ")

    println("Alaska and Arctic Canada North losing $(round(region_fits[varname = At(v), rgi=At(1), parameter= At(p), error = At(false)], digits=1)) ± $(round(region_fits[varname = At(v), rgi=At(1), parameter= At(p), error = At(true)], digits=1)) and $(round(region_fits[varname = At(v), rgi=At(3), parameter= At(p), error = At(false)], digits=1)) ± $(round(region_fits[varname = At(v), rgi=At(3), parameter= At(p), error = At(true)], digits=1)) km3 yr-1, respectively")

    println("and Antarctic glaciers gaining $(round(region_fits[varname = At(v), rgi=At(19), parameter= At(p), error = At(false)], digits=1)) ± $(round(region_fits[varname = At(v), rgi=At(19), parameter= At(p), error = At(true)], digits=1)) km3 yr-1")


    volume_to_mass_conversion_factor = round.(Int, region_fits[varname = At("dm"), parameter= At("trend"), error = At(false)] ./ region_fits[varname = At("dv"), parameter= At("trend"), error = At(false)] .*1000)
    println("Our results indicate a global average volume-to-mass conversion factor of $(volume_to_mass_conversion_factor[rgi = At(99)]) kg m⁻³")

    rgis = setdiff(drgi, [99, 13, 14, 15])
    drgis = dims(volume_to_mass_conversion_factor[At(rgis)], :rgi)
    index_min =argmin(volume_to_mass_conversion_factor[At(rgis)])
    GGA.rginum2label(drgis[index_min])
    index_max =argmax(volume_to_mass_conversion_factor[At(rgis)])
    println("and varies regionally from $(volume_to_mass_conversion_factor[At(rgis)][index_min]) kg m⁻³ in $(GGA.rginum2label(drgis[index_min])) to $(volume_to_mass_conversion_factor[At(rgis)][index_max]) kg m⁻³ in $(GGA.rginum2label(drgis[index_max]))")

    volume_to_mass_conversion_amplitude = round.(Int,region_fits[varname = At("fac"), parameter= At("amplitude"), error = At(false)] ./ region_fits[varname = At("dv"), parameter= At("amplitude"), error = At(false)] .* 100)
    println("Globally, firn air content evolution accounts for $(volume_to_mass_conversion_amplitude[rgi = At(99)])%")

    rgis = setdiff(drgi, [99, 13, 14, 15])
    drgis = dims(volume_to_mass_conversion_amplitude[At(rgis)], :rgi)
    index_min =argmin(volume_to_mass_conversion_amplitude[At(rgis)])
    GGA.rginum2label(drgis[index_min])
    index_max =argmax(volume_to_mass_conversion_amplitude[At(rgis)])
    println("ranges regionally from $(volume_to_mass_conversion_amplitude[At(rgis)][index_min])% in $(GGA.rginum2label(drgis[index_min])) to $(volume_to_mass_conversion_amplitude[At(rgis)][index_max])% in $(GGA.rginum2label(drgis[index_max]))")

    error_in_amp = round.(Int, ((region_fits[varname=At("dv"), parameter=At("amplitude"), error=At(false)] .* 0.85) .- region_fits[varname=At("dm"), parameter=At("amplitude"), error=At(false)])./ region_fits[varname=At("dm"), parameter=At("amplitude"), error=At(false)] * 100)
    rgi = 98
    println("error in amplitude: $(error_in_amp[rgi = At(rgi)])%")

    index_min =argmin(error_in_amp[At(rgis)])
    GGA.rginum2label(drgis[index_min])
    index_max =argmax(error_in_amp[At(rgis)])
    println("ranges regionally from $(error_in_amp[At(rgis)][index_min])% in $(GGA.rginum2label(drgis[index_min])) to $(error_in_amp[At(rgis)][index_max])% in $(GGA.rginum2label(drgis[index_max]))")



    # FIGURE 2 numbers
    println("\n--------------------- FIGURE 2 numbers ---------------------")
    rgi = 99;
    for v = ["acc", "fac", "dm", "discharge", "runoff", "ec", "rain", "refreeze"]
        println("$(v):\n    trend: $(round(Int,region_fits[varname = At(v), rgi=At(rgi), parameter= At("trend"), error = At(false)])) ± $(round(Int,region_fits[varname = At(v), rgi=At(rgi), parameter= At("trend"), error = At(true)])) Gt yr-1\n    acceleration:$(round(region_fits[varname = At(v), rgi=At(rgi), parameter= At("acceleration"), error = At(false)], digits=1)) ± $(round(region_fits[varname = At(v), rgi=At(rgi), parameter= At("acceleration"), error = At(true)], digits=1)) Gt yr-2 ")

        if v == "dm" || v == "runoff"
            println("    endorheic fraction: [$(round((1-endorheic_scale_correction[At(v), At(rgi)])*100, digits=1))%]")
        end
    end

    # check for mass closure
    dm0 = region_fits[varname = At("acc"), rgi=At(rgi), parameter= At("trend"), error = At(false)] - region_fits[varname = At("runoff"), rgi=At(rgi), parameter= At("trend"), error = At(false)] - region_fits[varname = At("discharge"), rgi=At(rgi), parameter= At("trend"), error = At(false)] - region_fits[varname = At("ec"), rgi=At(rgi), parameter= At("trend"), error = At(false)] 
    println("\nThere is a $(round(Int, (dm0 - region_fits[varname = At("dm"), rgi=At(rgi), parameter= At("trend"), error = At(false)]))) Gt yr-1 mass closure error [I ASSUME THAT THIS COMES FROM LINEAR FITTING TO EACH VARIABLE INDIVIDUALLY]")

    v = "net_acc";
    p = "trend";
    rgi = 99;
    println("We find that Earth’s glaciers accumulate a net $(round(Int, region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)])) ± $(round(Int, region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)])) Gt yr⁻¹")

    v = "dm";
    p = "trend";
    rgi = 99;
    println("while an additional $(round(Int, region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)])) ± $(round(Int, region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(true)])) Gt yr⁻¹ of unsustainable loss has")

    v = "gsi";
    p = "trend";
    rgis = [3,4,5,6,7,9];

    val8 = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)];
        println("Globally, glaciers have a GSI of $(round(val8, digits=2)), indicating that the sustainable freshwater flux is augmented by an additional $(round(Int,(1-val8)*100))% from net glacier mass loss ")

    println("$(GGA.rginum2label.(rgis)) —are in the poorest health, with GSIs between $(round.(extrema(region_fits[varname = At(v), rgi = At(rgis), parameter= At(p), error = At(false)]), digits=2))")


    rgi = 4;
    val9 = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)];
    println("The highest GSI of $(round(val9, digits=2)) is observed in Arctic Canada South, implying that net accumulation would need to increase by $(round(Int, (val9-1)*100))%")


    rgi = 3;
    val10 = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)];
    println("Arctic Canada North, with a GSI of $(round(val10, digits=2)), would similarly require a require a $(round(1-val10, digits=2)*100)% increase in net accumulation to achieve balance.")


    rgi = 19;
    val11 = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)];
    println("Antarctic peripheral glaciers have the lowest GSI of $(round(val11, digits=2)), indicating near-sustainable conditions")


    rgi = 98;
    val12 = region_fits[varname = At(v), rgi=At(rgi), parameter= At(p), error = At(false)];
    println("High Mountain Asia has a GSI of $(round(val12, digits=2)), much closer to balance conditions than")


    rgis = [8, 10, 11, 12, 98, 16, 17, 18]
    val13 = region_fits[varname = At(v), rgi=At(rgis), parameter= At(p), error = At(false)];
    println("$(GGA.rginum2label.(rgis)) have GSIs in a narrow range from $(round.(extrema(val13), digits=2))")


    xvar = "dm_grace"
    yvar = "dm"
    rgis =rgis = setdiff(drgis, [5, 19, 13, 14, 15, 99])
    x = region_fits_grace[varname=At(xvar), rgi=At(rgis), parameter=At("trend"), error=At(false)]
    yscale = endorheic_scale_correction[At("dm"), At(rgis)]
    y = eachslice(runs_rgi_fits_grace[varname=At(yvar), rgi=At(rgis), parameter=At("trend")], dims=:run)
    rmse = [sqrt(mean((x .- (y .* yscale)) .^ 2)) for y in y]
    rmse1 = rmse[run = At(ensemble_reference_file)]

    # without endorheic correction
    yscale0 = 1  
    rmse = [sqrt(mean((x .- (y .* yscale0)) .^ 2)) for y in y]
    rmse2 = rmse[run = At(ensemble_reference_file)]

    # assuming a fac correction of 0.85
    yscale_fac = 0.85
    yvar = "dv"
    y = eachslice(runs_rgi_fits_grace[varname=At(yvar), rgi=At(rgis), parameter=At("trend")], dims=:run)
        rmse = [sqrt(mean((x .- (y .* yscale_fac .* yscale)) .^ 2)) for y in y]
    rmse3 = rmse[run = At(ensemble_reference_file)]

    println("RMSE with GRACE: $(round(rmse1, digits=1)) Gt yr⁻¹, RMSE without endorheic correction: $(round(rmse2, digits=1)) Gt yr⁻¹, RMSE without fac correction: $(round(rmse3, digits=1)) Gt yr⁻¹")
end

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!! PLOTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
begin
    # set grace Antarctic and Greenland to NaNs othersise ice sheet will be plotted with glaciers
    regions["dm_grace"][At([5, 19]), :, :] .= NaN;

    # plot altimetry results for RGI regions
    exclude_regions = [13, 14, 15, 99]
    rgi_regions = setdiff(collect(dims(runs_rgi["dm_altim"], :rgi)), exclude_regions)

    variables = ["dm_altim"];
    outfile = joinpath(paths[:figures], "regional_$(join(variables,"_")).png");
    f, region_order, ylims = GGA.plot_multiregion_dvdm(
        regions;
        variables = ["dm_altim"], # last variable is plotted last
        units = "Gt yr⁻¹",
        rgi_regions,
        showlines = false,
        fontsize = 15,
        colormap = :Dark2_4,
        daterange = DateTime(2000, 1, 1):Month(1):DateTime(2025, 1, 1),
        );
    display(f)
    CairoMakie.save(outfile, f);

    # plot altimetry results with GRACE data for RGI regions
    variables = ["dm_grace", "dm"];
    
    exclude_regions = [5, 19, 13, 14, 15, 99] 

    outfile = joinpath(paths[:figures], "regional_$(join(variables,"_")).png");
    f, _, ylims = GGA.plot_multiregion_dvdm(
        regions;
        variables, # last variable is plotted last
        units = "Gt yr⁻¹",
        rgi_regions = setdiff(collect(dims(runs_rgi["dm_altim"], :rgi)), exclude_regions),
        fontsize=15,
        colormap=:Dark2_4,
        region_order = Dim{:rgi}(setdiff(region_order,  exclude_regions)),
        daterange=DateTime(2000, 1, 1):Month(1):DateTime(2025, 1, 1),
    );
    display(f)
    CairoMakie.save(outfile, f);

    # plot altimetry results with gemb data for RGI regions
    variables = ["dm_altim", "dm"]
    outfile = joinpath(paths[:figures], "regional_$(join(variables,"_")).png");
    f, _, ylims = GGA.plot_multiregion_dvdm(
        regions;
        variables, # last variable is plotted last
        units = "Gt yr⁻¹",
        rgi_regions,
        fontsize=15,
        colormap=:Dark2_4,
        region_order,
        all_error_bounds = true,
        daterange=DateTime(2000, 1, 1):Month(1):DateTime(2025, 1, 1),
        numbersinylabel = true,
    );
    display(f)
    CairoMakie.save(outfile, f)

    # plot altimetry results with gemb data for RGI regions
    variables = ["dm_glambie", "dm"]
    outfile = joinpath(paths[:figures], "regional_$(join(variables,"_")).png");
    f, _, ylims = GGA.plot_multiregion_dvdm(
        regions;
        variables, # last variable is plotted last
        units = "Gt yr⁻¹",
        rgi_regions,
        fontsize=15,
        colormap=:Dark2_4,
        region_order,
        all_error_bounds = true,
        daterange=dates_glambie_overlap,
    );
    display(f)
    CairoMakie.save(outfile, f)


    # plot altimetry results with gemb data for RGI regions
    variables = ["dm"];
    outfile = joinpath(paths[:figures], "regional_$(join(variables,"_")).png");
    f, _, ylims = GGA.plot_multiregion_dvdm(
        regions;
        variables, # last variable is plotted last
        units = "Gt yr⁻¹",
        rgi_regions,
        fontsize=15,
        colormap=:Dark2_4,
        region_order,
        daterange=DateTime(1980, 1, 1):Month(1):DateTime(2025, 1, 1),
    );
    display(f)
    CairoMakie.save(outfile, f);

    # GRACE rgi regions    
    colormap = GGA.resample(:Dark2_4, length(region_order))
    colormap = DimArray(colormap.colors, (region_order,))
    
    title = "glacier mass change [Gt/yr]"
    xvar = "dm_grace"
    
    for yvar = ["dm", "dm_altim"]
       
        f = Figure()
     
        ax = Axis(f[1, 1]; 
            ylabel = "this study [Gt/yr]",
            xlabel = "GRACE /FO [Gt/yr]",
            title=title
        )

        for rgi in region_order
            if (rgi == 5) || (rgi == 19)
                continue
            end
            x = region_fits_grace[At(xvar), At(rgi), At("trend"), At(false)]
            y = region_fits_grace[At(yvar), At(rgi), At("trend"), At(false)] 
            yscale = endorheic_scale_correction[At("dm"), At(rgi)]
            scatter!(x, y .* yscale; label=GGA.rginum2label(rgi), markersize=15, color=colormap[At(rgi)])
        end

        for rgi = 98
            scatter!(region_fits_grace[At(xvar), At(rgi), At("trend"), At(false)], region_fits_grace[At(yvar), At(rgi), At("trend"), At(false)]; label=nothing, markersize=15, color=colormap[At(rgi)], strokecolor=:black, strokewidth=3)

            scatter!(-20, -75; label=nothing, markersize=15, color=(colormap[At(rgi)], 0.0), strokecolor=:black, strokewidth=3)

            txt = "including endorheic"
            text!(-17, -75, text=txt, align=(:left, :center))
            
            #text!(grace_dm_gt[At(rgi)], dm_gt[At(false), At(true), At(rgi)], text=txt, align=(:left, :top))
        end

        xlims!(ax, -80, 20)
        ylims!(ax, -80, 20)
        lines!(ax, [-80, 20], [-80, 20], color=:black, linestyle=:dash)

        rgis = setdiff(region_order, [5, 19, 13, 14, 15, 99])
        x = region_fits_grace[At(xvar), At(rgis), At("trend"), At(false)];
        y = region_fits_grace[At(yvar), At(rgis), At("trend"), At(false)];
        yscale = endorheic_scale_correction[At("dm"), At(rgis)]
        rmse = sqrt(mean((x .- (y .* yscale)) .^ 2))

        text!(-75, 15, text="RMSE = $(round(rmse, digits=1)) Gt yr⁻¹", align=(:left, :top))

        leg = Legend(f[1, 2], ax)
        display(f)

        output_dir = joinpath(paths[:figures], "regional_$(xvar)_vs_$(yvar).png")
        save(output_dir, f)
    end

    xvar = "dm_altim"

    for yvar = ["dm"]

        f = Figure()

        ax = Axis(f[1, 1];
            ylabel="modeled [Gt/yr]",
            xlabel="altimetry [Gt/yr]",
            title=title
        )

        for rgi in region_order
            x = region_fits_grace[At(xvar), At(rgi), At("trend"), At(false)]
            y = region_fits_grace[At(yvar), At(rgi), At("trend"), At(false)]
            yscale = 1
            scatter!(x, y .* yscale; label=GGA.rginum2label(rgi), markersize=15, color=colormap[At(rgi)])
        end


        xlims!(ax, -80, 20)
        ylims!(ax, -80, 20)
        lines!(ax, [-80, 20], [-80, 20], color=:black, linestyle=:dash)

        rgis = setdiff(region_order, [5, 19])
        x = region_fits_grace[At(xvar), At(rgis), At("trend"), At(false)]
        y = region_fits_grace[At(yvar), At(rgis), At("trend"), At(false)]
        yscale = 1
        rmse = sqrt(mean((x .- (y .* yscale)) .^ 2))

        text!(-75, 15, text="RMSE = $(round(rmse, digits=1)) Gt/yr", align=(:left, :top))

        leg = Legend(f[1, 2], ax)
        display(f)

        output_dir = joinpath(paths[:figures], "regional_$(xvar)_vs_$(yvar).png")
        save(output_dir, f)
    end


    # GRACE rgi regions for ALL runs
    xvar = "dm_grace"
    for yvar = ["dm", "dm_altim"]
    #yvar = "dm"

        f = Figure()
     
        ax = Axis(f[1, 1]; 
            ylabel = "this study [Gt/yr]",
            xlabel = "GRACE /FO [Gt/yr]",
            title=title
        )

        for rgi in region_order
            if (rgi == 5) || (rgi == 19)
                continue
            end
            x = collect(region_fits_grace[At(xvar), At(rgi), At("trend"), At(false)])
            y = collect(runs_rgi_fits_grace[At(yvar), :, At(rgi), At("trend"), At(false)])
            yscale = endorheic_scale_correction[At("dm"), At(rgi)]
            scatter!(x .* ones(length(y)), y .* yscale; label=GGA.rginum2label(rgi), markersize=15, color=colormap[At(rgi)])
        end


        for rgi in region_order
            if (rgi == 5) || (rgi == 19)
                continue
            end
            x = region_fits_grace[At(xvar), At(rgi), At("trend"), At(false)]
            y = region_fits_grace[At(yvar), At(rgi), At("trend"), At(false)] 
            yscale = endorheic_scale_correction[At("dm"), At(rgi)]
            scatter!(x, y .* yscale; label=nothing, marker = '○', markersize=15, color=:black)
        end

        xlims!(ax, -80, 20)
        ylims!(ax, -80, 20)
        lines!(ax, [-80, 20], [-80, 20], color=:black, linestyle=:dash)

        rgis = setdiff(region_order, [5, 19])
        
        x = region_fits_grace[varname=At(xvar), rgi =At(rgis), parameter=At("trend"), error=At(false)];
        yscale = endorheic_scale_correction[At("dm"), At(rgis)]
        y = eachslice(runs_rgi_fits_grace[varname=At(yvar), rgi =At(rgis), parameter=At("trend")], dims = :run)
        rmse = [sqrt(mean((x .- (y.* yscale)).^2)) for y in y]

        println("ensemble runs with minimum RMSE: dm vs dm_grace")
        display(rmse[minimum(rmse) .== rmse])
     
        
        text!(-75, 15, text="RMSE = $(round(minimum(rmse), digits=1)) to $(round(maximum(rmse), digits=1)) Gt/yr", align=(:left, :top))

        leg = Legend(f[1, 2], ax)
        display(f)

        output_dir = joinpath(paths[:figures], "regional_$(xvar)_vs_$(yvar)_allruns.png")
        save(output_dir, f)
    end
end


# Plot global firn air content anomaly over time and save figure
begin

    center_on_dates = DateTime(2000, 2, 1):Month(1):DateTime(2001, 2, 1)
    runs_rgi = GGA.runs2rgi(path2runs_synthesized);
    runs_rgi = GGA.runs_center!(runs_rgi, center_on_dates);

    # extract reference run and calculate 95% confidence interval
    regions = GGA.runs_ref_and_err(runs_rgi, ensemble_reference_file; error_quantile, error_scaling);

    rgi = 99
    var = "fac"
    date_range = DateTime(2000, 1, 1) .. DateTime(2025, 1, 1)
    
    ts = regions[var][date = date_range, rgi = At(rgi), error = At(false)]
    ts_error = regions[var][date = date_range, rgi = At(rgi), error = At(true)]
    x = GGA.decimalyear.(collect(val(dims(ts, :date))))


    f = GGA._publication_figure(columns=1, rows=1);
    ax = Axis(f[1, 1]; ylabel="volume [km³]", title="Global firn air content anomaly")
    low = ts .- ts_error
    high = ts .+ ts_error

    CairoMakie.band!(ax, x, collect(low), collect(high); color=(:gray, 0.4))

    lines!(ax, x, ts.data)
    xlims!(ax, minimum(x) - 30/365, maximum(x) + 30/365)
    display(f)
    save(joinpath(GGA.pathlocal.figures, "Figure3.png"), f)
end


end