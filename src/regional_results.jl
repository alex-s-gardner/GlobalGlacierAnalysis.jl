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
begin
    import GlobalGlacierAnalysis as GGA
    using DataFrames
    using FileIO
    using CairoMakie
    using DimensionalData
    using Statistics
    using Dates
    using DataInterpolations
    using LsqFit
    using ColorSchemes
    using GeoDataFrames
    using StatsBase
    using FlexiJoins
    import GeometryOps as GO
    using NCDatasets
    using CSV

    error_quantile=0.95
    error_scaling=1.5

    # to include in uncertainty
    paths = GGA.pathlocal
    reference_run = "binned/2deg/glacier_dh_best_nmad5_v01_filled_ac_p2_synthesized.jld2"
    glacier_flux_path = joinpath(paths[:project_dir], "gardner2025_glacier_summary.nc")

    project_id = "v01"
    surface_mask = ["glacier", "glacier_rgi7"]
    dem_id = ["best", "cop30_v2"]
    curvature_correct = [false, true]
    amplitude_correct = [true]
    binning_method = ["median", "nmad5", "nmad3"]
    fill_param = [1, 2, 3, 4]
    binned_folder = ["/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"]

    discharge_fractional_error = 0.15
    dates2plot = [DateTime(1980, 1, 1),DateTime(2025, 1, 15)] # also dates used for trend fitting
    dates4trend = [DateTime(2000, 3, 1), DateTime(2024, 12, 15)] # use these date for comparing individual model runs to GRACE (only need for rmse anlysis)
    dates_firstdecade = [DateTime(2000, 1, 1), DateTime(2010, 1, 1)] 
    dates_lastdecade = [DateTime(2014, 1, 1), DateTime(2024, 1, 1)] 
    dates_altim_grace_overlap = [DateTime(2002, 5, 15), DateTime(2024, 12, 15)] # use these date for comparing 
    dates_glambie_overlap = [DateTime(2000, 6, 1), DateTime(2023, 1, 1)] # use these date for comparing individual model runs to GRACE (only need for rmse anlysis)
    center_on_dates = DateTime(2002, 5, 1):Month(1):DateTime(2005, 12, 15)

    
    path2reference = joinpath(paths[:data_dir], reference_run)
    glacier_summary_file = GGA.pathlocal[:glacier_summary]

    # path2perglacier = replace(path2reference, ".jld2" => "_perglacier.jld2")
    path2discharge = joinpath(paths[:data_dir], "GlacierOutlines/GlacierDischarge/global_glacier_discharge.jld2")
    path2rgi_regions = paths[:rgi6_regions_shp]

    path2river_flux = joinpath(paths[:river], "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.arrow")

    param_nt = (;project_id = [project_id], surface_mask, dem_id, curvature_correct, amplitude_correct, binning_method, fill_param, binned_folder)
    params = GGA.ntpermutations(param_nt)

    # only include files that exist
    path2runs = String[]
    for param in params
        binned_aligned_file = GGA.binned_filled_filepath(; param...)
        if isfile(binned_aligned_file)
            push!(path2runs, binned_aligned_file)
        end
    end
    path2runs = replace.(path2runs, "aligned.jld2" => "synthesized.jld2")

    ocean_area_km2 = 362.5 * 1E6
end

# pre-process data
begin #[10 min]
    # load discharge for each RGI [<1s]
    discharge_rgi = GGA.discharge_rgi(path2discharge, path2rgi_regions; fractional_error = discharge_fractional_error);

    # load results for all runs for each RGI [18s]
    runs_rgi = GGA.runs2rgi(path2runs);

    # fit trends and acceleration for each model and RGI [33s]
    runs_rgi_fits = GGA.rgi_trends(runs_rgi, discharge_rgi, dates4trend);

    # center data
    runs_rgi = GGA.runs_center!(runs_rgi, center_on_dates);
    GGA.check_for_all_nans(runs_rgi)

    # extract reference run and calculate 95% confidence interval
    regions = GGA.runs_ref_and_err(runs_rgi, path2reference; error_quantile, error_scaling);
    GGA.check_for_all_nans(regions)

    region_fits = GGA.region_fit_ref_and_err(runs_rgi_fits, path2reference; error_quantile, error_scaling, discharge = discharge_rgi)
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
    (dvarname, drgi, dparameters, derr) = dims(region_fits)
    region_fits = cat(region_fits, zeros(Dim{:varname}(["discharge"]), drgi, dparameters, derr), dims=:varname)
    region_fits[At("discharge"), :, At("trend"), At(false)] = discharge_rgi[:, At(false)]
    region_fits[At("discharge"), :, At("trend"), At(true)] = discharge_rgi[:, At(true)]
   
    # deteremine endorheic and non-endorheic glacier runoff [3 min]
    # rivers have only been calculated for rgi6 so this will throw an error if using glacier_rgi7
    if occursin("glacier_rgi7", reference_run)
        error("rivers have only been calculated for rgi6 so this will throw an error if using glacier_rgi7")
    else
        dm_gt_sinktype = GGA.rgi_endorheic(path2river_flux, glacier_summary_file; dates4trend)
        endorheic_scale_correction = 1 .- dm_gt_sinktype[:,:, At(true)] ./ dm_gt_sinktype[:,:, At(false)]
    end

    # compute trends over grace-altim overlap period
    runs_rgi_fits_grace = GGA.rgi_trends(runs_rgi, discharge_rgi, dates_altim_grace_overlap);
    dvarname2 = Dim{:varname}([collect(dvarname)..., "dm_grace"])
    region_fits_grace = zeros(dvarname2, drgi, dparameters, derr)
    region_fits_grace[At(collect(dvarname)), :, : ,:] = GGA.region_fit_ref_and_err(runs_rgi_fits_grace, path2reference; error_quantile, error_scaling, discharge = discharge_rgi)

    # regions 12 and 99 contain NaN so exclude them in trend calculation
    region_fits_grace[At("dm_grace"), At(setdiff(drgi, [12, 99])), : , At(false)] = GGA.rgi_trends(regions["dm_grace"][At(setdiff(drgi, [12, 99])),:, At(false)], dates_altim_grace_overlap);

    # compute trends over the GlaMBIE period
    dvarname2 = Dim{:varname}([collect(dvarname)..., "dm_glambie"])
    region_fits_glambie = zeros(dvarname2, drgi, dparameters, derr)
    runs_rgi_fits_glambie = GGA.rgi_trends(runs_rgi, discharge_rgi, dates_glambie_overlap);
    region_fits_glambie[At(collect(dvarname)), :, :, :] = GGA.region_fit_ref_and_err(runs_rgi_fits_glambie, path2reference; error_quantile, error_scaling, discharge = discharge_rgi)
    region_fits_glambie[At("dm_glambie"), :, : , At(false)] = GGA.rgi_trends(regions["dm_glambie"][:,:, At(false)], dates_glambie_overlap);

    runs_rgi_noaltim = filter(p -> p[1] in setdiff(keys(runs_rgi), ["dm_altim", "dv_altim"]), runs_rgi)
    runs_rgi_fits_firstdecade = GGA.rgi_trends(runs_rgi_noaltim, discharge_rgi, dates_firstdecade);
    region_fits_firstdecade = GGA.region_fit_ref_and_err(runs_rgi_fits_firstdecade, path2reference; error_quantile, error_scaling, discharge = discharge_rgi);
    runs_rgi_fits_lastdecade = GGA.rgi_trends(runs_rgi_noaltim, discharge_rgi, dates_lastdecade);
    region_fits_lastdecade = GGA.region_fit_ref_and_err(runs_rgi_fits_lastdecade, path2reference; error_quantile, error_scaling, discharge = discharge_rgi);
end;

# create table with rgi results and error bars
begin
    exclude_regions = [13, 14, 15]
    varnames = ["acc", "runoff", "dm_altim", "dm", "ec", "dv_altim", "fac", "net_acc", "gsi"]
    varnames = ["dm_altim", "dm", "acc", "runoff", "ec", "fac", "net_acc", "gsi", "discharge"]
    params = ["trend", "acceleration", "amplitude", "phase"]
    regional_results = GGA.error_bar_table(region_fits, varnames, params, setdiff(drgi, exclude_regions));

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
    CSV.write(output_file, regional_results_out)
end

begin
    exclude_regions = [13, 14, 15]
    
    varnames = ["gsi"]
    params = ["trend"]

    println("----------- from $(dates4trend[1]) to $(dates4trend[2]) -----------")
    regional_results = GGA.error_bar_table(region_fits, varnames, params, setdiff(drgi, exclude_regions); digits=2);
    GGA.show_error_bar_table(regional_results; cols2display = 1)

end

# dump summary conclusions to console
begin
    drgi = dims(region_fits, :rgi);
    rgis = setdiff(drgi, [99, 13,14,15])

    varname = "dm";
    
    println("Total glacier contibution to rivers and oceans: $(round(Int,region_fits[At("runoff"), At(99), At("trend"), At(false)] + region_fits[At("discharge"), At(99), At("trend"), At(false)]))) ± $(round(Int,sqrt(region_fits[At("dm"), At(99), At("trend"), At(true)].^2 + region_fits[At("discharge"), At(99), At("trend"), At(true)].^2))) Gt/yr")

    println("Total glacier loss: $(round(Int,region_fits[At("dm"), At(99), At("trend"), At(false)])) ± $(round(Int,region_fits[At("dm"), At(99), At("trend"), At(true)])) Gt/yr")

    println("Total glacier runoff acceleration: $(round(region_fits[At("runoff"), At(99), At("acceleration"), At(false)], digits=1)) ± $(round(region_fits[At("runoff"), At(99), At("acceleration"), At(true)], digits=1)) Gt/yr²")

    println("Total glacier accumulation acceleration: $(round(region_fits[At("acc"), At(99), At("acceleration"), At(false)], digits=1)) ± $(round(region_fits[At("acc"), At(99), At("acceleration"), At(true)], digits=1)) Gt/yr²")

    println("Total glacier loss: $(round(Int, region_fits[At("dm"), At(99), At("trend"), At(false)])) ± $(round(Int, region_fits[At("dm"), At(99), At("trend"), At(true)])) Gt/yr")

    println("Total glacier loss acceleration: $(round(region_fits[At("dm"), At(99), At("acceleration"), At(false)], digits=1)) ± $(round(region_fits[At("dm"), At(99), At("acceleration"), At(true)], digits=1)) Gt/yr²")
   
    rgi_subset = [1, 5, 3, 98, 2, 17]
    fact_of_total = round(Int, sum(region_fits[At(varname),At(rgi_subset),At("trend"),At(false)]) ./ sum(region_fits[At(varname),At(rgis),At("trend"),At(false)]) *100)
    println("$(join(GGA.rginum2label.(rgi_subset), ", ")) comprise $(fact_of_total)% of total mass change")

    varname = "dm";
    sig_index = region_fits[At(varname),:,At("acceleration"),At(false)] .<= -region_fits[At(varname),:,At("acceleration"),At(true)];
    println("Regions with significant acceleration in $(varname):")
    for rgi in drgi[sig_index]
        println("   $(GGA.rgi2label[GGA.rginum2txt[rgi]])")
    end
    println("")

    sig_index = region_fits[At(varname), :, At("acceleration"), At(false)] .>= region_fits[At(varname), :, At("acceleration"), At(true)];
    println("Regions with significant decceleration in $(varname):")
    for rgi in drgi[sig_index]
        println("   $(GGA.rgi2label[GGA.rginum2txt[rgi]])")
    end
    println("")

    # fraction of glacier loss that is endorheic
    for varname in ["dm", "runoff"]
        println("-----------------------------------$varname-----------------------------------------")
        dm_endhoric = round(Int, (1 - endorheic_scale_correction[At(varname), At(99)]) .* region_fits[At(varname), At(99), At("trend"), At(false)])
        dm_ocean = round(Int, (endorheic_scale_correction[At(varname), At(99)]) .* region_fits[At(varname), At(99), At("trend"), At(false)])
        dm_endhoric_frac = round((1 - endorheic_scale_correction[At(varname), At(99)])*100, digits=1)
        println("$(100 - dm_endhoric_frac)% ($(dm_ocean) Gt/yr, $(round((dm_ocean / ocean_area_km2* 1E6), digits=2)) mm/yr SLR) of glacier $varname is ocean terminating")
        println("$(dm_endhoric_frac)% ($(dm_endhoric) Gt/yr) of glacier $varname is endhoric")

        hma_frac_of_total_endorheic = dm_gt_sinktype[At(varname), At(98), At(true)] / dm_gt_sinktype[At(varname), At(99), At(true)];
        hma_endorheic = round(Int, (1 - endorheic_scale_correction[At(varname), At(98)]) .* region_fits[At(varname), At(98), At("trend"), At(false)])
        println("HMA comprises $(round(Int,hma_frac_of_total_endorheic*100))% ($hma_endorheic Gt/yr) of total endorheic $varname")
        
        na_frac_of_total_endorheic = dm_gt_sinktype[At(varname), At(12), At(true)] / dm_gt_sinktype[At(varname), At(99), At(true)];
        na_endorheic = round(dm_gt_sinktype[At(varname), At(12), At(true)], digits = 1)
        println("$(round(Int,(1 - endorheic_scale_correction[At(varname), At(12)])*100))% ($na_endorheic Gt/yr) of North Asia is comprised endorheic $varname")
        println("North Asia comprises $(round(Int,na_frac_of_total_endorheic*100))% ($na_endorheic Gt/yr) of total endorheic $varname")

        dm_total = round(Int,region_fits[At(varname), At(99), At("trend"), At(false)]);
        dm_ocean = round(Int,dm_total * endorheic_scale_correction[At(varname), At(99)]);
        println("Of the total glacier $varname of $(dm_total) Gt/yr, $(dm_ocean) Gt/yr flows to the ocean")
    end

    total_runoff = round(Int, region_fits[At("runoff"), At(99), At("trend"), At(false)])
    total_runoff_err = round(Int, region_fits[At("runoff"), At(99), At("trend"), At(true)])
    total_discharge = round(Int, region_fits[At("discharge"), At(99), At("trend"), At(false)])
    total_discharge_err = round(Int, region_fits[At("discharge"), At(99), At("trend"), At(true)])

    println("Total fresh water flux: $(total_runoff + total_discharge) ± $(round(Int, hypot(total_runoff_err, total_discharge_err))) Gt/yr")
   
    println("$(round(Int, total_discharge/(total_runoff + total_discharge)*100))% ($(total_discharge) ± $(total_discharge_err)) Gt/yr came from discharge")
   
    ant_discharge = round(Int, region_fits[At("discharge"), At(19), At("trend"), At(false)])
    ant_discharge_err = round(Int, region_fits[At("discharge"), At(19), At("trend"), At(true)])
    ant_frac = round(Int, ant_discharge / total_discharge * 100)
    println("$(ant_frac)% ($(ant_discharge) ± $(ant_discharge_err) Gt/yr) of all dischage occured in the Antarctica")

    arctic_island_discharge = round(Int,sum(region_fits[At("discharge"), 2..9, At("trend"), At(false)]))
    arctic_island_discharge_err = round(Int, arctic_island_discharge * discharge_fractional_error)
    arctic_island_frac = round(Int, arctic_island_discharge / total_discharge * 100)
    println("$(arctic_island_frac)% ($(arctic_island_discharge) ± $(arctic_island_discharge_err) Gt/yr) of all dischage occured in the Arctic islands ")

    patagonia_discharge = round(Int,sum(region_fits[At("discharge"), At(17), At("trend"), At(false)]))
    patagonia_discharge_err = round(Int, patagonia_discharge * discharge_fractional_error)
    patagonia_frac = round(Int,patagonia_discharge / total_discharge * 100)
    println("$(patagonia_frac)% ($(patagonia_discharge) ± $(patagonia_discharge_err) Gt/yr) of all dischage occured in Patagonia")

    alaska_discharge = round(Int,sum(region_fits[At("discharge"), At(1), At("trend"), At(false)]))
    alaska_discharge_err = round(Int, alaska_discharge * discharge_fractional_error)
    alaska_frac = round(Int, alaska_discharge / total_discharge * 100)
    println("$(alaska_frac)% ($(alaska_discharge) ± $(alaska_discharge_err) Gt/yr) of all dischage occured in Alaska")


    println("$(round(Int, total_runoff/(total_runoff + total_discharge)*100))% ($(total_runoff) ± $(total_runoff_err)) Gt/yr came from runoff")



    rgis = setdiff(drgi, [5, 19, 13, 14, 15, 99])
    for yvar = ["dm", "dm_altim"]
    #yvar = "dm"
    xvar = "dm_grace"
        println("----------------------------------------------------------------------------")
        println("                  $xvar vs $yvar for rgi regions")

        x = region_fits_grace[At(xvar), At(rgis), At("trend"), At(false)];
        y = region_fits_grace[At(yvar), At(rgis), At("trend"), At(false)];
        yscale = endorheic_scale_correction[At("dm"), At(rgis)]
        rmse = sqrt(mean((x .- (y .* yscale)) .^ 2))
        println("RMSE = $(round(rmse, digits=1)) Gt/yr")

        rmse_incl_endorheic = sqrt(mean((x .- y) .^ 2))
        println("RMSE including endorheic = $(round(rmse_incl_endorheic, digits=1)) Gt/yr")

        rmse_static_fac = sqrt(mean((x .- (y .* yscale) .* 0.85) .^ 2))
        println("RMSE using static fac = $(round(rmse_static_fac, digits=1)) Gt/yr")
        println("----------------------------------------------------------------------------\n")
    end

    println("\n----------------------------------------------------------------------------")
    println("                  Firn Air Content")

    println("Total glacier fac trend: $(round(region_fits[At("fac"), At(99), At("trend"), At(false)], digits=1)) ± $(round(region_fits[At("fac"), At(99), At("trend"), At(true)], digits=1)) km3/yr")
    println("Total glacier fac trend Alaska: $(round(region_fits[At("fac"), At(1), At("trend"), At(false)], digits=1)) ± $(round(region_fits[At("fac"), At(1), At("trend"), At(true)], digits=1)) km3/yr")
    println("Total glacier fac trend Canada Arctic North: $(round(region_fits[At("fac"), At(3), At("trend"), At(false)], digits=1)) ± $(round(region_fits[At("fac"), At(3), At("trend"), At(true)], digits=1)) km3/yr")
    println("Total glacier fac trend Antarctic: $(round(region_fits[At("fac"), At(19), At("trend"), At(false)], digits=1)) ± $(round(region_fits[At("fac"), At(19), At("trend"), At(true)], digits=1)) km3/yr")

    
    rgis = setdiff(drgi, [13, 14, 15])
    for rgi in rgis
        println("Volume to mass conversion $(GGA.rginum2label(rgi)): $(round(Int, region_fits[At("dm"), At(rgi), At("trend"), At(false)] / (region_fits[At("dv"), At(rgi), At("trend"), At(false)]) * 1000)) kg m⁻³")
    end

     for rgi in rgis
        println("seasonal volume change from FAC $(GGA.rginum2label(rgi)): $(GGA.rginum2label(rgi)): $(round(Int,region_fits[At("fac"), At(rgi), At("amplitude"), At(false)]/(region_fits[At("dv"), At(rgi), At("amplitude"), At(false)]) * 100))%")
    end


    # compute amplitude and trend difference if using static fac 
    amp_static_fac_frac = round.((region_fits[At("dv"), :, At("amplitude"), At(false)] * 0.85) ./ region_fits[At("dm"), :, At("amplitude"), At(false)] .- 1, digits=2)
    println()
    println("Amplitude would be `x` times larger if static fac was applied")
    show(stdout, "text/plain", amp_static_fac_frac)
    println("\n----------------------------------------------------------------------------\n")

    println("\n----------------------------------------------------------------------------")
    println("                  rates for GlamBIE overlap period (2000-2024)")

   

    for rgi in rgis
        println("$(GGA.rginum2label(rgi)): $(round(Int, region_fits_glambie[At("dm"), At(rgi), At("trend"), At(false)])) ± $(round(Int, region_fits_glambie[At("dm"), At(rgi), At("trend"), At(true)])) Gt/yr")
    end

    println("\n                  acceleration for GlamBIE overlap period (2000-2024)")

    rgi = 99
    println("$(GGA.rginum2label(rgi)): $(round(region_fits_glambie[At("dm"), At(rgi), At("acceleration"), At(false)], digits=1)) ± $(round(region_fits_glambie[At("dm"), At(rgi), At("acceleration"), At(true)], digits=1)) Gt/yr²")
    println("----------------------------------------------------------------------------\n")


    println("net accumulation $(GGA.rginum2label(99)): $(round(Int, region_fits[At("net_acc"), At(99), At("trend"), At(false)])) ± $(round(Int, region_fits[At("net_acc"), At(99), At("trend"), At(true)])) Gt/yr")
    println("mass loss $(GGA.rginum2label(99)): $(round(Int, region_fits[At("dm"), At(99), At("trend"), At(false)])) ± $(round(Int, region_fits[At("dm"), At(99), At("trend"), At(true)])) Gt/yr")

end

# create plots
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
        
        x = region_fits_grace[At(xvar), At(rgis), At("trend"), At(false)];
        yscale = endorheic_scale_correction[At("dm"), At(rgis)]
        y = eachslice(runs_rgi_fits_grace[At(yvar), :, At(rgis), At("trend"), At(false)], dims = :run)
        rmse = [sqrt(mean((x .- (y.* yscale)).^2)) for y in y]
     
        
        text!(-75, 15, text="RMSE = $(round(minimum(rmse), digits=1)) to $(round(maximum(rmse), digits=1)) Gt/yr", align=(:left, :top))

        leg = Legend(f[1, 2], ax)
        display(f)

        output_dir = joinpath(paths[:figures], "regional_$(xvar)_vs_$(yvar)_allruns.png")
        save(output_dir, f)
    end
end