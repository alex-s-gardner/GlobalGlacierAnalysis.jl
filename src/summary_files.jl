"""
    summary_files.jl

Generate summary files for glacier mass change analysis.

This script processes glacier elevation change data and glacier surface mass balance model
outputs to create summary files at both geotile and individual glacier levels. It calculates
mass change trends and associated uncertainties, combining data from multiple model runs to
produce robust estimates.

The workflow:
1. Loads and processes glacier elevation change data from multiple model runs
2. Computes mean values and uncertainties across different model configurations
3. Aggregates results by geotile and individual glacier
4. Performs validation checks on global mass change trends
5. Saves results to NetCDF files for further analysis

Key outputs:
- Glacier-level mass change time series with uncertainty estimates
- Geotile-level aggregated statistics
- NetCDF files containing all variables with proper units and metadata
"""

begin
    import GlobalGlacierAnalysis as GGA
    using FileIO
    using DataFrames
    import GeometryOps as GO
    using NCDatasets
    using DimensionalData
    using CairoMakie
    using Dates
    using Statistics
    using ProgressMeter
    using Unitful
    using NonlinearSolve
    using LsqFit
    using GlobalGlacierAnalysis.MyUnits
    paths = GGA.pathlocal

    Unitful.register(MyUnits)

    geotile_summary_file = joinpath(paths[:project_dir], "gardner2025_geotile_summary.nc")
    glacier_summary_file = joinpath(paths[:project_dir], "gardner2025_glacier_summary.nc")

    reference_run = "binned/2deg/glacier_dh_best_nmad5_v01_filled_ac_p2_synthesized.jld2"
    reference_period = (DateTime(2000, 1, 1), DateTime(2000, 12, 31))

    project_id = ["v01"]
    surface_mask = ["glacier", "glacier_rgi7"]
    dem_id = ["best", "cop30_v2"]
    curvature_correct = [false, true]
    amplitude_correct = [true]
    binning_method = ["median", "nmad5", "nmad3"]
    paramater_set = [1, 2, 3, 4]
    binned_folder = ["/mnt/bylot-r3/data/binned/2deg", "/mnt/bylot-r3/data/binned_unfiltered/2deg"]

    path2reference = joinpath(paths[:data_dir], reference_run)   

    param_nt = (;project_id, surface_mask, dem_id, curvature_correct, amplitude_correct, binning_method, paramater_set, binned_folder)
    params = GGA.ntpermutations(param_nt)

    # only include files that exist
    path2runs = String[]
    for param in params
        binned_aligned_file = GGA.binned_aligned_filepath(; param...)
        if isfile(binned_aligned_file)
            push!(path2runs, binned_aligned_file)
        end
    end
    path2runs = replace.(path2runs, "aligned.jld2" => "synthesized.jld2")

    ## calculate individual glacier hypsometry
    geomfile = Dict()
    for sm in surface_mask
        if sm == "glacier"
            geomfile[sm] = joinpath(paths.data_dir, "GlacierOutlines/rgi60/rgi60_Global.gpkg")
        elseif sm == "glacier_rgi7"
            geomfile[sm] = joinpath(paths.data_dir, "GlacierOutlines/RGI2000-v7.0-G-global-fix/rgi70_Global.gpkg")
        end
    end

    glacier_geotile_hyps_fn = Dict()
    for sm in surface_mask
        glacier_geotile_hyps_fn[sm] = replace(geomfile[sm], ".gpkg" => "geotile_hyps.jld2")
    end
end

@time begin #[10 min]
    sm = "glacier"
    binned_synthesized_dv_files = replace.(path2runs,".jld2" => "_gembfit_dv.jld2")
    binned_synthesized_dv_ref= replace("binned/2deg/glacier_dh_best_nmad5_v01_filled_ac_p2_synthesized.jld2",".jld2" => "_gembfit_dv.jld2")
    # only data with a matching surface mask can be downscaled
    binned_synthesized_dv_files = filter(x -> occursin("$(sm)_dh", x), binned_synthesized_dv_files)
    index_ref = findfirst(occursin.(binned_synthesized_dv_ref, binned_synthesized_dv_files));

    path2ref = binned_synthesized_dv_files[index_ref]
    path2files = setdiff(binned_synthesized_dv_files, [binned_synthesized_dv_files[index_ref]]);

    # return a DimensionalArray with the mean and error for all runs relative to the reference run, for all modeled variables 
    geotiles0 = GGA.geotiles_mean_error(path2ref, path2files; reference_period, p = 0.95) # returns variables in units of kg/m2 [90s]

    #TODO: it would be good to save this to disk and a gpkg file of trends and amplitudes currently in synthesis_plots_gis.jl... synthesis_plots_gis.jl should be deprecated at some point

    # load glaciers
    glaciers0 = load(glacier_geotile_hyps_fn[sm], "glaciers")
    glaciers0.area_km2 = sum.(glaciers0.area_km2)

    # glaciers have been split by geotile by model variables have been averaged by geotile groups so this is no longer needed... combine split glaciers
    glaciers = DataFrame()
    rgiids = unique(glaciers0.RGIId)
    glaciers[!, :rgiid] = rgiids
    glaciers[!, :area_km2] .= 0.0 * u"km^2"
    glaciers[!, :geotile] .= ""
    glaciers[!, :rgi] .= 0.0
    glaciers[!, :latitude] .= 0.0
    glaciers[!, :longitude] .= 0.0

    Threads.@threads for r in eachrow(glaciers)
        index = findall(isequal.(Ref(r.rgiid), glaciers0.RGIId))
        _, ind = findmax(glaciers0.area_km2[index])
        ind = index[ind]
        r.area_km2 = sum(glaciers0.area_km2[index])*u"km^2"
        r.geotile = glaciers0[ind, :geotile]
        r.latitude = glaciers0[ind, :CenLat]
        r.longitude = glaciers0[ind, :CenLon]
        r.rgi = parse(Float64, glaciers0[ind, :O1Region])
    end


    # sanity check: 2000-2024 dm trend should be -316 Gt/yr
    begin
        foo = geotiles0
        
        varname = "dm"
        d = DimensionalData.metadata(foo)
        dgeotile = dims(foo, :geotile)
        area = [d["area"][gt] for gt in dgeotile]
        v0 = uconvert.(u"Gt",dropdims(sum(foo[varname = At(varname), error=At(false)] .* area, dims = :geotile), dims = :geotile))
        
        dTi = dims(foo, :Ti)
        decyear = GGA.decimalyear.(dTi)
        decyear .-= mean(decyear)
        
        dates4trend = (DateTime(2000, 1, 1),DateTime(2024, 12, 31))
        fit = curve_fit(GGA.offset_trend_seasonal2, decyear[dates4trend[1]..dates4trend[2]], ustrip.(v0)[dates4trend[1]..dates4trend[2]], GGA.p_offset_trend_seasonal)
        
        p = lines(v0)
        display(p)
        println("trend = $(round(fit.param[2], digits=2)) Gt/yr")
    end

    glacier_out = GGA.geotiles_mean_error_glaciers(glaciers, geotiles0) #[2 min]

    # sanity check: 2000-2024 dm trend should be -316 Gt/yr
    begin
        foo = glacier_out
        
        varname = "dm"
        v0 = dropdims(sum(foo[varname = At(varname), error=At(false)], dims = :rgiid), dims = :rgiid)
        
        dTi = dims(foo, :Ti)
        decyear = GGA.decimalyear.(dTi)
        decyear .-= mean(decyear)
        
        dates4trend = (DateTime(2000, 1, 1),DateTime(2024, 12, 31))
        fit = curve_fit(GGA.offset_trend_seasonal2, decyear[dates4trend[1]..dates4trend[2]], ustrip.(v0)[dates4trend[1]..dates4trend[2]], GGA.p_offset_trend_seasonal)
        
        lines!(v0)
        display(p)
        println("trend = $(round(fit.param[2], digits=2)) Gt/yr")
    end

    foo = glacier_out
    v0 = dropdims(sum(foo[varname = At("dm"), error=At(false)], dims = :rgiid), dims = :rgiid)
    err0 = dropdims(sum(foo[varname = At("dm"), error=At(true)], dims = :rgiid), dims = :rgiid)
    p = lines(v0.- err0)
    lines!(v0.+ err0)
    lines!(v0)
    display(p)

   # From DimArray to DimStack for saveing
    vars = glacier_out[error= At(false)]
    var_error = glacier_out[error= At(true)]
    dvarerror = Dim{:varname}(parent(val(dims(vars, :varname))) .* "_error")
    var_error = set(var_error, :varname => dvarerror)
    da = cat(vars, var_error, dims=(:varname))
    dstack = DimStack(glacier_out[error= At(false)]; layersfrom = :varname)
   
    # save to netcdf
    NCDataset(glacier_summary_file, "c") do ds

        data_dims = dims(dstack)

        # First all dimensions need to be defined
        for dim in data_dims
            dname = string(DimensionalData.name(dim));
            defDim(ds, dname, length(dim))
        end

        # now add the variables
        for dim in data_dims
            dname = string(DimensionalData.name(dim));
            d = defVar(ds, dname, val(val(dim)), (dname,))
            if DateTime <: eltype(dim)
                d.attrib["cf_role"] = "timeseries_id";
            end
        end

        for vaname in keys(dstack)
            v = defVar(ds, "$vaname", ustrip.(parent(dstack[vaname])), string.(DimensionalData.name.(data_dims)))
            v.attrib["units"] = string(Unitful.unit(dstack[vaname][1]))
        end

        # add area_km2 [this was added but not yet tested... you might get an error that needs to be fixed]
        defVar(ds, "glacier_area", glaciers.area_km2, string.(DimensionalData.name.(data_dims[1:1])))
        v.attrib["units"] = "km^2"

        # add latitude and longitude [This is a hack until I can get geometry into a NetCDF]
        defVar(ds, "latitude", glaciers.latitude, string.(DimensionalData.name.(data_dims[1:1])))
        defVar(ds, "longitude", glaciers.longitude, string.(DimensionalData.name.(data_dims[1:1])))

        # add global attributes
        ds.attrib["title"] = "cumulative glacier mass flux outputs from Gardner et al. 2025"
        ds.attrib["version"] = "beta - " * Dates.format(now(), "yyyy-mm-dd")
        ds.attrib["featureType"] = "timeSeries";

        println("glacier output saved to $glacier_summary_file")
    end;
end
