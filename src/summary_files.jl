# summary_files.jl
#
# Generate summary files for glacier mass change analysis.
#
# This script processes glacier elevation change data and glacier surface mass balance model
# outputs to create summary files at both geotile and individual glacier levels. It calculates
# mass change trends and associated uncertainties, combining data from multiple model runs to
# produce robust estimates.
#
# The workflow:
# 1. Loads and processes glacier elevation change data from multiple model runs
# 2. Computes mean values and uncertainties across different model configurations
# 3. Aggregates results by geotile and individual glacier
# 4. Performs validation checks on global mass change trends
# 5. Saves results to NetCDF files for further analysis
#
# Key outputs:
# - Glacier-level mass change time series with uncertainty estimates
# - Geotile-level aggregated statistics
# - NetCDF files containing all variables with proper units and metadata

#function glacier_summary_files(binned_synthesized_dv_files, binned_synthesized_dv_file_ref; error_quantile = 0.95, error_scaling = 2.0, reference_period = (DateTime(2000, 1, 1), DateTime(2000, 12, 31)), surface_mask = "glacier")
     import GlobalGlacierAnalysis as GGA
     using Unitful
     using DataFrames
     using DimensionalData
     using NCDatasets

     Unitful.register(GGA.MyUnits)
   


    reference_period = (DateTime(2000, 1, 1), DateTime(2000, 12, 31))
    surface_mask = "glacier"
    error_quantile = 0.95
    error_scaling = 2.0

#function glacier_summary_file(binned_synthesized_dv_files, binned_synthesized_dv_file_ref; error_quantile = 0.95, error_scaling = 2.0, reference_period = (DateTime(2000, 1, 1), DateTime(2000, 12, 31)), surface_mask = "glacier")

    path2ref = binned_synthesized_dv_file_ref
    path2files = setdiff(binned_synthesized_dv_files, [binned_synthesized_dv_file_ref])
    paths = GGA.pathlocal
    if length(path2files) == length(binned_synthesized_dv_files)
        error("Reference run not found in binned_synthesized_dv_files")
    end

    # return a DimensionalArray with the mean and error for all runs relative to the reference run, for all modeled variables 
    geotiles0 = GGA.geotiles_mean_error(path2ref, path2files; error_quantile, error_scaling, reference_period) # returns variables in units of kg/m2 [90s]

    #TODO: it would be good to save this to disk and a gpkg file of trends and amplitudes currently in synthesis_plots_gis.jl... synthesis_plots_gis.jl should be deprecated at some point
    glacier_geotile_hyps_fn = replace(GGA.pathlocal[Symbol("$(surface_mask)_individual")], ".gpkg" => "geotile_hyps.jld2")
    glaciers0 = GGA.load(glacier_geotile_hyps_fn, "glacier")

    # load glaciers
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

    geotiles = GGA._geotile_load_align(; surface_mask, geotile_width)
    geotiles[!, :rgi] .= 0
    
    for rgi in 1:19
        geotiles[geotiles[:,"rgi$rgi"].>0, :rgi] .= rgi
    end

    Threads.@threads for r in eachrow(glaciers)
        index = findall(isequal.(Ref(r.rgiid), glaciers0.RGIId))
        _, ind = findmax(glaciers0.area_km2[index])
        ind = index[ind]
        r.area_km2 = sum(glaciers0.area_km2[index])*u"km^2"
        r.geotile = glaciers0[ind, :geotile]
        r.latitude = glaciers0[ind, :CenLat]
        r.longitude = glaciers0[ind, :CenLon]
        r.rgi = geotiles[findfirst(isequal(r.geotile), geotiles.id), :rgi]
    end
    
    vars2downscale = nothing
    
    drgiid = Dim{:rgiid}(glaciers.rgiid)
   (dvarname, dgeotile, dTi, derror) = dims(geotiles0)
    
    if !isnothing(vars2downscale)
        dvarname = Dim{:varname}(vars2downscale)
    end

    glacier_out = GGA.geotiles_mean_error_glaciers(glaciers, geotiles0) #[2 min]
   
   # From DimArray to DimStack for saveing
    lvars = glacier_out[error= At(false)]
    var_error = glacier_out[error= At(true)]
    dvarerror = Dim{:varname}(parent(val(dims(vars, :varname))) .* "_error")
    var_error = set(var_error, :varname => dvarerror)
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

    return dstack
end
