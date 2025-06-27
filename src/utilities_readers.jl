"""
    glacier_discharge(; datadir=pathlocal[:data_dir]) -> DataFrame

Load and combine glacier discharge data from multiple sources.

# Arguments
- `datadir`: Base directory containing glacier data files (default: pathlocal[:data_dir])

# Returns
- DataFrame containing glacier discharge information with columns:
  - `latitude`: Glacier center latitude
  - `longitude`: Glacier center longitude
  - `discharge_gtyr`: Discharge rate in gigatons per year
  - `discharge_err_gtyr`: Uncertainty in discharge rate
  - `frontal_ablation_gtyr`: Frontal ablation rate in gigatons per year

# Description
Combines discharge data from Kochtitzky (2022) for the Northern Hemisphere and 
Fuerst (2023) for Patagonia. Matches Patagonian glaciers with RGI database entries
to obtain coordinates, and standardizes data format across sources.
"""
function glacier_discharge()
    
    # Kochtitzky NH discharge and terminus retreate 
    nothern_hemisphere = CSV.read(pathlocal[:discharge_nh], DataFrame; header=14, skipto=16)

    df = DataFrame(Shapefile.Table(pathlocal[:rgi6_southern_andes]))

    patagonia = CSV.read(pathlocal[:discharge_npi], DataFrame; header=25, skipto=27)
    patagonia = vcat(patagonia, CSV.read(pathlocal[:discharge_spi], DataFrame; header=25, skipto=27))

    # find matching glaciers in RGI 
    patagonia[!, :"Latitude"] .= NaN
    patagonia[!, :"Longitude"] .= NaN
    df.Name[ismissing.(df.Name)] .= "junk"
    patagonia.Names = replace.(patagonia.Names, "Témpano" => "Tempano")
    patagonia.Names[Base.contains.(patagonia.Names, "Grey")] .= "Grey + Dickson"
    patagonia.Names[Base.contains.(patagonia.Names, "Upsala")] .= "Upsala + Cono"
    patagonia.Names[Base.contains.(patagonia.Names, "O'Higgins")] .= "OHiggins"

    for r in eachrow(patagonia)
        index = findfirst((r.Names .== df.Name) .| (r.Names .== df.RGIId))
        if isnothing(index)
            println("cold not find match for: $(r.Names)")

            continue
        end
        r.Latitude = df.CenLat[index]
        r.Longitude = df.CenLon[index]
    end

    df = DataFrame()
    df[!, :latitude] = vcat(nothern_hemisphere.lat, patagonia.Latitude)
    df[!, :longitude] = vcat(nothern_hemisphere.lon, patagonia.Longitude)
    df[!, :discharge_gtyr] = vcat(nothern_hemisphere[:, "2010_2020_mean_discharge_gt_per_year"], (patagonia[:, "Calving FG"] .* δice / 1000))
    df[!, :discharge_err_gtyr] = vcat(nothern_hemisphere[:, "2010_2020_mean_flux_err_gt"], (patagonia[:, "Unc. Calving FG"] .* δice / 1000))
    df[!, :frontal_ablation_gtyr] = vcat(nothern_hemisphere[:, "Frontal_ablation_2010_to_2020_gt_per_yr_mean"], (patagonia[:, "Frontal ablation Minowa (2000-2019)"] .* δice / 1000))
    df[!, :frontal_ablation_gtyr] = vcat(nothern_hemisphere[:, "Frontal_ablation_2010_to_2020_gt_per_yr_mean"], (patagonia[:, "Unc. Frontal ablation Minowa (2000-2019)"] .* δice / 1000))

    return df
end

"""
    read_zemp2019(; datadir=setpaths().zemp_2019) -> NamedTuple

Read and process glacier mass balance data from Zemp et al. 2019.

# Arguments
- `datadir`: Directory containing Zemp 2019 data files (default: from setpaths())

# Returns
- NamedTuple containing:
  - `dm_gt`: DimArray of cumulative mass change in Gt by RGI region and date
  - `err_gt`: DimArray of uncertainty values in Gt by RGI region and date
  - `all`: Dictionary of raw data by RGI region

# Description
Reads CSV files containing regional glacier mass balance data from Zemp et al. 2019,
extracts RGI region identifiers from filenames, and organizes the data into DimArrays
with dimensions for RGI regions and dates (1950-2016).
"""
function read_zemp2019(; datadir=setpaths().zemp_2019)
    # Read Zemp 2019 data
    fn_startswith = "Zemp_etal_results_region_"
    fn_endswith = ".csv"
    files = allfiles(datadir; fn_startswith, fn_endswith)

    all = Dict()
    for file in files
        # pull RGI number form file name
        foo = file[end-10:(end-7)]
        s = findfirst('_', foo)
        e = findlast('_', foo)
        rgi = "rgi$(foo[(s+1):(e-1)])"
        data = DataFrame((CSV.File(file; header=28, skipto=29, stripwhitespace=true)))
        push!(all, rgi => data)
    end

    date = DateTime.(1950:2016)
    rgi = "rgi" .* string.(1:19)

    drgi = Dim{:rgi}(rgi)
    ddate = Dim{:date}(date)

    dm_gt = DimArray(fill(NaN, (length(drgi), length(ddate))), (drgi, ddate))
    err_gt = DimArray(fill(NaN, (length(drgi), length(ddate))), (drgi, ddate))

    for r in rgi
        d = DateTime.(all[r].Year)
        dm_gt[At(r), At(d)] = cumsum(all[r].INT_Gt)
        err_gt[At(r), At(d)] = all[r].sig_Total_Gt
    end

    return (; dm_gt, err_gt, all)
end


"""
    read_marzeion2020(; datadir=setpaths().marzeion_2020) -> NamedTuple

Read and process glacier mass change data from Marzeion et al. 2020.

# Arguments
- `datadir`: Directory containing Marzeion 2020 data files (default: from setpaths())

# Returns
- NamedTuple containing:
  - `dm_gt`: DimArray of mass change in Gt with dimensions for RGI regions, dates, 
             climate models, glacier models, and scenarios

# Description
Reads NetCDF data containing global glacier mass change projections from Marzeion et al. 2020,
and organizes the data into a multidimensional DimArray with appropriate dimensions.
"""
function read_marzeion2020(; datadir=setpaths().marzeion_2020)

    foo = Dataset(datadir)
    rgi = "rgi" .* string.(Int.(foo["Region"]))
    date = DateTime.(foo["Time"])
    climate_model = getindex.(collect(foo["Climate_Model"].attrib), 2)
    glacier_model = getindex.(collect(foo["Glacier_Model"].attrib), 2)
    scenario = getindex.(collect(foo["Scenario"].attrib), 2)

    drgi = Dim{:rgi}(rgi)
    ddate = Dim{:date}(date)
    dclimate_model = Dim{:climate_model}(climate_model)
    dglacier_model = Dim{:glacier_model}(glacier_model)
    dscenario = Dim{:scenario}(scenario)

    foo = collect(foo["Mass"].var)
    dm_gt = DimArray(foo, (drgi, ddate, dclimate_model, dglacier_model, dscenario))

    return (; dm_gt)
end

"""
    read_marzeion2012(; datadir=setpaths().marzeion_2012) -> NamedTuple

Read and process glacier mass change data from Marzeion et al. 2012.

# Arguments
- `datadir`: Directory containing Marzeion 2012 data files (default: from setpaths())

# Returns
- NamedTuple containing:
  - `dm_gt`: DimArray of mass change in Gt with dimensions for RGI regions, dates, 
             climate models, and scenarios
  - `err_gt`: DimArray of mass change errors in Gt with the same dimensions

# Description
Reads MATLAB data containing global glacier mass change projections from Marzeion et al. 2012,
and organizes the data into multidimensional DimArrays with appropriate dimensions.
"""
function read_marzeion2012(; datadir=setpaths().marzeion_2012)
    foo = matread(datadir)

    scenario = collect(keys(foo))

    # for s in scenario
    s = first(scenario)
    date = vec(DateTime.(foo[s]["yr"]))
    climate_model = vec(foo[s]["model"])
    rgi = vec("rgi" .* string.(Int.(foo[s]["rgiid"])))
    gt2sle = foo[s]["gt2sle"] .* -1

    drgi = Dim{:rgi}(rgi)
    ddate = Dim{:date}(date)
    dclimate_model = Dim{:climate_model}(climate_model)
    dscenario = Dim{:scenario}(scenario)


    dm_gt = fill(NaN, (length(rgi), length(date), length(climate_model), length(scenario)))
    err_gt = fill(NaN, (length(rgi), length(date), length(climate_model), length(scenario)))

    dm_gt = DimArray(dm_gt, (drgi, ddate, dclimate_model, dscenario))
    err_gt = DimArray(err_gt, (drgi, ddate, dclimate_model, dscenario))

    for (i, s) in enumerate(scenario)
        climate_model0 = vec(foo[s]["model"])
        foo1 = foo[s]["mb"] ./ gt2sle
        foo1 = permutedims(foo1, (2, 1, 3))
        dm_gt[:, :, At(climate_model0), i] = foo1

        foo1 = foo[s]["mb_err"] ./ gt2sle
        foo1 = permutedims(foo1, (2, 1, 3))
        err_gt[:, :, At(climate_model0), i] = foo1
    end

    foo = (; dm_gt, err_gt)

    return foo
end


"""
    read_hock2019(; datadir=setpaths().hock_2019) -> NamedTuple

Load and organize glacier mass change projections from Hock et al. 2019.

# Arguments
- `datadir`: Path to the Hock 2019 dataset file (default: from setpaths())

# Returns
- NamedTuple containing `dm_gt`: DimArray of glacier volume change with dimensions for 
  RGI regions, dates, climate models, glacier models, and scenarios

# Description
Reads NetCDF data containing global glacier mass change projections from Hock et al. 2019,
and organizes the data into a multidimensional DimArray with appropriate dimensions.
The function handles the complex combination of scenarios, glacier models, and climate models.
"""
function read_hock2019(; datadir=setpaths().hock_2019)
    #datadir= setpaths().hock_2019
    foo = Dataset(datadir)
    rgi = "rgi" .* string.(Int.(foo["region"]))
    date = DateTime.(foo["time"])

    # from metadata
    scenario = ["A1B", "RCP26", "RCP45", "RCP60", "RCP85"]
    glacier_model = ["BM", "GieOer", "GloGEM", "GloGEMebal", "HYOGA2", "RHB", "VASlangen"]
    climate_model = ["ACCESS1-0", "BCC-CSM1-1", "BNU-ESM", "CCSM4", "CGCM3-1(T63)", "CNRM-CM3", "CNRM-CM5", "CSIRO-Mk3-0", "CSIRO-Mk3-6-0", "CanESM2", "ECHAM5-MPI-OM", "GFDL-CM2-0", "GFDL-CM2-1", "GFDL-CM3", "GFDL-ESM2G", "GFDL-ESM2M", "GISS-E2-R", "HadGEM2-AO", "HadGEM2-CC", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-LR", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC-ESM-CHEM", "MIROC5", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-CGCM3", "NorESM1-M", "NorESM1-ME", "PCM", "UKMO-HadCM3"]

    drgi = Dim{:rgi}(rgi)
    ddate = Dim{:date}(date)
    dclimate_model = Dim{:climate_model}(climate_model)
    dglacier_model = Dim{:glacier_model}(glacier_model)
    dscenario = Dim{:scenario}(scenario)

    coord = (drgi, ddate, dclimate_model, dglacier_model, dscenario)
    dm_gt = DimArray(fill(NaN, length.(coord)), coord)

    for i in 1:length(dclimate_model)
        for j in 1:length(dglacier_model)
            for k in 1:length(dscenario)
                idx = (foo["scenario"] .== k) .& (foo["glaciermodel"] .== j) .& (foo["forcingmodel"] .== i)
                f = foo["volume"].var[:, :, findall(idx)]
                if isempty(f)
                    continue
                else
                    dm_gt[:, :, i, j, k][:] = f
                end
            end
        end
    end

    return (; dm_gt)
end

"""
    read_ipccar6(; datadir=setpaths()[:ipcc_ar6], start2007=false)

Read IPCC AR6 glacier mass change data from CSV files.

# Arguments
- `datadir`: Directory containing IPCC AR6 data files (default: from setpaths())
- `start2007`: Whether to read data starting from 2007 (default: false)

# Returns
A NamedTuple containing:
- `dm_gt`: DimArray of glacier mass change values in Gt by RGI region, date, and scenario
- `err_gt`: DimArray of corresponding error values

# Description
Reads glacier mass change projections from IPCC AR6 Figure 9.21 data files.
The function processes both the main data files and corresponding error files,
organizing them into multidimensional arrays indexed by RGI region, date, and scenario.

rsync -r devon:/mnt/bylot-r3/data/binned/2deg/figures /Users/gardnera/Research/20_01_GlobalGlacierChange/version 2/
"""


function read_ipccar6(; datadir=setpaths()[:ipcc_ar6], start2007=false)
    # Data from IPCC AR6 Figure 9.21
    fn_endswith = ".csv"
    files = allfiles(datadir; fn_endswith)

    f2007 = contains.(files, Ref("2007"))
    if start2007
        files = files[f2007]
    else
        files = files[.!f2007]
    end

    err_files = contains.(files, Ref("error"))
    files_err = files[err_files]
    files = files[.!err_files]
    all = Dict()

    if start2007
        scenarios = [f[end-13:end-9] for f in files]
    else
        scenarios = [f[end-8:end-4] for f in files]
    end

    dscenario = Dim{:scenario}(scenarios)

    data = DataFrame((CSV.File(files[1]; header=1, skipto=2, stripwhitespace=true)))
    ridx = contains.(names(data), Ref("RGI"))

    rgi = "rgi" .* [n[5:end] for n in names(data)[ridx]]
    drgi = Dim{:rgi}(rgi)

    date = DateTime.(data[:, 1])
    ddate = Dim{:date}(date)

    dm_gt = DimArray(fill(NaN, length(drgi), length(ddate), length(dscenario)), (drgi, ddate, dscenario))
    err_gt = copy(dm_gt)

    for k in eachindex(files)
        # pull RGI number form file name
        data = DataFrame((CSV.File(files[k]; header=1, skipto=2, stripwhitespace=true)))
        for j = 1:19
            dm_gt[j, :, k] = data[:, j+1]
        end

        data = DataFrame((CSV.File(files_err[k]; header=1, skipto=2, stripwhitespace=true)))
        for j = 1:19
            err_gt[j, :, k] = data[:, j+1]
        end
    end

    foo = (; dm_gt, err_gt)

    return foo
end

"""
    wimberly2024(; filename = pathlocal[:wimberly_2024])

Loads and processes glacier model data from Wimberly et al. 2024.

# Arguments
- `filename`: Path to the CSV file containing the data (default: pathlocal[:wimberly_2024])

# Returns
- A DimArray containing glacier model runoff data in units of [km^3] organized by date, 
  GCM, SSP, basin, and model type. The function parses the combined GCM_SSP_Basin column 
  into separate dimensions and structures the data from three glacier models 
  (GloGEM, PyGEM, OGGM) into a DimArray.
"""
function wimberly2024(; filename=setpaths().wimberly_2024)

    df = CSV.read(filename, DataFrame)

    # split the GCM_SSP_Basin column into three columns
    df = transform(df, [:GCM_SSP_Basin] => ByRow(x -> split(x, "_")) => [:GCM, :SSP, :Basin])

    ddate = Ti(sort!(unique(df.Date)))
    dgcm = Dim{:GCM}(unique(df.GCM))
    dssp = Dim{:SSP}(unique(df.SSP))
    dbasin = Dim{:Basin}(unique(df.Basin))
    dmodel = Dim{:Model}(["GloGEM", "PyGEM", "OGGM"])

    gdf = groupby(df, [:GCM, :SSP, :Basin])

    da = zeros(ddate, dgcm, dssp, dbasin, dmodel)

    for gcm in dgcm
        for ssp in dssp
            for basin in dbasin
                for model in dmodel
                    da[At(gdf[(gcm, ssp, basin)][:, :Date]), At(gcm), At(ssp), At(basin), At(model)] = gdf[(gcm, ssp, basin)][:, model]
                end
            end
        end
    end

    return da
end


"""
    read_grace_rgi(; path2file=setpaths()[:grace_rgi])

Read GRACE (Gravity Recovery and Climate Experiment) data for RGI (Randolph Glacier Inventory) regions.

This function loads GRACE data from a MATLAB file and converts region codes from letter-based 
abbreviations to standard RGI numerical identifiers (e.g., "ALA" to "rgi1").

# Arguments
- `path2file`: Path to the GRACE RGI data file (MATLAB format)

# Returns
- A dictionary mapping RGI region identifiers to their corresponding GRACE data

# Note
The function converts region codes from the original letter-based abbreviations (e.g., "ALA" for Alaska)
to the standard RGI numerical format (e.g., "rgi1").
"""
function read_grace_rgi(;path2file=setpaths()[:grace_rgi])

    grace0 = matread(path2file)

    # rgi leter to digit mapping
    old2new = Dict(
        "PAT" => "rgi18", 
        "GRE" => "rgi5", 
        "NAS" => "rgi10", 
        "ICE" => "rgi6", 
        "TRP" => "rgi16",
        "SAW" => "rgi14",
        "SAE" => "rgi15",
        "CDN" => "rgi3",
        "CAS" => "rgi13",
        "AND" => "rgi17",
        "CDS" => "rgi4",
        "ANT" => "rgi19",
        "CEU" => "rgi11",
        "ALA" => "rgi1",
        "SVB" => "rgi7",
        "WNA" => "rgi2", 
        "NEZ" => "rgi18",
        "RAI" => "rgi9",
        "SCA" => "rgi8"
    )


    grace = Dict(any(keys(old2new).== key) ? (old2new[key]) => val : (key) => val for (key, val) in grace0)
    return grace
end

"""
    grace_masschange(; path2file=setpaths()[:grace_rgi])

Load and organize GRACE (Gravity Recovery and Climate Experiment) mass change data for glacier regions.

This function reads GRACE data for RGI (Randolph Glacier Inventory) regions and organizes it into a
dimensional array structure with dimensions for region, date, and error bounds.

# Arguments
- `path2file`: Path to the GRACE RGI data file (MATLAB format). Defaults to the path from `setpaths()[:grace_rgi]`.

# Returns
- A DimensionalData array with dimensions for RGI regions (1-19, plus 98 for HMA and 99 for Global),
  dates, and error bounds (false=value, true=error).

# Note
The error values are multiplied by 2 to represent 95% confidence intervals.
"""
function grace_masschange(; path2file=setpaths()[:grace_rgi])
    grace_raw = read_grace_rgi(; path2file)

    # organize GRACE grace data into a DimArray
    drgi = Dim{:rgi}(collect([1:19..., 98:99...]))
    ddate = Dim{:date}(vec(datenum2date.(grace_raw["rgi1"]["dM_gt_mdl_fill_date"])))
    derror = Dim{:error}([false, true])

    grace = fill(NaN, drgi, ddate, derror)

    for rgi in drgi
        #rgi = drgi[1]
        if rgi == 98
            rgi_id = "HMA"
        elseif rgi == 99
            rgi_id = "Global"
        else
            rgi_id = "rgi$(rgi)"
        end

        if haskey(grace_raw, rgi_id)
            grace[At(rgi), :, At(false)] = vec(grace_raw[rgi_id]["dM_gt_mdl_fill"])

            ## multiply by 2 to get 95% confidence interval
            grace[At(rgi), :, At(true)] = vec(grace_raw[rgi_id]["dM_sigma_gt_mdl_fill"]) * 2
        end
    end

    return grace
end

"""
    glambie2024(; path2glambie=paths = setpaths()[:glambie_2024]) -> DimArray

Load and process the GlaMBIE 2024 glacier mass balance dataset.

# Keywords
- `path2glambie`: Path to the GlaMBIE 2024 CSV file (default: from setpaths())

# Returns
- A DimensionalData.DimArray containing cumulative glacier mass balance in Gt
  with dimensions for RGI region (1-19, 98, 99), date (2000-2024), and error flag

# Notes
- Dates are shifted by up to 6 months from original data to align across hemispheres
- Missing values are linearly interpolated
- Region 98 (HMA) is calculated as the sum of regions 13-15
"""
function glambie2024(; path2glambie=paths = setpaths()[:glambie_2024])
    df = CSV.read(path2glambie, DataFrame)

    drgi = Dim{:rgi}(vcat(collect(1:19), 98, 99))
    derror = Dim{:error}([false, true])
    @warn "GlaMBIE dates are shifted by up to 6 months from the original data to align dates across hemispheres"
    ddate = Dim{:date}(collect(decimalyear2datetime.(2000:2024)); name="GlaMBIE 2024 cumulative glacier mass balance [GLAMBIE]")

    glambie_Gt = fill(NaN, drgi, ddate, derror)

    for (i, rgi) in enumerate(drgi)
        if rgi == 98
            continue
        elseif rgi == 99
            i = i - 1
        end
        scol = (i - 1) * 3 + 1
        date1 = round.(df[:, scol])
        valid = .!ismissing.(date1)
        date1 = collect(decimalyear2datetime.(date1[valid]))

        glambie_Gt[At(rgi), At(date1), At(false)] = coalesce.(df[valid, scol+1], NaN)
        glambie_Gt[At(rgi), At(date1), At(true)] = coalesce.(df[valid, scol+2], NaN)
    end


    glambie_Gt[At(98), :, :] .= 0
    for rgi in 13:15
        glambie_Gt[At(98), :, At(false)] .+= glambie_Gt[At(rgi), :, At(false)]
        glambie_Gt[At(98), :, At(true)] .+= glambie_Gt[At(rgi), :, At(true)] .^ 2
    end

    glambie_Gt[At(98), :, At(true)] = sqrt.(glambie_Gt[At(98), :, At(true)])

    @warn "GlaMBIE data is linearly interpolated over missing date values"
    for rgi in drgi
        for err_flag in derror
            valid = .!isnan.(glambie_Gt[At(rgi), :, At(err_flag)])
            if !all(valid)
                fit = curve_fit(model1_trend, collect(decimalyear.(ddate[valid])), collect(glambie_Gt[At(rgi), :, At(err_flag)][valid]), p_trend)
                glambie_Gt[At(rgi), .!valid, At(err_flag)] = model1_trend(decimalyear.(collect(ddate[.!valid])), fit.param)
            end
        end
    end

    return (glambie_Gt)
end
