begin
    progress = true

    using ProgressMeter
    using Altim
    import GeoFormatTypes as GFT
    using Proj
    import GeometryOps as GO
    import GeoInterface as GI
    using CairoMakie
    using LibGEOS
    using Rasters, ArchGDAL
    using JLD2
    using FileIO
    #using Altim
    using ArchGDAL
    using GeoDataFrames
    using DataFrames
    #using RangeExtractor
    using GeometryBasics
    using DimensionalData.Lookups
    import DimensionalData as DD
    using Colors
    using BinStatistics
    using LibGEOS # was not able to use GeometryOps for union
    using Tyler
    using Extents
    using StaticArrays
    using Downloads
    using NCDatasets
    using Unitful
    using Dates
    using ColorSchemes
    using Format
    # Load local configuration paths
    paths = Altim.pathlocal

    # exclude rivers at high latitudes and near the dateline when calculating population
    latlim = 85.0

    # path to flux data
    rivers_path = joinpath(paths.data_dir, "rivers/MERIT_Hydro_v07_Basins_v01")
    glacier_rivers_path = joinpath(rivers_path, "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.arrow")

    glacier_summary_file = joinpath("/mnt/bylot-r3/data/project_data", "gardner2025_glacier_summary.nc")
    glacier_summary_gmax_file = replace(glacier_summary_file, ".nc" => "_gmax.gpkg")

    population_file = joinpath("/mnt/bylot-r3/data/project_data", "gardner2025_gmax_buffer_population.nc")
end

begin
    # build a vrt of global high resolution population data
    gpw_dir = "/mnt/bylot-r3/data/sedac-popdensity-yeargrid5yr-v4.11";
    gpw_file = joinpath(gpw_dir, "gpw_v4_population_density_rev11_2020_30_sec_2020.tif");
    if !isfile(gpw_file)
        isdir(gpw_dir) || mkdir(gpw_dir)
        url = "https://data.ghg.center/sedac-popdensity-yeargrid5yr-v4.11/gpw_v4_population_density_rev11_2020_30_sec_2020.tif"
        Downloads.download(url, gpw_file)

        # NOTE: gpw can be petty crude in places like Asia, this could be impoved by merging with the hrsl dataset. 
        # https://dataforgood-fb-data.s3.amazonaws.com/hrsl-cogs/hrsl_general/hrsl_general-latest.vrt
        # The best way to do this would be to resample the hrsl data to 30 arc-second resolution and merge with the gpw data... the hrsl dataset is not global and the native 30 meter resolution is too overkill. 
    end;

    # load rivers
    rivers = GeoDataFrames.read(glacier_rivers_path);

    # exclude rivers at high latitudes and near the dateline
    midlatitude_poly = GI.Polygon(GI.LinearRing(
        [
        (-180.0, -latlim),
        (-180.0, latlim),
        (180.0, latlim),
        (180.0, -latlim),
        (-180.0, -latlim)
    ]
    ));

    rivers[!, :midlatitude] = GO.intersects.(getindex.(collect.(GI.getpoint.(rivers.geometry)), 1), Ref(midlatitude_poly));

    # exclude rivers at high latitudes and near the dateline
    rivers = rivers[rivers[:, :midlatitude], :];


    # load in river glacier runoff data
    glacier_gmax = GeoDataFrames.read(glacier_summary_gmax_file)[:, [:COMID, :gmax_avg, :runoff_max_avg]];
    glacier_gmax[!, :runoff_max_avg] = coalesce.(glacier_gmax.runoff_max_avg, 0.0)
    # add gmax to rivers
    if issorted(glacier_gmax.COMID)
        index = [searchsortedfirst(glacier_gmax[:, :COMID], comid) for comid in rivers[!, :COMID]];
        rivers[!, :gmax_avg] = glacier_gmax[index, :gmax_avg];
        rivers[!, :runoff_max_avg] = glacier_gmax[index, :runoff_max_avg];
    else
        index = [findfirst(==(comid), glacier_gmax.COMID) for comid in rivers[!, :COMID]];
        rivers[!, :gmax_avg] = glacier_gmax[index, :gmax_avg];
        rivers[!, :runoff_max_avg] = glacier_gmax[index, :runoff_max_avg];
    end

    country_polygons = GeoDataFrames.read(paths[:countries]);
    rename!(country_polygons, :SOVEREIGNT => :country);
    rename!(country_polygons, :CONTINENT => :continent);
    country_polygons.geometry = GI.convert.(Ref(LibGEOS), country_polygons.geometry);


    replace!(rivers[!, :country], "United States of America" => "U.S.A.")
    replace!(rivers[!, :country], "United Republic of Tanzania" => "Tanzania")

    replace!(country_polygons[!, :country], "United States of America" => "U.S.A.")
    replace!(country_polygons[!, :country], "United Republic of Tanzania" => "Tanzania")
end;

if !isfile(population_file)
    # load population raster
    gpw_ras0 = Raster(gpw_file)
    #gpw_ras is small enought that if can be read into memory
    gpw_ras0 = Rasters.read(gpw_ras0)
    gpw_ras0 = replace_missing(gpw_ras0, missingval=0)
    # gpw is in units of number of persons per square kilometer, convert to total number of persons
    gpw_ras = (Rasters.cellarea(gpw_ras0) / 1E6) .* gpw_ras0;

    # compute population for each country as a function of maximum glacier fraction and buffer radius
    gmax =[25, 50]# [25, 50]
    buffer_radii = [1_000, 5_000, 10_000, 20_000, 30_000, 40_000, 50_000]
    runoff = [0, 1.0, 2.5, 5.0, 10.0, 100.0, 1000.0] # runoff threshold [m³/s]

    # subset to just Canada for testing 
    #country_polygons = country_polygons[country_polygons[:, :country] .== "Russia", :]

    
    # 94s for all countires, gmax = 50, buffer_radii = 30km
    @time population = Altim.compute_population(gpw_ras, rivers, country_polygons, gmax, runoff, buffer_radii; progress=true)

    name = "population"
    global_attributes = Dict(
        "title" => "population living within a buffer distance of a river with a gmax threshold",
        "version" => "beta - " * Dates.format(now(), "yyyy-mm-dd"),
        )
    fn = Altim.dimarray2netcdf(population, population_file; name, units=nothing, global_attributes)
else
    population = Altim.netcdf2dimarray(population_file; varname=nothing)
end


# low values
gmax_low = 50
buffer_low = 10_000
runoff_low = 10

data = population[gmax=At(gmax_low), buffer=At(buffer_low), runoff=At(runoff_low)];
data = data[data .> 0]
dcountry = dims(data, :country)
dscenario = Dim{:scenario}(["low", "mid", "high"])
data = zeros(dcountry, dscenario)

data[At(val(dcountry)), At("low")] = population[country=At(val(dcountry)), gmax=At(gmax_low), buffer=At(buffer_low), runoff=At(runoff_low)]
println("low: $(round(Int, sum(population[gmax=At(gmax_low), buffer=At(buffer_low), runoff=At(runoff_low)])/1E6)) million")

# mid values
gmax_mid = 50
buffer_mid = 30_000
runoff_mid = 10

data[At(val(dcountry)), At("mid")] = population[country=At(val(dcountry)), gmax=At(gmax_mid), buffer=At(buffer_mid), runoff=At(runoff_mid)]
println("mid: $(round(Int, sum(population[gmax=At(gmax_mid), buffer=At(buffer_mid), runoff=At(runoff_mid)])/1E6)) million")

# high values
gmax_high = 25
buffer_high = 50_000
runoff_high = 1

data[At(val(dcountry)), At("high")] = population[country=At(val(dcountry)), gmax=At(gmax_high), buffer=At(buffer_high), runoff=At(runoff_high)]
println("high: $(round(Int, sum(population[gmax=At(gmax_high), buffer=At(buffer_high), runoff=At(runoff_high)])/1E6)) million")

sort_idx = sortperm(parent(data[scenario=At("mid")]); rev = false)
data = data[sort_idx,:]


begin
    fig = Figure(size = (400, 700))
    xstart = 10
    yticklables = parent(val(dims(data,:country)))
    ncountry = length(data[scenario=At("mid")])

    ax = Axis(
        fig[1,1],
        yticks = (1:ncountry,yticklables),
        xticks= ([1 * 10^i for i = 2:2:8], ["10²", "10⁴", "10⁶", "10⁸"]),
        xscale = log10,
        xminorticks = [10 * 10^i for i = 2:2:8],
        xminorticksvisible = true,
        xgridvisible = false,
        xlabel = "population",
    )
 
    xlims!(ax, low=xstart, high=3E8)
    ylims!(ax, (0.5, ncountry+0.5))

    for (i, d) in enumerate(eachslice(data, dims=:country))
        low = d[scenario = At("low")]
        high = d[scenario = At("high")]
        mid = d[scenario = At("mid")]
    
        scatter!(ax, tuple.(low, i); color = :black, markersize = 5)
        scatter!(ax, tuple.(high, i); color = :black, markersize = 5)
        lines!(ax, [low, high], [i,i]; color = :black)
        scatter!(ax, tuple.(mid, i); color = ColorSchemes.BuPu_3[end], markersize = 20)

        txt_color = :gray51
        if mid > 1E6
            text!(xstart * 1.4, i, text="$(string(round(mid/1E6, digits=1))) [$(string(round(low/1E6, digits=1))), $(string(round(high/1E6, digits=1)))]M", align=(:left, :center), color=txt_color)
        elseif mid > 1E3
            text!(xstart * 1.4, i, text="$(string(round(Int, mid/1E3))) [$(string(round(Int, low/1E3))), $(string(round(Int, high/1E3)))]k", align=(:left, :center), color=txt_color)
        else
            text!(xstart * 1.4, i, text="$(round(Int, mid)) [$(round(Int, low)), $(round(Int, high))]", align=(:left, :center), color=txt_color)
        end
    end

    display(fig)

    outfile = joinpath(paths[:figures], "country_population.png")
    save(outfile, fig)
end

# create a table of results
countries = collect(val(dims(population, :country)));
country2continent = Dict(country => country_polygons.continent[findfirst(country_polygons.country .== country)] for country in countries);

df = DataFrame(country=countries)
df[!, :continent] = getindex.(Ref(country2continent), df[!, :country]);
index = falses(nrow(df))

for gmax in [gmax_low, gmax_mid, gmax_high]
    df[!, "gmax ≥ $gmax"] = parent(population[gmax=At(gmax), buffer=At(buffer), runoff=At(runoff)])
    index .|= df[!, "gmax ≥ $gmax"] .> 0
end
df = df[index,:]

sort!(df, ["continent", "gmax ≥ 50"])

df[!, :geometry] .= [country_polygons.geometry[3]];

for dr in eachrow(df)
    index = findall(==(dr.country), country_polygons.country);
    area = GO.area.(country_polygons.geometry[index]);
    dr.geometry = country_polygons.geometry[index[argmax(area)]] 
end


with_theme(theme_dark()) do
    data = dropdims(sum(population[gmax=At(gmax), runoff=At(runoff)], dims=:country), dims=:country)

    fig = Figure(fontcolor=:green, fontsize = 24)
    
    ax = Axis(fig[1, 1]; 
        backgroundcolor=RGBA{Float32}(0.0f0, 0.0f0, 0.0f0, 0.0f0), 
        yticklabelcolor=(:white, 0.5), 
        xticklabelcolor=(:white, 0.5),
        ylabelcolor=(:white, 0.5),
        xlabelcolor=(:white, 0.5),
    )

    x = collect(val(dims(data, :buffer)))/1E3
    x = vcat(0, x)
    y = vcat(0, collect(data)) * 1E-6
    fill_between!(ax, x, 0, y; color=(:grey, 0.4))
    lines!(ax, x, y; color=(:white, 0.6))
    CairoMakie.xlims!(ax, [0, 50])
    ylims!(low=0)

    ax.ylabel =  "population [million]"
    ax.xlabel = "distance [km] from river with gmax ≥ $gmax%"
   
    outfile = joinpath(paths[:figures], "river_buffer_population.png")

    fig.scene.backgroundcolor = RGBA{Float32}(0.0f0, 0.0f0, 0.0f0, 0.0f0);
    ax.scene.backgroundcolor = RGBA{Float32}(0.0f0, 0.0f0, 0.0f0, 0.0f0);

    save(outfile, fig)
    display(fig)
end


runoff = 10;
buffer = 30_000;
gmax = 50;
pop = population[gmax=At(gmax), buffer=At(buffer), runoff=At(runoff)];
pop = pop[pop.>0];
sort_idx = sortperm(parent(pop); rev=true);
continents = [country_polygons.continent[findfirst(country_polygons.country .== country)] for country in collect(val(dims(pop, :country)))]

total_pop = sum(pop)
println("total population: $(round(Int, total_pop/1E6)) million: gmax ≥ $(gmax), buffer = $(maximum(buffer)) km, runoff > $(runoff) m³/s")

for continent in unique(continents)
    index = continents .== continent
    println("$(continent): $(round(Int, sum(pop[index])/total_pop *100))%")
end


runoff = 1;
buffer = 50_000;
gmax = 25;
pop = population[gmax=At(gmax), buffer=At(buffer), runoff=At(runoff)];
pop = pop[pop.>0];
sort_idx = sortperm(parent(pop); rev=true);
continents = [country_polygons.continent[findfirst(country_polygons.country .== country)] for country in collect(val(dims(pop, :country)))]

total_pop = sum(pop)
println("total population: $(round(Int, total_pop/1E6)) million: gmax ≥ $(gmax), buffer = $(maximum(buffer)) km, runoff > $(runoff) m³/s")
for continent in unique(continents)
    index = continents .== continent
    println("$(continent): $(round(Int, sum(pop[index])/total_pop *100))%")
end


subset = rivers[rivers.gmax_avg .>= gmax, :];
total_river_length = sum(subset.lengthkm);
println("for gmax ≥ $gmax");
println("total river length: $(round(Int, total_river_length/1E3)) thousand km");
println("for buffer = $buffer km");
println("total population: $(round(Int, sum(pop))) million");

total_length = sum(subset.lengthkm)
gdf = groupby(subset, :continent)
println("total river length: $(round(Int, total_length/1E3)) thousand km")
for g in gdf
    continent_length = sum(g.lengthkm)
    if continent_length .> 1E3
        println("$(g.continent[1]): $(round(Int, sum(g.lengthkm)/total_length*100)) % of total river length")
    end
end 
