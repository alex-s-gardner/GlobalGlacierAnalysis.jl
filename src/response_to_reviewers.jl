begin
    import CSV
    using DataFrames
    import GlobalGlacierAnalysis as GGA
    using DimensionalData
    using Dates
    using CairoMakie
    using Statistics
    using ProgressMeter

    path2wgms_mb = "/mnt/bylot-r3/data/glacier_mb/WGMS/DOI-WGMS-FoG-2025-02b/data/mass_balance.csv"
    path2wgms_glacier = "/mnt/bylot-r3/data/glacier_mb/WGMS/DOI-WGMS-FoG-2025-02b/data/glacier.csv"
    path2glacier = GGA.pathlocal[Symbol("glacier_individual")]
end;


begin
    wgms_mb = CSV.read(path2wgms_mb, DataFrame)

    # add glacier info to wgms_mb
    wgms_glacier = CSV.read(path2wgms_glacier, DataFrame)
    ia, ib = GGA.intersectindices(collect(wgms_mb.glacier_id), collect(wgms_glacier.id))
    wgms_mb[!, :latitude] = wgms_glacier[ib, :latitude]
    wgms_mb[!, :longitude] = wgms_glacier[ib, :longitude]
    wgms_mb[!, :geometry] = GGA.GI.Point.(wgms_mb.longitude, wgms_mb.latitude)

    # load in study glacier mass balance data
    path2outfile = joinpath(GGA.pathlocal[:project_dir], "Gardner2025_glacier_2deg.nc")
    ds = GGA.netcdf2dimstack(path2outfile)
    glaciers = GGA.GeoDataFrames.read(path2glacier)
end;


begin
    function add_runoff_and_dm!(wgms_mb, ds)
        wgms_mb[!, :runoff] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
        wgms_mb[!, :dm] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
        wgms_mb[!, :smb] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
        wgms_mb[!, :runoff_winter] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
        wgms_mb[!, :dm_winter] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
        wgms_mb[!, :smb_winter] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
        wgms_mb[!, :runoff_summer] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
        wgms_mb[!, :dm_summer] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
        wgms_mb[!, :smb_summer] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
        wgms_mb[!, :geotile_centroid] = Vector{Union{Missing,GGA.GI.Point}}(missing, nrow(wgms_mb))
        wgms_mb[!, :geotile_glacier_area] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
        wgms_mb[!, :geotile_area] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))

        # create an interval grid
        # create global raster for mask coverage calculation
        geotile_width_out = abs(dims(ds, :X)[1] - dims(ds, :X)[2])
        lond = X(GGA.DimensionalData.Dimensions.Lookups.Sampled(-180+geotile_width_out/2:geotile_width_out:180-geotile_width_out/2, sampling=GGA.DimensionalData.Dimensions.Lookups.Intervals(GGA.DimensionalData.Dimensions.Lookups.Center()); metadata=Dict("long_name" => "longitude", "units" => "degrees")))
        latd = Y(GGA.DimensionalData.Dimensions.Lookups.Sampled(-90+geotile_width_out/2:geotile_width_out:90-geotile_width_out/2, sampling=GGA.DimensionalData.Dimensions.Lookups.Intervals(GGA.DimensionalData.Dimensions.Lookups.Center()); metadata=Dict("long_name" => "latitude", "units" => "degrees")))

        # Calculate cell area for each cell [m^2]
        cell_area = GGA.Rasters.cellarea(GGA.Rasters.Raster(zeros(lond, latd); crs=GGA.GeoFormatTypes.EPSG(4326)))

        for r = eachrow(wgms_mb)
            #r = eachrow((wgms_mb))[1523]
            if !ismissing(r.begin_date) && !ismissing(r.end_date)
                r.geotile_glacier_area = ds.area[X=Near(r.longitude), Y=Near(r.latitude)]::Float64
                r.geotile_area = cell_area[X=Near(r.longitude), Y=Near(r.latitude)]::Float64
                runoff = ds.runoff[X=Near(r.longitude), Y=Near(r.latitude), Ti=Touches(r.begin_date, r.end_date)]::DimVector{Float64}
                runoff_cumulative = cumsum(runoff) # m.w.e.
                dm = ds.dm[X=Near(r.longitude), Y=Near(r.latitude), Ti=Touches(r.begin_date, r.end_date)]::DimVector{Float64}
                dm_cumulative = cumsum(dm)
                smb = ds.smb[X=Near(r.longitude), Y=Near(r.latitude), Ti=Touches(r.begin_date, r.end_date)]::DimVector{Float64}
                smb_cumulative = cumsum(smb)

                if length(runoff_cumulative) > 0
                    r.runoff = (runoff_cumulative[end] - runoff_cumulative[1]) / r.geotile_glacier_area / 1000 # m.w.e.
                    r.dm = (dm_cumulative[end] - dm_cumulative[1]) / r.geotile_glacier_area / 1000 # m.w.e.
                    r.smb = (smb_cumulative[end] - smb_cumulative[1]) / r.geotile_glacier_area / 1000 # m.w.e.
                    r.geotile_centroid = GGA.GI.Point(refdims(runoff)[1][1], refdims(runoff)[2][1])


                    if !ismissing(r.midseason_date)
                        r.midseason_date = DateTime(r.midseason_date)

                        runoff = ds.runoff[X=Near(r.longitude), Y=Near(r.latitude), Ti=Touches(r.begin_date, r.midseason_date)]::DimVector{Float64}
                        runoff_cumulative = cumsum(runoff) # m.w.e.
                        dm = ds.dm[X=Near(r.longitude), Y=Near(r.latitude), Ti=Touches(r.begin_date, r.midseason_date)]::DimVector{Float64}
                        dm_cumulative = cumsum(dm)
                        smb = ds.smb[X=Near(r.longitude), Y=Near(r.latitude), Ti=Touches(r.begin_date, r.midseason_date)]::DimVector{Float64}
                        smb_cumulative = cumsum(smb)

                        if length(runoff_cumulative) > 0
                            r.runoff_winter = (runoff_cumulative[end] - runoff_cumulative[1]) / r.geotile_glacier_area / 1000 # m.w.e.
                            r.dm_winter = (dm_cumulative[end] - dm_cumulative[1]) / r.geotile_glacier_area / 1000 # m.w.e.
                            r.smb_winter = (smb_cumulative[end] - smb_cumulative[1]) / r.geotile_glacier_area / 1000 # m.w.e.
                        end
                        runoff = ds.runoff[X=Near(r.longitude), Y=Near(r.latitude), Ti=Touches(r.midseason_date, r.end_date)]::DimVector{Float64}
                        runoff_cumulative = cumsum(runoff) # m.w.e.
                        dm = ds.dm[X=Near(r.longitude), Y=Near(r.latitude), Ti=Touches(r.midseason_date, r.end_date)]::DimVector{Float64}
                        dm_cumulative = cumsum(dm)
                        smb = ds.smb[X=Near(r.longitude), Y=Near(r.latitude), Ti=Touches(r.midseason_date, r.end_date)]::DimVector{Float64}
                        smb_cumulative = cumsum(smb)

                        if length(runoff_cumulative) > 0
                            r.runoff_summer = (runoff_cumulative[end] - runoff_cumulative[1]) / r.geotile_glacier_area / 1000 # m.w.e.
                            r.dm_summer = (dm_cumulative[end] - dm_cumulative[1]) / r.geotile_glacier_area / 1000 # m.w.e.
                            r.smb_summer = (smb_cumulative[end] - smb_cumulative[1]) / r.geotile_glacier_area / 1000 # m.w.e.
                        end

                    end
                end
            end
        end

        wgms_mb[!, :synth_minus_obs] = wgms_mb[:, :smb] .- wgms_mb[:, :annual_balance]
        wgms_mb[!, :synth_minus_obs_winter] = wgms_mb[:, :smb_winter] .- wgms_mb[:, :winter_balance]
        wgms_mb[!, :synth_minus_obs_summer] = wgms_mb[:, :smb_summer] .- wgms_mb[:, :summer_balance]
        wgms_mb[!, :both_valid] = (.!ismissing.(wgms_mb.annual_balance) .& .!ismissing.(wgms_mb.dm)) .& .!isnan.(wgms_mb.synth_minus_obs)
        wgms_mb[!, :both_valid_winter] = (.!ismissing.(wgms_mb.winter_balance) .& .!ismissing.(wgms_mb.dm_winter)) .& .!isnan.(wgms_mb.synth_minus_obs_winter)
        wgms_mb[!, :both_valid_summer] = (.!ismissing.(wgms_mb.summer_balance) .& .!ismissing.(wgms_mb.dm_summer)) .& .!isnan.(wgms_mb.synth_minus_obs_summer)
    end

    add_runoff_and_dm!(wgms_mb, ds)
end;

index = wgms_mb.both_valid .& (wgms_mb.begin_date .> DateTime(2000, 1, 1));

# add mising glacier areas
@showprogress desc = "Adding missing areas..." Threads.@threads for r in eachrow(wgms_mb)[index.&ismissing.(wgms_mb.area)]
    ind = findfirst(GGA.GO.intersects.(glaciers.geom, Ref(r.geometry)))
    if !isnothing(ind)
        r.area = round(glaciers[ind, :Area] * 1.E6)
    end
end;


# write unique validation glaciers to gpkg
begin
    gdf = groupby(wgms_mb[index, :], :glacier_name);
    df2 = combine(gdf, :geometry => first => :geometry, :geometry => length => :count);

    # save points locations
    GGA.GeoDataFrames.write(joinpath(GGA.pathlocal[:project_dir], "wgms_mb_validation_points.gpkg"), df2; crs=GGA.EPSG(4326))
end;

# save with glacier geometry
begin
    rename!(df2, :geometry => :point_location)
    df2[!, :geometry] .= repeat([glaciers.geom[1]], nrow(df2))

    @showprogress desc = "Adding glacier geometry to validation points..." Threads.@threads for r in eachrow(df2)
        ind = findfirst(GGA.GO.intersects.(glaciers.geom, Ref(r.point_location)))
        if !isnothing(ind)
            r.geometry = glaciers.geom[ind]
        end
    end

    # save with glacier geometry
    no_intersection = df2[:, :geometry] .== Ref(glaciers.geom[1])
    GGA.GeoDataFrames.write(joinpath(GGA.pathlocal[:project_dir], "wgms_mb_validation.gpkg"), df2[.!no_intersection, Not(:point_location)]; crs=GGA.EPSG(4326))
end;

# Peyto Glacier investigation
#=
index_peyto = (wgms_mb.glacier_id .== 57) .& (wgms_mb.year .> 1999) 

#index_peyto_rgi = findfirst(GGA.GO.intersects.(glaciers.geom, Ref(wgms_mb.geometry[findlast(index_peyto)])))
#glaciers[index_peyto_rgi, :]

area = ds.area[X=Near(wgms_mb.longitude[findlast(index_peyto)]), Y=Near(wgms_mb.latitude[findlast(index_peyto)])]::Float64
dm = ds.dm[X=Near(wgms_mb.longitude[findlast(index_peyto)]), Y=Near(wgms_mb.latitude[findlast(index_peyto)]), Ti=Touches(DateTime(2000, 1, 1), DateTime(2025, 1, 1))]::DimVector{Float64}
dm_cumulative = cumsum(dm) ./ area/1000

f = Figure();
ax1 = Axis(f[1, 1]);
scatterlines!(wgms_mb[index_peyto, :year], cumsum(wgms_mb[index_peyto, :annual_balance]); label="in situ");
scatterlines!(GGA.decimalyear.(dims(dm_cumulative, :Ti).val), dm_cumulative.data; label="modeled");
f
=#


# print total number of valid annual_balance observations in the WGMS dataset
println("total number of valid annual_balance observations in the WGMS dataset: $(sum(.!ismissing.(wgms_mb.annual_balance)))")

# print total number of valid observations in the WGMS dataset
println("total number of valid observations in the WGMS dataset: $(sum(index))")

# glacier area observed each year
index_area = index .& .!ismissing.(wgms_mb.area);
gdf = groupby(wgms_mb[index_area, :], :year);
df2 = combine(gdf, :area => sum => :area_sum, :annual_balance => std => :annual_balance_std, :synth_minus_obs => std => :synth_minus_obs_std);
wgms_annual_coverage = df2.area_sum ./ sum(ds.area[:]) *100;
println("average annual coverage: $(round(mean(wgms_annual_coverage), digits=2))%")

f = GGA._publication_figure(columns=1, rows=1);
ax1 = CairoMakie.Axis(f[1, 1]; ylabel="global glacier area observed [%]");
lines!(df2.year, wgms_annual_coverage);
fill_between!(ax1, df2.year, wgms_annual_coverage, 0);
xlims!(ax1, minimum(df2.year), maximum(df2.year));
ylims!(ax1, (0, nothing));
f
fname = joinpath(GGA.pathlocal.figures, "wgms_mb_annual_coverage.png");
CairoMakie.save(fname, f)


# calculate variogram of annual balance as a function of distance from observation
begin
    wgms_mb[!, :geometry] = GGA.GI.Point.(wgms_mb.longitude, wgms_mb.latitude)
    mincount = 3
    # calculate std of glaciers as a function of distance from observation

    delta_distance = 10_000 .* cumsum(collect(1:2:16)) 
    ddistance = Dim{:distance}(delta_distance) # meters
    index_variogram = .!ismissing.(wgms_mb.annual_balance) .& (wgms_mb.year .> 1999)
    variogram = zeros(ddistance)

    @showprogress desc = "Calculating variogram of annual balance as a function of distance from observation..." Threads.@threads for buffer_radius in ddistance
        intermediate_year = Vector{Float64}()
        gdf = groupby(wgms_mb[index_variogram, :], :year)

        buffer_ind = findfirst(delta_distance .== buffer_radius)
        if buffer_ind == 1
            buffer_inner = 0
        else
            buffer_inner = delta_distance[buffer_ind-1]
        end

        for df2 in gdf
        #df2 = gdf[10]
            intermediate_pt = Vector{Float64}()

            for r in eachrow(df2)
                # this needs to be a donut shape
                cap1 = GGA.UnitSphericalCap(r.geometry, buffer_radius)
                cap2 = GGA.UnitSphericalCap(r.geometry, buffer_inner)
                polygon1 = GGA.to_latlong_polygon(cap1, 30)
                polygon2 = GGA.to_latlong_polygon(cap2, 30)

                index_within_radius = GGA.GO.intersects.(Ref(polygon1), df2.geometry) .& .!GGA.GO.intersects.(Ref(polygon2), df2.geometry)
                if sum(index_within_radius) > mincount
                    absdiff = abs.(df2.annual_balance[index_within_radius] .- r.annual_balance)
                    mean_absdiff = mean(absdiff[absdiff .!= 0])
                    push!(intermediate_pt, mean_absdiff)
                end
            end
            if length(intermediate_pt) > 0
                push!(intermediate_year, mean(intermediate_pt))
            end
        end
        variogram[distance=At(buffer_radius)] = mean(intermediate_year)
    end
end;

gdf = groupby(wgms_mb[index, :], :geotile_centroid);
df2 = combine(gdf, :synth_minus_obs => std => :synth_minus_obs_std, :synth_minus_obs => length => :nobs, :geotile_area => first => :geotile_area, :synth_minus_obs => mean => :synth_minus_obs_mean);

radius = round(Int,sqrt(mean(df2.geotile_area[index_count]) / (pi)) / 1000)
println("average radius of geotiles with at least 65 observations: $(radius) km")
synth_minus_obs_mad = round(mean(abs.(wgms_mb[index, :synth_minus_obs])), digits=3)
println("mean absolute difference between in situ and this study: $(synth_minus_obs_mad)")

f = GGA._publication_figure(columns=1, rows=1);
ax1 = CairoMakie.Axis(f[1, 1]; ylabel="mean absolute difference [m w.e. yr⁻¹]", xlabel="distance to other observations [km]");
scatterlines!(ax1, ddistance.val ./ 1000, variogram.data);
plot!(ax1, Point(Float64(radius), synth_minus_obs_mad), color=:red);
f
fname = joinpath(GGA.pathlocal.figures, "wgms_mb_variogram.png");
CairoMakie.save(fname, f)

# now take average difference in all geotiles



index_area = index .& .!ismissing.(wgms_mb.area);
gdf = groupby(wgms_mb[index_area, :], :year);
df2 = combine(gdf, :area => sum => :area_sum, :annual_balance => std => :annual_balance_std, :synth_minus_obs => std => :synth_minus_obs_std);
wgms_annual_coverage = df2.area_sum ./ sum(ds.area[:]) * 100;
println("average annual coverage: $(round(mean(wgms_annual_coverage), digits=2))%")



begin
    bound_upper = 3;
    bound_lower = -5;
    bins = bound_lower:0.25:bound_upper
    df = GGA.binstats(wgms_mb[index, :], [:annual_balance, :dm], [bins, bins], :both_valid; missing_bins=true, col_function=sum)
    x = unique(df[:, 1])
    y = unique(df[:, 2])
    density = reshape(df.both_valid_sum, length(x), length(y))

    f = GGA._publication_figure(columns=1, rows=1);
    ax1 = CairoMakie.Axis(f[1, 1]; ylabel="in situ [m w.e. yr⁻¹]", xlabel="this study [m w.e. yr⁻¹]")
    hm = heatmap!(x, y, density)
    Colorbar(f[1, 2], hm)
    lines!(ax1, [bound_lower, bound_upper], [bound_lower, bound_upper], color=:black)
    display(f)
end;


hist(wgms_mb[index, :synth_minus_obs])
hist(wgms_mb[wgms_mb[:, :both_valid_summer], :synth_minus_obs_summer])
hist(wgms_mb[wgms_mb[:, :both_valid_winter], :synth_minus_obs_winter])





# print the nuber of geotiles with at least 1 observation
index_count = df2.nobs .> 0;
println("number of geotiles with at least 1 observation: in $(nrow(df2)) geotiles")

# print number of goetiles with at least 20 observations
index_count = df2.nobs .> 20;
println("number of geotiles with at least 20 observations: in $(sum(index_count)) geotiles, repesenting $(round(Int,sum(df2.nobs[index_count]) ./ sum(df2.nobs) * 100))% of all observations")

# print number of goetiles with at least 50 observations
index_count = df2.nobs .> 65;
println("number of geotiles with at least 65 observations: in $(sum(index_count)) geotiles, repesenting $(round(Int,sum(df2.nobs[index_count]) ./ sum(df2.nobs) * 100))% of all observations")



lines(cumsum(sort(df2.nobs; rev=true))./sum(df2.nobs)*100;  ylabel="% of observations", xlabel="number of geotiles")


index = .!ismissing.(wgms_mb.annual_balance) .& .!isnan.(wgms_mb.annual_balance) .& (wgms_mb.begin_date .< DateTime(2000, 1, 1));

if ismissing(gdf[1])
    gdf = gdf[2:end]
end

for mincount in 5:15
    geotile_stats = DataFrame(geotile_centroid=Vector{Union{Missing,GGA.GI.Point}}(missing, length(gdf)), in_situ_std=Vector{Union{Missing,Float64}}(missing, length(gdf)), nyears=Vector{Union{Missing,Int64}}(missing, length(gdf)))
  
    for (i, k) in enumerate(keys(gdf))
        #(i, k) = first(enumerate(keys(gdf)))

        geotile_stats[i, :geotile_centroid] = getindex(k, 1)
        gdf2 = groupby(gdf[i], :year)
        df2 = combine(gdf2, :annual_balance => std => :in_situ_std, :annual_balance => length => :nobs)

        if any(df2.nobs .>= mincount)
            geotile_stats[i, :in_situ_std] = mean(df2.in_situ_std[df2.nobs.>=mincount])
            geotile_stats[i, :nyears] = sum(df2.nobs.>=mincount)
        end
    end

    intertile_std = mean(geotile_stats.in_situ_std[.!ismissing.(geotile_stats.in_situ_std)])
    println("mincount: $mincount, intertile_std: $intertile_std")
end


#TODO



