using Altim, GeoTiles, DataFrames, Extents, Rasters, FastGeoProjections, Proj
using Shapefile, GeoInterface, GeometryOps, Plots, Arrow, BinStatistics, Statistics, NearestNeighbors
using Interpolations, ArchGDAL, CSV, Dates, MLJ, MLJLinearModels
import GeoFormatTypes as GFT;
import GeoDataFrames as GDF

force_remake = false;
region0 = [:RGI01, :RGI02, :RGI03, :RGI04, :RGI06, :RGI07, :RGI09, :RGI17, :RGI98, :RGI18, :RGI11]
region0 = [:RGI98]

ref_dem_grid_res = 100;
# hyps binning method
function central(x)
    # filter input
    dev = abs.(x .- median(x));
    mdev = median(dev);
    valid = dev .< (mdev * 1.4826 * 5)
    mean(x[valid])
end
gedi_regions = [:RGI02, :RGI98]

begin
    fname = "data/icesat_campaigns.csv";
    csv_reader = CSV.File(fname)
    icesat_campaigns = DataFrame(csv_reader)

    d = Dates.Date[]
    ucycle = Int16[]
    for (cnt, r) in enumerate(eachrow(icesat_campaigns))
        push!(d, r.start_date)
        push!(d, r.end_date)
        push!(ucycle, cnt)
        push!(ucycle, cnt)
    end
    d = decimalyear.(d)
    icesat_campaigns = linear_interpolation(Interpolations.deduplicate_knots!(d), ucycle)
end



product_paths =  project_paths(; project_id=:v01);

for region in region0
    println("region: $region")
    data_folder = "/mnt/devon-r2/data/dhdt_sampling/$region"

    region_ext, epsg = region_extent(region)
    if region == :RGI01
        region_ext = Extent(X=region_ext.X, Y=(55.0, region_ext.Y[2]))
    elseif region == :RGI02
        region_ext = Extent(X=region_ext.X, Y=(region_ext.Y[1], 55.0))
    end
    
    if any(region .== gedi_regions)
        mission0 = [:icesat, :icesat2, :gedi, :hugonnet]
    else
        mission0 = [:icesat, :icesat2, :hugonnet];
    end

    dem_file = joinpath(data_folder, "$(region)_hugonnet_all_ice.arrow")
    altim_files = joinpath.(Ref(data_folder), "$(Ref(region))_$(mission0)_all.arrow")

    if any(!(isfile.(altim_files)) || !isfile(dem_file)) || force_remake
        if !isdir(data_folder)
            mkdir(data_folder)
        end

        local_paths = setpaths()

        # dhdt files
        dhdt_dir = "2000-01-01_2020-01-01"
        f = joinpath(local_paths.hugonnet_v1_stacks, dhdt_dir)
        f = Altim.allfiles(f, subfolders=true, fn_endswith="_dhdt.tif")
        f = Raster.(f)

        # determine raster extents in lat lon
        f_exts = Extent{(:Y, :X)}[]
        f_epsg = EPSG{1}[]
        for k in f
            ext0 = extent(k)
            crs0 = Rasters.crs(k)
            epsg = GFT.convert(EPSG, Proj.CRS(crs0))
            trans = FastGeoProjections.Transformation(epsg, EPSG(4326))
            corner_x = [ext0.X[1], ext0.X[1], ext0.X[2], ext0.X[2]]
            corner_y = [ext0.Y[1], ext0.Y[2], ext0.Y[2], ext0.Y[1]]
            lat, lon = trans(corner_x, corner_y)

            lat_minmax = extrema(lat)
            lon_minmax = extrema(lon)

            push!(f_epsg, epsg)
            push!(f_exts, Extent(Y=lat_minmax, X=lon_minmax))
        end
        idx = Extents.intersects.(Ref(region_ext), f_exts)
        f_exts = f_exts[idx]
        f_epsg = f_epsg[idx]
        f = f[idx]
    end

    for mission in mission0
        println("   mission: $mission")
        altim_file = joinpath.(data_folder, "$(region)_$(mission)_all.arrow")

        if !isfile(altim_file) || force_remake
            if mission == :hugonnet
                altim = GeoTiles.readall(product_paths[mission].geotile; extent=region_ext, suffix=".arrow")
                select!.(altim, :height)
                rename!.(altim, :height => :dhdt)
            else
                df = GeoTiles.listtiles(product_paths[mission].geotile; extent=region_ext, suffix=".arrow")

                println("       extract ref dhdt at point locations")
                for d in eachrow(df)
                    outfile = GeoTile.path2tile(product_paths[mission].geotile, DataFrames.metadata(gt, "geotile_id"), ".dhdt_ref")

                    if !isfile(outfile) || force_remake
                        t1 = time();
                        # find intersecting rasters
                        inds = findall(Extents.overlaps.(Ref(d.extent), f_exts))
                        if isempty(inds)
                            # if there are not overlapping Rasters then skip
                            continue
                        end
                        gt = GeoTiles.read(d.path2file)
                        # initialize geotile dataframe
                        dhdt = DataFrame(:dhdt => [fill(NaN32, length(h)) for h in gt.height]) 
                        metadata!(dhdt, "geotile_id", DataFrames.metadata(gt, "geotile_id"))
                        
                        #=
                        # THIS IS NOT REALLY NEEDED AND JUST ADD COMPLEXITY TO THE DATA AND CODE
                        if mission == :icesat
                            dhdt[!, :cycle] .= Ref(Int16[])
                            for (i, r) in enumerate(eachrow(gt))
                                dhdt[i, :cycle] = round.(icesat_campaigns(decimalyear.(r.datetime)))
                            end
                            dhdt[!, :rgt] = [g.granule_info.info.rgt for g in eachrow(gt)]
                        elseif mission == :gedi
                            dhdt[!, :cycle] .= Int16(0);
                            dhdt[!, :rgt] = [g.granule_info.info.track for g in eachrow(gt)]
                        else
                            dhdt[!, :cycle] = [Int16.(g.granule_info.info.cycle) for g in eachrow(gt)]
                            dhdt[!, :rgt] = [g.granule_info.info.rgt for g in eachrow(gt)]
                        end

                        if mission == :icesat2
                            dhdt[!, :strong_beam] = [g.strong_beam for g in eachrow(gt)]
                        end
                        =#

                        for i in inds
                            trans = FastGeoProjections.Transformation(EPSG(4326), f_epsg[i])

                            for (j, r) in enumerate(eachrow(gt))
                                if !isempty(r.latitude) 
                                    x, y = trans(r.latitude, r.longitude)
                                    pnts = collect((p[1], p[2]) for p in zip(x,y))
                                    a = collect(extract(f[i], pnts))
                                    a = coalesce.(getindex.(a, 2), NaN32)
                                    
                                    if !(typeof(a) <: Vector)
                                        a = [a]
                                    end

                                    replace!(a, f[i].missingval => NaN32)
                                    ind = isnan.(dhdt[j, :dhdt])
                                    dhdt[j, :dhdt][ind] = a[ind]
                                end
                            end
                        end

                        GeoTiles.save(product_paths[mission].geotile, ".dhdt_ref", dhdt; filetype=:arrow)

                        t2 = time();
                        println("           $(DataFrames.metadata(gt, "geotile_id")): $(round(t2-t1, digits = 2))s")
                    end
                    
                    altim = GeoTiles.readall(product_paths[mission].geotile, suffix = ".dhdt_ref");
                end
            end

            # add heights, datetime, height_ref, lat, longitude, latitude, landice
            altim0 = DataFrame[];
            for gt in altim
                id = DataFrames.metadata(gt, "geotile_id")
                fpath = GeoTiles.path2tile(product_paths[mission].geotile, id, ".arrow")
                fpath_h = GeoTiles.path2tile(product_paths[mission].geotile, id, ".cop30_v2")
                fpath_masks = GeoTiles.path2tile(product_paths[mission].geotile, id, ".masks")
                if isfile(fpath_h) .& isfile(fpath_masks) .& isfile(fpath)
                    h = GeoTiles.read(fpath_h)
                    m = GeoTiles.read(fpath_masks)
                    a = GeoTiles.read(fpath)
                    if (length(h.height) == nrow(gt)) && (length(m.landice) == nrow(gt))
                        gt[!, :height] = a[:, :height]
                        gt[!, :datetime] = a[:, :datetime]
                        gt[!, :height_reference] = h[:, :height]
                        gt[!, :latitude] = h[:, :latitude]
                        gt[!, :longitude] = h[:, :longitude]
                        gt[!, :landice] = m[:, :landice]
                        #@time gt[!, :altim_valid] = [isempty(r.height) ? Bool[] : (abs.(r.height .- r.height_reference) .< 200) for r in eachrow(gt)]

                        # filter data to only include data that intersects landice
                        ind = map(any, gt.landice);
                        gt0 = gt[ind,:];
                        gt0[!, :altim_valid] = [isempty(r.height) ? Bool[] : (abs.(r.height .- r.height_reference) .< 200) for r in eachrow(gt0)]

                        gt0[!, :dhdt_count] =  [sum(.!isnan.(k)) for k in gt0.dhdt]

                        push!(altim0, gt0)
                    else
                        # these tiles typcically fall outside of "ice domain"
                        @warn("has different number of columns: $fpath")
                    end
                else
                    if !isfile(fpath_h)
                        @warn("$(id) is missing reference heights")
                    end
                    if !isfile(fpath_h)
                            @warn("$(id) is missing fpath_masks")
                    end

                    if !isfile(fpath)
                        @warn("$(id) is missing altimetry file")
                    end
                end        
            end

            # hcat altim sampling and remove rows without dhdt data
            altim = reduce(vcat, altim0)
            altim = Altim.dataframe_reduce(altim)
            idx = Altim.within.(Ref(region_ext), altim.longitude, altim.latitude)
            altim = altim[idx, :];
            Arrow.write(altim_file, altim::DataFrame);
        end
    end

    # compute elevation change from static dh/dt field
    if !isfile(dem_file) || force_remake
        println("       extract all ref dhdt")
        p = setpaths();
        shp = allfiles(p.RGI_dissolved; subfolders=true, fn_endswith=".shp")
        polygons = Shapefile.Handle.(shp)
        poly_ext = GeoInterface.bbox.(polygons)

        # initialize
        dem = (height=Vector{Int16}[], dhdt=Vector{Float32}[], latitude=Vector{Float32}[], longitude=Vector{Float32}[])

        for (i, k) in enumerate(f)
            t1 = time()
            ind = findall(Extents.intersects.(Ref(f_exts[i]), poly_ext))
            A = Raster(zeros(UInt16, k.dims[1], k.dims[2]); crs=f_epsg[i]);
            M = Raster(zeros(Bool, k.dims[1], k.dims[2]); crs=f_epsg[i])

            for j in ind
                # reprojection takes up most of the time
                poly = GeometryOps.reproject(polygons[j]; target_crs=f_epsg[i])
                A = Rasters.rasterize!(count, A, poly.geom; threaded=false, shape=:polygon, verbose = false, progress = false)
                M = M .| (A .> 0)
            end

            if !any(M)
                continue
            end

            geometry = DimPoints(k)[M]
            x = round.(getindex.(geometry, 1), digits=4)
            y = round.(getindex.(geometry, 2), digits=4)

            trans = FastGeoProjections.Transformation(f_epsg[i], EPSG(4326))
            lat, lon = trans(x,y)

            push!(dem.dhdt, k[M])
            push!(dem.latitude, lat)
            push!(dem.longitude, lon)
            # now get elevations
            h = itslive_extract(lat, lon, [:h]; always_xy=false)

            push!(dem.height,h.h)
            t2 = time()
            fn = splitpath(f[i].metadata["filepath"])[end]
            println("$(fn): $(round(t2-t1, digits = 2))s")
        end

        dem = (dhdt=reduce(vcat, dem[:dhdt]), height=reduce(vcat,dem[:height]), lat=reduce(vcat,dem[:lat]), lon=reduce(vcat,dem[:lon]))

        notvalid = abs.(dem.dhdt) .> 100
        dem.dhdt[notvalid] .= NaN32

        dem = DataFrame(dem)
        # remove duplicate points from overlaping projections

        # remove points that are closer than 75m from one another [VERY VERY SLOW]
        trans = FastGeoProjections.Transformation(EPSG(4326), EPSG(3413))
        x, y = trans(dem.lat, dem.lon)
        pts = vcat(x', y')
        kdtree = KDTree(pts; leafsize=10)
        idx = inrange(kdtree, pts, 75, true)

        idx = unique(getindex.(idx, 1))
        dem = dem[idx,:]
        idx = Altim.within.(Ref(region_ext), dem.lon, dem.lat)
        dem = dem[idx, :]
        # save file
        Arrow.write(dem_file, dem::DataFrame);
    end

    dem_file = joinpath(data_folder, "$(region)_hugonnet_all_ice.arrow")
    dem = DataFrame(Arrow.Table(dem_file))
    for mission in mission0
        # make dh/dt vs elevation plots
        altim_file = joinpath.(data_folder, "$(region)_$(mission)_all.arrow");
        altim = DataFrame(Arrow.Table(altim_file));
        altim[!, :decyear] = decimalyear.(altim.datetime)


        # add subcycles
        if (mission == :icesat)
            altim[!, :subcycle] = Int16.(round.(icesat_campaigns(altim.decyear)))
        else
            subcycle_dt = (91 / 365.25) / 3
            date_min, date_max = extrema(altim.decyear)
            subcycle_dates = date_min:subcycle_dt:date_max+subcycle_dt
            subcycle_dates = repeat(subcycle_dates, inner=2)
            subcycle_dates = subcycle_dates[2:end-1]
            subcycle = repeat(Int16.(round.((1:length(subcycle_dates)/2))), inner=2)
            subcycles = linear_interpolation(Interpolations.deduplicate_knots!(subcycle_dates), subcycle)
            altim[!, :subcycle] = Int16.(round.(subcycles(altim.decyear))) 
        end

        #TODO: need to duplicate variable for BinStatistics, this should be fixed in BinStatistics
        hbins = 0.0:100.0:10000.0;
        dem[!, :height2] = dem[!, :height];
        ha = BinStatistics.binstats(dem, :height, hbins, :height2);
        ha[!, :bincenter] = BinStatistics.bincenter.(ha.height);
        ha[!, :area_km2] = ha.nrow .* (ref_dem_grid_res/1E3).^2;

        # find appropiate hypsometry elevation limits 
        xmax = findfirst((cumsum(ha.area_km2) ./ sum(ha.area_km2)) .> .995)+1;
        xmin = findfirst((reverse(cumsum(reverse(ha.area_km2))) ./ sum(ha.area_km2)) .< 0.995)-1;
        hbins = hbins[xmin:xmax+1];
        ha = ha[xmin:xmax,:];

        # calculate dhdt as a function of elevation
        valid = .!isnan.(dem.dhdt);
        n_ref = sum(valid)
        bs = BinStatistics.binstats(dem[valid, :], :height, hbins, :dhdt; col_function=[central, Altim.mad], missing_bins=true);
        bs[!, :bincenter] = BinStatistics.bincenter.(bs.height);
        sort!(bs, :bincenter);

        # find and fill missing bins
        valid = .!(ismissing.(bs.dhdt_central) .| isnan.(bs.dhdt_central));
        interp_vars = [:dhdt_central, :dhdt_mad];
        
        for v in interp_vars
            interp_linear = linear_interpolation(Interpolations.deduplicate_knots!(bs.bincenter[valid]), bs[:, v][valid], extrapolation_bc=Flat())
            a = bs[:, v]
            a[.!valid] = interp_linear(bs.bincenter[.!valid])
            bs[:, v] = a
        end

        # identify locations where altimetry has a valid elevation
        altim = altim[altim.altim_valid, :];

        # save single subcycle of data as geojson for plotting in QGIS
        if false
            cycle = 2

            gt = altim[altim.cycle.==cycle, :]

            gt[!, :geometry] = ArchGDAL.createpoint.(gt.longitude, gt.latitude)
            geojson_file = joinpath(data_folder, "$(region)_$(mission)_cycle$(cycle).geojson")
            if mission == :icesat2
                GDF.write(geojson_file, gt[1:10:end,:], crs=GFT.EPSG(4326))
            else
                GDF.write(geojson_file, gt, crs=GFT.EPSG(4326))
            end
        end

        # calcualte volume change for each subcycle
        usubcycle = sort(unique(altim.subcycle));
        subcycle_stats = DataFrame(
            subcycle=usubcycle,
            count=zeros(Int64, size(usubcycle)),
            exclude=falses(size(usubcycle)),
            dvdt_ref_km3yr=zeros(size(usubcycle)),
            dv_km3=zeros(size(usubcycle)),
            decyear_obs=zeros(size(usubcycle)),
            area_km2=fill(sum(ha.area_km2),size(usubcycle))
        );

        altim[!, :dh] = altim.height .- altim.height_reference;
        ylim = (-7, 1);
        xlim = extrema(hbins);
        println(xlim)

        p = plot(size=(800, 400); legend=false, xlabel="elevation [m]", ylabel="height change [m yr⁻¹]", margin=5Plots.mm, ylim=ylim, xlim=xlim);
        
        valid0 = .!isnan.(altim.dhdt)
        alpha = 0.1
        for r in eachrow(subcycle_stats)
            ind = (altim.subcycle .== r.subcycle) .& valid0
            gt = altim[ind, :]

            df = BinStatistics.binstats(gt, :height, hbins, :dhdt; col_function=[central, Altim.mad], missing_bins=true)
            df[!, :dh_central] = BinStatistics.binstats(gt, :height, hbins, :dh; col_function=[central, Altim.mad], missing_bins=true)[:, :dh_central]
            df[!, :decyear_central] = BinStatistics.binstats(gt, :height, hbins, :decyear; col_function=[mean, Altim.mad], missing_bins=true)[:, :decyear_mean]

            df[!, :bincenter] = BinStatistics.bincenter.(df.height)
            sort!(df, :bincenter)
            # linearly interpolate missing values
            valid = .!(ismissing.(df.dhdt_central) .| isnan.(df.dhdt_central));
            interp_vars = [:dh_central, :dhdt_mad, :dhdt_central, :decyear_central]

            if sum(valid) > 1
                for v in interp_vars
                    interp_linear = linear_interpolation(df.bincenter[valid], df[:, v][valid], extrapolation_bc=Flat())
                    a = df[:, v]
                    a[.!valid] = interp_linear(df.bincenter[.!valid])
                    df[:, v] = a
                end
            elseif any(valid)
                for v in interp_vars
                    df[:, v] .= df[:, v][valid]
                end
            end

            r.dvdt_ref_km3yr = sum(ha.area_km2 .* (df.dhdt_central / 1E3)) 
            r.count = nrow(gt)

            if r.subcycle == 2
                println("       subcycle 2: N = $(r.count)")
            end

            r.dv_km3 = sum(ha.area_km2 .* df.dh_central / 1E3)
            r.decyear_obs = sum(ha.area_km2 .* df.decyear_central) ./ sum(ha.area_km2)
            plot!(p, 
                df.bincenter, 
                df.dhdt_central; 
                ribbon=df.dhdt_mad, 
                fillalpha=alpha, 
                label="subcycle $(r.subcycle) [$(Int16(round(r.dvdt_ref_km3yr))) km³ yr⁻¹]"
                )
        end

        dvdt_ref = sum(ha.area_km2 .* (bs.dhdt_central / 1E3));
        subcycle_stats[!, :dvdt_error_km3yr] = subcycle_stats.dvdt_ref_km3yr .- dvdt_ref;

        fn = joinpath(data_folder, "$(region)_$(mission)_subcycle_hyps.arrow");
        Arrow.write(fn, subcycle_stats::DataFrame);

        # this needs to be -66/0.85 for Alaska  
        plot!(p, bs.bincenter, bs.dhdt_central; ribbon=bs.dhdt_mad, fillalpha=alpha, label="truth [$(Int16(round(dvdt_ref))) km³ yr⁻¹]", linewidth=2, thickness_scaling=1, color=:black)

        #annotate!(100, 0, text("truth: $(Int16(round(dmdt))) Gt yr⁻¹\n per-cycle sampling: $(Int16(round(mean(dmdt_cycle)))) ± $(round(std(dmdt_cycle), digits =1))", :black, :left, 10))

        #plot!(p, foreground_color_legend=nothing, legend=:outerright)
        bar_color = RGB(68 / 255, 152 / 255, 242 / 255);
        bar!(twinx(), ha.bincenter, ha.area_km2, ylim=(0, 2 * maximum(ha.area_km2)), label="area sq. km.", xlim=xlim, ylabel="glacier area km²", legend=false, yguidefontcolor=bar_color, y_foreground_color_axis=bar_color, y_foreground_color_text=bar_color, y_foreground_color_border=bar_color); #, foreground_color_legend=nothing, legend=:outerright)
        mad_error = (Altim.mad(subcycle_stats.dvdt_error_km3yr))
        subcycle_stats.exclude = abs.(subcycle_stats.dvdt_error_km3yr) .> (Altim.mad(subcycle_stats.dvdt_error_km3yr) * 1.4826 * 3)

        println("       truth: $(Int16(round(dvdt_ref))) km³ yr⁻¹, N = $(n_ref)")
        println("       subcycle sampling: MAD error $(round(mad_error, digits =1))")

        plot_file = joinpath(data_folder, "$(region)_$(mission)_subcycle_hyps.png")
        savefig(p, plot_file)

        p = plot(size=(800, 400); legend=false, xlabel="elevation [m]", ylabel="height change [m yr⁻¹]", margin=5Plots.mm, ylim=ylim, xlim=xlim)
        plot!(p, bs.bincenter, bs.dhdt_central; ribbon=bs.dhdt_mad, fillalpha=alpha, label="truth [$(Int16(round(dvdt_ref))) km³ yr⁻¹]", linewidth=2, thickness_scaling=1, color=:black);
        bar!(twinx(), ha.bincenter, ha.area_km2, ylim=(0, 2 * maximum(ha.area_km2)), label="area sq. km.", xlim=xlim, ylabel="glacier area km²", legend=false, yguidefontcolor=bar_color, y_foreground_color_axis=bar_color, y_foreground_color_text=bar_color, y_foreground_color_border=bar_color);
        plot_file = joinpath(data_folder, "$(region)_reference_hyps.png")
        savefig(p, plot_file)
    end

    begin
        # make timeseries
        df = DataFrame[];
        for mission in mission0
            fn = joinpath(data_folder, "$(region)_$(mission)_subcycle_hyps.arrow");
            foo = DataFrame(Arrow.Table(fn));
            foo[.!foo.exclude, :]
            foo[!,:mission] .= mission;
            push!(df, foo);
        end

        df = reduce(vcat, df)
        
        # create design matrix for full model
        #regressor = HuberRegression(fit_intercept=false, scale_penalty_with_samples=false)
        regressor = LinearRegression(fit_intercept=false);
        model = model(t; t_intercept=2010) = hcat(ones(size(t)), (t .- t_intercept), cos.(2 * pi * t), sin.(2 * pi * t));

        t0 = df.decyear_obs;
        dv0 = df.dv_km3;
  
        X0 = model(t0);
        θ = fit(regressor, X0, dv0);

        println("   Fit to all:") 
        println("       trend = $(Int64(round(θ[2]))) km³ yr⁻¹")
        println("       amp = $(Int64(round(sqrt(θ[3].^2 + θ[4].^2)))) km³")
        
        min_t, max_t = extrema(t0);
        t0 = floor(min_t):0.01:ceil(max_t);
        dv0 = model(t0) * θ;

        alpha = 0.5;
        p = plot(size=(800, 400), legend=false, ylabel="volume change [km³ yr⁻¹]", margin=5Plots.mm);

        p = plot!(p,t0, dv0 .- dv0[1], color=:gray);

        for (i, mission) in enumerate(mission0[1:end])
            ind = df.mission .== mission
            dv0 = model(df[ind, :decyear_obs]) * θ
            valid = abs.(df[ind, :dv_km3] .- dv0) .< sqrt(θ[3].^2 + θ[4].^2)
            ind[ind] = valid
            p = plot!(p, df[ind,:decyear_obs], df[ind,:dv_km3] .- dv0[1], ribbon=std(df[ind,:dvdt_error_km3yr]*10), fillalpha=alpha);
            ucycle0 = p.series_list[i].plotattributes[:linecolor];
            p = plot!(p, df[ind,:decyear_obs], df[ind,:dv_km3] .- dv0[1], color=ucycle0, seriestype=:scatter, markersize=3);
        end

        plot_file = joinpath(data_folder, "$(region)_timeseries.png")
        savefig(p, plot_file)
    end
end