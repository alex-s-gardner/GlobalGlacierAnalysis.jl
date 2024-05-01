using Altim
using Interpolations
using OffsetArrays
using ImageFiltering
using GeoArrays
using Test
using Plots
using GeoInterface
using AWS
using Statistics
using Proj

@testset "GeoArray crop" begin
    # dem with slope of 1 in x and 1 in y
    zx = collect(1.:1.:100.);
    zy = collect(1.:2.:200.);
    dem = zx*ones(1, length(zy)) .+ ones(length(zx),1)*(zy')

    ga = GeoArray(dem)
    bbox!(ga, (min_x=.5, min_y=.5, max_x=100.5, max_y=200.5))  
    epsg!(ga, 32625)

    @warn("need to use compat GeoArrays 0.8.4 once released")
    # test cropping in GeoArrays [remove later]
    extent = (min_x=-100., min_y=-100., max_x=200., max_y=200.)
    bbox_overlap(GeoArrays.bbox(ga), extent)
    @test ga == crop(ga, extent)
end

@testset "slope extraction from DEM" begin
    # points
    x = collect(10:20)
    y = collect(10:20)

    # dem with slope of 1 in x and 1 in y
    zx = collect(1.0:1.0:100.0) # for GeoArrays x is rows
    zy = collect(1.0:2.0:200.0) # for GeoArrays y is columns
    dem = zx * ones(1, length(zy)) .+ ones(length(zx), 1) * (zy')

    ga = GeoArray(dem)
    bbox!(ga, (min_x=0.5, min_y=0.5, max_x=100.5, max_y=200.5))
    epsg!(ga, 32625)

    # set slope kernal as used in dem_height 
    dx = centered([-1. -2. -1.; 0. 0. 0.; 1. 2. 1.]) ./ (1*8) # for GeoArrays x is rows
    dy = centered([-1. 0. 1.; -2. 0. 2.; -1. 0. 1.]) ./ (1*8) # for GeoArrays y is columns
    dk = [dx,dy]

    k = centered(ones(3, 3)./(3*3))

    interp = Constant()

    z, dhdx, dhdy = pointextract(x, y, "EPSG:32625", ga; 
        nodatavalue = NaN, filter_kernel = k, derivative_kernel = dk, 
        interp = interp, replace_nodatavalue_withnans = false)

    @test(all(round.(z, digits=10) .== [19.0, 22.0, 23.0, 26.0, 27.0, 30.0, 31.0, 34.0, 35.0, 38.0, 39.0]))
    @test(all(round.(dhdx,digits = 10) .== 1.))
    @test(all(round.(dhdy,digits = 10) .== 2.))
end

@testset "pointextract" begin
    # points
    x = collect(10:20)
    y = collect(10:20)

    # dem with slope of 1 in x and 1 in y
    dhdx0 = 1.0
    dhdy0 = 2.0
    zx = collect(1.0:dhdx0:(100.0*dhdx0)) # for GeoArrays x is rows
    zy = collect(1.0:2.0:(100.0*dhdy0)) # for GeoArrays y is columns
    dem = zx * ones(1, length(zy)) .+ ones(length(zx), 1) * (zy')

    ga = GeoArray(dem)
    bbox!(ga, (min_x=0.5, min_y=0.5, max_x=100.5, max_y=200.5))
    epsg!(ga, 32625)

    # extraction kernel
    k = centered(ones(3, 3) ./ (3 * 3)) # for GeoArrays x is rows, y is columns

    z = pointextract(x, y, "EPSG:32625", ga;
    nodatavalue=NaN, filter_kernel=k, derivative_kernel=nothing,
    interp=Linear(), replace_nodatavalue_withnans=false)
        
    @test(all(round.(z, digits=10) .== collect(19.5:2.:39.5)))

    # set slope kernal
    dx = centered([-1.0 -2.0 -1.0; 0.0 0.0 0.0; 1.0 2.0 1.0]) ./ (1 * 8) # for GeoArrays x is rows
    dy = centered([-1.0 0.0 1.0; -2.0 0.0 2.0; -1.0 0.0 1.0]) ./ (1 * 8) # for GeoArrays y is columns
    dk = [dx, dy]

    z, dhdx, dhdy = pointextract(x, y, "EPSG:32625", ga; 
        nodatavalue = NaN, filter_kernel = k, derivative_kernel = dk, 
        interp=Constant(), replace_nodatavalue_withnans=false)

    @test(all(round.(z, digits=10) .== [19.0, 22.0, 23.0, 26.0, 27.0, 30.0, 31.0, 34.0, 35.0, 38.0, 39.0]))
    @test(all(round.(dhdx, digits=10) .== dhdx0))
    @test(all(round.(dhdy,digits = 10) .== dhdy0))
end

@testset "track_offset" begin
    # points
    x = collect(30:170);
    y = collect(10:150);

    # Create a dem with 4 gaussian peaks
    dem = Kernel.gaussian([25, 25]);
    dem = hcat(dem.parent, dem.parent*5);
    dem = vcat(dem, dem*2);

    ga = GeoArray(dem);
    bbox!(ga, (min_x=0.5, min_y=0.5, max_x=202.5, max_y=202.5));
    epsg!(ga, 32625);

    # extraction kernel
    k = centered(ones(3, 3) ./ (3 * 3)) # for GeoArrays x is rows, y is columns

    # Notes: if line is moved too far away then slopes will not be representative
    # this is an inherent limit of the approach
    x0 = x .+ 4.;
    y0 = y .+ -3.;


    z0 = pointextract(x0, y0, "EPSG:32625", ga; 
        nodatavalue = NaN, filter_kernel = k, derivative_kernel = nothing, 
        interp=Linear(), replace_nodatavalue_withnans=false);


    # set slope kernal
    dx = centered([-1.0 -2.0 -1.0; 0.0 0.0 0.0; 1.0 2.0 1.0]) ./ (1 * 8) # for GeoArrays x is rows
    dy = centered([-1.0 0.0 1.0; -2.0 0.0 2.0; -1.0 0.0 1.0]) ./ (1 * 8) # for GeoArrays y is columns
    dk = [dx, dy]

    z, dhdx, dhdy = pointextract(x, y, "EPSG:32625", ga;
        nodatavalue = NaN, filter_kernel = k, derivative_kernel = dk, 
        interp=Linear(), replace_nodatavalue_withnans=false);

    dx, dy, dz = Altim.track_offset(dhdx, dhdy,
        z0 .- z;
        interations=3,
        iter_thresh=7)

    z1, dhdx, dhdy = pointextract(x0 .- round(dx), y0 .- (round(dy)), "EPSG:32625", ga;
        nodatavalue=NaN, filter_kernel=k, derivative_kernel=dk,
        interp=Linear(), replace_nodatavalue_withnans=false);

    @test(z1 == z) 
end


@testset "track_offset for real DEM" begin

    #TODO: setup automatic download if file does not exist
    #aws_config = AWSConfig(;creds=nothing)
    #aws_config.region = "eu-central-1"
    #a = S3.list_objects("copernicus-dem-30m/"; aws_config)

    #;aws s3 cp --no-sign-request s3://copernicus-dem-30m/Copernicus_DSM_COG_10_N51_00_W125_00_DEM/Copernicus_DSM_COG_10_N51_00_W125_00_DEM.tif /net/devon/mnt/devon2-r1/devon0/gardnera/Documents/GitHub/Altim/data/

    dem = GeoArrays.read("data/Copernicus_DSM_COG_10_N51_00_W125_00_DEM.tif")
    # plot(dem)

    # axis dimentions
    lonDim, latDim = GeoArrays.ranges(dem)

    # extract along a line
    lon = collect(lonDim)[100:2100]
    lat = collect(latDim)[700:2700]

    # find local utm epsg and tranform from WGS84
    epsg = utm_epsg(mean(lon),mean(lat); always_xy=true)
    trans = Proj.Transformation("EPSG:4326", epsg, always_xy=true)

    xy = trans.(lon, lat)

    # plot line on map
    #scatter!(lon,lat)

    lonlat = [(lon[i],lat[i]) for i in eachindex(lon)]
    ind = CartesianIndex.(Tuple.(indices.(Ref(dem), lonlat)))

    z_og = dem[ind]
    #plot(z_og)

    k = centered(ones(1, 1) ./ (1))

    # set slope kernal
    dx = centered([-1.0 -2.0 -1.0; 0.0 0.0 0.0; 1.0 2.0 1.0]) ./ (1 * 8) # for GeoArrays x is rows
    dy = centered([-1.0 0.0 1.0; -2.0 0.0 2.0; -1.0 0.0 1.0]) ./ (1 * 8) # for GeoArrays y is columns
    dk = [dx, dy]

    z0, dhdlon0, dhdlat0 = pointextract(lon, lat, "EPSG:4326", dem; 
    nodatavalue = NaN, filter_kernel = k, derivative_kernel = dk, 
    interp = Linear(), replace_nodatavalue_withnans = false)
    # plot!(z0)
    @test .!any(abs.(z_og .- z0) .> 0.00001)

    dhdlon = copy(dem)
    dhdlon.A[:,:,1] = imfilter(dem.A[:,:,1],dx);
    #plot(dhdlon0)
    #plot!(dhdlon[ind])=
    @test .!any(abs.(dhdlon[ind] .- dhdlon0) .> 0.00001)

    #plot(dhdlon)

    dhdlat = copy(dem)
    dhdlat.A[:, :, 1] = imfilter(dem.A[:, :, 1], dy);
    #plot(dhdlat0)
    #plot!(dhdlat[ind])

    @test .!any(abs.(dhdlat[ind] .- dhdlat0) .> 0.00001)

    #plot(dhdlat)

    x_offset = 10
    y_offset = -10

    x = getindex.(xy, 1) .+ x_offset
    y = getindex.(xy, 2) .+ y_offset

    lonlat_offset = inv(trans).(x, y)

    z0, dhdlon, dhdlat = pointextract(getindex.(lonlat_offset, 1), getindex.(lonlat_offset, 2), "EPSG:4326", dem;
        nodatavalue=NaN, filter_kernel=k, derivative_kernel=dk,
        interp=Linear(), replace_nodatavalue_withnans=false)
    dhdlon = dhdlon ./ dem.f.linear[1, 1]
    dhdlat = dhdlat ./ dem.f.linear[2, 2]
    # plot(z)
    # plot!(z0)

    dhdx = similar(dhdlon)
    dhdy = similar(dhdlat)
    for i = eachindex(lon)
        cxy = trans(lon[i], lat[i])
        pt_dx = inv(trans).([cxy[1] - 0.5, cxy[1] + 0.5], [cxy[2], cxy[2]])
        pt_dy = inv(trans).([cxy[1], cxy[1]], [cxy[2] - 0.5, cxy[2] + 0.5])

        dhdx[i], dhdy[i] = (
            (dhdlon[i] * (pt_dx[2][1] - pt_dx[1][1])) + (dhdlat[i] * (pt_dx[2][2] - pt_dx[1][2])),
            (dhdlon[i] * (pt_dy[2][1] - pt_dy[1][1])) + (dhdlat[i] * (pt_dy[2][2] - pt_dy[1][2])))
    end
            
    dh = z_og .- z0
    dx, dy, dz = Altim.track_offset(dhdx, dhdy,
        dh;
        interations=3,
        iter_thresh=7);

    # check if good within 10% or 10 cm
    @test .!any(abs.(x_offset .+ dx) .> max((abs(x_offset) * 0.1),.1))
    @test .!any(abs.(y_offset .+ dy) .> max((abs(y_offset) * 0.1), .1))

    # now test offset_geo_fast
    dx, dy, dz, epsg = Altim.offset_geo_fast(lon, lat, dh, dhdlon, dhdlat)

    @test .!any(abs.(x_offset .+ dx) .> max((abs(x_offset) * 0.15),.1))
    @test .!any(abs.(y_offset .+ dy) .> max((abs(y_offset) * 0.1), .1))
    @test .!any(abs.(dz) .> .1)
    @test "EPSG:32610" == epsg

    # now test offset_geo
     dx, dy, dz, epsg = Altim.offset_geo(lon, lat, dh, dhdlon, dhdlat)

    @test .!any(abs.(x_offset .+ dx) .> max((abs(x_offset) * 0.15),.1))
    @test .!any(abs.(y_offset .+ dy) .> max((abs(y_offset) * 0.1), .1))
    @test .!any(abs.(dz) .> .1)
    @test "EPSG:32610" == epsg
end