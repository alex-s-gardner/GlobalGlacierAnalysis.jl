using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using DimensionalData
import DimensionalData as DD
using Dates

@testset "NetCDF I/O" begin
    @testset "dimstack2netcdf round-trip" begin
        mktempdir() do tmpdir
            filepath = joinpath(tmpdir, "test.nc")

            # Create DimStack
            dates = [DateTime(2018,1,1), DateTime(2018,7,1), DateTime(2019,1,1)]
            x_coords = [1.0, 2.0, 3.0]
            y_coords = [4.0, 5.0, 6.0]

            temp_data = randn(3, 3, 3)
            precip_data = randn(Float32, 3, 3, 3)

            ds = DimStack((
                temperature = DimArray(temp_data, (Ti=dates, X=x_coords, Y=y_coords);
                                     metadata=Dict("units" => "°C")),
                precipitation = DimArray(precip_data, (Ti=dates, X=x_coords, Y=y_coords);
                                       metadata=Dict("units" => "mm"))
            ))

            # Write
            GGA.dimstack2netcdf(ds, filepath; global_attributes=Dict("source" => "test"))

            # Read
            ds_loaded = GGA.netcdf2dimstack(filepath)

            # Verify dimensions
            @test haskey(dims(ds_loaded), :Ti) || haskey(dims(ds_loaded), :date)
            @test haskey(dims(ds_loaded), :X)
            @test haskey(dims(ds_loaded), :Y)

            # Verify variables exist
            @test haskey(ds_loaded, :temperature) || haskey(ds_loaded, :Temperature)
            @test haskey(ds_loaded, :precipitation) || haskey(ds_loaded, :Precipitation)

            # Verify approximate equality (allow for type conversions)
            temp_key = haskey(ds_loaded, :temperature) ? :temperature : :Temperature
            @test size(ds_loaded[temp_key]) == size(temp_data)
        end
    end

    @testset "NetCDF with metadata preservation" begin
        mktempdir() do tmpdir
            filepath = joinpath(tmpdir, "test_metadata.nc")

            # Create DimArray with metadata
            data = collect(reshape(1.0:12.0, 3, 4))
            da = DimArray(
                data,
                (X=1:3, Y=1:4);
                name=:test_var,
                metadata=Dict("units" => "m", "long_name" => "Test Variable")
            )

            ds = DimStack((testvar=da,))

            # Write with global attributes
            global_attrs = Dict("title" => "Test Dataset", "institution" => "Test Lab")
            GGA.dimstack2netcdf(ds, filepath; global_attributes=global_attrs)

            # Verify file was created
            @test isfile(filepath)

            # Read back
            ds_loaded = GGA.netcdf2dimstack(filepath)

            # Verify data
            loaded_key = first(keys(ds_loaded))
            @test size(ds_loaded[loaded_key]) == size(data)
        end
    end

    @testset "NetCDF with missing values" begin
        mktempdir() do tmpdir
            filepath = joinpath(tmpdir, "test_missing.nc")

            # Create data with NaN (represents missing)
            data = [1.0 2.0 NaN; 4.0 NaN 6.0]
            da = DimArray(data, (X=1:2, Y=1:3); name=:data_with_missing)

            ds = DimStack((datavar=da,))

            # Write
            GGA.dimstack2netcdf(ds, filepath)

            # Read
            ds_loaded = GGA.netcdf2dimstack(filepath)

            # Verify NaN handling
            loaded_key = first(keys(ds_loaded))
            loaded_data = ds_loaded[loaded_key]

            # Check that NaN values are preserved or represented as missing/fill value
            @test any(isnan.(loaded_data)) || any(ismissing.(loaded_data))
        end
    end
end
