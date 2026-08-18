using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using DimensionalData
import DimensionalData as DD
using Dates

@testset "DimArray Operations" begin
    @testset "DimArray creation and indexing" begin
        # Create basic DimArray
        data = rand(10, 20)
        da = DimArray(data, (Ti=1:10, X=1:20))

        @test size(da) == (10, 20)
        @test dims(da, Ti) isa DD.Dim{:Ti}
        @test dims(da, X) isa DD.Dim{:X}

        # Index with At()
        val = da[Ti=DD.At(5), X=DD.At(10)]
        @test val == data[5, 10]

        # Index with range
        sub = da[Ti=1:5]
        @test size(sub) == (5, 20)

        # Check dimension labels preserved
        @test haskey(dims(sub), :Ti)
        @test haskey(dims(sub), :X)
    end

    @testset "DimArray with DateTime" begin
        dates = DateTime(2018,1,1):Month(1):DateTime(2018,12,1)
        data = randn(length(dates), 5)

        da_time = DimArray(data, (date=dates, height=1:5))

        # Index by date
        jan_data = da_time[date=DD.At(DateTime(2018,1,1))]
        @test length(jan_data) == 5

        # Verify dimension types
        @test dims(da_time, :date) isa DD.Dim{:date}
    end

    @testset "DimStack operations" begin
        # Create two DimArrays
        temp_data = randn(10, 10)
        press_data = randn(10, 10)

        temp_array = DimArray(temp_data, (X=1:10, Y=1:10); name=:temperature)
        press_array = DimArray(press_data, (X=1:10, Y=1:10); name=:pressure)

        # Create DimStack
        stack = DimStack((temperature=temp_array, pressure=press_array))

        # Access layers
        @test stack[:temperature] == temp_array
        @test stack[:pressure] == press_array

        # Check dimensions
        @test size(stack) == (10, 10)

        # Convert to Dict
        dict_stack = Dict(stack)
        @test haskey(dict_stack, :temperature)
        @test haskey(dict_stack, :pressure)
    end

    @testset "Concatenation along new dimension" begin
        da1 = DimArray(ones(5, 5), (X=1:5, Y=1:5))
        da2 = DimArray(ones(5, 5) .* 2, (X=1:5, Y=1:5))

        # Concatenate along new :run dimension
        combined = cat(da1, da2, dims=DD.Dim{:run}([:run1, :run2]))

        @test size(combined) == (5, 5, 2)
        @test combined[X=1, Y=1, run=DD.At(:run1)] == 1.0
        @test combined[X=1, Y=1, run=DD.At(:run2)] == 2.0
    end

    @testset "Reduction operations" begin
        data = reshape(1:20, 4, 5)
        da = DimArray(data, (Ti=1:4, X=1:5))

        # Mean over Ti dimension
        mean_ti = dropdims(mean(da, dims=Ti), dims=Ti)
        @test size(mean_ti) == (5,)

        # Sum over X dimension
        sum_x = dropdims(sum(da, dims=X), dims=X)
        @test size(sum_x) == (4,)
    end
end
