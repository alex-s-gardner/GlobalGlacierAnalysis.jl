using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA

@testset "GlobalGlacierAnalysis.jl Tests" begin
    # Include test helpers
    include("test_helpers.jl")

    # Unit tests - Pure mathematical functions
    @testset "Unit Tests" begin
        include("unit/test_coordinates.jl")
        include("unit/test_time.jl")
        include("unit/test_statistics.jl")
        include("unit/test_array_ops.jl")
        include("unit/test_extents.jl")
        include("unit/test_terrain.jl")
        include("unit/test_models.jl")
        include("unit/test_project.jl")
        include("unit/test_build_archive.jl")
        include("unit/test_sliderule.jl")
    end

    # Algorithm tests - Complex operations with synthetic data
    @testset "Algorithm Tests" begin
        include("algorithms/test_binning.jl")
        include("algorithms/test_lowlevel_binning.jl")
        include("algorithms/test_gemb.jl")
        include("algorithms/test_routing.jl")
        include("algorithms/test_synthesis.jl")
        include("algorithms/test_postprocessing.jl")
    end

    # I/O tests - Data persistence and round-trip
    @testset "I/O Tests" begin
        include("io/test_dimarray.jl")
        include("io/test_jld2.jl")
        include("io/test_netcdf.jl")
        include("io/test_readers.jl")
    end

    # Integration tests - End-to-end workflows
    @testset "Integration Tests" begin
        include("integration/test_geotile_workflow.jl")
        include("integration/test_aggregation.jl")
        include("integration/test_properties.jl")
    end
end
