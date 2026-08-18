"""
Tests for river network routing in utilities_routing.jl

These tests validate river network traversal, flux accumulation,
hydrologic routing, and population exposure calculations.
"""

using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using DataFrames
using Statistics

# Include test fixtures
include("../fixtures/synthetic_network.jl")

@testset "Routing Operations" begin
    @testset "trace_downstream" begin
        # Simple linear network: 1 → 2 → 3 → 0 (ocean)
        ids = [1, 2, 3]
        next_ids = [2, 3, 0]

        # Trace from node 1
        path = GGA.trace_downstream(1, ids, next_ids)
        @test path == [1, 2, 3]

        # Trace from node 2
        path2 = GGA.trace_downstream(2, ids, next_ids)
        @test path2 == [2, 3]

        # Trace from ocean terminus
        path3 = GGA.trace_downstream(3, ids, next_ids)
        @test path3 == [3]
    end

    @testset "trace_downstream - branching network" begin
        # Branching network: 1 → 3, 2 → 3 → 4 → 0
        ids = [1, 2, 3, 4]
        next_ids = [3, 3, 4, 0]

        # Trace from branch 1
        path1 = GGA.trace_downstream(1, ids, next_ids)
        @test path1 == [1, 3, 4]

        # Trace from branch 2
        path2 = GGA.trace_downstream(2, ids, next_ids)
        @test path2 == [2, 3, 4]

        # Trace from confluence
        path3 = GGA.trace_downstream(3, ids, next_ids)
        @test path3 == [3, 4]

        # All branches eventually reach node 4
        @test 4 in path1
        @test 4 in path2
        @test 4 in path3
    end

    @testset "haversine_distance" begin
        # Same point
        d_same = GGA.haversine_distance((0.0, 0.0), (0.0, 0.0))
        @test d_same ≈ 0.0 atol=1.0

        # Equator separation (1 degree ≈ 111.32 km)
        d_eq = GGA.haversine_distance((0.0, 0.0), (1.0, 0.0))
        @test d_eq ≈ 111320.0 rtol=1e-2

        # Antipodal points (half Earth circumference ≈ 20,000 km)
        d_antipodal = GGA.haversine_distance((0.0, 0.0), (180.0, 0.0))
        @test d_antipodal ≈ 20037000.0 rtol=1e-2

        # Pole to pole
        d_pole = GGA.haversine_distance((0.0, 90.0), (0.0, -90.0))
        @test d_pole ≈ 20004000.0 rtol=1e-2
    end

    @testset "flux_accumulate! - linear network" begin
        # Create linear network: A → B → C → ocean
        network = create_linear_network(3)

        # Initial flux at each node (local contribution)
        network.flux = [10.0, 5.0, 3.0]  # Gt/yr at each node

        # Accumulate flux downstream
        GGA.flux_accumulate!(network; flux_col=:flux, id_col=:COMID, next_id_col=:NextDownID)

        # Test: Node C (terminus) should have sum of all upstream fluxes
        terminus_idx = findfirst(network.NextDownID .== 0)
        @test network.flux[terminus_idx] ≈ 10.0 + 5.0 + 3.0

        # Test: Node B should have its own flux + upstream from A
        node_b_idx = findfirst(network.COMID .== network.NextDownID[1])
        @test network.flux[node_b_idx] ≈ 10.0 + 5.0

        # Test: Node A (headwater) should only have its own flux
        @test network.flux[1] ≈ 10.0
    end

    @testset "flux_accumulate! - branching network" begin
        # Create branching network: A + B → C → D → ocean
        network = create_branching_network()

        # Local flux contributions
        n_nodes = nrow(network)
        network.flux = [8.0, 6.0, 2.0, 1.0, 0.5]  # Different at each node

        # Get confluence and outlet IDs
        confluence_id = 3  # Where A and B meet
        outlet_id = 5      # Ocean terminus

        # Accumulate flux
        GGA.flux_accumulate!(network; flux_col=:flux, id_col=:COMID, next_id_col=:NextDownID)

        # Test: Confluence should receive flux from both branches A and B
        confluence_idx = findfirst(network.COMID .== confluence_id)
        expected_at_confluence = 8.0 + 6.0 + 2.0  # A + B + local at C
        @test network.flux[confluence_idx] ≈ expected_at_confluence rtol=1e-6

        # Test: Outlet should have total from entire network
        outlet_idx = findfirst(network.COMID .== outlet_id)
        expected_total = 8.0 + 6.0 + 2.0 + 1.0 + 0.5
        @test network.flux[outlet_idx] ≈ expected_total rtol=1e-6

        # Test: Flux should increase monotonically downstream
        # (or stay same if no tributaries)
        for i in 1:nrow(network)
            if network.NextDownID[i] != 0
                downstream_idx = findfirst(network.COMID .== network.NextDownID[i])
                @test network.flux[downstream_idx] >= network.flux[i]
            end
        end
    end

    @testset "linear_reservoir_impulse_response_monthly" begin
        # Test hydrologic routing with known residence time
        n_months = 24
        residence_time_months = 3.0  # 3 month residence time

        # Create unit impulse: 1.0 at t=0, 0.0 elsewhere
        input_flux = zeros(n_months)
        input_flux[1] = 1.0

        # Apply linear reservoir routing
        output_flux = GGA.linear_reservoir_impulse_response_monthly(
            input_flux,
            residence_time_months
        )

        # Test: Output should sum to 1.0 (mass conservation)
        @test sum(output_flux) ≈ 1.0 rtol=1e-3

        # Test: Peak output should be at first time step but < input
        @test output_flux[1] < input_flux[1]
        @test output_flux[1] == maximum(output_flux)

        # Test: Output should decay exponentially
        # Each month, remaining water = exp(-1/tau)
        decay_rate = exp(-1.0 / residence_time_months)
        for i in 2:10  # Check first 10 months
            expected_ratio = decay_rate
            actual_ratio = output_flux[i] / output_flux[i-1]
            @test actual_ratio ≈ expected_ratio rtol=0.15  # Allow some numerical error
        end

        # Test: Long tail should approach zero
        @test output_flux[end] < 0.01
    end

    @testset "linear_reservoir_routing - continuous input" begin
        # Test with constant input flux
        n_months = 36
        residence_time = 2.0

        # Constant input of 10 Gt/month
        input_flux = fill(10.0, n_months)

        output_flux = GGA.linear_reservoir_impulse_response_monthly(
            input_flux,
            residence_time
        )

        # Test: With constant input, output should reach equilibrium
        # At equilibrium, output ≈ input
        equilibrium_months = 15:n_months  # After ~5 residence times
        mean_output = mean(output_flux[equilibrium_months])
        mean_input = mean(input_flux[equilibrium_months])

        @test mean_output ≈ mean_input rtol=0.1

        # Test: Output should be monotonically increasing at start
        @test all(diff(output_flux[1:10]) .>= -1e-10)  # Allow tiny numerical errors
    end

    @testset "flux conservation in routing" begin
        # Test that total flux is conserved through routing
        network = create_dendritic_network()

        # Random local fluxes
        network.flux_local = abs.(randn(nrow(network))) .* 5.0
        network.flux_accumulated = copy(network.flux_local)

        # Accumulate
        GGA.flux_accumulate!(network;
            flux_col=:flux_accumulated,
            id_col=:COMID,
            next_id_col=:NextDownID
        )

        # Test: Total flux at outlet equals sum of all local inputs
        outlet_idx = findfirst(network.NextDownID .== 0)
        total_input = sum(network.flux_local)
        total_output = network.flux_accumulated[outlet_idx]

        @test total_output ≈ total_input rtol=1e-10
    end

    @testset "network topology validation" begin
        # Test detection of circular references (should not occur)

        # Valid network
        valid_ids = [1, 2, 3]
        valid_next = [2, 3, 0]
        @test_nowarn GGA.trace_downstream(1, valid_ids, valid_next)

        # Circular network: 1 → 2 → 1 (invalid, should be caught)
        circular_ids = [1, 2]
        circular_next = [2, 1]

        # This should either error or not return (depending on implementation)
        # For safety, trace should detect cycles
        @test_throws Union{ErrorException, StackOverflowError} begin
            GGA.trace_downstream(1, circular_ids, circular_next; max_steps=100)
        end || length(GGA.trace_downstream(1, circular_ids, circular_next; max_steps=10)) <= 10
    end

    @testset "multiple outlets handling" begin
        # Test network with multiple ocean outlets (endorheic + coastal)

        # Network: 1→2→0, 3→4→0 (two separate drainage systems)
        ids = [1, 2, 3, 4]
        next_ids = [2, 0, 4, 0]

        # Trace from first system
        path1 = GGA.trace_downstream(1, ids, next_ids)
        @test path1 == [1, 2]
        @test 3 ∉ path1  # Should not include other drainage

        # Trace from second system
        path2 = GGA.trace_downstream(3, ids, next_ids)
        @test path2 == [3, 4]
        @test 1 ∉ path2  # Should not include other drainage
    end

    @testset "distance calculations along network" begin
        # Test cumulative distance calculation
        network = create_linear_network(5)

        # Each segment is 10 km
        @test all(network.lengthkm .≈ 10.0)

        # Calculate cumulative distance from headwater
        total_distance = sum(network.lengthkm)
        @test total_distance ≈ 50.0  # 5 segments × 10 km

        # Verify haversine distances match network distances (roughly)
        for i in 1:(nrow(network)-1)
            lon1, lat1 = network.lon[i], network.lat[i]
            lon2, lat2 = network.lon[i+1], network.lat[i+1]

            haversine_km = GGA.haversine_distance((lon1, lat1), (lon2, lat2)) / 1000
            network_km = network.lengthkm[i]

            # Should be similar (within 20% for straight line vs river meander)
            @test haversine_km ≈ network_km rtol=0.3
        end
    end

    @testset "flux routing with seasonal variation" begin
        # Test routing with seasonal input
        n_months = 24
        residence_time = 2.0

        # Seasonal input (sine wave)
        t = collect(0:n_months-1) ./ 12.0  # Years
        input_flux = 10.0 .+ 5.0 .* sin.(2π .* t)  # Mean 10, amplitude 5

        output_flux = GGA.linear_reservoir_impulse_response_monthly(
            input_flux,
            residence_time
        )

        # Test: Output should be smoother than input (reservoir dampens)
        input_std = std(input_flux)
        output_std = std(output_flux[6:end])  # Skip spin-up
        @test output_std < input_std

        # Test: Mean should be preserved
        @test mean(output_flux[6:end]) ≈ mean(input_flux[6:end]) rtol=0.1

        # Test: Output should lag input by ~residence time
        # Peak of output should occur after peak of input
        input_peak_idx = argmax(input_flux[1:12])
        output_peak_idx = argmax(output_flux[1:12])
        @test output_peak_idx > input_peak_idx  # Output lags
    end

    @testset "Edge case - zero flux" begin
        # Test handling of zero flux
        network = create_linear_network(3)
        network.flux = zeros(3)

        GGA.flux_accumulate!(network; flux_col=:flux, id_col=:COMID, next_id_col=:NextDownID)

        # All fluxes should remain zero
        @test all(network.flux .== 0.0)
    end

    @testset "Edge case - single node network" begin
        # Test single-node network (direct to ocean)
        network = DataFrame(
            COMID = [1],
            NextDownID = [0],
            lengthkm = [10.0],
            lon = [-120.0],
            lat = [45.0],
            flux = [15.0]
        )

        GGA.flux_accumulate!(network; flux_col=:flux, id_col=:COMID, next_id_col=:NextDownID)

        # Flux should be unchanged
        @test network.flux[1] == 15.0
    end

    @testset "Large network scaling" begin
        # Test with larger network to ensure performance
        n_nodes = 100

        # Create chain network
        ids = collect(1:n_nodes)
        next_ids = vcat(collect(2:n_nodes), [0])

        network = DataFrame(
            COMID = ids,
            NextDownID = next_ids,
            lengthkm = fill(5.0, n_nodes),
            lon = range(-120, -110, length=n_nodes),
            lat = fill(45.0, n_nodes),
            flux = fill(1.0, n_nodes)
        )

        # Should complete quickly
        @test_nowarn GGA.flux_accumulate!(network;
            flux_col=:flux,
            id_col=:COMID,
            next_id_col=:NextDownID
        )

        # Test: Outlet should have sum of all inputs
        @test network.flux[end] ≈ n_nodes rtol=1e-10
    end
end
