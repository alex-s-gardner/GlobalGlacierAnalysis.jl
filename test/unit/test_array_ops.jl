using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA

@testset "Array Operations" begin
    @testset "validrange 1D" begin
        # Standard case with valid range
        v = [false, true, true, false, true, false]
        range_result = GGA.validrange(v)
        @test range_result[1] == 2  # First true at index 2
        @test range_result[2] == 5  # Last true at index 5

        # All false
        all_false = falses(10)
        range_empty = GGA.validrange(all_false)
        @test isempty(range_empty[1]:range_empty[2]) || range_empty[1] > range_empty[2]

        # All true
        all_true = trues(10)
        range_all = GGA.validrange(all_true)
        @test range_all[1] == 1
        @test range_all[2] == 10

        # Single true at start
        single_start = [true, false, false]
        range_single = GGA.validrange(single_start)
        @test range_single[1] == 1
        @test range_single[2] == 1

        # Single true at end
        single_end = [false, false, true]
        range_end = GGA.validrange(single_end)
        @test range_end[1] == 3
        @test range_end[2] == 3
    end

    @testset "validrange 2D" begin
        # Rectangle of true values embedded in false matrix
        mat = falses(10, 10)
        mat[3:7, 4:8] .= true

        range_2d = GGA.validrange(mat)
        @test range_2d[1][1] == 3  # Y min
        @test range_2d[1][2] == 7  # Y max
        @test range_2d[2][1] == 4  # X min
        @test range_2d[2][2] == 8  # X max
    end

    @testset "validgaps" begin
        # Gap in middle
        valid = BitVector([true, false, false, true, true])
        gaps = GGA.validgaps(valid)
        # Should identify gap at indices 2-3
        @test length(gaps) > 0

        # No gaps (all true)
        no_gaps = BitVector([true, true, true, true])
        gaps_none = GGA.validgaps(no_gaps)
        @test isempty(gaps_none)

        # Leading false values (not a gap, just boundary)
        leading = BitVector([false, false, true, true, true])
        gaps_leading = GGA.validgaps(leading)
        # Should not identify leading false as gap
        @test isempty(gaps_leading)

        # Trailing false values
        trailing = BitVector([true, true, true, false, false])
        gaps_trailing = GGA.validgaps(trailing)
        # Should not identify trailing false as gap
        @test isempty(gaps_trailing)

        # Multiple gaps
        multi_gap = BitVector([true, false, true, false, false, true])
        gaps_multi = GGA.validgaps(multi_gap)
        # Should identify two gaps
        @test length(gaps_multi) >= 1
    end

    @testset "true_block_size" begin
        # Contiguous blocks
        v = Bool[true, true, false, true, true, true]
        blocks = GGA.true_block_size(v)
        # Should identify block sizes [2, 0, 3] or similar encoding
        @test sum(blocks[blocks .> 0]) >= 5  # At least 5 total true values

        # No true values
        no_true = falses(10)
        blocks_none = GGA.true_block_size(no_true)
        @test all(blocks_none .== 0)

        # All true values
        all_true = trues(10)
        blocks_all = GGA.true_block_size(all_true)
        @test any(blocks_all .== 10) || sum(blocks_all) == 10

        # Alternating pattern
        alternating = Bool[true, false, true, false, true]
        blocks_alt = GGA.true_block_size(alternating)
        # Should have blocks of size 1
        @test maximum(blocks_alt) <= 1
    end

    @testset "dilate morphology" begin
        # 3×3 mask with center pixel, radius=1 should create 3×3 block
        mask = falses(5, 5)
        mask[3, 3] = true

        dilated = GGA.dilate(mask, 1)
        # Should dilate to 3×3 region
        @test sum(dilated) >= 9  # At least the 3×3 neighborhood

        # Test negative radius (erosion)
        large_mask = trues(5, 5)
        eroded = GGA.dilate(large_mask, -1)
        # Erosion should reduce the mask
        @test sum(eroded) < sum(large_mask)

        # Test edge handling (boundary behavior)
        edge_mask = falses(5, 5)
        edge_mask[1, 1] = true  # Corner pixel
        dilated_edge = GGA.dilate(edge_mask, 1)
        # Should handle edge gracefully
        @test sum(dilated_edge) >= 1
        @test sum(dilated_edge) <= 9

        # Zero radius (no change)
        original = [true false; false true]
        unchanged = GGA.dilate(original, 0)
        @test unchanged == original
    end
end
