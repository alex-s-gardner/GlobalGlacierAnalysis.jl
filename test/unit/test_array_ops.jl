using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA

@testset "Array Operations" begin
    # `validrange` returns one range per dimension, as a Tuple -- so a 1D input yields a 1-tuple
    # whose only element is the range, not a (first, last) pair.
    @testset "validrange 1D" begin
        # Standard case with valid range
        v = [false, true, true, false, true, false]
        range_result = GGA.validrange(v)
        @test length(range_result) == 1
        @test range_result[1] == 2:5  # first true at 2, last at 5

        # All false -> empty range, not an error
        range_empty = GGA.validrange(falses(10))
        @test isempty(range_empty[1])

        # All true
        @test GGA.validrange(trues(10))[1] == 1:10

        # Single true at start
        @test GGA.validrange([true, false, false])[1] == 1:1

        # Single true at end
        @test GGA.validrange([false, false, true])[1] == 3:3
    end

    @testset "validrange 2D" begin
        # Rectangle of true values embedded in false matrix
        mat = falses(10, 10)
        mat[3:7, 4:8] .= true

        range_2d = GGA.validrange(mat)
        @test length(range_2d) == 2
        @test range_2d[1] == 3:7  # rows
        @test range_2d[2] == 4:8  # columns

        # the returned ranges select exactly the true block back out
        @test all(mat[range_2d...])
        @test sum(mat[range_2d...]) == sum(mat)

        # all false in 2D
        empty_2d = GGA.validrange(falses(4, 4))
        @test isempty(empty_2d[1]) && isempty(empty_2d[2])
    end

    # `validgaps` returns a same-length mask, true only at invalid points *between* the first and
    # last valid point. Leading/trailing invalid runs are boundary, not gaps -- so the mask comes
    # back all-false rather than empty.
    @testset "validgaps" begin
        # Gap in middle
        gaps = GGA.validgaps(BitVector([true, false, false, true, true]))
        @test gaps == BitVector([false, true, true, false, false])
        @test count(gaps) == 2

        # No gaps (all true)
        gaps_none = GGA.validgaps(BitVector([true, true, true, true]))
        @test length(gaps_none) == 4
        @test !any(gaps_none)

        # Leading false values are boundary, not a gap
        @test !any(GGA.validgaps(BitVector([false, false, true, true, true])))

        # Trailing false values are boundary, not a gap
        @test !any(GGA.validgaps(BitVector([true, true, true, false, false])))

        # Multiple gaps
        gaps_multi = GGA.validgaps(BitVector([true, false, true, false, false, true]))
        @test gaps_multi == BitVector([false, true, false, true, true, false])
        @test count(gaps_multi) == 3

        # No valid points at all -> nothing is a gap
        @test !any(GGA.validgaps(falses(5)))
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

        # Test negative radius (erosion). Neighborhoods are clamped to the array bounds rather
        # than zero-padded, so a fully-true array has no zero neighbour anywhere and correctly
        # does not erode -- erosion has to be checked against a block with a real edge inside.
        @test sum(GGA.dilate(trues(5, 5), -1)) == 25

        block = falses(5, 5)
        block[2:4, 2:4] .= true
        eroded = GGA.dilate(block, -1)
        @test sum(eroded) < sum(block)
        @test sum(eroded) == 1        # only the centre survives
        @test eroded[3, 3]

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
