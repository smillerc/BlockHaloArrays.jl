using .Threads
using BlockHaloArrays
using ThreadPools
using Test

@testset "Simple" begin

    dims = (30, 20)
    B = rand(dims...)
    nhalo = 2

    nblocks = 6

    A1 = BlockHaloArray(dims, nhalo, nblocks; T=Float64);
    A2 = BlockHaloArray(B, nhalo, nblocks)

    @test length(A1.blocks) == nblocks
    @test length(A2.blocks) == nblocks


    @test A1.blockdims == (3, 2)

    # for a (3,2) block layout, here are the block id's
    #
    #    *---* *---* *---*
    #    | 4 | | 5 | | 6 |
    # j  *---* *---* *---*
    #    *---* *---* *---*
    #    | 1 | | 2 | | 3 |
    #    *---* *---* *---*
    #            i

    @test A1.neighbor_blocks[1][:ilo] == -1
    @test A1.neighbor_blocks[1][:ihi] == 2
    @test A1.neighbor_blocks[1][:jlo] == -1
    @test A1.neighbor_blocks[1][:jhi] == 4
    @test A1.neighbor_blocks[1][:ilojhi] == -1
    @test A1.neighbor_blocks[1][:ihijhi] == 5
    @test A1.neighbor_blocks[1][:ilojlo] == -1
    @test A1.neighbor_blocks[1][:ihijlo] == -1
    
    @test A1.neighbor_blocks[2][:ilo] == 1
    @test A1.neighbor_blocks[2][:ihi] == 3
    @test A1.neighbor_blocks[2][:jlo] == -1
    @test A1.neighbor_blocks[2][:jhi] == 5
    @test A1.neighbor_blocks[2][:ilojhi] == 4
    @test A1.neighbor_blocks[2][:ihijhi] == 6
    @test A1.neighbor_blocks[2][:ilojlo] == -1
    @test A1.neighbor_blocks[2][:ihijlo] == -1

    function init(A)
        dom = domainview(A, threadid())
        fill!(dom, threadid())
    end

    @sync for tid in 1:nblocks
        ThreadPools.@tspawnat tid init(A1)
    end
    
    update_halo!(A1)

    hi_domn_start, hi_domn_end, hi_halo_start, hi_halo_end = hi_indices()
    lo_halo_start, lo_halo_end, lo_domn_start, lo_domn_end = lo_indices()
    
end
