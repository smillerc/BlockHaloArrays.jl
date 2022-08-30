using .Threads
using BlockHaloArrays
using ThreadPools
using Test
using EllipsisNotation

@testset "Simple" begin

    dims = (30, 20)
    B = rand(dims...)
    nhalo = 2

    nblocks = 6

    A1 = BlockHaloArray(dims, nhalo, nblocks; T=Float64);
    A2 = BlockHaloArray(B, nhalo, nblocks)

    @test length(A1.blocks) == nblocks
    @test length(A2.blocks) == nblocks


    @test A1.block_layout == (3, 2)

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
    
    sync_halo!(A1)

    blockid = 2
    hi_domn_start, hi_domn_end, hi_halo_start, hi_halo_end = BlockHaloArrays.hi_indices(A1.blocks[blockid], A1.nhalo)
    lo_halo_start, lo_halo_end, lo_domn_start, lo_domn_end = BlockHaloArrays.lo_indices(A1.blocks[blockid], A1.nhalo)

    ihi_domn_start, jhi_domn_start = hi_domn_start
    ihi_domn_end, jhi_domn_end = hi_domn_end
    ihi_halo_start, jhi_halo_start = hi_halo_start
    ihi_halo_end, jhi_halo_end = hi_halo_end
    ilo_halo_start, jlo_halo_start = lo_halo_start
    ilo_halo_end, jlo_halo_end = lo_halo_end
    ilo_domn_start, jlo_domn_start = lo_domn_start
    ilo_domn_end, jlo_domn_end = lo_domn_end

    current_block = A1.blocks[blockid]
    @test all(current_block[ilo_halo_start:ilo_halo_end, jlo_domn_start:jhi_domn_end] .== 1) # ilo
    @test all(current_block[ihi_halo_start:ihi_halo_end, jlo_domn_start:jhi_domn_end] .== 3) # ihi
    @test all(current_block[ilo_domn_start:ihi_domn_end, jlo_halo_start:jlo_halo_end] .== 0) # jlo
    @test all(current_block[ilo_domn_start:ihi_domn_end, jhi_halo_start:jhi_halo_end] .== 5) # jhi
    
end

@testset "Different Halo Dimensions" begin

    dims = (4, 30, 20)
    B = rand(dims...)
    nhalo = 2
    halodims = (2, 3)

    nblocks = 6

    A1 = BlockHaloArray(dims, halodims, nhalo, nblocks; T=Float64);
    A2 = BlockHaloArray(B, halodims, nhalo, nblocks)

    @test length(A1.blocks) == nblocks
    @test length(A2.blocks) == nblocks

    @test A1.block_layout == (3, 2)

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
    
    sync_halo!(A1)

    blockid = 2
    hi_domn_start, hi_domn_end, hi_halo_start, hi_halo_end = BlockHaloArrays.hi_indices(A1.blocks[blockid], A1.nhalo)
    lo_halo_start, lo_halo_end, lo_domn_start, lo_domn_end = BlockHaloArrays.lo_indices(A1.blocks[blockid], A1.nhalo)

    _, ihi_domn_start, jhi_domn_start = hi_domn_start
    _, ihi_domn_end, jhi_domn_end = hi_domn_end
    _, ihi_halo_start, jhi_halo_start = hi_halo_start
    _, ihi_halo_end, jhi_halo_end = hi_halo_end
    _, ilo_halo_start, jlo_halo_start = lo_halo_start
    _, ilo_halo_end, jlo_halo_end = lo_halo_end
    _, ilo_domn_start, jlo_domn_start = lo_domn_start
    _, ilo_domn_end, jlo_domn_end = lo_domn_end

    current_block = A1.blocks[blockid]
    @test all(current_block[:, ilo_halo_start:ilo_halo_end, jlo_domn_start:jhi_domn_end] .== 1) # ilo
    @test all(current_block[:, ihi_halo_start:ihi_halo_end, jlo_domn_start:jhi_domn_end] .== 3) # ihi
    @test all(current_block[:, ilo_domn_start:ihi_domn_end, jlo_halo_start:jlo_halo_end] .== 0) # jlo
    @test all(current_block[:, ilo_domn_start:ihi_domn_end, jhi_halo_start:jhi_halo_end] .== 5) # jhi
    
end
