using .Threads
using BlockHaloArrays
using ThreadPools
using Test
using EllipsisNotation


function init(A)
    dom = domainview(A, threadid())
    fill!(dom, threadid())
end

@testset "Simple" begin

    dims = (30, 20)
    B = rand(dims...)
    nhalo = 2

    nblocks = 6

    A1 = BlockHaloArray(dims, nhalo, nblocks; T=Float64)
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

    # negative IDs mean it crosses over the domain boundaries, e.g. periodic
    @test A1.neighbor_blocks[1][:ilo] == -3
    @test A1.neighbor_blocks[1][:ihi] == 2
    @test A1.neighbor_blocks[1][:jlo] == -4
    @test A1.neighbor_blocks[1][:jhi] == 4
    @test A1.neighbor_blocks[1][:ilojhi] == -6
    @test A1.neighbor_blocks[1][:ihijhi] == 5
    @test A1.neighbor_blocks[1][:ilojlo] == -6
    @test A1.neighbor_blocks[1][:ihijlo] == -5

    @test A1.neighbor_blocks[2][:ilo] == 1
    @test A1.neighbor_blocks[2][:ihi] == 3
    @test A1.neighbor_blocks[2][:jlo] == -5
    @test A1.neighbor_blocks[2][:jhi] == 5
    @test A1.neighbor_blocks[2][:ilojhi] == 4
    @test A1.neighbor_blocks[2][:ihijhi] == 6
    @test A1.neighbor_blocks[2][:ilojlo] == -4
    @test A1.neighbor_blocks[2][:ihijlo] == -6

    function init(A)
        dom = domainview(A, threadid())
        fill!(dom, threadid())
    end

    @sync for tid in 1:nblocks
        ThreadPools.@tspawnat tid init(A1)
    end

    include_periodic_bc = false
    sync_halo!(A1, include_periodic_bc)

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

    A1 = BlockHaloArray(dims, halodims, nhalo, nblocks; T=Float64)
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

    # negative IDs mean it crosses over the domain boundaries, e.g. periodic
    @test A1.neighbor_blocks[1][:ilo] == -3
    @test A1.neighbor_blocks[1][:ihi] == 2
    @test A1.neighbor_blocks[1][:jlo] == -4
    @test A1.neighbor_blocks[1][:jhi] == 4
    @test A1.neighbor_blocks[1][:ilojhi] == -6
    @test A1.neighbor_blocks[1][:ihijhi] == 5
    @test A1.neighbor_blocks[1][:ilojlo] == -6
    @test A1.neighbor_blocks[1][:ihijlo] == -5

    @test A1.neighbor_blocks[2][:ilo] == 1
    @test A1.neighbor_blocks[2][:ihi] == 3
    @test A1.neighbor_blocks[2][:jlo] == -5
    @test A1.neighbor_blocks[2][:jhi] == 5
    @test A1.neighbor_blocks[2][:ilojhi] == 4
    @test A1.neighbor_blocks[2][:ihijhi] == 6
    @test A1.neighbor_blocks[2][:ilojlo] == -4
    @test A1.neighbor_blocks[2][:ihijlo] == -6

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

@testset "Different Number Types" begin
    dims = (5, 10, 20)
    nhalo = 2
    nblocks = 4

    A_f16 = BlockHaloArray(dims, nhalo, nblocks; T=Float16)
    A_f32 = BlockHaloArray(dims, nhalo, nblocks; T=Float32)
    A_f64 = BlockHaloArray(dims, nhalo, nblocks; T=Float64)
end

@testset "1D Array, 1D halo dims" begin
    dims = (20,)
    nhalo = 2
    nblocks = 4
    A = BlockHaloArray(dims, nhalo, nblocks; T=Float64)

    @test A.neighbor_blocks[1][:ilo] == -4
    @test A.neighbor_blocks[1][:ihi] == 2
    @test A.neighbor_blocks[2][:ilo] == 1
    @test A.neighbor_blocks[2][:ihi] == 3
    @test A.neighbor_blocks[3][:ilo] == 2
    @test A.neighbor_blocks[3][:ihi] == 4
    @test A.neighbor_blocks[4][:ilo] == 3
    @test A.neighbor_blocks[4][:ihi] == -1

    @sync for tid in 1:nblocks
        ThreadPools.@tspawnat tid init(A)
    end

    sync_halo!(A)

    blockid = 2
    hi_domn_start, hi_domn_end, hi_halo_start, hi_halo_end = BlockHaloArrays.hi_indices(A.blocks[blockid], A.nhalo)
    lo_halo_start, lo_halo_end, lo_domn_start, lo_domn_end = BlockHaloArrays.lo_indices(A.blocks[blockid], A.nhalo)

    ihi_halo_start = first(hi_halo_start)
    ihi_halo_end = first(hi_halo_end)
    ilo_halo_start = first(lo_halo_start)
    ilo_halo_end = first(lo_halo_end)

    current_block = A.blocks[blockid]
    @test all(current_block[ilo_halo_start:ilo_halo_end] .== 1) # ilo
    @test all(current_block[ihi_halo_start:ihi_halo_end] .== 3) # ihi
end

@testset "2D Array, 1D halo dims" begin
    dims = (10, 20)
    halodims = (2,)
    nhalo = 2
    nblocks = 4
    A = BlockHaloArray(dims, nhalo, nblocks; T=Float64)
end

@testset "2D Array, 2D halo dims" begin
    dims = (10, 20)
    nhalo = 2
    nblocks = 4
    A = BlockHaloArray(dims, nhalo, nblocks; T=Float64)
end

@testset "3D Array, 1D halo dims" begin
    dims = (10, 20, 30)
    halodims = (2,)
    nhalo = 2
    nblocks = 4
    A = BlockHaloArray(dims, nhalo, nblocks; T=Float64)
end

@testset "3D Array, 2D halo dims" begin
    dims = (10, 20, 30)
    halodims = (2, 3)
    nhalo = 2
    nblocks = 4
    A = BlockHaloArray(dims, nhalo, nblocks; T=Float64)
end

@testset "3D Array, 3D halo dims" begin
    dims = (10, 20, 30)
    nhalo = 2
    nblocks = 6
    A = BlockHaloArray(dims, nhalo, nblocks; T=Float64)

    @sync for tid in 1:nblocks
        ThreadPools.@tspawnat tid init(A)
    end

    sync_halo!(A)

end

@testset "Flat index testing" begin

    dims = (4, 50, 50)
    global_idx = (1, 26, 26)
    halodims = (2, 3)
    A = BlockHaloArray(dims, halodims, 2, 2)

    block_idx = ntuple(i ->
            BlockHaloArrays.get_block_idx(global_idx[A.halodims[i]],
                A._cummulative_blocksize_per_dim,
                A.halodims[i],
                A.halodims),
        length(A.halodims))

    local_idx = ntuple(i ->
            BlockHaloArrays.get_local_idx(global_idx,
                block_idx,
                A.nhalo,
                A._cummulative_blocksize_per_dim,
                i,
                A.halodims),
        length(A.globaldims))

    @test block_idx == (2, 1)
    @test local_idx == (1, 3, 28)
end