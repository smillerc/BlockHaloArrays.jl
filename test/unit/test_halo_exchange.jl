@testitem "2D Halo, 2D Array" begin
    include("common.jl")

    dims = (50, 50)
    nhalo = 2
    nblocks = 6

    A = BlockHaloArray(dims, nhalo, nblocks; T=Int64)
    @test length(A.blocks) == nblocks
    @test A.block_layout == (3, 2)

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

    blk_neighbor = [
        Dict(:ilo => -3, :ihi => 2, :jlo => -4, :jhi => 4, :ilojhi => -6, :ihijhi => 5, :ilojlo => -6, :ihijlo => -5), # 1
        Dict(:ilo => 1, :ihi => 3, :jlo => -5, :jhi => 5, :ilojhi => 4, :ihijhi => 6, :ilojlo => -4, :ihijlo => -6), # 2
        Dict(:ilo => 2, :ihi => -1, :jlo => -6, :jhi => 6, :ilojhi => 5, :ihijhi => -4, :ilojlo => -5, :ihijlo => -4), # 3
        Dict(:ilo => -6, :ihi => 5, :jlo => 1, :jhi => -1, :ilojhi => -3, :ihijhi => -2, :ilojlo => -3, :ihijlo => 2), # 4
        Dict(:ilo => 4, :ihi => 6, :jlo => 2, :jhi => -2, :ilojhi => -1, :ihijhi => -3, :ilojlo => 1, :ihijlo => 3), # 5
        Dict(:ilo => 5, :ihi => -4, :jlo => 3, :jhi => -3, :ilojhi => -2, :ihijhi => -1, :ilojlo => 2, :ihijlo => -1), # 6
    ]

    for blockid in 1:nblocks
        for key in keys(first(blk_neighbor))
            @test A.neighbor_blocks[blockid][key] == blk_neighbor[blockid][key]
        end
    end

    @sync for tid in 1:nblocks
        @tspawnat tid init(A, tid)
    end

    include_periodic_bc = true
    updatehalo!(A, include_periodic_bc)

    for blockid in 1:nblocks
        dv = domainview(A, blockid)
        @test all(dv .== blockid)
    end

    for blockid in 1:nblocks
        hi_domn_start, hi_domn_end, hi_halo_start, hi_halo_end = BlockHaloArrays.hi_indices(A.blocks[blockid], A.nhalo)
        lo_halo_start, lo_halo_end, lo_domn_start, lo_domn_end = BlockHaloArrays.lo_indices(A.blocks[blockid], A.nhalo)

        ihi_domn_start, jhi_domn_start = hi_domn_start
        ihi_domn_end, jhi_domn_end = hi_domn_end
        ihi_halo_start, jhi_halo_start = hi_halo_start
        ihi_halo_end, jhi_halo_end = hi_halo_end
        ilo_halo_start, jlo_halo_start = lo_halo_start
        ilo_halo_end, jlo_halo_end = lo_halo_end
        ilo_domn_start, jlo_domn_start = lo_domn_start
        ilo_domn_end, jlo_domn_end = lo_domn_end

        current_block = A.blocks[blockid]
        ilo = A.neighbor_blocks[blockid][:ilo] |> Float64 |> abs
        ihi = A.neighbor_blocks[blockid][:ihi] |> Float64 |> abs
        jlo = A.neighbor_blocks[blockid][:jlo] |> Float64 |> abs
        jhi = A.neighbor_blocks[blockid][:jhi] |> Float64 |> abs

        @test all(current_block[ilo_halo_start:ilo_halo_end, jlo_domn_start:jhi_domn_end] .== ilo)
        @test all(current_block[ihi_halo_start:ihi_halo_end, jlo_domn_start:jhi_domn_end] .== ihi)
        @test all(current_block[ilo_domn_start:ihi_domn_end, jlo_halo_start:jlo_halo_end] .== jlo)
        @test all(current_block[ilo_domn_start:ihi_domn_end, jhi_halo_start:jhi_halo_end] .== jhi)
    end

end

@testitem "2D Halo, 3D Array" begin
    include("common.jl")

    dims = (4, 30, 20)
    halodims = (2, 3)
    nhalo = 2
    nblocks = 6

    A = BlockHaloArray(dims, halodims, nhalo, nblocks; T=Int64)
    @test length(A.blocks) == nblocks
    @test A.block_layout == (3, 2)

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
    blk_neighbor = [
        Dict(:ilo => -3, :ihi => 2, :jlo => -4, :jhi => 4, :ilojhi => -6, :ihijhi => 5, :ilojlo => -6, :ihijlo => -5), # 1
        Dict(:ilo => 1, :ihi => 3, :jlo => -5, :jhi => 5, :ilojhi => 4, :ihijhi => 6, :ilojlo => -4, :ihijlo => -6), # 2
        Dict(:ilo => 2, :ihi => -1, :jlo => -6, :jhi => 6, :ilojhi => 5, :ihijhi => -4, :ilojlo => -5, :ihijlo => -4), # 3
        Dict(:ilo => -6, :ihi => 5, :jlo => 1, :jhi => -1, :ilojhi => -3, :ihijhi => -2, :ilojlo => -3, :ihijlo => 2), # 4
        Dict(:ilo => 4, :ihi => 6, :jlo => 2, :jhi => -2, :ilojhi => -1, :ihijhi => -3, :ilojlo => 1, :ihijlo => 3), # 5
        Dict(:ilo => 5, :ihi => -4, :jlo => 3, :jhi => -3, :ilojhi => -2, :ihijhi => -1, :ilojlo => 2, :ihijlo => -1), # 6
    ]

    for blockid in 1:nblocks
        for key in keys(first(blk_neighbor))
            @test A.neighbor_blocks[blockid][key] == blk_neighbor[blockid][key]
        end
    end

    @sync for tid in 1:nblocks
        @tspawnat tid init(A, tid)
    end

    for blockid in 1:nblocks
        dv = domainview(A, blockid)
        @test all(dv .== blockid)
    end

    include_periodic_bc = true
    updatehalo!(A, include_periodic_bc)

    for blockid in 1:nblocks
        hi_domn_start, hi_domn_end, hi_halo_start, hi_halo_end = BlockHaloArrays.hi_indices(A.blocks[blockid], A.nhalo, A.halodims)
        lo_halo_start, lo_halo_end, lo_domn_start, lo_domn_end = BlockHaloArrays.lo_indices(A.blocks[blockid], A.nhalo, A.halodims)

        _, ihi_domn_start, jhi_domn_start = hi_domn_start
        _, ihi_domn_end, jhi_domn_end = hi_domn_end
        _, ihi_halo_start, jhi_halo_start = hi_halo_start
        _, ihi_halo_end, jhi_halo_end = hi_halo_end
        _, ilo_halo_start, jlo_halo_start = lo_halo_start
        _, ilo_halo_end, jlo_halo_end = lo_halo_end
        _, ilo_domn_start, jlo_domn_start = lo_domn_start
        _, ilo_domn_end, jlo_domn_end = lo_domn_end

        current_block = A.blocks[blockid]
        ilo = A.neighbor_blocks[blockid][:ilo] |> Float64 |> abs
        ihi = A.neighbor_blocks[blockid][:ihi] |> Float64 |> abs
        jlo = A.neighbor_blocks[blockid][:jlo] |> Float64 |> abs
        jhi = A.neighbor_blocks[blockid][:jhi] |> Float64 |> abs

        @test all(current_block[:, ilo_halo_start:ilo_halo_end, jlo_domn_start:jhi_domn_end] .== ilo)
        @test all(current_block[:, ihi_halo_start:ihi_halo_end, jlo_domn_start:jhi_domn_end] .== ihi)
        @test all(current_block[:, ilo_domn_start:ihi_domn_end, jlo_halo_start:jlo_halo_end] .== jlo)
        @test all(current_block[:, ilo_domn_start:ihi_domn_end, jhi_halo_start:jhi_halo_end] .== jhi)
    end

end

@testitem "1D Array, 1D halo dims" begin
    include("common.jl")

    dims = (20,)
    nhalo = 2
    nblocks = 4
    A = BlockHaloArray(dims, nhalo, nblocks; T=Int32)

    @test A.neighbor_blocks[1][:ilo] == -4
    @test A.neighbor_blocks[1][:ihi] == 2
    @test A.neighbor_blocks[2][:ilo] == 1
    @test A.neighbor_blocks[2][:ihi] == 3
    @test A.neighbor_blocks[3][:ilo] == 2
    @test A.neighbor_blocks[3][:ihi] == 4
    @test A.neighbor_blocks[4][:ilo] == 3
    @test A.neighbor_blocks[4][:ihi] == -1

    @sync for tid in 1:nblocks
        @tspawnat tid init(A, tid)
    end

    include_periodic_bc = true
    updatehalo!(A, include_periodic_bc)

    for blockid in 1:nblocks
        hi_domn_start, hi_domn_end, hi_halo_start, hi_halo_end = BlockHaloArrays.hi_indices(A.blocks[blockid], A.nhalo)
        lo_halo_start, lo_halo_end, lo_domn_start, lo_domn_end = BlockHaloArrays.lo_indices(A.blocks[blockid], A.nhalo)

        ihi_halo_start = first(hi_halo_start)
        ihi_halo_end = first(hi_halo_end)
        ilo_halo_start = first(lo_halo_start)
        ilo_halo_end = first(lo_halo_end)

        current_block = A.blocks[blockid]
        ilo = A.neighbor_blocks[blockid][:ilo] |> Float64 |> abs
        ihi = A.neighbor_blocks[blockid][:ihi] |> Float64 |> abs

        @test all(current_block[ilo_halo_start:ilo_halo_end] .== ilo)
        @test all(current_block[ihi_halo_start:ihi_halo_end] .== ihi)
    end

end

@testitem "3D Array, 3D halo dims" begin
    include("common.jl")

    dims = (10, 20, 30)
    nhalo = 2
    nblocks = 6
    A = BlockHaloArray(dims, nhalo, nblocks; T=Float64)

    @sync for tid in 1:nblocks
        @tspawnat tid init(A, tid)
    end

    include_periodic_bc = true
    updatehalo!(A, include_periodic_bc)

    for blockid in 1:nblocks
        hi_domn_start, hi_domn_end, hi_halo_start, hi_halo_end = BlockHaloArrays.hi_indices(A.blocks[blockid], A.nhalo)
        lo_halo_start, lo_halo_end, lo_domn_start, lo_domn_end = BlockHaloArrays.lo_indices(A.blocks[blockid], A.nhalo)

        ihi_domn_start, jhi_domn_start, khi_domn_start = hi_domn_start
        ihi_domn_end, jhi_domn_end, khi_domn_end = hi_domn_end
        ihi_halo_start, jhi_halo_start, khi_halo_start = hi_halo_start
        ihi_halo_end, jhi_halo_end, khi_halo_end = hi_halo_end
        ilo_halo_start, jlo_halo_start, klo_halo_start = lo_halo_start
        ilo_halo_end, jlo_halo_end, klo_halo_end = lo_halo_end
        ilo_domn_start, jlo_domn_start, klo_domn_start = lo_domn_start
        ilo_domn_end, jlo_domn_end, klo_domn_end = lo_domn_end

        current_block = A.blocks[blockid]
        ilo = A.neighbor_blocks[blockid][:ilo] |> Float64 |> abs
        ihi = A.neighbor_blocks[blockid][:ihi] |> Float64 |> abs
        jlo = A.neighbor_blocks[blockid][:jlo] |> Float64 |> abs
        jhi = A.neighbor_blocks[blockid][:jhi] |> Float64 |> abs
        klo = A.neighbor_blocks[blockid][:klo] |> Float64 |> abs
        khi = A.neighbor_blocks[blockid][:khi] |> Float64 |> abs

        @test all(current_block[ilo_halo_start:ilo_halo_end, jlo_domn_start:jhi_domn_end, klo_domn_start:khi_domn_end] .== ilo)
        @test all(current_block[ihi_halo_start:ihi_halo_end, jlo_domn_start:jhi_domn_end, klo_domn_start:khi_domn_end] .== ihi)

        @test all(current_block[ilo_domn_start:ihi_domn_end, jlo_halo_start:jlo_halo_end, klo_domn_start:khi_domn_end] .== jlo)
        @test all(current_block[ilo_domn_start:ihi_domn_end, jhi_halo_start:jhi_halo_end, klo_domn_start:khi_domn_end] .== jhi)

        @test all(current_block[ilo_domn_start:ihi_domn_end, jlo_domn_start:jhi_domn_end, klo_halo_start:klo_halo_end] .== klo)
        @test all(current_block[ilo_domn_start:ihi_domn_end, jlo_domn_start:jhi_domn_end, khi_halo_start:khi_halo_end] .== khi)
    end

end
