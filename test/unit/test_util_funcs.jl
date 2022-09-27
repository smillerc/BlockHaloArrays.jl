@testitem "halo/donorview" begin
    include("common.jl")
    dims = (50, 50)
    nhalo = 2
    nblocks = 4
    A = BlockHaloArray(dims, nhalo, nblocks)

    @test size(donorview(A, 1, :ilo)) == (2, 25)
    @test size(haloview(A, 1, :ilo)) == (2, 25)
    @test size(donorview(A, 1, :jhi)) == (25, 2)
    @test size(haloview(A, 1, :jhi)) == (25, 2)
end

@testitem "nblocks" begin
    include("common.jl")
    A = BlockHaloArray((50,50), 2, 4)
    @test nblocks(A) == 4
end

@testitem "eltype" begin
    include("common.jl")
    A = BlockHaloArray((50,50), 2, 4, T=Int32)
    @test eltype(A) == Int32
end

@testitem "size" begin
    include("common.jl")
    A = BlockHaloArray((50,50), 2, 4, T=Int32)
    @test size(A) == (50, 50)
end

@testitem "onboundary" begin
    include("common.jl")

    dims = (50, 50)
    nhalo = 2
    nblocks = 6

    A = BlockHaloArray(dims, nhalo, nblocks; T=Int64)

    # for a (3,2) block layout, here are the block id's
    #
    #    *---* *---* *---*
    #    | 4 | | 5 | | 6 |
    # j  *---* *---* *---*
    #    *---* *---* *---*
    #    | 1 | | 2 | | 3 |
    #    *---* *---* *---*
    #            i

    @test onboundary(A, 1, :ilo) == true
    @test onboundary(A, 1, :ihi) == false
    @test onboundary(A, 1, :jlo) == true
    @test onboundary(A, 1, :jhi) == false

    @test onboundary(A, 2, :ilo) == false
    @test onboundary(A, 2, :ihi) == false
    @test onboundary(A, 2, :jlo) == true
    @test onboundary(A, 2, :jhi) == false

    @test onboundary(A, 3, :ilo) == false
    @test onboundary(A, 3, :ihi) == true
    @test onboundary(A, 3, :jlo) == true
    @test onboundary(A, 3, :jhi) == false

    @test onboundary(A, 4, :ilo) == true
    @test onboundary(A, 4, :ihi) == false
    @test onboundary(A, 4, :jlo) == false
    @test onboundary(A, 4, :jhi) == true

    @test onboundary(A, 5, :ilo) == false
    @test onboundary(A, 5, :ihi) == false
    @test onboundary(A, 5, :jlo) == false
    @test onboundary(A, 5, :jhi) == true

    @test onboundary(A, 6, :ilo) == false
    @test onboundary(A, 6, :ihi) == true
    @test onboundary(A, 6, :jlo) == false
    @test onboundary(A, 6, :jhi) == true
end
