@testitem "halo/donorview" begin
    include("common.jl")
    dims = (50, 50)
    nhalo = 2
    n_blocks = 4
    A = BlockHaloArray(dims, nhalo, n_blocks)

    @test size(donorview(A, 1, :ilo)) == (2, 25)
    @test size(haloview(A, 1, :ilo)) == (2, 25)
    @test size(donorview(A, 1, :jhi)) == (25, 2)
    @test size(haloview(A, 1, :jhi)) == (25, 2)
end

@testitem "domainview" begin
    include("common.jl")
    nhalo = 3
    n_blocks = 4
    dims = (40,40)
    A = BlockHaloArray(dims, nhalo, n_blocks)

    blockid = 1
    @test size(domainview(A, blockid)) == (20, 20)
    @test size(domainview(A, blockid, 1)) == (22, 22)
    @test size(domainview(A, blockid, 2)) == (24, 24)
    @test size(domainview(A, blockid, 3)) == (26, 26)
    @test_throws ErrorException domainview(A, blockid, 4)
    @test_throws ErrorException domainview(A, blockid, -1)
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
    n_blocks = 6

    A = BlockHaloArray(dims, nhalo, n_blocks; T=Int64)

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

@testitem "1D Indexing" begin
    include("common.jl")
    x = rand(50)
    nhalo = 2
    n_blocks = 4
    A = BlockHaloArray(x, nhalo, n_blocks)

    @test all(A[1] .== A.blocks[1])
    @test_nowarn A[4] .= 2.5
end

@testitem "n-D Indexing" begin
    include("common.jl")
    x = rand(40, 40)
    y = rand(4, 40, 40)

    nhalo = 2
    n_blocks = 4
    A = BlockHaloArray(x, nhalo, n_blocks)

    halodims = (2, 3)
    B = BlockHaloArray(y, halodims, nhalo, n_blocks)

    @test A[1] == A.blocks[1]
    @test size(A[1]) == (24, 24)
    @test A[1][3,3] == x[1,1]
    @test_nowarn A[1][3,3] = 1.2

    @test B[1][1,3,3] == y[1,1,1]
    @test B[1][:,3,3] == y[:,1,1]
    @test_nowarn B[1][1,3,3] =6.7
end

@testitem "globalindices" begin

    x = rand(40, 40)
    nhalo = 2
    n_blocks = 4
    A = BlockHaloArray(x, nhalo, n_blocks)

    @test globalindices(A, 1, (3, 3)) == (1, 1)
    @test globalindices(A, 4, (3, 3)) == (21, 21)

    @show globalindices(A, 4, (22, 22))
end

@testitem "copy!" begin
    include("common.jl")

    dims = (10, 10)
    nhalo = 2
    n_blocks = 2
    AA = rand(dims...)
    AA2 = zeros(dims...)

    AA_mismatch = rand(12,10)
    BHA = BlockHaloArray(dims, nhalo, n_blocks)

    @test_throws ErrorException copy!(AA_mismatch, BHA)
    @test_throws ErrorException copy!(BHA, AA_mismatch)

    copy!(BHA, AA) # copy from the abstract array to the block halo array
    @test AA == flatten(BHA)

    copy!(AA2, BHA) # copy from the abstract array to the block halo array
    @test AA2 == flatten(BHA)

end
