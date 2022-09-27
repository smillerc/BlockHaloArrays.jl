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
