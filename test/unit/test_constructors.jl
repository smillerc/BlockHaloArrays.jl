
@testitem "Copy Constructor" begin
    dims = (50, 50)
    nhalo = 2
    nblocks = 6
    B = rand(dims...)
    A = BlockHaloArray(B, nhalo, nblocks)
end

@testitem "Check tile dims" begin

    dims = (60, 10)
    nhalo = 2
    nblocks = 6
    tiledims = (3,3) # mismatched partitioning

    @test_throws ErrorException("Invalid tile_dims; the number of blocks is not consistent") BlockHaloArray(dims, nhalo, nblocks; tile_dims = tiledims)
end


@testitem "Provide tile dims" begin
    dims = (60, 10)
    nhalo = 2
    nblocks = 6
    tiledims = (6,1)
    A = BlockHaloArray(dims, nhalo, nblocks; tile_dims = tiledims)

    @test A.block_layout == tiledims

    @test A.global_blockranges[1] == (1:10, 1:10)
    @test A.global_blockranges[2] == (11:20, 1:10)
    @test A.global_blockranges[3] == (21:30, 1:10)
    @test A.global_blockranges[4] == (31:40, 1:10)
    @test A.global_blockranges[5] == (41:50, 1:10)
    @test A.global_blockranges[6] == (51:60, 1:10)

    for block in eachindex(A.blocks)
        @test size(A[block]) == (14, 14) # 10 + 2nhalo
    end
end

@testitem "2D Halo, 2D Array, Different Types" begin
    include("common.jl")

    dims = (50, 50)
    nhalo = 2
    nblocks = 6

    for T in [
        Int16, Int32, Int64, Int128,
        Float16, Float32, Float64,
        ComplexF64
    ]
        A = BlockHaloArray(dims, nhalo, nblocks; T=T)

        @test eltype(A) == T
    end
end

@testitem "2D Array, 1D halo dims" begin
    dims = (10, 20)
    halodims = (2,)
    nhalo = 2
    nblocks = 4
    A = BlockHaloArray(dims, nhalo, nblocks)
end

@testitem "2D Array, 2D halo dims" begin
    dims = (10, 20)
    nhalo = 2
    nblocks = 4
    A = BlockHaloArray(dims, nhalo, nblocks)
end

@testitem "3D Array, 1D halo dims" begin
    dims = (10, 20, 30)
    halodims = (2,)
    nhalo = 2
    nblocks = 4
    A = BlockHaloArray(dims, nhalo, nblocks)
end

@testitem "3D Array, 2D halo dims" begin
    dims = (10, 20, 30)
    halodims = (2, 3)
    nhalo = 2
    nblocks = 4
    A = BlockHaloArray(dims, nhalo, nblocks)
end
