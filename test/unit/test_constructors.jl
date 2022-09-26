
@testitem "Copy Constructor" begin
    dims = (50, 50)
    nhalo = 2
    nblocks = 6
    B = rand(dims...)
    A = BlockHaloArray(B, nhalo, nblocks)
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
