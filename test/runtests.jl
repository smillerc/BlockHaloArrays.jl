using .Threads
using BlockHaloArrays
using Test

@testset "BlockHaloArrays.jl" begin

    dims = (512, 256)
    nhalo = 2
    A = BlockHaloArray(dims, nhalo; nblocks=nthreads(), T=Float64)

    @test length(A.blocks) == nthreads()
end
