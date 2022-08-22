module BlockHaloArrays

using Base.Threads, Base.Iterators, LinearAlgebra
using ThreadPools, ThreadPinning, NumaAllocators

export BlockHaloArray
export flatten, repartition!

"""
"""
struct BlockHaloArray{AA<:Array{T,N1} where {T,N1},N}
    blocks::Vector{AA}
    blockdims::NTuple{N,Int64}
    global_blockranges::Array{NTuple{N,UnitRange{Int64}},N}
    nhalo::Int64
    loop_limits::Vector{Vector{Int64}}
    globaldims::NTuple{N,Int64}
end

include("partitioning.jl")
include("halo_exchange.jl")

"""

# Arguments
 - `dims::NTuple{N,Int}`: Array dimensions
 - `nhalo::Integer`: Number of halo entries (equal in all dimensions)

# Keyword Arguments
 - `nblocks::Integer`: Number of blocks to divide the array into; default is nthreads()
 - `T`:: Array number type; default is Float64 
"""
function BlockHaloArray(dims::NTuple{N,Int}, nhalo::Integer, nblocks=nthreads(); T=Float64) where {N}

    if nblocks > nthreads()
        @error "Unable to partition; nblocks > nthreads"
    end

    blocks = Vector{Array{Float64,N}}(undef, nblocks)

    tile_dims = block_layout(nblocks, N) |> Tuple

    block_sizes = split_count.(dims, collect(tile_dims))
    block_ranges = get_block_ranges(dims, nblocks)

    nneighbors = 3^N - 1
    neighbors = Vector{Vector{Int64}}(undef, nblocks)

    CI = CartesianIndices(tile_dims)
    LI = LinearIndices(tile_dims)

    threadid_to_numa_mapping = map_threadid_to_numa()

    for threadid in LI
        block_dim = size.(block_ranges[threadid]) |> Iterators.flatten |> collect |> Tuple
        numa_id = threadid_to_numa_mapping[threadid]
        # blocks[threadid] = zeros(block_dim .+ 2nhalo)
        # blocks[threadid] = Array{T}(undef, block_dim .+ 2nhalo)

        # allocate on the thread's numa node
        blocks[threadid] = Array{T}(numa(numa_id), block_dim .+ 2nhalo)
    end


    loop_limits = Vector{Vector{Int64}}(undef, nblocks)
    for I in LI
        loop_limits[I] = [(first(ax) + nhalo, last(ax) - nhalo)
                          for ax in axes(blocks[I])] |> flatten |> collect
    end
    # for I in CI
    # 	ijk = Tuple(I)
    # 	@show ijk
    # 	neigh = [(l-1, l+1) for l in ijk] |> Tuple
    # 	@show neigh
    # end

    A = BlockHaloArray(blocks, tile_dims, block_ranges, nhalo, loop_limits, dims)#, neighbors, LI)

    # testing NUMA first-touch policy
    @sync for tid in 1:nblocks
        ThreadPools.@tspawnat tid fill!(A.blocks[tid], 0)
    end

    A
end

function BlockHaloArray(A::AbstractArray{T,N}, nhalo::Integer, nblocks=nthreads()) where {T,N}

    dims = size(A)
    A_blocked = BlockHaloArray(dims, nhalo, nblocks, T=T)

    block_ranges = get_block_ranges(dims, nblocks)

    for tid in LinearIndices(block_ranges)
        domain_indices = UnitRange.((axes(A.blocks[tid]) .|> first) .+ nhalo,
            (axes(A.blocks[tid]) .|> last) .- nhalo)

        domain_view = view(A_blocked.blocks[tid], domain_indices...)

        A_view = view(A, block_ranges[tid]...)
        copy!(domain_view, A_view)
    end

    return A_blocked
end

"""Create a dictionary mapping the threadid to the NUMA node."""
function map_threadid_to_numa()
    buf = []
    for (id, node) in enumerate(cpuids_per_numa())
        for tid in node
            push!(buf, (tid + 1, id - 1)) # +/- 1 is due to difference in 0 and 1-based indexing with lscpu
        end
    end
    Dict(buf)
end

import Base.eltype, Base.size, Base.axes

eltype(A::BlockHaloArray) = eltype(first(A.blocks))

size(A::BlockHaloArray) = size.(A.blocks)
size(A::BlockHaloArray, dim) = size.(A.blocks, dim)

axes(A::BlockHaloArray) = axes.(A.blocks)
axes(A::BlockHaloArray, dim) = axes.(A.blocks, dim)

nblocks(A::BlockHaloArray) = length(A.blocks)

end
