module BlockHaloArrays

using Base.Threads, Base.Iterators, LinearAlgebra
using ThreadPools, ThreadPinning, NumaAllocators

export BlockHaloArray

"""
"""
struct BlockHaloArray{AA<:Array{T,N1} where {T,N1},N}
    blocks::Vector{AA}
    blockdims::NTuple{N,Int64}
    global_blockranges::Array{NTuple{N,UnitRange{Int64}},N}
    nhalo::Int64
    loop_limits::Vector{Vector{Int64}}
end

include("partitioning.jl")
include("halo_exchange.jl")

"""

# Arguments
 - `dims::NTuple{N,Int}`: Array dimensions
 - `nhalo::Int`: Number of halo entries (equal in all dimensions)

# Keyword Arguments
 - `nblocks::Int`: Number of blocks to divide the array into; default is nthreads()
 - `T`:: Array number type; default is Float64 
"""
function BlockHaloArray(dims::NTuple{N,Int}, nhalo::Int; nblocks=nthreads(), T=Float64) where {N}

    blocks = Vector{Array{Float64,N}}(undef, nblocks)

    tile_dims = num_tiles(nblocks, N) |> Tuple

    block_sizes = split_count.(dims, collect(tile_dims))
    block_ranges = get_block_ranges(dims, nblocks)

    nneighbors = 3^N - 1
    neighbors = Vector{Vector{Int64}}(undef, nblocks)

    CI = CartesianIndices(tile_dims)
    LI = LinearIndices(tile_dims)

    threadid_to_numa_mapping = map_threadid_to_numa()

    for threadid in LI
        block_dim = size.(block_ranges[threadid]) |> flatten |> collect |> Tuple
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

    A = BlockHaloArray(blocks, tile_dims, block_ranges, nhalo, loop_limits)#, neighbors, LI)
    _first_touch(A) # for NUMA first-touch policy
    A
end

function _first_touch(A)
    @sync for tid in 1:nthreads()
        t = ThreadPools.@tspawnat tid init_blocks(A, tid)
    end
end

function init_blocks(A, tid)
    fill!(A.blocks[tid], 0)
end


function map_threadid_to_numa()
    buf = []
    for (id, node) in enumerate(cpuids_per_numa())
        for tid in node
            push!(buf, (tid + 1, id - 1)) # +/- 1 is due to difference in 0 and 1-based indexing with lscpu
        end
    end
    Dict(buf)
end

end
