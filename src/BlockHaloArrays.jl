module BlockHaloArrays

import Base.eltype, Base.size, Base.axes

using Base.Threads, Base.Iterators, LinearAlgebra
using ThreadPools, ThreadPinning, NumaAllocators
using EllipsisNotation

export BlockHaloArray
export flatten, repartition!, update_halo!
export domainview



abstract type AbstractBlockHaloArray end

"""
A blocked array structure that stores thread-specific data in blocks. This facilitates a micro domain-decompisition
for shared-memory applications. Each thread operates on it's own block of data. This provides better performance
scaling than multi-threaded loops

# Fields
 - `blocks::Vector{AA}`: 
 - `blockdims::NTuple{N,Int}`: dimensions of each block
 - `global_blockranges::Array{NTuple{N,UnitRange{Int}},N}`: Indexing/ranges of each block from the global perspective
 - `nhalo::Int`: Number of halo regions, e.g. 2 entries along each dimension
 - `loop_limits::Vector{Vector{Int}}`: Looping limits for convienence e.g. `[ilo,ihi,jlo,jhi]`
 - `globaldims::NTuple{N,Int}`: Dimensions of the array if it were a simple `Array{T,N}`, e.g. `(20,20)`

"""
struct BlockHaloArray{T, N, AA<:Array{T,N}} <: AbstractBlockHaloArray
    blocks::Vector{AA}
    blockdims::NTuple{N,Int}
    global_blockranges::Array{NTuple{N,UnitRange{Int}},N}
    nhalo::Int
    loop_limits::Vector{Vector{Int}}
    globaldims::NTuple{N,Int}
    neighbor_blocks::Vector{Dict{Symbol,Int}}
end

"""
An MPI-aware BlockHaloArray. The only difference between this and the
plain `BlockHaloArray` is the addition of send/receive buffers that help
MPI communication.

# Fields
 - `blocks::Vector{AA}`: 
 - `blockdims::NTuple{N,Int}`: dimensions of each block
 - `global_blockranges::Array{NTuple{N,UnitRange{Int}},N}`: Indexing/ranges of each block from the global perspective
 - `nhalo::Int`: Number of halo regions, e.g. 2 entries along each dimension
 - `loop_limits::Vector{Vector{Int}}`: Looping limits for convienence e.g. `[ilo,ihi,jlo,jhi]`
 - `globaldims::NTuple{N,Int}`: Dimensions of the array if it were a simple `Array{T,N}`, e.g. `(20,20)`
 - `_global_halo_send_buf::Vector{Array{T,N}}`: Buffers used to send across MPI ranks
 - `_global_halo_recv_buf::Vector{Array{T,N}}`: Buffers used to receive across MPI ranks

"""
struct MPIBlockHaloArray{T, N, AA<:Array{T,N}} <: AbstractBlockHaloArray
    blocks::Vector{AA}
    blockdims::NTuple{N,Int}
    global_blockranges::Array{NTuple{N,UnitRange{Int}},N}
    nhalo::Int
    loop_limits::Vector{Vector{Int}}
    globaldims::NTuple{N,Int}
    neighbor_blocks::Vector{Dict{Symbol,Int}}
    _global_halo_send_buf::Vector{Array{T,N}}
    _global_halo_recv_buf::Vector{Array{T,N}}
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
function BlockHaloArray(dims::NTuple{N,Int}, nhalo::Integer, nblocks=nthreads(); T=Float64, use_numa=true) where {N}

    if nblocks > nthreads()
        @error "Unable to partition; nblocks > nthreads"
    end

    blocks = Vector{Array{Float64,N}}(undef, nblocks)

    tile_dims = block_layout(nblocks, N) |> Tuple

    # block_sizes = split_count.(dims, collect(tile_dims))
    block_ranges = get_block_ranges(dims, nblocks)

    # nneighbors = 3^N - 1

    CI = CartesianIndices(tile_dims)
    LI = LinearIndices(tile_dims)

    threadid_to_numa_mapping = map_threadid_to_numa()

    for threadid in LI
        block_dim = size.(block_ranges[threadid]) |> Iterators.flatten |> collect |> Tuple
        numa_id = threadid_to_numa_mapping[threadid]
        # blocks[threadid] = zeros(block_dim .+ 2nhalo)
        # blocks[threadid] = Array{T}(undef, block_dim .+ 2nhalo)

        # allocate on the thread's numa node
        if use_numa
            blocks[threadid] = Array{T}(numa(numa_id), block_dim .+ 2nhalo)
        else
            blocks[threadid] = Array{T}(undef, block_dim .+ 2nhalo)
        end
    end

    loop_limits = Vector{Vector{Int}}(undef, nblocks)
    for I in LI
        loop_limits[I] = [(first(ax) + nhalo, last(ax) - nhalo)
                          for ax in axes(blocks[I])] |> flatten |> collect
    end

    neighbors = get_2d_neighbor_blocks(tile_dims)
    A = BlockHaloArray(blocks, tile_dims, block_ranges, nhalo, loop_limits, dims, neighbors)

    # # testing NUMA first-touch policy
    if use_numa
        @sync for tid in 1:nblocks
            ThreadPools.@tspawnat tid fill!(A.blocks[tid], 0)
        end
    end

    A
end

function BlockHaloArray(A::AbstractArray{T,N}, nhalo::Integer, nblocks=nthreads()) where {T,N}

    dims = size(A)
    A_blocked = BlockHaloArray(dims, nhalo, nblocks, T=T)

    block_ranges = get_block_ranges(dims, nblocks)

    for tid in LinearIndices(block_ranges)
        domain_indices = UnitRange.((axes(A_blocked.blocks[tid]) .|> first) .+ nhalo,
            (axes(A_blocked.blocks[tid]) .|> last) .- nhalo)

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


eltype(A::AbstractBlockHaloArray) = eltype(first(A.blocks))

size(A::AbstractBlockHaloArray) = size.(A.blocks)
size(A::AbstractBlockHaloArray, dim) = size.(A.blocks, dim)

axes(A::AbstractBlockHaloArray) = axes.(A.blocks)
axes(A::AbstractBlockHaloArray, dim) = axes.(A.blocks, dim)

nblocks(A::AbstractBlockHaloArray) = length(A.blocks)

"""
    domainview(A::BlockHaloArray, blockid) -> SubArray

Get the SubArray view of the domain region of a block in the BlockHaloArray.
"""
function domainview(A::BlockHaloArray, blockid::Integer)
    _, _, lo_dom_start, _ = lo_indices(A.blocks[blockid], A.nhalo)
    _, hi_dom_end, _, _ = hi_indices(A.blocks[blockid], A.nhalo)

    idx_range = UnitRange.(lo_dom_start, hi_dom_end)

    return @views A.blocks[blockid][.., idx_range...]
end

end
