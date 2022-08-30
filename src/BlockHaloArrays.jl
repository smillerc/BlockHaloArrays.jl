module BlockHaloArrays

import Base.eltype, Base.size, Base.axes

using Base.Threads, Base.Iterators, LinearAlgebra
using ThreadPools, ThreadPinning, NumaAllocators
using EllipsisNotation

export BlockHaloArray
export flatten, repartition!, sync_halo!
export domainview
export onboundary



abstract type AbstractBlockHaloArray end

"""
A blocked array structure that stores thread-specific data in blocks. This facilitates a micro domain-decompisition
for shared-memory applications. Each thread operates on it's own block of data. This provides better performance
scaling than multi-threaded loops

# Fields
 - `blocks::Vector{AA}`: 
 - `block_layout::NTuple{D,Int}`: number of blocks along each dimension
 - `global_blockranges::Array{NTuple{D,UnitRange{Int}},D}`: Indexing/ranges of each block from the global perspective
 - `nhalo::Int`: Number of halo regions, e.g. 2 entries along each dimension
 - `loop_limits::Vector{Vector{Int}}`: Looping limits for convienence e.g. `[ilo,ihi,jlo,jhi]`
 - `globaldims::NTuple{D,Int}`: Dimensions of the array if it were a simple `Array{T,D}`, e.g. `(20,20)`

"""
struct BlockHaloArray{T,N,NH,NBL,AA<:Array{T,N}} <: AbstractBlockHaloArray
    blocks::Vector{AA}
    block_layout::NTuple{NBL,Int}
    halodims::NTuple{NH,Int}
    global_blockranges::Array{NTuple{N,UnitRange{Int}}}
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
 - `block_layout::NTuple{D,Int}`: number of blocks along each dimension
 - `global_blockranges::Array{NTuple{D,UnitRange{Int}},D}`: Indexing/ranges of each block from the global perspective
 - `nhalo::Int`: Number of halo regions, e.g. 2 entries along each dimension
 - `loop_limits::Vector{Vector{Int}}`: Looping limits for convienence e.g. `[ilo,ihi,jlo,jhi]`
 - `globaldims::NTuple{D,Int}`: Dimensions of the array if it were a simple `Array{T,D}`, e.g. `(20,20)`
 - `_global_halo_send_buf::Vector{Array{T,D}}`: Buffers used to send across MPI ranks
 - `_global_halo_recv_buf::Vector{Array{T,D}}`: Buffers used to receive across MPI ranks

"""
struct MPIBlockHaloArray{T,N,NH,NBL,AA<:Array{T,N}} <: AbstractBlockHaloArray
    blocks::Vector{AA}
    block_layout::NTuple{NBL,Int}
    halodims::NTuple{NH,Int}
    global_blockranges::Array{NTuple{N,UnitRange{Int}}}
    nhalo::Int
    loop_limits::Vector{Vector{Int}}
    globaldims::NTuple{N,Int}
    neighbor_blocks::Vector{Dict{Symbol,Int}}
    _global_halo_send_buf::Vector{Array{T,NH}}
    _global_halo_recv_buf::Vector{Array{T,NH}}
end

include("partitioning.jl")
include("halo_exchange.jl")

"""
Construct a BlockHaloArray

# Arguments
 - `dims::NTuple{N,Int}`: Array dimensions
 - `nhalo::Integer`: Number of halo entries (equal in all dimensions)

# Keyword Arguments
 - `nblocks::Integer`: Number of blocks to divide the array into; default is nthreads()
 - `T`:: Array number type; default is Float64 
"""
function BlockHaloArray(dims::NTuple{N,Int}, halodims::NTuple{N2,Int}, nhalo::Integer, nblocks=nthreads(); T=Float64, use_numa=true) where {N,N2}

    alldims = Tuple(1:length(dims))
    non_halo_dims = Tuple([i for (i, v) in enumerate(dims) if !(i in halodims)])

    for dim in halodims
        if !(dim in alldims)
            error("Invalid halo dimension: $(dim) not in any of the array axes $alldims")
        end
    end

    if !isempty(non_halo_dims)
        if !(all(non_halo_dims .< halodims))
            error("The axes for halo exchange must be the outmost-most of the array: halodims=$halodims, non-halo dims=$non_halo_dims")
        end
    end

    if any(halodims .> length(dims))
        error("Some (or all) of the given halo_dims $(halodims) are incompatible with the dimensionality of the given array A")
    end

    mismatched_blocks = false
    if nblocks > Threads.nthreads()
        @warn "nblocks ($nblocks) > nthreads ($(nthreads()))"
        mismatched_blocks = true
    end

    blocks = Vector{Array{Float64,N}}(undef, nblocks)
    halo_only_sizes = Tuple([v for (i, v) in enumerate(dims) if i in halodims])
    tile_dims = block_layout(nblocks, length(halodims)) |> Tuple
    halo_only_dims = Tuple([v for (i, v) in enumerate(dims) if i in halodims])
    non_halo_dim_sizes = Tuple([v for (i, v) in enumerate(dims) if !(i in halodims)])

    if halodims == Tuple(1:length(dims))
        block_ranges = get_block_ranges(dims, nblocks)
    else
        block_ranges_halo_only = get_block_ranges(halo_only_dims, nblocks)
        block_ranges = update_block_ranges_with_non_halo_dims(block_ranges_halo_only, dims, halodims)
    end

    block_sizes = [collect(flatten(size.(block))) for block in block_ranges]
    
    # pad the dimensions that will include halo regions
    for block in block_sizes
        for i in eachindex(block)
            if i in halodims
                block[i] += 2nhalo
            end
        end
    end

    for threadid in eachindex(blocks)       
        # allocate on the thread's numa node
        if use_numa
            try
                threadid_to_numa_mapping = map_threadid_to_numa()
                numa_id = threadid_to_numa_mapping[threadid]
                blocks[threadid] = Array{T}(numa(numa_id), block_sizes[threadid]...)
            catch
                if threadid == firstindex(blocks) 
                    @warn "Unable to allocate blocks on the thread-local NUMA node"
                end
                blocks[threadid] = Array{T}(undef, block_sizes[threadid]...)
            end
        else
            blocks[threadid] = Array{T}(undef, block_sizes[threadid]...)
        end
    end

    loop_limits = Vector{Vector{Int}}(undef, nblocks)
    for i in eachindex(loop_limits)
        loop_limits[i] = [(first(ax) + nhalo, last(ax) - nhalo)
                          for ax in axes(blocks[i])] |> flatten |> collect
    end

    neighbors = get_neighbor_blocks(tile_dims)

    A = BlockHaloArray(blocks, tile_dims, halodims, block_ranges, nhalo, loop_limits, dims, neighbors)

    # # testing NUMA first-touch policy
    if !mismatched_blocks
        if use_numa
            @sync for tid in 1:nblocks
                ThreadPools.@tspawnat tid fill!(A.blocks[tid], 0)
            end
        end
    else
        for tid in 1:nblocks
            fill!(A.blocks[tid], 0)
        end
    end

    A
end

function BlockHaloArray(A::AbstractArray{T,N}, nhalo::Integer, nblocks=nthreads()) where {T,N}
    halodims = Tuple(1:length(size(A)))
    BlockHaloArray(A, halodims, nhalo, nblocks)
end

"""

Construct a BlockHaloArray from a normal Array
"""
function BlockHaloArray(A::AbstractArray{T,N}, halodims::NTuple{N2,Integer}, nhalo::Integer, nblocks=nthreads()) where {T,N,N2}
    dims = size(A)
    A_blocked = BlockHaloArray(dims, halodims, nhalo, nblocks, T=T)
    block_ranges = A_blocked.global_blockranges

    for tid in LinearIndices(block_ranges)
        domain_view = domainview(A_blocked, tid)
        A_view = view(A, block_ranges[tid]...)
        copy!(domain_view, A_view)
    end

    return A_blocked
end

function BlockHaloArray(dims::NTuple{N,Int}, nhalo::Integer, nblocks=nthreads(); T=Float64, use_numa=true) where {N}

    halodims = Tuple(1:length(dims))
    BlockHaloArray(dims, halodims, nhalo, nblocks; T=T, use_numa=use_numa)
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
Determine if the block is on the boundary

# Arguments
 - `A::AbstractBlockHaloArray`: The array in question
 - `blockid::Integer`: The id of the block. This is normally associated with the thread id
 - `boundary::Symbol`: Which boundary are we checking for? Examples include :ilo, :jhi, etc... 
"""
function onboundary(A::AbstractBlockHaloArray, blockid::Integer, boundary::Symbol)
    try
        if A.neighbor_blocks[blockid][boundary] == -1
            return true
        else
            return false
        end
    catch
        return false
    end
end

"""
    domainview(A::BlockHaloArray, blockid) -> SubArray

Get the SubArray view of the domain region of a block in the BlockHaloArray.
"""
function domainview(A::BlockHaloArray, blockid::Integer)
    _, _, lo_dom_start, _ = lo_indices(A.blocks[blockid], A.nhalo)
    _, hi_dom_end, _, _ = hi_indices(A.blocks[blockid], A.nhalo)

    idx_range = UnitRange.(lo_dom_start, hi_dom_end)
    idx_range_vec = collect(idx_range)
    for dim in eachindex(idx_range_vec)
        if !(dim in A.halodims)
            idx_range_vec[dim] = axes(A.blocks[blockid], dim)
        end
    end
    return @views A.blocks[blockid][.., Tuple(idx_range_vec)...]
end

end # module