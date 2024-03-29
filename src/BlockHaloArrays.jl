module BlockHaloArrays

using Base.Threads, Base.Iterators, LinearAlgebra
import Base.eltype, Base.size, Base.axes, Base.copy!
import Base.eachindex, Base.first, Base.firstindex, Base.last, Base.lastindex

using NumaAllocators
using EllipsisNotation
using OffsetArrays

export BlockHaloArray
export flatten, repartition!, updatehalo!, updateblockhalo!
export globalindices, nblocks
export domainview
export onboundary
export donorview, haloview
export blockindex, localindex
export @tspawnat
export copy_domain!, copy_halo!

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
struct BlockHaloArray{NB,T,N,NH,NBL,AA<:Array{T,N},SA,NL} <: AbstractBlockHaloArray
    blocks::NTuple{NB,AA}
    block_layout::NTuple{NBL,Int}
    halodims::NTuple{NH,Int}
    global_blockranges::Vector{NTuple{N,UnitRange{Int}}}
    nhalo::Int
    loop_limits::NTuple{NB, NTuple{NL, Int}}
    globaldims::NTuple{N,Int}
    neighbor_blocks::Vector{Dict{Symbol,Int}}

    # private vars
    _cummulative_blocksize_per_dim::OffsetVector{Vector{Int64},Vector{Vector{Int64}}}
    _linear_indices::LinearIndices{NBL,NTuple{NBL,Base.OneTo{Int}}}
    _donor_views::Vector{Dict{Symbol,SA}}
    _halo_views::Vector{Dict{Symbol,SA}}
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
    global_blockranges::Vector{NTuple{N,UnitRange{Int}}}
    nhalo::Int
    loop_limits::NTuple{NBL, NTuple{4, Int}}
    globaldims::NTuple{N,Int}
    neighbor_blocks::Vector{Dict{Symbol,Int}}
    _global_halo_send_buf::Vector{Array{T,NH}}
    _global_halo_recv_buf::Vector{Array{T,NH}}
end

include("spawn.jl")
include("partitioning.jl")
include("halo_exchange.jl")
include("indexing.jl")

"""
Construct a BlockHaloArray

# Arguments
 - `dims::NTuple{N,Int}`: Array dimensions
 - `nhalo::Integer`: Number of halo entries (equal in all dimensions)

# Keyword Arguments
 - `nblocks::Integer`: Number of blocks to divide the array into; default is nthreads()
 - `T`:: Array number type; default is Float64
"""
function BlockHaloArray(dims::NTuple{N,Int}, halodims::NTuple{N2,Int}, nhalo::Integer, nblocks=Base.Threads.nthreads(); T=Float64, use_numa=true, tile_dims=nothing) where {N,N2}

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

    blocks = Vector{Array{T,N}}(undef, nblocks)
    halo_only_sizes = Tuple([v for (i, v) in enumerate(dims) if i in halodims])

    if tile_dims === nothing
        tile_dims = block_layout(nblocks, length(halodims)) |> Tuple
    else
        if prod(tile_dims) != nblocks
            error("Invalid tile_dims; the number of blocks is not consistent")
        end
    end
    halo_only_dims = Tuple([v for (i, v) in enumerate(dims) if i in halodims])
    nhalodims = length(halo_only_dims)
    non_halo_dim_sizes = Tuple([v for (i, v) in enumerate(dims) if !(i in halodims)])

    if halodims == Tuple(1:length(dims))
        block_ranges = get_block_ranges(dims, tile_dims)
    else
        block_ranges_halo_only = get_block_ranges(halo_only_dims, tile_dims)
        block_ranges = update_block_ranges_with_non_halo_dims(block_ranges_halo_only, dims, halodims)
    end

    block_sizes = [collect(flatten(size.(block))) for block in block_ranges]

    # this is the cumulative block size along each dimension, indexed via [tile_dimension][1:blocks_in_tile_dim]
    cummulative_blocksize_per_dim = OffsetArray(
        cumsum.(collect(split_count.(dims[[i for i in halodims]], tile_dims))),
        UnitRange(first(halodims), last(halodims))
    )

    LI = LinearIndices(tile_dims) # used to convert the cartesian block id into the linear one (for the get/set index methods)

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
                # if threadid == firstindex(blocks)
                #     @warn "Unable to allocate blocks on the thread-local NUMA node"
                # end
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

    loop_lims = Tuple(Tuple.(loop_limits))
    neighbors = get_neighbor_blocks(tile_dims)

    # Pass undef'ed vector to the constructor, since the views
    # need to be done _after_ the type is created.
    _halo_views = gen_halo_views(blocks, nhalo, halodims)
    HaloViewTypes = typeof(first(_halo_views)[:ilo])
    halo_views = Vector{Dict{Symbol, HaloViewTypes}}(undef, nblocks)

    _donor_views = gen_donor_views(blocks, nhalo, halodims)
    DonorViewTypes = typeof(first(_donor_views)[:ilo])
    donor_views = Vector{Dict{Symbol, DonorViewTypes}}(undef, nblocks)

    A = BlockHaloArray(tuple(blocks...), tile_dims, halodims, vec(block_ranges), nhalo,
        loop_lims, dims, neighbors,
        cummulative_blocksize_per_dim, LI,
        donor_views, halo_views)

    # Get the halo reciever and domain donor views of each block.
    # These views are used during the halo sync process
    halo_views = gen_halo_views(A)
    donor_views = gen_donor_views(A)

    for i in eachindex(A.blocks)
        A._halo_views[i] = halo_views[i]
        A._donor_views[i] = donor_views[i]
    end

    # # testing NUMA first-touch policy
    if !mismatched_blocks
        if use_numa
            @sync for tid in 1:nblocks
                @tspawnat tid fill!(A.blocks[tid], 0)
            end
        end
    else
        for tid in 1:nblocks
            fill!(A.blocks[tid], 0)
        end
    end

    return A
end

function BlockHaloArray(A::AbstractArray{T,N}, nhalo::Integer, nblocks=Base.Threads.nthreads(); use_numa=true, tile_dims=nothing) where {T,N}
    halodims = Tuple(1:length(size(A)))
    BlockHaloArray(A, halodims, nhalo, nblocks; use_numa=use_numa, tile_dims=tile_dims)
end

"""
Construct a BlockHaloArray from a normal Array
"""
function BlockHaloArray(A::AbstractArray{T,N}, halodims::NTuple{N2,Integer}, nhalo::Integer, nblocks=Base.Threads.nthreads(); use_numa=true, tile_dims=nothing) where {T,N,N2}
    dims = size(A)
    A_blocked = BlockHaloArray(dims, halodims, nhalo, nblocks; T=T, use_numa=use_numa, tile_dims=tile_dims)
    block_ranges = A_blocked.global_blockranges

    for tid in LinearIndices(block_ranges)
        domain_view = domainview(A_blocked, tid)
        A_view = view(A, block_ranges[tid]...)
        copy!(domain_view, A_view)
    end

    return A_blocked
end

function BlockHaloArray(dims::NTuple{N,Int}, nhalo::Integer, nblocks=Base.Threads.nthreads(); T=Float64, use_numa=true, tile_dims=nothing) where {N}

    halodims = Tuple(1:length(dims))
    BlockHaloArray(dims, halodims, nhalo, nblocks; T=T, use_numa=use_numa, tile_dims=tile_dims)
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

size(A::AbstractBlockHaloArray) = A.globaldims
size(A::AbstractBlockHaloArray, dim) = size.(A.blocks, dim)

axes(A::AbstractBlockHaloArray) = axes.(A.blocks)
axes(A::AbstractBlockHaloArray, dim) = axes.(A.blocks, dim)

first(A::AbstractBlockHaloArray) = first(A.blocks)
firstindex(A::AbstractBlockHaloArray) = firstindex(A.blocks)

last(A::AbstractBlockHaloArray) = last(A.blocks)
lastindex(A::AbstractBlockHaloArray) = lastindex(A.blocks)

eachindex(A::AbstractBlockHaloArray) = eachindex(A.blocks)

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
        if A.neighbor_blocks[blockid][boundary] < 0
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
function domainview(A::BlockHaloArray, blockid::Integer, offset)

    if blockid < 1 || blockid > nblocks(A)
        error("Invalid blockid, must be 1 <= blockid <= $(nblocks(A))")
    end

    if offset > A.nhalo
        error("offset must be <= nhalo")
    end

    if offset < 0
        error("offset must be > 0")
    end

    function f(i, halodims, offset)
        if i in halodims
            return offset
        else
            return 0
        end
    end
    idx_offset = ntuple(i -> f(i, A.halodims, offset), length(A.globaldims))

    _, _, lo_dom_start, _ = lo_indices(A.blocks[blockid], A.nhalo, A.halodims)
    _, hi_dom_end, _, _ = hi_indices(A.blocks[blockid], A.nhalo, A.halodims)

    lo_dom_start = lo_dom_start .- idx_offset
    hi_dom_end = hi_dom_end .+ idx_offset

    idx_range = UnitRange.(lo_dom_start, hi_dom_end)
    return @views A.blocks[blockid][.., idx_range...]
end

domainview(A::BlockHaloArray, blockid::Integer) = domainview(A, blockid, 0)

"""Get the halo region at a particular location, e.g. `:ilo` for block `blockid`"""
haloview(A::BlockHaloArray, blockid::Integer, location::Symbol) = A._halo_views[blockid][location]

"""Get the domain donor region at a particular location, e.g. `:ilo` for block `blockid`"""
donorview(A::BlockHaloArray, blockid::Integer, location::Symbol) = A._donor_views[blockid][location]


"""
    copy!(dst, src) -> dst

Copy from an AbstractArray into a BlockHaloArray. The global dimensions of the
BlockHaloArray must be the same as the AbstractArray
"""
function copy!(BHA::BlockHaloArray, AA::AbstractArray)

    if eltype(AA) != eltype(BHA)
        @warn "Mismatched AbstractArray and BlockHaloArray eltypes"
    end

    if size(AA) != BHA.globaldims
        error("Mismatched array dimensions: BlockHaloArray.globaldims $(BHA.globaldims) != size(AbstractArray) $(size(AA))")
    end

    @sync for block_id in 1:nblocks(BHA)
        @tspawnat block_id _array_to_block!(BHA, AA, block_id)
    end
end

function copy!(BHA1::BlockHaloArray, BHA2::BlockHaloArray)

    if eltype(BHA1) != eltype(BHA2)
        @warn "Mismatched AbstractArray and BlockHaloArray eltypes"
    end

    if size(BHA1) != size(BHA2)
        error("Mismatched array dimensions: BlockHaloArray.globaldims $(BHA.globaldims) != size(AbstractArray) $(size(AA))")
    end

    @sync for block_id in 1:nblocks(BHA1)
        @tspawnat block_id begin
            copy!(BHA1[block_id], BHA2[block_id])
        end
    end
end

function copy_domain!(BHA1::BlockHaloArray, BHA2::BlockHaloArray)

    if eltype(BHA1) != eltype(BHA2)
        @warn "Mismatched AbstractArray and BlockHaloArray eltypes"
    end

    if size(BHA1) != size(BHA2)
        error("Mismatched array dimensions: BlockHaloArray.globaldims $(BHA.globaldims) != size(AbstractArray) $(size(AA))")
    end

    @sync for block_id in 1:nblocks(BHA1)
        @tspawnat block_id begin
            B1 = domainview(BHA1, block_id)
            B2 = domainview(BHA2, block_id)
            copy!(B1, B2)
        end
    end
end

function copy_halo!(dst::BlockHaloArray, src::BlockHaloArray)
    if size(src) != size(dst)
        error("Mismatched array dimensions")
    end

    dst_hv = dst._halo_views
    src_hv = src._halo_views
    for blk in 1:nblocks(dst)
        for loc in keys(src_hv[blk])
            copy!(dst_hv[blk][loc], src_hv[blk][loc])
        end
    end

end

function _array_to_block!(BHA::BlockHaloArray, AA::AbstractArray, block_id::Int)
    dv = domainview(BHA, block_id)
    av = view(AA, BHA.global_blockranges[block_id]...)
    copy!(dv, av) # update the block
end

"""
    copy!(dst, src) -> dst

Copy from a BlockHaloArray into an AbstractArray. The global dimensions of the
BlockHaloArray must be the same as the AbstractArray
"""
function copy!(AA::AbstractArray, BHA::BlockHaloArray)
    if size(AA) != BHA.globaldims
        error("Mismatched array dimensions: BlockHaloArray.globaldims $(BHA.globaldims) != size(AbstractArray) $(size(AA))")
    end

    @sync for block_id in 1:nblocks(BHA)
        @tspawnat block_id _block_to_array!(AA, BHA, block_id)
    end
end

function _block_to_array!(AA::AbstractArray, BHA::BlockHaloArray, block_id::Int)
    dv = domainview(BHA, block_id)
    av = view(AA, BHA.global_blockranges[block_id]...)
    copy!(av, dv) # update the abstract array
end

end # module
