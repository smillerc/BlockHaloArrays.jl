
import Base.Iterators.flatten

"""
    block_layout(nblocks, N)

Determine the block layout based on the number of total blocks `nblocks`
and the dimensionality `N` of the domain

# Arguments
 - `nblocks::Integer`
 - `N::Integer`
"""
function block_layout(nblocks::Integer, N::Integer)
    @assert(1 <= N <= 3)
    if N == 1
        return [nblocks]
    elseif N == 2
        return _block_layout_2d(nblocks)
    elseif N == 3
        return _block_layout_3d(nblocks)
    end
end

"""
    split_count(N::Integer, n::Integer)

Return a vector of `n` integers which are approximately equally sized and sum to `N`.
This borrows from https://juliaparallel.org/MPI.jl/latest/examples/06-scatterv/
"""
function split_count(N::Integer, n::Integer)
    q, r = divrem(N, n)
    return [i <= r ? q + 1 : q for i = 1:n]
end

"""Return all common denominators of n"""
function denominators(n::Integer)
    denominators = Vector{Int}(undef, 0)
    for i in 1:n
        if mod(n, i) == 0
            push!(denominators, i)
        end
    end
    return denominators
end

"""Find the optimal block layout in 2D given total number of tiles `n`"""
function _block_layout_2d(n)
    # find all common denominators of the total number of images
    denoms = denominators(n)

    # find all combinations of common denominators
    # whose product equals the total number of images
    dim1 = Vector{Int}(undef, 0)
    dim2 = Vector{Int}(undef, 0)
    for j in eachindex(denoms)
        for i in eachindex(denoms)
            if denoms[i] * denoms[j] == n
                push!(dim1, denoms[i])
                push!(dim2, denoms[j])
            end
        end
    end
    # pick the set of common denominators with the minimal norm
    # between two elements -- rectangle closest to a square
    n_tiles = [dim1[1], dim2[1]]
    for i in 2:length(dim1)
        n1 = norm([dim1[i], dim2[i]] .- sqrt(n))
        n2 = norm(n_tiles .- sqrt(n))
        if n1 < n2
            n_tiles = [dim1[i], dim2[i]]
        end
    end
    return n_tiles
end

"""Find the optimal block layout in 3D given total number of tiles `n`"""
function _block_layout_3d(n)
    # find all common denominators of the total number of images
    denoms = denominators(n)

    # find all combinations of common denominators
    # whose product equals the total number of images
    dim1 = Vector{Int}(undef, 0)
    dim2 = Vector{Int}(undef, 0)
    dim3 = Vector{Int}(undef, 0)
    for k in eachindex(denoms)
        for j in eachindex(denoms)
            for i in eachindex(denoms)
                if denoms[i] * denoms[j] * denoms[k] == n
                    push!(dim1, denoms[i])
                    push!(dim2, denoms[j])
                    push!(dim3, denoms[k])
                end
            end
        end
    end
    # pick the set of common denominators with the minimal norm
    # between two elements -- rectangle closest to a square
    n_tiles = [dim1[1], dim2[1], dim3[1]]
    for i in 2:length(dim1)
        n1 = norm([dim1[i], dim2[i], dim3[i]] .- sqrt(n))
        n2 = norm(n_tiles .- sqrt(n))
        if n1 < n2
            n_tiles = [dim1[i], dim2[i], dim3[i]]
        end
    end
    return n_tiles
end


"""
    get_block_ranges(dims::NTuple{N, Int}, nblocks::Integer) -> Array{NTuple{ndims,UnitRange{Int}}}

Get the ranges of each block (in the global space)

# Arguments
 - `dims::NTuple{N, Int}`: Dimensions of the flat array to split into blocks
 - `nblocks::Integer`: Number of total blocks
"""
function get_block_ranges(dims::NTuple{N,Int}, nblocks::Integer) where {N}
    ndims = length(dims)
    tile_dims = block_layout(nblocks, ndims) |> Tuple
    block_sizes = split_count.(dims, collect(tile_dims))

    end_idx = cumsum.(block_sizes)
    start_idx = copy(end_idx)
    for i in eachindex(end_idx)
        start_idx[i] = end_idx[i] - block_sizes[i] .+ 1
    end

    blks = Array{NTuple{ndims,UnitRange{Int}}}(undef, tile_dims)

    for block_I in CartesianIndices(tile_dims)
        blks[block_I] = [start_idx[dim][x]:end_idx[dim][x]
                         for (dim, x) in enumerate(Tuple(block_I))] |> Tuple
    end

    blks
end


"""
    repartition(A::AbstractBlockHaloArray, nblocks) -> BlockHaloArray

Repartition the BlockHaloArray into a different block layout
"""
function repartition!(A::AbstractBlockHaloArray, nblocks::Integer)

    if nblocks == 1
        @warn "Trying to repartition!() into 1 block is prohibited"
        return A
    elseif nblocks == length(A.blocks)
        @warn "New block layout is the same as the original"
        return A
    end

    Aflat = flatten(A)
    return BlockHaloArray(Aflat, A.nhalo, nblocks)
end

"""
    flatten(A::AbstractBlockHaloArray) -> Array

Return a flattened version of a BlockHaloArray. This is a copy, since a view
of the current block structure isn't possible.
"""
function flatten(A::AbstractBlockHaloArray)

    A_flat = zeros(eltype(A), A.globaldims)
    block_ranges = A.global_blockranges

    for tid in LinearIndices(block_ranges)
        domain_indices = UnitRange.((axes(A.blocks[tid]) .|> first) .+ A.nhalo,
            (axes(A.blocks[tid]) .|> last) .- A.nhalo)

        domain_view = view(A.blocks[tid], domain_indices...)
        A_flat_view = view(A_flat, block_ranges[tid]...)
        copy!(A_flat_view, domain_view)
    end

    return A_flat
end


function globalsize(A::AbstractBlockHaloArray, nohalo=true)
    ndims = length(size(first(A.blocks)))
    global_dims = zeros(Int, ndims)

    for block in A.blocks
        for dim in 1:ndims
            global_dims[dim] += size(block, dim)
        end
    end

    if nohalo
        global_dims .-= length(A.blocks) * ndims * A.nhalo
    end

    global_dims |> Tuple
end