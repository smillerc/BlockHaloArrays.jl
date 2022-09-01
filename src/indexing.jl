
"""
Get the cartesian block index within the block layout. 
For example, in a (2,3) block tiling scheme, based on a given
global index `global_i`, and the dimension-specific cummulative `block_sizes`,
find the 2D index. The `dim` argument is used to ignore dimensions
that are not include in the `halodims`.
"""
function get_block_idx(global_i, block_sizes, dim, halodims)
    block_i = 1
    if !(dim in halodims)
        return 0
    end
    for (i, val) in enumerate(block_sizes[dim])
        if global_i <= val
            block_i = i - 1
            break
        end
    end

    return block_i + 1
end


"""For a given single dimension, find the block-local index"""
function _get_single_dim_local_index(global_i, block_i, nhalo, block_size)
    if block_i == 1
        return global_i + nhalo
    else
        dblock = block_size[block_i] - block_size[block_i-1]
        return global_i - dblock + nhalo
    end
end

"""
Get the block-local index based on the given global index. The `block_i` value
is the dimension-specific block cartesian index. Dimensions not included in the 
`halodims` are a special case where the global index is the same as the local index.
"""
function get_local_idx(global_idx, block_i, nhalo, block_sizes, dim, halodims)
    if !(dim in halodims)
        return global_idx[dim]
    else
        hdim = dim - (length(global_idx) - length(halodims))
        li = _get_single_dim_local_index(global_idx[dim],
            block_i[hdim], nhalo,
            block_sizes[dim])
        return li
    end
end

"""
Overload the getindex, or [], method for a BlockHaloArray. This is a 'flat' iterator of sorts, 
since the actual data within the BlockHaloArray is a series of separate blocks. Iterating through
the array in this manner will be slower due to the additional arithmetic need to find the global
to local index conversion for each block.
"""
function Base.getindex(A::BlockHaloArray{T,N,NH,NBL,AA}, global_idx::Vararg{Int,N}) where {T,N,NH,NBL,AA}
    block_idx = ntuple(i ->
            get_block_idx(global_idx[A.halodims[i]],
                A._cummulative_blocksize_per_dim,
                A.halodims[i],
                A.halodims),
        length(A.halodims))

    local_idx = ntuple(i ->
            get_local_idx(global_idx,
                block_idx,
                A.nhalo,
                A._cummulative_blocksize_per_dim,
                i,
                A.halodims),
        length(A.globaldims))

    block_id::Int = A._linear_indices[block_idx...]
    return A.blocks[block_id][local_idx...]
end

function Base.getindex(A::BlockHaloArray, block_idx::Int)
    return A.blocks[block_idx]
end

"""
Overload the setindex, or A[I] = ... , method for a BlockHaloArray. This is a 'flat' iterator of sorts, 
since the actual data within the BlockHaloArray is a series of separate blocks. Iterating through
the array in this manner will be slower due to the additional arithmetic need to find the global
to local index conversion for each block.
"""
function Base.setindex!(A::BlockHaloArray, v, global_idx::Vararg{Int,N}) where {N}
    block_idx = ntuple(i ->
            get_block_idx(global_idx[A.halodims[i]],
                A._cummulative_blocksize_per_dim,
                A.halodims[i],
                A.halodims),
        length(A.halodims))

    local_idx = ntuple(i ->
            get_local_idx(global_idx,
                block_idx,
                A.nhalo,
                A._cummulative_blocksize_per_dim,
                i,
                A.halodims),
        length(A.globaldims))

    block_id::Int = A._linear_indices[block_idx...]
    
    setindex!(A.blocks[block_id], v, local_idx...)
end

function Base.setindex!(A::BlockHaloArray, v, block_idx::Int)
    setindex!(A.blocks, v, block_idx)
end


# A = BlockHaloArray((4, 50, 50), (2, 3), 2, 2; T=Float64);
# A.block_layout
# for blockid in eachindex(A.blocks)
#     fill!(A.blocks[blockid], blockid)
# end


# global_idx = (2, 26, 50)

# @code_warntype A[2, 50, 50]

# # @code_warntype A[51]
# A[1, 2, 50]

# @allocated A[1, 2, 50]
