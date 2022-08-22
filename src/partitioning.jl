
function num_tiles(ntiles, dim)
	if dim == 1
		return [ntiles]
	elseif dim == 2
		return num_2d_tiles(ntiles)
	elseif dim == 3
		return num_3d_tiles(ntiles)
	else
		return dim
	end
end

function split_count(N::Integer, n::Integer)
    q,r = divrem(N, n)
    return [i <= r ? q+1 : q for i = 1:n]
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

"""Returns the optimal number of tiles in (i,j) given total number of tiles n"""
function num_2d_tiles(n)
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
    num_2d_tiles = [dim1[1], dim2[1]]
    for i in 2:length(dim1)
        n1 = norm([dim1[i], dim2[i]] .- sqrt(n))
        n2 = norm(num_2d_tiles .- sqrt(n))
        if n1 < n2 
            num_2d_tiles = [dim1[i], dim2[i]]
        end
    end
    return num_2d_tiles
end

"""Returns the optimal number of tiles in (i,j,k) given total number of tiles n"""
function num_3d_tiles(n)
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
    num_3d_tiles = [dim1[1], dim2[1], dim3[1]]
    for i in 2:length(dim1)
        n1 = norm([dim1[i], dim2[i], dim3[i]] .- sqrt(n))
        n2 = norm(num_3d_tiles .- sqrt(n))
        if n1 < n2 
            num_3d_tiles = [dim1[i], dim2[i], dim3[i]]
        end
    end
    return num_3d_tiles
end

function get_block_ranges(Adims, nblocks)

	ndims = length(Adims)
	tile_dims = num_tiles(nblocks, ndims) |> Tuple
	block_sizes = split_count.(Adims, collect(tile_dims))
	
	end_idx = cumsum.(block_sizes)
	start_idx = copy(end_idx)
	for i in eachindex(end_idx)
		start_idx[i] = end_idx[i] - block_sizes[i] .+ 1
	end

	blks = Array{NTuple{ndims, UnitRange{Int64}}}(undef, tile_dims)
	
	for block_I in CartesianIndices(tile_dims)
		blks[block_I] = [start_idx[dim][x]:end_idx[dim][x] 
			for (dim, x) in enumerate(Tuple(block_I))] |> Tuple
	end
	blks	
end