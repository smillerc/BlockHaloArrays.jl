
function get_1D_halo_regions(dims, nhalo; T=Float64)

    #  7  6  5  |  7  6  5
    #  8     4  |  8     4
    #  1  2  3  |  1  2  3
    #  ---------|---------
    #  7  6  5  |  7  6  5
    #  8     4  |  8     4
    #  1  2  3  |  1  2  3

    return [
        zeros(T, (nhalo)),   # 1
        zeros(T, (nhalo)), # 2
    ]
end

function get_2D_halo_regions(dims, nhalo; T=Float64)

    #  7  6  5  |  7  6  5
    #  8     4  |  8     4
    #  1  2  3  |  1  2  3
    #  ---------|---------
    #  7  6  5  |  7  6  5
    #  8     4  |  8     4
    #  1  2  3  |  1  2  3

    #   j
    #   |
    #   *--- i

    ni, nj = dims
    return [
        zeros(T, (nhalo, nhalo)),   # 1
        zeros(T, (ni, nhalo)), # 2
        zeros(T, (nhalo, nhalo)),   # 3
        zeros(T, (nhalo, nj)), # 4
        zeros(T, (nhalo, nhalo)),   # 5
        zeros(T, (ni, nhalo)), # 6
        zeros(T, (nhalo, nj)), # 7
        zeros(T, (nhalo, nhalo)),   # 8
    ]
end

function get_3D_halo_regions(dims, nhalo; T=Float64)

    #   front   |   middle   |    back
    #  7  6  5  |  16 15 14  |  24 23 22
    #  8  9  4  |  17    13  |  25 26 21
    #  1  2  3  |  10 11 12  |  18 19 20 

    #     k
    #    /
    #   *--- i
    #   |
    #   j

    ni, nj, nk = dims
    return [
        zeros(T, (nhalo, nhalo, nhalo)), # 1 
        zeros(T, (ni, nhalo, nhalo)),    # 2 
        zeros(T, (nhalo, nhalo, nhalo)), # 3 
        zeros(T, (nhalo, nj, nhalo)),    # 4 
        zeros(T, (nhalo, nhalo, nhalo)), # 5 
        zeros(T, (ni, nhalo, nhalo)),    # 6 
        zeros(T, (nhalo, nhalo, nhalo)), # 7 
        zeros(T, (nhalo, nj, nhalo)),    # 8 
        zeros(T, (ni, nj, nhalo)),       # 9 
        zeros(T, (nhalo, nhalo, nk)),    # 10
        zeros(T, (nhalo, nhalo, nk)),    # 11
        zeros(T, (nhalo, nhalo, nk)),    # 12
        zeros(T, (nhalo, nhalo, nk)),    # 13
        zeros(T, (nhalo, nhalo, nk)),    # 14
        zeros(T, (nhalo, nhalo, nk)),    # 15
        zeros(T, (nhalo, nhalo, nk)),    # 16
        zeros(T, (nhalo, nhalo, nk)),    # 17
        zeros(T, (nhalo, nhalo, nhalo)), # 18
        zeros(T, (ni, nhalo, nhalo)),    # 19
        zeros(T, (nhalo, nhalo, nhalo)), # 20
        zeros(T, (nhalo, nj, nhalo)),    # 21
        zeros(T, (nhalo, nhalo, nhalo)), # 22
        zeros(T, (ni, nhalo, nhalo)),    # 23
        zeros(T, (nhalo, nhalo, nhalo)), # 24
        zeros(T, (nhalo, nj, nhalo)),    # 25
        zeros(T, (ni, nj, nhalo)),       # 26
    ]
end

function halo_exhange_map_1d()
    return Dict(:ilo => :ihi, :ihi => :ilo)
end

function halo_exhange_map_2d()
    return Dict(:ilo => :ihi, :ihi => :ilo, :jlo => :jhi, :jhi => :jlo,
        :ilojlo => :ihijhi, :ihijlo => :ilojhi, :ihijhi => :ilojlo, :ilojhi => :ihijlo)
end

"""
    domain_donor_ranges_2d(block, nhalo) -> Dict

Get the ranges for each donor region in the doman that sends data to neighbor block halo regions
"""
function domain_donor_ranges_2d(block, nhalo)

    # domn == domain area
    _, _, lo_domn_start, lo_domn_end = lo_indices(block, nhalo)
    hi_domn_start, hi_domn_end, _, _ = hi_indices(block, nhalo)

    ilo_domn_start, jlo_domn_start = lo_domn_start
    ilo_domn_end, jlo_domn_end = lo_domn_end
    ihi_domn_start, jhi_domn_start = hi_domn_start
    ihi_domn_end, jhi_domn_end = hi_domn_end

    # key -> neighbor id, value -> array index ranges
    return Dict(
        :ilo => (ilo_domn_start:ilo_domn_end, jlo_domn_start:jhi_domn_end),
        :jlo => (ilo_domn_start:ihi_domn_end, jlo_domn_start:jlo_domn_end),
        :ihi => (ihi_domn_start:ihi_domn_end, jlo_domn_start:jhi_domn_end),
        :jhi => (ilo_domn_start:ihi_domn_end, jhi_domn_start:jhi_domn_end),
        :ilojlo => (ilo_domn_start:ilo_domn_end, jlo_domn_start:jlo_domn_end),
        :ihijlo => (ihi_domn_start:ihi_domn_end, jlo_domn_start:jlo_domn_end),
        :ilojhi => (ilo_domn_start:ilo_domn_end, jhi_domn_start:jhi_domn_end),
        :ihijhi => (ihi_domn_start:ihi_domn_end, jhi_domn_start:jhi_domn_end),
    )
end

"""
    halo_reciever_ranges_2d(block, nhalo) -> Dict

Get the ranges for each halo region that recieves data from neighbor blocks
"""
function halo_reciever_ranges_2d(block, nhalo)

    # domn == domain area
    lo_halo_start, lo_halo_end, lo_domn_start, _ = lo_indices(block, nhalo)
    _, hi_domn_end, hi_halo_start, hi_halo_end = hi_indices(block, nhalo)

    ilo_halo_start, jlo_halo_start = lo_halo_start
    ilo_halo_end, jlo_halo_end = lo_halo_end
    ilo_domn_start, jlo_domn_start = lo_domn_start
    ihi_domn_end, jhi_domn_end = hi_domn_end
    ihi_halo_start, jhi_halo_start = hi_halo_start
    ihi_halo_end, jhi_halo_end = hi_halo_end

    # key -> neighbor id, value -> array index ranges
    return Dict(
        :ilo => (ilo_halo_start:ilo_halo_end, jlo_domn_start:jhi_domn_end),
        :ihi => (ihi_halo_start:ihi_halo_end, jlo_domn_start:jhi_domn_end),
        :jlo => (ilo_domn_start:ihi_domn_end, jlo_halo_start:jlo_halo_end),
        :jhi => (ilo_domn_start:ihi_domn_end, jhi_halo_start:jhi_halo_end),
        :ilojlo => (ilo_halo_start:ilo_halo_end, jlo_halo_start:jlo_halo_end),
        :ilojhi => (ilo_halo_start:ilo_halo_end, jhi_halo_start:jhi_halo_end),
        :ihijlo => (ihi_halo_start:ihi_halo_end, jlo_halo_start:jlo_halo_end),
        :ihijhi => (ihi_halo_start:ihi_halo_end, jhi_halo_start:jhi_halo_end),
    )
end

function halo_exhange_map_3d()
    #   front   |   middle   |    back
    #  7  6  5  |  16 15 14  |  24 23 22
    #  8  9  4  |  17    13  |  25 26 21
    #  1  2  3  |  10 11 12  |  18 19 20 

    return Dict(
        1 => 22,
        2 => 23,
        3 => 24,
        4 => 25,
        5 => 18,
        6 => 19,
        7 => 20,
        8 => 21,
        9 => 26,
        10 => 14,
        11 => 15,
        12 => 16,
        13 => 17,
        14 => 10,
        15 => 11,
        16 => 12,
        17 => 13,
        18 => 5,
        19 => 6,
        20 => 7,
        21 => 8,
        22 => 1,
        23 => 2,
        24 => 3,
        25 => 4,
        26 => 9,
    )
end

function halo_exchange(A::BlockHaloArray)
end

function _halo_exchange_1d(A::BlockHaloArray)
end

function _halo_exchange_2d(A::BlockHaloArray)
    # mapping = halo_exchange_map_2d()
    @sync for tid in 1:length(A.blocks)
        ThreadPools.@tspawnat tid neighbor_exhange(A, tid)
    end
end

function update_halo!(A::BlockHaloArray)
    _halo_exchange_2d(A::BlockHaloArray)
end

@inline function valid_neighbor(id::Integer)
    if id > 0
        return true
    else
        return false
    end
end

"""Copy the current block's domain data to the neighbor block's halo data"""
function neighbor_exhange(A::BlockHaloArray, block_id)
    current_block = @views A.blocks[block_id]

    domain_ranges = domain_donor_ranges_2d(current_block, A.nhalo)
    halo_ranges = halo_reciever_ranges_2d(current_block, A.nhalo)

    for (dom_id, halo_id) in halo_exhange_map_2d()
        neighbor_id = A.neighbor_blocks[block_id][halo_id]
        if valid_neighbor(neighbor_id)
            donor = @views current_block[.., domain_ranges[dom_id]...]
            halo = @views A.blocks[neighbor_id][.., halo_ranges[halo_id]...]
            copy!(halo, donor)
        end
    end
    return nothing
end

function _halo_exchange_3d(A::BlockHaloArray)
end

function halo_exchange(A::MPIBlockHaloArray)
end

"""Get the high-side indices for an array `A` given a number of halo entries `nhalo`"""
function hi_indices(A::AbstractArray, nhalo::Integer)
    hi_halo_end = last.(axes(A))               # end index of the halo region
    hi_halo_start = hi_halo_end .- nhalo .+ 1  # start index of the halo region
    hi_domn_end = hi_halo_start .- 1           # end index of the inner domain
    hi_domn_start = hi_domn_end .- nhalo .+ 1  # start index of the inner domain

    return (hi_domn_start, hi_domn_end, hi_halo_start, hi_halo_end)
end

"""Get the low-side indices for an array `A` given a number of halo entries `nhalo`"""
function lo_indices(A::AbstractArray, nhalo::Integer)
    lo_halo_start = first.(axes(A))            # start index of the halo region
    lo_halo_end = lo_halo_start .+ nhalo .- 1  # end index of the halo region
    lo_domn_start = lo_halo_end .+ 1           # start index of the inner domain
    lo_domn_end = lo_domn_start .+ nhalo .- 1  # end index of the inner domain

    return (lo_halo_start, lo_halo_end, lo_domn_start, lo_domn_end)
end

function get_2d_neighbor_blocks(tile_dims)
    nblocks = prod(tile_dims)
    nneighbors = 8

    CI = CartesianIndices(tile_dims)
    LI = LinearIndices(tile_dims)
    block_neighbors = Vector{Dict{Symbol,Int}}(undef, nblocks)

    neighbor_block_sym = [:ilojlo, :jlo, :ihijlo, :ihi, :ihijhi, :jhi, :ilojhi, :ilo]
    for block_idx in CI
        i, j = Tuple(block_idx)

        # set up the block neighbors based on cartesian indexing
        neighbor_indices = [
            (i - 1, j - 1),
            (i, j - 1),
            (i + 1, j - 1),
            (i + 1, j),
            (i + 1, j + 1),
            (i, j + 1),
            (i - 1, j + 1),
            (i - 1, j)
        ]

        block_neighbor_set = Vector{Tuple{Symbol,Int}}(undef, nneighbors)
        
        # if the index is out of bounds, it's invalid. Otherwise, save the 1D index of the block 
        # that is the proper neighbor
        for (neighbor_id, idx) in enumerate(neighbor_indices)

            neighbor_symbol = neighbor_block_sym[neighbor_id]
            valid_neighbor = checkbounds(Bool, CI, idx...)
            if valid_neighbor
                neighbor_block = LI[idx...]
            else
                neighbor_block = -1
            end
            block_neighbor_set[neighbor_id] = (neighbor_symbol, neighbor_block)
        end
        block_neighbors[LI[i, j]] = Dict(block_neighbor_set)
    end

    return block_neighbors
end