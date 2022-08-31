
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

function halo_exhange_map_3d()
    return Dict(
        :ilo => :ihi,
        :ihi => :ilo,
        :jlo => :jhi,
        :jhi => :jlo,
        :klo => :khi,
        :khi => :klo,
        :ilojlo => :ihijhi,
        :ihijlo => :ilojhi,
        :ilojhi => :ihijlo,
        :ihijhi => :ilojlo,
        :ilokhi => :ihiklo,
        :iloklo => :ihikhi,
        :ihikhi => :iloklo,
        :ihiklo => :ilokhi,
        :jhikhi => :jloklo,
        :jlokhi => :jhiklo,
        :jhiklo => :jlokhi,
        :jloklo => :jhikhi,
        :ilojhikhi => :ihijloklo,
        :ihijhikhi => :ilojloklo,
        :ilojlokhi => :ihijhiklo,
        :ihijlokhi => :ilojhiklo,
        :ilojhiklo => :ihijlokhi,
        :ihijhiklo => :ilojlokhi,
        :ilojloklo => :ihijhikhi,
        :ihijloklo => :ilojhikhi,
    )

end
"""
    domain_donor_ranges_1d(block, nhalo) -> Dict

Get the ranges for each donor region in the doman that sends data to neighbor block halo regions
"""
function domain_donor_ranges_1d(block, nhalo, halodims)

    # domn == domain area
    _, _, lo_domn_start, lo_domn_end = lo_indices(block, nhalo)
    hi_domn_start, hi_domn_end, _, _ = hi_indices(block, nhalo)

    ilo_domn_start = [v for (i, v) in enumerate(lo_domn_start) if i in halodims] |> first
    ilo_domn_end = [v for (i, v) in enumerate(lo_domn_end) if i in halodims] |> first
    ihi_domn_start = [v for (i, v) in enumerate(hi_domn_start) if i in halodims] |> first
    ihi_domn_end = [v for (i, v) in enumerate(hi_domn_end) if i in halodims] |> first

    # key -> neighbor id, value -> array index ranges
    return Dict(
        :ilo => (ilo_domn_start:ilo_domn_end,),
        :ihi => (ihi_domn_start:ihi_domn_end,),
    )
end

"""
    domain_donor_ranges_2d(block, nhalo) -> Dict

Get the ranges for each donor region in the doman that sends data to neighbor block halo regions
"""
function domain_donor_ranges_2d(block, nhalo, halodims)

    # domn == domain area
    _, _, lo_domn_start, lo_domn_end = lo_indices(block, nhalo)
    hi_domn_start, hi_domn_end, _, _ = hi_indices(block, nhalo)

    ilo_domn_start, jlo_domn_start = [v for (i, v) in enumerate(lo_domn_start) if i in halodims]
    ilo_domn_end, jlo_domn_end = [v for (i, v) in enumerate(lo_domn_end) if i in halodims]
    ihi_domn_start, jhi_domn_start = [v for (i, v) in enumerate(hi_domn_start) if i in halodims]
    ihi_domn_end, jhi_domn_end = [v for (i, v) in enumerate(hi_domn_end) if i in halodims]

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
    domain_donor_ranges_3d(block, nhalo) -> Dict

Get the ranges for each donor region in the doman that sends data to neighbor block halo regions
"""
function domain_donor_ranges_3d(block, nhalo, halodims)

    # dom == domain area
    _, _, lo_dom_start, lo_dom_end = lo_indices(block, nhalo)
    hi_dom_start, hi_dom_end, _, _ = hi_indices(block, nhalo)

    ilo_dom_start, jlo_dom_start, klo_dom_start = [v for (i, v) in enumerate(lo_dom_start) if i in halodims]
    ilo_dom_end, jlo_dom_end, klo_dom_end = [v for (i, v) in enumerate(lo_dom_end) if i in halodims]
    ihi_dom_start, jhi_dom_start, khi_dom_start = [v for (i, v) in enumerate(hi_dom_start) if i in halodims]
    ihi_dom_end, jhi_dom_end, khi_dom_end = [v for (i, v) in enumerate(hi_dom_end) if i in halodims]

    # key -> neighbor id, value -> array index ranges
    return Dict(
        :ilo => (ilo_dom_start:ilo_dom_end, jlo_dom_start:jhi_dom_end, klo_dom_start:khi_dom_end),
        :ihi => (ihi_dom_start:ihi_dom_end, jlo_dom_start:jhi_dom_end, klo_dom_start:khi_dom_end),
        :jlo => (ilo_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end, klo_dom_start:khi_dom_end),
        :jhi => (ilo_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end, klo_dom_start:khi_dom_end),
        :klo => (ilo_dom_start:ihi_dom_end, jlo_dom_start:jhi_dom_end, klo_dom_start:klo_dom_end),
        :khi => (ilo_dom_start:ihi_dom_end, jlo_dom_start:jhi_dom_end, khi_dom_start:khi_dom_end),
        :ilojlo => (ilo_dom_start:ilo_dom_end, jlo_dom_start:jlo_dom_end, klo_dom_start:khi_dom_end),
        :ihijlo => (ihi_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end, klo_dom_start:khi_dom_end),
        :ilojhi => (ilo_dom_start:ilo_dom_end, jhi_dom_start:jhi_dom_end, klo_dom_start:khi_dom_end),
        :ihijhi => (ihi_dom_start:ihi_dom_end, jhi_dom_start:jhi_dom_end, klo_dom_start:khi_dom_end),
        :ilokhi => (ilo_dom_start:ilo_dom_end, jlo_dom_start:jhi_dom_end, khi_dom_start:khi_dom_end),
        :iloklo => (ilo_dom_start:ilo_dom_end, jlo_dom_start:jhi_dom_end, klo_dom_start:klo_dom_end),
        :ihikhi => (ihi_dom_start:ihi_dom_end, jlo_dom_start:jhi_dom_end, khi_dom_start:khi_dom_end),
        :ihiklo => (ihi_dom_start:ihi_dom_end, jlo_dom_start:jhi_dom_end, klo_dom_start:klo_dom_end),
        :jhikhi => (ilo_dom_start:ilo_dom_end, jhi_dom_start:jhi_dom_end, khi_dom_start:khi_dom_end),
        :jlokhi => (ilo_dom_start:ilo_dom_end, jlo_dom_start:jlo_dom_end, khi_dom_start:khi_dom_end),
        :jhiklo => (ilo_dom_start:ilo_dom_end, jhi_dom_start:jhi_dom_end, klo_dom_start:klo_dom_end),
        :jloklo => (ilo_dom_start:ilo_dom_end, jlo_dom_start:jlo_dom_end, klo_dom_start:klo_dom_end),
        :ilojhikhi => (ilo_dom_start:ilo_dom_end, jhi_dom_start:jhi_dom_end, khi_dom_start:khi_dom_end),
        :ihijhikhi => (ihi_dom_start:ihi_dom_end, jhi_dom_start:jhi_dom_end, khi_dom_start:khi_dom_end),
        :ilojlokhi => (ilo_dom_start:ilo_dom_end, jlo_dom_start:jlo_dom_end, khi_dom_start:khi_dom_end),
        :ihijlokhi => (ihi_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end, khi_dom_start:khi_dom_end),
        :ilojhiklo => (ilo_dom_start:ilo_dom_end, jhi_dom_start:jhi_dom_end, klo_dom_start:klo_dom_end),
        :ihijhiklo => (ihi_dom_start:ihi_dom_end, jhi_dom_start:jhi_dom_end, klo_dom_start:klo_dom_end),
        :ilojloklo => (ilo_dom_start:ilo_dom_end, jlo_dom_start:jlo_dom_end, klo_dom_start:klo_dom_end),
        :ihijloklo => (ihi_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end, klo_dom_start:klo_dom_end),
    )
end

"""
    halo_reciever_ranges_1d(block, nhalo) -> Dict

Get the ranges for each halo region that recieves data from neighbor blocks
"""
function halo_reciever_ranges_1d(block, nhalo, halodims)

    # domn == domain area
    lo_halo_start, lo_halo_end, _, _ = lo_indices(block, nhalo)
    _, _, hi_halo_start, hi_halo_end = hi_indices(block, nhalo)

    ilo_halo_start = [v for (i, v) in enumerate(lo_halo_start) if i in halodims] |> first
    ilo_halo_end = [v for (i, v) in enumerate(lo_halo_end) if i in halodims] |> first
    ihi_halo_start = [v for (i, v) in enumerate(hi_halo_start) if i in halodims] |> first
    ihi_halo_end = [v for (i, v) in enumerate(hi_halo_end) if i in halodims] |> first

    # key -> neighbor id, value -> array index ranges
    return Dict(
        :ilo => (ilo_halo_start:ilo_halo_end,),
        :ihi => (ihi_halo_start:ihi_halo_end,),
    )
end

"""
    halo_reciever_ranges_2d(block, nhalo) -> Dict

Get the ranges for each halo region that recieves data from neighbor blocks
"""
function halo_reciever_ranges_2d(block, nhalo, halodims)

    # domn == domain area
    lo_halo_start, lo_halo_end, lo_domn_start, _ = lo_indices(block, nhalo)
    _, hi_domn_end, hi_halo_start, hi_halo_end = hi_indices(block, nhalo)

    ilo_halo_start, jlo_halo_start = [v for (i, v) in enumerate(lo_halo_start) if i in halodims]
    ilo_halo_end, jlo_halo_end = [v for (i, v) in enumerate(lo_halo_end) if i in halodims]
    ilo_domn_start, jlo_domn_start = [v for (i, v) in enumerate(lo_domn_start) if i in halodims]
    ihi_domn_end, jhi_domn_end = [v for (i, v) in enumerate(hi_domn_end) if i in halodims]
    ihi_halo_start, jhi_halo_start = [v for (i, v) in enumerate(hi_halo_start) if i in halodims]
    ihi_halo_end, jhi_halo_end = [v for (i, v) in enumerate(hi_halo_end) if i in halodims]

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

"""
    halo_reciever_ranges_3d(block, nhalo) -> Dict

Get the ranges for each halo region that recieves data from neighbor blocks
"""
function halo_reciever_ranges_3d(block, nhalo, halodims)

    # domn == domain area
    lo_halo_start, lo_halo_end, lo_domn_start, _ = lo_indices(block, nhalo)
    _, hi_domn_end, hi_halo_start, hi_halo_end = hi_indices(block, nhalo)

    ilo_halo_start, jlo_halo_start, klo_halo_start  = [v for (i, v) in enumerate(lo_halo_start) if i in halodims]
    ilo_halo_end, jlo_halo_end, klo_halo_end  = [v for (i, v) in enumerate(lo_halo_end) if i in halodims]
    ilo_domn_start, jlo_domn_start, klo_domn_start  = [v for (i, v) in enumerate(lo_domn_start) if i in halodims]
    ihi_domn_end, jhi_domn_end, khi_domn_end  = [v for (i, v) in enumerate(hi_domn_end) if i in halodims]
    ihi_halo_start, jhi_halo_start, khi_halo_start  = [v for (i, v) in enumerate(hi_halo_start) if i in halodims]
    ihi_halo_end, jhi_halo_end, khi_halo_end  = [v for (i, v) in enumerate(hi_halo_end) if i in halodims]

    # key -> neighbor id, value -> array index ranges
    return Dict(
        :ilo =>       (ilo_halo_start:ilo_halo_end,jlo_domn_start:jhi_domn_end,klo_domn_start:khi_domn_end),
        :ihi =>       (ihi_halo_start:ihi_halo_end,jlo_domn_start:jhi_domn_end,klo_domn_start:khi_domn_end),
        :jlo =>       (ilo_domn_start:ihi_domn_end,jlo_halo_start:jlo_halo_end,klo_domn_start:khi_domn_end),
        :jhi =>       (ilo_domn_start:ihi_domn_end,jhi_halo_start:jhi_halo_end,klo_domn_start:khi_domn_end),
        :klo =>       (ilo_domn_start:ihi_domn_end,jlo_domn_start:jhi_domn_end,klo_halo_start:klo_halo_end),
        :khi =>       (ilo_domn_start:ihi_domn_end,jlo_domn_start:jhi_domn_end,khi_halo_start:khi_halo_end),
        :ihijhi =>    (ihi_halo_start:ihi_halo_end,jhi_halo_start:jhi_halo_end,klo_domn_start:khi_domn_end),
        :ihijlo =>    (ihi_halo_start:ihi_halo_end,jlo_halo_start:jlo_halo_end,klo_domn_start:khi_domn_end),
        :ilojhi =>    (ilo_halo_start:ilo_halo_end,jhi_halo_start:jhi_halo_end,klo_domn_start:khi_domn_end),
        :ilojlo =>    (ilo_halo_start:ilo_halo_end,jlo_halo_start:jlo_halo_end,klo_domn_start:khi_domn_end),
        :ihikhi =>    (ihi_halo_start:ihi_halo_end,jlo_domn_start:jhi_domn_end,khi_halo_start:khi_halo_end),
        :ihiklo =>    (ihi_halo_start:ihi_halo_end,jlo_domn_start:jhi_domn_end,klo_halo_start:klo_halo_end),
        :iloklo =>    (ilo_halo_start:ilo_halo_end,jlo_domn_start:jhi_domn_end,klo_halo_start:klo_halo_end),
        :ilokhi =>    (ilo_halo_start:ilo_halo_end,jlo_domn_start:jhi_domn_end,khi_halo_start:khi_halo_end),
        :jhiklo =>    (ilo_domn_start:ihi_domn_end,jhi_halo_start:jhi_halo_end,klo_halo_start:klo_halo_end),
        :jloklo =>    (ilo_domn_start:ihi_domn_end,jlo_halo_start:jlo_halo_end,klo_halo_start:klo_halo_end),
        :jhikhi =>    (ilo_domn_start:ihi_domn_end,jhi_halo_start:jhi_halo_end,khi_halo_start:khi_halo_end),
        :jlokhi =>    (ilo_domn_start:ihi_domn_end,jlo_halo_start:jlo_halo_end,khi_halo_start:khi_halo_end),
        :ihijloklo => (ihi_halo_start:ihi_halo_end,jlo_halo_start:jlo_halo_end,klo_halo_start:klo_halo_end),
        :ilojhiklo => (ilo_halo_start:ilo_halo_end,jhi_halo_start:jhi_halo_end,klo_halo_start:klo_halo_end),
        :ilojloklo => (ilo_halo_start:ilo_halo_end,jlo_halo_start:jlo_halo_end,klo_halo_start:klo_halo_end),
        :ihijhiklo => (ihi_halo_start:ihi_halo_end,jhi_halo_start:jhi_halo_end,klo_halo_start:klo_halo_end),
        :ihijlokhi => (ihi_halo_start:ihi_halo_end,jlo_halo_start:jlo_halo_end,khi_halo_start:khi_halo_end),
        :ilojhikhi => (ilo_halo_start:ilo_halo_end,jhi_halo_start:jhi_halo_end,khi_halo_start:khi_halo_end),
        :ilojlokhi => (ilo_halo_start:ilo_halo_end,jlo_halo_start:jlo_halo_end,khi_halo_start:khi_halo_end),
        :ihijhikhi => (ihi_halo_start:ihi_halo_end,jhi_halo_start:jhi_halo_end,khi_halo_start:khi_halo_end),
    )
end

"""

Synchronize the halo regions within each block. This spawns tasks so that
each thread/block copies from it's neighbor block's domain into the current block's halo region.

# Arguments
A::BlockHaloArray, include_periodic_bc=false
"""
function sync_halo!(A::BlockHaloArray, include_periodic_bc=false)
    @sync for tid in 1:length(A.blocks)
        ThreadPools.@tspawnat tid _neighbor_exhange(A, tid, include_periodic_bc)
    end
end

@inline function valid_neighbor(id::Integer)
    if id > 0
        return true
    else
        return false
    end
end

"""Copy data from the neighbor block into the current block's halo region"""
function _neighbor_exhange(A::BlockHaloArray, block_id::Integer, include_periodic_bc=false)
    current_block = @views A.blocks[block_id]

    exchange_map, domain_ranges, halo_ranges = get_neighbor_mapping(A, block_id)

    for (dom_id, halo_id) in exchange_map

        # The convention is that periodic neighbors block ids are < 0 as a hint to the 
        # user and code. 
        if include_periodic_bc
            neighbor_id = abs(A.neighbor_blocks[block_id][halo_id])
        else
            neighbor_id = A.neighbor_blocks[block_id][halo_id]
        end

        if valid_neighbor(neighbor_id)
            donor = @views A.blocks[neighbor_id][.., domain_ranges[dom_id]...]
            halo = @views current_block[.., halo_ranges[halo_id]...]

            copy!(halo, donor) # update the halo region
        end
    end
    return nothing
end

"""Get the reciever and donor indices for 1D halo exchanges"""
function get_neighbor_mapping(A::BlockHaloArray{T,N,1,NBL,AA}, block_id::Integer) where {T,N,NBL,AA}
    donor_ranges = domain_donor_ranges_1d(A.blocks[block_id], A.nhalo, A.halodims)
    halo_ranges = halo_reciever_ranges_1d(A.blocks[block_id], A.nhalo, A.halodims)
    exchange_map = halo_exhange_map_1d()

    return exchange_map, donor_ranges, halo_ranges
end

"""Get the reciever and donor indices for 2D halo exchanges"""
function get_neighbor_mapping(A::BlockHaloArray{T,N,2,NBL,AA}, block_id::Integer) where {T,N,NBL,AA}
    donor_ranges = domain_donor_ranges_2d(A.blocks[block_id], A.nhalo, A.halodims)
    halo_ranges = halo_reciever_ranges_2d(A.blocks[block_id], A.nhalo, A.halodims)
    exchange_map = halo_exhange_map_2d()

    return exchange_map, donor_ranges, halo_ranges
end

"""Get the reciever and donor indices for 3D halo exchanges"""
function get_neighbor_mapping(A::BlockHaloArray{T,N,3,NBL,AA}, block_id::Integer) where {T,N,NBL,AA}
    donor_ranges = domain_donor_ranges_3d(A.blocks[block_id], A.nhalo, A.halodims)
    halo_ranges = halo_reciever_ranges_3d(A.blocks[block_id], A.nhalo, A.halodims)
    exchange_map = halo_exhange_map_3d()

    return exchange_map, donor_ranges, halo_ranges
end

"""Get the upper indices for an array `A` given a number of halo entries `nhalo`"""
function hi_indices(A::AbstractArray, nhalo::Integer)
    hi_halo_end = last.(axes(A))               # end index of the halo region
    hi_halo_start = hi_halo_end .- nhalo .+ 1  # start index of the halo region
    hi_domn_end = hi_halo_start .- 1           # end index of the inner domain
    hi_domn_start = hi_domn_end .- nhalo .+ 1  # start index of the inner domain

    return (hi_domn_start, hi_domn_end, hi_halo_start, hi_halo_end)
end

"""Get the upper indices for an array `A` given a number of halo entries `nhalo`. This 
version properly accounts for non-halo dimensions"""
function hi_indices(A::AbstractArray, nhalo::Integer, halodims::NTuple{N,Integer}) where {N}

    indices = collect.(hi_indices(A, nhalo))

    for index in indices
        for dim in eachindex(index)
            if !(dim in halodims)
                index[dim] = last(axes(A, dim))
            end
        end
    end
    return Tuple.(indices)
end

"""Get the lower indices for an array `A` given a number of halo entries `nhalo`"""
function lo_indices(A::AbstractArray, nhalo::Integer)
    lo_halo_start = first.(axes(A))            # start index of the halo region
    lo_halo_end = lo_halo_start .+ nhalo .- 1  # end index of the halo region
    lo_domn_start = lo_halo_end .+ 1           # start index of the inner domain
    lo_domn_end = lo_domn_start .+ nhalo .- 1  # end index of the inner domain

    return (lo_halo_start, lo_halo_end, lo_domn_start, lo_domn_end)
end

"""Get the lower indices for an array `A` given a number of halo entries `nhalo`. This 
version properly accounts for non-halo dimensions"""
function lo_indices(A::AbstractArray, nhalo::Integer, halodims::NTuple{N,Integer}) where {N}

    indices = collect.(lo_indices(A, nhalo))

    for index in indices
        for dim in eachindex(index)
            if !(dim in halodims)
                index[dim] = first(axes(A, dim))
            end
        end
    end
    return Tuple.(indices)
end

"""Get the neighbor block id's for a 1D decomposition"""
function get_neighbor_blocks_no_periodic(tile_dims::NTuple{1,Integer})
    nblocks = prod(tile_dims)
    nneighbors = 2

    CI = CartesianIndices(tile_dims)
    LI = LinearIndices(tile_dims)
    block_neighbors = Vector{Dict{Symbol,Int}}(undef, nblocks)

    neighbor_block_sym = [:ilo, :ihi]
    for i in LI

        # set up the block neighbors based on cartesian indexing
        neighbor_indices = [(i - 1), (i + 1)]

        block_neighbor_set = Vector{Tuple{Symbol,Int}}(undef, nneighbors)

        # if the index is out of bounds, it's invalid. 
        # Otherwise, save the 1D index of the block 
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
        block_neighbors[i] = Dict(block_neighbor_set)
    end

    return block_neighbors
end

"""Get the neighbor block id's for a 2D decomposition"""
function get_neighbor_blocks_no_periodic(tile_dims::NTuple{2,Integer})
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
        block_neighbors[LI[i, j]] = Dict(block_neighbor_set) # e.g. (:ilo => 4, :jhi => ...), where 4 is the ilo neighbor
    end

    return block_neighbors
end

"""Get the neighbor block id's for a 3D decomposition"""
function get_neighbor_blocks_no_periodic(tile_dims::NTuple{3,Integer})
    nblocks = prod(tile_dims)
    nneighbors = 3^length(tile_dims) - 1

    CI = CartesianIndices(tile_dims)
    LI = LinearIndices(tile_dims)
    block_neighbors = Vector{Dict{Symbol,Int}}(undef, nblocks)

    neighbor_block_sym = [
        :ihi,
        :ihijhi,
        :ihijhikhi,
        :ihijhiklo,
        :ihijlo,
        :ihijlokhi,
        :ihijloklo,
        :ihikhi,
        :ihiklo,
        :ilo,
        :ilojhi,
        :ilojhikhi,
        :ilojhiklo,
        :ilojlo,
        :ilojlokhi,
        :ilojloklo,
        :ilokhi,
        :iloklo,
        :jhi,
        :jhikhi,
        :jhiklo,
        :jlo,
        :jlokhi,
        :jloklo,
        :khi,
        :klo,
    ]

    for block_idx in CI
        i, j, k = Tuple(block_idx)

        # set up the block neighbors based on cartesian indexing
        neighbor_indices = [
            (i + 1, j, k),
            (i + 1, j + 1, k),
            (i + 1, j + 1, k + 1),
            (i + 1, j + 1, k - 1),
            (i + 1, j - 1, k),
            (i + 1, j - 1, k + 1),
            (i + 1, j - 1, k - 1),
            (i + 1, j, k + 1),
            (i + 1, j, k - 1),
            (i - 1, j, k),
            (i - 1, j + 1, k),
            (i - 1, j + 1, k + 1),
            (i - 1, j + 1, k - 1),
            (i - 1, j - 1, k),
            (i - 1, j - 1, k + 1),
            (i - 1, j - 1, k - 1),
            (i - 1, j, k + 1),
            (i - 1, j, k - 1),
            (i, j + 1, k),
            (i, j + 1, k + 1),
            (i, j + 1, k - 1),
            (i, j - 1, k),
            (i, j - 1, k + 1),
            (i, j - 1, k - 1),
            (i, j, k + 1),
            (i, j, k - 1),
        ]

        block_neighbor_set = Vector{Tuple{Symbol,Int}}(undef, nneighbors)

        # if the index is out of bounds, it's invalid. 
        # Otherwise, save the 1D index of the block 
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
        block_neighbors[LI[i, j, k]] = Dict(block_neighbor_set)
    end

    return block_neighbors
end

"""Get the neighbor block id's for a 1D decomposition (with periodic edges)"""
function get_neighbor_blocks(tile_dims::NTuple{1,Integer})
    nblocks = prod(tile_dims)
    nneighbors = 2

    CI = CartesianIndices(tile_dims)
    LI = LinearIndices(tile_dims)
    block_neighbors = Vector{Dict{Symbol,Int}}(undef, nblocks)

    neighbor_block_sym = [:ilo, :ihi]
    for i in LI

        # set up the block neighbors based on cartesian indexing
        neighbor_indices = [[i - 1], [i + 1]]

        block_neighbor_set = Vector{Tuple{Symbol,Int}}(undef, nneighbors)

        # if the index is out of bounds, it's invalid. 
        # Otherwise, save the 1D index of the block 
        # that is the proper neighbor
        for (neighbor_id, idx) in enumerate(neighbor_indices)
            periodic = false
            neighbor_symbol = neighbor_block_sym[neighbor_id]
            for (i, dim) in enumerate(idx)
                if dim < 1
                    periodic = true
                    idx[i] = tile_dims[i]
                elseif dim > tile_dims[i]
                    periodic = true
                    idx[i] = 1
                end
            end
            neighbor_block = LI[CI[idx...]]
            if periodic
                neighbor_block *= -1
            end

            block_neighbor_set[neighbor_id] = (neighbor_symbol, neighbor_block)
        end
        block_neighbors[i] = Dict(block_neighbor_set)
    end

    return block_neighbors
end

"""Get the neighbor block id's for a 2D decomposition (with periodic edges)"""
function get_neighbor_blocks(tile_dims::NTuple{2,Integer})
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
            [i - 1, j - 1],
            [i, j - 1],
            [i + 1, j - 1],
            [i + 1, j],
            [i + 1, j + 1],
            [i, j + 1],
            [i - 1, j + 1],
            [i - 1, j]
        ]

        block_neighbor_set = Vector{Tuple{Symbol,Int}}(undef, nneighbors)
        for (neighbor_id, idx) in enumerate(neighbor_indices)
            periodic = false
            neighbor_symbol = neighbor_block_sym[neighbor_id]
            for (i, dim) in enumerate(idx)
                if dim < 1
                    periodic = true
                    idx[i] = tile_dims[i]
                elseif dim > tile_dims[i]
                    periodic = true
                    idx[i] = 1
                end
            end
            neighbor_block = LI[CI[idx...]]
            if periodic
                neighbor_block *= -1
            end

            block_neighbor_set[neighbor_id] = (neighbor_symbol, neighbor_block)
        end
        block_neighbors[LI[i, j]] = Dict(block_neighbor_set) # e.g. (:ilo => 4, :jhi => ...), where 4 is the ilo neighbor
    end

    return block_neighbors
end

"""Get the neighbor block id's for a 3D decomposition (with periodic edges)"""
function get_neighbor_blocks(tile_dims::NTuple{3,Integer})
    nblocks = prod(tile_dims)
    nneighbors = 3^length(tile_dims) - 1

    CI = CartesianIndices(tile_dims)
    LI = LinearIndices(tile_dims)
    block_neighbors = Vector{Dict{Symbol,Int}}(undef, nblocks)

    neighbor_block_sym = [
        :ihi,
        :ihijhi,
        :ihijhikhi,
        :ihijhiklo,
        :ihijlo,
        :ihijlokhi,
        :ihijloklo,
        :ihikhi,
        :ihiklo,
        :ilo,
        :ilojhi,
        :ilojhikhi,
        :ilojhiklo,
        :ilojlo,
        :ilojlokhi,
        :ilojloklo,
        :ilokhi,
        :iloklo,
        :jhi,
        :jhikhi,
        :jhiklo,
        :jlo,
        :jlokhi,
        :jloklo,
        :khi,
        :klo,
    ]

    for block_idx in CI
        i, j, k = Tuple(block_idx)

        # set up the block neighbors based on cartesian indexing
        neighbor_indices = [
            [i + 1, j, k],
            [i + 1, j + 1, k],
            [i + 1, j + 1, k + 1],
            [i + 1, j + 1, k - 1],
            [i + 1, j - 1, k],
            [i + 1, j - 1, k + 1],
            [i + 1, j - 1, k - 1],
            [i + 1, j, k + 1],
            [i + 1, j, k - 1],
            [i - 1, j, k],
            [i - 1, j + 1, k],
            [i - 1, j + 1, k + 1],
            [i - 1, j + 1, k - 1],
            [i - 1, j - 1, k],
            [i - 1, j - 1, k + 1],
            [i - 1, j - 1, k - 1],
            [i - 1, j, k + 1],
            [i - 1, j, k - 1],
            [i, j + 1, k],
            [i, j + 1, k + 1],
            [i, j + 1, k - 1],
            [i, j - 1, k],
            [i, j - 1, k + 1],
            [i, j - 1, k - 1],
            [i, j, k + 1],
            [i, j, k - 1],
        ]

        block_neighbor_set = Vector{Tuple{Symbol,Int}}(undef, nneighbors)

        # if the index is out of bounds, it's invalid. 
        # Otherwise, save the 1D index of the block 
        # that is the proper neighbor
        for (neighbor_id, idx) in enumerate(neighbor_indices)
            periodic = false
            neighbor_symbol = neighbor_block_sym[neighbor_id]
            for (i, dim) in enumerate(idx)
                if dim < 1
                    periodic = true
                    idx[i] = tile_dims[i]
                elseif dim > tile_dims[i]
                    periodic = true
                    idx[i] = 1
                end
            end
            neighbor_block = LI[CI[idx...]]
            if periodic
                neighbor_block *= -1
            end

            block_neighbor_set[neighbor_id] = (neighbor_symbol, neighbor_block)
        end
        block_neighbors[LI[i, j, k]] = Dict(block_neighbor_set)
    end

    return block_neighbors
end

