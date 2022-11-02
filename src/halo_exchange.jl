const halo_1d_exchange_map = Dict(:ilo => :ihi, :ihi => :ilo)

const halo_2d_exchange_map = Dict(
    :ilo => :ihi, :ihi => :ilo, :jlo => :jhi, :jhi => :jlo,
    :ilojlo => :ihijhi, :ihijlo => :ilojhi, :ihijhi => :ilojlo, :ilojhi => :ihijlo
)

const halo_3d_exchange_map = Dict(
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

const hmap = (halo_1d_exchange_map, halo_2d_exchange_map, halo_3d_exchange_map)

"""1D halo exchange mapping, e.g. donor => reciever block ID"""
function halo_exchange_map(::NTuple{1,Int})
    return Dict(:ilo => :ihi, :ihi => :ilo)
end

"""2D halo exchange mapping, e.g. donor => reciever block ID"""
function halo_exchange_map(::NTuple{2,Int})
    return Dict(:ilo => :ihi, :ihi => :ilo, :jlo => :jhi, :jhi => :jlo,
        :ilojlo => :ihijhi, :ihijlo => :ilojhi, :ihijhi => :ilojlo, :ilojhi => :ihijlo)
end

"""3D halo exchange mapping, e.g. donor => reciever block ID"""
function halo_exchange_map(::NTuple{3,Int})
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
    domain_donor_ranges(block, nhal, halodims::NTuple{1, Int}) -> Dict

Get the 1D ranges for each donor region in the doman that sends data to neighbor block halo regions
"""
function domain_donor_ranges(block, nhalo, halodims::NTuple{1,Int})

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
    domain_donor_ranges(block, nhal, halodims::NTuple{2, Int}) -> Dict

Get the 2D ranges for each donor region in the doman that sends data to neighbor block halo regions
"""
function domain_donor_ranges(block, nhalo, halodims::NTuple{2,Int})

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
    domain_donor_ranges(block, nhal, halodims::NTuple{3, Int}) -> Dict

Get the 3D ranges for each donor region in the doman that sends data to neighbor block halo regions
"""
function domain_donor_ranges(block, nhalo, halodims::NTuple{3,Int})

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
        :jhikhi => (ilo_dom_start:ihi_dom_end, jhi_dom_start:jhi_dom_end, khi_dom_start:khi_dom_end),
        :jlokhi => (ilo_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end, khi_dom_start:khi_dom_end),
        :jhiklo => (ilo_dom_start:ihi_dom_end, jhi_dom_start:jhi_dom_end, klo_dom_start:klo_dom_end),
        :jloklo => (ilo_dom_start:ihi_dom_end, jlo_dom_start:jlo_dom_end, klo_dom_start:klo_dom_end),
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
    halo_reciever_ranges(block, nhalo, halodims::NTuple{1, Int}) -> Dict

Get the 1D ranges for each halo region that recieves data from neighbor blocks
"""
function halo_reciever_ranges(block, nhalo, halodims::NTuple{1,Int})

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
    halo_reciever_ranges(block, nhalo, halodims::NTuple{2, Int}) -> Dict

Get the 2D ranges for each halo region that recieves data from neighbor blocks
"""
function halo_reciever_ranges(block, nhalo, halodims::NTuple{2,Int})

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
    halo_reciever_ranges(block, nhalo, halodims::NTuple{3, Int}) -> Dict

Get the 3D ranges for each halo region that recieves data from neighbor blocks
"""
function halo_reciever_ranges(block, nhalo, halodims::NTuple{3,Int})

    # domn == domain area
    lo_halo_start, lo_halo_end, lo_domn_start, _ = lo_indices(block, nhalo)
    _, hi_domn_end, hi_halo_start, hi_halo_end = hi_indices(block, nhalo)

    ilo_halo_start, jlo_halo_start, klo_halo_start = [v for (i, v) in enumerate(lo_halo_start) if i in halodims]
    ilo_halo_end, jlo_halo_end, klo_halo_end = [v for (i, v) in enumerate(lo_halo_end) if i in halodims]
    ilo_domn_start, jlo_domn_start, klo_domn_start = [v for (i, v) in enumerate(lo_domn_start) if i in halodims]
    ihi_domn_end, jhi_domn_end, khi_domn_end = [v for (i, v) in enumerate(hi_domn_end) if i in halodims]
    ihi_halo_start, jhi_halo_start, khi_halo_start = [v for (i, v) in enumerate(hi_halo_start) if i in halodims]
    ihi_halo_end, jhi_halo_end, khi_halo_end = [v for (i, v) in enumerate(hi_halo_end) if i in halodims]

    # key -> neighbor id, value -> array index ranges
    return Dict(
        :ilo => (ilo_halo_start:ilo_halo_end, jlo_domn_start:jhi_domn_end, klo_domn_start:khi_domn_end),
        :ihi => (ihi_halo_start:ihi_halo_end, jlo_domn_start:jhi_domn_end, klo_domn_start:khi_domn_end),
        :jlo => (ilo_domn_start:ihi_domn_end, jlo_halo_start:jlo_halo_end, klo_domn_start:khi_domn_end),
        :jhi => (ilo_domn_start:ihi_domn_end, jhi_halo_start:jhi_halo_end, klo_domn_start:khi_domn_end),
        :klo => (ilo_domn_start:ihi_domn_end, jlo_domn_start:jhi_domn_end, klo_halo_start:klo_halo_end),
        :khi => (ilo_domn_start:ihi_domn_end, jlo_domn_start:jhi_domn_end, khi_halo_start:khi_halo_end),
        :ihijhi => (ihi_halo_start:ihi_halo_end, jhi_halo_start:jhi_halo_end, klo_domn_start:khi_domn_end),
        :ihijlo => (ihi_halo_start:ihi_halo_end, jlo_halo_start:jlo_halo_end, klo_domn_start:khi_domn_end),
        :ilojhi => (ilo_halo_start:ilo_halo_end, jhi_halo_start:jhi_halo_end, klo_domn_start:khi_domn_end),
        :ilojlo => (ilo_halo_start:ilo_halo_end, jlo_halo_start:jlo_halo_end, klo_domn_start:khi_domn_end),
        :ihikhi => (ihi_halo_start:ihi_halo_end, jlo_domn_start:jhi_domn_end, khi_halo_start:khi_halo_end),
        :ihiklo => (ihi_halo_start:ihi_halo_end, jlo_domn_start:jhi_domn_end, klo_halo_start:klo_halo_end),
        :iloklo => (ilo_halo_start:ilo_halo_end, jlo_domn_start:jhi_domn_end, klo_halo_start:klo_halo_end),
        :ilokhi => (ilo_halo_start:ilo_halo_end, jlo_domn_start:jhi_domn_end, khi_halo_start:khi_halo_end),
        :jloklo => (ilo_domn_start:ihi_domn_end, jlo_halo_start:jlo_halo_end, klo_halo_start:klo_halo_end),
        :jhikhi => (ilo_domn_start:ihi_domn_end, jhi_halo_start:jhi_halo_end, khi_halo_start:khi_halo_end),
        :jhiklo => (ilo_domn_start:ihi_domn_end, jhi_halo_start:jhi_halo_end, klo_halo_start:klo_halo_end),
        :jlokhi => (ilo_domn_start:ihi_domn_end, jlo_halo_start:jlo_halo_end, khi_halo_start:khi_halo_end),
        :ihijloklo => (ihi_halo_start:ihi_halo_end, jlo_halo_start:jlo_halo_end, klo_halo_start:klo_halo_end),
        :ilojhiklo => (ilo_halo_start:ilo_halo_end, jhi_halo_start:jhi_halo_end, klo_halo_start:klo_halo_end),
        :ilojloklo => (ilo_halo_start:ilo_halo_end, jlo_halo_start:jlo_halo_end, klo_halo_start:klo_halo_end),
        :ihijhiklo => (ihi_halo_start:ihi_halo_end, jhi_halo_start:jhi_halo_end, klo_halo_start:klo_halo_end),
        :ihijlokhi => (ihi_halo_start:ihi_halo_end, jlo_halo_start:jlo_halo_end, khi_halo_start:khi_halo_end),
        :ilojhikhi => (ilo_halo_start:ilo_halo_end, jhi_halo_start:jhi_halo_end, khi_halo_start:khi_halo_end),
        :ilojlokhi => (ilo_halo_start:ilo_halo_end, jlo_halo_start:jlo_halo_end, khi_halo_start:khi_halo_end),
        :ihijhikhi => (ihi_halo_start:ihi_halo_end, jhi_halo_start:jhi_halo_end, khi_halo_start:khi_halo_end),
    )
end

"""
    updatehalo!(A, include_periodic_bc=false)

Synchronize the halo regions within each block. This spawns tasks so that
each thread/block copies from it's neighbor block's domain into the current block's halo region.

# Arguments
 - `A`: `BlockHaloArray`
 - `include_periodic_bc`: Update the halo regions that are on periodic boundaries
"""
function updatehalo!(A::BlockHaloArray, include_periodic_bc)
    @sync for tid in eachindex(A.blocks)
        @tspawnat tid updateblockhalo!(A, tid, include_periodic_bc)
    end
end

function updatehalo!(A::BlockHaloArray)
    @sync for tid in eachindex(A.blocks)
        @tspawnat tid updateblockhalo!(A, tid)
    end
end

@inline function valid_neighbor(id::Integer)
    if id > 0
        return true
    else
        return false
    end
end

"""
    updateblockhalo!(A::BlockHaloArray, block_id::Integer, include_periodic_bc=false)

Copy data from the neighbor block into the current block's halo region

# Arguments
 - `A`: `BlockHaloArray`
 - `block_id::Integer`: Block index
 - `include_periodic_bc`: Update the halo regions that are on periodic boundaries
"""
updateblockhalo!(A, current_block_id) = updateblockhalo!(A, current_block_id, false)

function updateblockhalo!(A::BlockHaloArray{T,N,NH,NBL,AA,SA}, current_block_id::Integer, include_periodic_bc::Bool) where {T,N,NH,NBL,AA,SA}
    exchange_map = hmap[NH]

    for (dom_id, halo_id) in exchange_map

        neighbor_block_id = A.neighbor_blocks[current_block_id][halo_id]
        
        # The convention is that periodic neighbors block ids are < 0 as a hint to the
        # user and code.
        if include_periodic_bc
            neighbor_block_id = abs(neighbor_block_id)
        end

        if valid_neighbor(neighbor_block_id)
            copy!(
                A._halo_views[current_block_id][halo_id],
                A._donor_views[neighbor_block_id][dom_id]
            ) # copy donor -> halo

            # for idx in eachindex(A._halo_views[current_block_id][halo_id])
            #     A._halo_views[current_block_id][halo_id][idx] = A._donor_views[neighbor_block_id][dom_id][idx]
            # end
        end
    end
    return nothing
end

"""Get the upper indices for an array `A` given a number of halo entries `nhalo`"""
function hi_indices(A, nhalo)
    hmod = 1 .* (nhalo .== 0)
    hi_halo_end = last.(axes(A))
    hi_halo_start = hi_halo_end .- nhalo .+ 1 .- hmod
    hi_domn_end = hi_halo_start .- 1 .+ hmod
    hi_domn_start = hi_domn_end .- nhalo .+ 1 .- hmod

    return (hi_domn_start, hi_domn_end, hi_halo_start, hi_halo_end)
end

"""Get the upper indices for an array `A` given a number of halo entries `nhalo`. This
version properly accounts for non-halo dimensions"""
function hi_indices(A::AbstractArray{T,N}, nhalo, halodims) where {N,T}

    function f(i, halodims, nhalo)
        if i in halodims
            return nhalo
        else
            return 0
        end
    end
    halo_tuple = ntuple(i -> f(i, halodims, nhalo), Val(N))
    hi_indices(A, halo_tuple)
end

"""Get the lower indices for an array `A` given a number of halo entries `nhalo`"""
function lo_indices(A, nhalo)
    lmod = 1 .* (nhalo .== 0)
    lo_halo_start = first.(axes(A))
    lo_halo_end = lo_halo_start .+ nhalo .- 1 .+ lmod
    lo_domn_start = lo_halo_end .+ 1 .- lmod
    lo_domn_end = lo_domn_start .+ nhalo .- 1 .+ lmod

    return (lo_halo_start, lo_halo_end, lo_domn_start, lo_domn_end)
end

"""Get the lower indices for an array `A` given a number of halo entries `nhalo`. This
version properly accounts for non-halo dimensions"""
function lo_indices(A::AbstractArray{T,N}, nhalo, halodims) where {N,T}

    function f(i, halodims, nhalo)
        if i in halodims
            return nhalo
        else
            return 0
        end
    end
    halo_tuple = ntuple(i -> f(i, halodims, nhalo), Val(N))
    lo_indices(A, halo_tuple)
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

"""
Generate the SubArray views of each halo region in `A::BlockHaloArray` that
"reciever" views. These are the regions updated in the halo exchange copy/update.
"""
function gen_halo_views(A)

    blk_halo_views = Vector{Dict{Symbol,SubArray}}(undef, nblocks(A))

    for blk_i in eachindex(A.blocks)
        halo_reciever = halo_reciever_ranges(A.blocks[blk_i], A.nhalo, A.halodims)
        blk_halo_views[blk_i] = Dict(
            k => @view A.blocks[blk_i][.., halo_reciever[k]...] for k in keys(halo_reciever)
        )
    end

    blk_halo_views
end

function gen_halo_views(blocks, nhalo, halodims)

    blk_halo_views = Vector{Dict{Symbol,SubArray}}(undef, length(blocks))

    for blk_i in eachindex(blocks)
        halo_reciever = halo_reciever_ranges(blocks[blk_i], nhalo, halodims)
        blk_halo_views[blk_i] = Dict(
            k => @view blocks[blk_i][.., halo_reciever[k]...] for k in keys(halo_reciever)
        )
    end

    blk_halo_views
end

"""
Generate the SubArray views of each domain region in `A::BlockHaloArray` that are called
"donor" views. These are the regions copied from in the halo exchange copy/update.
"""
function gen_donor_views(A)

    blk_donor_views = Vector{Dict{Symbol,SubArray}}(undef, nblocks(A))

    for blk_i in eachindex(A.blocks)
        donor = domain_donor_ranges(A.blocks[blk_i], A.nhalo, A.halodims)

        blk_donor_views[blk_i] = Dict(
            k => @view A.blocks[blk_i][.., donor[k]...] for k in keys(donor)
        )
    end

    blk_donor_views
end

function gen_donor_views(blocks, nhalo, halodims)

    blk_donor_views = Vector{Dict{Symbol,SubArray}}(undef, length(blocks))

    for blk_i in eachindex(blocks)
        donor = domain_donor_ranges(blocks[blk_i], nhalo, halodims)

        blk_donor_views[blk_i] = Dict(
            k => @view blocks[blk_i][.., donor[k]...] for k in keys(donor)
        )
    end

    blk_donor_views
end
