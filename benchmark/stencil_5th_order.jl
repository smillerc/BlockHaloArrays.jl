import Pkg

Pkg.activate("..")

using .Threads
using ThreadPools
using ThreadPinning: pinthreads
using BlockHaloArrays
using LoopVectorization

abstract type AbstractReconstructionScheme end
abstract type AbstractMUSCLReconstruction <: AbstractReconstructionScheme end
abstract type AbstractLimiter end

const global L = 1
const global R = 2
const global SMALL_NUM = 1e-30

struct MLP5Reconstruction <: AbstractMUSCLReconstruction end


function reconstruct_with_blocks!(RS::MLP5Reconstruction, U, i_edge, j_edge, blockid, ϵ=eps(Float64))

    # i_edge & j_edge holds the i+1/2 & j+1/2 L and R values
    # these limits must be the no-halo versions. If the true array size is (10,10) with
    # 2 halo cells, the limits must be ilo, ihi = 3, 8 and jlo, jhi = 3, 8
    ilo, ihi, jlo, jhi = first(U).loop_limits[blockid]

    for (ϕ̄, ϕ_recon) in zip(U, i_edge) # (ρ, u, v, p)
        @turbo for j = jlo:jhi
            for i = ilo-1:ihi
                Δ_minus_half = ϕ̄[blockid][i, j] - ϕ̄[blockid][i-1, j]
                Δ_plus_three_half = ϕ̄[blockid][i+2, j] - ϕ̄[blockid][i+1, j]
                Δ_plus_half = ϕ̄[blockid][i+1, j] - ϕ̄[blockid][i, j]
                Δ_plus_five_half = ϕ̄[blockid][i+3, j] - ϕ̄[blockid][i+2, j]
                Δ_minus_three_half = ϕ̄[blockid][i-1, j] - ϕ̄[blockid][i-2, j]

                Δ_minus_half = Δ_minus_half * (abs(Δ_minus_half) >= ϵ)
                Δ_plus_three_half = Δ_plus_three_half * (abs(Δ_plus_three_half) >= ϵ)
                Δ_plus_half = Δ_plus_half * (abs(Δ_plus_half) >= ϵ)
                Δ_minus_three_half = Δ_minus_three_half * (abs(Δ_minus_three_half) >= ϵ)
                Δ_plus_five_half = Δ_plus_five_half * (abs(Δ_plus_five_half) >= ϵ)

                rL = Δ_plus_half / (Δ_minus_half + SMALL_NUM)
                rL_m1 = Δ_minus_half / (Δ_minus_three_half + SMALL_NUM)
                rL_p1 = Δ_plus_three_half / (Δ_plus_half + SMALL_NUM)

                rR = Δ_minus_half / (Δ_plus_half + SMALL_NUM)
                rR_p1 = Δ_plus_half / (Δ_plus_three_half + SMALL_NUM)
                rR_p2 = Δ_plus_three_half / (Δ_plus_five_half + SMALL_NUM)

                β_L = ((-2 / (rL_m1 + SMALL_NUM)) + 11.0 + 24rL - 3rL * rL_p1) / 30.0
                β_R = ((-2 / (rR_p2 + SMALL_NUM)) + 11.0 + 24rR_p1 - 3rR_p1 * rR) / 30.0

                tanθ_ij = abs(ϕ̄[blockid][i, j+1] - ϕ̄[blockid][i, j-1]) / (abs(ϕ̄[blockid][i+1, j] - ϕ̄[blockid][i-1, j]) + SMALL_NUM)
                tanθ_ij_p1 = abs(ϕ̄[blockid][i+1, j+1] - ϕ̄[blockid][i+1, j-1]) / (abs(ϕ̄[blockid][i+2, j] - ϕ̄[blockid][i, j]) + SMALL_NUM)
                α_L, α_R = get_alpha(rL, rR_p1, tanθ_ij, tanθ_ij_p1)

                ϕ_recon[blockid][1, i, j] = ϕ̄[blockid][i, j] + 0.5max(0, min(α_L * rL, α_L, β_L)) * Δ_minus_half
                ϕ_recon[blockid][2, i, j] = ϕ̄[blockid][i+1, j] - 0.5max(0, min(α_R * rR_p1, α_R, β_R)) * Δ_plus_three_half
            end
        end
    end

    for (ϕ̄, ϕ_recon) in zip(U, j_edge) # (ρ, u, v, p)
        @turbo for j = jlo-1:jhi
            for i = ilo:ihi
                Δ_minus_half = ϕ̄[blockid][i, j] - ϕ̄[blockid][i, j-1]
                Δ_plus_three_half = ϕ̄[blockid][i, j+2] - ϕ̄[blockid][i, j+1]
                Δ_plus_half = ϕ̄[blockid][i, j+1] - ϕ̄[blockid][i, j]
                Δ_minus_three_half = ϕ̄[blockid][i, j-1] - ϕ̄[blockid][i, j-2]
                Δ_plus_five_half = ϕ̄[blockid][i, j+3] - ϕ̄[blockid][i, j+3]

                Δ_minus_half = Δ_minus_half * (abs(Δ_minus_half) >= ϵ)
                Δ_plus_three_half = Δ_plus_three_half * (abs(Δ_plus_three_half) >= ϵ)
                Δ_plus_half = Δ_plus_half * (abs(Δ_plus_half) >= ϵ)
                Δ_minus_three_half = Δ_minus_three_half * (abs(Δ_minus_three_half) >= ϵ)
                Δ_plus_five_half = Δ_plus_five_half * (abs(Δ_plus_five_half) >= ϵ)

                rL = Δ_plus_half / (Δ_minus_half + SMALL_NUM)
                rL_m1 = Δ_minus_half / (Δ_minus_three_half + SMALL_NUM)
                rL_p1 = Δ_plus_three_half / (Δ_plus_half + SMALL_NUM)

                rR = Δ_minus_half / (Δ_plus_half + SMALL_NUM)
                rR_p1 = Δ_plus_half / (Δ_plus_three_half + SMALL_NUM)
                rR_p2 = Δ_plus_three_half / (Δ_plus_five_half + SMALL_NUM)

                β_L = ((-2 / (rL_m1 + SMALL_NUM)) + 11.0 + 24rL - 3rL * rL_p1) / 30.0
                β_R = ((-2 / (rR_p2 + SMALL_NUM)) + 11.0 + 24rR_p1 - 3rR_p1 * rR) / 30.0

                tanθ_ij = abs(ϕ̄[blockid][i+1, j] - ϕ̄[blockid][i-1, j]) / (abs(ϕ̄[blockid][i, j+1] - ϕ̄[blockid][i, j-1]) + SMALL_NUM)
                tanθ_ij_p1 = abs(ϕ̄[blockid][i+1, j+1] - ϕ̄[blockid][i-1, j+1]) / (abs(ϕ̄[blockid][i, j+2] - ϕ̄[blockid][i, j]) + SMALL_NUM)
                α_L, α_R = get_alpha(rL, rR_p1, tanθ_ij, tanθ_ij_p1)

                ϕ_recon[blockid][1, i, j] = ϕ̄[blockid][i, j] + 0.5max(0, min(α_L * rL, α_L, β_L)) * Δ_minus_half
                ϕ_recon[blockid][2, i, j] = ϕ̄[blockid][i, j+1] - 0.5max(0, min(α_R * rR_p1, α_R, β_R)) * Δ_plus_three_half
            end
        end
    end

    return nothing
end

function reconstruct_flat!(RS::MLP5Reconstruction, U, i_edge, j_edge, looplimits, ϵ=eps(Float64))

    # i_edge & j_edge holds the i+1/2 & j+1/2 L and R values
    # these limits must be the no-halo versions. If the true array size is (10,10) with
    # 2 halo cells, the limits must be ilo, ihi = 3, 8 and jlo, jhi = 3, 8
    ilo, ihi, jlo, jhi = looplimits

    for (ϕ̄, ϕ_recon) in zip(U, i_edge) # (ρ, u, v, p)
        @tturbo for j = jlo:jhi
            for i = ilo-1:ihi
                Δ_minus_half = ϕ̄[i, j] - ϕ̄[i-1, j]
                Δ_plus_three_half = ϕ̄[i+2, j] - ϕ̄[i+1, j]
                Δ_plus_half = ϕ̄[i+1, j] - ϕ̄[i, j]
                Δ_plus_five_half = ϕ̄[i+3, j] - ϕ̄[i+2, j]
                Δ_minus_three_half = ϕ̄[i-1, j] - ϕ̄[i-2, j]

                Δ_minus_half = Δ_minus_half * (abs(Δ_minus_half) >= ϵ)
                Δ_plus_three_half = Δ_plus_three_half * (abs(Δ_plus_three_half) >= ϵ)
                Δ_plus_half = Δ_plus_half * (abs(Δ_plus_half) >= ϵ)
                Δ_minus_three_half = Δ_minus_three_half * (abs(Δ_minus_three_half) >= ϵ)
                Δ_plus_five_half = Δ_plus_five_half * (abs(Δ_plus_five_half) >= ϵ)

                rL = Δ_plus_half / (Δ_minus_half + SMALL_NUM)
                rL_m1 = Δ_minus_half / (Δ_minus_three_half + SMALL_NUM)
                rL_p1 = Δ_plus_three_half / (Δ_plus_half + SMALL_NUM)

                rR = Δ_minus_half / (Δ_plus_half + SMALL_NUM)
                rR_p1 = Δ_plus_half / (Δ_plus_three_half + SMALL_NUM)
                rR_p2 = Δ_plus_three_half / (Δ_plus_five_half + SMALL_NUM)

                β_L = ((-2 / (rL_m1 + SMALL_NUM)) + 11.0 + 24rL - 3rL * rL_p1) / 30.0
                β_R = ((-2 / (rR_p2 + SMALL_NUM)) + 11.0 + 24rR_p1 - 3rR_p1 * rR) / 30.0

                tanθ_ij = abs(ϕ̄[i, j+1] - ϕ̄[i, j-1]) / (abs(ϕ̄[i+1, j] - ϕ̄[i-1, j]) + SMALL_NUM)
                tanθ_ij_p1 = abs(ϕ̄[i+1, j+1] - ϕ̄[i+1, j-1]) / (abs(ϕ̄[i+2, j] - ϕ̄[i, j]) + SMALL_NUM)
                α_L, α_R = get_alpha(rL, rR_p1, tanθ_ij, tanθ_ij_p1)

                ϕ_recon[1, i, j] = ϕ̄[i, j] + 0.5max(0, min(α_L * rL, α_L, β_L)) * Δ_minus_half
                ϕ_recon[2, i, j] = ϕ̄[i+1, j] - 0.5max(0, min(α_R * rR_p1, α_R, β_R)) * Δ_plus_three_half
            end
        end
    end

    for (ϕ̄, ϕ_recon) in zip(U, j_edge) # (ρ, u, v, p)
        @tturbo for j = jlo-1:jhi
            for i = ilo:ihi
                Δ_minus_half = ϕ̄[i, j] - ϕ̄[i, j-1]
                Δ_plus_three_half = ϕ̄[i, j+2] - ϕ̄[i, j+1]
                Δ_plus_half = ϕ̄[i, j+1] - ϕ̄[i, j]
                Δ_minus_three_half = ϕ̄[i, j-1] - ϕ̄[i, j-2]
                Δ_plus_five_half = ϕ̄[i, j+3] - ϕ̄[i, j+3]

                Δ_minus_half = Δ_minus_half * (abs(Δ_minus_half) >= ϵ)
                Δ_plus_three_half = Δ_plus_three_half * (abs(Δ_plus_three_half) >= ϵ)
                Δ_plus_half = Δ_plus_half * (abs(Δ_plus_half) >= ϵ)
                Δ_minus_three_half = Δ_minus_three_half * (abs(Δ_minus_three_half) >= ϵ)
                Δ_plus_five_half = Δ_plus_five_half * (abs(Δ_plus_five_half) >= ϵ)

                rL = Δ_plus_half / (Δ_minus_half + SMALL_NUM)
                rL_m1 = Δ_minus_half / (Δ_minus_three_half + SMALL_NUM)
                rL_p1 = Δ_plus_three_half / (Δ_plus_half + SMALL_NUM)

                rR = Δ_minus_half / (Δ_plus_half + SMALL_NUM)
                rR_p1 = Δ_plus_half / (Δ_plus_three_half + SMALL_NUM)
                rR_p2 = Δ_plus_three_half / (Δ_plus_five_half + SMALL_NUM)

                β_L = ((-2 / (rL_m1 + SMALL_NUM)) + 11.0 + 24rL - 3rL * rL_p1) / 30.0
                β_R = ((-2 / (rR_p2 + SMALL_NUM)) + 11.0 + 24rR_p1 - 3rR_p1 * rR) / 30.0

                tanθ_ij = abs(ϕ̄[i+1, j] - ϕ̄[i-1, j]) / (abs(ϕ̄[i, j+1] - ϕ̄[i, j-1]) + SMALL_NUM)
                tanθ_ij_p1 = abs(ϕ̄[i+1, j+1] - ϕ̄[i-1, j+1]) / (abs(ϕ̄[i, j+2] - ϕ̄[i, j]) + SMALL_NUM)
                α_L, α_R = get_alpha(rL, rR_p1, tanθ_ij, tanθ_ij_p1)

                ϕ_recon[1, i, j] = ϕ̄[i, j] + 0.5max(0, min(α_L * rL, α_L, β_L)) * Δ_minus_half
                ϕ_recon[2, i, j] = ϕ̄[i, j+1] - 0.5max(0, min(α_R * rR_p1, α_R, β_R)) * Δ_plus_three_half
            end
        end
    end

    return nothing
end

@inline function get_alpha(rL, rR_p1, tanθ_ij, tanθ_ij_p1)
    α_L_term = ((2max(1, rL) *
                 (1 + max(0, (tanθ_ij_p1 / (rR_p1 + SMALL_NUM))))) / (1 + tanθ_ij))

    α_R_term = ((2max(1, rR_p1) *
                 (1 + max(0, (tanθ_ij / (rL + SMALL_NUM))))) / (1 + tanθ_ij_p1))

    # This is the g(x) = max(1, min(2, alpha)) function
    α_L = max(1, min(2, α_L_term))
    α_R = max(1, min(2, α_R_term))

    return α_L, α_R
end

function run()

    @show nthreads()

    pinthreads(:compact)

    recon = MLP5Reconstruction()

    dims = (512, 1024) .* 2
    @show dims
    nhalo = 3
    nblocks = nthreads()
    T = Float64
    @show T

    # -------------------------------------------------------------------------------------
    # Block version
    # -------------------------------------------------------------------------------------

    ρ_flat = rand(T, dims...)
    u_flat = rand(T, dims...)
    v_flat = rand(T, dims...)
    p_flat = rand(T, dims...)

    U = (ρ=BlockHaloArray(ρ_flat, nhalo, nblocks),
        u=BlockHaloArray(u_flat, nhalo, nblocks),
        v=BlockHaloArray(v_flat, nhalo, nblocks),
        p=BlockHaloArray(p_flat, nhalo, nblocks))

    edgedims = (2, dims...)
    halodims = (2, 3)

    i_edges = (ρ=BlockHaloArray(zeros(T, edgedims...), halodims, nhalo, nblocks),
        u=BlockHaloArray(zeros(T, edgedims...), halodims, nhalo, nblocks),
        v=BlockHaloArray(zeros(T, edgedims...), halodims, nhalo, nblocks),
        p=BlockHaloArray(zeros(T, edgedims...), halodims, nhalo, nblocks))

    j_edges = deepcopy(i_edges)

    # cache it
    reconstruct_with_blocks!(recon, U, i_edges, j_edges, 1)

    max_iter = 10
    global t_elaps_blocked = 0.0
    for _ in 1:max_iter
        global t_elaps_blocked += @elapsed begin
            @sync for thread_id in 1:nblocks
                @tspawnat thread_id reconstruct_with_blocks!(recon, U, i_edges, j_edges, thread_id)
            end
        end
    end
    t_block = t_elaps_blocked / max_iter


    # -------------------------------------------------------------------------------------
    # Flat version
    # -------------------------------------------------------------------------------------

    ρ_edge_flat = zeros(T, edgedims)
    u_edge_flat = zeros(T, edgedims)
    v_edge_flat = zeros(T, edgedims)
    p_edge_flat = zeros(T, edgedims)

    U_flat = (ρ=ρ_flat,
        u=u_flat,
        v=v_flat,
        p=p_flat)

    i_edges_flat = (ρ=ρ_edge_flat,
        u=u_edge_flat,
        v=v_edge_flat,
        p=p_edge_flat)

    j_edges_flat = deepcopy(i_edges_flat)


    ilo = first(axes(U_flat.ρ, 1)) + nhalo
    jlo = first(axes(U_flat.ρ, 2)) + nhalo
    ihi = last(axes(U_flat.ρ, 1)) - nhalo
    jhi = last(axes(U_flat.ρ, 2)) - nhalo
    looplimits = (ilo, ihi, jlo, jhi)

    # cache it
    reconstruct_flat!(recon, U_flat, i_edges_flat, j_edges_flat, looplimits)

    global t_elaps_flat = 0.0
    for _ in 1:max_iter
        global t_elaps_flat += @elapsed begin
            reconstruct_flat!(recon, U_flat, i_edges_flat, j_edges_flat, looplimits)
        end
    end
    t_flat = t_elaps_flat / max_iter

    @show t_flat, t_block

    open("results_5th_order.csv", "a") do io
        println(io, "$(nthreads()),$t_flat,$t_block")
    end

    nothing

end

run()