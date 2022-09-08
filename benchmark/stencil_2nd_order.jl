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

struct MinMod <: AbstractLimiter end
@inline limit(::MinMod, r) = max(0, min(r, 1))

struct VanLeer <: AbstractLimiter end
@inline limit(::VanLeer, r) = (r + abs(r)) / (1 + abs(r))

struct SuperBee <: AbstractLimiter end
@inline limit(::SuperBee, r) = max(0, min(2r, 1), min(r, 2))


struct MUSCLReconstruction{L<:AbstractLimiter} <: AbstractMUSCLReconstruction
    limiter::L
end

function reconstruct_with_blocks!(RS::MUSCLReconstruction, U, i_edge, j_edge, blockid)

    # i_edge & j_edge holds the i+1/2 & j+1/2 L and R values
    # these limits must be the no-halo versions. If the true array size is (10,10) with
    # 2 halo cells, the limits must be ilo, ihi = 3, 8 and jlo, jhi = 3, 8
    ilo, ihi, jlo, jhi = first(U).loop_limits[blockid]
    ϵ = eps(Float64)

    for (ϕ̄, ϕ_recon) in zip(U, i_edge) # (ρ, u, v, p)
        @turbo for j = jlo:jhi
            for i = ilo-1:ihi
                ϕ̄ᵢⱼ = ϕ̄[blockid][i, j]
                ϕ̄ᵢ₊₁ⱼ = ϕ̄[blockid][i+1, j]

                Δ_minus_half = ϕ̄ᵢⱼ - ϕ̄[blockid][i-1, j]
                Δ_plus_half = ϕ̄ᵢ₊₁ⱼ - ϕ̄ᵢⱼ
                Δ_plus_three_half = ϕ̄[blockid][i+2, j] - ϕ̄ᵢ₊₁ⱼ

                # Machine epsilon checks
                Δ_minus_half = Δ_minus_half * (abs(Δ_minus_half) >= ϵ)
                Δ_plus_half = Δ_plus_half * (abs(Δ_plus_half) >= ϵ)
                Δ_plus_three_half = Δ_plus_three_half * (abs(Δ_plus_three_half) >= ϵ)

                rL = Δ_plus_half / (Δ_minus_half + 1e-30)
                rR = Δ_plus_half / (Δ_plus_three_half + 1e-30)

                lim_L = limit(RS.limiter, rL)
                lim_R = limit(RS.limiter, rR)

                ϕ_recon[blockid][1, i, j] = ϕ̄ᵢⱼ + 0.5lim_L * Δ_minus_half
                ϕ_recon[blockid][2, i, j] = ϕ̄ᵢ₊₁ⱼ - 0.5lim_R * Δ_plus_three_half
            end
        end
    end

    for (ϕ̄, ϕ_recon) in zip(U, j_edge) # (ρ, u, v, p)
        @turbo for j = jlo-1:jhi
            for i = ilo:ihi
                ϕ̄ᵢⱼ = ϕ̄[blockid][i, j]
                ϕ̄ᵢⱼ₊₁ = ϕ̄[blockid][i, j+1]

                Δ_minus_half = ϕ̄ᵢⱼ - ϕ̄[blockid][i, j-1]
                Δ_plus_half = ϕ̄ᵢⱼ₊₁ - ϕ̄ᵢⱼ
                Δ_plus_three_half = ϕ̄[blockid][i, j+2] - ϕ̄ᵢⱼ₊₁

                # Machine epsilon checks
                Δ_minus_half = Δ_minus_half * (abs(Δ_minus_half) >= ϵ)
                Δ_plus_half = Δ_plus_half * (abs(Δ_plus_half) >= ϵ)
                Δ_plus_three_half = Δ_plus_three_half * (abs(Δ_plus_three_half) >= ϵ)

                rL = Δ_plus_half / (Δ_minus_half + 1e-30)
                rR = Δ_plus_half / (Δ_plus_three_half + 1e-30)

                lim_L = limit(RS.limiter, rL)
                lim_R = limit(RS.limiter, rR)

                ϕ_recon[blockid][1, i, j] = ϕ̄ᵢⱼ + 0.5lim_L * Δ_minus_half
                ϕ_recon[blockid][2, i, j] = ϕ̄ᵢⱼ₊₁ - 0.5lim_R * Δ_plus_three_half
            end
        end
    end


    return nothing
end

function reconstruct_flat!(RS::MUSCLReconstruction, U, i_edge, j_edge, looplimits)

    # i_edge & j_edge holds the i+1/2 & j+1/2 L and R values
    # these limits must be the no-halo versions. If the true array size is (10,10) with
    # 2 halo cells, the limits must be ilo, ihi = 3, 8 and jlo, jhi = 3, 8
    ilo, ihi, jlo, jhi = looplimits
    # @show ilo, ihi, jlo, jhi
    ϵ = eps(Float64)

    for (ϕ̄, ϕ_recon) in zip(U, i_edge) # (ρ, u, v, p)
        @tturbo for j = jlo:jhi
            for i = ilo-1:ihi
                ϕ̄ᵢⱼ = ϕ̄[i, j]
                ϕ̄ᵢ₊₁ⱼ = ϕ̄[i+1, j]

                Δ_minus_half = ϕ̄ᵢⱼ - ϕ̄[i-1, j]
                Δ_plus_half = ϕ̄ᵢ₊₁ⱼ - ϕ̄ᵢⱼ
                Δ_plus_three_half = ϕ̄[i+2, j] - ϕ̄ᵢ₊₁ⱼ

                # Machine epsilon checks
                Δ_minus_half = Δ_minus_half * (abs(Δ_minus_half) >= ϵ)
                Δ_plus_half = Δ_plus_half * (abs(Δ_plus_half) >= ϵ)
                Δ_plus_three_half = Δ_plus_three_half * (abs(Δ_plus_three_half) >= ϵ)

                rL = Δ_plus_half / (Δ_minus_half + 1e-30)
                rR = Δ_plus_half / (Δ_plus_three_half + 1e-30)

                lim_L = limit(RS.limiter, rL)
                lim_R = limit(RS.limiter, rR)

                ϕ_recon[1, i, j] = ϕ̄ᵢⱼ + 0.5lim_L * Δ_minus_half
                ϕ_recon[2, i, j] = ϕ̄ᵢ₊₁ⱼ - 0.5lim_R * Δ_plus_three_half
            end
        end
    end

    for (ϕ̄, ϕ_recon) in zip(U, j_edge) # (ρ, u, v, p)
        @tturbo for j = jlo-1:jhi
            for i = ilo:ihi
                ϕ̄ᵢⱼ = ϕ̄[i, j]
                ϕ̄ᵢⱼ₊₁ = ϕ̄[i, j+1]

                Δ_minus_half = ϕ̄ᵢⱼ - ϕ̄[i, j-1]
                Δ_plus_half = ϕ̄ᵢⱼ₊₁ - ϕ̄ᵢⱼ
                Δ_plus_three_half = ϕ̄[i, j+2] - ϕ̄ᵢⱼ₊₁

                # Machine epsilon checks
                Δ_minus_half = Δ_minus_half * (abs(Δ_minus_half) >= ϵ)
                Δ_plus_half = Δ_plus_half * (abs(Δ_plus_half) >= ϵ)
                Δ_plus_three_half = Δ_plus_three_half * (abs(Δ_plus_three_half) >= ϵ)

                rL = Δ_plus_half / (Δ_minus_half + 1e-30)
                rR = Δ_plus_half / (Δ_plus_three_half + 1e-30)

                lim_L = limit(RS.limiter, rL)
                lim_R = limit(RS.limiter, rR)

                ϕ_recon[1, i, j] = ϕ̄ᵢⱼ + 0.5lim_L * Δ_minus_half
                ϕ_recon[2, i, j] = ϕ̄ᵢⱼ₊₁ - 0.5lim_R * Δ_plus_three_half
            end
        end
    end


    return nothing
end

function run()

    @show nthreads()

    pinthreads(:compact)

    recon = MUSCLReconstruction(MinMod())

    dims = (512, 1024) .* 10
    @show dims
    nhalo = 2
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
    
    t_block = @elapsed begin
        @sync for thread_id in 1:nblocks
            @tspawnat thread_id reconstruct_with_blocks!(recon, U, i_edges, j_edges, thread_id)
        end
    end

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

    t_flat = @elapsed begin
        reconstruct_flat!(recon, U_flat, i_edges_flat, j_edges_flat, looplimits)
    end

    @show t_flat, t_block, t_flat / t_block

    nothing

end

run()