using .Threads
using ThreadPinning: pinthreads
using BlockHaloArrays
using BlockHaloArrays: @tspawnat
using LoopVectorization
using TimerOutputs

abstract type AbstractReconstructionScheme end
abstract type AbstractMUSCLReconstruction <: AbstractReconstructionScheme end
abstract type AbstractLimiter end

const to = TimerOutput()
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

struct MLP5Reconstruction <: AbstractMUSCLReconstruction end

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

function reconstruct_with_blocks!(::MLP5Reconstruction, U, i_edge, j_edge, blockid, ϵ=eps(Float64))

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

function reconstruct_with_blocks_v2!(::MLP5Reconstruction, U, i_edge, j_edge, blockid, ϵ=eps(Float64))

    # i_edge & j_edge holds the i+1/2 & j+1/2 L and R values
    # these limits must be the no-halo versions. If the true array size is (10,10) with
    # 2 halo cells, the limits must be ilo, ihi = 3, 8 and jlo, jhi = 3, 8
    
    limits = Tuple(first(U).loop_limits[blockid])

    for (ϕ̄, ϕ_recon_i, ϕ_recon_j) in zip(U, i_edge, j_edge) # (ρ, u, v, p)
        A = @views ϕ̄[blockid]
        ϕ_recon_i_edge = @views ϕ_recon_i[blockid]
        ϕ_recon_j_edge = @views ϕ_recon_j[blockid]
        reconstruct_mlp5!(A, ϕ_recon_i_edge, ϕ_recon_j_edge, limits, ϵ)
    end
    return nothing
end

function reconstruct_mlp5!(A, i_edge, j_edge, limits::NTuple{4, Int}, ϵ=eps(Float64))

    # i_edge & j_edge holds the i+1/2 & j+1/2 L and R values
    # these limits must be the no-halo versions. If the true array size is (10,10) with
    # 2 halo cells, the limits must be ilo, ihi = 3, 8 and jlo, jhi = 3, 8
    ilo, ihi, jlo, jhi = limits

    @turbo for j = jlo:jhi
        for i = ilo-1:ihi
            Δ_minus_half = A[i, j] - A[i-1, j]
            Δ_plus_three_half = A[i+2, j] - A[i+1, j]
            Δ_plus_half = A[i+1, j] - A[i, j]
            Δ_plus_five_half = A[i+3, j] - A[i+2, j]
            Δ_minus_three_half = A[i-1, j] - A[i-2, j]

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

            tanθ_ij = abs(A[i, j+1] - A[i, j-1]) / (abs(A[i+1, j] - A[i-1, j]) + SMALL_NUM)
            tanθ_ij_p1 = abs(A[i+1, j+1] - A[i+1, j-1]) / (abs(A[i+2, j] - A[i, j]) + SMALL_NUM)
            α_L, α_R = get_alpha(rL, rR_p1, tanθ_ij, tanθ_ij_p1)

            i_edge[1, i, j] = A[i, j] + 0.5max(0, min(α_L * rL, α_L, β_L)) * Δ_minus_half
            i_edge[2, i, j] = A[i+1, j] - 0.5max(0, min(α_R * rR_p1, α_R, β_R)) * Δ_plus_three_half
        end
    end

    @turbo for j = jlo-1:jhi
        for i = ilo:ihi
            Δ_minus_half = A[i, j] - A[i, j-1]
            Δ_plus_three_half = A[i, j+2] - A[i, j+1]
            Δ_plus_half = A[i, j+1] - A[i, j]
            Δ_minus_three_half = A[i, j-1] - A[i, j-2]
            Δ_plus_five_half = A[i, j+3] - A[i, j+3]

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

            tanθ_ij = abs(A[i+1, j] - A[i-1, j]) / (abs(A[i, j+1] - A[i, j-1]) + SMALL_NUM)
            tanθ_ij_p1 = abs(A[i+1, j+1] - A[i-1, j+1]) / (abs(A[i, j+2] - A[i, j]) + SMALL_NUM)
            α_L, α_R = get_alpha(rL, rR_p1, tanθ_ij, tanθ_ij_p1)

            j_edge[1, i, j] = A[i, j] + 0.5max(0, min(α_L * rL, α_L, β_L)) * Δ_minus_half
            j_edge[2, i, j] = A[i, j+1] - 0.5max(0, min(α_R * rR_p1, α_R, β_R)) * Δ_plus_three_half
        end
    end

    return nothing
end

function reconstruct_flat!(::MLP5Reconstruction, U, i_edge, j_edge, looplimits, ϵ=eps(Float64))

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

function run_stencil_benchmark(recon, nhalo, T, results_file)

    @show nthreads()

    pinthreads(:compact)

    dims = (512, 1024) .* 10
    nblocks = nthreads()
    @show dims
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
    for _ in 1:max_iter
        @timeit to "reconstruct_with_blocks!" begin
            @sync for thread_id in 1:nblocks
                @tspawnat thread_id reconstruct_with_blocks!(recon, U, i_edges, j_edges, thread_id)
            end
        end
    end

    # cache it
    reconstruct_with_blocks_v2!(recon, U, i_edges, j_edges, 1)

    max_iter = 10
    for _ in 1:max_iter
        @timeit to "reconstruct_with_blocks_v2!" begin
            @sync for thread_id in 1:nblocks
                @tspawnat thread_id reconstruct_with_blocks_v2!(recon, U, i_edges, j_edges, thread_id)
            end
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

    for _ in 1:max_iter
        @timeit to "reconstruct_flat!" begin
            reconstruct_flat!(recon, U_flat, i_edges_flat, j_edges_flat, looplimits)
        end
    end

    show(to)
    println()

    t_flat = TimerOutputs.time(to["reconstruct_flat!"])
    t_block = TimerOutputs.time(to["reconstruct_with_blocks!"])

    println("t_flat / t_block = $(t_flat / t_block)")

    open(results_file, "a") do io
        println(io, "$(nthreads()),$t_flat,$t_block")
    end

    nothing
end