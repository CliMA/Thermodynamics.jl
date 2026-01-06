using Test

import RootSolvers as RS

import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import ClimaParams as CP

include("TestedProfiles.jl")

# Determine array type.
# - If ARGS[1] is "CuArray" or "Array", honor it.
# - Otherwise, default to GPU when CUDA is available and functional.
arg = get(ARGS, 1, "")
if arg == "Array"
    ArrayType = Array
elseif arg == "CuArray"
    import CUDA
    ArrayType = CUDA.CuArray
    CUDA.allowscalar(false)
else
    ArrayType = try
        import CUDA
        if CUDA.functional()
            CUDA.allowscalar(false)
            CUDA.CuArray
        else
            Array
        end
    catch
        Array
    end
end

# Pre-allocate parameter sets for different floating-point types.
# Use ClimaParams to construct the TOML dict explicitly to avoid relying on extensions
# being loaded implicitly in this standalone test runner.
const param_set_Float64 = TP.ThermodynamicsParameters(CP.create_toml_dict(Float64))
const param_set_Float32 = TP.ThermodynamicsParameters(CP.create_toml_dict(Float32))

"""
    parameter_set(::Type{FT})

Get the appropriate parameter set for a given floating-point type.

# Arguments
- `FT`: Floating-point type (Float32 or Float64)

# Returns
- Thermodynamic parameter set for the specified type
"""
parameter_set(::Type{Float64}) = param_set_Float64
parameter_set(::Type{Float32}) = param_set_Float32

@show ArrayType

"""
    sat_adjust!(dst_T, dst_ql, dst_qi, inputs..., maxiter, tol, param_set)

CPU/GPU implementation for testing the functional `saturation_adjustment` API.

For each point, this kernel runs `saturation_adjustment` for multiple independent-variable
formulations and stores `(T, q_liq, q_ice)` for each.
"""
@inline function sat_adjust_point!(
    param_set,
    dst_T,
    dst_ql,
    dst_qi,
    ρ,
    p,
    p_ρ,
    q_tot,
    e_int_ρ,
    e_int_p,
    h_p,
    θ_ρ,
    θ_p,
    i,
    maxiter,
    tol,
)
    # 1) ρeq
    (T, ql, qi) = TD.saturation_adjustment(
        RS.NewtonsMethod,
        param_set,
        TD.ρeq(),
        ρ[i],
        e_int_ρ[i],
        q_tot[i],
        maxiter,
        tol,
    )
    dst_T[1, i] = T
    dst_ql[1, i] = ql
    dst_qi[1, i] = qi

    # 2) peq
    (T, ql, qi) = TD.saturation_adjustment(
        RS.SecantMethod,
        param_set,
        TD.peq(),
        p[i],
        e_int_p[i],
        q_tot[i],
        maxiter,
        tol,
    )
    dst_T[2, i] = T
    dst_ql[2, i] = ql
    dst_qi[2, i] = qi

    # 3) phq
    (T, ql, qi) = TD.saturation_adjustment(
        RS.SecantMethod,
        param_set,
        TD.phq(),
        p[i],
        h_p[i],
        q_tot[i],
        maxiter,
        tol,
    )
    dst_T[3, i] = T
    dst_ql[3, i] = ql
    dst_qi[3, i] = qi

    # 4) pρq
    (T, ql, qi) = TD.saturation_adjustment(
        RS.SecantMethod,
        param_set,
        TD.pρq(),
        p_ρ[i],
        ρ[i],
        q_tot[i],
        maxiter,
        tol,
    )
    dst_T[4, i] = T
    dst_ql[4, i] = ql
    dst_qi[4, i] = qi

    # 5) ρθ_liq_ice_q
    (T, ql, qi) = TD.saturation_adjustment(
        RS.SecantMethod,
        param_set,
        TD.ρθ_liq_ice_q(),
        ρ[i],
        θ_ρ[i],
        q_tot[i],
        maxiter,
        tol,
    )
    dst_T[5, i] = T
    dst_ql[5, i] = ql
    dst_qi[5, i] = qi

    # 6) pθ_liq_ice_q
    (T, ql, qi) = TD.saturation_adjustment(
        RS.SecantMethod,
        param_set,
        TD.pθ_liq_ice_q(),
        p[i],
        θ_p[i],
        q_tot[i],
        maxiter,
        tol,
    )
    dst_T[6, i] = T
    dst_ql[6, i] = ql
    dst_qi[6, i] = qi

    return nothing
end

function sat_adjust_cpu!(
    param_set,
    dst_T,
    dst_ql,
    dst_qi,
    ρ,
    p,
    p_ρ,
    q_tot,
    e_int_ρ,
    e_int_p,
    h_p,
    θ_ρ,
    θ_p,
    maxiter,
    tol,
)
    n = length(q_tot)
    @inbounds for i in 1:n
        sat_adjust_point!(
            param_set,
            dst_T,
            dst_ql,
            dst_qi,
            ρ,
            p,
            p_ρ,
            q_tot,
            e_int_ρ,
            e_int_p,
            h_p,
            θ_ρ,
            θ_p,
            i,
            maxiter,
            tol,
        )
    end
    return nothing
end

if get(ARGS, 1, "Array") == "CuArray"
    import CUDA
    # Define CUDA-kernel code only when CUDA is actually available.
    # (`CUDA.@cuda` is a macro and must not be expanded unless CUDA is loaded.)
    @eval begin
        function sat_adjust_cuda!(
            param_set,
            dst_T,
            dst_ql,
            dst_qi,
            ρ,
            p,
            p_ρ,
            q_tot,
            e_int_ρ,
            e_int_p,
            h_p,
            θ_ρ,
            θ_p,
            maxiter,
            tol,
        )
            function kernel!(
                param_set,
                dst_T,
                dst_ql,
                dst_qi,
                ρ,
                p,
                p_ρ,
                q_tot,
                e_int_ρ,
                e_int_p,
                h_p,
                θ_ρ,
                θ_p,
                maxiter,
                tol,
            )
                i = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
                if i <= length(q_tot)
                    @inbounds sat_adjust_point!(
                        param_set,
                        dst_T,
                        dst_ql,
                        dst_qi,
                        ρ,
                        p,
                        p_ρ,
                        q_tot,
                        e_int_ρ,
                        e_int_p,
                        h_p,
                        θ_ρ,
                        θ_p,
                        i,
                        maxiter,
                        tol,
                    )
                end
                return nothing
            end

            threads = 256
            blocks = cld(length(q_tot), threads)
            CUDA.@cuda threads = threads blocks = blocks kernel!(
                param_set,
                dst_T,
                dst_ql,
                dst_qi,
                ρ,
                p,
                p_ρ,
                q_tot,
                e_int_ρ,
                e_int_p,
                h_p,
                θ_ρ,
                θ_p,
                maxiter,
                tol,
            )
            CUDA.synchronize()
            return nothing
        end
    end
end

@testset "Thermodynamics - kernels" begin
    FT = Float32  # Use Float32 for GPU compatibility
    param_set = parameter_set(FT)

    # Create equilibrium-consistent test inputs on CPU, then move to device.
    profiles = TestedProfiles.PhaseEquilProfiles(param_set, Array{FT})
    n = min(length(profiles.z), 256)
    sl = 1:n
    T0 = profiles.T[sl]
    ρ0 = profiles.ρ[sl]
    q0 = profiles.q_tot[sl]
    q_liq0 = profiles.q_liq[sl]
    q_ice0 = profiles.q_ice[sl]

    # ρ-based targets
    e_int_ρ = TD.internal_energy_sat.(Ref(param_set), T0, ρ0, q0)
    p_ρ = TD.air_pressure.(Ref(param_set), T0, ρ0, q0, q_liq0, q_ice0)
    θ_ρ = TD.liquid_ice_pottemp.(Ref(param_set), T0, ρ0, q0, q_liq0, q_ice0)

    # p-based targets (match peq/phq/pθ internals: ρ(T) = air_density(T,p,q_tot))
    p0 = profiles.p[sl]
    ρ_p = TD.air_density.(Ref(param_set), T0, p0, q0)
    (q_liq_p, q_ice_p) =
        TD.condensate_partition.(Ref(param_set), T0, ρ_p, q0) |> x -> (first.(x), last.(x))
    e_int_p = TD.internal_energy_sat.(Ref(param_set), T0, ρ_p, q0)
    h_p = TD.enthalpy_sat.(Ref(param_set), T0, ρ_p, q0)
    θ_p =
        TD.liquid_ice_pottemp_given_pressure.(Ref(param_set), T0, p0, q0, q_liq_p, q_ice_p)

    # Move arrays to device
    d_ρ = ArrayType(ρ0)
    d_p = ArrayType(p0)
    d_p_ρ = ArrayType(p_ρ)
    d_q = ArrayType(q0)
    d_e_int_ρ = ArrayType(e_int_ρ)
    d_e_int_p = ArrayType(e_int_p)
    d_h_p = ArrayType(h_p)
    d_θ_ρ = ArrayType(θ_ρ)
    d_θ_p = ArrayType(θ_p)

    d_T = ArrayType(zeros(FT, 6, n))
    d_ql = ArrayType(zeros(FT, 6, n))
    d_qi = ArrayType(zeros(FT, 6, n))

    if ArrayType === Array
        sat_adjust_cpu!(
            param_set,
            d_T,
            d_ql,
            d_qi,
            d_ρ,
            d_p,
            d_p_ρ,
            d_q,
            d_e_int_ρ,
            d_e_int_p,
            d_h_p,
            d_θ_ρ,
            d_θ_p,
            80,
            FT(1e-10),
        )
    else
        sat_adjust_cuda!(
            param_set,
            d_T,
            d_ql,
            d_qi,
            d_ρ,
            d_p,
            d_p_ρ,
            d_q,
            d_e_int_ρ,
            d_e_int_p,
            d_h_p,
            d_θ_ρ,
            d_θ_p,
            80,
            FT(1e-10),
        )
    end

    # Compare device results with CPU reference, and check correctness invariants.
    T_gpu = Array(d_T)
    ql_gpu = Array(d_ql)
    qi_gpu = Array(d_qi)

    # Helper for tight-but-robust phase partition comparison.
    approx_tight(a, b) = isapprox(a, b; rtol = FT(100) * eps(FT), atol = FT(100) * eps(FT))

    for i in 1:n
        # CPU reference for each formulation
        (T_ref, ql_ref, qi_ref) = TD.saturation_adjustment(
            RS.NewtonsMethod,
            param_set,
            TD.ρeq(),
            ρ0[i],
            e_int_ρ[i],
            q0[i],
            80,
            FT(1e-10),
        )
        @test isapprox(T_gpu[1, i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_gpu[1, i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_gpu[1, i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (ρeq): internal energy and phase partition at saturation
        @test isapprox(
            TD.internal_energy_sat(param_set, T_gpu[1, i], ρ0[i], q0[i]),
            e_int_ρ[i];
            rtol = FT(5e-6),
        )
        let (ql_chk, qi_chk) = TD.condensate_partition(param_set, T_gpu[1, i], ρ0[i], q0[i])
            @test approx_tight(ql_gpu[1, i], ql_chk)
            @test approx_tight(qi_gpu[1, i], qi_chk)
        end

        (T_ref, ql_ref, qi_ref) = TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.peq(),
            p0[i],
            e_int_p[i],
            q0[i],
            80,
            FT(1e-10),
        )
        @test isapprox(T_gpu[2, i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_gpu[2, i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_gpu[2, i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (peq): energy and partition using ρ(T,p,q_tot)
        let ρ_eff = TD.air_density(param_set, T_gpu[2, i], p0[i], q0[i])
            @test isapprox(
                TD.internal_energy_sat(param_set, T_gpu[2, i], ρ_eff, q0[i]),
                e_int_p[i];
                rtol = FT(5e-6),
            )
            let (ql_chk, qi_chk) =
                    TD.condensate_partition(param_set, T_gpu[2, i], ρ_eff, q0[i])
                @test approx_tight(ql_gpu[2, i], ql_chk)
                @test approx_tight(qi_gpu[2, i], qi_chk)
            end
        end

        (T_ref, ql_ref, qi_ref) = TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.phq(),
            p0[i],
            h_p[i],
            q0[i],
            80,
            FT(1e-10),
        )
        @test isapprox(T_gpu[3, i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_gpu[3, i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_gpu[3, i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (phq): enthalpy and partition using ρ(T,p,q_tot)
        let ρ_eff = TD.air_density(param_set, T_gpu[3, i], p0[i], q0[i])
            @test isapprox(
                TD.enthalpy_sat(param_set, T_gpu[3, i], ρ_eff, q0[i]),
                h_p[i];
                rtol = FT(5e-6),
            )
            let (ql_chk, qi_chk) =
                    TD.condensate_partition(param_set, T_gpu[3, i], ρ_eff, q0[i])
                @test approx_tight(ql_gpu[3, i], ql_chk)
                @test approx_tight(qi_gpu[3, i], qi_chk)
            end
        end

        (T_ref, ql_ref, qi_ref) = TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.pρq(),
            p_ρ[i],
            ρ0[i],
            q0[i],
            80,
            FT(1e-10),
        )
        @test isapprox(T_gpu[4, i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_gpu[4, i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_gpu[4, i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (pρq): pressure target and partition at saturation
        @test isapprox(
            TD.air_pressure(
                param_set,
                T_gpu[4, i],
                ρ0[i],
                q0[i],
                ql_gpu[4, i],
                qi_gpu[4, i],
            ),
            p_ρ[i];
            rtol = FT(5e-6),
        )
        let (ql_chk, qi_chk) = TD.condensate_partition(param_set, T_gpu[4, i], ρ0[i], q0[i])
            @test approx_tight(ql_gpu[4, i], ql_chk)
            @test approx_tight(qi_gpu[4, i], qi_chk)
        end

        (T_ref, ql_ref, qi_ref) = TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.ρθ_liq_ice_q(),
            ρ0[i],
            θ_ρ[i],
            q0[i],
            80,
            FT(1e-10),
        )
        @test isapprox(T_gpu[5, i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_gpu[5, i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_gpu[5, i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (ρθ_liq_ice_q): θ target and partition at saturation
        @test isapprox(
            TD.liquid_ice_pottemp(
                param_set,
                T_gpu[5, i],
                ρ0[i],
                q0[i],
                ql_gpu[5, i],
                qi_gpu[5, i],
            ),
            θ_ρ[i];
            rtol = FT(5e-6),
        )
        let (ql_chk, qi_chk) = TD.condensate_partition(param_set, T_gpu[5, i], ρ0[i], q0[i])
            @test approx_tight(ql_gpu[5, i], ql_chk)
            @test approx_tight(qi_gpu[5, i], qi_chk)
        end

        (T_ref, ql_ref, qi_ref) = TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.pθ_liq_ice_q(),
            p0[i],
            θ_p[i],
            q0[i],
            80,
            FT(1e-10),
        )
        @test isapprox(T_gpu[6, i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_gpu[6, i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_gpu[6, i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (pθ_liq_ice_q): θ target and partition using ρ(T,p,q_tot)
        let ρ_eff = TD.air_density(param_set, T_gpu[6, i], p0[i], q0[i])
            let (ql_chk, qi_chk) =
                    TD.condensate_partition(param_set, T_gpu[6, i], ρ_eff, q0[i])
                @test approx_tight(ql_gpu[6, i], ql_chk)
                @test approx_tight(qi_gpu[6, i], qi_chk)
                @test isapprox(
                    TD.liquid_ice_pottemp_given_pressure(
                        param_set,
                        T_gpu[6, i],
                        p0[i],
                        q0[i],
                        ql_chk,
                        qi_chk,
                    ),
                    θ_p[i];
                    rtol = FT(5e-6),
                )
            end
        end

        # Basic invariants for all formulations
        @test all(ql_gpu[:, i] .>= 0)
        @test all(qi_gpu[:, i] .>= 0)
    end

    @testset "Broadcasting on ArrayType ($ArrayType)" begin
        # Broadcasting checks: run a few key pure functions on the target ArrayType
        # and compare against CPU. This catches GPU broadcast/kernel lowering issues.
        dT0 = ArrayType(T0)
        dρ0 = ArrayType(ρ0)
        dq0 = ArrayType(q0)
        dql0 = ArrayType(q_liq0)
        dqi0 = ArrayType(q_ice0)

        # air_pressure broadcast
        p_dev = TD.air_pressure.(Ref(param_set), dT0, dρ0, dq0, dql0, dqi0)
        p_cpu = TD.air_pressure.(Ref(param_set), T0, ρ0, q0, q_liq0, q_ice0)
        @test all(Array(p_dev) .≈ p_cpu)

        # internal_energy broadcast
        e_dev = TD.internal_energy.(Ref(param_set), dT0, dq0, dql0, dqi0)
        e_cpu = TD.internal_energy.(Ref(param_set), T0, q0, q_liq0, q_ice0)
        @test all(Array(e_dev) .≈ e_cpu)

        # liquid_ice_pottemp broadcast (ρ-based)
        θ_dev = TD.liquid_ice_pottemp.(Ref(param_set), dT0, dρ0, dq0, dql0, dqi0)
        θ_cpu = TD.liquid_ice_pottemp.(Ref(param_set), T0, ρ0, q0, q_liq0, q_ice0)
        @test all(Array(θ_dev) .≈ θ_cpu)
    end
end
