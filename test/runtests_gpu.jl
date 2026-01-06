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

    # p-based targets (match pe/ph/pθ internals: ρ(T) = air_density(T,p,q_tot))
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

    maxiter = 80
    tol = FT(1e-10)

    # 1) ρeq
    # 1) ρeq
    results = broadcast(d_ρ, d_e_int_ρ, d_q) do ρ, e_int, q
        TD.saturation_adjustment(
            RS.NewtonsMethod,
            param_set,
            TD.ρe(),
            ρ,
            e_int,
            q,
            maxiter,
            tol,
        )
    end
    d_T[1, :] .= getproperty.(results, :T)
    d_ql[1, :] .= getproperty.(results, :q_liq)
    d_qi[1, :] .= getproperty.(results, :q_ice)

    # 2) pe
    results = broadcast(d_p, d_e_int_p, d_q) do p, e_int, q
        TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.pe(),
            p,
            e_int,
            q,
            maxiter,
            tol,
        )
    end
    d_T[2, :] .= getproperty.(results, :T)
    d_ql[2, :] .= getproperty.(results, :q_liq)
    d_qi[2, :] .= getproperty.(results, :q_ice)

    # 3) ph
    results = broadcast(d_p, d_h_p, d_q) do p, h, q
        TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.ph(),
            p,
            h,
            q,
            maxiter,
            tol,
        )
    end
    d_T[3, :] .= getproperty.(results, :T)
    d_ql[3, :] .= getproperty.(results, :q_liq)
    d_qi[3, :] .= getproperty.(results, :q_ice)

    # 4) pρ
    results = broadcast(d_p_ρ, d_ρ, d_q) do p, ρ, q
        TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.pρ(),
            p,
            ρ,
            q,
            maxiter,
            tol,
        )
    end
    d_T[4, :] .= getproperty.(results, :T)
    d_ql[4, :] .= getproperty.(results, :q_liq)
    d_qi[4, :] .= getproperty.(results, :q_ice)

    # 5) ρθ_li
    results = broadcast(d_ρ, d_θ_ρ, d_q) do ρ, θ, q
        TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.ρθ_li(),
            ρ,
            θ,
            q,
            maxiter,
            tol,
        )
    end
    d_T[5, :] .= getproperty.(results, :T)
    d_ql[5, :] .= getproperty.(results, :q_liq)
    d_qi[5, :] .= getproperty.(results, :q_ice)

    # 6) pθ_li
    results = broadcast(d_p, d_θ_p, d_q) do p, θ, q
        TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.pθ_li(),
            p,
            θ,
            q,
            maxiter,
            tol,
        )
    end
    d_T[6, :] .= getproperty.(results, :T)
    d_ql[6, :] .= getproperty.(results, :q_liq)
    d_qi[6, :] .= getproperty.(results, :q_ice)

    # Compare device results with CPU reference, and check correctness invariants.
    T_gpu = Array(d_T)
    ql_gpu = Array(d_ql)
    qi_gpu = Array(d_qi)

    # Helper for tight-but-robust phase partition comparison.
    approx_tight(a, b) = isapprox(a, b; rtol = FT(100) * eps(FT), atol = FT(100) * eps(FT))

    for i in 1:n
        # CPU reference for each formulation
        res = TD.saturation_adjustment(
            RS.NewtonsMethod,
            param_set,
            TD.ρe(),
            ρ0[i],
            e_int_ρ[i],
            q0[i],
            80,
            FT(1e-10),
        )
        T_ref = res.T
        ql_ref = res.q_liq
        qi_ref = res.q_ice
        @test isapprox(T_gpu[1, i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_gpu[1, i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_gpu[1, i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (ρe): internal energy and phase partition at saturation
        @test isapprox(
            TD.internal_energy_sat(param_set, T_gpu[1, i], ρ0[i], q0[i]),
            e_int_ρ[i];
            rtol = FT(5e-6),
        )
        let (ql_chk, qi_chk) = TD.condensate_partition(param_set, T_gpu[1, i], ρ0[i], q0[i])
            @test approx_tight(ql_gpu[1, i], ql_chk)
            @test approx_tight(qi_gpu[1, i], qi_chk)
        end

        res = TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.pe(),
            p0[i],
            e_int_p[i],
            q0[i],
            80,
            FT(1e-10),
        )
        T_ref = res.T
        ql_ref = res.q_liq
        qi_ref = res.q_ice
        @test isapprox(T_gpu[2, i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_gpu[2, i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_gpu[2, i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (pe): energy and partition using ρ(T,p,q_tot)
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

        res = TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.ph(),
            p0[i],
            h_p[i],
            q0[i],
            80,
            FT(1e-10),
        )
        T_ref = res.T
        ql_ref = res.q_liq
        qi_ref = res.q_ice
        @test isapprox(T_gpu[3, i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_gpu[3, i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_gpu[3, i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (ph): enthalpy and partition using ρ(T,p,q_tot)
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

        res = TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.pρ(),
            p_ρ[i],
            ρ0[i],
            q0[i],
            80,
            FT(1e-10),
        )
        T_ref = res.T
        ql_ref = res.q_liq
        qi_ref = res.q_ice
        @test isapprox(T_gpu[4, i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_gpu[4, i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_gpu[4, i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (pρ): pressure target and partition at saturation
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

        res = TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.ρθ_li(),
            ρ0[i],
            θ_ρ[i],
            q0[i],
            80,
            FT(1e-10),
        )
        T_ref = res.T
        ql_ref = res.q_liq
        qi_ref = res.q_ice
        @test isapprox(T_gpu[5, i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_gpu[5, i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_gpu[5, i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (ρθ_li): θ target and partition at saturation
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

        res = TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.pθ_li(),
            p0[i],
            θ_p[i],
            q0[i],
            80,
            FT(1e-10),
        )
        T_ref = res.T
        ql_ref = res.q_liq
        qi_ref = res.q_ice
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
