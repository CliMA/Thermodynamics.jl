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

if ArrayType == Array
    @info "Running GPU tests on CPU (ArrayType = Array)"
else
    @info "Running GPU tests on GPU (ArrayType = $ArrayType)"
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

    maxiter = 80
    tol = FT(1e-10)

    # GPU tests ensure that the GPU path produces the same physical outputs as the 
    # CPU reference implementation for the same states.
    # Use default methods (no method type parameter) to avoid dynamic dispatch issues on GPU.

    # 1) ρeq - using NewtonsMethod (via type dispatch for GPU compatibility)
    gpu_results_ρe =
        TD.saturation_adjustment.(
            Ref(RS.NewtonsMethod),
            Ref(param_set),
            Ref(TD.ρe()),
            d_ρ,
            d_e_int_ρ,
            d_q,
            Ref(maxiter),
            Ref(tol),
        )
    T_ρe = map(x -> x.T, gpu_results_ρe)
    ql_ρe = map(x -> x.q_liq, gpu_results_ρe)
    qi_ρe = map(x -> x.q_ice, gpu_results_ρe)

    # 2) pe - using SecantMethod (via type dispatch for GPU compatibility)
    gpu_results_pe =
        TD.saturation_adjustment.(
            Ref(RS.SecantMethod),
            Ref(param_set),
            Ref(TD.pe()),
            d_p,
            d_e_int_p,
            d_q,
            Ref(maxiter),
            Ref(tol),
        )
    T_pe = map(x -> x.T, gpu_results_pe)
    ql_pe = map(x -> x.q_liq, gpu_results_pe)
    qi_pe = map(x -> x.q_ice, gpu_results_pe)

    # 3) ph - using SecantMethod (via type dispatch for GPU compatibility)
    gpu_results_ph =
        TD.saturation_adjustment.(
            Ref(RS.SecantMethod),
            Ref(param_set),
            Ref(TD.ph()),
            d_p,
            d_h_p,
            d_q,
            Ref(maxiter),
            Ref(tol),
        )
    T_ph = map(x -> x.T, gpu_results_ph)
    ql_ph = map(x -> x.q_liq, gpu_results_ph)
    qi_ph = map(x -> x.q_ice, gpu_results_ph)

    # 4) pρ - using SecantMethod (via type dispatch for GPU compatibility)
    gpu_results_pρ =
        TD.saturation_adjustment.(
            Ref(RS.SecantMethod),
            Ref(param_set),
            Ref(TD.pρ()),
            d_p_ρ,
            d_ρ,
            d_q,
            Ref(maxiter),
            Ref(tol),
        )
    T_pρ = map(x -> x.T, gpu_results_pρ)
    ql_pρ = map(x -> x.q_liq, gpu_results_pρ)
    qi_pρ = map(x -> x.q_ice, gpu_results_pρ)

    # 5) ρθ_li - using SecantMethod (via type dispatch for GPU compatibility)
    gpu_results_ρθ_li =
        TD.saturation_adjustment.(
            Ref(RS.SecantMethod),
            Ref(param_set),
            Ref(TD.ρθ_li()),
            d_ρ,
            d_θ_ρ,
            d_q,
            Ref(maxiter),
            Ref(tol),
        )
    T_ρθ_li = map(x -> x.T, gpu_results_ρθ_li)
    ql_ρθ_li = map(x -> x.q_liq, gpu_results_ρθ_li)
    qi_ρθ_li = map(x -> x.q_ice, gpu_results_ρθ_li)

    # 6) pθ_li - using SecantMethod (via type dispatch for GPU compatibility)
    gpu_results_pθ_li =
        TD.saturation_adjustment.(
            Ref(RS.SecantMethod),
            Ref(param_set),
            Ref(TD.pθ_li()),
            d_p,
            d_θ_p,
            d_q,
            Ref(maxiter),
            Ref(tol),
        )
    T_pθ_li = map(x -> x.T, gpu_results_pθ_li)
    ql_pθ_li = map(x -> x.q_liq, gpu_results_pθ_li)
    qi_pθ_li = map(x -> x.q_ice, gpu_results_pθ_li)

    # Convert GPU results to CPU arrays for comparison
    T_ρe_cpu = Array(T_ρe)
    ql_ρe_cpu = Array(ql_ρe)
    qi_ρe_cpu = Array(qi_ρe)
    T_pe_cpu = Array(T_pe)
    ql_pe_cpu = Array(ql_pe)
    qi_pe_cpu = Array(qi_pe)
    T_ph_cpu = Array(T_ph)
    ql_ph_cpu = Array(ql_ph)
    qi_ph_cpu = Array(qi_ph)
    T_pρ_cpu = Array(T_pρ)
    ql_pρ_cpu = Array(ql_pρ)
    qi_pρ_cpu = Array(qi_pρ)
    T_ρθ_li_cpu = Array(T_ρθ_li)
    ql_ρθ_li_cpu = Array(ql_ρθ_li)
    qi_ρθ_li_cpu = Array(qi_ρθ_li)
    T_pθ_li_cpu = Array(T_pθ_li)
    ql_pθ_li_cpu = Array(ql_pθ_li)
    qi_pθ_li_cpu = Array(qi_pθ_li)

    # Helper for tight-but-robust phase partition comparison.
    approx_tight(a, b) = isapprox(a, b; rtol = FT(100) * eps(FT), atol = FT(100) * eps(FT))

    for i in 1:n
        # CPU reference for each formulation (using explicit methods for CPU)
        res = TD.saturation_adjustment(
            RS.NewtonsMethod,
            param_set,
            TD.ρe(),
            ρ0[i],
            e_int_ρ[i],
            q0[i],
            maxiter,
            tol,
        )
        T_ref = res.T
        ql_ref = res.q_liq
        qi_ref = res.q_ice
        @test isapprox(T_ρe_cpu[i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_ρe_cpu[i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_ρe_cpu[i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (ρe): internal energy and phase partition at saturation
        @test isapprox(
            TD.internal_energy_sat(param_set, T_ρe_cpu[i], ρ0[i], q0[i]),
            e_int_ρ[i];
            rtol = FT(5e-6),
        )
        let (ql_chk, qi_chk) = TD.condensate_partition(param_set, T_ρe_cpu[i], ρ0[i], q0[i])
            @test approx_tight(ql_ρe_cpu[i], ql_chk)
            @test approx_tight(qi_ρe_cpu[i], qi_chk)
        end

        res = TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.pe(),
            p0[i],
            e_int_p[i],
            q0[i],
            maxiter,
            tol,
        )
        T_ref = res.T
        ql_ref = res.q_liq
        qi_ref = res.q_ice
        @test isapprox(T_pe_cpu[i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_pe_cpu[i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_pe_cpu[i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (pe): energy and partition using ρ(T,p,q_tot)
        let ρ_eff = TD.air_density(param_set, T_pe_cpu[i], p0[i], q0[i])
            @test isapprox(
                TD.internal_energy_sat(param_set, T_pe_cpu[i], ρ_eff, q0[i]),
                e_int_p[i];
                rtol = FT(5e-6),
            )
            let (ql_chk, qi_chk) =
                    TD.condensate_partition(param_set, T_pe_cpu[i], ρ_eff, q0[i])
                @test approx_tight(ql_pe_cpu[i], ql_chk)
                @test approx_tight(qi_pe_cpu[i], qi_chk)
            end
        end

        res = TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.ph(),
            p0[i],
            h_p[i],
            q0[i],
            maxiter,
            tol,
        )
        T_ref = res.T
        ql_ref = res.q_liq
        qi_ref = res.q_ice
        @test isapprox(T_ph_cpu[i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_ph_cpu[i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_ph_cpu[i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (ph): enthalpy and partition using ρ(T,p,q_tot)
        let ρ_eff = TD.air_density(param_set, T_ph_cpu[i], p0[i], q0[i])
            @test isapprox(
                TD.enthalpy_sat(param_set, T_ph_cpu[i], ρ_eff, q0[i]),
                h_p[i];
                rtol = FT(5e-6),
            )
            let (ql_chk, qi_chk) =
                    TD.condensate_partition(param_set, T_ph_cpu[i], ρ_eff, q0[i])
                @test approx_tight(ql_ph_cpu[i], ql_chk)
                @test approx_tight(qi_ph_cpu[i], qi_chk)
            end
        end

        res = TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.pρ(),
            p_ρ[i],
            ρ0[i],
            q0[i],
            maxiter,
            tol,
        )
        T_ref = res.T
        ql_ref = res.q_liq
        qi_ref = res.q_ice
        @test isapprox(T_pρ_cpu[i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_pρ_cpu[i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_pρ_cpu[i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (pρ): pressure target and partition at saturation
        @test isapprox(
            TD.air_pressure(
                param_set,
                T_pρ_cpu[i],
                ρ0[i],
                q0[i],
                ql_pρ_cpu[i],
                qi_pρ_cpu[i],
            ),
            p_ρ[i];
            rtol = FT(5e-6),
        )
        let (ql_chk, qi_chk) = TD.condensate_partition(param_set, T_pρ_cpu[i], ρ0[i], q0[i])
            @test approx_tight(ql_pρ_cpu[i], ql_chk)
            @test approx_tight(qi_pρ_cpu[i], qi_chk)
        end

        res = TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.ρθ_li(),
            ρ0[i],
            θ_ρ[i],
            q0[i],
            maxiter,
            tol,
        )
        T_ref = res.T
        ql_ref = res.q_liq
        qi_ref = res.q_ice
        @test isapprox(T_ρθ_li_cpu[i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_ρθ_li_cpu[i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_ρθ_li_cpu[i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (ρθ_li): θ target and partition at saturation
        @test isapprox(
            TD.liquid_ice_pottemp(
                param_set,
                T_ρθ_li_cpu[i],
                ρ0[i],
                q0[i],
                ql_ρθ_li_cpu[i],
                qi_ρθ_li_cpu[i],
            ),
            θ_ρ[i];
            rtol = FT(5e-6),
        )
        let (ql_chk, qi_chk) =
                TD.condensate_partition(param_set, T_ρθ_li_cpu[i], ρ0[i], q0[i])
            @test approx_tight(ql_ρθ_li_cpu[i], ql_chk)
            @test approx_tight(qi_ρθ_li_cpu[i], qi_chk)
        end

        res = TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.pθ_li(),
            p0[i],
            θ_p[i],
            q0[i],
            maxiter,
            tol,
        )
        T_ref = res.T
        ql_ref = res.q_liq
        qi_ref = res.q_ice
        @test isapprox(T_pθ_li_cpu[i], T_ref; rtol = FT(5e-6))
        @test isapprox(ql_pθ_li_cpu[i], ql_ref; rtol = FT(1e-6), atol = FT(1e-12))
        @test isapprox(qi_pθ_li_cpu[i], qi_ref; rtol = FT(1e-6), atol = FT(1e-12))
        # Invariance checks (pθ_liq_ice_q): θ target and partition using ρ(T,p,q_tot)
        let ρ_eff = TD.air_density(param_set, T_pθ_li_cpu[i], p0[i], q0[i])
            let (ql_chk, qi_chk) =
                    TD.condensate_partition(param_set, T_pθ_li_cpu[i], ρ_eff, q0[i])
                @test approx_tight(ql_pθ_li_cpu[i], ql_chk)
                @test approx_tight(qi_pθ_li_cpu[i], qi_chk)
                @test isapprox(
                    TD.liquid_ice_pottemp_given_pressure(
                        param_set,
                        T_pθ_li_cpu[i],
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
        @test ql_ρe_cpu[i] >= 0
        @test qi_ρe_cpu[i] >= 0
        @test ql_pe_cpu[i] >= 0
        @test qi_pe_cpu[i] >= 0
        @test ql_ph_cpu[i] >= 0
        @test qi_ph_cpu[i] >= 0
        @test ql_pρ_cpu[i] >= 0
        @test qi_pρ_cpu[i] >= 0
        @test ql_ρθ_li_cpu[i] >= 0
        @test qi_ρθ_li_cpu[i] >= 0
        @test ql_pθ_li_cpu[i] >= 0
        @test qi_pθ_li_cpu[i] >= 0
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
