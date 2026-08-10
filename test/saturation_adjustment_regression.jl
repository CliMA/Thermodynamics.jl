"""
# Saturation adjustment regression tests

These pin behaviour that was previously wrong, so that it cannot silently regress. Each
testset names the failure mode it guards against.

The states used here are built to be *exactly* in equilibrium: the condensate is computed
from the same (T, ρ) or (T, p) the solver will be handed, so the solver's answer can be
compared against a known temperature rather than against a second call to the same code.
"""

# Build an exactly-equilibrium state at fixed density.
function equilibrium_state_ρ(param_set, T, p, q_tot)
    (q_liq, q_ice) = TD._condensate_partition_from_p(param_set, T, p, q_tot)
    ρ = TD.air_density(param_set, T, p, q_tot, q_liq, q_ice)
    # Re-partition at that density so the state is self-consistent in (T, ρ)
    (q_liq, q_ice) = TD.condensate_partition(param_set, T, ρ, q_tot)
    ρ = TD.air_density(param_set, T, p, q_tot, q_liq, q_ice)
    return (; ρ, q_liq, q_ice)
end

@testset "Thermodynamics - saturation adjustment regressions" begin
    for FT in (Float32, Float64)
        param_set = FT == Float64 ? param_set_Float64 : param_set_Float32

        # Cloud-heavy states spanning warm, mixed-phase and cold conditions. These carry
        # far more condensate than the TestedProfiles columns, which is the regime where
        # the solver used to fail.
        cloudy_states = [
            (T = FT(300), p = FT(100000), q_tot = FT(0.030)),
            (T = FT(295), p = FT(90000), q_tot = FT(0.025)),
            (T = FT(285), p = FT(90000), q_tot = FT(0.015)),
            (T = FT(275), p = FT(85000), q_tot = FT(0.008)),
            (T = FT(265), p = FT(80000), q_tot = FT(0.006)),
            (T = FT(255), p = FT(70000), q_tot = FT(0.004)),
            (T = FT(240), p = FT(60000), q_tot = FT(0.002)),
        ]

        @testset "Converged solves recover the state exactly ($FT)" begin
            # Before the pressure-based formulations took their partition from the pressure,
            # pe/ph/pθ_li converged to a temperature biased by ~0.16 K on states like these,
            # and no amount of iteration removed it. Tolerances here are tight on purpose.
            rtol = FT === Float32 ? FT(1e-4) : FT(1e-9)
            atol = FT === Float32 ? FT(2e-2) : FT(1e-6)
            maxiter = 100
            tol = FT === Float32 ? FT(1e-6) : FT(1e-12)

            for st in cloudy_states
                (; ρ, q_liq, q_ice) =
                    equilibrium_state_ρ(param_set, st.T, st.p, st.q_tot)
                e_int = TD.internal_energy(param_set, st.T, st.q_tot, q_liq, q_ice)
                h = TD.enthalpy(param_set, st.T, st.q_tot, q_liq, q_ice)
                θ_p = TD.liquid_ice_pottemp_given_pressure(
                    param_set,
                    st.T,
                    st.p,
                    st.q_tot,
                    q_liq,
                    q_ice,
                )
                θ_ρ = TD.liquid_ice_pottemp(
                    param_set,
                    st.T,
                    ρ,
                    st.q_tot,
                    q_liq,
                    q_ice,
                )
                p_ρ = TD.air_pressure(param_set, st.T, ρ, st.q_tot, q_liq, q_ice)

                cases = (
                    (TD.ρe(), ρ, e_int),
                    (TD.pe(), st.p, e_int),
                    (TD.ph(), st.p, h),
                    (TD.pρ(), p_ρ, ρ),
                    (TD.pθ_li(), st.p, θ_p),
                    (TD.ρθ_li(), ρ, θ_ρ),
                )
                for (indep_vars, var₁, var₂) in cases
                    sol = TD.saturation_adjustment(
                        RS.NewtonsMethod,
                        param_set,
                        indep_vars,
                        var₁,
                        var₂,
                        st.q_tot,
                        maxiter,
                        tol,
                    )
                    @test sol.converged
                    @test isapprox(sol.T, st.T; rtol = rtol, atol = atol)
                end
            end
        end

        @testset "Fixed iterations converge monotonically ($FT)" begin
            # The fixed-iteration solver could previously enter a limit cycle: the error did
            # not shrink with maxiter, and at the worst states it stayed tens of K off no
            # matter how many iterations were allowed. Guard the property that actually
            # matters now — more iterations never make the answer worse, and enough
            # iterations reach the exact answer.
            for st in cloudy_states
                (; ρ, q_liq, q_ice) =
                    equilibrium_state_ρ(param_set, st.T, st.p, st.q_tot)
                e_int = TD.internal_energy(param_set, st.T, st.q_tot, q_liq, q_ice)

                errs = map((2, 3, 4, 6, 10)) do m
                    sol = TD.saturation_adjustment(
                        param_set,
                        TD.ρe(),
                        ρ,
                        e_int,
                        st.q_tot;
                        maxiter = m,
                    )
                    abs(sol.T - st.T)
                end
                # Non-increasing to within round-off
                slack = FT === Float32 ? FT(1e-2) : FT(1e-8)
                for i in 2:length(errs)
                    @test errs[i] <= errs[i - 1] + slack
                end
                # And converged by ten iterations
                @test errs[end] < (FT === Float32 ? FT(1e-2) : FT(1e-6))
            end
        end

        @testset "Default iteration count matches its documented envelope ($FT)" begin
            # What the default `maxiter = 2` promises is set by the gap between the
            # unsaturated first guess and the solution — the warming from condensing the
            # excess vapor — not by the absolute humidity. The envelope below is the one
            # stated in the `saturation_adjustment` docstring, with headroom; the states
            # here span gaps from 6 K to 25 K, well beyond the sub-4 K range that the
            # tested profiles occupy.
            envelope(gap) =
                gap <= 12 ? FT(0.05) : gap <= 18 ? FT(0.5) : FT(5)

            for st in cloudy_states
                (; ρ, q_liq, q_ice) =
                    equilibrium_state_ρ(param_set, st.T, st.p, st.q_tot)
                e_int = TD.internal_energy(param_set, st.T, st.q_tot, q_liq, q_ice)
                T_unsat = TD.air_temperature(param_set, e_int, st.q_tot)
                gap = abs(st.T - T_unsat)
                sol = TD.saturation_adjustment(param_set, TD.ρe(), ρ, e_int, st.q_tot)
                @test abs(sol.T - st.T) < envelope(gap)
            end
        end

        @testset "Tested profiles are well inside that envelope ($FT)" begin
            # The realistic columns the package ships with sit at gaps below 4 K, where two
            # iterations are accurate to a few millikelvin. This is the regime a
            # time-stepping model actually operates in.
            profiles = TestedProfiles.EquilMoistProfiles(param_set, Array{FT})
            (; T, ρ, q_tot) = profiles
            idxs = unique(round.(Int, range(1, length(T), length = 100)))
            worst_gap = zero(FT)
            worst_err = zero(FT)
            for i in idxs
                e_int = TD.internal_energy_sat(param_set, T[i], ρ[i], q_tot[i])
                T_unsat = TD.air_temperature(param_set, e_int, q_tot[i])
                worst_gap = max(worst_gap, abs(T[i] - T_unsat))
                sol = TD.saturation_adjustment(
                    param_set,
                    TD.ρe(),
                    ρ[i],
                    e_int,
                    q_tot[i],
                )
                worst_err = max(worst_err, abs(sol.T - T[i]))
            end
            @test worst_gap < FT(5)
            @test worst_err < FT(0.01)
        end

        @testset "Saturation excess derivative matches AD ($FT)" begin
            # ∂q_vap_sat_∂T omitted the liquid-fraction term, which made it wrong by up to
            # 7.4% through the mixed-phase band while remaining exact above freezing.
            ρ = FT(0.9)
            rtol = FT === Float32 ? FT(1e-2) : FT(1e-6)
            for T in FT.((235, 245, 255, 265, 271, 274, 290))
                ad = ForwardDiff.derivative(
                    t -> TD.q_vap_saturation(param_set, t, ρ),
                    T,
                )
                @test TD.∂q_vap_sat_∂T(param_set, T, ρ) ≈ ad rtol = rtol
            end
        end

        @testset "Non-integer pow_icenuc is safe below T_icenuc ($FT)" begin
            # A negative base raised to a non-integer power threw a DomainError, which is
            # unrecoverable inside a GPU kernel. Only the default pow_icenuc = 1 avoided it.
            base = FT == Float64 ? param_set_Float64 : param_set_Float32
            fields = fieldnames(typeof(base))
            for n in FT.((0.5, 1.0, 1.5))
                ps = TP.ThermodynamicsParameters{FT}(;
                    (f => getfield(base, f) for f in fields)...,
                    pow_icenuc = n,
                )
                for T in FT.((150, 200, 233, 250, 273, 300))
                    λ = TD.liquid_fraction_ramp(ps, T)
                    @test isfinite(λ)
                    @test 0 <= λ <= 1
                    @test isfinite(TD.∂λ_∂T_ramp(ps, T))
                end
            end
        end

        @testset "Temperature guards ($FT)" begin
            # Saturation vapor pressure used to throw for negative temperatures, which the
            # solver could reach from a bad iterate.
            @test TD.saturation_vapor_pressure(param_set, zero(FT), TD.Liquid()) ≈ zero(FT)
            @test isfinite(TD.saturation_vapor_pressure(param_set, eps(FT), TD.Liquid()))
        end
    end

    @testset "Argument promotion keeps AD type-stable" begin
        # The promoting fallbacks were shadowed by the ::APS methods and never ran, so
        # differentiating with respect to a latent heat left the zero-temperature guard
        # comparing a float against a dual and the return type inferred as Any.
        param_set = param_set_Float64
        D = ForwardDiff.Dual{Nothing, Float64, 1}
        @test only(
            Base.return_types(
                TD.saturation_vapor_pressure_calc,
                (typeof(param_set), Float64, D, Float64),
            ),
        ) === D
        @test only(
            Base.return_types(
                TD.latent_heat_generic,
                (typeof(param_set), Float64, D, Float64),
            ),
        ) === D
    end

    @testset "solution_type reaches the solver" begin
        # solution_type() had no call sites, so DataCollection could never collect and
        # always reported zeros.
        param_set = param_set_Float64
        (; ρ, q_liq, q_ice) = equilibrium_state_ρ(param_set, 290.0, 85000.0, 0.019)
        e_int = TD.internal_energy(param_set, 290.0, 0.019, q_liq, q_ice)

        TD.DataCollection.reset_stats()
        @test TD.DataCollection.get_data().call_counter == 0

        @eval TD solution_type() = $(RS.VerboseSolution)()
        try
            for _ in 1:3
                TD.saturation_adjustment(
                    RS.NewtonsMethod,
                    param_set,
                    TD.ρe(),
                    ρ,
                    e_int,
                    0.019,
                    20,
                    1e-6,
                )
            end
            @test TD.DataCollection.get_data().call_counter == 3
        finally
            @eval TD solution_type() = $(RS.CompactSolution)()
            TD.DataCollection.reset_stats()
        end
    end
end
