"""
# Default Behavior Accuracy Test Suite

This file contains tests for saturation adjustment accuracy and convergence.
"""

using Test
import ForwardDiff
using Thermodynamics
import Thermodynamics as TD
import Thermodynamics.Parameters as TP
using Thermodynamics.TestedProfiles

# Include common parameters
include("common_parameters.jl")

@testset "Thermodynamics - default behavior accuracy" begin
    # Input arguments should be accurate within machine precision
    # Temperature is approximated via saturation adjustment, and should be within a physical tolerance

    or(a, b) = a || b
    for ArrayType in array_types
        FT = eltype(ArrayType)
        param_set = TP.ThermodynamicsParameters(FT)
        _cp_d = TP.cp_d(param_set)
        _eint_v0 = TP.e_int_v0(param_set)

        profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
        (; T, p, e_int, ρ, θ_liq_ice, phase_type) = profiles
        (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = profiles

        RH_sat_mask = or.(RH .> 1, RH .≈ 1)
        RH_unsat_mask = .!or.(RH .> 1, RH .≈ 1)
        ts = PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot)
        @test all(saturated.(param_set, ts[RH_sat_mask]))
        @test !any(saturated.(param_set, ts[RH_unsat_mask]))

        # Test Clausius Clapeyron relation
        k = findfirst(q -> q > 0.01, q_tot) # test for one value with q_tot above some threshhold
        ts_sol = TD.PhaseEquil_ρTq(param_set, ρ[k], T[k], q_tot[k])

        function q_vap_sat(_T::FT) where {FT}
            _ρ = TD.air_density(param_set, ts_sol)
            _q_tot = TD.total_specific_humidity(param_set, ts_sol)
            _phase_type = PhaseEquil{FT}
            _q_pt = PhasePartition_equil(
                param_set,
                _T,
                oftype(_T, _ρ),
                oftype(_T, _q_tot),
                _phase_type,
            )
            return TD.q_vap_saturation(
                param_set,
                _T,
                oftype(_T, _ρ),
                _phase_type,
                _q_pt,
            )
        end

        function ∂q_vap_sat_∂T_vs_T(_T::FT) where {FT}
            _ρ = TD.air_density(param_set, ts_sol)
            _q_tot = TD.total_specific_humidity(param_set, ts_sol)
            _phase_type = PhaseEquil{FT}
            _λ = TD.liquid_fraction(param_set, ts_sol)
            _q_pt = TD.PhasePartition_equil(
                param_set,
                _T,
                oftype(_T, _ρ),
                oftype(_T, _q_tot),
                _phase_type,
            )
            _q_vap_sat =
                TD.q_vap_saturation(param_set, _T, _ρ, _phase_type, _q_pt)
            return TD.∂q_vap_sat_∂T(
                param_set,
                oftype(_T, _λ),
                _T,
                oftype(_T, _q_vap_sat),
            )
        end

        ∂q_vap_sat_∂T_fd = _T -> ForwardDiff.derivative(q_vap_sat, _T)
        @test all(
            isapprox.(
                log.(∂q_vap_sat_∂T_fd.(T)),
                log.(∂q_vap_sat_∂T_vs_T.(T));
                rtol = 2e-2,
            ),
        )

        # PhaseEquil (freezing)
        _T_freeze = FT(TP.T_freeze(param_set))
        e_int_upper =
            internal_energy_sat.(
                param_set,
                Ref(_T_freeze + sqrt(eps(FT))),
                ρ,
                q_tot,
                phase_type,
            )
        e_int_lower =
            internal_energy_sat.(
                param_set,
                Ref(_T_freeze - sqrt(eps(FT))),
                ρ,
                q_tot,
                phase_type,
            )
        _e_int = (e_int_upper .+ e_int_lower) / 2
        ts = PhaseEquil_ρeq.(param_set, ρ, _e_int, q_tot)
        @test all(
            isapprox.(
                air_temperature.(param_set, ts),
                Ref(_T_freeze),
                atol = atol_temperature,
            ),
        )

        # Args needs to be in sync with PhaseEquil:
        ts =
            PhaseEquil_ρeq.(
                param_set,
                ρ,
                _e_int,
                q_tot,
                8,
                FT(rtol_temperature),
                RS.SecantMethod,
            )
        @test all(
            isapprox.(
                air_temperature.(param_set, ts),
                Ref(_T_freeze),
                atol = atol_temperature,
            ),
        )

        # PhaseEquil
        ts_exact = PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot, 100, FT(1e-6))
        ts = PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot)
        @test all(
            isapprox.(
                T,
                air_temperature.(param_set, ts),
                atol = atol_temperature,
            ),
        )

        # Should be machine accurate (because ts contains `e_int`,`ρ`,`q_tot`):
        @test all(compare_moisture.(param_set, ts, ts_exact))
        @test all(
            internal_energy.(param_set, ts) .≈
            internal_energy.(param_set, ts_exact),
        )
        @test all(
            air_density.(param_set, ts) .≈ air_density.(param_set, ts_exact),
        )
        # Approximate (temperature must be computed via saturation adjustment):
        @test all(
            isapprox.(
                air_pressure.(param_set, ts),
                air_pressure.(param_set, ts_exact),
                rtol = rtol_pressure,
            ),
        )
        @test all(
            isapprox.(
                air_temperature.(param_set, ts),
                air_temperature.(param_set, ts_exact),
                atol = atol_temperature,
            ),
        )

        dry_mask = abs.(q_tot .- 0) .< eps(FT)
        q_dry = q_pt[dry_mask]
        @test all(
            condensate.(q_pt) .==
            getproperty.(q_pt, :liq) .+ getproperty.(q_pt, :ice),
        )
        @test all(has_condensate.(q_dry) .== false)

        e_tot = total_energy.(param_set, ts, e_kin, e_pot)
        _cp_d = FT(TP.cp_d(param_set))
        @test all(
            specific_enthalpy.(param_set, ts) .≈
            e_int .+
            gas_constant_air.(param_set, ts) .* air_temperature.(param_set, ts),
        )
        @test all(
            specific_enthalpy.(param_set, ts) .≈
            e_int .+
            air_pressure.(param_set, ts) ./ air_density.(param_set, ts),
        )
        @test all(
            total_specific_enthalpy.(param_set, ts, e_tot) .≈
            specific_enthalpy.(param_set, ts) .+ e_kin .+ e_pot,
        )
        @test all(
            moist_static_energy.(param_set, ts, e_pot) .≈
            specific_enthalpy.(param_set, ts) .+ e_pot,
        )
        @test all(
            moist_static_energy.(param_set, ts, e_pot) .≈
            total_specific_enthalpy.(param_set, ts, e_tot) .- e_kin,
        )
        @test all(
            virtual_dry_static_energy.(param_set, ts, e_pot) .≈
            _cp_d .* virtual_temperature.(param_set, ts) .+ e_pot,
        )

        # PhaseEquil
        ts_exact =
            PhaseEquil_ρeq.(
                param_set,
                ρ,
                e_int,
                q_tot,
                100,
                FT(1e-6),
                RS.SecantMethod,
            )
        ts =
            PhaseEquil_ρeq.(
                param_set,
                ρ,
                e_int,
                q_tot,
                35,
                FT(rtol_temperature),
                RS.SecantMethod,
            ) # Needs to be in sync with default
        # Should be machine accurate (because ts contains `e_int`,`ρ`,`q_tot`):
        @test all(compare_moisture.(param_set, ts, ts_exact))
        @test all(
            internal_energy.(param_set, ts) .≈
            internal_energy.(param_set, ts_exact),
        )
        @test all(
            air_density.(param_set, ts) .≈ air_density.(param_set, ts_exact),
        )
        # Approximate (temperature must be computed via saturation adjustment):
        @test all(
            isapprox.(
                air_pressure.(param_set, ts),
                air_pressure.(param_set, ts_exact),
                rtol = rtol_pressure,
            ),
        )
        @test all(
            isapprox.(
                air_temperature.(param_set, ts),
                air_temperature.(param_set, ts_exact),
                atol = atol_temperature,
            ),
        )

        # PhaseEquil_ρθq
        ts_exact = PhaseEquil_ρθq.(param_set, ρ, θ_liq_ice, q_tot, 45, FT(1e-6))
        ts = PhaseEquil_ρθq.(param_set, ρ, θ_liq_ice, q_tot)
        # Should be machine accurate:
        @test all(
            air_density.(param_set, ts) .≈ air_density.(param_set, ts_exact),
        )
        @test all(compare_moisture.(param_set, ts, ts_exact))
        # Approximate (temperature must be computed via saturation adjustment):
        @test all(
            abs.(
                internal_energy.(param_set, ts) .-
                internal_energy.(param_set, ts_exact)
            ) .<=
            (atol_energy_temperature .+ rtol_humidity * _eint_v0 .* q_tot),
        )
        @test all(
            isapprox.(
                liquid_ice_pottemp.(param_set, ts),
                liquid_ice_pottemp.(param_set, ts_exact),
                atol = atol_temperature,
            ),
        )
        @test all(
            isapprox.(
                air_temperature.(param_set, ts),
                air_temperature.(param_set, ts_exact),
                atol = atol_temperature,
            ),
        )

        # PhaseEquil_pθq
        ts_exact = PhaseEquil_pθq.(param_set, p, θ_liq_ice, q_tot, 40, FT(1e-6))
        ts = PhaseEquil_pθq.(param_set, p, θ_liq_ice, q_tot)

        ts =
            PhaseEquil_pθq.(
                param_set,
                p,
                θ_liq_ice,
                q_tot,
                40,
                FT(rtol_temperature),
                RS.RegulaFalsiMethod,
            )
        # Should be machine accurate:
        @test all(compare_moisture.(param_set, ts, ts_exact))
        # Approximate (temperature must be computed via saturation adjustment):
        @test all(
            isapprox.(
                air_density.(param_set, ts),
                air_density.(param_set, ts_exact),
                rtol = rtol_density,
            ),
        )
        @test all(
            abs.(
                internal_energy.(param_set, ts) .-
                internal_energy.(param_set, ts_exact)
            ) .<=
            (atol_energy_temperature .+ rtol_humidity * _eint_v0 .* q_tot),
        )
        @test all(
            isapprox.(
                liquid_ice_pottemp.(param_set, ts),
                liquid_ice_pottemp.(param_set, ts_exact),
                atol = atol_temperature,
            ),
        )
        @test all(
            isapprox.(
                air_temperature.(param_set, ts),
                air_temperature.(param_set, ts_exact),
                atol = atol_temperature,
            ),
        )

        # PhaseEquil_pθq (freezing)
        _T_freeze = FT(TP.T_freeze(param_set))

        function θ_liq_ice_closure(T, p, q_tot)
            q_pt_closure = TD.PhasePartition_equil_given_p(
                param_set,
                T,
                p,
                q_tot,
                phase_type,
            )
            θ_liq_ice_close = TD.liquid_ice_pottemp_given_pressure(
                param_set,
                T,
                p,
                q_pt_closure,
            )
            return θ_liq_ice_close
        end

        θ_liq_ice_upper =
            θ_liq_ice_closure.(Ref(_T_freeze + sqrt(eps(FT))), p, q_tot)
        θ_liq_ice_lower =
            θ_liq_ice_closure.(Ref(_T_freeze - sqrt(eps(FT))), p, q_tot)
        θ_liq_ice_mid = (θ_liq_ice_upper .+ θ_liq_ice_lower) ./ 2

        ts_lower = PhaseEquil_pθq.(param_set, p, θ_liq_ice_lower, q_tot)
        ts_upper = PhaseEquil_pθq.(param_set, p, θ_liq_ice_upper, q_tot)
        ts_mid = PhaseEquil_pθq.(param_set, p, θ_liq_ice_mid, q_tot)

        @test count(air_temperature.(param_set, ts_lower) .== Ref(_T_freeze)) ≥
              217
        @test count(air_temperature.(param_set, ts_upper) .== Ref(_T_freeze)) ≥
              217
        @test count(air_temperature.(param_set, ts_mid) .== Ref(_T_freeze)) ≥
              1296
        # we should do this instead, but we're failing because some inputs are bad
        # E.g. p ~ 110_000 Pa, q_tot ~ 0.16, which results in negative θ_liq_ice
        # This means that we should probably update our tested profiles.
        # @test all(air_temperature.(param_set, ts_lower) .== Ref(_T_freeze))
        # @test all(air_temperature.(param_set, ts_upper) .== Ref(_T_freeze))
        # @test all(air_temperature.(param_set, ts_mid) .== Ref(_T_freeze))

        # @show ρ, θ_liq_ice, q_pt
        # PhaseNonEquil_ρθq
        ts_exact =
            PhaseNonEquil_ρθq.(param_set, ρ, θ_liq_ice, q_pt, 40, FT(1e-6))
        ts = PhaseNonEquil_ρθq.(param_set, ρ, θ_liq_ice, q_pt)
        # Should be machine accurate:
        @test all(compare_moisture.(param_set, ts, ts_exact))
        @test all(
            air_density.(param_set, ts) .≈ air_density.(param_set, ts_exact),
        )
        # Approximate (temperature must be computed via non-linear solve):
        @test all(
            abs.(
                internal_energy.(param_set, ts) .-
                internal_energy.(param_set, ts_exact)
            ) .<= (
                atol_energy_temperature .+
                rtol_humidity * _eint_v0 .* getproperty.(q_pt, :tot)
            ),
        )
        # Potential temperature comparison - use a more appropriate tolerance
        pot_temp_errors =
            abs.(
                liquid_ice_pottemp.(param_set, ts) .-
                liquid_ice_pottemp.(param_set, ts_exact)
            )
        # For potential temperature, use a tolerance based on the iterative nature of the calculation
        pot_temp_tolerances = 2.25 * atol_temperature  # Larger tol to allow for condensate effects
        failing_pot_temp = findall(pot_temp_errors .> pot_temp_tolerances)
        # if !isempty(failing_pot_temp)
        #     @info "Potential temperature test failures" n=length(failing_pot_temp)
        #     @show maximum(pot_temp_errors)
        #     @show maximum(pot_temp_tolerances)
        # end
        @test all(pot_temp_errors .<= pot_temp_tolerances)
        @test all(
            isapprox.(
                air_temperature.(param_set, ts),
                air_temperature.(param_set, ts_exact),
                atol = atol_temperature,
            ),
        )

    end

end 