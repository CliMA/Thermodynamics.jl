"""
# Physical property tests

Two kinds of test live here:

1. **Direct coverage** of exported functions that the rest of the suite only reaches
   indirectly, or not at all. Several of these have no caller anywhere in `src/`, so
   without a direct test an error in them would go unnoticed.

2. **Property tests**: thermodynamic identities, orderings, and monotonicity that must hold
   regardless of the parameter values. These constrain the physics rather than pinning
   numbers computed from the same code, so they stay meaningful as parameters change.
"""

@testset "Thermodynamics - physical properties" begin
    for FT in (Float32, Float64)
        param_set = FT == Float64 ? param_set_Float64 : param_set_Float32

        T_freeze = TP.T_freeze(param_set)
        T_triple = TP.T_triple(param_set)
        T_icenuc = TP.T_icenuc(param_set)
        R_d = TP.R_d(param_set)
        R_v = TP.R_v(param_set)
        cp_d = TP.cp_d(param_set)
        cv_d = TP.cv_d(param_set)

        @testset "Speed of sound ($FT)" begin
            # Dry limit: c = sqrt(γ R_d T) with γ = cp_d/cv_d
            T = FT(300)
            @test TD.soundspeed_air(param_set, T) ≈ sqrt(cp_d / cv_d * R_d * T)
            # ~347 m/s in dry air at 300 K
            @test FT(340) < TD.soundspeed_air(param_set, T) < FT(355)
            # Moist air is faster than dry air at the same temperature
            @test TD.soundspeed_air(param_set, T, FT(0.02)) >
                  TD.soundspeed_air(param_set, T)
            # Monotone in temperature
            @test TD.soundspeed_air(param_set, FT(310)) >
                  TD.soundspeed_air(param_set, FT(290))
        end

        @testset "Supersaturation ($FT)" begin
            T = FT(280)
            ρ = FT(1)
            q_sat = TD.q_vap_saturation(param_set, T, ρ, TD.Liquid())
            # Exactly saturated air has zero supersaturation
            @test TD.supersaturation(param_set, q_sat, ρ, T, TD.Liquid()) ≈ FT(0) atol =
                sqrt(eps(FT))
            # 10% excess vapor gives S = 0.1
            @test TD.supersaturation(param_set, FT(1.1) * q_sat, ρ, T, TD.Liquid()) ≈
                  FT(0.1) rtol = FT(1e-4)
            # Subsaturated air is negative, and bounded below by -1 as q_vap → 0
            @test TD.supersaturation(param_set, FT(0.5) * q_sat, ρ, T, TD.Liquid()) < 0
            @test TD.supersaturation(param_set, zero(FT), ρ, T, TD.Liquid()) ≈ -FT(1)
            # Below freezing, air saturated over liquid is supersaturated over ice
            T_cold = FT(250)
            q_sat_liq = TD.q_vap_saturation(param_set, T_cold, ρ, TD.Liquid())
            @test TD.supersaturation(param_set, q_sat_liq, ρ, T_cold, TD.Ice()) > 0
            # The five-argument form agrees when handed the same saturation vapor pressure
            p_v_sat = TD.saturation_vapor_pressure(param_set, T, TD.Liquid())
            @test TD.supersaturation(param_set, q_sat, ρ, T, p_v_sat) ≈
                  TD.supersaturation(param_set, q_sat, ρ, T, TD.Liquid())
        end

        @testset "Volumetric vapor mixing ratio ($FT)" begin
            q_tot = FT(0.02)
            # No condensate: r_vol = (R_v/R_d) * q_vap / (1 - q_tot)
            expected = TP.Rv_over_Rd(param_set) * q_tot / (1 - q_tot)
            @test TD.vol_vapor_mixing_ratio(param_set, q_tot) ≈ expected
            # Dry air has no vapor
            @test TD.vol_vapor_mixing_ratio(param_set, zero(FT)) ≈ FT(0)
            # Condensate does not count toward the vapor mixing ratio
            @test TD.vol_vapor_mixing_ratio(param_set, q_tot, FT(0.005)) <
                  TD.vol_vapor_mixing_ratio(param_set, q_tot)
        end

        @testset "Exner function and potential temperatures ($FT)" begin
            T = FT(280)
            ρ = FT(1)
            p = TD.air_pressure(param_set, T, ρ)
            p0 = TP.p_ref_theta(param_set)

            # Π = (p/p₀)^(R_m/cp_m), and the two forms must agree
            @test TD.exner(param_set, T, ρ) ≈ (p / p0)^(R_d / cp_d)
            @test TD.exner(param_set, T, ρ) ≈
                  TD.exner_given_pressure(param_set, p)
            # θ = T/Π, so θ = T at the reference pressure
            @test TD.potential_temperature(param_set, T, ρ) ≈ T / TD.exner(param_set, T, ρ)
            ρ0 = TD.air_density(param_set, T, p0)
            @test TD.potential_temperature(param_set, T, ρ0) ≈ T
            # Below the reference pressure, θ exceeds T
            @test TD.potential_temperature(param_set, T, ρ / 2) > T

            # θ_v = θ R_m/R_d, so moist air has the larger virtual potential temperature,
            # and the two coincide for dry air
            @test TD.virtual_pottemp(param_set, T, ρ) ≈
                  TD.potential_temperature(param_set, T, ρ)
            @test TD.virtual_pottemp(param_set, T, ρ, FT(0.02)) >
                  TD.potential_temperature(param_set, T, ρ, FT(0.02))
        end

        @testset "Saturation vapor pressure ordering ($FT)" begin
            # Liquid and ice coincide at the triple point
            @test TD.saturation_vapor_pressure(param_set, T_triple, TD.Liquid()) ≈
                  TD.saturation_vapor_pressure(param_set, T_triple, TD.Ice())
            # Below it, ice has the lower saturation vapor pressure. This is what drives the
            # Wegener-Bergeron-Findeisen process, and a swapped Liquid()/Ice() dispatch
            # would reverse it.
            for T in FT.((200, 220, 240, 260, 270))
                @test TD.saturation_vapor_pressure(param_set, T, TD.Ice()) <
                      TD.saturation_vapor_pressure(param_set, T, TD.Liquid())
            end
            # Both are strictly increasing in temperature
            for T in FT.((200, 240, 273, 300))
                for phase in (TD.Liquid(), TD.Ice())
                    @test TD.saturation_vapor_pressure(param_set, T + 1, phase) >
                          TD.saturation_vapor_pressure(param_set, T, phase)
                end
            end
            # The mixture is bracketed by the pure phases, and reduces to them at λ = 1, 0
            T = FT(250)
            p_liq = TD.saturation_vapor_pressure(param_set, T, TD.Liquid())
            p_ice = TD.saturation_vapor_pressure(param_set, T, TD.Ice())
            @test TD.saturation_vapor_pressure_mixture(param_set, T, one(FT)) ≈ p_liq
            @test TD.saturation_vapor_pressure_mixture(param_set, T, zero(FT)) ≈ p_ice
            for λ in FT.((0.25, 0.5, 0.75))
                p_mix = TD.saturation_vapor_pressure_mixture(param_set, T, λ)
                @test p_ice < p_mix < p_liq
            end
        end

        @testset "Clausius-Clapeyron ($FT)" begin
            # d(ln p_v^*)/dT = L/(R_v T²) for each pure phase, checked by central differences
            # at several temperatures rather than one, so that a compensating constant error
            # cannot hide.
            δ = FT(0.05)
            for T in FT.((240, 260, 280, 300))
                for (phase, L) in (
                    (TD.Liquid(), TD.latent_heat_vapor(param_set, T)),
                    (TD.Ice(), TD.latent_heat_sublim(param_set, T)),
                )
                    p_plus = TD.saturation_vapor_pressure(param_set, T + δ, phase)
                    p_minus = TD.saturation_vapor_pressure(param_set, T - δ, phase)
                    p_mid = TD.saturation_vapor_pressure(param_set, T, phase)
                    dlnp_dT = (log(p_plus) - log(p_minus)) / (2δ)
                    @test dlnp_dT ≈ L / (R_v * T^2) rtol = FT(1e-3)
                end
            end

            # Across the mixed-phase band the effective latent heat is the liquid-fraction
            # weighted mean, and the saturation curve additionally migrates between the two
            # pure curves as λ changes with temperature. Compare the analytic derivative
            # against differences of q_vap_saturation itself, which includes both effects.
            ρ = FT(0.9)
            for T in FT.((235, 245, 255, 265, 272))
                q_plus = TD.q_vap_saturation(param_set, T + δ, ρ)
                q_minus = TD.q_vap_saturation(param_set, T - δ, ρ)
                fd = (q_plus - q_minus) / (2δ)
                @test TD.∂q_vap_sat_∂T(param_set, T, ρ) ≈ fd rtol = FT(1e-2)
            end
        end

        @testset "Monotonicity and bounds ($FT)" begin
            ρ = FT(1)
            # Saturation specific humidity increases with temperature, decreases with density
            for T in FT.((240, 270, 300))
                @test TD.q_vap_saturation(param_set, T + 1, ρ) >
                      TD.q_vap_saturation(param_set, T, ρ)
                @test TD.q_vap_saturation(param_set, T, 2ρ) <
                      TD.q_vap_saturation(param_set, T, ρ)
            end
            # Liquid fraction is a fraction, and both parameterizations respect that
            for T in FT.((200, 230, 250, 265, 273, 280, 320))
                @test 0 <= TD.liquid_fraction_ramp(param_set, T) <= 1
                @test 0 <= TD.liquid_fraction(param_set, T, zero(FT), zero(FT)) <= 1
                @test 0 <= TD.liquid_fraction(param_set, T, FT(1e-3), FT(1e-3)) <= 1
            end
            # Both are non-decreasing in temperature
            for T in FT.((200, 240, 260, 270, 274))
                @test TD.liquid_fraction_ramp(param_set, T + 1) >=
                      TD.liquid_fraction_ramp(param_set, T)
            end
            # Relative humidity stays within [0, 1] and reaches 1 at saturation
            T = FT(285)
            p = FT(90000)
            q_sat = TD.q_vap_saturation_from_pressure(param_set, FT(0), p, T)
            @test TD.relative_humidity(param_set, T, p, q_sat) ≈ FT(1) rtol = FT(1e-3)
            for q in FT.((0, 0.001, 0.005, 0.01, 0.05))
                RH = TD.relative_humidity(param_set, T, p, q)
                @test 0 <= RH <= 1
            end
        end

        @testset "Liquid fraction: condensate-free ramp ($FT)" begin
            # The four-argument form falls back to a 0.2 K ramp ending at T_freeze when no
            # condensate is present. This branch is otherwise only reached indirectly.
            @test TD.liquid_fraction(param_set, T_freeze, zero(FT), zero(FT)) ≈ FT(1)
            @test TD.liquid_fraction(param_set, T_freeze - FT(0.2), zero(FT), zero(FT)) ≈
                  FT(0)
            @test TD.liquid_fraction(param_set, T_freeze - FT(0.1), zero(FT), zero(FT)) ≈
                  FT(0.5) rtol = FT(1e-3)
            @test TD.liquid_fraction(param_set, T_freeze + 10, zero(FT), zero(FT)) ≈ FT(1)
            @test TD.liquid_fraction(param_set, T_freeze - 10, zero(FT), zero(FT)) ≈ FT(0)
            # With condensate present the fraction is set by the condensate itself
            @test TD.liquid_fraction(param_set, FT(250), FT(3e-3), FT(1e-3)) ≈ FT(0.75)
            # The ramp parameterization is pinned at its endpoints
            @test TD.liquid_fraction_ramp(param_set, T_freeze) ≈ FT(1)
            @test TD.liquid_fraction_ramp(param_set, T_icenuc) ≈ FT(0)
        end

        @testset "Dry limit ($FT)" begin
            # With no water at all, the mixture properties collapse to the dry-air values
            z = zero(FT)
            T = FT(290)
            @test TD.cp_m(param_set, z, z, z) ≈ cp_d
            @test TD.cv_m(param_set, z, z, z) ≈ cv_d
            @test TD.gas_constant_air(param_set, z, z, z) ≈ R_d
            # Dry internal energy carries a -R_d T_0 offset relative to cv_d (T - T_0), so
            # that enthalpy comes out as cp_d (T - T_0) with no offset.
            T_0 = TP.T_0(param_set)
            @test TD.internal_energy(param_set, T, z, z, z) ≈
                  cv_d * (T - T_0) - R_d * T_0
            @test TD.internal_energy(param_set, T, z, z, z) ≈
                  TD.internal_energy_dry(param_set, T)
            @test TD.enthalpy(param_set, T, z, z, z) ≈ cp_d * (T - T_0)
            # h = e_int + R_m T in the dry limit
            @test TD.enthalpy(param_set, T, z, z, z) ≈
                  TD.internal_energy(param_set, T, z, z, z) + R_d * T
            @test TD.virtual_temperature(param_set, T, z, z, z) ≈ T
            @test TD.vapor_specific_humidity(z, z, z) ≈ z
            # Saturation adjustment of dry air returns the dry temperature and no condensate
            ρ = FT(1)
            e_int = TD.internal_energy(param_set, T, z, z, z)
            sol = TD.saturation_adjustment(param_set, TD.ρe(), ρ, e_int, z)
            @test sol.T ≈ T rtol = FT(1e-5)
            @test sol.q_liq ≈ z
            @test sol.q_ice ≈ z
        end

        @testset "Round trips ($FT)" begin
            # Build a genuinely cloudy equilibrium state
            T = FT(285)
            p = FT(90000)
            q_tot = FT(0.015)
            (q_liq, q_ice) = TD._condensate_partition_from_p(param_set, T, p, q_tot)
            ρ = TD.air_density(param_set, T, p, q_tot, q_liq, q_ice)
            rtol = FT === Float32 ? FT(1e-3) : FT(1e-6)

            # enthalpy ↔ air_temperature(::ph)
            h = TD.enthalpy(param_set, T, q_tot, q_liq, q_ice)
            @test TD.air_temperature(param_set, TD.ph(), h, q_tot, q_liq, q_ice) ≈ T rtol =
                rtol
            # internal energy ↔ air_temperature(::ρe)
            e_int = TD.internal_energy(param_set, T, q_tot, q_liq, q_ice)
            @test TD.air_temperature(param_set, TD.ρe(), e_int, q_tot, q_liq, q_ice) ≈ T rtol =
                rtol
            # pressure/density ↔ air_temperature(::pρ)
            @test TD.air_temperature(param_set, TD.pρ(), p, ρ, q_tot, q_liq, q_ice) ≈ T rtol =
                rtol
            # θ_li at fixed pressure is an exact inversion
            θ_p = TD.liquid_ice_pottemp_given_pressure(
                param_set,
                T,
                p,
                q_tot,
                q_liq,
                q_ice,
            )
            @test TD.air_temperature(
                param_set,
                TD.pθ_li(),
                p,
                θ_p,
                q_tot,
                q_liq,
                q_ice,
            ) ≈ T rtol = rtol

            # air_temperature(::ρθ_li) is a second-order Taylor approximation rather than an
            # exact inverse, so quantify its truncation error instead of demanding equality.
            θ_ρ = TD.liquid_ice_pottemp(param_set, T, ρ, q_tot, q_liq, q_ice)
            T_taylor =
                TD.air_temperature(param_set, TD.ρθ_li(), ρ, θ_ρ, q_tot, q_liq, q_ice)
            @test abs(T_taylor - T) < FT(1)

            # Vapor pressure ↔ specific humidity
            p_vap = TD.partial_pressure_vapor(param_set, p, q_tot, q_liq, q_ice)
            q_vap = TD.vapor_specific_humidity(q_tot, q_liq, q_ice)
            @test TD.q_vap_from_p_vap(param_set, T, ρ, p_vap) ≈ q_vap rtol = rtol

            # Specific humidity ↔ mixing ratio
            r = TD.specific_humidity_to_mixing_ratio(q_vap, q_tot)
            @test r * (1 - q_tot) ≈ q_vap rtol = rtol
        end
    end
end
