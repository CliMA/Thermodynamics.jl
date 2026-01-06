"""
# Correctness Test Suite (functional API)

This file contains tests for fundamental thermodynamic relations and physical laws,
using the non-deprecated functional API (no `PhasePartition`/state types).
"""

@testset "Thermodynamics - correctness (functional)" begin
    for FT in (Float32, Float64)
        param_set = FT == Float64 ? param_set_Float64 : param_set_Float32

        @testset "Ideal gas law consistency ($FT)" begin
            T = FT(300)
            ρ = FT(1.2)
            q_tot = FT(0.01)
            q_liq = FT(0)
            q_ice = FT(0)
            p = TD.air_pressure(param_set, T, ρ, q_tot, q_liq, q_ice)
            R_m = TD.gas_constant_air(param_set, q_tot, q_liq, q_ice)
            @test p ≈ R_m * ρ * T
            @test TD.air_density(param_set, T, p, q_tot, q_liq, q_ice) ≈ ρ
        end

        @testset "Mixture cp/cv are mass-fraction consistent ($FT)" begin
            q_tot = FT(0.02)
            q_liq = FT(0.003)
            q_ice = FT(0.001)
            q_vap = TD.vapor_specific_humidity(q_tot, q_liq, q_ice)

            cp_expected =
                (1 - q_tot) * TP.cp_d(param_set) +
                q_vap * TP.cp_v(param_set) +
                q_liq * TP.cp_l(param_set) +
                q_ice * TP.cp_i(param_set)
            cv_expected =
                (1 - q_tot) * TP.cv_d(param_set) +
                q_vap * TP.cv_v(param_set) +
                q_liq * TP.cv_l(param_set) +
                q_ice * TP.cv_i(param_set)

            @test TD.cp_m(param_set, q_tot, q_liq, q_ice) ≈ cp_expected
            @test TD.cv_m(param_set, q_tot, q_liq, q_ice) ≈ cv_expected
        end

        @testset "Energy/temperature inversions are consistent ($FT)" begin
            T0 = FT(287)
            q_tot = FT(0.015)
            q_liq = FT(0.0)
            q_ice = FT(0.0)
            e_int = TD.internal_energy(param_set, T0, q_tot, q_liq, q_ice)
            T = TD.air_temperature(param_set, TD.ρeq(), e_int, q_tot, q_liq, q_ice)
            @test isapprox(T, T0; atol = FT(atol_temperature), rtol = FT(0))
        end

        @testset "cp_m - cv_m = R_m ($FT)" begin
            # Calorically-perfect ideal-gas mixture identity, with condensed phases contributing no R.
            q_tot = FT(0.02)
            q_liq = FT(0.003)
            q_ice = FT(0.001)
            R_m = TD.gas_constant_air(param_set, q_tot, q_liq, q_ice)
            @test TD.cp_m(param_set, q_tot, q_liq, q_ice) -
                  TD.cv_m(param_set, q_tot, q_liq, q_ice) ≈
                  R_m
        end

        @testset "h = e_int + R_m T ($FT)" begin
            # Ideal-gas mixture identity (condensed-phase specific volume neglected).
            T = FT(295)
            q_tot = FT(0.017)
            q_liq = FT(0.002)
            q_ice = FT(0.001)
            e_int = TD.internal_energy(param_set, T, q_tot, q_liq, q_ice)
            h = TD.enthalpy(param_set, T, q_tot, q_liq, q_ice)
            R_m = TD.gas_constant_air(param_set, q_tot, q_liq, q_ice)
            @test h ≈ e_int + R_m * T
        end

        @testset "Latent heats are consistent at T0 ($FT)" begin
            T_0 = TP.T_0(param_set)
            @test TD.latent_heat_vapor(param_set, T_0) ≈ TP.LH_v0(param_set)
            @test TD.latent_heat_sublim(param_set, T_0) ≈ TP.LH_s0(param_set)
            @test TD.latent_heat_fusion(param_set, T_0) ≈ TP.LH_f0(param_set)
        end

        @testset "Latent heat identity: L_s = L_v + L_f ($FT)" begin
            for T in (FT(240), FT(273.15), FT(300))
                @test TD.latent_heat_sublim(param_set, T) ≈
                      TD.latent_heat_vapor(param_set, T) +
                      TD.latent_heat_fusion(param_set, T)
            end
        end

        @testset "Saturation vapor pressure at triple point ($FT)" begin
            T_tr = TP.T_triple(param_set)
            p_tr = TP.press_triple(param_set)
            @test TD.saturation_vapor_pressure(param_set, T_tr, TD.Liquid()) ≈ p_tr
            @test TD.saturation_vapor_pressure(param_set, T_tr, TD.Ice()) ≈ p_tr
        end

        @testset "Clausius-Clapeyron (d ln e_s / dT) ($FT)" begin
            R_v = TP.R_v(param_set)
            T = FT(280)
            δ = FT(0.1)
            e_sat = T_ -> TD.saturation_vapor_pressure(param_set, T_, TD.Liquid())
            dlog_es_dT_fd = (log(e_sat(T + δ)) - log(e_sat(T - δ))) / (2δ)
            L = TD.latent_heat_vapor(param_set, T)
            dlog_es_dT_cc = L / (R_v * T^2)
            @test isapprox(dlog_es_dT_fd, dlog_es_dT_cc; rtol = FT(rtol_temperature_fd))
        end

        @testset "Entropy reduces to dry-air form in the dry limit ($FT)" begin
            T1 = FT(280)
            T2 = FT(310)
            p1 = FT(9e4)
            p2 = FT(8e4)

            # Full entropy reduces to dry-air entropy when q_tot=q_liq=q_ice=0
            @test TD.entropy(param_set, p1, T1) ≈ TD.entropy_dry(param_set, p1, T1)

            # Difference form matches the ideal-gas entropy change for dry air
            Δs = TD.entropy_dry(param_set, p2, T2) - TD.entropy_dry(param_set, p1, T1)
            Δs_expected =
                TP.cp_d(param_set) * log(T2 / T1) - TP.R_d(param_set) * log(p2 / p1)
            @test isapprox(Δs, Δs_expected; rtol = FT(5e-6), atol = FT(5e-6))
        end

        @testset "Saturation equilibrium partitioning constraints ($FT)" begin
            T = FT(270)
            ρ = FT(1.1)

            q_sat = TD.q_vap_saturation(param_set, T, ρ)

            # Unsaturated: no condensate
            q_tot = FT(0.8) * q_sat
            (q_liq, q_ice) = TD.condensate_partition(param_set, T, ρ, q_tot)
            @test q_liq == 0
            @test q_ice == 0

            # Saturated: vapor is at saturation and condensate is nonnegative
            q_tot = FT(1.3) * q_sat
            (q_liq, q_ice) = TD.condensate_partition(param_set, T, ρ, q_tot)
            @test q_liq ≥ 0
            @test q_ice ≥ 0
            @test q_liq + q_ice ≤ q_tot + FT(100) * eps(FT)

            q_vap = TD.vapor_specific_humidity(q_tot, q_liq, q_ice)
            q_vap_sat = TD.q_vap_saturation(param_set, T, ρ, q_liq, q_ice)
            @test isapprox(q_vap, q_vap_sat; rtol = FT(1e-6), atol = FT(50) * eps(FT))
        end
    end
end
