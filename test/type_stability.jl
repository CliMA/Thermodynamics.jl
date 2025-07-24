"""
# Type Stability Test Suite

This file contains tests for performance across different floating-point types.
"""

using Test
using Thermodynamics
import Thermodynamics as TD
import Thermodynamics.Parameters as TP
using Thermodynamics.TestedProfiles

# Saturation adjustment tolerance (relative change of temperature between consecutive iterations)
rtol_temperature = 1e-4

@testset "Thermodynamics - Type Stability" begin

    # NOTE: `Float32` saturation adjustment tends to have more difficulty
    # with converging to the same tolerances as `Float64`, so they're relaxed here.
    ArrayType = Array{Float32}
    FT = eltype(ArrayType)
    param_set = TP.ThermodynamicsParameters(FT)

    profiles = TestedProfiles.PhaseDryProfiles(param_set, ArrayType)
    (; T, p, e_int, ρ, θ_liq_ice, phase_type) = profiles
    (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = profiles
    h = e_int + p ./ ρ

    θ_dry = dry_pottemp.(param_set, T, ρ)
    ts_dry = PhaseDry.(param_set, e_int, ρ)
    ts_dry_ρp = PhaseDry_ρp.(param_set, ρ, p)
    ts_dry_pT = PhaseDry_pT.(param_set, p, T)
    ts_dry_ρθ = PhaseDry_ρθ.(param_set, ρ, θ_dry)
    ts_dry_pθ = PhaseDry_pθ.(param_set, p, θ_dry)
    ts_dry_ph = PhaseDry_ph.(param_set, p, h)

    profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
    (; T, p, e_int, ρ, θ_liq_ice, phase_type) = profiles
    (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = profiles

    ts_eq =
        PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot, 15, FT(rtol_temperature))
    e_tot = total_energy.(param_set, ts_eq, e_kin, e_pot)

    ts_T =
        PhaseEquil_ρTq.(
            param_set,
            air_density.(param_set, ts_eq),
            air_temperature.(param_set, ts_eq),
            q_tot,
        )
    ts_Tp =
        PhaseEquil_pTq.(
            param_set,
            air_pressure.(param_set, ts_eq),
            air_temperature.(param_set, ts_eq),
            q_tot,
        )
    ts_ρp =
        PhaseEquil_ρpq.(
            param_set,
            air_density.(param_set, ts_eq),
            air_pressure.(param_set, ts_eq),
            q_tot,
        )

    @test all(
        air_temperature.(param_set, ts_T) .≈ air_temperature.(param_set, ts_Tp),
    )
    @test all(air_pressure.(param_set, ts_T) .≈ air_pressure.(param_set, ts_Tp))
    @test all(
        total_specific_humidity.(param_set, ts_T) .≈
        total_specific_humidity.(param_set, ts_Tp),
    )

    ts_neq = PhaseNonEquil.(param_set, e_int, ρ, q_pt)
    ts_ρT_neq = PhaseNonEquil_ρTq.(param_set, ρ, T, q_pt)
    ts_pT_neq = PhaseNonEquil_pTq.(param_set, p, T, q_pt)

    ts_θ_liq_ice_eq =
        PhaseEquil_ρθq.(
            param_set,
            ρ,
            θ_liq_ice,
            q_tot,
            40,
            FT(rtol_temperature),
        )
    ts_θ_liq_ice_eq_p =
        PhaseEquil_pθq.(
            param_set,
            p,
            θ_liq_ice,
            q_tot,
            40,
            FT(rtol_temperature),
        )
    ts_θ_liq_ice_neq = PhaseNonEquil_ρθq.(param_set, ρ, θ_liq_ice, q_pt)
    ts_θ_liq_ice_neq_p = PhaseNonEquil_pθq.(param_set, p, θ_liq_ice, q_pt)

    for ts in (
        ts_dry,
        ts_dry_ρp,
        ts_dry_pT,
        ts_dry_ρθ,
        ts_dry_pθ,
        ts_dry_ph,
        ts_eq,
        ts_T,
        ts_Tp,
        ts_ρp,
        ts_neq,
        ts_ρT_neq,
        ts_pT_neq,
        ts_θ_liq_ice_eq,
        ts_θ_liq_ice_eq_p,
        ts_θ_liq_ice_neq,
        ts_θ_liq_ice_neq_p,
    )
        @test typeof.(soundspeed_air.(param_set, ts)) == typeof.(e_int)
        @test typeof.(gas_constant_air.(param_set, ts)) == typeof.(e_int)
        @test typeof.(specific_enthalpy.(param_set, ts)) == typeof.(e_int)
        @test typeof.(vapor_specific_humidity.(param_set, ts)) == typeof.(e_int)
        @test typeof.(relative_humidity.(param_set, ts)) == typeof.(e_int)
        @test typeof.(air_pressure.(param_set, ts)) == typeof.(e_int)
        @test typeof.(air_density.(param_set, ts)) == typeof.(e_int)
        @test typeof.(total_specific_humidity.(param_set, ts)) == typeof.(e_int)
        @test typeof.(liquid_specific_humidity.(param_set, ts)) ==
              typeof.(e_int)
        @test typeof.(ice_specific_humidity.(param_set, ts)) == typeof.(e_int)
        @test typeof.(cp_m.(param_set, ts)) == typeof.(e_int)
        @test typeof.(cv_m.(param_set, ts)) == typeof.(e_int)
        @test typeof.(air_temperature.(param_set, ts)) == typeof.(e_int)
        @test typeof.(internal_energy_sat.(param_set, ts)) == typeof.(e_int)
        @test typeof.(internal_energy.(param_set, ts)) == typeof.(e_int)
        @test typeof.(internal_energy_dry.(param_set, ts)) == typeof.(e_int)
        @test typeof.(internal_energy_vapor.(param_set, ts)) == typeof.(e_int)
        @test typeof.(internal_energy_liquid.(param_set, ts)) == typeof.(e_int)
        @test typeof.(internal_energy_ice.(param_set, ts)) == typeof.(e_int)
        @test typeof.(latent_heat_vapor.(param_set, ts)) == typeof.(e_int)
        @test typeof.(latent_heat_sublim.(param_set, ts)) == typeof.(e_int)
        @test typeof.(latent_heat_fusion.(param_set, ts)) == typeof.(e_int)
        @test typeof.(q_vap_saturation.(param_set, ts)) == typeof.(e_int)
        @test typeof.(q_vap_saturation_liquid.(param_set, ts)) == typeof.(e_int)
        @test typeof.(q_vap_saturation_ice.(param_set, ts)) == typeof.(e_int)
        @test typeof.(saturation_excess.(param_set, ts)) == typeof.(e_int)
        @test typeof.(liquid_fraction.(param_set, ts)) == typeof.(e_int)
        @test typeof.(liquid_ice_pottemp.(param_set, ts)) == typeof.(e_int)
        @test typeof.(dry_pottemp.(param_set, ts)) == typeof.(e_int)
        @test typeof.(exner.(param_set, ts)) == typeof.(e_int)
        @test typeof.(liquid_ice_pottemp_sat.(param_set, ts)) == typeof.(e_int)
        @test typeof.(specific_volume.(param_set, ts)) == typeof.(e_int)
        @test typeof.(supersaturation.(param_set, ts, Ice())) == typeof.(e_int)
        @test typeof.(supersaturation.(param_set, ts, Liquid())) ==
              typeof.(e_int)
        @test typeof.(virtual_pottemp.(param_set, ts)) == typeof.(e_int)
        @test typeof.(specific_entropy.(param_set, ts)) == typeof.(e_int)
        @test eltype.(gas_constants.(param_set, ts)) == typeof.(e_int)

        @test typeof.(total_specific_enthalpy.(param_set, ts, e_tot)) ==
              typeof.(e_int)
        @test typeof.(moist_static_energy.(param_set, ts, e_pot)) ==
              typeof.(e_int)
        @test typeof.(getproperty.(PhasePartition.(param_set, ts), :tot)) ==
              typeof.(e_int)
        @test typeof.(virtual_dry_static_energy.(param_set, ts, e_pot)) ==
              typeof.(e_int)
    end
    @test typeof(
        q_vap_from_RH_liquid(param_set, FT(100000), FT(275), FT(0.8)),
    ) == FT
end 