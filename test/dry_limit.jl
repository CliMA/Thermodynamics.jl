"""
# Dry Limit Test Suite

This file contains tests for behavior when approaching dry air conditions.
"""

@testset "Thermodynamics - Dry Limit" begin

    ArrayType = Array{Float64}
    FT = eltype(ArrayType)
    param_set = TP.ThermodynamicsParameters(FT)
    profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
    (; T, p, e_int, ρ, θ_liq_ice, phase_type) = profiles
    (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = profiles

    # PhasePartition test is noisy, so do this only once:
    ts_dry = PhaseDry(param_set, first(e_int), first(ρ))
    ts_eq =
        PhaseEquil_ρeq(param_set, first(ρ), first(e_int), typeof(first(ρ))(0))
    @test PhasePartition(param_set, ts_eq).tot ≈
          PhasePartition(param_set, ts_dry).tot
    @test PhasePartition(param_set, ts_eq).liq ≈
          PhasePartition(param_set, ts_dry).liq
    @test PhasePartition(param_set, ts_eq).ice ≈
          PhasePartition(param_set, ts_dry).ice

    @test mixing_ratios(param_set, ts_eq).tot ≈
          mixing_ratios(param_set, ts_dry).tot
    @test mixing_ratios(param_set, ts_eq).liq ≈
          mixing_ratios(param_set, ts_dry).liq
    @test mixing_ratios(param_set, ts_eq).ice ≈
          mixing_ratios(param_set, ts_dry).ice

    ts_dry = PhaseDry.(param_set, e_int, ρ)
    ts_eq = PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot .* 0)

    @test all(
        gas_constant_air.(param_set, ts_eq) .≈
        gas_constant_air.(param_set, ts_dry),
    )
    @test all(
        relative_humidity.(param_set, ts_eq) .≈
        relative_humidity.(param_set, ts_dry),
    )
    @test all(
        air_pressure.(param_set, ts_eq) .≈ air_pressure.(param_set, ts_dry),
    )
    @test all(air_density.(param_set, ts_eq) .≈ air_density.(param_set, ts_dry))
    @test all(
        specific_volume.(param_set, ts_eq) .≈
        specific_volume.(param_set, ts_dry),
    )
    @test all(
        total_specific_humidity.(param_set, ts_eq) .≈
        total_specific_humidity.(param_set, ts_dry),
    )
    @test all(
        liquid_specific_humidity.(param_set, ts_eq) .≈
        liquid_specific_humidity.(param_set, ts_dry),
    )
    @test all(
        ice_specific_humidity.(param_set, ts_eq) .≈
        ice_specific_humidity.(param_set, ts_dry),
    )
    @test all(cp_m.(param_set, ts_eq) .≈ cp_m.(param_set, ts_dry))
    @test all(cv_m.(param_set, ts_eq) .≈ cv_m.(param_set, ts_dry))
    @test all(
        air_temperature.(param_set, ts_eq) .≈
        air_temperature.(param_set, ts_dry),
    )
    @test all(
        internal_energy.(param_set, ts_eq) .≈
        internal_energy.(param_set, ts_dry),
    )
    @test all(
        internal_energy_sat.(param_set, ts_eq) .≈
        internal_energy_sat.(param_set, ts_dry),
    )
    @test all(
        internal_energy_dry.(param_set, ts_eq) .≈
        internal_energy_dry.(param_set, ts_dry),
    )
    @test all(
        internal_energy_vapor.(param_set, ts_eq) .≈
        internal_energy_vapor.(param_set, ts_dry),
    )
    @test all(
        internal_energy_liquid.(param_set, ts_eq) .≈
        internal_energy_liquid.(param_set, ts_dry),
    )
    @test all(
        internal_energy_ice.(param_set, ts_eq) .≈
        internal_energy_ice.(param_set, ts_dry),
    )
    @test all(
        soundspeed_air.(param_set, ts_eq) .≈ soundspeed_air.(param_set, ts_dry),
    )
    @test all(
        supersaturation.(param_set, ts_eq, Ice()) .≈
        supersaturation.(param_set, ts_dry, Ice()),
    )
    @test all(
        supersaturation.(param_set, ts_eq, Liquid()) .≈
        supersaturation.(param_set, ts_dry, Liquid()),
    )
    @test all(
        latent_heat_vapor.(param_set, ts_eq) .≈
        latent_heat_vapor.(param_set, ts_dry),
    )
    @test all(
        latent_heat_sublim.(param_set, ts_eq) .≈
        latent_heat_sublim.(param_set, ts_dry),
    )
    @test all(
        latent_heat_fusion.(param_set, ts_eq) .≈
        latent_heat_fusion.(param_set, ts_dry),
    )
    @test all(
        q_vap_saturation.(param_set, ts_eq) .≈
        q_vap_saturation.(param_set, ts_dry),
    )
    @test all(
        q_vap_saturation_liquid.(param_set, ts_eq) .≈
        q_vap_saturation_liquid.(param_set, ts_dry),
    )
    @test all(
        q_vap_saturation_ice.(param_set, ts_eq) .≈
        q_vap_saturation_ice.(param_set, ts_dry),
    )
    @test all(
        saturation_excess.(param_set, ts_eq) .≈
        saturation_excess.(param_set, ts_dry),
    )
    @test all(
        liquid_fraction.(param_set, ts_eq) .≈
        liquid_fraction.(param_set, ts_dry),
    )
    @test all(
        liquid_ice_pottemp.(param_set, ts_eq) .≈
        liquid_ice_pottemp.(param_set, ts_dry),
    )
    @test all(dry_pottemp.(param_set, ts_eq) .≈ dry_pottemp.(param_set, ts_dry))
    @test all(
        virtual_pottemp.(param_set, ts_eq) .≈
        virtual_pottemp.(param_set, ts_dry),
    )
    @test all(
        specific_entropy.(param_set, ts_eq) .≈
        specific_entropy.(param_set, ts_dry),
    )
    @test all(
        liquid_ice_pottemp_sat.(param_set, ts_eq) .≈
        liquid_ice_pottemp_sat.(param_set, ts_dry),
    )
    @test all(exner.(param_set, ts_eq) .≈ exner.(param_set, ts_dry))

    @test all(
        saturation_vapor_pressure.(param_set, ts_eq, Ice()) .≈
        saturation_vapor_pressure.(param_set, ts_dry, Ice()),
    )
    @test all(
        saturation_vapor_pressure.(param_set, ts_eq, Liquid()) .≈
        saturation_vapor_pressure.(param_set, ts_dry, Liquid()),
    )
    @test all(
        first.(gas_constants.(param_set, ts_eq)) ≈
        first.(gas_constants.(param_set, ts_dry)),
    )
    @test all(
        last.(gas_constants.(param_set, ts_eq)) ≈
        last.(gas_constants.(param_set, ts_dry)),
    )

end
