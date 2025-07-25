"""
# Dry Limit Test Suite

This file contains tests for behavior when approaching dry air conditions.
"""

@testset "Thermodynamics - Dry Limit" begin
    ArrayType = Array{Float64}
    FT = eltype(ArrayType)
    param_set = FT == Float64 ? param_set_Float64 : param_set_Float32

    profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
    (; e_int, ρ, q_tot) = profiles

    @testset "PhasePartition and Mixing Ratios" begin
        ts_dry = PhaseDry(param_set, first(e_int), first(ρ))
        ts_eq = PhaseEquil_ρeq(
            param_set,
            first(ρ),
            first(e_int),
            typeof(first(ρ))(0),
        )
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
    end

    @testset "Thermodynamic Properties" begin
        ts_dry = PhaseDry.(param_set, e_int, ρ)
        ts_eq = PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot .* 0)

        for (func, args) in (
            (gas_constant_air, ()),
            (relative_humidity, ()),
            (air_pressure, ()),
            (air_density, ()),
            (specific_volume, ()),
            (total_specific_humidity, ()),
            (liquid_specific_humidity, ()),
            (ice_specific_humidity, ()),
            (cp_m, ()),
            (cv_m, ()),
            (air_temperature, ()),
            (internal_energy, ()),
            (internal_energy_sat, ()),
            (internal_energy_dry, ()),
            (internal_energy_vapor, ()),
            (internal_energy_liquid, ()),
            (internal_energy_ice, ()),
            (soundspeed_air, ()),
            (latent_heat_vapor, ()),
            (latent_heat_sublim, ()),
            (latent_heat_fusion, ()),
            (q_vap_saturation, ()),
            (q_vap_saturation_liquid, ()),
            (q_vap_saturation_ice, ()),
            (saturation_excess, ()),
            (liquid_fraction, ()),
            (liquid_ice_pottemp, ()),
            (dry_pottemp, ()),
            (virtual_pottemp, ()),
            (specific_entropy, ()),
            (liquid_ice_pottemp_sat, ()),
            (exner, ()),
        )
            @test all(
                func.(param_set, ts_eq, args...) .≈
                func.(param_set, ts_dry, args...),
            )
        end

        # Test supersaturation separately with different phase arguments
        @test all(
            supersaturation.(param_set, ts_eq, Ice()) .≈
            supersaturation.(param_set, ts_dry, Ice()),
        )
        @test all(
            supersaturation.(param_set, ts_eq, Liquid()) .≈
            supersaturation.(param_set, ts_dry, Liquid()),
        )
        @test all(
            saturation_vapor_pressure.(param_set, ts_eq, Ice()) .≈
            saturation_vapor_pressure.(param_set, ts_dry, Ice()),
        )
        @test all(
            saturation_vapor_pressure.(param_set, ts_eq, Liquid()) .≈
            saturation_vapor_pressure.(param_set, ts_dry, Liquid()),
        )
        @test all(
            first.(gas_constants.(param_set, ts_eq)) .≈
            first.(gas_constants.(param_set, ts_dry)),
        )
        @test all(
            last.(gas_constants.(param_set, ts_eq)) .≈
            last.(gas_constants.(param_set, ts_dry)),
        )
    end
end
