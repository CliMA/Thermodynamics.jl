"""
# Functional type stability tests

These tests focus on the functional API 
"""

@testset "Thermodynamics - type stability (functional)" begin
    for FT in (Float32, Float64)
        param_set = FT == Float64 ? param_set_Float64 : param_set_Float32

        T = FT(300)
        ρ = FT(1.1)
        p = TD.air_pressure(param_set, T, ρ)
        q_tot = FT(0.02)
        q_liq = FT(0.003)
        q_ice = FT(0.001)
        e_int = TD.internal_energy(param_set, T, q_tot, q_liq, q_ice)
        h = TD.enthalpy(param_set, T, q_tot, q_liq, q_ice)

        @testset "Scalar outputs ($FT)" begin
            @test @inferred(TD.air_pressure(param_set, T, ρ, q_tot, q_liq, q_ice)) isa FT
            @test @inferred(TD.air_density(param_set, T, p, q_tot, q_liq, q_ice)) isa FT
            @test @inferred(TD.gas_constant_air(param_set, q_tot, q_liq, q_ice)) isa FT
            @test @inferred(TD.cp_m(param_set, q_tot, q_liq, q_ice)) isa FT
            @test @inferred(TD.cv_m(param_set, q_tot, q_liq, q_ice)) isa FT
            @test @inferred(TD.internal_energy(param_set, T, q_tot, q_liq, q_ice)) isa FT
            @test @inferred(TD.enthalpy(param_set, T, q_tot, q_liq, q_ice)) isa FT
            @test @inferred(TD.saturation_vapor_pressure(param_set, T, TD.Liquid())) isa FT
            @test @inferred(TD.saturation_vapor_pressure(param_set, T, TD.Ice())) isa FT
            @test @inferred(TD.q_vap_saturation(param_set, T, ρ)) isa FT
            @test @inferred(TD.q_vap_saturation_from_pressure(param_set, q_tot, p, T)) isa
                  FT
            @test @inferred(TD.relative_humidity(param_set, T, p, q_tot, q_liq, q_ice)) isa
                  FT
            @test @inferred(
                TD.air_temperature(param_set, TD.ρe(), e_int, q_tot, q_liq, q_ice)
            ) isa FT
            @test @inferred(
                TD.air_temperature(param_set, TD.ph(), h, q_tot, q_liq, q_ice)
            ) isa FT
            @test @inferred(
                TD.air_temperature(param_set, TD.pρ(), p, ρ, q_tot, q_liq, q_ice)
            ) isa FT
        end
    end
end
