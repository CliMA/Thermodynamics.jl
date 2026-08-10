"""
# Functional type stability tests

These tests focus on the functional API 
"""

function test_zero_allocations(param_set, T, ρ, q_tot, q_liq, q_ice, e_int)
    # Warm up inside the function
    TD.air_pressure(param_set, T, ρ, q_tot, q_liq, q_ice)
    TD.internal_energy(param_set, T, q_tot, q_liq, q_ice)
    TD.air_temperature(param_set, TD.ρe(), e_int, q_tot, q_liq, q_ice)
    TD.saturation_adjustment(param_set, TD.ρe(), ρ, e_int, q_tot)

    # Return allocations
    alloc1 = @allocated TD.air_pressure(param_set, T, ρ, q_tot, q_liq, q_ice)
    alloc2 = @allocated TD.internal_energy(param_set, T, q_tot, q_liq, q_ice)
    alloc3 = @allocated TD.air_temperature(param_set, TD.ρe(), e_int, q_tot, q_liq, q_ice)
    alloc4 = @allocated TD.saturation_adjustment(param_set, TD.ρe(), ρ, e_int, q_tot)
    return alloc1, alloc2, alloc3, alloc4
end

"""
    saturation_adjustment_allocations(param_set, ρ, p, e_int, h, θ_p, θ_ρ, q_tot)

Allocations of the convenience `saturation_adjustment` method for each of the six
formulations. These are the branch-free kernels targeted at GPUs, so all six must be
allocation-free, not just `ρe`.
"""
function saturation_adjustment_allocations(
    param_set,
    ρ,
    p,
    e_int,
    h,
    θ_p,
    θ_ρ,
    q_tot,
)
    cases = (
        (TD.ρe(), ρ, e_int),
        (TD.pe(), p, e_int),
        (TD.ph(), p, h),
        (TD.pρ(), p, ρ),
        (TD.pθ_li(), p, θ_p),
        (TD.ρθ_li(), ρ, θ_ρ),
    )
    # Warm up every method before measuring
    for (indep_vars, var₁, var₂) in cases
        TD.saturation_adjustment(param_set, indep_vars, var₁, var₂, q_tot)
    end
    return map(cases) do (indep_vars, var₁, var₂)
        @allocated TD.saturation_adjustment(param_set, indep_vars, var₁, var₂, q_tot)
    end
end

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
            @test @inferred(TD.relative_humidity(param_set, T, p)) isa FT
            @test @inferred(TD.relative_humidity(param_set, T, p, q_tot)) isa FT
            @test @inferred(TD.relative_humidity(param_set, T, p, q_tot, q_liq)) isa FT
            @test @inferred(
                TD.air_temperature(param_set, TD.ρe(), e_int, q_tot, q_liq, q_ice)
            ) isa FT
            @test @inferred(
                TD.air_temperature(param_set, TD.ph(), h, q_tot, q_liq, q_ice)
            ) isa FT
            @test @inferred(
                TD.air_temperature(param_set, TD.pρ(), p, ρ, q_tot, q_liq, q_ice)
            ) isa FT

            # Saturation-adjustment surface: all six formulations, plus the phase
            # partitioning and analytic derivatives that the solvers call every iteration.
            θ_p = TD.liquid_ice_pottemp_given_pressure(
                param_set,
                T,
                p,
                q_tot,
                q_liq,
                q_ice,
            )
            θ_ρ = TD.liquid_ice_pottemp(param_set, T, ρ, q_tot, q_liq, q_ice)
            for (indep_vars, var₁, var₂) in (
                (TD.ρe(), ρ, e_int),
                (TD.pe(), p, e_int),
                (TD.ph(), p, h),
                (TD.pρ(), p, ρ),
                (TD.pθ_li(), p, θ_p),
                (TD.ρθ_li(), ρ, θ_ρ),
            )
                sol = @inferred TD.saturation_adjustment(
                    param_set,
                    indep_vars,
                    var₁,
                    var₂,
                    q_tot,
                )
                @test sol.T isa FT
            end
            @test @inferred(TD.condensate_partition(param_set, T, ρ, q_tot)) isa
                  Tuple{FT, FT}
            @test @inferred(TD.liquid_fraction(param_set, T, q_liq, q_ice)) isa FT
            @test @inferred(TD.liquid_fraction_ramp(param_set, T)) isa FT
            @test @inferred(TD.∂q_vap_sat_∂T(param_set, T, ρ)) isa FT
            @test @inferred(TD.∂e_int_∂T_sat_ρ(param_set, T, ρ, q_tot)) isa FT
            @test @inferred(TD.∂e_int_∂T_sat_p(param_set, T, p, q_tot)) isa FT
            @test @inferred(TD.∂h_∂T_sat_p(param_set, T, p, q_tot)) isa FT
            @test @inferred(TD.∂θ_li_∂T_sat_ρ(param_set, T, ρ, q_tot)) isa FT
            @test @inferred(TD.∂θ_li_∂T_sat_p(param_set, T, p, q_tot)) isa FT
            @test @inferred(TD.∂p_∂T_sat_ρ(param_set, T, ρ, q_tot)) isa FT
            @test @inferred(TD.latent_heat_mixed(param_set, T, FT(0.5))) isa FT
            @test @inferred(TD.entropy(param_set, p, T, q_tot, q_liq, q_ice)) isa FT
            @test @inferred(TD.soundspeed_air(param_set, T, q_tot, q_liq, q_ice)) isa FT
            @test @inferred(TD.exner(param_set, T, ρ, q_tot, q_liq, q_ice)) isa FT
        end
        @testset "Zero allocations ($FT)" begin
            allocs = test_zero_allocations(param_set, T, ρ, q_tot, q_liq, q_ice, e_int)
            @test allocs[1] == 0
            @test allocs[2] == 0
            @test allocs[3] == 0
            @test allocs[4] == 0
        end
        @testset "Zero allocations, all formulations ($FT)" begin
            θ_p = TD.liquid_ice_pottemp_given_pressure(
                param_set,
                T,
                p,
                q_tot,
                q_liq,
                q_ice,
            )
            θ_ρ = TD.liquid_ice_pottemp(param_set, T, ρ, q_tot, q_liq, q_ice)
            allocs = saturation_adjustment_allocations(
                param_set,
                ρ,
                p,
                e_int,
                h,
                θ_p,
                θ_ρ,
                q_tot,
            )
            for a in allocs
                @test a == 0
            end
        end
    end
end
