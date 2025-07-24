"""
# Exceptions Test Suite

This file contains tests for error handling on failed convergence.
"""

@testset "Thermodynamics - Exceptions on Failed Convergence" begin

    ArrayType = Array{Float64}
    FT = eltype(ArrayType)
    param_set = TP.ThermodynamicsParameters(FT)
    profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
    (; T, p, e_int, ρ, θ_liq_ice, phase_type) = profiles
    (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = profiles

    @test_throws ErrorException TD.saturation_adjustment.(
        RS.NewtonsMethod,
        param_set,
        e_int,
        ρ,
        q_tot,
        Ref(phase_type),
        2,
        FT(1e-10),
    )

    @test_throws ErrorException TD.saturation_adjustment.(
        RS.SecantMethod,
        param_set,
        e_int,
        ρ,
        q_tot,
        Ref(phase_type),
        2,
        FT(1e-10),
    )

    @test_throws ErrorException TD.saturation_adjustment_given_peq.(
        RS.SecantMethod,
        param_set,
        p,
        e_int,
        q_tot,
        Ref(phase_type),
        2,
        FT(1e-10),
    )

    T_virt = T # should not matter: testing for non-convergence
    @test_throws ErrorException TD.temperature_and_humidity_given_TᵥρRH.(
        param_set,
        T_virt,
        ρ,
        RH,
        Ref(phase_type),
        2,
        RS.ResidualTolerance(FT(1e-10)),
    )

    @test_throws ErrorException TD.air_temperature_given_ρθq_nonlinear.(
        param_set,
        ρ,
        θ_liq_ice,
        2,
        RS.ResidualTolerance(FT(1e-10)),
        q_pt,
    )

    @test_throws ErrorException TD.saturation_adjustment_given_ρθq.(
        param_set,
        ρ,
        θ_liq_ice,
        q_tot,
        Ref(phase_type),
        2,
        RS.ResidualTolerance(FT(1e-10)),
    )

    @test_throws ErrorException TD.saturation_adjustment_given_pθq.(
        RS.SecantMethod,
        param_set,
        p,
        θ_liq_ice,
        q_tot,
        Ref(phase_type),
        2,
        FT(1e-10),
    )

    @test_throws ErrorException TD.saturation_adjustment_given_pθq.(
        RS.NewtonsMethodAD,
        param_set,
        p,
        θ_liq_ice,
        q_tot,
        Ref(phase_type),
        2,
        FT(1e-10),
    )

    @test_throws ErrorException TD.saturation_adjustment_ρpq.(
        RS.NewtonsMethodAD,
        param_set,
        ρ,
        p,
        q_tot,
        Ref(phase_type),
        2,
        FT(1e-10),
    )

end
