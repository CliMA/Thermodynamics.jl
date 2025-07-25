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
    (; q_tot, q_pt, RH) = profiles

    maxiter = 2
    tol = FT(1e-10)
    T_virt = T

    @testset "Saturation Adjustment" begin
        @test_throws ErrorException TD.saturation_adjustment.(
            RS.NewtonsMethod,
            param_set,
            e_int,
            ρ,
            q_tot,
            Ref(phase_type),
            maxiter,
            tol,
        )

        @test_throws ErrorException TD.saturation_adjustment.(
            RS.SecantMethod,
            param_set,
            e_int,
            ρ,
            q_tot,
            Ref(phase_type),
            maxiter,
            tol,
        )

        @test_throws ErrorException TD.saturation_adjustment_given_peq.(
            RS.SecantMethod,
            param_set,
            p,
            e_int,
            q_tot,
            Ref(phase_type),
            maxiter,
            tol,
        )

        @test_throws ErrorException TD.saturation_adjustment_given_ρθq.(
            param_set,
            ρ,
            θ_liq_ice,
            q_tot,
            Ref(phase_type),
            maxiter,
            RS.ResidualTolerance(tol),
        )

        @test_throws ErrorException TD.saturation_adjustment_given_pθq.(
            RS.SecantMethod,
            param_set,
            p,
            θ_liq_ice,
            q_tot,
            Ref(phase_type),
            maxiter,
            tol,
        )

        @test_throws ErrorException TD.saturation_adjustment_given_pθq.(
            RS.NewtonsMethodAD,
            param_set,
            p,
            θ_liq_ice,
            q_tot,
            Ref(phase_type),
            maxiter,
            tol,
        )

        @test_throws ErrorException TD.saturation_adjustment_ρpq.(
            RS.NewtonsMethodAD,
            param_set,
            ρ,
            p,
            q_tot,
            Ref(phase_type),
            maxiter,
            tol,
        )
    end

    @testset "Temperature and Humidity" begin
        @test_throws ErrorException TD.temperature_and_humidity_given_TᵥρRH.(
            param_set,
            T_virt,
            ρ,
            RH,
            Ref(phase_type),
            maxiter,
            RS.ResidualTolerance(tol),
        )

        @test_throws ErrorException TD.air_temperature_given_ρθq_nonlinear.(
            param_set,
            ρ,
            θ_liq_ice,
            maxiter,
            RS.ResidualTolerance(tol),
            q_pt,
        )
    end
end
