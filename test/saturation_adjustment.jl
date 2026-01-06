"""
Tests for saturation adjustment (functional API)

Focus:
- convergence (typical cases)
- internal consistency (returned state reproduces the target thermo variable)
"""

@testset "Thermodynamics - saturation_adjustment (functional)" begin
    for FT in (Float32, Float64)
        param_set = FT == Float64 ? param_set_Float64 : param_set_Float32

        # Choose a representative density and temperature
        ρ0 = FT(1.0)
        T0 = FT(290)
        p0_dry = TD.air_pressure(param_set, T0, ρ0)

        @testset "Unsaturated returns zero condensate ($FT)" begin
            q_tot = FT(1e-3)
            e_int = TD.internal_energy(param_set, T0, q_tot)
            (T, q_liq, q_ice) = TD.saturation_adjustment(
                RS.SecantMethod,
                param_set,
                TD.ρe(),
                ρ0,
                e_int,
                q_tot,
                50,
                FT(1e-10),
            )
            @test isapprox(T, T0; atol = FT(atol_temperature), rtol = FT(0))
            @test q_liq == zero(FT)
            @test q_ice == zero(FT)
        end

        @testset "Saturated equilibrium is a fixed point ($FT)" begin
            # Make a saturated case by picking q_tot slightly above saturation at (T0, ρ0)
            q_tot = FT(1.2) * TD.q_vap_saturation(param_set, T0, ρ0)
            e_int_sat = TD.internal_energy_sat(param_set, T0, ρ0, q_tot)
            (q_liq0, q_ice0) = TD.condensate_partition(param_set, T0, ρ0, q_tot)
            p0 = TD.air_pressure(param_set, T0, ρ0, q_tot, q_liq0, q_ice0)
            h_sat = TD.enthalpy_sat(param_set, T0, ρ0, q_tot)
            θ_liq_ice_p = TD.liquid_ice_pottemp_given_pressure(
                param_set,
                T0,
                p0,
                q_tot,
                q_liq0,
                q_ice0,
            )
            θ_liq_ice_ρ = TD.liquid_ice_pottemp(param_set, T0, ρ0, q_tot, q_liq0, q_ice0)

            # ρeq
            let (T, q_liq, q_ice) = TD.saturation_adjustment(
                    RS.NewtonsMethod,
                    param_set,
                    TD.ρe(),
                    ρ0,
                    e_int_sat,
                    q_tot,
                    50,
                    FT(1e-10),
                )
                @test isapprox(T, T0; atol = FT(atol_temperature), rtol = FT(0))
                @test isapprox(q_liq, q_liq0; atol = FT(0), rtol = FT(1e-6))
                @test isapprox(q_ice, q_ice0; atol = FT(0), rtol = FT(1e-6))
                @test isapprox(
                    TD.internal_energy_sat(param_set, T, ρ0, q_tot),
                    e_int_sat;
                    rtol = FT(1e-6),
                )
            end

            # peq
            let (T, q_liq, q_ice) = TD.saturation_adjustment(
                    RS.SecantMethod,
                    param_set,
                    TD.pe(),
                    p0,
                    e_int_sat,
                    q_tot,
                    80,
                    FT(1e-10),
                )
                ρ = TD.air_density(param_set, T, p0, q_tot)
                @test isapprox(T, T0; atol = FT(atol_temperature), rtol = FT(0))
                @test isapprox(
                    TD.internal_energy_sat(param_set, T, ρ, q_tot),
                    e_int_sat;
                    rtol = FT(1e-6),
                )
                @test q_liq + q_ice ≈ TD.saturation_excess(param_set, T, ρ, q_tot)
            end

            # phq
            let (T, q_liq, q_ice) = TD.saturation_adjustment(
                    RS.SecantMethod,
                    param_set,
                    TD.ph(),
                    p0,
                    h_sat,
                    q_tot,
                    80,
                    FT(1e-10),
                )
                ρ = TD.air_density(param_set, T, p0, q_tot)
                @test isapprox(T, T0; atol = FT(atol_temperature), rtol = FT(0))
                @test isapprox(
                    TD.enthalpy_sat(param_set, T, ρ, q_tot),
                    h_sat;
                    rtol = FT(1e-6),
                )
                @test q_liq + q_ice ≈ TD.saturation_excess(param_set, T, ρ, q_tot)
            end

            # pρq (solve for T at saturation such that computed pressure matches)
            let (T, q_liq, q_ice) = TD.saturation_adjustment(
                    RS.SecantMethod,
                    param_set,
                    TD.pρ(),
                    p0,
                    ρ0,
                    q_tot,
                    80,
                    FT(1e-10),
                )
                @test isapprox(T, T0; atol = FT(atol_temperature), rtol = FT(0))
                @test isapprox(
                    TD.air_pressure(param_set, T, ρ0, q_tot, q_liq, q_ice),
                    p0;
                    rtol = FT(1e-6),
                )
            end

            # pθ_liq_ice_q
            let (T, q_liq, q_ice) = TD.saturation_adjustment(
                    RS.SecantMethod,
                    param_set,
                    TD.pθ_li(),
                    p0,
                    θ_liq_ice_p,
                    q_tot,
                    80,
                    FT(1e-10),
                )
                @test isapprox(T, T0; atol = FT(atol_temperature), rtol = FT(0))
                θ_chk = TD.liquid_ice_pottemp_given_pressure(
                    param_set,
                    T,
                    p0,
                    q_tot,
                    q_liq,
                    q_ice,
                )
                @test isapprox(θ_chk, θ_liq_ice_p; rtol = FT(1e-6))
            end

            # ρθ_liq_ice_q
            let (T, q_liq, q_ice) = TD.saturation_adjustment(
                    RS.SecantMethod,
                    param_set,
                    TD.ρθ_li(),
                    ρ0,
                    θ_liq_ice_ρ,
                    q_tot,
                    80,
                    FT(1e-10),
                )
                @test isapprox(T, T0; atol = FT(atol_temperature), rtol = FT(0))
                θ_chk = TD.liquid_ice_pottemp(param_set, T, ρ0, q_tot, q_liq, q_ice)
                @test isapprox(θ_chk, θ_liq_ice_ρ; rtol = FT(1e-6))
            end
        end
    end
end
