"""
Tests for default convenience methods of saturation_adjustment.

This suite verifies that the convenience methods with default settings (without 
explicit solver method type) work correctly across all formulations:

- ρe: Uses RS.NewtonsMethod with forced_fixed_iters=true
- Other formulations: Use RS.SecantMethod
"""

@testset "Thermodynamics - saturation_adjustment default methods" begin
    for FT in (Float32, Float64)
        param_set = FT == Float64 ? param_set_Float64 : param_set_Float32

        # Default error tolerance for all methods (0.1 K)
        default_atol_temperature = 0.1

        approx_tight(a, b) =
            isapprox(a, b; rtol = FT(100) * eps(FT), atol = FT(100) * eps(FT))

        function check_partition(Tsol, ρ_eff, q0, ql, qi)
            (ql_exp, qi_exp) = TD.condensate_partition(param_set, Tsol, ρ_eff, q0)
            @test approx_tight(ql, ql_exp)
            @test approx_tight(qi, qi_exp)
            return nothing
        end

        @testset "ρe default method accuracy ($FT)" begin
            profiles = TestedProfiles.PhaseEquilProfiles(param_set, Array{FT})
            (; T, p, ρ, q_tot) = profiles

            # Sample ~120 points across the full profile grid
            idxs = unique(round.(Int, range(1, length(T), length = 250)))

            function targets(i::Int)
                T0 = T[i]
                p0 = p[i]
                ρ0 = ρ[i]
                q0 = q_tot[i]

                # ρ-based targets
                (q_liq_ρ, q_ice_ρ) = TD.condensate_partition(param_set, T0, ρ0, q0)
                e_int_ρ = TD.internal_energy_sat(param_set, T0, ρ0, q0)

                return (; T0, p0, ρ0, q0, e_int_ρ)
            end

            # Test default ρe method (forced_fixed_iters=true)
            max_err = FT(0)
            chunk_errs = FT[]
            for i in idxs
                inp = targets(i)

                # Use default convenience method (no explicit solver type)
                res = TD.saturation_adjustment(
                    param_set,
                    TD.ρe(),
                    inp.ρ0,
                    inp.e_int_ρ,
                    inp.q0,
                )

                err = abs(res.T - inp.T0)
                max_err = max(max_err, err)
                push!(chunk_errs, err)

                if err > FT(default_atol_temperature)
                    @info "Outlier detected" FT i err T_ref = inp.T0 p = inp.p0 q_tot =
                        inp.q0 ρ = inp.ρ0
                end

                check_partition(res.T, inp.ρ0, inp.q0, res.q_liq, res.q_ice)
            end
            @info "Max error ρe default method" FT max_err

            # Allow a few outliers to exceed strict tolerance
            @test max_err < 5 * default_atol_temperature

            # Check percentage of points with error < default_atol_temperature
            n_strict = count(e -> e < FT(default_atol_temperature), chunk_errs)
            pct_strict = n_strict / length(chunk_errs) * 100
            @info "Percentage of points with error < $(default_atol_temperature) K" FT pct_strict
            @test pct_strict > 98
        end

        @testset "Other formulations default methods ($FT)" begin
            profiles = TestedProfiles.PhaseEquilProfiles(param_set, Array{FT})
            (; T, p, ρ, q_tot) = profiles

            # Sample fewer points for other formulations (faster test)
            idxs = unique(round.(Int, range(1, length(T), length = 100)))

            function targets(i::Int)
                T0 = T[i]
                p0 = p[i]
                ρ0 = ρ[i]
                q0 = q_tot[i]

                # ρ-based targets
                (q_liq_ρ, q_ice_ρ) = TD.condensate_partition(param_set, T0, ρ0, q0)
                e_int_ρ = TD.internal_energy_sat(param_set, T0, ρ0, q0)
                p_ρ = TD.air_pressure(param_set, T0, ρ0, q0, q_liq_ρ, q_ice_ρ)
                θ_ρ = TD.liquid_ice_pottemp(param_set, T0, ρ0, q0, q_liq_ρ, q_ice_ρ)

                # p-based targets
                ρ_p = TD.air_density(param_set, T0, p0, q0)
                (q_liq_p, q_ice_p) = TD.condensate_partition(param_set, T0, ρ_p, q0)
                e_int_p = TD.internal_energy_sat(param_set, T0, ρ_p, q0)
                h_p = TD.enthalpy_sat(param_set, T0, ρ_p, q0)
                θ_p =
                    TD.liquid_ice_pottemp_given_pressure(
                        param_set,
                        T0,
                        p0,
                        q0,
                        q_liq_p,
                        q_ice_p,
                    )

                return (;
                    T0,
                    p0,
                    ρ0,
                    q0,
                    e_int_ρ,
                    p_ρ,
                    θ_ρ,
                    e_int_p,
                    h_p,
                    θ_p,
                )
            end

            # Test pe default method
            for i in idxs
                inp = targets(i)
                let (; T, q_liq, q_ice) = TD.saturation_adjustment(
                        param_set,
                        TD.pe(),
                        inp.p0,
                        inp.e_int_p,
                        inp.q0,
                    )
                    @test isfinite(T)
                    @test isapprox(
                        T,
                        inp.T0;
                        atol = FT(default_atol_temperature),
                        rtol = FT(0),
                    )
                    ρ_eff = TD.air_density(param_set, T, inp.p0, inp.q0)
                    check_partition(T, ρ_eff, inp.q0, q_liq, q_ice)
                end
            end

            # Test ph default method
            for i in idxs
                inp = targets(i)
                let (; T, q_liq, q_ice) = TD.saturation_adjustment(
                        param_set,
                        TD.ph(),
                        inp.p0,
                        inp.h_p,
                        inp.q0,
                    )
                    @test isfinite(T)
                    @test isapprox(
                        T,
                        inp.T0;
                        atol = FT(default_atol_temperature),
                        rtol = FT(0),
                    )
                    ρ_eff = TD.air_density(param_set, T, inp.p0, inp.q0)
                    check_partition(T, ρ_eff, inp.q0, q_liq, q_ice)
                end
            end

            # Test pρ default method
            for i in idxs
                inp = targets(i)
                let (; T, q_liq, q_ice) = TD.saturation_adjustment(
                        param_set,
                        TD.pρ(),
                        inp.p_ρ,
                        inp.ρ0,
                        inp.q0,
                    )
                    @test isfinite(T)
                    @test isapprox(
                        T,
                        inp.T0;
                        atol = FT(default_atol_temperature),
                        rtol = FT(0),
                    )
                    check_partition(T, inp.ρ0, inp.q0, q_liq, q_ice)
                end
            end

            # Test pθ_li default method
            for i in idxs
                inp = targets(i)
                let (; T, q_liq, q_ice) = TD.saturation_adjustment(
                        param_set,
                        TD.pθ_li(),
                        inp.p0,
                        inp.θ_p,
                        inp.q0,
                    )
                    @test isfinite(T)
                    @test isapprox(
                        T,
                        inp.T0;
                        atol = FT(default_atol_temperature),
                        rtol = FT(0),
                    )
                    ρ_eff = TD.air_density(param_set, T, inp.p0, inp.q0)
                    check_partition(T, ρ_eff, inp.q0, q_liq, q_ice)
                end
            end

            # Test ρθ_li default method
            for i in idxs
                inp = targets(i)
                let (; T, q_liq, q_ice) = TD.saturation_adjustment(
                        param_set,
                        TD.ρθ_li(),
                        inp.ρ0,
                        inp.θ_ρ,
                        inp.q0,
                    )
                    @test isfinite(T)
                    @test isapprox(
                        T,
                        inp.T0;
                        atol = FT(default_atol_temperature),
                        rtol = FT(0),
                    )
                    check_partition(T, inp.ρ0, inp.q0, q_liq, q_ice)
                end
            end
        end
    end
end
