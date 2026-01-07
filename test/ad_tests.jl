"""
Automatic differentiation tests for Thermodynamics.jl
"""

using ForwardDiff
using Thermodynamics: ρe, pe, ph


@testset "Thermodynamics - AD compatibility" begin
    # Use Float64 parameter set from test_common.jl
    FT = Float64
    param_set = param_set_Float64

    # Test conditions
    T_test = FT(280)  # Temperature [K]
    ρ_test = FT(1.2)  # Density [kg/m³]
    p_test = FT(101325)  # Pressure [Pa]

    # Finite difference step size
    ε = sqrt(eps(FT))

    # Helper for AD vs FD
    function check_derivative(f, x; rtol = 5e-2, atol = 1e-8)
        ad = ForwardDiff.derivative(f, x)
        fd = (f(x + ε) - f(x - ε)) / (2ε)
        if abs(ad) < atol && abs(fd) < atol
            @test true
        else
            @test isapprox(ad, fd; rtol, atol)
        end
    end

    @testset "Basic thermodynamic derivatives" begin
        # Test dp/de_int at fixed (ρ, q_tot)
        q_tot_dry = FT(0.0)
        q_tot_moist = FT(0.01)

        for q_tot in [q_tot_dry, q_tot_moist]
            @testset "dp/de_int with q_tot=$q_tot" begin
                # Partition for consistency
                (q_liq, q_ice) =
                    TD.condensate_partition(param_set, T_test, ρ_test, q_tot)

                # Function: given e_int, compute p
                # Note: assumes q_liq and q_ice are fixed (frozen composition derivative)
                function p_from_e_int(e_int)
                    T = TD.air_temperature(param_set, e_int, q_tot, q_liq, q_ice)
                    return TD.air_pressure(
                        param_set,
                        T,
                        ρ_test,
                        q_tot,
                        q_liq,
                        q_ice,
                    )
                end

                # Reference e_int
                e_int_ref =
                    TD.internal_energy(param_set, T_test, q_tot, q_liq, q_ice)

                check_derivative(p_from_e_int, e_int_ref; rtol = 2e-3)
            end
        end
    end

    @testset "Saturation adjustment derivatives (ρe)" begin
        maxiter = 100
        tol = FT(1e-14)

        # Unsaturated case
        q_tot_unsat = FT(0.001)
        e_int_unsat = TD.internal_energy(param_set, T_test, q_tot_unsat)

        # Saturated case (Liquid)
        T_sat_liq = FT(290)
        q_vap_sat_liq = TD.q_vap_saturation(param_set, T_sat_liq, ρ_test)
        q_tot_sat_liq = q_vap_sat_liq * FT(1.5)
        (q_liq_sat_liq, q_ice_sat_liq) = TD.condensate_partition(
            param_set,
            T_sat_liq,
            ρ_test,
            q_tot_sat_liq,
        )
        e_int_sat_liq = TD.internal_energy(
            param_set,
            T_sat_liq,
            q_tot_sat_liq,
            q_liq_sat_liq,
            q_ice_sat_liq,
        )

        # Saturated case (Ice)
        T_sat_ice = FT(230)
        q_vap_sat_ice = TD.q_vap_saturation(param_set, T_sat_ice, ρ_test)
        q_tot_sat_ice = q_vap_sat_ice * FT(1.5)
        (q_liq_sat_ice, q_ice_sat_ice) = TD.condensate_partition(
            param_set,
            T_sat_ice,
            ρ_test,
            q_tot_sat_ice,
        )
        e_int_sat_ice = TD.internal_energy(
            param_set,
            T_sat_ice,
            q_tot_sat_ice,
            q_liq_sat_ice,
            q_ice_sat_ice,
        )

        for (name, e_int, q_tot) in [
            ("unsaturated", e_int_unsat, q_tot_unsat),
            ("saturated_liq", e_int_sat_liq, q_tot_sat_liq),
            ("saturated_ice", e_int_sat_ice, q_tot_sat_ice),
        ]
            @testset "$name (ρe)" begin
                # Closures for derivatives
                get_sol(ρ, e, q) = TD.saturation_adjustment(
                    RS.SecantMethod,
                    param_set,
                    ρe(),
                    ρ,
                    e,
                    q,
                    maxiter,
                    tol,
                )

                # Wrt e_int
                check_derivative(e -> get_sol(ρ_test, e, q_tot).T, e_int)
                check_derivative(e -> get_sol(ρ_test, e, q_tot).q_liq, e_int)
                check_derivative(e -> get_sol(ρ_test, e, q_tot).q_ice, e_int)

                # Wrt ρ
                check_derivative(ρ -> get_sol(ρ, e_int, q_tot).T, ρ_test)
                check_derivative(ρ -> get_sol(ρ, e_int, q_tot).q_liq, ρ_test)
                check_derivative(ρ -> get_sol(ρ, e_int, q_tot).q_ice, ρ_test)

                # Wrt q_tot
                check_derivative(q -> get_sol(ρ_test, e_int, q).T, q_tot)
                check_derivative(q -> get_sol(ρ_test, e_int, q).q_liq, q_tot)
                check_derivative(q -> get_sol(ρ_test, e_int, q).q_ice, q_tot)
            end
        end
    end

    @testset "Saturation adjustment derivatives (pe)" begin
        maxiter = 100
        tol = FT(1e-14)

        # Unsaturated
        q_tot_unsat = FT(0.001)
        e_int_unsat = TD.internal_energy(param_set, T_test, q_tot_unsat)

        # Saturated
        T_sat = FT(290)
        ρ_sat = TD.air_density(param_set, T_sat, p_test, q_tot_unsat)
        q_vap_sat = TD.q_vap_saturation(param_set, T_sat, ρ_sat)
        q_tot_sat = q_vap_sat * FT(1.5)
        (q_liq_sat, q_ice_sat) =
            TD.condensate_partition(param_set, T_sat, ρ_sat, q_tot_sat)
        e_int_sat = TD.internal_energy(
            param_set,
            T_sat,
            q_tot_sat,
            q_liq_sat,
            q_ice_sat,
        )

        for (name, e_int, q_tot) in [
            ("unsaturated", e_int_unsat, q_tot_unsat),
            ("saturated", e_int_sat, q_tot_sat),
        ]
            @testset "$name (pe)" begin
                get_sol(p, e, q) = TD.saturation_adjustment(
                    RS.SecantMethod,
                    param_set,
                    pe(),
                    p,
                    e,
                    q,
                    maxiter,
                    tol,
                )

                # Wrt e_int
                check_derivative(e -> get_sol(p_test, e, q_tot).T, e_int)
                check_derivative(e -> get_sol(p_test, e, q_tot).q_liq, e_int)
                check_derivative(e -> get_sol(p_test, e, q_tot).q_ice, e_int)

                # Wrt p
                check_derivative(p -> get_sol(p, e_int, q_tot).T, p_test)
                check_derivative(p -> get_sol(p, e_int, q_tot).q_liq, p_test)
                check_derivative(p -> get_sol(p, e_int, q_tot).q_ice, p_test)

                # Wrt q_tot
                check_derivative(q -> get_sol(p_test, e_int, q).T, q_tot)
                check_derivative(q -> get_sol(p_test, e_int, q).q_liq, q_tot)
                check_derivative(q -> get_sol(p_test, e_int, q).q_ice, q_tot)
            end
        end
    end

    @testset "Saturation adjustment derivatives (ph)" begin
        maxiter = 100
        tol = FT(1e-14)

        # Unsaturated
        q_tot_unsat = FT(0.005)
        h_unsat = TD.enthalpy(param_set, T_test, q_tot_unsat)

        # Saturated (re-using pe values, just conversion to h needed)
        # Note: Need self-consistent h.
        T_sat = FT(290)
        ρ_sat = TD.air_density(param_set, T_sat, p_test, q_tot_unsat)
        q_vap_sat = TD.q_vap_saturation(param_set, T_sat, ρ_sat)
        q_tot_sat = q_vap_sat * FT(1.5)
        (q_liq_sat, q_ice_sat) =
            TD.condensate_partition(param_set, T_sat, ρ_sat, q_tot_sat)
        h_sat =
            TD.enthalpy(param_set, T_sat, q_tot_sat, q_liq_sat, q_ice_sat)

        for (name, h, q_tot) in [
            ("unsaturated", h_unsat, q_tot_unsat),
            ("saturated", h_sat, q_tot_sat),
        ]
            @testset "$name (ph)" begin
                get_sol(p, h, q) = TD.saturation_adjustment(
                    RS.SecantMethod,
                    param_set,
                    ph(),
                    p,
                    h,
                    q,
                    maxiter,
                    tol,
                )

                # Wrt h
                check_derivative(h_val -> get_sol(p_test, h_val, q_tot).T, h)
                check_derivative(h_val -> get_sol(p_test, h_val, q_tot).q_liq, h)
                check_derivative(h_val -> get_sol(p_test, h_val, q_tot).q_ice, h)

                # Wrt p
                check_derivative(p -> get_sol(p, h, q_tot).T, p_test)
                check_derivative(p -> get_sol(p, h, q_tot).q_liq, p_test)
                check_derivative(p -> get_sol(p, h, q_tot).q_ice, p_test)

                # Wrt q_tot
                check_derivative(q -> get_sol(p_test, h, q).T, q_tot)
                check_derivative(q -> get_sol(p_test, h, q).q_liq, q_tot)
                check_derivative(q -> get_sol(p_test, h, q).q_ice, q_tot)
            end
        end
    end
end
