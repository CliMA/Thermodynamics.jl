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
    
    @testset "Basic thermodynamic derivatives" begin
        # Test dp/de_int at fixed (ρ, q_tot)
        q_tot_dry = FT(0.0)
        q_tot_moist = FT(0.01)
        
        for q_tot in [q_tot_dry, q_tot_moist]
            @testset "dp/de_int with q_tot=$q_tot" begin
                # Partition for consistency
                (q_liq, q_ice) = TD.condensate_partition(param_set, T_test, ρ_test, q_tot)

                # Function: given e_int, compute p
                # Note: assumes q_liq and q_ice are fixed (frozen composition derivative)
                function p_from_e_int(e_int)
                    T = TD.air_temperature(param_set, e_int, q_tot, q_liq, q_ice)
                    return TD.air_pressure(param_set, T, ρ_test, q_tot, q_liq, q_ice)
                end
                
                # Reference e_int
                e_int_ref = TD.internal_energy(param_set, T_test, q_tot, q_liq, q_ice)
                
                # ForwardDiff derivative
                dp_de_ad = ForwardDiff.derivative(p_from_e_int, e_int_ref)
                
                # Finite difference derivative
                dp_de_fd = (p_from_e_int(e_int_ref + ε) - p_from_e_int(e_int_ref - ε)) / (2ε)
                
                @test isapprox(dp_de_ad, dp_de_fd, rtol=2e-3)
            end
        end
    end
    
    @testset "Saturation adjustment derivatives (ρe)" begin
        maxiter = 50
        tol = FT(1e-6)
        
        # Unsaturated case: low q_tot
        q_tot_unsat = FT(0.001)
        e_int_unsat = TD.internal_energy(param_set, T_test, q_tot_unsat)
        
        # Saturated case (Liquid): high q_tot (at saturation)
        T_sat_liq = FT(290)
        q_vap_sat_liq = TD.q_vap_saturation(param_set, T_sat_liq, ρ_test)
        q_tot_sat_liq = q_vap_sat_liq * FT(1.5)  # Supersaturated
        # Compute equilibrium partition first to get consistent e_int
        (q_liq_sat_liq, q_ice_sat_liq) = TD.condensate_partition(param_set, T_sat_liq, ρ_test, q_tot_sat_liq)
        e_int_sat_liq = TD.internal_energy(param_set, T_sat_liq, q_tot_sat_liq, q_liq_sat_liq, q_ice_sat_liq)

        # Saturated case (Ice): cold temperature
        T_sat_ice = FT(230)
        q_vap_sat_ice = TD.q_vap_saturation(param_set, T_sat_ice, ρ_test)
        q_tot_sat_ice = q_vap_sat_ice * FT(1.5)
        # Compute equilibrium partition first
        (q_liq_sat_ice, q_ice_sat_ice) = TD.condensate_partition(param_set, T_sat_ice, ρ_test, q_tot_sat_ice)
        e_int_sat_ice = TD.internal_energy(param_set, T_sat_ice, q_tot_sat_ice, q_liq_sat_ice, q_ice_sat_ice)
        
        for (name, e_int, q_tot) in [
            ("unsaturated", e_int_unsat, q_tot_unsat),
            ("saturated_liq", e_int_sat_liq, q_tot_sat_liq),
            ("saturated_ice", e_int_sat_ice, q_tot_sat_ice),
        ]
            @testset "dq_liq/de_int and dq_ice/de_int ($name, ρe)" begin
                function q_liq_from_e_int(e_int)
                    sol = TD.saturation_adjustment(
                        RS.SecantMethod,
                        param_set,
                        ρe(),
                        ρ_test,
                        e_int,
                        q_tot,
                        maxiter,
                        tol,
                    )
                    return sol.q_liq
                end
                
                function q_ice_from_e_int(e_int)
                    sol = TD.saturation_adjustment(
                        RS.SecantMethod,
                        param_set,
                        ρe(),
                        ρ_test,
                        e_int,
                        q_tot,
                        maxiter,
                        tol,
                    )
                    return sol.q_ice
                end
                
                # ForwardDiff derivatives
                dq_liq_de_ad = ForwardDiff.derivative(q_liq_from_e_int, e_int)
                dq_ice_de_ad = ForwardDiff.derivative(q_ice_from_e_int, e_int)
                
                # Finite difference derivatives
                dq_liq_de_fd = (q_liq_from_e_int(e_int + ε) - q_liq_from_e_int(e_int - ε)) / (2ε)
                dq_ice_de_fd = (q_ice_from_e_int(e_int + ε) - q_ice_from_e_int(e_int - ε)) / (2ε)
                
                # Note: unsaturated case should have zero derivatives
                if name == "unsaturated"
                    @test abs(dq_liq_de_ad) < 1e-10
                    @test abs(dq_ice_de_ad) < 1e-10
                    @test abs(dq_liq_de_fd) < 1e-10
                    @test abs(dq_ice_de_fd) < 1e-10
                else
                    # Saturated case: compare AD to FD
                    @test isapprox(dq_liq_de_ad, dq_liq_de_fd, rtol=1e-2, atol=1e-10)
                    @test isapprox(dq_ice_de_ad, dq_ice_de_fd, rtol=1e-2, atol=1e-10)
                    
                    # Verify we actually have condensate in the reference solution
                    sol = TD.saturation_adjustment(RS.SecantMethod, param_set, ρe(), ρ_test, e_int, q_tot, maxiter, tol)
                    if name == "saturated_liq"
                        @test sol.q_liq > 0
                        @test sol.q_ice == 0
                    elseif name == "saturated_ice"
                        @test sol.q_ice > 0
                        @test sol.q_liq == 0
                    end
                end
            end
        end
    end
    
    @testset "Saturation adjustment derivatives (pe)" begin
        maxiter = 50
        tol = FT(1e-6)
        
        # Unsaturated case
        q_tot_unsat = FT(0.001)
        e_int_unsat = TD.internal_energy(param_set, T_test, q_tot_unsat)
        
        # Saturated case
        T_sat = FT(290)
        ρ_sat = TD.air_density(param_set, T_sat, p_test, q_tot_unsat) # approx rho
        q_vap_sat = TD.q_vap_saturation(param_set, T_sat, ρ_sat)
        q_tot_sat = q_vap_sat * FT(1.5)
        # Partition for correct e_int
        (q_liq_sat, q_ice_sat) = TD.condensate_partition(param_set, T_sat, ρ_sat, q_tot_sat)
        e_int_sat = TD.internal_energy(param_set, T_sat, q_tot_sat, q_liq_sat, q_ice_sat)
        
        for (name, e_int, q_tot) in [
            ("unsaturated", e_int_unsat, q_tot_unsat),
            ("saturated", e_int_sat, q_tot_sat),
        ]
            @testset "dT/de_int ($name, pe)" begin
                function T_from_e_int(e_int)
                    sol = TD.saturation_adjustment(
                        RS.SecantMethod,
                        param_set,
                        pe(),
                        p_test,
                        e_int,
                        q_tot,
                        maxiter,
                        tol,
                    )
                    return sol.T
                end
                
                # ForwardDiff derivative
                dT_de_ad = ForwardDiff.derivative(T_from_e_int, e_int)
                
                # Finite difference derivative
                dT_de_fd = (T_from_e_int(e_int + ε) - T_from_e_int(e_int - ε)) / (2ε)
                
                @test isapprox(dT_de_ad, dT_de_fd, rtol=2e-3)
            end
        end
    end
    
    @testset "Saturation adjustment derivatives (ph)" begin
        maxiter = 50
        tol = FT(1e-6)
        
        # Test with enthalpy
        q_tot = FT(0.005)
        h_test = TD.enthalpy(param_set, T_test, q_tot)
        
        @testset "dT/dh (ph)" begin
            function T_from_h(h)
                sol = TD.saturation_adjustment(
                    RS.SecantMethod,
                    param_set,
                    ph(),
                    p_test,
                    h,
                    q_tot,
                    maxiter,
                    tol,
                )
                return sol.T
            end
            
            # ForwardDiff derivative
            dT_dh_ad = ForwardDiff.derivative(T_from_h, h_test)
            
            # Finite difference derivative
            dT_dh_fd = (T_from_h(h_test + ε) - T_from_h(h_test - ε)) / (2ε)
            
            @test isapprox(dT_dh_ad, dT_dh_fd, rtol=1e-3)
        end
    end
end
