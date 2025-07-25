"""
# Isentropic Processes Test Suite

This file contains tests for adiabatic processes and ideal gas law consistency.
"""

using Random

@testset "Thermodynamics - isentropic processes" begin
    for ArrayType in array_types
        FT = eltype(ArrayType)
        param_set = TP.ThermodynamicsParameters(FT)

        # Extract thermodynamic parameters using the common function
        (
            _R_d,
            _Rv_over_Rd,
            _R_v,
            _cp_d,
            _cp_v,
            _cp_l,
            _cp_i,
            _cv_d,
            _cv_v,
            _cv_l,
            _cv_i,
            _T_0,
            _e_int_v0,
            _e_int_i0,
            _LH_v0,
            _LH_s0,
            _LH_f0,
            _press_triple,
            _T_triple,
            _T_freeze,
            _T_icenuc,
            _T_min,
            _T_max,
            _p_ref_theta,
            _kappa_d,
        ) = extract_thermodynamic_parameters(param_set)

        profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
        (; T, p, e_int, ρ, θ_liq_ice, phase_type) = profiles
        (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = profiles

        # Test ideal gas law consistency across all profiles
        T_idgl = TD.air_temperature_given_ρp.(param_set, p, ρ, q_pt)
        @test all(T .≈ T_idgl)

        Φ = FT(1)
        Random.seed!(15)
        perturbation = FT(0.1) * rand(FT, length(T))

        # Test adiabatic processes with perturbed ambient conditions
        T∞, p∞ = T .* perturbation, p .* perturbation
        @test air_temperature.(param_set, p, θ_liq_ice, DryAdiabaticProcess()) ≈
              (p ./ _p_ref_theta) .^ (_R_d / _cp_d) .* θ_liq_ice
        @test TD.air_pressure_given_θ.(
            param_set,
            θ_liq_ice,
            Φ,
            DryAdiabaticProcess(),
        ) ≈ _p_ref_theta .* (1 .- Φ ./ (θ_liq_ice .* _cp_d)) .^ (_cp_d / _R_d)
        @test air_pressure.(param_set, T, T∞, p∞, DryAdiabaticProcess()) ≈
              p∞ .* (T ./ T∞) .^ (FT(1) / _kappa_d)
    end
end
