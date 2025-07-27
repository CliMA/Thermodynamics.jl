"""
Tests for dry adiabatic processes.
"""

using Random

@testset "Thermodynamics - isentropic processes" begin
    for ArrayType in array_types
        FT = eltype(ArrayType)
        param_set = FT == Float64 ? param_set_Float64 : param_set_Float32
        _R_d = TP.R_d(param_set)
        _cp_d = TP.cp_d(param_set)
        _p_ref_theta = TP.p_ref_theta(param_set)
        _kappa_d = TP.kappa_d(param_set)

        profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
        (; T, p, ρ, θ_liq_ice, q_pt) = profiles

        @testset "Ideal Gas Law" begin
            T_idgl = TD.air_temperature_given_ρp.(param_set, p, ρ, q_pt)
            @test all(T .≈ T_idgl)
        end

        @testset "Dry Adiabatic Processes" begin
            Φ = FT(1)
            Random.seed!(15)
            perturbation = FT(0.1) * rand(FT, length(T))

            T∞, p∞ = T .* perturbation, p .* perturbation
            @test air_temperature.(
                param_set,
                p,
                θ_liq_ice,
                DryAdiabaticProcess(),
            ) ≈ (p ./ _p_ref_theta) .^ (_R_d / _cp_d) .* θ_liq_ice
            @test TD.air_pressure_given_θ.(
                param_set,
                θ_liq_ice,
                Φ,
                DryAdiabaticProcess(),
            ) ≈
                  _p_ref_theta .*
                  (1 .- Φ ./ (θ_liq_ice .* _cp_d)) .^ (_cp_d / _R_d)
            @test air_pressure.(param_set, T, T∞, p∞, DryAdiabaticProcess()) ≈
                  p∞ .* (T ./ T∞) .^ (FT(1) / _kappa_d)
        end
    end
end
