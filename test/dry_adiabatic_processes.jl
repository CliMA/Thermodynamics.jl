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

        profiles = TestedProfiles.EquilMoistProfiles(param_set, ArrayType)
        (; T, p, ρ, θ_li, q_tot, q_liq, q_ice) = profiles

        @testset "Ideal Gas Law" begin
            T_idgl =
                TD.air_temperature.(
                    Ref(param_set),
                    Ref(TD.pρ()),
                    p,
                    ρ,
                    q_tot,
                    q_liq,
                    q_ice,
                )
            @test all(T .≈ T_idgl)
        end

        @testset "Dry Adiabatic Processes" begin
            Φ = FT(1)
            # A local RNG, so that this file does not perturb the global stream that later
            # test files draw from. The perturbation is centered on 1 (±5%): scaling by a
            # factor in [0, 0.1) as before would drive T∞ towards zero, where (T/T∞)^(1/κ_d)
            # overflows in Float32, and would exercise a physically meaningless regime.
            rng = MersenneTwister(15)
            perturbation = FT(1) .+ FT(0.1) .* (rand(rng, FT, length(T)) .- FT(0.5))

            T∞, p∞ = T .* perturbation, p .* perturbation
            @test air_temperature.(
                param_set,
                DryAdiabaticProcess(),
                p,
                θ_li,
            ) ≈ (p ./ _p_ref_theta) .^ (_R_d / _cp_d) .* θ_li
            @test TD.air_pressure_given_θ.(
                param_set,
                DryAdiabaticProcess(),
                θ_li,
                Φ,
            ) ≈
                  _p_ref_theta .*
                  (1 .- Φ ./ (θ_li .* _cp_d)) .^ (_cp_d / _R_d)
            @test air_pressure.(param_set, DryAdiabaticProcess(), T, T∞, p∞) ≈
                  p∞ .* (T ./ T∞) .^ (FT(1) / _kappa_d)
        end
    end
end
