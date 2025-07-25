import Random

"""
# Data-driven Test Suite

This file contains data-driven tests for constructor robustness and specific known test cases.
"""

"""
    sample_range(; param_set, e_int_range, ρ_range, q_tot_range, n_samples)

Generate random thermodynamic states within specified ranges to test constructor robustness.
This helps identify potential convergence issues with the saturation adjustment algorithm.

# Arguments
- `param_set`: Thermodynamic parameter set
- `e_int_range`: Range for internal energy values
- `ρ_range`: Range for density values  
- `q_tot_range`: Range for total specific humidity values
- `n_samples`: Number of random samples to generate

# Returns
- Nothing, but creates thermodynamic states to test constructor behavior
"""
function sample_range(; param_set, e_int_range, ρ_range, q_tot_range, n_samples)
    for i in 1:n_samples
        e_int = Random.rand(e_int_range)
        ρ = Random.rand(ρ_range)
        q_tot = Random.rand(q_tot_range)
        ts = PhaseEquil_ρeq(param_set, ρ, e_int, q_tot, 4)
    end
end

@testset "Thermodynamics - Data-driven tests" begin
    for ArrayType in array_types
        FT = eltype(ArrayType)
        param_set = FT == Float64 ? param_set_Float64 : param_set_Float32
        Random.seed!(1234)

        @testset "Constructor robustness - $FT" begin
            # Realistic atmospheric ranges for internal energy, density, and humidity
            e_int_range = (FT(28311.801716):FT(30981.514836))  # Internal energy range (J/kg)
            ρ_range = (FT(1.124755):FT(1.129586))              # Density range (kg/m³)
            q_tot_range = (FT(0.011897):FT(0.013305))          # Total humidity range (kg/kg)
            n_samples = 1_000                                  # Reduced from 11k for faster CI
            sample_range(;
                param_set,
                e_int_range,
                ρ_range,
                q_tot_range,
                n_samples,
            )
        end

        @testset "pθq data-driven tests - $FT" begin
            # Test cases that should work correctly - realistic atmospheric conditions
            pθq = [
                (;
                    p = 82307.30319719888,
                    θ_liq_ice = 296.7074342326652,
                    q_tot = 0.01894019026929829,
                ),
                (;
                    p = 88357.42589002676,
                    θ_liq_ice = 296.7074342326652,
                    q_tot = 0.01894019026929829,
                ),
                (;
                    p = 81090.35731696963,
                    θ_liq_ice = 296.7074342326652,
                    q_tot = 0.01894019026929829,
                ),
                (;
                    p = 75461.41343839701,
                    θ_liq_ice = 296.7074342326652,
                    q_tot = 0.01894019026929829,
                ),
                (;
                    p = 71158.23329080557,
                    θ_liq_ice = 296.7074342326652,
                    q_tot = 0.01894019026929829,
                ),
                (;
                    p = 68032.73180723468,
                    θ_liq_ice = 296.7074342326652,
                    q_tot = 0.01894019026929829,
                ),
                (;
                    p = 59906.16044902247,
                    θ_liq_ice = 296.7074342326652,
                    q_tot = 0.01894019026929829,
                ),
                (;
                    p = 51740.77281055945,
                    θ_liq_ice = 296.7074342326652,
                    q_tot = 0.01894019026929829,
                ),
                (;
                    p = 50056.7894012541,
                    θ_liq_ice = 296.7074342326652,
                    q_tot = 0.01894019026929829,
                ),
                (;
                    p = 48806.70567178993,
                    θ_liq_ice = 296.7074342326652,
                    q_tot = 0.01894019026929829,
                ),
                (;
                    p = 27564.73889538213,
                    θ_liq_ice = 296.7074342326652,
                    q_tot = 0.01894019026929829,
                ),
            ]

            config = ()
            for case in pθq
                @test begin
                    ts = PhaseEquil_pθq(
                        param_set,
                        FT(case.p),
                        FT(case.θ_liq_ice),
                        FT(case.q_tot),
                        config...,
                    )
                    air_pressure(param_set, ts) == FT(case.p)
                end
            end
        end
    end
end
