using Test
using Thermodynamics
import Thermodynamics.Parameters as TP

@testset "liquid_fraction_ramp" begin
    FT = Float64
    param_set = param_set_Float64

    T_freeze = TP.T_freeze(param_set)
    T_icenuc = TP.T_icenuc(param_set)

    # Test strictly above freezing
    @test Thermodynamics.liquid_fraction_ramp(param_set, T_freeze + 1) == FT(1)

    # Test at freezing
    # Based on the implementation: T > T_freeze -> 1. T <= T_freeze -> calculation.
    # If T = T_freeze, it falls into supercooled liquid branch (or effectively 1 if calculation yields 1).
    # lambda = ((T - Ti)/(Tf - Ti))^n. If T=Tf, lambda = 1^n = 1.
    @test Thermodynamics.liquid_fraction_ramp(param_set, T_freeze) ≈ FT(1)

    # Test strictly below nucleation
    @test Thermodynamics.liquid_fraction_ramp(param_set, T_icenuc - 1) == FT(0)

    # Test at nucleation
    # T > Ti check. If T = Ti, fails T > Ti. So returns 0.
    @test Thermodynamics.liquid_fraction_ramp(param_set, T_icenuc) == FT(0)

    # Test intermediate value
    T_mid = (T_freeze + T_icenuc) / 2
    λ = Thermodynamics.liquid_fraction_ramp(param_set, T_mid)
    @test 0 < λ < 1

    # Manual calculation check
    n = TP.pow_icenuc(param_set)
    expected_λ = ((T_mid - T_icenuc) / (T_freeze - T_icenuc))^n
    @test λ ≈ expected_λ
end
