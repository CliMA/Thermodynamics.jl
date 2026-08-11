using Test
using Thermodynamics
import Thermodynamics.Parameters as TP

# Run in both precisions: the ramp raises a base to `pow_icenuc` and compares against the
# freezing and nucleation temperatures, both of which are more fragile in Float32.
@testset "liquid_fraction_ramp" begin
    for FT in (Float32, Float64)
        param_set = FT == Float64 ? param_set_Float64 : param_set_Float32

        T_freeze = TP.T_freeze(param_set)
        T_icenuc = TP.T_icenuc(param_set)

        @testset "Endpoints and interior ($FT)" begin
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

        @testset "Far outside the ramp ($FT)" begin
            # The ramp argument is clamped before exponentiation, so temperatures well
            # outside the interval must stay finite and in range rather than throwing.
            for T in FT.((1, 50, 150, 200, 320, 400))
                λ = Thermodynamics.liquid_fraction_ramp(param_set, T)
                @test isfinite(λ)
                @test 0 <= λ <= 1
            end
            @test Thermodynamics.liquid_fraction_ramp(param_set, FT(400)) == FT(1)
            @test Thermodynamics.liquid_fraction_ramp(param_set, FT(1)) == FT(0)
        end

        @testset "Derivative is consistent with the ramp ($FT)" begin
            # ∂λ/∂T vanishes outside the ramp and matches finite differences inside it.
            @test Thermodynamics.∂λ_∂T_ramp(param_set, T_freeze + 10) == FT(0)
            @test Thermodynamics.∂λ_∂T_ramp(param_set, T_icenuc - 10) == FT(0)

            δ = FT(0.05)
            for T in FT.((240, 250, 260, 270))
                fd =
                    (
                        Thermodynamics.liquid_fraction_ramp(param_set, T + δ) -
                        Thermodynamics.liquid_fraction_ramp(param_set, T - δ)
                    ) / (2δ)
                @test Thermodynamics.∂λ_∂T_ramp(param_set, T) ≈ fd rtol = FT(1e-2)
            end
        end
    end
end
