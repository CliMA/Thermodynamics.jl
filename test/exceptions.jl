using Test
import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import ClimaParams as CP
import RootSolvers as RS

@testset "Thermodynamics - exceptions/convergence" begin
    for FT in (Float32, Float64)
        toml_dict = CP.create_toml_dict(FT)
        param_set = TP.ThermodynamicsParameters(toml_dict)

        # Construct a state that is definitely saturated and requires iteration
        ρ = FT(1.1)
        T_true = FT(300)
        q_tot = FT(0.025)
        e_int = TD.internal_energy_sat(param_set, T_true, ρ, q_tot)

        @testset "Non-convergence is reported, not thrown ($FT)" begin
            # Test that providing insufficient maxiter results in converged = false
            # instead of throwing an error.
            res = TD.saturation_adjustment(
                RS.SecantMethod,
                param_set,
                TD.ρe(),
                ρ,
                e_int,
                q_tot,
                1, # maxiter too small
                FT(1e-6),
            )

            @test res.converged == false
            @test isfinite(res.T)
            # Even unconverged, the returned state must be usable: a finite temperature
            # and a non-negative phase partition.
            @test res.q_liq >= 0
            @test res.q_ice >= 0

            # Verify that with enough iterations it does converge
            res_converged = TD.saturation_adjustment(
                RS.SecantMethod,
                param_set,
                TD.ρe(),
                ρ,
                e_int,
                q_tot,
                20,
                FT(1e-6),
            )
            @test res_converged.converged == true
            @test isapprox(res_converged.T, T_true; rtol = FT(1e-4))
        end

        @testset "Extreme inputs stay finite ($FT)" begin
            # States far outside the atmospheric range should still return something
            # usable rather than NaN, Inf, or an exception. The fixed-iteration solvers in
            # particular have no bracketing to fall back on.
            for (T_x, p_x, q_x) in (
                (FT(180), FT(1000), FT(1e-6)),     # cold and very thin
                (FT(320), FT(101325), FT(0.05)),   # hot and very moist
                (FT(240), FT(20000), FT(0.0)),     # completely dry
            )
                ρ_x = TD.air_density(param_set, T_x, p_x, q_x)
                e_x = TD.internal_energy_sat(param_set, T_x, ρ_x, q_x)
                sol = TD.saturation_adjustment(param_set, TD.ρe(), ρ_x, e_x, q_x)
                @test isfinite(sol.T)
                @test sol.T > 0
                @test isfinite(sol.q_liq) && sol.q_liq >= 0
                @test isfinite(sol.q_ice) && sol.q_ice >= 0
                @test sol.q_liq + sol.q_ice <= q_x + sqrt(eps(FT))
            end
        end
    end

    @testset "Unsupported formulations are a MethodError" begin
        # The formulation argument is dispatched on, so an unsupported one must fail
        # loudly at the call site rather than silently selecting a fallback.
        param_set = TP.ThermodynamicsParameters(Float64)
        @test_throws MethodError TD.saturation_adjustment(
            param_set,
            :not_a_formulation,
            1.0,
            50000.0,
            0.01,
        )
        @test_throws MethodError TD.air_temperature(
            param_set,
            :not_a_formulation,
            50000.0,
            0.01,
        )
    end
end
