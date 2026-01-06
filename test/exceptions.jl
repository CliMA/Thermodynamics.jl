"""
# Exceptions Test Suite (functional API)

These tests ensure the functional `saturation_adjustment` API can surface
non-convergence via `Thermodynamics.error_on_non_convergence()`.
"""

@testset "Thermodynamics - Exceptions on Failed Convergence (functional)" begin
    for FT in (Float32, Float64)
        param_set = FT == Float64 ? param_set_Float64 : param_set_Float32

        # Configure to throw on non-convergence, and ensure we restore defaults afterwards.
        Thermodynamics.error_on_non_convergence() = true
        Thermodynamics.print_warning() = true
        try
            ρ = FT(1.0)
            T = FT(290)
            # Force a saturated case (harder for tiny maxiter)
            q_tot = FT(1.2) * TD.q_vap_saturation(param_set, T, ρ)
            e_int_sat = TD.internal_energy_sat(param_set, T, ρ, q_tot)

            @test_throws ErrorException TD.saturation_adjustment(
                RS.NewtonsMethod,
                param_set,
                TD.ρeq(),
                ρ,
                e_int_sat,
                q_tot,
                1, # maxiter
                FT(1e-12),
            )

            @test_throws ErrorException TD.saturation_adjustment(
                RS.SecantMethod,
                param_set,
                TD.ρeq(),
                ρ,
                e_int_sat,
                q_tot,
                1, # maxiter
                FT(1e-12),
            )
        finally
            Thermodynamics.error_on_non_convergence() = false
            Thermodynamics.print_warning() = false
        end
    end
end
