using Test
import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import ClimaParams as CP
import RootSolvers as RS

@testset "Thermodynamics - exceptions/convergence" begin
    FT = Float64
    toml_dict = CP.create_toml_dict(FT)
    param_set = TP.ThermodynamicsParameters(toml_dict)

    # Construct a state that is definitely saturated and requires iteration
    ρ = FT(1.1)
    T_true = FT(300)
    q_tot = FT(0.025)
    e_int = TD.internal_energy_sat(param_set, T_true, ρ, q_tot)

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
        1e-6,
    )

    @test res.converged == false
    @test isfinite(res.T)

    # Verify that with enough iterations it does converge
    res_converged = TD.saturation_adjustment(
        RS.SecantMethod,
        param_set,
        TD.ρe(),
        ρ,
        e_int,
        q_tot,
        20,
        1e-6,
    )
    @test res_converged.converged == true
    @test isapprox(res_converged.T, T_true; rtol = 1e-5)
end
