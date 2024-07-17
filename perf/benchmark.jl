include("common.jl")

@testset "Thermodynamics - Performance pθq constructor" begin

    trial = BenchmarkTools.@benchmark thermo_state_ρeq(
        $param_set,
        $ρ,
        $e_int,
        $q_tot,
    )
    show(stdout, MIME("text/plain"), trial)
    trial = BenchmarkTools.@benchmark thermo_state_pθq(
        $param_set,
        $p,
        $θ_liq_ice,
        $q_tot,
    )
    show(stdout, MIME("text/plain"), trial)
    trial =
        BenchmarkTools.@benchmark thermo_state_pTq($param_set, $p, $T, $q_tot)
    show(stdout, MIME("text/plain"), trial)

end
