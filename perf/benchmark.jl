include("common.jl")

@testset "Thermodynamics - Performance pθq constructor" begin

    trial = BenchmarkTools.@benchmark $thermo_state_ρeq()
    show(stdout, MIME("text/plain"), trial)
    trial = BenchmarkTools.@benchmark $thermo_state_pθq()
    show(stdout, MIME("text/plain"), trial)
    trial = BenchmarkTools.@benchmark $thermo_state_pTq()
    show(stdout, MIME("text/plain"), trial)

end
