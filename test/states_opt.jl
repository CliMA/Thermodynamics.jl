if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using Test
import Thermodynamics
const TD = Thermodynamics

import UnPack
import BenchmarkTools

import RootSolvers
const RS = RootSolvers

import CLIMAParameters
const CP = CLIMAParameters

struct EarthParameterSet <: CP.AbstractEarthParameterSet end
const param_set = EarthParameterSet()

@testset "Thermodynamics - Performance ρeq constructor" begin
    ArrayType = Array{Float64}
    FT = eltype(ArrayType)
    profiles = TD.TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
    UnPack.@unpack e_int, ρ, q_tot = profiles

    function thermo_state_ρeq_newton()
        ts =
            TD.PhaseEquil_dev_only.(
                param_set,
                ρ,
                e_int,
                q_tot;
                sat_adjust_method = RS.NewtonsMethod,
            )
    end

    function thermo_state_ρeq_regula_falsi()
        ts =
            TD.PhaseEquil_dev_only.(
                param_set,
                ρ,
                e_int,
                q_tot;
                sat_adjust_method = RS.RegulaFalsiMethod, maxiter = 20,
            )
    end

    function thermo_state_ρeq_secant()
        ts =
            TD.PhaseEquil_dev_only.(
                param_set,
                ρ,
                e_int,
                q_tot;
                sat_adjust_method = RS.SecantMethod, maxiter = 50,
            )
    end

    trial = BenchmarkTools.@benchmark $thermo_state_ρeq_newton()
    show(stdout, MIME("text/plain"), trial)
    trial = BenchmarkTools.@benchmark $thermo_state_ρeq_regula_falsi()
    show(stdout, MIME("text/plain"), trial)
    trial = BenchmarkTools.@benchmark $thermo_state_ρeq_secant()
    show(stdout, MIME("text/plain"), trial)

end

@testset "Thermodynamics - Performance pθq constructor" begin
    ArrayType = Array{Float64}
    profiles = TD.TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
    UnPack.@unpack p, θ_liq_ice, q_tot = profiles

    function thermo_state_pθq()
        ts = TD.PhaseEquil_pθq.(param_set, p, θ_liq_ice, q_tot)
    end
    trial = BenchmarkTools.@benchmark $thermo_state_pθq()
    show(stdout, MIME("text/plain"), trial)
end
