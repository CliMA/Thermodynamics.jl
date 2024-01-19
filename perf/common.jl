using Test

import BenchmarkTools
import RootSolvers as RS

import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import CLIMAParameters as CP

const FT = Float64
const param_set = TP.ThermodynamicsParameters(FT)

ArrayType = Array{FT}
profiles = TD.TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
(; e_int, T, ρ, p, θ_liq_ice, q_tot) = profiles

function thermo_state_ρeq()
    ts = TD.PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot)
end

function thermo_state_pθq()
    ts = TD.PhaseEquil_pθq.(param_set, p, θ_liq_ice, q_tot)
end

function thermo_state_pTq()
    ts = TD.PhaseEquil_pTq.(param_set, p, T, q_tot)
end
