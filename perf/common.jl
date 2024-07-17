using Test

import BenchmarkTools
import RootSolvers as RS

import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import ClimaParams as CP

const FT = Float64
const param_set = TP.ThermodynamicsParameters(FT)

function thermo_state_ρeq(param_set, _ρ, _e_int, _q_tot)
    for i in 1:length(_q_tot)
        (ρ, e_int, q_tot) = (_ρ[i], _e_int[i], _q_tot[i])
        ts = TD.PhaseEquil_ρeq(param_set, ρ, e_int, q_tot)
    end
    return nothing
end

function thermo_state_pθq(param_set, _p, _θ_liq_ice, _q_tot)
    for i in 1:length(_q_tot)
        (p, θ_liq_ice, q_tot) = (_p[i], _θ_liq_ice[i], _q_tot[i])
        ts = TD.PhaseEquil_pθq(param_set, p, θ_liq_ice, q_tot)
    end
    return nothing
end

function thermo_state_pTq(param_set, _p, _T, _q_tot)
    for i in 1:length(_q_tot)
        (p, T, q_tot) = (_p[i], _T[i], _q_tot[i])
        ts = TD.PhaseEquil_pTq(param_set, p, T, q_tot)
    end
    return nothing
end

ArrayType = Array{FT}
profiles = TD.TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
(; e_int, T, ρ, p, θ_liq_ice, q_tot) = profiles
