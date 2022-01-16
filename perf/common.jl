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

ArrayType = Array{Float64}
profiles = TD.TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
UnPack.@unpack e_int, T, ρ, p, θ_liq_ice, q_tot = profiles

function thermo_state_ρeq()
    ts = TD.PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot)
end

function thermo_state_pθq()
    ts = TD.PhaseEquil_pθq.(param_set, p, θ_liq_ice, q_tot)
end

function thermo_state_pTq()
    ts = TD.PhaseEquil_pTq.(param_set, p, T, q_tot)
end
