using Test
import Thermodynamics
const TD = Thermodynamics

import UnPack
import BenchmarkTools

import RootSolvers
const RS = RootSolvers

using TOML
planet_parse = TOML.parsefile(joinpath("test", "planet_parameters.toml"))
param_dict = Dict{Symbol, Float64}()
for (key, val) in planet_parse
    # In the future - we will use the full names,
    # param_dict[Symbol(key)] = val["value"]

    # for now we use the aliases
    param_dict[Symbol(val["alias"])] = val["value"]
end
full_parameter_set = (; param_dict...)

param_set = TD.ThermodynamicsParameters(full_parameter_set)

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
