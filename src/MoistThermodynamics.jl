"""
    MoistThermodynamics

Moist thermodynamic functions, e.g., for air pressure (atmosphere equation
of state), latent heats of phase transitions, saturation vapor pressures, and
saturation specific humidities.
"""
module MoistThermodynamics

using DocStringExtensions

@inline q_pt_0(::Type{FT}) where {FT} = PhasePartition{FT}(FT(0), FT(0), FT(0))

using RootSolvers

include("Parameters.jl")
using .Parameters

include("states.jl")
include("isentropic.jl")
include("relations.jl")
include("test_moist_thermo.jl")

end
