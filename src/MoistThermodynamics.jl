"""
    MoistThermodynamics

Moist thermodynamic functions, e.g., for air pressure (atmosphere equation
of state), latent heats of phase transitions, saturation vapor pressures, and
saturation specific humidities.
"""
module MoistThermodynamics

using DocStringExtensions
using RootSolvers
using CLIMAParameters
using CLIMAParameters.Planet

@inline q_pt_0(::Type{FT}) where {FT} = PhasePartition{FT}(FT(0), FT(0), FT(0))

include("states.jl")
include("isentropic.jl")
include("relations.jl")

end
