"""
    Thermodynamics

Moist thermodynamic functions for atmospheric modeling, including air pressure,
latent heats, saturation vapor pressures, and saturation specific humidities.

## Parameter Sets

Many functions defined in this module rely on ClimaParams.jl.
ClimaParams.jl defines several functions (e.g., many planet
parameters). For example, to compute the mole-mass ratio of 
dry air and water vapor:

```julia
import ClimaParams as CP
import Thermodynamics.Parameters as TP
FT = Float64
param_set = TP.ThermodynamicsParameters(FT)
Rv_over_Rd = TP.Rv_over_Rd(param_set)
```

Because these parameters are widely used throughout this module,
`param_set` is an argument for many Thermodynamics functions.

## Saturation adjustment

One entry point is [`saturation_adjustment`](@ref). It accepts:
- `method`: A root-solving method type from RootSolvers.jl (e.g., `RS.NewtonsMethod`, `RS.SecantMethod`).
- `param_set`: The thermodynamics parameter set.
- `formulation`: The thermodynamic formulation (e.g., [`ρe()`](@ref), [`ph()`](@ref)).
- Thermodynamic state variables appropriate for the formulation.

Supported methods in RootSolvers.jl:
- `NewtonsMethod`: Newton method with analytic gradients (recommended for `ρe`).
- `NewtonsMethodAD`: Newton method with automatic differentiation.
- `SecantMethod`: Secant method (derivative-free).
- `BrentsMethod`: Brent's method (hybrid root-finding).

## Thermodynamic functions

Many thermodynamic functions (e.g., [`air_temperature`](@ref), [`air_pressure`](@ref)) are
dispatched over the independent variables used to compute them. The available
independent variable types (subtypes of [`IndepVars`](@ref)) are:

- [`ρe`](@ref): Density, internal energy, and specific humidities
- [`pe`](@ref): Pressure, internal energy, and specific humidities
- [`ph`](@ref): Pressure, enthalpy, and specific humidities
- [`pρ`](@ref): Pressure, density, and specific humidities
- [`pθ_li`](@ref): Pressure, liquid-ice potential temperature, and specific humidities
- [`ρθ_li`](@ref): Density, liquid-ice potential temperature, and specific humidities
"""
module Thermodynamics

import RootSolvers as RS

include("Parameters.jl")
import .Parameters
const TP = Parameters
const APS = TP.AbstractThermodynamicsParameters

include("ThermoTypes.jl")

include("depr_PhasePartitionTypes.jl")

@inline solution_type() = RS.CompactSolution()
include("DataCollection.jl")
import .DataCollection

include("aux_functions.jl")


include("air_properties.jl")
include("air_humidities.jl")
include("air_energies.jl")
include("air_temperatures.jl")
include("air_pressure_and_density.jl")
include("air_saturation_functions.jl")
include("saturation_adjustment.jl")
include("air_entropies.jl")
include("air_dry_adiabatic.jl")
include("TemperatureProfiles.jl")

# Soon to be removed
include("depr_air_temperatures.jl")
include("depr_saturation_adjustment.jl")
include("depr_air_states.jl")
include("depr_state_methods.jl")
include("depr_phase_partition_methods.jl")

Base.broadcastable(dap::DryAdiabaticProcess) = tuple(dap)
Base.broadcastable(phase::Phase) = tuple(phase)

end # module Thermodynamics
