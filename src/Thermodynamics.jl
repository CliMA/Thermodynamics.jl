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
_Rv_over_Rd = TP.Rv_over_Rd(param_set)
```

Because these parameters are widely used throughout this module,
`param_set` is an argument for many Thermodynamics functions.
    
## Saturation adjustment

Saturation adjustment functions accept:
- `sat_adjust_method`: numerical method type (from RootSolvers.jl)
- A function returning the numerical method instance (e.g.,
    `sa_numerical_method_ρpq` returns an instance of the numerical
    method for the `ρ-p-q_tot` formulation)

Supported methods in RootSolvers.jl:
- `NewtonsMethod`: Newton method with analytic gradients
- `NewtonsMethodAD`: Newton method with autodiff
- `SecantMethod`: Secant method
- `RegulaFalsiMethod`: Regula-Falsi method
"""
module Thermodynamics

import DocStringExtensions
const DSE = DocStringExtensions

import RootSolvers
const RS = RootSolvers

import KernelAbstractions
const KA = KernelAbstractions

include("Parameters.jl")
import .Parameters
const TP = Parameters
const APS = TP.AbstractThermodynamicsParameters

include("ThermodynamicTypes.jl")

# For printing literal strings on the gpu
include("printing.jl")

# Allow users to print warning and throw errors on non-convergence
# by importing:
# ```julia
# import Thermodynamics
# Thermodynamics.error_on_non_convergence() = true
# Thermodynamics.print_warning() = true
# ```
# By default, we don't print warnings and don't throw errors on non-convergence.
@inline error_on_non_convergence() = false
@inline print_warning() = false

# Default phase partition for dry air
@inline function q_pt_0(ps::APS)
    FT = eltype(ps)
    return PhasePartition(FT(0), FT(0), FT(0))
end

@inline solution_type() = RS.CompactSolution()
include("DataCollection.jl")
import .DataCollection

include("aux_functions.jl")
include("air_states.jl")

include("air_properties.jl")
include("air_phase_partition.jl")
include("air_humidities.jl")
include("air_energies.jl")
include("air_temperatures.jl")
include("air_pressure_and_density.jl")
include("air_saturation_functions.jl")
include("saturation_adjustment.jl")
include("air_entropies.jl")
include("air_dry_adiabatic.jl")
include("config_numerical_method.jl")
include("TemperatureProfiles.jl")
include("TestedProfiles.jl")

# State methods
include("air_state_methods.jl")

Base.broadcastable(dap::DryAdiabaticProcess) = tuple(dap)
Base.broadcastable(phase::Phase) = tuple(phase)

end #module Thermodynamics.jl
