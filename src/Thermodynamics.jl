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

Saturation adjustment functions accept:
- `sat_adjust_method`: numerical method type (from RootSolvers.jl)
- A function returning the numerical method instance (e.g.,
    `sa_numerical_method_ρpq` returns an instance of the numerical
    method for the `ρ-p-q_tot` formulation)

Supported methods in RootSolvers.jl:
- `NewtonsMethod`: Newton method with analytic gradients
- `NewtonsMethodAD`: Newton method with autodiff
- `SecantMethod`: Secant method
- `BrentsMethod`: Brent's method (hybrid root-finding)
"""
module Thermodynamics

import RootSolvers as RS
import KernelAbstractions as KA

include("Parameters.jl")
import .Parameters
const TP = Parameters
const APS = TP.AbstractThermodynamicsParameters

include("ThermoTypes.jl")

include("depr_PhasePartitionTypes.jl")

# For printing literal strings on the gpu
include("printing.jl")

# Allow users to print warnings and throw errors on non-convergence.
#
# Recommended configuration API:
# ```julia
# import Thermodynamics
# Thermodynamics.set_error_on_non_convergence!(true)
# Thermodynamics.set_print_warning!(true)
# ```
#
# (Advanced) Users can still override these methods directly:
# ```julia
# import Thermodynamics
# Thermodynamics.error_on_non_convergence() = true
# Thermodynamics.print_warning() = true
# ```
#
# By default, we don't print warnings and don't throw errors on non-convergence.
const _error_on_non_convergence_ref = Ref(false)
const _print_warning_ref = Ref(false)

@inline error_on_non_convergence() = false # _error_on_non_convergence_ref[]
@inline print_warning() = false # _print_warning_ref[]

set_error_on_non_convergence!(x::Bool) = (_error_on_non_convergence_ref[] = x; nothing)
set_print_warning!(x::Bool) = (_print_warning_ref[] = x; nothing)

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
include("config_sa_method.jl")
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
