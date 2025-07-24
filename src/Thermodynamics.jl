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

# For printing literal strings on the gpu
include("printing.jl")

# Allow users to skip error on non-convergence
# by importing:
# ```julia
# import Thermodynamics
# Thermodynamics.error_on_non_convergence() = false
# ```
# Error on convergence must be the default
# behavior because this can result in printing
# very large logs resulting in CI to seemingly hang.
@inline error_on_non_convergence() = true

# Allow users to skip printing warnings on non-convergence
@inline print_warning() = true

@inline q_pt_0(::Type{FT}) where {FT} = PhasePartition(FT(0), FT(0), FT(0))

@inline solution_type() = RS.CompactSolution()
include("DataCollection.jl")
import .DataCollection

include("aux_functions.jl")
include("states.jl")

include("phase_partition.jl")
include("material_properties_air.jl")
include("energies_and_enthalpies.jl")
include("humidities.jl")
include("temperatures.jl")
include("pressure_and_density.jl")
include("saturation_functions.jl")
include("saturation_adjustment.jl")
include("entropies.jl")
include("dry_adiabatic_relations.jl")
include("config_numerical_method.jl")
include("TemperatureProfiles.jl")
include("TestedProfiles.jl")

Base.broadcastable(dap::DryAdiabaticProcess) = tuple(dap)
Base.broadcastable(phase::Phase) = tuple(phase)

end #module Thermodynamics.jl
