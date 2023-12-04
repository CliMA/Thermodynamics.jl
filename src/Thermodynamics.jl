"""
    Thermodynamics

Moist thermodynamic functions, e.g., for air pressure (atmosphere equation
of state), latent heats of phase transitions, saturation vapor pressures, and
saturation specific humidities.


## AbstractParameterSet's

Many functions defined in this module rely on ClimaParams.jl.
ClimaParams.jl defines several functions (e.g., many planet
parameters). For example, to compute the mole-mass ratio:

```julia
import ClimaParams as CP
import Thermodynamics.Parameters as TP
FT = Float64
param_set = TP.ThermodynamicsParameters(FT)
_molmass_ratio = TP.molmass_ratio(param_set)
```

Because these parameters are widely used throughout this module,
`param_set` is an argument for many Thermodynamics functions.

## Numerical methods for saturation adjustment

Saturation adjustment function are designed to accept
 - `sat_adjust_method` a type used to dispatch which numerical method to use

and a function to return an instance of the numerical method. For example:

 - `sa_numerical_method_ρpq` returns an instance of the numerical
    method. One of these functions must be defined for the particular
    numerical method and the particular formulation (`ρ-p-q_tot` in this case).

The currently supported numerical methods, in RootSolvers.jl, are:
 - `NewtonsMethod` uses Newton method with analytic gradients
 - `NewtonsMethodAD` uses Newton method with autodiff
 - `SecantMethod` uses Secant method
 - `RegulaFalsiMethod` uses Regula-Falsi method
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

include("states.jl")
include("relations.jl")
include("isentropic.jl")
include("config_numerical_method.jl")
include("TemperatureProfiles.jl")
include("TestedProfiles.jl")

Base.broadcastable(dap::DryAdiabaticProcess) = tuple(dap)
Base.broadcastable(phase::Phase) = tuple(phase)

# For backwards compatibility with package extensions
if !isdefined(Base, :get_extension)
    include(joinpath("..", "ext", "CreateParametersExt.jl"))
end

end #module Thermodynamics.jl
