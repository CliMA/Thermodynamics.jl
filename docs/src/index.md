# Thermodynamics

## Principal design

This package uses an abstraction that leverages the idea of a thermodynamic state:

 - given two (or more) independent intrinsic thermodynamic properties, we can establish a thermodynamic state and
 - given a thermodynamic state, we can compute any thermodynamic property

 
## Simple example

For our example, we first load packages, and create a set of thermodynamic parameters, using a convenience constructor, offered through [ClimaParams.jl](https://github.com/CliMA/ClimaParams.jl). Then we create a thermodynamic state using density, liquid-ice potential temperature, and total specific humidity. Finally, we compute air temperature from the thermodynamic state.

```@example
using ClimaParams # needed in environment to load convenience parameter struct wrappers
import Thermodynamics as TD
params = TD.Parameters.ThermodynamicsParameters(Float64);

ts = TD.PhaseEquil_ρθq(params, 1.0, 374.0, 0.01);
T = TD.air_temperature(params, ts)
```

And that's it. See a full list of different thermodynamic state constructors, in case you want to create a thermodynamic state with different variables, [here](https://clima.github.io/Thermodynamics.jl/dev/API/#Thermodynamic-State-Constructors).

See a full list of quantities that you can compute from a thermodynamic state, see our thermodynamic state-compatible methods [here](https://clima.github.io/Thermodynamics.jl/dev/API/#Thermodynamic-state-methods).

## How to guide

Thermodynamics.jl provides all thermodynamic functions needed for the
atmosphere and functions shared across model components. The functions are
general for a moist atmosphere that includes suspended cloud condensate in
the working fluid; the special case of a dry atmosphere is obtained for zero
specific humidities (or simply by omitting the optional specific humidity
arguments in the functions that are needed for a dry atmosphere). The
general formulation assumes that there are tracers for specific humidity
`q`, partitioned into

 - `q.tot` total water specific humidity
 - `q.liq` liquid specific humidity
 - `q.ice` ice specific humidity

to characterize the thermodynamic state and composition of moist air.

There are several types of functions:

1. Equation of state (ideal gas law):
    * `air_pressure`
2. Specific gas constant and isobaric and isochoric specific heats of moist
   air:
    * `gas_constant_air`
    * `cp_m`
    * `cv_m`
3. Specific latent heats of vaporization, fusion, and sublimation:
    * `latent_heat_vapor`
    * `latent_heat_fusion`
    * `latent_heat_sublim`
4. Saturation vapor pressure and specific humidity over liquid and ice:
    * `sat_vapor_press_liquid`
    * `sat_vapor_press_ice`
    * `sat_shum`
5. Functions computing energies and inverting them to obtain temperatures
    * `total_energy`
    * `internal_energy`
    * `air_temperature`
6. Functions to compute temperatures and partitioning of water into phases in
   thermodynamic equilibrium (when Gibbs' phase rule implies that the entire
   thermodynamic state of moist air, including the liquid and ice specific
   humidities, can be calculated from the 3 thermodynamic state variables, such
   as energy, pressure, and total specific humidity)
    * `liquid_fraction` (fraction of condensate that is liquid)
    * `saturation_adjustment` (compute temperature from energy, density, and
      total specific humidity)
7. Auxiliary functions for diagnostic purposes, e.g., other thermodynamic
quantities
    * `liquid_ice_pottemp` (liquid-ice potential temperature)

A moist dynamical core that assumes equilibrium thermodynamics can be
obtained from a dry dynamical core with total energy as a prognostic
variable by including a tracer for the total specific humidity `q.tot`,
using the functions, e.g., for the energies in the module, and computing
the temperature `T` and the liquid and ice specific humidities (`q.liq` and
`q.ice`) from the internal energy `e_int` by saturation adjustment.

## Dycore pseudo code

Here, we outline how users might use Thermodynamics inside a circulation model.
Users are encouraged to first establish a thermodynamic state with one of our
[Thermodynamic State Constructors](@ref). For example, we would construct
a moist thermodynamic state using

```julia
ts = PhaseEquil_ρeq(param_set, ρ, e_int, q_tot)
```

here, `ρ` is the density of the moist air, and the internal energy `e_int =
e_tot - e_kin - geopotential` is the total energy `e_tot` minus kinetic energy
`e_kin` and potential energy `geopotential` (all energies per unit mass). Once
we've established a thermodynamic state, we can call [Thermodynamic state
methods](@ref) that support thermodynamic states:

```julia
T = air_temperature(param_set, ts)
q = PhasePartition(param_set, ts)
```

No changes to the "right-hand sides" of the dynamical equations are needed
for a moist dynamical core that supports clouds, as long as they do not
precipitate. Additional source-sink terms arise from precipitation.

Schematically, the workflow in such a core would look as follows:
```julia
# initialize
geopotential = grav * z
q_tot        = ...
ρ            = ...

(u, v, w)    = ...
e_kin        = 0.5 * (u^2 + v^2 + w^2)

e_tot        = total_energy(param_set, e_kin, geopotential, T, q_tot)

do timestep  # timestepping loop

  # advance dynamical variables by a timestep (temperature typically
  # appears in terms on the rhs, such as radiative transfer)
  advance(u, v, w, ρ, e_tot, q_tot)

  # compute internal energy from dynamic variables
  e_int = e_tot - 0.5 * (u^2 + v^2 + w^2) - geopotential

  # compute temperature, pressure and condensate specific humidities,
  ts = PhaseEquil_ρeq(param_set, ρ, e_int, q_tot)
  T = air_temperature(param_set, ts)
  q = PhasePartition(param_set, ts)
  p = air_pressure(param_set, ts)

end
```

For a dynamical core that additionally uses the liquid and ice specific
humidities `q.liq` and `q.ice` as prognostic variables, and thus explicitly
allows the presence of non-equilibrium phases such as supercooled water,
the saturation adjustment in the above workflow is replaced calling a
non-equilibrium moist thermodynamic state:

```julia
q_tot, q_liq, q_ice = ...
ts = PhaseNonEquil(param_set, e_int, ρ, PhasePartition(q_tot, q_liq, q_ice))
T = air_temperature(param_set, ts)
p = air_pressure(param_set, ts)
```
