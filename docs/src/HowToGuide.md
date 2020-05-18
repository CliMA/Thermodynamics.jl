# How-to-guide

## Usage

Users are encouraged to first establish a thermodynamic state with one of our [Thermodynamic State Constructors](@ref). For example, we would construct a moist thermodynamic state using

```julia
ts = PhaseEquil(param_set, e_int, ρ, q_tot);
```

here, `ρ` is the density of the moist air, and the internal energy `e_int = e_tot - e_kin - geopotential` is the total energy `e_tot` minus kinetic energy `e_kin` and potential energy `geopotential` (all energies per unit mass). Once we've established a thermodynamic state, we can call [Thermodynamic state methods](@ref) that support thermodynamic states:

```julia
T = air_temperature(ts);
q = PhasePartition(ts);
```

No changes to the "right-hand sides" of the dynamical equations are needed for a moist dynamical core that supports clouds, as long as they do not precipitate. Additional source-sink terms arise from precipitation.

Schematically, the workflow in such a core would look as follows:
```julia
# initialize
geopotential = grav * z
q_tot          = ...
ρ            = ...

(u, v, w)    = ...
e_kin           = 0.5 * (u^2 + v^2 + w^2)

e_tot        = total_energy(e_kin, geopotential, T, q_tot)

do timestep   # timestepping loop

  # advance dynamical variables by a timestep (temperature typically
  # appears in terms on the rhs, such as radiative transfer)
  advance(u, v, w, ρ, e_tot, q_tot)

  # compute internal energy from dynamic variables
  e_int = e_tot - 0.5 * (u^2 + v^2 + w^2) - geopotential

  # compute temperature, pressure and condensate specific humidities,
  ts = PhaseEquil(param_set, e_int, ρ, q_tot);
  T = air_temperature(ts);
  q = PhasePartition(ts);
  p = air_pressure(ts);

end
```

For a dynamical core that additionally uses the liquid and ice specific humidities `q.liq` and `q.ice` as prognostic variables, and thus explicitly allows the presence of non-equilibrium phases such as supercooled water, the saturation adjustment in the above workflow is replaced calling a non-equilibrium moist thermodynamic state:
```julia
q_tot, q_liq, q_ice = ...
ts = PhaseNonEquil(param_set, e_int, ρ, PhasePartition(q_tot, q_liq, q_ice));
T = air_temperature(ts);
p = air_pressure(ts);
```

## Extending

If Thermodynamics.jl does not have a particular thermodynamic constructor that is needed, one can implement a new one in `Thermodynamics/src/states.jl`. In this constructor, one must add whichever arguments they wish to offer as inputs, then translate this thermodynamic state into one of:

 - `PhaseDry` a dry thermodynamic state, uniquely determined by two independent thermodynamic properties
 - `PhaseEquil` a moist thermodynamic state in thermodynamic equilibrium, uniquely determined by three independent thermodynamic properties
 - `PhaseNonEquil` a moist thermodynamic state in thermodynamic non-equilibrium, uniquely determined by four independent thermodynamic properties

For example, to add a thermodynamic state constructor that accepts temperature, density and total specific humidity, we could add the following code to states:

```
"""
    TemperatureSHumEquil_given_density(param_set, T, ρ, q_tot)

Constructs a [`PhaseEquil`](@ref) thermodynamic state from temperature.

 - `param_set` parameter set, used to dispatch planet parameter function calls
 - `T` temperature
 - `ρ` density
 - `q_tot` total specific humidity
"""
function TemperatureSHumEquil(
    param_set::APS,
    T::FT,
    ρ::FT,
    q_tot::FT,
) where {FT <: Real}
    q = PhasePartition_equil(param_set, T, ρ, q_tot)
    e_int = internal_energy(param_set, T, q)
    return PhaseEquil{FT, typeof(param_set)}(param_set, e_int, ρ, q_tot, T)
end
```
