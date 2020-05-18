# Thermodynamics

Thermodynamics.jl provides all thermodynamic functions needed for the atmosphere and functions shared across model components. The functions are general for a moist atmosphere that includes suspended cloud condensate in the working fluid; the special case of a dry atmosphere is obtained for zero specific humidities (or simply by omitting the optional specific humidity arguments in the functions that are needed for a dry atmosphere). The general formulation assumes that there are tracers for specific humidity `q`, partitioned into

 - `q.tot` total water specific humidity
 - `q.liq` liquid specific humidity
 - `q.ice` ice specific humidity

to characterize the thermodynamic state and composition of moist air.

There are several types of functions:

1. Equation of state (ideal gas law):
    * `air_pressure`
2. Specific gas constant and isobaric and isochoric specific heats of moist air:
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
6. Functions to compute temperatures and partitioning of water into phases in thermodynamic equilibrium (when Gibbs' phase rule implies that the entire thermodynamic state of moist air, including the liquid and ice specific humidities, can be calculated from the 3 thermodynamic state variables, such as energy, pressure, and total specific humidity)
    * `liquid_fraction` (fraction of condensate that is liquid)
    * `saturation_adjustment` (compute temperature from energy, density, and total specific humidity)
7. Auxiliary functions for diagnostic purposes, e.g., other thermodynamic quantities
    * `liquid_ice_pottemp` (liquid-ice potential temperature)

A moist dynamical core that assumes equilibrium thermodynamics can be obtained from a dry dynamical core with total energy as a prognostic variable by including a tracer for the total specific humidity `q.tot`, using the functions, e.g., for the energies in the module, and computing the temperature `T` and the liquid and ice specific humidities (`q.liq` and `q.ice`) from the internal energy `e_int` by saturation adjustment.
