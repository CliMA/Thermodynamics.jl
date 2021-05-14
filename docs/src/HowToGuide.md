# How to guide

If Thermodynamics.jl does not have a particular thermodynamic
constructor that is needed, you can implement a new one in
`Thermodynamics/states.jl`. In this constructor, you must
add whichever arguments you wish to offer as inputs, then translate this
thermodynamic state into one of:

 - `PhaseDry` a dry thermodynamic state, uniquely determined by two
   independent thermodynamic properties
 - `PhaseEquil` a moist thermodynamic state in thermodynamic equilibrium,
   uniquely determined by three independent thermodynamic properties
 - `PhaseNonEquil` a moist thermodynamic state in thermodynamic
   non-equilibrium, uniquely determined by four independent thermodynamic
   properties
