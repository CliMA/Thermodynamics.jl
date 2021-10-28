# Thermodynamics

```@meta
CurrentModule = Thermodynamics
```
```@docs
Thermodynamics
```

## Thermodynamic State Constructors

```@docs
PhasePartition
ThermodynamicState
PhaseDry
PhaseDry_pT
PhaseDry_pθ
PhaseDry_pe
PhaseDry_ρθ
PhaseDry_ρT
PhaseDry_ρp
PhaseEquil
PhaseEquil_ρeq
PhaseEquil_ρTq
PhaseEquil_pTq
PhaseEquil_pθq
PhaseEquil_peq
PhaseEquil_ρθq
PhaseEquil_ρpq
PhaseNonEquil
PhaseNonEquil_ρTq
PhaseNonEquil_pTq
PhaseNonEquil_ρθq
PhaseNonEquil_pθq
PhaseNonEquil_peq
PhaseNonEquil_ρpq
```

## Types

```@docs
Ice
Liquid
```

## Thermodynamic state methods

```@docs
air_density
air_pressure
air_temperature
condensate
cp_m
cv_m
dry_pottemp
exner
gas_constant_air
gas_constants
has_condensate
ice_specific_humidity
internal_energy
internal_energy_dry
internal_energy_vapor
internal_energy_liquid
internal_energy_ice
internal_energy_sat
latent_heat_fusion
latent_heat_liq_ice
latent_heat_sublim
latent_heat_vapor
liquid_fraction
liquid_ice_pottemp
liquid_ice_pottemp_sat
liquid_specific_humidity
mixing_ratios
vol_vapor_mixing_ratio
moist_static_energy
q_vap_saturation
q_vap_saturation_liquid
q_vap_saturation_ice
relative_humidity
saturated
saturation_adjustment
saturation_excess
saturation_vapor_pressure
soundspeed_air
specific_enthalpy
specific_volume
supersaturation
total_energy
total_specific_enthalpy
total_specific_humidity
vapor_specific_humidity
virtual_pottemp
virtual_temperature
specific_entropy
```

## Internal methods

```@docs
shum_to_mixing_ratio
q_vap_saturation_generic
q_vap_saturation_from_pressure
PhasePartition_equil
PhasePartition_equil_given_p
```

## Dispatch types

```@docs
DryAdiabaticProcess
```

## Temperature Profiles
```@docs
TemperatureProfiles.IsothermalProfile
TemperatureProfiles.TemperatureProfile
TemperatureProfiles.DryAdiabaticProfile
TemperatureProfiles.DecayingTemperatureProfile
```

## Tested profiles

```@docs
Thermodynamics.TestedProfiles
```
