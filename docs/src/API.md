# Thermodynamics

```@meta
CurrentModule = Thermodynamics
```
```@docs
Thermodynamics
```

## Thermodynamics Parameters

```@docs
Parameters.ThermodynamicsParameters
```

## Thermodynamic State Constructors

```@docs
PhasePartition
ThermodynamicState
PhaseDry
PhaseDry_ρe
PhaseDry_pT
PhaseDry_pθ
PhaseDry_pe
PhaseDry_ph
PhaseDry_ρθ
PhaseDry_ρT
PhaseDry_ρp
PhaseEquil
PhaseEquil_ρeq
PhaseEquil_ρTq
PhaseEquil_pTq
PhaseEquil_pθq
PhaseEquil_peq
PhaseEquil_phq
PhaseEquil_ρθq
PhaseEquil_ρpq
PhaseNonEquil
PhaseNonEquil_ρTq
PhaseNonEquil_pTq
PhaseNonEquil_ρθq
PhaseNonEquil_pθq
PhaseNonEquil_peq
PhaseNonEquil_phq
PhaseNonEquil_ρpq
```

## Types

```@docs
Ice
Liquid
```

## Thermodynamic Functions

```@docs
air_density
air_pressure
air_temperature
air_pressure_given_θ
air_temperature_given_ρp
condensate
condensate_specific_humidity
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
moist_static_energy
virtual_dry_static_energy
q_vap_from_RH_liquid
q_vap_saturation
q_vap_saturation_liquid
q_vap_saturation_ice
q_vap_from_p_vap
q_vap_saturation_from_density
partial_pressure_vapor
partial_pressure_dry
relative_humidity
saturated
saturation_adjustment
saturation_excess
saturation_vapor_pressure
soundspeed_air
specific_enthalpy
specific_enthalpy_dry
specific_enthalpy_vapor
specific_enthalpy_liquid
specific_enthalpy_ice
specific_volume
supersaturation
total_energy
total_specific_enthalpy
total_specific_humidity
vapor_pressure_deficit
vapor_specific_humidity
vol_vapor_mixing_ratio
virtual_pottemp
virtual_temperature
specific_entropy

```

## Internal Methods

```@docs
shum_to_mixing_ratio
q_vap_saturation_generic
q_vap_saturation_from_pressure
PhasePartition_equil
PhasePartition_equil_given_p
```

## Dispatch Types

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

## Tested Profiles

```@docs
Thermodynamics.TestedProfiles
```

## Data Collection

```@docs
Thermodynamics.DataCollection
```
