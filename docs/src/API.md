# Thermodynamics

```@meta
CurrentModule = Thermodynamics
```

```@docs
Thermodynamics
```

## Thermodynamics Parameters

```@docs
Parameters
Parameters.ThermodynamicsParameters
```

## Types

```@docs
Phase
Liquid
Ice
PhasePartition
IndepVars
ρe
pe
ph
pρ
pθ_li
ρθ_li
DryAdiabaticProcess
```

## Thermodynamic Functions

### Air Properties

```@docs
gas_constant_air
cp_m
cv_m
soundspeed_air
```

### Air Humidities

```@docs
liquid_specific_humidity
ice_specific_humidity
vapor_specific_humidity
condensate_specific_humidity
mixing_ratios
vol_vapor_mixing_ratio
partial_pressure_dry
partial_pressure_vapor
vapor_pressure_deficit
relative_humidity
condensate_partition
has_condensate
liquid_fraction
specific_humidity_to_mixing_ratio
```

### Air Energies

```@docs
internal_energy
internal_energy_dry
internal_energy_vapor
internal_energy_liquid
internal_energy_ice
enthalpy
enthalpy_dry
enthalpy_vapor
enthalpy_liquid
enthalpy_ice
total_energy
total_enthalpy
dry_static_energy
vapor_static_energy
moist_static_energy
virtual_dry_static_energy
internal_energy_sat
latent_heat_fusion
latent_heat_sublim
latent_heat_vapor
humidity_weighted_latent_heat
latent_heat_mixed
```

### Air Temperatures

```@docs
air_temperature
potential_temperature
potential_temperature_given_pressure
liquid_ice_pottemp
liquid_ice_pottemp_given_pressure
liquid_ice_pottemp_sat
virtual_pottemp
virtual_temperature
```

### Air Pressure and Density

```@docs
air_pressure
air_density
exner
exner_given_pressure
air_pressure_given_θ
specific_volume
```

### Air Saturation Functions

```@docs
saturation_vapor_pressure
q_vap_saturation
supersaturation
saturation_excess
q_vap_from_RH
q_vap_from_p_vap
```

### Saturation Adjustment

```@docs
saturation_adjustment
```

### Air Entropies

```@docs
entropy
entropy_dry
entropy_vapor
```

## Temperature Profiles

```@docs
TemperatureProfiles.IsothermalProfile
TemperatureProfiles.TemperatureProfile
TemperatureProfiles.DryAdiabaticProfile
TemperatureProfiles.DecayingTemperatureProfile
```

## Data Collection

```@docs
Thermodynamics.DataCollection
```

## Deprecated Functions

These functions are deprecated and will be removed in a future release.

```@docs
air_temperature_given_hq
air_temperature_given_pρq
air_temperature_given_pθq
air_temperature_given_ρθq
air_temperature_given_ρθq_nonlinear
saturated
total_specific_humidity
q_vap_from_RH_liquid
temperature_and_humidity_given_TᵥρRH
```

## Thermodynamic State Constructors (Deprecated)

```@docs
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

## Internal Methods

```@docs
q_vap_saturation_from_pressure
PhasePartition_equil
PhasePartition_equil_given_p
```
