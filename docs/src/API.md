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
vapor_specific_humidity
condensate_specific_humidity
vol_vapor_mixing_ratio
partial_pressure_dry
partial_pressure_vapor
vapor_pressure_deficit
relative_humidity
condensate_partition
has_condensate
liquid_fraction
liquid_fraction_ramp
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
latent_heat_fusion
latent_heat_sublim
latent_heat_vapor
latent_heat_mixed
humidity_weighted_latent_heat
```

### Air Temperatures

```@docs
air_temperature
potential_temperature
potential_temperature_given_pressure
liquid_ice_pottemp
liquid_ice_pottemp_given_pressure
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
```

### Air Saturation Functions

```@docs
saturation_vapor_pressure
q_vap_saturation
supersaturation
saturation_excess
q_vap_from_RH
q_vap_from_p_vap
q_vap_saturation_from_pressure
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
TemperatureProfiles.TemperatureProfile
TemperatureProfiles.DecayingTemperatureProfile
TemperatureProfiles.IsothermalProfile
TemperatureProfiles.DryAdiabaticProfile
```

## Data Collection

```@docs
Thermodynamics.DataCollection
```

## Deprecated Functions

These functions are deprecated and will be removed in a future release.

### Backward Compatibility Wrappers

These wrappers exist for backward compatibility with older versions of the package.

```@docs
specific_enthalpy
specific_enthalpy_dry
specific_enthalpy_vapor
specific_enthalpy_liquid
specific_enthalpy_ice
dry_pottemp
total_specific_enthalpy
q_vap_saturation_generic
latent_heat_liq_ice
```

### Other Deprecated Functions

```@docs
air_temperature_given_hq
air_temperature_given_pρq
air_temperature_given_pθq
air_temperature_given_ρθq
air_temperature_given_ρθq_nonlinear
saturated
total_specific_humidity
liquid_specific_humidity
ice_specific_humidity
mixing_ratios
specific_volume
q_vap_from_RH_liquid
temperature_and_humidity_given_TᵥρRH
liquid_ice_pottemp_sat
PhasePartition_equil
PhasePartition_equil_given_p
```

## Thermodynamic State Constructors (Deprecated)

```@docs
ThermodynamicState
PhasePartition
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
