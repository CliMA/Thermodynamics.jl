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
latent_heat
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
∂q_vap_sat_∂T
supersaturation
saturation_excess
q_vap_from_RH
q_vap_from_p_vap
q_vap_saturation_from_pressure
```

### Saturation Adjustment

```@docs
saturation_adjustment
∂e_int_∂T_sat_ρ
∂e_int_∂T_sat_p
∂θ_li_∂T_sat_ρ
∂θ_li_∂T_sat_p
∂p_∂T_sat_ρ
∂h_∂T_sat_p
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

