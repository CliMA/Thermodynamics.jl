# Tested Profiles

## Overview

Thermodynamics.jl is tested using a comprehensive set of thermodynamic profiles specified in `src/TestedProfiles.jl`. These profiles represent a wide range of atmospheric conditions and are used to validate the thermodynamic calculations.

## Purpose

The tested profiles serve several important functions:

- **Validation**: Ensure thermodynamic calculations are correct across diverse conditions
- **Coverage**: Test the full range of atmospheric temperatures and humidities
- **Robustness**: Verify numerical stability under various thermodynamic states
- **Benchmarking**: Provide consistent test cases for performance evaluation

!!! note "Related Testing Tools"
    For additional testing scenarios, see [Temperature Profiles](TemperatureProfiles.md)
    for pre-defined atmospheric temperature profiles.

## Profile Characteristics

The tested profiles cover:

- **Altitude range**: 0 to 25 km (surface to upper troposphere)
- **Temperature range**: 150 K to 340 K (cryogenic to hot conditions)
- **Humidity range**: 0 to 102% relative humidity (subsaturated to supersaturated)
- **Density range**: Full atmospheric density variation with height

## Dry Phase Profiles

Dry phase profiles test thermodynamic calculations without moisture, ensuring the dry air thermodynamics are correctly implemented.

```@example
import Thermodynamics as TD
import Plots
import ClimaParams as CP
import Thermodynamics.Parameters as TP
FT = Float64
param_set = TP.ThermodynamicsParameters(FT)

profiles = TD.TestedProfiles.PhaseDryProfiles(param_set, Array{FT});
(;T, ρ, z) = profiles
p1 = Plots.scatter(ρ, z./10^3, xlabel="Density [kg/m^3]", ylabel="z [km]", title="Density");
p2 = Plots.scatter(T, z./10^3, xlabel="Temperature [K]", ylabel="z [km]", title="Temperature");
Plots.plot(p1, p2, layout=(1,2))
Plots.savefig("tested_profiles_dry.svg");
```
![](tested_profiles_dry.svg)

### Key Features
- **No moisture**: All humidity components are zero
- **Pressure variation**: Follows hydrostatic balance
- **Temperature variation**: Uses decaying temperature profile
- **Density calculation**: Computed from ideal gas law

## Moist Phase Profiles (Thermodynamic Equilibrium)

Moist phase profiles test thermodynamic calculations with moisture in thermodynamic equilibrium, including saturation adjustment.

```@example
import Thermodynamics as TD
import Plots
import ClimaParams as CP
import Thermodynamics.Parameters as TP
FT = Float64
param_set = TP.ThermodynamicsParameters(FT)

profiles = TD.TestedProfiles.PhaseEquilProfiles(param_set, Array{FT});
(;T, ρ, q_tot, z) = profiles
p1 = Plots.scatter(ρ, z./10^3, xlabel="Density [kg/m^3]", ylabel="z [km]", title="Density");
p2 = Plots.scatter(T, z./10^3, xlabel="Temperature [K]", ylabel="z [km]", title="Temperature");
p3 = Plots.scatter(q_tot*1000, z./10^3, xlabel="Total specific\nhumidity [g/kg]", ylabel="z [km]", title="Total specific\nhumidity");
Plots.plot(p1, p2, p3, layout=(1,3))
Plots.savefig("tested_profiles_virt_temp.svg")
```
![](tested_profiles_virt_temp.svg)

### Key Features
- **Moisture included**: Total specific humidity varies with height
- **Saturation adjustment**: Phase partitioning determined by equilibrium
- **Wide humidity range**: From subsaturated to supersaturated conditions

## Profile Generation

The profiles are generated using:

1. **Altitude grid**: 50 points from 0 to 25 km
2. **Temperature profile**: Decaying temperature from 340 K at surface to 150 K at 25 km
3. **Pressure**: Computed using hydrostatic balance
4. **Humidity**: Relative humidity from 0% to 102% (including supersaturation)

## Usage in Testing

These profiles are used in the test suite to validate:

- **State constructors**: All thermodynamic state creation methods
- **Property calculations**: Temperature, pressure, density, humidity
- **Energy calculations**: Internal energy, enthalpy, potential temperature
- **Phase transitions**: Saturation adjustment and phase partitioning
- **Numerical stability**: Convergence and accuracy across parameter space

## Integration with Development

The tested profiles are automatically used in:

- **Unit tests**: Validate individual function correctness
- **Integration tests**: Verify end-to-end thermodynamic calculations
- **Performance tests**: Benchmark computational efficiency

!!! tip "Development Workflow"
    When adding new thermodynamic functionality, ensure it works correctly
    with these tested profiles. The profiles provide comprehensive coverage
    of atmospheric conditions and help catch numerical issues early.
