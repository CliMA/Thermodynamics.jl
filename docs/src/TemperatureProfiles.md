# Atmospheric Temperature Profiles

```@meta
CurrentModule = Thermodynamics.TemperatureProfiles
```

## Overview

The `Thermodynamics.TemperatureProfiles` module provides pre-defined atmospheric temperature profiles for testing and validation purposes, and to be used as reference states in atmosphere models. These profiles represent idealized atmospheric conditions commonly used in climate modeling and thermodynamic testing.

## Available Profiles

Three temperature profiles are currently available:

1. **IsothermalProfile**: Constant temperature with height
2. **DecayingTemperatureProfile**: Temperature decreases smoothly with height
3. **DryAdiabaticProfile**: Temperature follows dry adiabatic lapse rate

## Usage

Each profile constructor requires two arguments:

- `param_set`: A thermodynamic parameter set from [ClimaParams.jl](https://github.com/CliMA/ClimaParams.jl)
- `z`: Altitude (height above surface)

The profiles return both temperature and pressure as a function of altitude.

### Isothermal Profile

A constant temperature profile, useful for testing thermodynamic calculations under uniform conditions.

```@example
import Thermodynamics as TD
import Plots
import ClimaParams as CP
import Thermodynamics.Parameters as TP
FT = Float64
param_set = TP.ThermodynamicsParameters(FT)
z = range(FT(0), stop = FT(2.5e4), length = 50);

isothermal = TD.TemperatureProfiles.IsothermalProfile(param_set, FT);
args = isothermal.(Ref(param_set), z);
T = first.(args);
p = last.(args);

p1 = Plots.plot(T, z./10^3, xlabel="Temperature [K]");
p2 = Plots.plot(p./10^3, z./10^3, xlabel="Pressure [kPa]");
Plots.plot(p1, p2, layout=(1,2), title="Isothermal", ylabel="z [km]");
Plots.savefig("isothermal.svg");
```

![](isothermal.svg)

### Decaying Temperature Profile

A profile where temperature decreases smoothly with height, representing typical atmospheric conditions.

```@example
import Thermodynamics as TD
import Plots
import ClimaParams as CP
import Thermodynamics.Parameters as TP
FT = Float64
param_set = TP.ThermodynamicsParameters(FT)
z = range(FT(0), stop = FT(2.5e4), length = 50);

decaying = TD.TemperatureProfiles.DecayingTemperatureProfile{FT}(param_set);
args = decaying.(Ref(param_set), z);
T = first.(args);
p = last.(args);

p1 = Plots.plot(T, z./10^3, xlabel="Temperature [K]");
p2 = Plots.plot(p./10^3, z./10^3, xlabel="Pressure [kPa]");
Plots.plot(p1, p2, layout=(1,2), ylabel="z [km]", title="Decaying");
Plots.savefig("decaying.svg")
```

![](decaying.svg)

### Dry Adiabatic Profile

A profile following the dry adiabatic lapse rate, where potential temperature is constant with height.

```@example
import Thermodynamics as TD
import Plots
import ClimaParams as CP
import Thermodynamics.Parameters as TP
FT = Float64
param_set = TP.ThermodynamicsParameters(FT)
z = range(FT(0), stop = FT(2.5e4), length = 50);

dry_adiabatic = TD.TemperatureProfiles.DryAdiabaticProfile{FT}(param_set);
args = dry_adiabatic.(Ref(param_set), z)
T = first.(args)
p = last.(args)
θ_dry = TD.potential_temperature_given_pressure.(Ref(param_set), T, p)

p1 = Plots.plot(T, z./10^3, xlabel="Temperature [K]");
p2 = Plots.plot(p./10^3, z./10^3, xlabel="Pressure [kPa]");
p3 = Plots.plot(θ_dry, z./10^3, xlabel="Potential temperature [K]");
Plots.plot(p1, p2, p3, layout=(1,3), ylabel="z [km]", title="Dry adiabatic")
Plots.savefig("dry_adiabatic.svg")
```

![](dry_adiabatic.svg)

## Applications

These temperature profiles are useful for:

- **Testing thermodynamic calculations** under controlled conditions
- **Benchmarking performance** across different atmospheric conditions
- **Educational purposes** to understand atmospheric thermodynamics

## Integration with Thermodynamic States

Temperature profiles can be combined with thermodynamic state calculations:

```julia
# Example: Create thermodynamic states along a temperature profile
z = range(0, 25e3, length=100)  # 0 to 25 km
profile = TD.TemperatureProfiles.DryAdiabaticProfile{Float64}(param_set)

# Get temperature and pressure at each altitude
T_pairs = profile.(Ref(param_set), z)
T = first.(T_pairs)
p = last.(T_pairs)

# Compute density from temperature and pressure using the functional API
# For dry air: ρ = p / (R_d * T)
R_d = 287.0  # Specific gas constant for dry air [J/(kg·K)]
ρ = [p[i] / (R_d * T[i]) for i in 1:length(z)]

# Now compute other properties as needed, e.g., potential temperature
θ = [TD.potential_temperature_given_pressure(param_set, T[i], p[i]) for i in 1:length(z)]
```

## Extending the Module

Additional temperature profiles can be added by:

1. **Creating a struct** containing the parameters needed for the profile
2. **Implementing a functor** that takes `(param_set, z)` and returns `(T, p)`

Example structure:

```julia
struct CustomProfile{FT} <: TemperatureProfile{FT}
    # Profile parameters
end

function (profile::CustomProfile{FT})(param_set, z) where {FT}
    # Calculate T and p based on z
    return (T, p)
end
```

!!! note "Implementation Details"
    All profiles implement the `TemperatureProfile` interface and return
    temperature-pressure pairs as a function of altitude. The profiles are
    designed to be composable with other thermodynamic calculations.
