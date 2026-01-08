<div align="center">
  <img src="docs/src/assets/logo.svg" alt="Thermodynamics.jl Logo" width="128" height="128">
</div>

# Thermodynamics.jl

The `Thermodynamics.jl` package implements the thermodynamic formulation of the [CliMA Earth System Model](https://clima.caltech.edu) (Yatunin et al., 2026). It provides a consistent, accurate, and efficient framework for moist thermodynamics based on the **Rankine-Kirchhoff approximations** (Romps, 2021), providing consistent and accurate thermodynamic functions for moist air including all phases of water (vapor, liquid, and ice).

|||
|-----------------------------:|:-------------------------------------------------|
| **Documentation**            | [![dev][docs-latest-img]][docs-latest-url]       |
| **Docs Build**               | [![docs build][docs-bld-img]][docs-bld-url]      |
| **GHA CI**                   | [![gha ci][gha-ci-img]][gha-ci-url]              |
| **Code Coverage**            | [![codecov][codecov-img]][codecov-url]           |
| **Downloads**                | [![Downloads][dlt-img]][dlt-url]                 |

[docs-bld-img]: https://github.com/CliMA/Thermodynamics.jl/actions/workflows/docs.yml/badge.svg
[docs-bld-url]: https://github.com/CliMA/Thermodynamics.jl/actions/workflows/docs.yml

[docs-latest-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-latest-url]: https://CliMA.github.io/Thermodynamics.jl/dev/

[gha-ci-img]: https://github.com/CliMA/Thermodynamics.jl/actions/workflows/ci.yml/badge.svg
[gha-ci-url]: https://github.com/CliMA/Thermodynamics.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/CliMA/Thermodynamics.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/Thermodynamics.jl

[dlt-img]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FThermodynamics&query=total_requests&label=Downloads
[dlt-url]: https://juliapkgstats.com/pkg/Thermodynamics

## Quick Start

### Installation

```julia
using Pkg
Pkg.add("Thermodynamics")
Pkg.add("ClimaParams")
```

### Basic Usage

Thermodynamics.jl provides a **functional, stateless API**. You import the package (`TD`) and pass a **parameter set** plus **thermodynamic variables** (e.g., density, internal energy, specific humidities) directly to functions.

```julia
import Thermodynamics as TD
# Use RootSolvers for the saturation adjustment method
import RootSolvers as RS
using ClimaParams

# 1. Create thermodynamic parameters
#    (requires a definition of the parameter set, e.g. from ClimaParams)
params = TD.Parameters.ThermodynamicsParameters(Float64)

# 2. Define your thermodynamic variables
ρ     = 1.1        # Density [kg/m³]
e_int = 200000.0   # Internal energy [J/kg]
q_tot = 0.015      # Total specific humidity [kg/kg]
q_liq = 0.005      # Liquid specific humidity [kg/kg]
q_ice = 0.001      # Ice specific humidity [kg/kg]

# 3. Compute properties directly
T = TD.air_temperature(params, e_int, q_tot, q_liq, q_ice)
p = TD.air_pressure(params, T, ρ, q_tot, q_liq, q_ice)
```

### Saturation Adjustment

To find the equilibrium temperature and phase partition from non-equilibrium variables (e.g., given `ρ`, `e_int`, `q_tot`), use `saturation_adjustment`:

```julia
# Solve for equilibrium (T, q_liq, q_ice) given (ρ, e_int, q_tot)
# using SecantMethod
sol = TD.saturation_adjustment(
    RS.SecantMethod,        # Root-solving method
    params,                 # Parameter set
    TD.ρe(),                # Formulation: Density & Internal Energy
    ρ, e_int, q_tot,        # Input variables
    10,                     # Max iterations
    1e-3                    # Relative tolerance
)

println("Equilibrium T: ", sol.T)
println("Liquid q: ",      sol.q_liq)
println("Ice q: ",         sol.q_ice)
println("Converged: ",     sol.converged)
```

## Key Features

### 🌟 **Comprehensive Thermodynamics**

- **Complete moist air thermodynamics** including all water phases (vapor, liquid, ice).
- **Stateless, functional API** for maximum flexibility and ease of integration.
- **Consistent formulation** assuming a **calorically perfect gas** mixture.

### ⚡ **High Performance**

- **Type-stable** and **GPU-compatible** (CUDA.jl, AMDGPU.jl, etc.).
- **AD-compatible** (ForwardDiff.jl, etc.) for differentiable physics.
- **Zero-allocation** design for core functions.

### 🔧 **Flexible Design**

- **Multiple formulations**: Solve for equilibrium from `(ρ, e_int)`, `(p, h)`, `(p, θ_liq_ice)`, etc.
- **Extensible parameters**: Easily adapt to different planetary atmospheres via `ClimaParams`.

## Core Design Principles

### **Functional & Stateless**

Functions in Thermodynamics.jl are pure and stateless. They take a `ThermodynamicsParameters` struct and the necessary thermodynamic variables (e.g., `T`, `ρ`, `q`...) as arguments. This design fits naturally into large-scale simulations (like `ClimaAtmos.jl`) where prognostic variables are managed externally.

### **Working Fluid**

The working fluid is **moist air** (dry air + water vapor + liquid water + ice, which may include precipitation). We treat it as a mixture of ideal gases and condensed phases, ensuring rigorous mass and energy conservation.

### **Consistent Formulation**

All quantities are derived from the **calorically perfect gas** assumption with constant specific heat capacities. This provides a consistent, closed set of equations for saturation vapor pressures (the so-called Rankine-Kirchhoff approximation), latent heats, and other derived quantities.

## Documentation

- **[Mathematical Formulation](https://clima.github.io/Thermodynamics.jl/dev/Formulation/)** - Theoretical background.
- **[API Reference](https://clima.github.io/Thermodynamics.jl/dev/API/)** - Detailed function documentation.
- **[How-To Guide](https://clima.github.io/Thermodynamics.jl/dev/HowToGuide/)** - Recipes and examples.

## Integration with Climate Models

Thermodynamics.jl is the thermodynamic core for the [CliMA](https://github.com/CliMA) ecosystem, including:

- [ClimaAtmos](https://github.com/CliMA/ClimaAtmos.jl)
- [ClimaLand](https://github.com/CliMA/ClimaLand.jl)
- [ClimaOcean](https://github.com/CliMA/ClimaOcean.jl)
- [ClimaCoupler](https://github.com/CliMA/ClimaCoupler.jl)
- [CloudMicrophysics](https://github.com/CliMA/CloudMicrophysics.jl)
- [SurfaceFluxes](https://github.com/CliMA/SurfaceFluxes.jl)
- [KinematicDriver](https://github.com/CliMA/KinematicDriver.jl)

## Getting Help

For questions, verify the [documentation](https://clima.github.io/Thermodynamics.jl/dev/) or open an issue on [GitHub](https://github.com/CliMA/Thermodynamics.jl).
