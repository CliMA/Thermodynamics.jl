<div align="center">
  <img src="docs/src/assets/logo.svg" alt="Thermodynamics.jl Logo" width="128" height="128">
</div>

# Thermodynamics.jl

A comprehensive Julia package for Earth system thermodynamics, providing consistent and accurate thermodynamic functions for moist air including all phases of water (vapor, liquid, and ice).

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
```julia
import Thermodynamics as TD
using ClimaParams

# Create thermodynamic parameters
params = TD.Parameters.ThermodynamicsParameters(Float64)

# Create a thermodynamic state
ρ = 1.0
e_int = -7.e4
q_tot = 0.01
ts = TD.PhaseEquil_ρeq(params, ρ, e_int, q_tot)

# Compute thermodynamic properties from state
T = TD.air_temperature(params, ts)
p = TD.air_pressure(params, ts)
q = TD.PhasePartition(params, ts)
```

## Key Features

### 🌟 **Comprehensive Thermodynamics**
- **Complete moist air thermodynamics** including all water phases (vapor, liquid, ice)
- **Consistent formulation** for use across all Earth system model components
- **Precipitation included** in the atmospheric working fluid for full thermodynamic consistency
- **Calorically perfect gas approximation** enabling closed-form expressions for saturation vapor pressure

### ⚡ **High Performance**
- **Type-stable implementations** for optimal Julia performance
- **GPU-compatible** implementations
- **Differentiable implementation** compatible with Julia's automatic differentiation capabilities
- **Efficient saturation adjustment** with Newton's method and analytical derivatives

### 🔧 **Flexible Design**
- **Multiple thermodynamic state constructors** for different use cases
- **Direct function access** for efficient single calculations
- **Equilibrium and non-equilibrium** phase partitioning
- **Extensible parameter system** for different planetary atmospheres

## Core Design Principles

### **Working Fluid Definition**
The working fluid includes **moist air with precipitation**, ensuring mass and energy conservation across all phases and thermodynamic consistency throughout the system.

### **Consistent Formulation**
All thermodynamic quantities are derived from a single fundamental approximation of **calorically perfect gases** with constant specific heat capacities, providing accuracy within 1-3% for atmospheric conditions.

### **Thermodynamic State Abstraction**
Optionally, given two (or more) independent intrinsic thermodynamic properties, we can establish a thermodynamic state from which any thermodynamic property can be computed.

## Documentation

- **[Mathematical Formulation](https://clima.github.io/Thermodynamics.jl/dev/Formulation/)** - Complete theoretical framework and equations
- **[API Reference](https://clima.github.io/Thermodynamics.jl/dev/API/)** - Complete function documentation
- **[How-To Guide](https://clima.github.io/Thermodynamics.jl/dev/HowToGuide/)** - Practical usage examples and patterns

## Integration with Climate Models

Thermodynamics.jl is used by several CliMA components:

- [ClimaAtmos](https://github.com/CliMA/ClimaAtmos.jl)
- [ClimaLand](https://github.com/CliMA/ClimaLand.jl)
- [ClimaOcean](https://github.com/CliMA/ClimaOcean.jl)
- [ClimaCoupler](https://github.com/CliMA/ClimaCoupler.jl)
- [CloudMicrophysics](https://github.com/CliMA/CloudMicrophysics.jl)
- [SurfaceFluxes](https://github.com/CliMA/SurfaceFluxes.jl)
- [KinematicDriver](https://github.com/CliMA/KinematicDriver.jl)

## Getting Help

For questions and issues, please check the [documentation](https://clima.github.io/Thermodynamics.jl/dev/) or open an issue on the [GitHub repository](https://github.com/CliMA/Thermodynamics.jl).
