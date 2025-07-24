# Thermodynamics.jl

A comprehensive Julia package for atmospheric thermodynamics, providing consistent and accurate thermodynamic functions for moist air including all phases of water (vapor, liquid, and ice).

## Table of Contents

1. [Quick Start](#quick-start)
2. [Documentation Overview](#documentation-overview)
3. [Key Features](#key-features)
4. [Core Design Principles](#core-design-principles)
5. [Getting Started](#getting-started)
6. [Usage Examples](#usage-examples)
7. [Integration with Climate Models](#integration-with-climate-models)

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

# Or compute directly from independent variables
q_liq = 0.005
q_ice = 0.0004
q = TD.PhasePartition(q_tot, q_liq, q_ice)
T = TD.air_temperature(params, e_int, q)  # From internal energy and humidity
p = TD.air_pressure(params, ρ, T, q)      # From density, temperature, and humidity
```

## Documentation Overview

### 📚 **Core Documentation**

- **[Mathematical Formulation](Formulation.md)** - Complete theoretical framework and equations
  - Fundamental assumptions and working fluid definition
  - Equation of state and heat capacities
  - Internal energies, enthalpies, and latent heats
  - Saturation vapor pressure and specific humidity
  - Saturation adjustment algorithms
  - Auxiliary thermodynamic functions

- **[API Reference](API.md)** - Complete function documentation
  - Thermodynamic state constructors
  - Equation of state functions
  - Energy and temperature functions
  - Saturation and phase equilibrium functions
  - Auxiliary diagnostic functions

### 🛠️ **User Guides**

- **[How-To Guide](HowToGuide.md)** - Practical usage examples and patterns
  - Installation and setup
  - Common use cases and workflows

### 🔬 **Advanced Topics**

- **[Temperature Profiles](TemperatureProfiles.md)** - Pre-defined atmospheric profiles to be used as reference states in atmosphere models and for testing

- **[Tested Profiles](TestedProfiles.md)** - Thermodynamic profiles used for testing of the package

- **[Clausius-Clapeyron Validation](Clausius_Clapeyron.md)** - Validation of analytical derivatives

### 👨‍💻 **Developer Resources**

- **[Saturation Adjustment Convergence](SaturationAdjustmentConvergence.md)** - Convergence testing for numerical methods

### 📚 **Published References and Background**

- **[References](References.md)** - Bibliography of theoretical foundations

## Key Features

### 🌟 **Comprehensive Thermodynamics**
- **Complete moist air thermodynamics** including all water phases (vapor, liquid, ice)
- **Consistent formulation** for use across all model components
- **Precipitation included** in the atmospheric working fluid for full thermodynamic consistency
- **Calorically perfect gas approximation** enabling closed-form (Rankine-Kirchhoff) expressions for saturation vapor pressure

### ⚡ **High Performance**
- **Type-stable implementations** for optimal Julia performance
- **GPU-compatible** implementations
- **Efficient saturation adjustment** with Newton's method and analytical derivatives for equilibrium thermodynamics formulations

### 🔧 **Flexible Design**
- **Multiple thermodynamic state constructors** for different use cases
- **Direct function access** for efficient single calculations
- **Equilibrium and non-equilibrium** phase partitioning
- **Extensible parameter system** for different planetary atmospheres
- **Comprehensive testing** and validation suite

## Core Design Principles

### **Thermodynamic State Abstraction**
The package leverages the fundamental principle that:
- Given two (or more) independent intrinsic thermodynamic properties, we can establish a thermodynamic state
- Given a thermodynamic state, we can compute any thermodynamic property

This abstraction provides a clean, consistent interface for all thermodynamic calculations.

### **Working Fluid Definition**
The working fluid includes **moist air with precipitation**, ensuring:
- Mass and energy conservation across all phases
- Thermodynamic consistency throughout the system
- Unified treatment of cloud and precipitation condensate

### **Consistent Formulation**
All thermodynamic quantities are derived from a single fundamental approximation:
- **Calorically perfect gases** with constant specific heat capacities
- **Closed-form expressions** for all thermodynamic quantities
- **Accuracy within 1-3%** for atmospheric conditions
- **Computational efficiency** without numerical integration or additional ad-hoc approximations

## Getting Started

### Installation
```julia
using Pkg
Pkg.add("Thermodynamics")
Pkg.add("ClimaParams")
```

## Usage Examples

### **Equilibrium Thermodynamics (Saturation Adjustment)**
```julia
import Thermodynamics as TD 
using ClimaParams

params = TD.Parameters.ThermodynamicsParameters(Float64)

# Create state with internal energy, density, and total humidity
ρ = 1.0
e_int = -7.0e4
q_tot = 0.01
ts = TD.PhaseEquil_ρeq(params, ρ, e_int, q_tot)

# Temperature and phase partitioning computed automatically
T = TD.air_temperature(params, ts)
q = TD.PhasePartition(params, ts)
```

### **Non-Equilibrium Thermodynamics**
```julia
import Thermodynamics as TD
using ClimaParams

params = TD.Parameters.ThermodynamicsParameters(Float64)

# Explicit phase partitioning
q_tot = 0.01
q_liq = 0.005
q_ice = 0.0003
q = TD.PhasePartition(q_tot, q_liq, q_ice)
ρ = 1.0
e_int = -7.0e4
ts = TD.PhaseNonEquil(params, e_int, ρ, q)

# Temperature computation from thermodynamic state
T = TD.air_temperature(params, ts)

# Alternative direct computation, avoiding the thermodynamic state
T = TD.air_temperature(params, e_int, q)
```

### **Saturation Calculations**
```julia
import Thermodynamics as TD
using ClimaParams

params = TD.Parameters.ThermodynamicsParameters(Float64)

# Create a thermodynamic state for testing
ρ = 1.0
e_int = -7.0e4
q_tot = 0.01
ts = TD.PhaseEquil_ρeq(params, ρ, e_int, q_tot)
T = TD.air_temperature(params, ts)
p = TD.air_pressure(params, ts)

# Saturation vapor pressure
p_v_sat = TD.saturation_vapor_pressure(params, T, TD.Liquid())

# Saturation specific humidity
q_v_sat = TD.q_vap_saturation(params, T, ρ, typeof(ts))

# Relative humidity
RH = TD.relative_humidity(params, ts)

# Alternative, avoiding the thermodynamic state
q_liq = 0.005
q_ice = 0.0003
RH_alt = TD.relative_humidity(params, T, p, typeof(ts), TD.PhasePartition(q_tot, q_liq, q_ice))
```

## Integration with Climate Models

### **Dynamical Core Integration**
The package is designed for seamless integration with atmospheric dynamical cores, schematically as follows:

```julia
# Initialize
import Thermodynamics as TD
using ClimaParams

FT = Float64
params = TD.Parameters.ThermodynamicsParameters(FT)

# Define physical constants and initial conditions
grav = 9.81  # m/s²
cv_d = 718.0  # J/(kg K)
z = 1000.0   # m
geopotential = grav * z
q_tot = 0.01
ρ = 1.0
T = 300.0

# Initial velocity components
u, v, w = 10.0, 5.0, 1.0
e_kin = 0.5 * (u^2 + v^2 + w^2)
e_tot = TD.total_energy(params, e_kin, geopotential, T, TD.PhasePartition(q_tot))

# Timestepping loop (simplified example)
for timestep in 1:10
    # Advance dynamical variables (simplified)
    u += FT(0.1)
    v += FT(-0.1)
    w += FT(0.01)
    e_tot += cv_d * FT(1)
    
    # Compute internal energy
    e_kin = 0.5 * (u^2 + v^2 + w^2)
    e_int = e_tot - e_kin - geopotential
    
    # Saturation adjustment
    ts = TD.PhaseEquil_ρeq(params, ρ, e_int, q_tot)
    T = TD.air_temperature(params, ts)
    q = TD.PhasePartition(params, ts)
    p = TD.air_pressure(params, ts)
end
```

## Next Steps

1. **Read the [Mathematical Formulation](Formulation.md)** for theoretical background
2. **Explore the [API Reference](API.md)** for complete function documentation
3. **Follow the [How-To Guide](HowToGuide.md)** for practical examples
4. **Check [Saturation Adjustment Convergence](SaturationAdjustmentConvergence.md)** for numerical method testing

---

!!! note "Citation"
    If you use Thermodynamics.jl in your research, please cite the relevant papers listed in the [References](References.md) section.

!!! tip "Getting Help"
    For questions and issues, please check the documentation or open an issue on the [GitHub repository](https://github.com/CliMA/Thermodynamics.jl).
