# Thermodynamics.jl

A comprehensive Julia package for atmospheric thermodynamics, providing consistent and accurate thermodynamic functions for moist air including all phases of water (vapor, liquid, and ice). `Thermodynamics.jl` implements the thermodynamic formulation of the CliMA Earth System Model ([Yatunin2026](@cite)). It provides a consistent, accurate, and efficient framework for moist thermodynamics based on the **Rankine-Kirchhoff approximations** ([Romps2021](@cite)).

## Table of Contents

1. [Quick Start](#quick-start)
2. [Documentation Overview](#documentation-overview)
3. [Key Features](#key-features)
4. [Core Design Principles](#core-design-principles)
5. [Usage Examples](#usage-examples)
6. [Integration with Climate Models](#integration-with-climate-models)

## Quick Start

### Installation

```julia
using Pkg
Pkg.add("Thermodynamics")
Pkg.add("ClimaParams")
```

### Basic Usage

```@example Index
import Thermodynamics as TD
# Use RootSolvers for the saturation adjustment method
import RootSolvers as RS
using ClimaParams

# 1. Create thermodynamic parameters
#    (requires a definition of the parameter set, e.g. from ClimaParams)
params = TD.Parameters.ThermodynamicsParameters(Float64)

# 2. Define your thermodynamic variables
ρ     = 1.1        # Density [kg/m³]
e_int = 130000.0   # Internal energy [J/kg] (approx. 290 K)
q_tot = 0.015      # Total specific humidity [kg/kg]
q_liq = 0.005      # Liquid specific humidity [kg/kg]
q_ice = 0.001      # Ice specific humidity [kg/kg]

# 3. Compute properties directly
T = TD.air_temperature(params, e_int, q_tot, q_liq, q_ice)
p = TD.air_pressure(params, T, ρ, q_tot, q_liq, q_ice)
(T=T, p=p)
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
  - Thermodynamics Parameters
  - Types
  - Thermodynamic Functions
  - Temperature Profiles
  - Data Collection

### 🛠️ **User Guides**

- **[How-To Guide](HowToGuide.md)** - Practical usage examples and patterns
  - Installation and setup
  - Common use cases and workflows

### 🔬 **Advanced Topics**

- **[Temperature Profiles](TemperatureProfiles.md)** - Pre-defined atmospheric profiles to be used as reference states in atmosphere models and for testing
- **[Tested Profiles](TestedProfiles.md)** - Thermodynamic profiles used for testing of the package (Internal Use)

### 👨‍💻 **Developer Resources**

- **[Saturation Adjustment Convergence](SaturationAdjustmentConvergence.md)** - Convergence testing for numerical methods

### 📚 **Published References and Background**

- **[References](References.md)** - Bibliography of theoretical foundations

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
- **Comprehensive testing** and validation suite.

## Core Design Principles

### **Functional & Stateless**

Functions in Thermodynamics.jl are pure and stateless. They take a `ThermodynamicsParameters` struct and the necessary thermodynamic variables (e.g., `T`, `ρ`, `q`...) as arguments. This design fits naturally into large-scale simulations (like `ClimaAtmos.jl`) where prognostic variables are managed externally.

### **Working Fluid**

The working fluid is **moist air** (dry air + water vapor + liquid water + ice, which may include precipitation). We treat it as a mixture of ideal gases and condensed phases, ensuring rigorous mass and energy conservation.

### **Consistent Formulation**

All quantities are derived from the **calorically perfect gas** assumption with constant specific heat capacities. This provides a consistent, closed consistency set of equations for saturation vapor pressures (the so-called Rankine-Kirchhoff approximation), latent heats, and other derived quantities.

## Usage Examples

### **Equilibrium Thermodynamics (Saturation Adjustment)**

To find the equilibrium temperature and phase partition from non-equilibrium variables (e.g., given `ρ`, `e_int`, `q_tot`), use `saturation_adjustment`:

```@example Index
# Define thermodynamic variables: internal energy, density, and total humidity
ρ = 1.0
e_int = -7.0e4
q_tot = 0.01

# Solve for equilibrium
# (using SecantMethod, but NewtonsMethod is also available for ρe)
sol = TD.saturation_adjustment(
    RS.SecantMethod, 
    params, 
    TD.ρe(), 
    ρ, e_int, q_tot, 
    10, 1e-4
)

println("T: ", sol.T)
println("q_liq: ", sol.q_liq)
println("q_ice: ", sol.q_ice)
```

### **Non-Equilibrium Thermodynamics**

Directly compute properties from thermodynamic variables without saturation adjustment:

```@example Index
# Define variables
q_tot = 0.01
q_liq = 0.005
q_ice = 0.0003
ρ = 1.0
e_int = -7.0e4

# Compute temperature directly
T = TD.air_temperature(params, e_int, q_tot, q_liq, q_ice)

# Compute pressure directly
p = TD.air_pressure(params, T, ρ, q_tot, q_liq, q_ice)
```

### **Saturation Calculations**

```@example Index
T = 280.0
ρ = 1.0
q_tot = 0.01

# Saturation vapor pressure
p_v_sat = TD.saturation_vapor_pressure(params, T, TD.Liquid())

# Saturation specific humidity
q_v_sat = TD.q_vap_saturation(params, T, ρ)

# Phase partition assuming equilibrium
(q_liq, q_ice) = TD.condensate_partition(params, T, ρ, q_tot)
```

## Integration with Climate Models

### **Dynamical Core Integration**

The package is designed for seamless integration with atmospheric dynamical cores, schematically as follows:

```@example Index
# Initialize
import Thermodynamics as TD
import RootSolvers as RS
using ClimaParams

FT = Float64
params = TD.Parameters.ThermodynamicsParameters(FT)

# Define physical constants and initial conditions
grav = TD.Parameters.grav(params)
cv_d = TD.Parameters.cv_d(params)
z = 1000.0   # m
geopotential = grav * z
q_tot = 0.01
ρ = 1.0
T = 300.0

# Initial velocity components
u, v, w = 10.0, 5.0, 1.0
e_kin = 0.5 * (u^2 + v^2 + w^2)

# Specific total energy e_tot = e_int + e_kin + geopotential
(q_liq, q_ice) = TD.condensate_partition(params, T, ρ, q_tot)
e_int = TD.internal_energy(params, T, q_tot, q_liq, q_ice)
e_tot = e_int + e_kin + geopotential

# Timestepping loop (simplified example)
for timestep in 1:10
    # Advance dynamical variables (simplified)
    # Note: In a real simulation, we would update the state vector.
    # Here we just re-use the variable names for demonstration.
    u_new = u + FT(0.1) 
    e_tot_new = e_tot + cv_d * FT(1)
    
    # Compute new internal energy
    e_kin_new = 0.5 * (u_new^2 + v^2 + w^2)
    e_int_new = e_tot_new - e_kin_new - geopotential
    
    # Saturation adjustment to get T, q_liq, q_ice
    local sol = TD.saturation_adjustment(
        RS.NewtonsMethod, 
        params, 
        TD.ρe(), 
        ρ, e_int_new, q_tot, 
        10, 1e-4
    )
    
    T_new = sol.T
    q_liq_new = sol.q_liq
    q_ice_new = sol.q_ice
    
    # Compute pressure for dynamics
    p_new = TD.air_pressure(params, T_new, ρ, q_tot, q_liq_new, q_ice_new)
    
    # Update state for next iteration (in real code)
    global u = u_new
    global e_tot = e_tot_new
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
