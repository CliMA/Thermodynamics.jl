# How-To Guide

This guide covers the essential aspects of using `Thermodynamics.jl`, specifically focusing on its **functional, stateless API**.

## Table of Contents

1. [Basic Setup](#basic-setup)
2. [Core Workflow](#core-workflow)
3. [Thermodynamic Calculations](#thermodynamic-calculations)
4. [Performance Considerations](#performance-considerations)
5. [Integration with Models](#integration-with-models)

## Basic Setup

### Installation and Import

```julia
using Pkg
Pkg.add("Thermodynamics")
Pkg.add("ClimaParams")
```

```@example HowToGuide
import Thermodynamics as TD
import RootSolvers as RS  # Needed for saturation adjustment
using ClimaParams
```

### Creating Parameters

```@example HowToGuide
# Create thermodynamic parameters for Float64 precision
params = TD.Parameters.ThermodynamicsParameters(Float64);

# For Float32 precision (useful for GPU computations)
params_f32 = TD.Parameters.ThermodynamicsParameters(Float32);
```

## Core Workflow

`Thermodynamics.jl` operates on **independent thermodynamic variables** directly, rather than wrapping them in a state configuration object. The general pattern is:

1. **Define your independent variables** (e.g., density `ρ`, internal energy `e_int`, specific humidities `q_...`).
2. **Pass them to a function** along with the parameter set `params`.

```@example HowToGuide
# Define variables
ρ = 1.0          # Density [kg/m³]
e_int = -7.0e4   # Internal energy [J/kg]
q_tot = 0.01     # Total specific humidity [kg/kg]
q_liq = 0.005    # Liquid specific humidity [kg/kg]
q_ice = 0.0      # Ice specific humidity [kg/kg]

# Compute a property
T = TD.air_temperature(params, e_int, q_tot, q_liq, q_ice)
p = TD.air_pressure(params, T, ρ, q_tot, q_liq, q_ice)

(T=T, p=p)
```

!!! tip "Temperature Profiles for Testing"
    For testing and validation, you can use pre-defined atmospheric temperature profiles from `TD.TemperatureProfiles`. See [Temperature Profiles](@ref) for available profiles and usage examples.

## Thermodynamic Calculations

### **1. Equilibrium Calculations (Saturation Adjustment)**

When you have conservative variables (e.g., `ρ`, `e_int`, `q_tot`) and need to find the temperature `T` and phase partition `(q_liq, q_ice)` that satisfy thermodynamic equilibrium (saturation), use [`saturation_adjustment`](@ref).

```@example HowToGuide
# Input variables
ρ = 1.0
e_int = -7.0e4
q_tot = 0.01

# Solve for equilibrium
# We select a root-finding method from RootSolvers.jl (e.g., SecantMethod, NewtonsMethod)
# We specify the formulation (TD.ρe() means inputs are Density & Internal Energy)
sol = TD.saturation_adjustment(
    RS.SecantMethod,        # Method
    params,                 # Parameter set
    TD.ρe(),                # Formulation
    ρ, e_int, q_tot,        # Inputs
    10,                     # Max iterations
    1e-3                    # Tolerance
)

T_equil = sol.T
q_liq_equil = sol.q_liq
q_ice_equil = sol.q_ice
converged = sol.converged

println("Equilibrium T: $T_equil")
```

**Supported formulations for saturation adjustment:**

- `TD.ρe()`: Inputs: `ρ`, `e_int`, `q_tot`
- `TD.pe()`: Inputs: `p`, `e_int`, `q_tot`
- `TD.ph()`: Inputs: `p`, `h`, `q_tot`
- `TD.pθ_li()`: Inputs: `p`, `θ_li`, `q_tot`

### **2. Non-Equilibrium Calculations (Explicit Phases)**

If you already know the phase partition (e.g., `q_liq` and `q_ice` are prognostic variables in your model), you can compute properties directly without iteration.

```@example HowToGuide
q_tot = 0.01
q_liq = 0.002
q_ice = 0.001

# Temperature from internal energy
e_int = -6.5e4
T = TD.air_temperature(params, e_int, q_tot, q_liq, q_ice)

# Temperature from enthalpy
h = -5.0e4
T_from_h = TD.air_temperature(params, TD.ph(), h, q_tot, q_liq, q_ice)

# Pressure
p = TD.air_pressure(params, T, ρ, q_tot, q_liq, q_ice)

# Sound speed
c_s = TD.soundspeed_air(params, T, q_tot, q_liq, q_ice)
```

### **3. Saturation Properties**

Compute saturation properties given thermodynamic variables.

```@example HowToGuide
T = 290.0
ρ = 1.1

# Saturation vapor pressure over liquid
p_v_sat_l = TD.saturation_vapor_pressure(params, T, TD.Liquid())

# Saturation vapor pressure over ice
p_v_sat_i = TD.saturation_vapor_pressure(params, T, TD.Ice())

# Saturation specific humidity
q_v_sat = TD.q_vap_saturation(params, T, ρ)

(p_v_sat_l=p_v_sat_l, p_v_sat_i=p_v_sat_i, q_v_sat=q_v_sat)
```

## Performance Considerations

### **Type Stability**

Ensure all inputs map to the same floating-point type (e.g., `Float64` or `Float32`). Mixed precision can cause allocations and slowdowns.

```@example HowToGuide
# Good: Consistent types
FT = Float64
params = TD.Parameters.ThermodynamicsParameters(FT)
val = TD.air_temperature(params, FT(-7.0e4), FT(0.01))

# Avoid: Mixed types (e.g. Float64 params with Float32 inputs)
```

### **Vectorized Operations**

`Thermodynamics.jl` functions broadcast efficiently over arrays.

```@example HowToGuide
# Arrays of thermodynamic variables
e_int_arr = [-7.0e4, -6.5e4, -6.0e4]
densities = [1.0, 1.1, 1.2]
q_tots    = [0.01, 0.012, 0.015]

# Broadcast the function call
T_arr = TD.air_temperature.(Ref(params), e_int_arr, q_tots)
```

### **GPU Compatibility**

The functional API is GPU-friendly. Usage is identical, provided arrays are on the GPU (e.g., `CuArray`) and parameters are initialized with the correct type (`Float32`).

```julia
# Use Float32 for GPU
FT = Float32
params_f32 = TD.Parameters.ThermodynamicsParameters(FT)

# On GPU (conceptual)
# T_gpu = TD.air_temperature.(Ref(params_f32), e_int_gpu, q_tot_gpu)
```

## Integration with Models

### **Prognostic Variable Management**

In a weather or climate model, you typically evolve a state vector. `Thermodynamics.jl` acts as a kernel to close the system of equations.

```@example HowToGuide
# Example: Computing pressure for the momentum equation
# ----------------------------------------------------
# Prognostic variables from the model state:
ρ = 1.0 # arbitrary value for example
e_tot = 20000.0 # arbitrary value for example
u, v, w = 10.0, 0.0, 0.0 # arbitrary values
q_tot, q_liq, q_ice = 0.01, 0.001, 0.0 # arbitrary values

# 1. Recover internal energy
e_kin = 0.5 * (u^2 + v^2 + w^2)
e_pot = 0.0 # arbitrary
e_int = e_tot - e_kin - e_pot

# 2. Recover temperature (non-equilibrium if q_liq/ice are prognostic)
T = TD.air_temperature(params, e_int, q_tot, q_liq, q_ice)

# 3. Compute pressure
p = TD.air_pressure(params, T, ρ, q_tot, q_liq, q_ice)
```

### **Handling Phase Changes**

If your model handles phase changes (microphysics), you might step `q_liq` and `q_ice` explicitly. If you assume instantaneous equilibrium (saturation adjustment), you use [`saturation_adjustment`](@ref) at the end or beginning of the step.

```@example HowToGuide
# Saturation Adjustment Step
# Re-define variables for example completeness
ρ = 1.0
e_int = -7.0e4
q_tot = 0.01
sol = TD.saturation_adjustment(RS.NewtonsMethod, params, TD.ρe(), ρ, e_int, q_tot, 15, 1e-4)

# Update thermodynamic variables
T_new = sol.T
q_liq_new = sol.q_liq
q_ice_new = sol.q_ice

(T_new, q_liq_new, q_ice_new)
```

## Common Pitfalls

### **Input Units**

- **SI Units**: All inputs must be in SI units (kg, m, s, K, Pa, J).
- **Specific Quantities**: Energies and humidities are *specific* (per kg of moist air).

### **Specific Humidity Definition**

- `q_tot`, `q_liq`, `q_ice` are mass fractions (kg/kg moist air).
- Ensure `q_liq + q_ice <= q_tot`. `q_vapor` is implicitly `q_tot - q_liq - q_ice`.
