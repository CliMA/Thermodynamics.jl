# How-To Guide

This guide covers the essential aspects of using `Thermodynamics.jl`, specifically focusing on its **functional, stateless API**.

## Table of Contents

1. [Basic Setup](@ref)
2. [Core Workflow](@ref)
3. [Thermodynamic Calculations](@ref)
4. [Performance Considerations](@ref)
5. [Integration with Models](@ref)
6. [Automatic Differentiation](@ref)
7. [Common Pitfalls](@ref)

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

### **1. Phase Equilibrium Calculations (Saturation Adjustment)**

When you have conservative variables (e.g., `ρ`, `e_int`, `q_tot`) and need to find the temperature `T` and phase partition `(q_liq, q_ice)` that satisfy phase equilibrium, use [`saturation_adjustment`](@ref).

```@example HowToGuide
# Input variables
ρ = 1.0
e_int = -7.0e4
q_tot = 0.01

# Solve for phase equilibrium
# We use the convenience method which handles defaults automatically
# For ρe(), this defaults to the optimized fixed-iteration solver
sol = TD.saturation_adjustment(
    params,                 # Parameter set
    TD.ρe(),                # Formulation
    ρ, e_int, q_tot         # Inputs
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

### **2. Phase Non-Equilibrium Calculations (Explicit Phases)**

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

#### **Optimized Saturation Adjustment**

For GPU performance, avoiding branch divergence is important. For the `ρe` formulation, use the fixed-iteration path by passing `forced_fixed_iters=true` as the last positional argument.

```julia
# GPU-optimized broadcasting using the convenience method
# This automatically uses the fast, branch-free fixed-iteration solver
sol = TD.saturation_adjustment.(
    Ref(params_f32),
    Ref(TD.ρe()),
    ρ_gpu, e_int_gpu, q_tot_gpu
)
```

For CPU single calls, you can omit the last two arguments to use the standard solver:

```julia
# CPU single-call (uses defaults)
sol = TD.saturation_adjustment(
    params_f32,
    TD.ρe(),
    ρ_val, e_int_val, q_tot_val
)
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

# 2. Recover temperature (if q_liq/ice are prognostic)
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
sol = TD.saturation_adjustment(params, TD.ρe(), ρ, e_int, q_tot)

# Update thermodynamic variables
T_new = sol.T
q_liq_new = sol.q_liq
q_ice_new = sol.q_ice

(T_new, q_liq_new, q_ice_new)
```

## Automatic Differentiation

`Thermodynamics.jl` is compatible with automatic differentiation (AD) tools such as [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl). This enables computing thermodynamic derivatives for sensitivity analysis, optimization, and adjoint-based methods.

### **Computing Derivatives Through Saturation Adjustment**

A common use case is computing how phase variables change with respect to input thermodynamic variables. For example, computing `∂q_liq/∂e_int` (the rate of change of liquid condensate with internal energy):

```@example HowToGuide
using ForwardDiff

# Setup: conditions at T = 290 K
T_sat = 290.0
ρ = 1.2

# Create supersaturated conditions (q_tot > q_vap_sat)
q_vap_sat = TD.q_vap_saturation(params, T_sat, ρ)
q_tot = q_vap_sat * 1.5  # 50% supersaturated

# Get equilibrium partition
(q_liq_eq, q_ice_eq) = TD.condensate_partition(params, T_sat, ρ, q_tot)

# Compute consistent internal energy
e_int = TD.internal_energy(params, T_sat, q_tot, q_liq_eq, q_ice_eq)

# Define a function from e_int → q_liq through saturation adjustment
function q_liq_from_e_int(e)
    sol = TD.saturation_adjustment(
        params,
        TD.ρe(),
        ρ, e, q_tot
    )
    return sol.q_liq
end

# Compute ∂q_liq/∂e_int using ForwardDiff
dq_liq_de_int = ForwardDiff.derivative(q_liq_from_e_int, e_int)

println("∂q_liq/∂e_int = ", dq_liq_de_int, " [kg/kg per J/kg]")
```

### **AD-Compatible Derivatives**

All core thermodynamic functions and saturation adjustment routines propagate dual numbers correctly:

```@example HowToGuide
# ∂T/∂e_int at fixed composition
function T_from_e_int(e)
    q_tot, q_liq, q_ice = 0.01, 0.002, 0.0
    TD.air_temperature(params, e, q_tot, q_liq, q_ice)
end

e_int_ref = -7.0e4
dT_de_int = ForwardDiff.derivative(T_from_e_int, e_int_ref)

println("∂T/∂e_int = ", dT_de_int, " [K per J/kg]")
```

!!! note "Root Solver Compatibility"
    When using AD through [`saturation_adjustment`](@ref), prefer `RS.SecantMethod` or `RS.NewtonsMethod` for smooth derivative propagation.

## Common Pitfalls

### **Input Units**

- **SI Units**: All inputs must be in SI units (kg, m, s, K, Pa, J).
- **Specific Quantities**: Energies and humidities are *specific* (per kg of moist air).

### **Specific Humidity Definition**

- `q_tot`, `q_liq`, `q_ice` are mass fractions (kg/kg moist air).
- Ensure `q_liq + q_ice <= q_tot`. `q_vapor` is implicitly `q_tot - q_liq - q_ice`.
