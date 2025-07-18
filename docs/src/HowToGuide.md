# How-To Guide

This guide covers the essential aspects of using Thermodynamics.jl, from basic usage to advanced patterns.

## Table of Contents

1. [Basic Setup](#basic-setup)
2. [Thermodynamic States](#thermodynamic-states)
3. [Common Usage Patterns](#common-usage-patterns)
4. [Performance Considerations](#performance-considerations)
5. [Integration with Models](#integration-with-models)
6. [Extending the Package](#extending-the-package)

## Basic Setup

### Installation and Import

```julia
using Pkg
Pkg.add("Thermodynamics")

import Thermodynamics as TD
```

### Creating Parameters

```julia
# Create thermodynamic parameters for Float64 precision
params = TD.Parameters.ThermodynamicsParameters(Float64)

# For Float32 precision (useful for GPU computations)
params_f32 = TD.Parameters.ThermodynamicsParameters(Float32)
```

## Thermodynamic States

Thermodynamics.jl uses a state-based approach where you create a thermodynamic state from independent variables, then compute any other thermodynamic property from that state.

!!! tip "Temperature Profiles"
    For testing and validation, you can use pre-defined atmospheric temperature profiles.
    See [Temperature Profiles](TemperatureProfiles.md) for available profiles and usage examples.

### Available State Constructors

#### **Equilibrium States** (Saturation Adjustment)
```julia
# From density, internal energy, and total humidity
ts = TD.PhaseEquil_ρeq(params, ρ, e_int, q_tot)

# From density, potential temperature, and total humidity  
ts = TD.PhaseEquil_ρθq(params, ρ, θ_liq_ice, q_tot)

# From pressure, internal energy, and total humidity
ts = TD.PhaseEquil_peq(params, p, e_int, q_tot)
```

#### **Non-Equilibrium States** (Explicit Phase Partitioning)
```julia
# With explicit liquid and ice partitioning
ts = TD.PhaseNonEquil(params, e_int, ρ, TD.PhasePartition(q_tot, q_liq, q_ice))

# Or create phase partition separately
q = TD.PhasePartition(q_tot, q_liq, q_ice)
ts = TD.PhaseNonEquil(params, e_int, ρ, q)
```

#### **Dry Air States**
```julia
# Dry air from pressure and temperature
ts = TD.PhaseDry_pT(params, p, T)
```

### Extracting Properties from States

```julia
# Basic thermodynamic properties
T = TD.air_temperature(params, ts)      # Temperature
p = TD.air_pressure(params, ts)         # Pressure
ρ = TD.air_density(params, ts)          # Density

# Phase partitioning
q = TD.PhasePartition(params, ts)       # Complete phase partition
q_tot = TD.total_specific_humidity(params, ts)  # Total humidity
q_liq = TD.liquid_specific_humidity(params, ts) # Liquid humidity
q_ice = TD.ice_specific_humidity(params, ts)    # Ice humidity
q_vap = TD.vapor_specific_humidity(params, ts)  # Vapor humidity

# Alternative: Extract directly from PhasePartition object
q = TD.PhasePartition(params, ts)
q_tot = q.tot  # Total specific humidity
q_liq = q.liq  # Liquid specific humidity  
q_ice = q.ice  # Ice specific humidity
q_vap = q.tot - q.liq - q.ice  # Vapor specific humidity (computed)

# Energy quantities
e_int = TD.internal_energy(params, ts)  # Internal energy
```

### Direct Function Usage

Thermodynamic functions can also be used directly without creating state objects, which can be more efficient because it avoids storing the state in memory:

```julia
# Direct phase partitioning (for non-equilibrium)
q = TD.PhasePartition(q_tot, q_liq, q_ice)         # Create phase partition directly

# Direct temperature calculations
T = TD.air_temperature(params, e_int, q)           # From internal energy and humidity
T = TD.air_temperature_from_enthalpy(params, h, q) # From enthalpy and humidity

# Direct pressure calculations  
p = TD.air_pressure(params, ρ, T, q)               # From density, temperature, humidity

# Direct humidity calculations
q_vap_sat = TD.q_vap_saturation(params, T, ρ, PhaseEquil{Float64})  # Saturation vapor humidity
p_v_sat = TD.saturation_vapor_pressure(params, T, Liquid())         # Saturation vapor pressure over liquid

# Direct energy calculations
e_int = TD.internal_energy(params, T, q)          # From temperature and humidity
h = TD.specific_enthalpy(params, T, q)            # From temperature and humidity

```

!!! tip "When to Use Each Approach"
    - **State-based approach**: Use when you need multiple properties from the same thermodynamic state and storing them is more efficient than re-computing them each time (e.g., because of iterative saturation adjustment)
    - **Direct functions**: Use for calculations when you already have the required variables or re-computing them is more efficient than storing them
    - **Saturation adjustment**: State constructors automatically handle saturation adjustment

## Performance Considerations

### **Type Stability**
```julia
# Good: Type-stable operations
params = TD.Parameters.ThermodynamicsParameters(Float64)
ts = TD.PhaseEquil_ρeq(params, ρ, e_int, q_tot)
T = TD.air_temperature(params, ts)

# Avoid: Mixed precision (can cause type instability)
params_f64 = TD.Parameters.ThermodynamicsParameters(Float64)
ts = TD.PhaseEquil_ρeq(params_f64, Float32(ρ), Float32(e_int), Float32(q_tot))
```

### **Vectorized Operations**
```julia
# For arrays of thermodynamic states
ρ_array = [1.0, 1.1, 1.2]
e_int_array = [2.0e5, 2.1e5, 2.2e5]
q_tot_array = [0.01, 0.012, 0.008]

# State-based approach 
ts_array = [TD.PhaseEquil_ρeq(params, ρ, e_int, q_tot) 
            for (ρ, e_int, q_tot) in zip(ρ_array, e_int_array, q_tot_array)]
T_array = [TD.air_temperature(params, ts) for ts in ts_array]

# Direct function approach 
T_array = [TD.air_temperature(params, e_int, TD.PhasePartition(q_tot)) 
           for (e_int, q_tot) in zip(e_int_array, q_tot_array)]
```

### **GPU Compatibility**
```julia
# Use Float32 for GPU computations
params_gpu = TD.Parameters.ThermodynamicsParameters(Float32)

# Ensure all inputs are Float32
ts_gpu = TD.PhaseEquil_ρeq(params_gpu, Float32(ρ), Float32(e_int), Float32(q_tot))
```

## Integration with Models

### **Dynamical Core Integration**

```julia
# Initialize prognostic variables
ρ = initial_density
e_tot = initial_total_energy
q_tot = initial_total_humidity

# Timestepping loop
for timestep in timesteps
    # Advance prognostic variables
    advance_dynamics!(ρ, e_tot, q_tot)
    
    # Extract internal energy
    e_kin = compute_kinetic_energy(u, v, w)
    e_pot = compute_potential_energy(z)
    e_int = e_tot - e_kin - e_pot
    
    # Saturation adjustment (state-based approach)
    ts = TD.PhaseEquil_ρeq(params, ρ, e_int, q_tot)
    T = TD.air_temperature(params, ts)
    q = TD.PhasePartition(params, ts)
    
    # Or direct function approach
    T = TD.air_temperature(params, e_int, TD.PhasePartition(q))  # No saturation adjustment
    ts = TD.PhaseEquil_ρeq(params, ρ, e_int, q_tot)
    T = TD.air_temperature(params, ts)
    q = TD.PhasePartition(params, ts)
    
    # Use temperature in physics (radiation, etc.)
    compute_physics!(T, q)
end
```

## Extending the Package

### **Adding New Thermodynamic State Constructors**

If Thermodynamics.jl doesn't have a constructor for your specific use case, you can implement one in `src/states.jl`. The constructor must translate your inputs into one of the fundamental state types:

```julia
# Example: Constructor from pressure, temperature, and humidity
function PhaseEquil_pTq(param_set::APS, p::FT, T::FT, q_tot::FT) where {FT}
    # Compute density from equation of state
    ρ = TD.air_density(param_set, T, p, TD.PhasePartition(q_tot))
    
    # Compute internal energy
    e_int = TD.internal_energy(param_set, T, TD.PhasePartition(q_tot))
    
    # Return equilibrium state
    return TD.PhaseEquil_ρeq(param_set, ρ, e_int, q_tot)
end
```

### **Available Base State Types**

- **`PhaseDry`**: Dry air state (2 independent variables)
- **`PhaseEquil`**: Moist air in thermodynamic equilibrium (3 independent variables)
- **`PhaseNonEquil`**: Moist air in non-equilibrium (3+ independent variables)

### **Best Practices for Extensions**

1. **Maintain type stability**: Use consistent floating-point types
2. **Follow naming conventions**: Use descriptive names with parameter types
3. **Add tests**: Include unit tests for new constructors using [Tested Profiles](TestedProfiles.md)
4. **Document thoroughly**: Add docstrings explaining the constructor's purpose

## Common Pitfalls and Solutions

### **Pitfall 1: Incorrect State Type**
```julia
# Wrong: Using equilibrium constructor for non-equilibrium conditions
ts = TD.PhaseEquil_ρeq(params, ρ, e_int, q_tot)  # Assumes saturation adjustment

# Right: Use non-equilibrium constructor when you have explicit partitioning
ts = TD.PhaseNonEquil(params, e_int, ρ, q)
```

### **Pitfall 2: Mixed Precision**
```julia
# Wrong: Mixing Float32 and Float64
params = TD.Parameters.ThermodynamicsParameters(Float64)
ts = TD.PhaseEquil_ρeq(params, Float32(ρ), Float32(e_int), Float32(q_tot))

# Right: Consistent precision
params = TD.Parameters.ThermodynamicsParameters(Float32)
ts = TD.PhaseEquil_ρeq(params, Float32(ρ), Float32(e_int), Float32(q_tot))
```

## Next Steps

1. **Explore the [API Reference](API.md)** for complete function documentation
2. **Read the [Mathematical Formulation](Formulation.md)** for theoretical background
3. **Check [Saturation Adjustment Convergence](SaturationAdjustmentConvergence.md)** for numerical method testing
4. **Explore the API Reference** for complete function documentation

---

!!! tip "Getting Help"
    For specific questions or issues, check the documentation or open an issue on the [GitHub repository](https://github.com/CliMA/Thermodynamics.jl).
