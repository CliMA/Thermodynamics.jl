# Mathematical Formulation

## Table of Contents

1. [Introduction](#1-introduction)
2. [Fundamental Assumptions](#2-fundamental-assumptions)
3. [Working Fluid and Equation of State](#3-working-fluid-and-equation-of-state)
   - [3.1 Mass Fractions and Notation](#31-mass-fractions-and-notation)
   - [3.2 Equation of State](#32-equation-of-state)
4. [Heat Capacities](#4-heat-capacities)
5. [Latent Heats](#5-latent-heats)
6. [Internal Energies](#6-internal-energies)
7. [Enthalpies](#7-enthalpies)
8. [Moist Static Energy](#8-moist-static-energy)
9. [Saturation Vapor Pressure](#9-saturation-vapor-pressure)
10. [Saturation Specific Humidity](#10-saturation-specific-humidity)
11. [Saturation Adjustment](#11-saturation-adjustment)
12. [Auxiliary Thermodynamic Functions](#12-auxiliary-thermodynamic-functions)
    - [12.1 Relative Humidity](#121-relative-humidity)
    - [12.2 Potential Temperature](#122-potential-temperature)
    - [12.3 Virtual Temperature and Virtual Potential Temperature](#123-virtual-temperature-and-virtual-potential-temperature)
    - [12.4 Liquid-Ice Potential Temperature](#124-liquid-ice-potential-temperature)
    - [12.5 Speed of Sound](#125-speed-of-sound)
13. [Summary and Implementation Guidelines](#13-summary-and-implementation-guidelines)
    - [13.1 Key Theoretical Framework](#131-key-theoretical-framework)
    - [13.2 Core Equations](#132-core-equations)
    - [13.3 Implementation Strategy](#133-implementation-strategy)
    - [13.4 Validation and Testing](#134-validation-and-testing)
    - [13.5 Extensions and Limitations](#135-extensions-and-limitations)

!!! note "Cross-References"
    This documentation is closely related to the [API documentation](API.md) which provides detailed function signatures and usage examples. The theoretical framework described here is implemented in the [`Thermodynamics.jl`](https://github.com/CliMA/Thermodynamics.jl) package.

## 1. Introduction

Here we introduce one consistent set of thermodynamic approximations for all model components. The key to thermodynamic consistency at reasonable accuracy is to take the specific heat capacities of the constituents of moist air (dry air, water vapor, liquid water, and ice) to be constant, i.e., to assume the gases to be **calorically perfect**. We discuss how to derive all other thermodynamic quantities that are needed on the basis of this one approximation ([Yatunin2026](@cite)[Romps2008](@cite),[Ambaum2020](@cite), [Romps2021](@cite)). This includes:

- Giving accurate and easily adaptable closed-form expressions for internal energies, enthalpies, specific latent heats, and saturation vapor pressures
- Showing how to construct consistent sets of thermodynamic equations that either (i) assume phase equilibrium and require only one prognostic water variable, or (ii) do not assume phase equilibrium (but do assume thermal equilibrium) and require prognostic variables for all water phases
- Showing how to obtain temperatures from energy variables under either phase equilibrium assumptions (by `saturation adjustment`) or phase non-equilibrium assumptions (by a closed-form expression for temperature).

The resulting thermodynamic functions are implemented in [`Thermodynamics.jl`](https://github.com/CliMA/Thermodynamics.jl).

Specific thermodynamic formulations often vary in how they approximate the relevant material properties. The formulation used in `Thermodynamics.jl` balances three criteria:

1. **Accuracy**: The formulation is accurate enough for atmospheric modeling (e.g., errors in saturation vapor pressure are within a few percent).
2. **Consistency**: The formulation is thermodynamically consistent (e.g., it conserves energy and satisfies the Clausius-Clapeyron relation).
3. **Simplicity and Efficiency**: The formulation leads to closed-form expressions that can be evaluated efficiently.

The implementation follows the thermodynamic formulation of the CliMA Earth System Model ([Yatunin2026](@cite)). It relies on the **Rankine-Kirchhoff approximations**, which provide a consistent framework for moist thermodynamics ([Romps2021](@cite)).

Specific choices of thermodynamic constants have been made in [ClimaParams.jl](https://github.com/CliMA/ClimaParams.jl) to maximize accuracy given the Rankine-Kirchhoff approximations ([Ambaum2020](@cite), [Yatunin2026](@cite)).

!!! note "Physical Motivation"
    The assumption of calorically perfect gases (constant specific heat capacities) is justified for atmospheric conditions because the error of approximating them as constant is less than 1% for dry air, the main constituent of moist air, and at most a few percent for the water phases. This approximation enables closed-form expressions for all thermodynamic quantities while maintaining sufficient accuracy for atmospheric modeling.

## 2. Fundamental Assumptions

Our thermodynamic framework is based on the following fundamental assumptions:

1. **Ideal Gas Law**: Dry air and water vapor behave as ideal gases
2. **Calorically Perfect Gases**: Specific heat capacities are constant (temperature-independent)
3. **Negligible Condensate Volume**: The specific volume of liquid water and ice is neglected relative to gas phases
4. **Thermal Equilibrium**: All phases have the same temperature
5. **Reference State**: All thermodynamic quantities are defined relative to a reference temperature `T₀`

!!! tip "Implementation Note"
    These assumptions enable closed-form expressions for all thermodynamic quantities, making the implementation computationally efficient while maintaining sufficient accuracy for atmospheric modeling.

## 3. Working Fluid and Equation of State

The working fluid is **moist air**. That is, it is an ideal mixture of dry air, water vapor, and condensed water (liquid and ice). Atmospheric models may choose to include precipitation (e.g., rain, snow, graupel) within the definitions of the condensed water specific humidities ($q_{liq}$, $q_{ice}$), or treat them as separate species. The thermodynamic formulation remains valid in either case, provided that the phases included in the working fluid are in thermal equilibrium.

!!! note "Modeling Choice: Precipitation Included in Working Fluid"
    In the CliMA Earth System Model, precipitation is included as part of the working fluid. This choice means:
    - Mass and energy conservation are maintained across all phases (including precipitation)
    - Thermodynamic consistency is preserved throughout the system
    - All condensed water phases (cloud and precipitation) are assumed to be in thermal equilibrium with the surrounding air (which is an approximation, see [Yatunin2026](@cite) for discussion)

Dry air and water vapor are taken to be ideal gases. The specific volume of all condensed phases (cloud liquid, cloud ice, and potentially precipitation) is neglected relative to that of the gas phases (it is a factor $10^3$ less than that of the gas phases). All phases are assumed to have the same temperature. However, the condensates do not need to be in phase equilibrium with the other fluid constituents; out-of-equilibrium phases such as supercooled liquid can exist.

!!! tip "Implementation Note"
    When precipitation is included in the working fluid, its mass contributes to the total water content $q_t$ and affects the density and specific heat capacities. This is implemented consistently throughout the `Thermodynamics.jl` package.

---

!!! info "Sidebar: Pros and Cons of Including Precipitation in the Working Fluid"
    **Pros:**
    - Ensures mass and energy conservation across all phases
    - Simplifies the thermodynamic framework and implementation
    - Enables consistent inclusion of microphysics schemes that treat hydrometeors as occupying a continuous spectrum (without artificial spectral gaps between cloud condensate and precipitation)

    **Cons:**
    - Assumes precipitation is in thermal equilibrium with air, which may not always be true (e.g., large raindrops or hail)
    
    See [Yatunin2026](@cite) for further discussion and justification of this modeling choice.

---

### 3.1 Mass Fractions and Notation

The density of the moist air is denoted by $\rho$. We use the following notation for the mass fractions of the moist air mixture (mass of a constituent divided by the total mass of the working fluid):

| Symbol                | Definition                                                                 |
|-----------------------|----------------------------------------------------------------------------|
| $q_d$                 | dry air mass fraction                                                      |
| $q_v$                 | water vapor specific humidity                                              |
| $q_l$                 | liquid water specific humidity (**includes both cloud liquid and precipitating liquid**) |
| $q_i$                 | ice specific humidity (**includes both cloud ice and precipitating ice/snow**) |
| $q_c = q_l + q_i$     | condensate specific humidity (all condensed water, cloud + precipitation)  |
| $q_t = q_v + q_c$     | total specific humidity (all water phases)                                 |

!!! note "Precipitation Inclusion"
    The liquid and ice specific humidities $q_l$ and $q_i$ include both cloud condensate and precipitation. This unified treatment ensures that all condensed water phases are thermodynamically consistent and can exchange mass and energy with the gas phases.

Because this enumerates all constituents of the working fluid, we have $q_t + q_d = 1$. In Earth's atmosphere, the water vapor specific humidity $q_v$ generally dominates the total specific humidity $q_t$ and is usually $O(10^{-2})$ or smaller; the condensate specific humidity is typically $O(10^{-4})$. Hence, water is a trace constituent of the atmosphere, and only a small fraction of atmospheric water is in condensed phases.

### 3.2 Equation of State

The pressure $p$ of the working fluid is the sum of the partial pressures of dry air and water vapor, both taken to be ideal gases. Neglecting the volume of the condensed phases (but not their masses), this gives $p = \rho (R_d q_d + R_v q_v) T$, where $R_d$ is the specific gas constant of dry air, and $R_v$ is the specific gas constant of water vapor.

!!! note "Precipitation and the Equation of State"
    Although precipitation is included in the working fluid, it does **not** contribute to the pressure because its specific volume is neglected. However, precipitation mass affects the density andspecific heat capacity.

This can also be written as

```math
\begin{equation}
    p = \rho R_m T,
\label{e:eos}
\end{equation}
```

where

```math
\begin{equation}
\begin{aligned}
    R_m(q) & = R_d (1 - q_t) + R_v q_v \\
           & = R_d \left[ 1 + (\varepsilon_{dv}-1)q_t - \varepsilon_{dv} q_c\right]
\end{aligned}
\label{e:Rm}
\end{equation}
```

is the **specific gas constant of moist air** (which varies with composition); here, we have used $q_d = 1-q_t$ and $q_v = q_t - q_c$ and introduced the ratio of the molar masses of dry air and water vapor $\varepsilon_{dv} = R_v/R_d$ ($\approx 1.61$).

!!! note "Mathematical Note"
    Although called a "constant," $R_m$ actually varies with the composition of the moist air. This terminology reflects its role as the effective gas constant in the equation of state, similar to how the molecular weight of a gas mixture varies with composition.

Equations \eqref{e:eos} and \eqref{e:Rm} constitute the equation of state of the working fluid. We use the notation $q=(q_t, q_l, q_i)$ for the tuple of specific humidities that determine the composition of moist air.

!!! example "Typical Values"
    For Earth's atmosphere at sea level:

    | Quantity | Value |
    |----------|-------|
    | $R_d$ | $287.0$ J/(kg·K) |
    | $R_v$ | $461.5$ J/(kg·K) |
    | $\varepsilon_{dv}$ | $1.61$ |
    | $R_m$ (typical moist air, $q_t = 0.01$) | $288.7$ J/(kg·K) |
    | $R_m$ (air with precipitation, $q_t = 0.015$) | $289.0$ J/(kg·K) |

## 4. Heat Capacities

The isochoric specific heat capacities of the constituents of moist air are:

| Symbol      | Definition                                   |
|-------------|----------------------------------------------|
| $c_{vd}$    | Isochoric specific heat capacity of dry air   |
| $c_{vv}$    | Isochoric specific heat capacity of water vapor |
| $c_{vl}$    | Isochoric specific heat capacity of liquid water |
| $c_{vi}$    | Isochoric specific heat capacity of ice       |

Our key thermodynamic approximation is to take these isochoric specific heat capacities to be constants, i.e., we take the gases to be **calorically perfect**. This is an approximation because they depend weakly on temperature. But for atmospheric conditions, the error of approximating them as constant is less than 1% for dry air, the main constituent of moist air, and at most a few percent for the water phases.

The difference between the isochoric and isobaric specific heat capacities is proportional to the specific volume. Consistent with taking the specific volume of liquid water and ice to be zero, we take the isochoric and isobaric specific heat capacities of the condensed phases to be equal. The isobaric specific heat capacities of the constituents then are:

| Symbol                      | Definition                                      |
|-----------------------------|-------------------------------------------------|
| $c_{pd} = c_{vd} + R_d$     | Isobaric specific heat capacity of dry air       |
| $c_{pv} = c_{vv} + R_v$     | Isobaric specific heat capacity of water vapor   |
| $c_{pl} = c_{vl}$           | Isobaric specific heat capacity of liquid water  |
| $c_{pi} = c_{vi}$           | Isobaric specific heat capacity of ice           |

!!! tip "Implementation Note"
    The relationship $c_p = c_v + R$ for ideal gases follows from the definition of enthalpy $h = u + pv$ and the ideal gas law $pv = RT$. For condensed phases, we neglect the specific volume term ($v$), so $c_p ≈ c_v$.

The corresponding specific heat capacities of moist air are the weighted sum of those of the constituents:

```math
\begin{equation}
\begin{aligned}
    c_{\cdot m}(q) & = (1-q_t) c_{\cdot d} + q_v c_{\cdot v} + q_l c_{\cdot l} + q_i c_{\cdot i} \\
    & = c_{\cdot d} + (c_{\cdot v} - c_{\cdot d})q_t + (c_{\cdot l} - c_{\cdot v})q_l + (c_{\cdot i} - c_{\cdot v})q_i
\end{aligned}
\label{e:SpecificHeat}
\end{equation}
```

where $\cdot$ stands for $v$ or $p$ and we have used $q_v = q_t - q_l - q_i$.

!!! note "Mathematical Note"
    The second form of the equation is computationally more efficient as it avoids computing $q_v$ explicitly. This rearrangement is used in the implementation to improve performance.

Straightforward substitution shows that the above relation between the specific heat capacities of the constituents also holds for moist air as a whole:

```math
\begin{equation}\label{e:SpecificHeatRelation}
    c_{pm}(q) = c_{vm}(q) + R_m(q).
\end{equation}
```

!!! example "Typical Values"
    For Earth's atmosphere at standard conditions:

    | Quantity | Value |
    |----------|-------|
    | $c_{vd}$ | $717.6$ J/(kg·K) |
    | $c_{vv}$ | $1410.0$ J/(kg·K) |
    | $c_{vl}$ | $4219.0$ J/(kg·K) |
    | $c_{vi}$ | $2106.0$ J/(kg·K) |
    | $c_{vm}$ (typical moist air, $q_t = 0.01$) | $720.0$ J/(kg·K) |
    | $c_{pm}$ (typical moist air, $q_t = 0.01$) | $1008.0$ J/(kg·K) |

!!! tip "Implementation Note"
    The specific heat capacities are implemented as weighted sums in the [`cp_m`](@ref) and [`cv_m`](@ref) functions. The implementation uses the rearranged form of equation \eqref{e:SpecificHeat} for computational efficiency.

## 5. Latent Heats

Kirchhoff's relation states that the specific latent enthalpy (heat) $L$ of a phase change depends on temperature $T$ through

```math
\begin{equation}
    \frac{dL}{dT} = \Delta c_p,
\end{equation}
```

where $\Delta c_p$ is the difference in isobaric specific heat capacities between the phase with the higher and lower specific volume. For the constant isobaric specific heat capacities that we assume, this can be integrated to give

```math
\begin{equation}
    L(T) = L_0 + \Delta c_p (T-T_0),
    \label{e:LHTemperature}
\end{equation}
```

where $T_0$ is a reference temperature and $L_0$ is the specific latent heat at $T_0$.

!!! note "Physical Interpretation"
    Kirchhoff's relation follows from the fact that the enthalpy difference between phases changes with temperature due to the different heat capacities of the phases. The latent heat represents the energy required to transform a unit mass from one phase to another at constant pressure.

!!! tip "Implementation Note"
    The linear temperature dependence of latent heats enables closed-form expressions for saturation vapor pressure and other thermodynamic quantities. This approximation is accurate to within a few percent for atmospheric temperature ranges.

For the phase transitions of water, this implies specifically:

| Formula | Description |
|---------|-------------|
| $L_v(T) = L_{v,0} + (c_{pv} - c_{pl}) (T - T_0)$ | Latent heat of vaporization |
| $L_f(T) = L_{f,0} + (c_{pl} - c_{pi}) (T - T_0)$ | Latent heat of fusion |
| $L_s(T) = L_{s,0} + (c_{pv} - c_{pi}) (T - T_0)$ | Latent heat of sublimation |

With $L_{s,0} = L_{v,0} + L_{f,0}$, this gives $L_s(T) = L_v(T) + L_f(T)$, as it should.

!!! example "Typical Values"
    At the reference temperature $T_0 = 273.15$ K:

    | Quantity | Value |
    |----------|-------|
    | $L_{v,0}$ | $2.501 \times 10^6$ J/kg (latent heat of vaporization) |
    | $L_{f,0}$ | $0.334 \times 10^6$ J/kg (latent heat of fusion) |
    | $L_{s,0}$ | $2.835 \times 10^6$ J/kg (latent heat of sublimation) |

    At $T = 300$ K:

    | Quantity | Value |
    |----------|-------|
    | $L_v$ | $2.430 \times 10^6$ J/kg |
    | $L_f$ | $0.334 \times 10^6$ J/kg |
    | $L_s$ | $2.764 \times 10^6$ J/kg |

!!! tip "Implementation Note"
    The latent heats are implemented in the [`latent_heat_vapor`](@ref), [`latent_heat_fusion`](@ref), and [`latent_heat_sublim`](@ref) functions. The weighted latent heat for mixed-phase conditions is computed in [`latent_heat_mixed`](@ref).

## 6. Internal Energies

The specific internal energies of the constituents of moist air can be written as

```math
\begin{equation}
\begin{aligned}
I_d(T) & = c_{vd} (T - T_0) - R_d T_0,  \\
I_v(T) & = c_{vv} (T - T_0) + I_{v,0},\\
I_l(T) & = c_{vl} (T - T_0), \\
I_i(T) & = c_{vi} (T - T_0) - I_{i,0}.
\end{aligned}
\label{e:InternalEnergies}
\end{equation}
```

Here, the reference specific internal energy $I_{v,0}$ is the difference in specific internal energy between vapor and liquid at the reference temperature $T_0$, and $I_{i,0}$ is the difference in specific internal energy between ice and liquid at $T_0$. We have included an arbitrary constant offset $- R_d T_0$ in the definition of the dry specific internal energy as that simplifies the corresponding specific enthalpies \eqref{e:Enthalpies}. The formulation is **reference-temperature invariant**, meaning that the physics is independent of the choice of the arbitrary reference temperature $T_0$ used to define energies, enthalpies, and entropies, provided that boundary conditions (as implemented in [SurfaceFluxes.jl](https://github.com/CliMA/SurfaceFluxes.jl)) also respect this invariance. Measurable thermodynamic variables such as temperature, pressure, etc. do not depend on a shift in the reference temperature $T_0$.

!!! note "Physical Interpretation"
    The internal energy represents the total energy of a substance excluding kinetic and potential energy. The reference energies $I_{v,0}$ and $I_{i,0}$ represent the energy differences between phases at the reference temperature, accounting for the fact that vapor has higher internal energy than liquid, and ice has lower internal energy than liquid.

The internal energy of moist air is the weighted sum of that of the constituents:

```math
\begin{equation}
\begin{aligned}
     I(T, q) & = (1-q_t) I_d(T) + q_v I_v(T) + q_l I_l(T) + q_i I_i(T)\\
          & = c_{vm}(q) (T - T_0)  + q_v I_{v,0} - q_i I_{i,0} - (1 - q_t) R_d T_0.
\end{aligned}
\label{e:TotalInternalEnergy}
\end{equation}
```

!!! tip "Implementation Note"
    The internal energy is implemented as a weighted sum of constituent energies in the [`internal_energy`](@ref) function. The constituent energies are computed separately in [`internal_energy_dry`](@ref), [`internal_energy_vapor`](@ref), [`internal_energy_liquid`](@ref), and [`internal_energy_ice`](@ref).

The internal energy can be inverted to obtain the temperature given $I$ and the specific humidities:

```math
\begin{equation}
    T = T_0 + \frac{I - (q_t - q_l - q_i) I_{v,0} + q_i I_{i,0} + (1 - q_t) R_d T_0}{c_{vm}(q)},
    \label{e:temperature}
\end{equation}
```

where we have used $q_v = q_t - q_l - q_i$. This allows one to recover temperature given internal energy and specific humidities as state variables.

!!! note "Mathematical Note"
    The temperature recovery equation \eqref{e:temperature} is crucial to compute temperature from internal energy and composition when internal (or total) energy is used as a prognostic variable. This inversion is possible because internal energy is a monotonic function of temperature for our assumptions.

The reference specific internal energies $I_{v,0}$ and $I_{i,0}$ are related to the reference specific latent heats $L_{v,0}$ and $L_{f,0}$, which indicate the enthalpy differences between the phases at $T_0$. The reference specific internal energies are obtained from the reference specific latent heats by subtracting the "$pv$" term, which is $p_k/\rho_k$ for the relevant partial pressure $p_k$ and specific volume $1/\rho_k$ of the phase $k$ (and hence is zero for the condensed phases, whose specific volume we neglect). This gives

```math
\begin{equation}
\begin{aligned}
     I_{v,0} &= L_{v, 0} - R_v T_0,\\
     I_{i,0} &= L_{f, 0}.
\end{aligned}
\label{e:RefInternalEnergies}
\end{equation}
```

!!! example "Typical Values"
    At the reference temperature $T_0 = 273.15$ K:

    | Quantity | Value |
    |----------|-------|
    | $I_{v,0}$ | $2.501 \times 10^6$ J/kg (vapor reference energy) |
    | $I_{i,0}$ | $0.334 \times 10^6$ J/kg (ice reference energy) |

    For typical moist air at $T = 300$ K with $q_t = 0.01$:

    | Quantity | Value |
    |----------|-------|
    | $I$ | $21.5 \times 10^3$ J/kg (total internal energy) |

## 7. Enthalpies

The specific enthalpies of the constituents of moist air are obtained by adding $p_k/\rho_k$ for phase $k$ to the corresponding specific internal energy \eqref{e:InternalEnergies}. Again neglecting the specific volumes of the condensed phases and using the relations \eqref{e:RefInternalEnergies} between reference specific energies and latent heats, this gives:

```math
\begin{equation}
\label{e:Enthalpies}
\begin{aligned}
    h_d(T) = I_d(T) + R_d T &= c_{pd}(T-T_0), \\
    h_v(T) = I_v(T) + R_v T &= c_{pv}(T-T_0) + L_{v,0}, \\
    h_l(T) = I_l(T) &= c_{pl}(T-T_0), \\
    h_i(T) = I_i(T) &= c_{pi}(T-T_0) - L_{f,0}.
\end{aligned}
\end{equation}
```

The formulation is **reference-temperature invariant**, meaning that the physics is independent of the choice of the arbitrary reference temperature $T_0$ used to define enthalpies and entropies, provided that boundary conditions (as implemented in [SurfaceFluxes.jl](https://github.com/CliMA/SurfaceFluxes.jl)) also respect this invariance.

!!! note "Physical Interpretation"
    Enthalpy represents the total energy of a substance including the work done against pressure. For ideal gases, enthalpy includes the $pv$ term, while for condensed phases this term is neglected due to their small specific volume $v$. Enthalpy is the relevant energy quantity for fluid transport, including at boundaries.

The enthalpy of moist air is the weighted sum of the constituent enthalpies:

```math
\begin{equation}
\begin{split}\label{e:enthalpy_definition}
    h(T, q)  &= (1-q_t) h_d + q_v h_v + q_l h_l + q_i h_i \\
        &= c_{pm}(q) (T-T_0) + q_v L_{v,0} - q_i L_{f,0}\\
        &= I(q, T) + R_m T,
\end{split}
\end{equation}
```

where the last equality used $c_{pm} = c_{vm} + R_m$ (Eq. \ref{e:SpecificHeatRelation}).

!!! tip "Implementation Note"
    The enthalpy is implemented in the [`enthalpy`](@ref) function. The relationship $h = I + R_m T$ is used for efficient computation, avoiding the need to compute individual constituent enthalpies.

The enthalpy is the relevant thermodynamic energy quantity in fluid transport. It arises in boundary conditions for energy fluxes and in the modeling of subgrid-scale (SGS) turbulent transport.

!!! example "Typical Values"
    For typical moist air at $T = 300$ K with $q_t = 0.01$, we have the specific enthalpy $h = 302.0 \times 10^3$ J/kg. This is significantly larger than the specific internal energy due to the $R_m T$ term.

## 8. Moist Static Energy

The sum of the specific enthalpy of moist air and the specific gravitational potential energy ``Φ`` is the moist static energy [Neelin1987](@cite)

```math
\begin{equation}\label{e:MSE}
\mathrm{MSE} = h + Φ.
\end{equation}
```

Moist static energy arises naturally as the static energy component that is transported in moist air. "Static" here refers to the fact that the (small) kinetic energy contribution to the total energy is neglected. The global integral of moist static energy is approximately conserved in adiabatic processes, even in the presence of reversible phase transitions and latent heat release. It is also approximately materially conserved ([Romps2015](@cite)).

## 9. Saturation Vapor Pressure

The Clausius-Clapeyron relation describes how the saturation vapor pressure $p_v^*$ of an ideal gas over a plane surface of condensate depends on temperature:

```math
\begin{equation}
\label{e:ClausiusClapeyron}
    \frac{d \log(p_v^*)}{dT} = \frac{L}{R_v T^2}.
\end{equation}
```

Here, $L$ is the specific latent heat of the phase transition, which is $L_v$ for the saturation vapor pressure over liquid, or $L_s$ for the saturation vapor pressure over ice.

Substituting the linear relation \eqref{e:LHTemperature} between latent heat and temperature, and taking $p_\mathrm{tr}$ to be the vapor pressure at the triple point (by definition equal to the saturation vapor pressures both over liquid and ice), the Clausius-Clapeyron relation can be integrated to give a closed-form expression for the saturation vapor pressure. This closed-form expression is known as the **Rankine-Kirchhoff approximation** ([Romps2021](@cite)), named after its independent discoverers. It is the unique solution consistent with our thermodynamic assumptions (constant heat capacities and ideal gas law):

```math
\begin{equation}
    p_v^* = p_{\mathrm{tr}} \left( \frac{T}{T_{\mathrm{tr}}} \right)^{\frac{\Delta c_p}{R_v}}
        \exp \left[ \frac{L_0 - \Delta c_p T_0}{R_v}
        \left( \frac{1}{T_{\mathrm{tr}}} - \frac{1}{T} \right) \right].
        \label{e:SatVaporPressure}
\end{equation}
```

!!! tip "Implementation Note"
    The saturation vapor pressure is implemented in the [`saturation_vapor_pressure`](@ref) function. The closed-form expression enables efficient computation without numerical integration.

With $L_0 = L_{v,0}$ or $L_0 = L_{s,0}$ and the corresponding heat capacity difference $\Delta c_p$, this gives saturation vapor pressures over liquid or ice that are accurate within 3% for temperatures between 200K and 330K ([Ambaum2020](@cite)). The accuracy of this approximation depends on the choice of thermodynamic constants; the values used in `Thermodynamics.jl` (specified in [ClimaParams.jl](https://github.com/CliMA/ClimaParams.jl)) are chosen to minimize errors in the Rankine-Kirchhoff approximation ([Yatunin2026](@cite)).

!!! example "Typical Values"
    At $T = 300$ K:

    | Quantity | Value |
    |----------|-------|
    | $p_v^*$ (liquid) | $3537$ Pa |
    | $p_v^*$ (ice)    | $286$ Pa  |

    At $T = 273.15$ K (triple point):

    | Quantity | Value |
    |----------|-------|
    | $p_v^*$ (liquid) = $p_v^*$ (ice) | $611$ Pa |

    The ratio of liquid to ice saturation vapor pressure at 300 K is approximately 12.4, reflecting the higher energy required for sublimation compared to vaporization.

To obtain the saturation vapor pressure over a mixture of liquid and ice (e.g., in mixed-phase clouds), using a weighted average of the relevant specific latent heats in the vapor pressure \eqref{e:SatVaporPressure} leads to a thermodynamically consistent formulation ([Pressel2015](@cite)). That is, if a fraction $\lambda_p$ of the condensate is liquid and the complement $1-\lambda_p$ is ice, calculating the saturation vapor pressure with a specific latent heat $\lambda_p L_v + (1-\lambda_p)L_s$ gives a thermodynamically consistent saturation vapor pressure over the mixture.

!!! note "Physical Interpretation"
    The weighted latent heat approach ensures thermodynamic consistency when computing saturation vapor pressure over mixed-phase conditions. This is important for modeling mixed-phase clouds where both liquid and ice coexist.

In phase equilibrium, the liquid fraction $\lambda_p(T)$ is a function of temperature alone. For temperatures above the freezing temperature $T_{\text{freeze}}$, $\lambda_p(T) = 1$. In strict phase equilibrium, $\lambda_p(T) = 0$ below freezing; that is, supercooled liquid does not exist. To parameterize the presence of supercooled liquid in phase equilibrium, a continuous function is used to represent a nonzero liquid fraction between the temperature of homogeneous ice nucleation $T_{\text{icenuc}}$ and the freezing temperature $T_{\text{freeze}}$. In `Thermodynamics.jl`, this is parameterized following [Kaul2015](@cite) as a power-law ramp interpolation, with some exponent $n$, such that:

```math
\lambda_p(T) = 
\begin{cases} 
0 & T \le T_{\text{icenuc}} \\
\left(\frac{T - T_{\text{icenuc}}}{T_{\text{freeze}} - T_{\text{icenuc}}}\right)^n & T_{\text{icenuc}} < T < T_{\text{freeze}} \\
1 & T \ge T_{\text{freeze}}
\end{cases}
```

This smooth transition represents the statistical presence of supercooled liquid in a phase equilibrium framework.

Out of phase equilibrium, when prognostic variables for liquid ($q_{liq}$) and ice ($q_{ice}$) specific humidities are available, the liquid fraction is defined simply as the ratio of liquid to total condensate:

```math
\lambda = \frac{q_{liq}}{q_{liq} + q_{ice}}
```

This definition makes no assumption about the temperature dependence of the phase partitioning and allows for phase non-equilibrium states (e.g., supercooled liquid at temperatures where the equilibrium function would force freezing).

!!! tip "Implementation Note"
    The [`liquid_fraction`](@ref) function in `Thermodynamics.jl` dispatches on the arguments provided.
    - [`liquid_fraction(param_set, T)`](@ref) computes the phase equilibrium temperature-dependent fraction.
    - [`liquid_fraction(param_set, T, q_liq, q_ice)`](@ref) computes the fraction from specific humidities. If no condensate is present (`q_liq + q_ice ≈ 0`), it falls back to a **slightly smoothed Heaviside function** (a linear ramp over $\pm 0.1$ K around freezing) to ensure differentiability of derived quantities such as saturation vapor pressure.

## 10. Saturation Specific Humidity

From the saturation vapor pressure $p_v^*$, the saturation specific humidity can be computed using the ideal gas law \eqref{e:eos}, giving the density of water vapor at saturation $\rho_v^* = p_v^*(T)/(R_v T)$, and hence the saturation specific humidity

```math
\begin{equation}
     q_v^* = \frac{\rho_v^*}{\rho} = \frac{p_v^*(T)}{\rho R_v T}.
\label{e:SatShum}
\end{equation}
```

## 11. Saturation Adjustment

The thermodynamic state of a moist air parcel in local phase equilibrium is uniquely defined by its density $\rho$, total specific humidity $q_t$, and internal energy $I$ (or temperature $T$ or other thermodynamic variables). Thus, a moist dynamical core that assumes phase equilibrium thermodynamics can be constructed from a dry dynamical core with internal (or total) energy as a prognostic variable by including only a tracer for the total specific humidity $q_t$, and calculating the temperature and condensate specific humidities from $\rho$, $q_t$, and $I$.

Obtaining the temperature and condensate specific humidities from the state variables $\rho$, $q_t$, and $I$ is the problem of finding the root $T$ of

```math
\begin{equation}
I^*(T; \rho, q_t) - I = 0,
\label{e:SatAdjustment}
\end{equation}
```

where ``I^*(T; \rho, q_t)`` is the internal energy at phase equilibrium. In an unsaturated equilibrium, there is no condensate, so ``I^*`` is the internal energy with ``q_l=q_i=0``. At saturation, the internal energy ``I^*`` depends on the vapor specific humidity, ``q_v = q_v^*(T, \rho)``, and on the saturation excess (total condensate)

```math
\begin{equation}
q_c^* = \max\bigl[q_t - q_v^*(T, \rho), 0\bigr],
\end{equation}
```

which is partitioned according to the liquid fraction ``λ_p`` into

```math
\begin{equation}
q_l^* = λ_p(T) q_c^* \quad \text{and} \quad q_i^* = \bigl[1-λ_p(T)\bigr]q_c^*.
\label{e:PhasePartition}
\end{equation}
```

In saturated conditions, finding the root of \eqref{e:SatAdjustment} is a nonlinear problem, which must be solved iteratively or approximately, in what is known as a saturation adjustment procedure.

A zeroth-order approximation of the temperature $T$ satisfying the saturation adjustment condition \eqref{e:SatAdjustment} is obtained by assuming unsaturated conditions. In that case, the expression \eqref{e:temperature} for temperature, with $q_l=q_i=0$, gives the unsaturated temperature

```math
\begin{equation}
    T_1 = T_0 + \frac{I - q_t I_{v,0}}{c_{vm}^*}.
\end{equation}
```

Here, the isochoric specific heat capacity in phase equilibrium, $c_{vm}^* = c_{vm}(q^*)$, is the specific heat capacity under phase equilibrium partitioning $q^*$ of the phases, which here, for unsaturated conditions, means $q^*=(q_t; q_l=0, q_i=0)$. If the total specific humidity $q_t$ is less than the saturation specific humidity at $T_1$ ($q_t \le q_v^*(T_1, \rho)$), the air is indeed unsaturated, and $T=T_1$ is the exact temperature consistent given $I$, $\rho$, and $q_t$.

If the air is saturated ($q_t > q_v^*(T_1, \rho)$), successively improved temperature estimates $T_{n+1}$ can be obtained from the temperature $T_n$ ($n=1,\dots$) by Newton's method. Linearizing the saturation internal energy $I^*(T; \rho, q_t)$ around the temperature $T_n$ gives

```math
\begin{equation}
    I^*(T; \rho, q_t) \approx I^*(T_n; \rho, q_t) + \left.\frac{\partial I^*(T; \rho; q_t)}{\partial T}\right|_{T_n} (T - T_n),
\end{equation}
```

and solving for the temperature $T$ gives the first-order Newton update

```math
\begin{equation}
    T_{n+1} = T_{n} - \frac{I^*(T_{n}; \rho, q_t) - I}{(\partial I^*/\partial T)|_{T_{n}}}.
\end{equation}
```

The derivative ``\partial I^*/\partial T|_{T_n}`` is obtained by differentiation of the internal energy \eqref{e:TotalInternalEnergy}. The implementation in `Thermodynamics.jl` includes the full derivative, including the temperature dependence of the liquid fraction $\lambda_p(T)$:

```math
\begin{equation}
     \left.\frac{\partial I^*}{\partial T}\right|_{T_n}
     = c_{vm}^* + (e_{vap} - e_{cond}) \left. \frac{\partial q_v^*}{\partial T}\right|_{T_n} + (e_{liq} - e_{ice}) q_{cond}^* \frac{\partial \lambda_p}{\partial T},
\end{equation}
```

where $q_{cond}^* = q_t - q_v^*$, and $e_{cond} = \lambda_p e_{liq} + (1-\lambda_p) e_{ice}$. The derivative of the saturation specific humidity is given by

```math
\begin{equation}
    \left. \frac{\partial q_v^*}{\partial T}\right|_{T_n} = q_v^*(T_n) \left( \frac{L}{R_v T_n^2} - \frac{1}{T_n} \right),
\end{equation}
```

which follows from differentiation of the ideal gas law for vapor and the Clausius-Clapeyron relation. Note the inclusion of the $-1/T$ term, which arises from the density dependence. The derivative of the liquid fraction $\partial \lambda_p / \partial T$ is non-zero in the supercooled liquid mixed-phase region.

The resulting successive Newton approximations $T_n$ generally converge quadratically. Because condensate specific humidities are usually small, $T_1$ provides a close initial estimate, and few iterations are needed. Even the first-order approximation $T\approx T_2$ often suffices. With the smoothed liquid fraction $\lambda_p$ ramp (see [`liquid_fraction`](@ref)), the derivative of $I^*$ with respect to temperature remains continuous across the phase transition, allowing the saturation adjustment to converge reliably without requiring special treatment or limiters.

Using saturation adjustment makes it possible to construct a moist dynamical core that has the total specific humidity $q_t$ as the only prognostic moisture variable. The price for this simplicity is the necessity to solve a nonlinear problem iteratively (or approximately) at each time step, and being confined to a phase equilibrium framework which cannot adequately account for non-equilibrium processes. Using explicit tracers for the condensates $q_l$ and $q_i$ in addition to $q_t$ avoids iterations at each time step and allows the inclusion of explicit non-equilibrium processes, such as those leading to the formation of supercooled liquid in mixed-phase clouds.

## 12. Auxiliary Thermodynamic Functions

Several auxiliary thermodynamic functions are commonly used.

### 12.1 Relative Humidity

The relative humidity is defined as the ratio of the partial pressure of water vapor $p_v$ to the saturation vapor pressure $p_v^*$,

```math
\mathrm{RH} = \frac{p_v}{p_v^*}.
```

Using the ideal gas law for water vapor, ``p_v = q_v \rho R_v T``, this can be written as

```math
\begin{equation}
    \mathrm{RH} = \frac{q_v \rho R_v T}{p_v^*},
\end{equation}
```

where ``p_v^*`` is the saturation vapor pressure \eqref{e:SatVaporPressure}. Over a mixture of ice and liquid, the saturation vapor pressure \eqref{e:SatVaporPressure} is evaluated with a specific latent heat ``L = λ_p L_v + (1-λ_p) L_s`` that is a weighted sum of those for vaporization and sublimation.

### 12.2 Potential Temperature

The potential temperature $\theta$ is the temperature an air mass would have if brought adiabatically from pressure $p$ and temperature $T$ to some reference pressure $p_0$ (typically taken to be mean sea level pressure):

```math
\begin{equation}
θ = \frac{T}{\Pi},
\label{e:PotTempPressT}
\end{equation}
```

where ``\Pi`` is known as the Exner function

```math
\begin{equation}
    \Pi  = \left( \frac{p}{p_0} \right)^κ \quad \text{with} \quad κ = \frac{R_m}{c_{pm}}.
\end{equation}
```

Note that the adiabatic exponent $\kappa$ takes the effect of  moisture on the effective gas "constant" and specific heat capacity of air into account.

### 12.3 Virtual Temperature and Virtual Potential Temperature

The virtual or density temperature $T_v$ is the temperature dry air would need to have to have the same density as moist air at the same pressure.
Using the ideal gas law $p/\rho = R_m T$, this implies $R_m T  = R_d T_v $, or

```math
\begin{equation} \label{e:virtual_temp}
T_v = \frac{R_m}{R_d} T.
\end{equation}
```

A virtual potential temperature can be defined analogously:

```math
\begin{equation}
θ_v = \frac{R_m}{R_d} θ.
\label{e:virtual_pottemp}
\end{equation}
```

!!! note "Relationship between virtual temperature and 'density temperature'"
    Some texts distinguish a "(condensate-ignoring) virtual temperature" and a "density temperature", and an analogous condensate-ignoring virtual potential temperature and density potential temperature.
    In those texts, the definition of density temperature incorporates condensate mass but their "condensate-ignoring virtual temperature" does not.
    We always take the mass of any condensate into account in the thermodynamics of moist air, so this distinction is irrelevant here.
    In other words, because the virtual temperature defined above incorporates the mass of condensate into ``R_m(q)`` (and virtual potential temperature is additionally defined in terms of ``c_{pm}(q)``) via the potential temperature exponent ``κ``), there is no distinction between our virtual temperature and a hypothetical "density temperature".

### 12.4 Liquid-Ice Potential Temperature

When the amount of condensate in air is small and the temperature $T$ is not too small (e.g., [Tripoli1981](@cite)), the (linearized) liquid-ice potential temperature,

```math
\begin{equation}\label{e:LiquidIcePottemp}
θ_{li} = θ \left( 1 - \frac{L_{v,0} q_l + L_{s, 0} q_i}{c_{pm} T} \right),
\end{equation}
```

is approximately materially conserved in adiabatic and reversible processes (including phase transitions). It is approximately the potential temperature \eqref{e:PotTempPressT} an air parcel would have if all liquid water in the parcel were evaporated and all ice sublimated. It is the limit of a more general expression for liquid-ice potential temperature for small $q_l$ and $q_i$ and taking the specific latent heats to be constant (e.g., [Bryan2004](@cite)). The liquid-ice potential temperature \eqref{e:LiquidIcePottemp} and variants thereof are sometimes used as variables in numerical models. We use it for diagnostic purposes, for comparison with other studies.

The liquid-ice potential temperature $\theta_{li}$ can be inverted for the temperature given pressure $p$ (and hence $\Pi$) and the specific humidities $q_t$, $q_l$, and $q_i$:

```math
\begin{equation}
    T = \Pi θ_{li} + \frac{L_{v, 0} q_l + L_{s, 0} q_i}{c_{pm}}.
\label{e:TempFromThetaLiGivenP}
\end{equation}
```

Alternatively, when density $\rho$ instead of pressure $p$ is given, the temperature can be obtained by Taylor expansion from the liquid-ice potential temperature $\theta_{li}$,

```math
\begin{equation}
T \approx T_u + \frac{L_{v,0} q_l + L_{s,0}q_i}{c_{vm}} - \frac{κ}{2} \frac{1}{T_u}\left(\frac{L_{v,0} q_l + L_{s,0} q_i}{c_{vm}}\right)^2
\end{equation}
\label{e:TempFromThetaLiGivenRho}
```

where

```math
\begin{equation}
   T_u =  \left( \frac{\rho R_m θ_{li}}{p_0} \right)^{R_m/c_{vm}} θ_{li}
\end{equation}
```

is the temperature that would correspond to $\theta_{li}$ in unsaturated conditions, i.e., when the condensate specific humidities $q_l$ and $q_i$ are zero. However, the specific heats $c_{vm}$ and $c_{pm}$ and the moist gas constant $R_m$ are evaluated with the given total and condensate specific humidities $q_t$, $q_l$, and $q_i$.

This expression for temperature as a function of liquid-ice potential temperature is obtained from \eqref{e:TempFromThetaLiGivenP} by substituting for pressure in the Exner function $\Pi$ from the ideal gas law, $p=\rho R_m T$, solving for temperature using a second-order Taylor expansion around $T_u$ for small condensate specific humidities, and using the relation $1-\kappa = c_{vm}/c_{pm}$, which follows from $c_{pm} - R_m = c_{vm}$. The relation for temperature \eqref{e:TempFromThetaLiGivenRho} holds to second order in condensate specific humidities $q_l$ and $q_i$. That is, the inversion relation \eqref{e:TempFromThetaLiGivenRho} holds to one higher order of accuracy than the definition of the liquid-ice potential temperature \eqref{e:LiquidIcePottemp} itself, which is only first-order accurate in the condensate specific humidities $q_l$ and $q_i$.

### 12.5 Speed of Sound

The speed of sound in (moist) unstratified air is

```math
\begin{equation}
 c_s = \left(\frac{c_{pm}}{c_{vm}} R_m T \right)^{1/2},
\label{e:soundspeed}
\end{equation}
```

with the appropriate gas constants for moist air. In the presence of stratification, additional terms arise ([Durran1999](@cite)).

!!! note "Physical Interpretation"
    The speed of sound represents the rate at which pressure perturbations propagate through the fluid. The expression accounts for the effect of moisture on the gas constant and specific heat capacities of the air mixture.

!!! example "Typical Values"
    For dry air at $T = 300$ K: $c_s \approx 347$ m/s
    For moist air with $q_t = 0.01$ at $T = 300$ K: $c_s ≈ 348$ m/s

!!! tip "Implementation Note"
    The speed of sound is implemented in the [`soundspeed_air`](@ref) function. It is used in numerical schemes that require the characteristic wave speeds for stability analysis.

## 13. Summary and Implementation Guidelines

### 13.1 Key Theoretical Framework

The thermodynamic framework presented here is based on a single fundamental approximation: **calorically perfect gases** with constant specific heat capacities. This approximation enables:

1. **Closed-form expressions** for all thermodynamic quantities
2. **Consistent thermodynamics** across all model components
3. **Computational efficiency** without numerical integration
4. **Accuracy within 1-3%** for atmospheric conditions

### 13.2 Core Equations

The most important equations for implementation are:

- **Equation of State**: $p = \rho R_m T$ with $R_m(q) = R_d (1 - q_t) + R_v q_v$
- **Internal Energy**: $I(T, q) = c_{vm}(q) (T - T_0) + q_v I_{v,0} - q_i I_{i,0} - (1 - q_t) R_d T_0$
- **Temperature Recovery**: $T = T_0 + (I - q_v I_{v,0} + q_i I_{i,0} + (1 - q_t) R_d T_0) / c_{vm}(q)$
- **Saturation Vapor Pressure**: Closed-form expression \eqref{e:SatVaporPressure}

### 13.3 Implementation Strategy

#### **Phase Equilibrium vs. Phase Non-Equilibrium Approaches**

1. **Phase Equilibrium Approach** (Saturation Adjustment):
    - Prognostic variables: $\rho$, $q_t$, $I$
    - Diagnostic variables: $T$, $q_l$, $q_i$
    - Requires iterative solution for temperature

2. **Phase Non-Equilibrium Approach**:
    - Prognostic variables: $\rho$, $q_t$, $q_l$, $q_i$, $I$
    - Diagnostic variable: $T$
    - Closed-form temperature recovery

#### **Numerical Considerations**

- **Saturation Adjustment**: Use Newton's method with analytical derivatives (a derivative-free secant method is also provided)
- **Temperature Recovery**: Monotonic function enables robust inversion
- **Mixed-Phase Conditions**: Weighted latent heats ensure consistency
- **Performance**: Use rearranged forms to avoid redundant computations

### 13.4 Validation and Testing

The implementation should be validated against:

1. **Energy Conservation**: Total energy should be conserved even in the presence of moist processes
2. **Temperature Recovery**: Inversion should recover original temperature
3. **Saturation Consistency**: Saturation conditions should be self-consistent
4. **Phase Transitions**: Latent heat release should match energy changes

!!! note "Testing Framework"
    Comprehensive testing is performed using [Tested Profiles](TestedProfiles.md)
    that cover the full range of atmospheric conditions.

`Thermodynamics.jl` contains a battery of unit tests that check for consistency of the thermodynamic formulation and convergence of saturation adjustement in a broad range of conditions.

### 13.5 Extensions and Limitations

#### **Current Limitations**

- Constant specific heat capacities (1-3% error)
- Neglected condensate volume

#### **Potential Extensions**

- Temperature-dependent heat capacities
- Non-ideal gases

!!! tip "Implementation Note"
    The `Thermodynamics.jl` package provides a complete implementation of this framework with comprehensive testing and validation. See the [API documentation](API.md) for detailed function signatures and usage examples.

!!! note "References"
    This formulation is described in [Yatunin2026](@cite), building on the theoretical framework described in [Romps2008](@cite), [Bott2008](@cite), and [Marquet2016](@cite), with additional developments for mixed-phase conditions from [Pressel2015](@cite).
