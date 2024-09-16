# Mathematical formulation

The thermodynamics of moist air is often subject to empirical approximations, which usually are opaque, internally inconsistent, and/or inconsistent across model components. For example, microphysical process models often use different approximations for thermodynamic quantities such as saturation vapor pressures than the dynamical core. The often bewildering array of approximations makes it difficult to achieve global conservation, e.g., of energy, and it complicates the use of models for other planetary atmospheres, with different thermodynamic parameters.

Here we introduce one consistent set of thermodynamic approximations for all model components. The key to thermodynamic consistency at reasonable accuracy is to take the specific heat capacities of the constituents of moist air (dry air, water vapor, liquid water, and ice) to be constant, i.e., to assume the gases to be calorically perfect. We discuss how to derive all other thermodynamic quantities that are needed on the basis of this one approximation ([Romps2008](@cite),[Marquet2015](@cite)). This includes:

 - Giving accurate and easily adaptable closed-form expressions for internal energies, enthalpies, specific latent heats, and saturation vapor pressures
 - Showing how to construct consistent sets of thermodynamic equations that either (i) assume equilibrium of the phases and require only one prognostic water variable, or (ii) do not assume equilibrium of the phases and require prognostic variables for all water phases
 - Showing how to obtain temperatures from energy variables under either phase equilibrium assumptions (by `saturation adjustment`) or phase non-equilibrium assumptions (by a closed-form expression for temperature).

The resulting thermodynamic functions are implemented in [Thermodynamics.jl](https://github.com/CliMA/Thermodynamics.jl).

## Working Fluid and Equation of State

The working fluid of the atmosphere model is moist, potentially cloudy air, considered to be an ideal mixture of dry air, water vapor, and condensed water (liquid and ice) in clouds. Dry air and water vapor are taken to be ideal gases. The specific volume of the cloud condensate is neglected relative to that of the gas phases (it is a  factor ``10^{3}`` less than that of the gas phases). All phases are assumed to have the same temperature, and are advected with the same velocity. The cloud condensates may be sedimenting relative to the gas phases, but slowly enough to be in thermal equilibrium with the surrounding fluid. However, the condensates do not need to be in thermodynamic equilibrium with the other fluid constituents; out-of-equilibrium phases such as supercooled liquid can exist. Falling condensate (precipitation) is not considered part of the working fluid because it generally cannot be assumed to be in thermal equilibrium with the surrounding fluid; it is treated separately.

The density of the moist air is denoted by ``ρ``. We use the following notation for the mass fractions of the moist air mixture (mass of a constituent divided by the total mass of the working fluid):

 - ``q_d``: dry air mass fraction
 - ``q_v``: water vapor specific humidity
 - ``q_l``: liquid water specific humidity
 - ``q_i``: ice specific humidity
 - ``q_c = q_l + q_i``: condensate specific humidity
 - ``q_t = q_v + q_c``: total specific humidity

Because this enumerates all constituents of the working fluid, we have ``q_t + q_d = 1``. In Earth's atmosphere, the water vapor specific humidity ``q_v`` generally dominates the total specific humidity ``q_t`` and is usually ``𝒪(10^{-2})`` or smaller; the condensate specific humidity is typically ``𝒪(10^{-4})``. Hence, water is a trace constituent of the atmosphere, and only a small fraction of atmospheric water is in condensed phases.

The pressure ``p`` of the working fluid is the sum of the partial pressures of dry air and water vapor, both taken to be ideal gases. Neglecting the volume of the condensed phases (but not their masses), this gives ``p = ρ (R_d q_d + R_v q_v) T``, where ``R_d`` is the specific gas constant of dry air, and ``R_v`` is the specific gas constant of water vapor.

This can also be written as

```math
\begin{equation}
    p = ρ R_m T,
\label{e:eos}
\end{equation}
```
where
```math
\begin{equation}
\begin{aligned}
    R_m(q) & = R_d (1 - q_t) + R_v q_v \\
        & = R_d \left[ 1 + (ε_{dv}-1)q_t - ε_{dv} q_c\right]
\end{aligned}
\label{e:Rm}
\end{equation}
```

is the specific gas ``constant`` of moist air (which is not a constant); here, we have used ``q_d = 1-q_t`` and ``q_v = q_t - q_c`` and introduced the ratio of the molar masses of dry air and water vapor ``ε_{dv} = R_v/R_d`` (≈ 1.61). Equations \eqref{e:eos} and \eqref{e:Rm} constitute the equation of state of the working fluid. We use the notation ``q=(q_t, q_l, q_i)`` for the tuple of specific humidities that determine the composition of moist air.

## Heat Capacities

The isochoric specific heat capacities of the constituents of moist air are:
 - ``c_{vd}``: Isochoric specific heat capacity of dry air;
 - ``c_{vv}``: Isochoric specific heat capacity of water vapor;
 - ``c_{vl}``: Isochoric specific heat capacity of liquid water;
 - ``c_{vi}``: Isochoric specific heat capacity of ice.

Our key thermodynamic approximation is to take these isochoric specific heat capacities to be constants, i.e., we take the gases to be calorically perfect. This is an approximation because they depend weakly on temperature. But for atmospheric conditions, the error of approximating them as constant is less than 1% for dry air, the main constituent of moist air, and at most a few percent for the water phases.

The difference between the isochoric and isobaric specific heat capacities is proportional to the specific volume. Consistent with taking the specific volume of liquid water and ice to be zero, we take the isochoric and isobaric specific heat capacities of the condensed phases to be equal. The isobaric specific heat capacities of the constituents then are:

 - ``c_{pd} = c_{vd} + R_d``: Isobaric specific heat capacity of dry air;
 - ``c_{pv} = c_{vv} + R_v``: Isobaric specific heat capacity of water vapor;
 - ``c_{pl} = c_{vl}``: Isobaric specific heat capacity of liquid water;
 - ``c_{pi} = c_{vi}``: Isobaric specific heat capacity of ice.

The corresponding specific heat capacities of moist air are the weighted sum of those of the constituents:
```math
\begin{equation}
\begin{aligned}
    c_{⋅ m}(q) & = (1-q_t) c_{⋅ d} + q_v c_{⋅ v} + q_l c_{⋅ l} + q_i c_{⋅ i} \\
    & = c_{⋅ d} + (c_{⋅ v} - c_{⋅ d})q_t + (c_{⋅ l} - c_{⋅ v})q_l + (c_{⋅ i} - c_{⋅ v})q_i
\end{aligned}
\label{e:SpecificHeat}
\end{equation}
```
where ``⋅`` stands for ``v`` or ``p`` and we have used ``q_v = q_t -q_l - q_i``. Straightforward substitution shows that the above relation between the specific heat capacities of the constituents also holds for moist air as a whole:

```math
\begin{equation}\label{e:SpecificHeatRelation}
    c_{pm}(q) = c_{vm}(q) + R_m(q).
\end{equation}
```
## Latent Heats

Kirchhoff's relation states that the specific latent enthalpy (heat) ``L`` of a phase change depends on temperature ``T`` through
```math
\begin{equation}
    \frac{dL}{dT} = Δ c_p,
\end{equation}
```
where ``Δ c_p`` is the difference in isobaric specific heat capacities between the phase with the higher and lower specific volume. For the constant isobaric specific heat capacities that we assume, this can be integrated to give
```math
\begin{equation}
    L(T) = L_0 + Δ c_p (T-T_0),
    \label{e:LHTemperature}
\end{equation}
```
where ``T_0`` is a reference temperature and ``L_0`` is the specific latent heat at ``T_0``.

For the phase transitions of water, this implies specifically:

 - ``L_v(T) = L_{v,0} + (c_{pv} - c_{pl}) (T - T_0)``: Latent heat of vaporization;
 - ``L_f(T) = L_{f,0} + (c_{pl} - c_{pi}) (T - T_0)``: Latent heat of fusion;
 - ``L_s(T) = L_{s,0} + (c_{pv} - c_{pi}) (T - T_0)``: Latent heat of sublimation.

With ``L_{s,0} = L_{v,0} + L_{f,0}``, this gives ``L_s(T) = L_v(T) + L_f(T)``, as it should.

## Internal Energies

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
Here, the reference specific internal energy ``I_{v,0}`` is the difference in specific internal energy between vapor and liquid at the reference temperature ``T_0``, and ``I_{i,0}`` is the difference in specific internal energy between ice and liquid at ``T_0``. The internal energy of moist air is the weighted sum of that of the constituents,
```math
\begin{equation}
\begin{aligned}
     I(T, q) & = (1-q_t) I_d(T) + q_v I_v(T) + q_l I_l(T) + q_i I_i(T)\\
          & = c_{vm}(q) (T - T_0)  + q_v I_{v,0} - q_i I_{i,0} - (1 - q_t) R_d T_0.
\end{aligned}
\label{e:totalInternalEnergy}
\end{equation}
```
The internal energy can be inverted to obtain the temperature given ``I`` and the specific humidities,
```math
\begin{equation}
    T = T_0 + \frac{I - (q_t - q_l - q_i) I_{v,0} + q_i I_{i,0} + (1 - q_t) R_d T_0}{c_{vm}(q)},
    \label{e:temperature}
\end{equation}
```
where we have used ``q_v = q_t - q_l - q_i``. This allows one to recover temperature given internal energy and specific humidities as state variables.

The reference specific internal energies ``I_{v,0}`` and ``I_{i,0}`` are related to the reference specific latent heats ``L_{v,0}`` and ``L_{f,0}``, which indicate the enthalpy differences between the phases at ``T_0``. The reference specific internal energies are obtained from the reference specific latent heats by subtracting the ``pV`` term, which is ``p_k/ρ_k`` for the relevant partial pressure ``p_k`` and specific volume ``1/ρ_k`` of the phase ``k`` (and hence is zero for the condensed phases, whose specific volume we neglect). This gives
```math
\begin{equation}
\begin{aligned}
     I_{v,0} &= L_{v, 0} - R_v T_0,\\
     I_{i,0} &= L_{f, 0}.
\end{aligned}
\label{e:RefInternalEnergies}
\end{equation}
```

## Enthalpies

The specific enthalpies of the constituents of moist air are obtained by adding ``p_k/ρ_k`` for phase ``k`` to the corresponding specific internal energy \eqref{e:InternalEnergies}. Again neglecting the specific volumes of the condensed phases and using the relations \eqref{e:RefInternalEnergies} between reference specific energies and latent heats, this gives:
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
where the last equality used ``c_{pm} = c_{vm} + R_m`` (Eq. \ref{e:SpecificHeatRelation}). The enthalpy is the relevant thermodynamic energy quantity in fluid transport. It arises in boundary conditions for energy fluxes and in the modeling of subgrid-scale (SGS) turbulent transport. For those purposes, we need gradients of the enthalpy, which can be written as
```math
\begin{equation}
\label{e:EnthalpyGradient}
    ∇h = c_{pm}(q) ∇T - h_d(T) ∇q_t
    + h_v(T) ∇q_v + h_l(T) ∇q_l + h_i(T) ∇q_i.
\end{equation}
```
This cleanly separates gradients involving temperature and gradients involving specific humidities.

## Moist Static Energy

The sum of the specific enthalpy of moist air and the specific gravitational potential energy ``Φ`` is the moist static energy [Neelin1987](@cite)
```math
\begin{equation}\label{e:MSE}
\mathrm{MSE} = h + Φ.
\end{equation}
```

Moist static energy arises naturally as the static energy component that is transported in moist air. ``Static`` here refers to the fact that the (small) kinetic energy contribution to the total energy is neglected (see section TODO). The global integral of moist static energy is approximately conserved in adiabatic processes, even in the presence of reversible phase transitions and latent heat release. It is also approximately materially conserved ([Romps2015](@cite)).

## Saturation Vapor Pressure

The Clausius-Clapeyron relation describes how the saturation vapor pressure ``p_v^*`` of an ideal gas over a plane surface of condensate depends on temperature:
```math
\begin{equation}
\label{e:Clausius_Clapeyron}
    \frac{d \log(p_v^*)}{dT} = \frac{L}{R_v T^2}.
\end{equation}
```

Here, ``L`` is the specific latent heat of the phase transition, which is ``L_v`` for the saturation vapor pressure over liquid, or ``L_s`` for the saturation vapor pressure over ice. Substituting the linear relation \eqref{e:LHTemperature} between latent heat and temperature, and taking ``p_\mathrm{tr}`` to be the vapor pressure at the triple point (by definition equal to the saturation vapor pressures both over liquid and ice), the Clausius-Clapeyron relation can be integrated to give a closed-form expression for the vapor pressure that is consistent with our thermodynamic assumptions:
```math
\begin{equation}
    p_v^* = p_{\mathrm{tr}} \left( \frac{T}{T_{\mathrm{tr}}} \right)^{\frac{\Delta c_p}{R_v}}
        \exp \left[ \frac{L_0 - \Delta c_p T_0}{R_v}
        \left( \frac{1}{T_{\mathrm{tr}}} - \frac{1}{T} \right) \right].
        \label{e:SatVaporPressure}
\end{equation}
```

With ``L_0 = L_{v,0}`` or ``L_0 = L_{s,0}`` and the corresponding heat capacity difference ``\Delta c_p``, this gives saturation vapor pressures over liquid or ice that are accurate within 3% for temperatures between 200K and 330K (with accuracy better than 1% for typical near-surface conditions).

To obtain the saturation vapor pressure over a mixture of liquid and ice (e.g., in mixed-phase clouds), using a weighted average of the relevant specific latent heats in the vapor pressure \eqref{e:SatVaporPressure} leads to a thermodynamically consistent formulation ([Pressel2015](@cite)). That is, if a fraction ``λ_p`` of the condensate is liquid and the complement ``1-λ_p`` is ice, calculating the saturation vapor pressure with a specific latent heat ``λ_p L_v + (1-λ_p)L_s`` gives a thermodynamically consistent saturation vapor pressure over the mixture.

In thermodynamic equilibrium, the liquid fraction
```math
\begin{equation}\label{e:liquid_fraction}
    λ_p(T) = \mathcal{H}(T-T_{\mathrm{freeze}})
\end{equation}
```
is a Heaviside function ``\mathcal{H}`` of temperature, being 0 below the freezing temperature ``T_{\mathrm{freeze}}`` and 1 above it. However, out of thermodynamic equilibrium, supercooled liquid can exist between the temperature of homogeneous ice nucleation ``T_{\mathrm{icenuc}}`` and the freezing temperature ``T_{\mathrm{freeze}}``. In most climate models, this is modeled by a continuous function
```math
\begin{equation}
    λ_p(T) =
    \begin{cases}
    0 & \text{for } T\le T_{\mathrm{icenuc}}\\
    0<λ_i(T)<1 & \text{for } T_{\mathrm{icenuc}} < T <  T_{\mathrm{freeze}}\\
    1   & \text{for } T\ge T_{\mathrm{freeze}},
    \end{cases}
\end{equation}
```
where ``λ_i`` interpolates between 0 at the temperature of homogeneous ice nucleation and 1 at the freezing temperature. However, it is important to recognize that this is merely an attempt to model out-of-equilibrium phases such as supercooled liquid within a thermodynamic equilibrium framework (where phase partitioning only depends on thermodynamic state variables but not on the history of air masses); this is not generally possible, and we will adopt alternative approaches.

## Saturation Specific Humidity

From the saturation vapor pressure ``p_v^*``, the saturation specific humidity can be computed using the ideal gas law \eqref{e:eos}, giving the density of water vapor at saturation ``ρ_v^* = p_v^*(T)/(R_v T)``, and hence the saturation specific humidity
```math
\begin{equation}
     q_v^* = \frac{ρ_v^*}{ρ} = \frac{p_v^*(T)}{ρ R_v T}.
\label{e:sat_shum}
\end{equation}
```

## Saturation Adjustment

Gibbs' phase rule states that in thermodynamic equilibrium, the temperature ``T`` and liquid and ice specific humidities ``q_l`` and ``q_i`` can be obtained from the three thermodynamic state variables density ``ρ``, total water specific humidity ``q_t``, and internal energy ``I``. Thus, a moist dynamical core that assumes equilibrium thermodynamics can be constructed from a dry dynamical core with total energy as a prognostic variable by including only a tracer for the total specific humidity ``q_t``, and calculating the temperature and condensate specific humidities from ``ρ``, ``q_t``, and ``I``.

Obtaining the temperature and condensate specific humidities from the state variables ``ρ``, ``q_t``, and ``I`` is the problem of finding the root ``T`` of
```math
\begin{equation}
I^*(T; ρ, q_t) - I = 0,
\label{e:SatAdjustment}
\end{equation}
```
where ``I^*(T; ρ, q_t)`` is the internal energy at equilibrium. In an unsaturated equilibrium, there is no condensate, so ``I^*`` is the internal energy with ``q_l=q_i=0``. At saturation, the internal energy ``I^*`` depends on the vapor specific humidity, ``q_v = q_v^*(T, ρ)``, and on the saturation excess (total condensate)
```math
\begin{equation}
q_c^* = \max\bigl[q_t - q_v^*(T, ρ), 0\bigr],
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

A zeroth-order approximation of the temperature ``T`` satisfying the saturation adjustment condition \eqref{e:SatAdjustment} is obtained by assuming unsaturated conditions. In that case, the expression \eqref{e:temperature} for temperature, with ``q_l=q_i=0``, gives the unsaturated temperature
```math
\begin{equation}
    T_1 = T_0 + \frac{I - q_t I_{v,0}}{c_{vm}^*}.
\end{equation}
```
Here, the isochoric specific heat capacity in equilibrium, ``c_{vm}^* = c_{vm}(q^*)``, is the specific heat capacity under equilibrium partitioning ``q^*`` of the phases, which here, for unsaturated conditions, means ``q^*=(q_t; q_l=0, q_i=0)``. If the total specific humidity ``q_t`` is less than the saturation specific humidity at ``T_1`` (``q_t \le q_v^*(T_1, ρ)``), the air is indeed unsaturated, and ``T=T_1`` is the exact temperature consistent given ``I``, ``ρ``, and ``q_t``.

If the air is saturated (``q_t > q_v^*(T_1, ρ)``), successively improved temperature estimates ``T_{n+1}`` can be obtained from the temperature ``T_n`` (``n=1,\dots``) by Newton's method, with analytical gradients. Linearizing the saturation internal energy ``I^*(T; ρ, q_t)`` around the temperature ``T_n`` gives
```math
\begin{equation}
    I^*(T; ρ, q_t) \approx I^*(T_n; ρ, q_t) + \left.\frac{\partial I^*(T; ρ; q_t)}{\partial T}\right|_{T_n} (T - T_n),
\end{equation}
```
and solving for the temperature ``T`` gives the first-order Newton update
```math
\begin{equation}
    T_{n+1} = T_{n} - \frac{I^*(T_{n}; ρ, q_t) - I}{(\partial I^*/\partial T)|_{T_{n}}}.
\end{equation}
```
The derivative ``\partial I^*/\partial T|_{T_n}`` is obtained by differentiation of the internal energy \eqref{e:total_internal_energy}, \hl{[add derivatives of phase partitioning function and make that function smooth]}
```math
\begin{multline}
     \left.\frac{\partial I^*(T; ρ, q_t)}{\partial T}\right|_{T_n}
     = c_{vm}^*(q_t, T_n) \\
     +  \left( I_{v,0} + [1-λ_p(T_n)]I_{i,0} + (T_n - T_0) \left. \frac{dc_{vm}^*}{dq_v^*}\right|_{T_n} \right) \left. \frac{\partial q_v^*(T; ρ, q_t)}{\partial T}\right|_{T_n},
\end{multline}
```
where ``c_{vm}^*(q_t, T_n) = c_{vm}[q^*(T_n)]`` is the isochoric specific heat in equilibrium at temperature ``T_n``, with ``q_v = q_v^*(T_n)`` and with the corresponding phase partitioning ``q^* = (q_t, q_l^*, q_i^*)`` according to \eqref{e:PhasePartition}. The derivative of the saturation specific humidity, ``\partial q_v^*(T;ρ, q_t)/\partial T``, is to be taken at a fixed density ``ρ`` and total specific humidity ``q_t``, like the other derivatives. We have neglected the singular derivative of ``λ_p`` at the freezing temperature ``T_{\mathrm{freeze}}``. The two remaining derivatives are that of the isochoric specific heat,
```math
\begin{equation}
    \left. \frac{dc_{vm}^*}{dq_v^*}\right|_{T_n} = c_{vv} - λ_p(T_n) c_{vl} - [1-λ_p(T_n)]c_{vi},
\end{equation}
```
obtained from the definition \eqref{e:SpecificHeat} of the specific heat of moist air, and that of the saturation specific humidity,
```math
\begin{equation}
    \left. \frac{\partial q_v^*(T; ρ, q_t)}{\partial T}\right|_{T_n} = q_v^*(T_n) \frac{L}{R_v T_n^2} \quad \text{with} \quad L = λ_p(T_n) L_v + [1-λ_p(T_n)] L_s,
\end{equation}
```
obtained from the Clausius-Clapeyron relation \eqref{e:Clausius_Clapeyron} and the relation \eqref{e:sat_shum} between specific humidity and vapor pressure.

The resulting successive Newton approximations ``T_n`` generally converge quadratically. Because condensate specific humidities are usually small, ``T_1`` provides a close initial estimate, and few iterations are needed. Even the first-order approximation ``T\approx T_2`` often suffices. However, convergence may not be achieved near the phase transition at the freezing temperature ``T_{\mathrm{freeze}}`` because the derivative of ``I^*`` with respect to temperature is discontinuous there. In that case, the number of iterations needs to be limited (2--3 iterations generally suffice).

Using saturation adjustment makes it possible to construct a moist dynamical core that has the total specific humidity ``q_t`` as the only prognostic moisture variable. The price for this simplicity is the necessity to solve a nonlinear problem iteratively (or approximately) at each time step, and being confined to an equilibrium thermodynamics framework which cannot adequately account for non-equilibrium processes. Using explicit tracers for the condensates ``q_l`` and ``q_i`` in addition to ``q_t`` avoids iterations at each time step and allows the inclusion of explicit non-equilibrium processes, such as those leading to the formation of supercooled liquid in mixed-phase clouds.

## Auxiliary Thermodynamic Functions

Several auxiliary thermodynamic functions are commonly used.

### Relative Humidity

The relative humidity is defined as the ratio of the partial pressure of water vapor ``p_v`` to the saturation vapor pressure ``p_v^*``,
```math
\mathrm{RH} = \frac{p_v}{p_v^*}.
```
Using the ideal gas law for water vapor, ``p_v = q_v ρ R_v T``, this can be written as
```math
\begin{equation}
    \mathrm{RH} = \frac{q_v ρ R_v T}{p_v^*},
\end{equation}
```
where ``p_v^*`` is the saturation vapor pressure \eqref{e:SatVaporPressure}. Over a mixture of ice and liquid, the saturation vapor pressure \eqref{e:SatVaporPressure} is evaluated with a specific latent heat ``L = λ_p L_v + (1-λ_p) L_s`` that is a weighted sum of those for vaporization and sublimation.

### Potential Temperature

The potential temperature ``θ`` is the temperature an air mass would have if brought adiabatically from pressure ``p`` and temperature ``T`` to some reference pressure ``p_0`` (typically taken to be mean sea level pressure):
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
Note that the adiabatic exponent ``κ`` takes the effect of  moisture on the effective gas ``constant`` and specific heat capacity of air into account.

### Virtual Temperature and Virtual Potential Temperature

The virtual or density temperature ``T_v`` is the temperature dry air would need to have to have the same density as moist air at the same pressure.
Using the ideal gas law ``p/ρ = R_m T``, this implies ``R_m T  = R_d T_v ``, or

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
    In other words, because the virtual temperature defined above incorporates the mass of condensate into ``R_m(q)`` (and virtual potential temperature is additionally defined in terms of ``c_pm(q)``) via the potential temperature exponent ``κ``), there is no distinction between our virtual temperature and a hypothetical "density temperature".

### Liquid-Ice Potential Temperature

When the amount of condensate in air is small and the temperature ``T`` is not too small (e.g., [Tripoli1981](@cite)), the (linearized) liquid-ice potential temperature,
```math
\begin{equation}\label{e:liquid_ice_pottemp}
θ_{li} = θ \left( 1 - \frac{L_{v,0} q_l + L_{s, 0} q_i}{c_{pm} T} \right),
\end{equation}
```
is approximately materially conserved in adiabatic and reversible processes (including phase transitions). It is approximately the potential temperature \eqref{e:PotTempPressT} an air parcel would have if all liquid water in the parcel were evaporated and all ice sublimated. It is the limit of a more general expression for liquid-ice potential temperature for small ``q_l`` and ``q_i`` and taking the specific latent heats to be constant (e.g., [Bryan2004](@cite)). The liquid-ice potential temperature \eqref{e:liquid_ice_pottemp} and variants thereof are sometimes used as variables in numerical models. We use it for diagnostic purposes, for comparison with other studies.

The liquid-ice potential temperature ``θ_{li}`` can be inverted for the temperature given pressure ``p`` (and hence ``\Pi``) and the specific humidities ``q_t``, ``q_l``, and ``q_i``:
```math
\begin{equation}
    T = \Pi θ_{li} + \frac{L_{v, 0} q_l + L_{s, 0} q_i}{c_{pm}}.
\label{e:TempFromThetaLiGivenP}
\end{equation}
```

Alternatively, when density ``ρ`` instead of pressure ``p`` is given, the temperature can be obtained by Taylor expansion from the liquid-ice potential temperature ``θ_{li}``,
```math
\begin{equation}
T \approx T_u + \frac{L_{v,0} q_l + L_{s,0}q_i}{c_{vm}} - \frac{κ}{2} \frac{1}{T_u}\left(\frac{L_{v,0} q_l + L_{s,0} q_i}{c_{vm}}\right)^2
\end{equation}
\label{e:TempFromThetaLiGivenRho}
```
where
```math
\begin{equation}
   T_u =  \left( \frac{ρ R_m θ_{li}}{p_0} \right)^{R_m/c_{vm}} θ_{li}
\end{equation}
```
is the temperature that would correspond to ``θ_{li}`` in unsaturated conditions, i.e., when the condensate specific humidities ``q_l`` and ``q_i`` are zero. However, the specific heats ``c_{vm}`` and ``c_{pm}`` and the moist gas constant ``R_m`` are evaluated with the given total and condensate specific humidities ``q_t``, ``q_l``, and ``q_i``.

TODO: there was a commented equation here, do we need it?

This expression for temperature as a function of liquid-ice potential temperature is obtained from \eqref{e:TempFromThetaLiGivenP} by substituting for pressure in the Exner function ``\Pi`` from the ideal gas law, ``p=ρ R_m T``, solving for temperature using a second-order Taylor expansion around ``T_u`` for small condensate specific humidities, and using the relation ``1-κ = c_{vm}/c_{pm}``, which follows from ``c_{pm} - R_m = c_{vm}``. The relation for temperature \eqref{e:TempFromThetaLiGivenRho} holds to second order in condensate specific humidities ``q_l`` and ``q_i``. That is, the inversion relation \eqref{e:TempFromThetaLiGivenRho} holds to one higher order of accuracy than the definition of the liquid-ice potential temperature \eqref{e:liquid_ice_pottemp} itself, which is only first-order accurate in the condensate specific humidities ``q_l`` and ``q_i``.

### Speed of Sound
The speed of sound in (moist) unstratified air is
```math
\begin{equation}
 c_s = \left(\frac{c_{pm}}{c_{vm}} R_m T \right)^{1/2},
\label{e:soundspeed}
\end{equation}
```
with the appropriate gas constants for moist air. In the presence of stratification, additional terms arise ([Durran1999](@cite)).
