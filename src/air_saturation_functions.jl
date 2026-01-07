export liquid_fraction
export has_condensate
export saturation_vapor_pressure
export q_vap_saturation
export q_vap_saturation_from_pressure
export supersaturation
export saturation_excess
export condensate_partition
export vapor_pressure_deficit

"""
    liquid_fraction(param_set, T)
    liquid_fraction(param_set, T, q_liq, q_ice)

The fraction of condensate that is liquid.

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T`: temperature [K]
 - (optional) `q_liq`: liquid specific humidity [kg/kg]
 - (optional) `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `λ`: liquid fraction [dimensionless], 0 ≤ λ ≤ 1

If `q_liq` and `q_ice` are provided, the liquid fraction is computed from them.
If `q_liq + q_ice` exceeds a small threshold (see [`has_condensate`](@ref)), `q_liq / (q_liq + q_ice)`
is returned. If there is effectively no condensate, a smooth temperature-dependent partitioning is used
(linear ramp from 0 to 1 over ±0.1 K around freezing).

If `q_liq` and `q_ice` are not provided, the liquid fraction is computed from
temperature using a power law interpolation between `T_icenuc` and `T_freeze`.

Edge cases:
- For `T > T_freeze`, this returns `1`; for `T ≤ T_icenuc`, it returns `0`.
- The temperature-only form uses a (generally broader) supercooled-liquid transition between `T_icenuc` and `T_freeze`,
  whereas the `(T, q_liq, q_ice)` form uses the narrow ±0.1 K transition *only* when `q_liq ≈ q_ice ≈ 0`.

# Reference
Kaul et al. (2015), "Sensitivities in large-eddy simulations of mixed-phase Arctic stratocumulus
clouds using a simple microphysics approach," *Monthly Weather Review*, **143**, 4393-4421,
doi:[10.1175/MWR-D-14-00319.1](https://doi.org/10.1175/MWR-D-14-00319.1).
"""
@inline function liquid_fraction(param_set::APS, T, q_liq, q_ice)
    FT = eltype(param_set)
    q_c = condensate_specific_humidity(q_liq, q_ice)

    # If no condensate, use sharp temperature dependent partitioning
    Tᶠ = TP.T_freeze(param_set)
    ΔT = FT(0.1) # Smooth over +/- 0.1 K
    # Linear ramp from 0 to 1 over [Tᶠ - ΔT, Tᶠ + ΔT]
    λ_T = clamp((T - (Tᶠ - ΔT)) / (2 * ΔT), zero(T), one(T))

    return ifelse(has_condensate(q_c), q_liq / q_c, λ_T)
end

@inline function liquid_fraction(param_set::APS, T)
    FT = eltype(param_set)

    # Interpolation between homogeneous nucleation and freezing temperatures
    Tᶠ = TP.T_freeze(param_set)   # freezing temperature
    Tⁱ = TP.T_icenuc(param_set)   # temperature of homogeneous ice nucleation
    n = TP.pow_icenuc(param_set)  # power law partial ice nucleation parameter
    λᵖ = ((T - Tⁱ) / (Tᶠ - Tⁱ))^n

    above_freezing = T > Tᶠ
    supercooled_liquid = (T ≤ Tᶠ) & (T > Tⁱ)

    return ifelse(
        above_freezing,
        one(T),
        ifelse(supercooled_liquid, λᵖ, zero(T)),
    )
end

"""
    has_condensate(q_c)

Bool indicating if condensate exists, i.e., q_c > eps.

We use a threshold of `eps` rather than `0` to avoid division by zero in functions
like `liquid_fraction` and to robustly handle numerical noise.
"""
@inline has_condensate(q_c) = q_c > ϵ_numerics(typeof(q_c))

"""
    saturation_vapor_pressure(param_set, T, ::Liquid)
    saturation_vapor_pressure(param_set, T, ::Ice)

The saturation vapor pressure over a plane surface of condensate.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `Liquid()` or `Ice()` to dispatch over the condensate type

# Returns
 - `p_v^*`: saturation vapor pressure [Pa]

The saturation vapor pressure is computed by integration of the Clausius-Clapeyron
relation, assuming constant specific heat capacities in the so-called Rankine-Kirchhoff
approximation.
"""
@inline function saturation_vapor_pressure(param_set::APS, T, ::Liquid)
    LH_v0 = TP.LH_v0(param_set)
    cp_v = TP.cp_v(param_set)
    cp_l = TP.cp_l(param_set)
    return saturation_vapor_pressure_calc(param_set, T, LH_v0, cp_v - cp_l)
end

@inline function saturation_vapor_pressure(param_set::APS, T, ::Ice)
    LH_s0 = TP.LH_s0(param_set)
    cp_v = TP.cp_v(param_set)
    cp_i = TP.cp_i(param_set)
    return saturation_vapor_pressure_calc(param_set, T, LH_s0, cp_v - cp_i)
end

"""
    saturation_vapor_pressure_calc(param_set, T, LH_0, Δcp)

Internal function. Computes the saturation vapor pressure using the Rankine-Kirchhoff
approximation.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `LH_0`: latent heat at reference temperature `T_0` [J/kg]
 - `Δcp`: difference in isobaric specific heat capacity between the two phases [J/(kg·K)]

# Returns
 - `p_v^*`: saturation vapor pressure [Pa]

The computed value is:

``p_v^*(T) = p_{tr} \\left( \\frac{T}{T_{tr}} \\right)^{\\Delta c_{p} / R_v} \\exp \\left[ \\frac{L_0 - \\Delta c_{p} T_0}{R_v} \\left( \\frac{1}{T_{tr}} - \\frac{1}{T} \\right) \\right]``

where ``T_{tr}`` is the triple point temperature, ``p_{tr}`` is the triple point pressure,
``T_0`` is the reference temperature, and ``R_v`` is the gas constant for water vapor.
"""
@inline function saturation_vapor_pressure_calc(param_set::APS, T, LH_0, Δcp)
    press_triple = TP.press_triple(param_set)
    R_v = TP.R_v(param_set)
    T_triple = TP.T_triple(param_set)
    T_0 = TP.T_0(param_set)

    return press_triple *
           # Use fast_power to compute (T / T_triple)^(Δcp / R_v) more efficiently
           fast_power(T / T_triple, Δcp / R_v) *
           exp((LH_0 - Δcp * T_0) / R_v * (1 / T_triple - 1 / T))
end

# Promote the arguments to a common type to allow AD with dual numbers
saturation_vapor_pressure_calc(param_set, T, LH_0, Δcp) =
    saturation_vapor_pressure_calc(param_set, promote(T, LH_0, Δcp)...)

"""
    saturation_vapor_pressure(param_set, T)
    saturation_vapor_pressure(param_set, T, q_liq, q_ice)

The saturation vapor pressure over liquid, ice, or a mixture of liquid and ice.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `q_liq`: (optional) liquid specific humidity [kg/kg]
 - `q_ice`: (optional) ice specific humidity [kg/kg]

# Returns
 - `p_v^*`: saturation vapor pressure [Pa]

If `q_liq` and `q_ice` are provided, the saturation vapor pressure is computed
from a weighted mean of the latent heats of vaporization and sublimation, with
the weights given by the liquid fraction (see [`liquid_fraction`](@ref)).
If `q_liq` and `q_ice` are 0, the saturation vapor pressure is that over liquid
above freezing and over ice below freezing.

Otherwise, the liquid fraction below freezing is computed from a temperature dependent
parameterization `liquid_fraction(param_set, T)`.

Edge case: the `saturation_vapor_pressure(param_set, T, q_liq, q_ice)` form includes a small smooth transition
around freezing when `q_liq == q_ice == 0` (via `liquid_fraction(param_set, T, q_liq, q_ice)`).
"""
@inline function saturation_vapor_pressure(param_set::APS, T, q_liq, q_ice)
    λ = liquid_fraction(param_set, T, q_liq, q_ice)
    return saturation_vapor_pressure_mixture(param_set, T, λ)
end

@inline function saturation_vapor_pressure(param_set::APS, T)
    λ = liquid_fraction(param_set, T)
    return saturation_vapor_pressure_mixture(param_set, T, λ)
end

"""
    saturation_vapor_pressure_mixture(param_set, T, λ)

Internal function. Compute the saturation vapor pressure over a mixture of liquid
and/or ice, given the temperature `T` and liquid fraction `λ`.

The computation uses a weighted mean of the temperature-dependent latent heats of
vaporization and sublimation, weighted by the liquid fraction, following Pressel 
et al. (JAMES, 2015).
"""
@inline function saturation_vapor_pressure_mixture(param_set::APS, T, λ)
    LH_v0 = TP.LH_v0(param_set)
    LH_s0 = TP.LH_s0(param_set)
    cp_v = TP.cp_v(param_set)
    cp_l = TP.cp_l(param_set)
    cp_i = TP.cp_i(param_set)

    # Effective latent heat at T_0 and effective difference in isobaric specific
    # heats of the mixture, which combined yield the effective (temperature-dependent)
    # latent heat
    LH_0 = λ * LH_v0 + (1 - λ) * LH_s0
    Δcp = λ * (cp_v - cp_l) + (1 - λ) * (cp_v - cp_i)

    # Saturation vapor pressure over possible mixture of liquid and ice
    return saturation_vapor_pressure_calc(param_set, T, LH_0, Δcp)
end

"""
    q_vap_saturation(param_set, T, ρ)
    q_vap_saturation(param_set, T, ρ, q_liq, q_ice)

The saturation specific humidity.

# Arguments
- `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
- `T`: temperature [K]
- `ρ`: air density [kg/m³]
- (optional) `q_liq`: liquid specific humidity [kg/kg]
- (optional) `q_ice`: ice specific humidity [kg/kg]

# Returns
- `q_v^*`: saturation specific humidity [kg/kg]

If `q_liq` and `q_ice` are provided, the saturation specific humidity is that over a
mixture of liquid and ice, computed in a thermodynamically consistent way from the
weighted sum of the latent heats of the respective phase transitions. That is, the 
saturation vapor pressure and from it the saturation specific humidity are computed 
from a weighted mean of the latent heats of vaporization and sublimation, with the 
weights given by the liquid fraction (see [`liquid_fraction`](@ref)). If `q_liq` and 
`q_ice` are 0, the saturation specific humidity is that over liquid above freezing and 
over ice below freezing.
 
Otherwise, the fraction of liquid is given by the temperature dependent
`liquid_fraction(param_set, T)`.

# Reference
Pressel et al. (2015), "Numerics and subgrid-scale modeling in large eddy simulations of
stratocumulus clouds," *Journal of Advances in Modeling Earth Systems*, **7**(3), 1199-1220,
doi:[10.1002/2014MS000376](https://doi.org/10.1002/2014MS000376).
"""
@inline function q_vap_saturation(param_set::APS, T, ρ, q_liq, q_ice)
    p_v_sat = saturation_vapor_pressure(param_set, T, q_liq, q_ice)
    return q_vap_from_p_vap(param_set, T, ρ, p_v_sat)
end

@inline function q_vap_saturation(param_set::APS, T, ρ)
    p_v_sat = saturation_vapor_pressure(param_set, T)
    return q_vap_from_p_vap(param_set, T, ρ, p_v_sat)
end

"""
    q_vap_saturation(param_set, T, ρ, phase::Phase)

The saturation specific humidity.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `ρ`: (moist-)air density [kg/m³]
 - `phase`: the phase to compute saturation over (either `Liquid()` or `Ice()`)

# Returns
 - `q_v^*`: saturation specific humidity [kg/kg]
"""
@inline function q_vap_saturation(param_set::APS, T, ρ, phase::Phase)
    p_v_sat = saturation_vapor_pressure(param_set, T, phase)
    return q_vap_from_p_vap(param_set, T, ρ, p_v_sat)
end

"""
    q_vap_saturation_from_pressure(param_set, q_tot, p, T)
    q_vap_saturation_from_pressure(param_set, q_tot, p, T, q_liq, q_ice)

The saturation specific humidity.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `q_tot`: total water specific humidity [kg/kg]
 - `p`: air pressure [Pa]
 - `T`: air temperature [K]
 - (optional) `q_liq`: liquid specific humidity [kg/kg]
 - (optional) `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `q_v^*`: saturation specific humidity [kg/kg]

If `q_liq` and `q_ice` are provided, the saturation vapor pressure is computed
from a weighted mean of the latent heats of vaporization and sublimation, with
the weights given by the liquid fraction (see [`liquid_fraction`](@ref)).
If `q_liq` and `q_ice` are 0, the saturation vapor pressure is that over liquid
above freezing and over ice below freezing.

Otherwise, the liquid fraction is computed from a temperature dependent parameterization
`liquid_fraction(param_set, T)`.

The saturation specific humidity is computed as:
``q_v^* = (R_d / R_v) * (1 - q_{tot}) * p_v^* / (p - p_v^*)``
where `p_v^*` is the saturation vapor pressure.

Edge case: this expression assumes `p > p_v^*(T)`; if `p ≤ p_v^*(T)` the denominator changes sign.
"""
@inline function q_vap_saturation_from_pressure(
    param_set::APS,
    q_tot,
    p,
    T,
    q_liq,
    q_ice,
)
    p_v_sat = saturation_vapor_pressure(param_set, T, q_liq, q_ice)
    return q_vap_saturation_from_pressure_calc(param_set, q_tot, p, p_v_sat)
end

@inline function q_vap_saturation_from_pressure(param_set::APS, q_tot, p, T)
    p_v_sat = saturation_vapor_pressure(param_set, T)
    return q_vap_saturation_from_pressure_calc(param_set, q_tot, p, p_v_sat)
end

"""
    q_vap_saturation_from_pressure_calc(param_set, q_tot, p, p_v_sat)

Internal function. Compute the saturation specific humidity from the saturation vapor pressure.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `q_tot`: total specific humidity [kg/kg]
 - `p`: air pressure [Pa]
 - `p_v_sat`: saturation vapor pressure [Pa]

# Returns
 - `q_v^*`: saturation specific humidity [kg/kg]

The saturation specific humidity is computed as:
``q_v^* = (R_d / R_v) * (1 - q_{tot}) * p_v^* / (p - p_v^*)``
and is set to 1 if `p - p_v_sat` is less than machine epsilon.
"""
@inline function q_vap_saturation_from_pressure_calc(
    param_set::APS,
    q_tot,
    p,
    p_v_sat,
)
    FT = eltype(param_set)
    R_v = TP.R_v(param_set)
    R_d = TP.R_d(param_set)
    q_v_sat = ifelse(
        p - p_v_sat ≥ ϵ_numerics(FT),
        R_d / R_v * (1 - q_tot) * p_v_sat / (p - p_v_sat),
        FT(1),
    )
    return q_v_sat
end

"""
    supersaturation(param_set, q_vap, ρ, T[, phase::Phase = Liquid()])

The supersaturation over water or ice.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `q_vap`: vapor specific humidity [kg/kg]
 - `ρ`: air density [kg/m³]
 - `T`: air temperature [K]
 - `phase`: (optional) liquid or ice phase to dispatch over (default: `Liquid()`)

# Returns
 - `S`: supersaturation [dimensionless], defined as `p_v/p_v^* - 1`
"""
@inline function supersaturation(
    param_set::APS,
    q_vap,
    ρ,
    T,
    phase::Phase = Liquid(),
)
    p_v_sat = saturation_vapor_pressure(param_set, T, phase)

    return supersaturation(param_set, q_vap, ρ, T, p_v_sat)
end

"""
    supersaturation(param_set, q_vap, ρ, T, p_v_sat)

The supersaturation given the saturation vapor pressure.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `q_vap`: vapor specific humidity [kg/kg]
 - `ρ`: air density [kg/m³]
 - `T`: temperature [K]
 - `p_v_sat`: saturation vapor pressure [Pa]

# Returns
 - `S`: supersaturation [dimensionless], defined as `p_v/p_v_sat - 1`
"""
@inline function supersaturation(param_set::APS, q_vap, ρ, T, p_v_sat)
    p_v = q_vap * (ρ * TP.R_v(param_set) * T)

    return p_v / p_v_sat - 1
end

"""
    saturation_excess(param_set, T, ρ, q_tot)
    saturation_excess(param_set, T, ρ, q_tot, q_liq, q_ice)

The saturation excess in equilibrium.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `ρ`: (moist-)air density [kg/m³]
 - `q_tot`: total specific humidity [kg/kg]
 - (optional) `q_liq`: liquid specific humidity [kg/kg]
 - (optional) `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `q_ex`: saturation excess [kg/kg]

The saturation excess is the difference between the total specific humidity `q_tot`
and the saturation specific humidity in equilibrium, and it is defined to be
nonzero only if this difference is positive: `q_ex = max(0, q_tot - q_v^*)`.
"""
@inline function saturation_excess(param_set::APS, T, ρ, q_tot, q_liq, q_ice)
    p_vap_sat = saturation_vapor_pressure(param_set, T, q_liq, q_ice)
    return saturation_excess(param_set, T, ρ, q_tot, p_vap_sat)
end

@inline function saturation_excess(param_set::APS, T, ρ, q_tot)
    p_vap_sat = saturation_vapor_pressure(param_set, T)
    return saturation_excess(param_set, T, ρ, q_tot, p_vap_sat)
end

"""
    saturation_excess(param_set, T, ρ, q_tot, p_vap_sat)

The saturation excess given the saturation vapor pressure.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `ρ`: air density [kg/m³]
 - `q_tot`: total specific humidity [kg/kg]
 - `p_vap_sat`: saturation vapor pressure [Pa]

# Returns
 - `q_ex`: saturation excess [kg/kg]

The saturation excess is the difference between the total specific humidity `q_tot`
and the saturation specific humidity, and it is defined to be nonzero only if
this difference is positive: `q_ex = max(0, q_tot - q_v^*)`.
"""
@inline function saturation_excess(param_set::APS, T, ρ, q_tot, p_vap_sat)
    q_vap_sat = q_vap_from_p_vap(param_set, T, ρ, p_vap_sat)
    return max(0, q_tot - q_vap_sat)
end

"""
    condensate_partition(param_set, T, ρ, q_tot)

Compute the equilibrium liquid and ice specific humidities from the condensate
in thermodynamic equilibrium, returning a tuple `(q_liq, q_ice)`.

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T`: temperature [K]
 - `ρ`: (moist-)air density [kg/m³]
 - `q_tot`: total specific humidity [kg/kg]

# Returns
 - `(q_liq, q_ice)`: tuple of liquid and ice specific humidities [kg/kg]

The condensate is partitioned into liquid and ice using the temperature-dependent
liquid fraction (see [`liquid_fraction`](@ref)) and the saturation excess (see 
[`saturation_excess`](@ref)).
"""
@inline function condensate_partition(param_set::APS, T, ρ, q_tot)
    λ = liquid_fraction(param_set, T)
    q_c = saturation_excess(param_set, T, ρ, q_tot)
    q_liq = λ * q_c
    q_ice = (1 - λ) * q_c
    return (q_liq, q_ice)
end

"""
    vapor_pressure_deficit(param_set, T, p, q_tot, q_liq, q_ice, phase::Phase)

The vapor pressure deficit (saturation vapor pressure minus actual vapor pressure, 
truncated to be non-negative) over a specific phase (Liquid or Ice).

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: air temperature [K]
 - `p`: air pressure [Pa]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]
 - `phase`: phase to compute saturation over (`Liquid()` or `Ice()`)

# Returns
 - `VPD`: vapor pressure deficit [Pa]
"""
@inline function vapor_pressure_deficit(
    param_set::APS,
    T,
    p,
    q_tot,
    q_liq,
    q_ice,
    phase::Phase,
)
    es = saturation_vapor_pressure(param_set, T, phase)
    ea = partial_pressure_vapor(param_set, p, q_tot, q_liq, q_ice)
    return ReLU(es - ea)
end

"""
    vapor_pressure_deficit(param_set, T, p, q_tot=0, q_liq=0, q_ice=0)

The vapor pressure deficit (saturation vapor pressure minus actual vapor pressure, 
truncated to be non-negative) over liquid water for temperatures above freezing 
and over ice for temperatures below freezing.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: air temperature [K]
 - `p`: air pressure [Pa]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `VPD`: vapor pressure deficit [Pa]

If the specific humidities are not given, the result is the saturation vapor pressure.
"""
@inline function vapor_pressure_deficit(
    param_set::APS,
    T,
    p,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    Tᶠ = TP.T_freeze(param_set)
    above_freezing = T > Tᶠ
    return ifelse(
        above_freezing,
        vapor_pressure_deficit(param_set, T, p, q_tot, q_liq, q_ice, Liquid()),
        vapor_pressure_deficit(param_set, T, p, q_tot, q_liq, q_ice, Ice()),
    )
end
