export liquid_fraction
export saturation_vapor_pressure
export q_vap_saturation_generic
export q_vap_saturation
export q_vap_saturation_from_pressure
export saturation_excess
export supersaturation
export has_condensate
export vapor_pressure_deficit

"""
    liquid_fraction(param_set, T)
    liquid_fraction(param_set, T, q_liq, q_ice)

The fraction of condensate that is liquid.

If `q_liq` and `q_ice` are provided, the liquid fraction is computed from them.
If they are nonzero, `q_liq / (q_liq + q_ice)` is returned. If they are both
zero, a sharp temperature dependent partitioning is used (linear ramp from 0 to 
1 over +/- 0.1 K around freezing), so that the liquid fraction is zero below
freezing and one above freezing.

If `q_liq` and `q_ice` are not provided, the liquid fraction is computed from
temperature using a power law interpolation between `T_icenuc` and `T_freeze`
(Kaul et al., 2015).
"""
@inline function liquid_fraction(
    param_set::APS,
    T,
    q_liq,
    q_ice,
)
    FT = eltype(param_set)
    q_c = q_liq + q_ice

    # If no condensate, use sharp temperature dependent partitioning
    Tᶠ = TP.T_freeze(param_set)
    ΔT = FT(0.1) # Smooth over +/- 0.1 K
    # Linear ramp from 0 to 1 over [Tᶠ - ΔT, Tᶠ + ΔT]
    λ_T = clamp((T - (Tᶠ - ΔT)) / (2 * ΔT), zero(FT), one(FT))

    return ifelse(has_condensate(q_c), q_liq / q_c, λ_T)
end

@inline function liquid_fraction(
    param_set::APS,
    T,
)
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
        one(FT),
        ifelse(supercooled_liquid, λᵖ, zero(FT)),
    )
end

"""
    has_condensate(q_c)

Bool indicating if condensate exists, i.e., q_c > eps.
"""
@inline has_condensate(q_c) = q_c > eps(typeof(q_c))

"""
    saturation_vapor_pressure(param_set, T, ::Liquid)
    saturation_vapor_pressure(param_set, T, ::Ice)

The saturation vapor pressure over a plane surface of condensate, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `Liquid()` or `Ice()` to dispatch over the condensate type

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

Computes the saturation vapor pressure using the Rankine-Kirchhoff approximation, 
given:

 - `param_set` an `AbstractParameterSet`
 - `T` temperature
 - `LH_0` latent heat at reference temperature `T_0`
 - `Δcp` difference in isobaric specific heat capacity between the two phases

The computed value is:

``p_v^*(T) = p_{tr} \\left( \\frac{T}{T_{tr}} \\right)^{\\Delta c_{p} / R_v} \\exp \\left[ \\frac{L_0 - \\Delta c_{p} T_0}{R_v} \\left( \\frac{1}{T_{tr}} - \\frac{1}{T} \\right) \\right]``

where ``T_{tr}`` is the triple point temperature, ``p_{tr}`` is the triple point pressure, ``T_0`` is the reference temperature, and ``R_v`` is the gas constant for water vapor.
"""
@inline function saturation_vapor_pressure_calc(param_set::APS, T, LH_0, Δcp)
    press_triple = TP.press_triple(param_set)
    R_v = TP.R_v(param_set)
    T_triple = TP.T_triple(param_set)
    T_0 = TP.T_0(param_set)

    return press_triple *
           # (T / T_triple)^(Δcp / R_v) *
           fast_power(T / T_triple, Δcp / R_v) *
           exp((LH_0 - Δcp * T_0) / R_v * (1 / T_triple - 1 / T))
end

# Promote the arguments to a common type to allow AD with dual numbers
saturation_vapor_pressure_calc(param_set, T, LH_0, Δcp) =
    saturation_vapor_pressure_calc(param_set, promote(T, LH_0, Δcp)...)

"""
    saturation_vapor_pressure(param_set, T)
    saturation_vapor_pressure(param_set, T, q_liq, q_ice)

The saturation vapor pressure over liquid, ice, or a mixture of liquid and ice, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If `q_liq` and `q_ice` are provided, the saturation vapor pressure is computed
from a weighted mean of the latent heats of vaporization and sublimation, with
the weights given by the liquid fraction `liquid_fraction(param_set, T, q_liq, q_ice)`.
If `q_liq` and `q_ice` are 0, the saturation vapor pressure is that over liquid
above freezing and over ice below freezing.

Otherwise, the liquid fraction below freezing is computed from a temperature dependent
parameterization `liquid_fraction(param_set, T)`.
"""
@inline function saturation_vapor_pressure(
    param_set::APS,
    T,
    q_liq,
    q_ice,
)
    λ = liquid_fraction(param_set, T, q_liq, q_ice)
    return saturation_vapor_pressure_mixture(param_set, T, λ)
end

@inline function saturation_vapor_pressure(
    param_set::APS,
    T,
)
    λ = liquid_fraction(param_set, T)
    return saturation_vapor_pressure_mixture(param_set, T, λ)
end

"""
    saturation_vapor_pressure_mixture(param_set, T, λ)

Compute the saturation vapor pressure over a mixture of liquid and ice,
given the temperature `T` and liquid fraction `λ`.

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
    #latent heat
    LH_0 = λ * LH_v0 + (1 - λ) * LH_s0
    Δcp = λ * (cp_v - cp_l) + (1 - λ) * (cp_v - cp_i)

    # Saturation vapor pressure over possible mixture of liquid and ice
    return saturation_vapor_pressure_calc(param_set, T, LH_0, Δcp)
end

"""
    q_vap_saturation(param_set, T, ρ)
    q_vap_saturation(param_set, T, ρ, q_liq, q_ice)

The saturation specific humidity, given

- `param_set`: an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
- `T`: temperature
- `ρ`: air density
- (optional) `q_liq`: liquid specific humidity
- (optional) `q_ice`: ice specific humidity

If `q_liq` and `q_ice` are provided, the saturation specific humidity is that over a
mixture of liquid and ice, computed in a thermodynamically consistent way from the
weighted sum of the latent heats of the respective phase transitions (Pressel et al.,
JAMES, 2015). That is, the saturation vapor pressure and from it the saturation specific
humidity are computed from a weighted mean of the latent heats of vaporization and
sublimation, with the weights given by the liquid fraction. If `q_liq` and `q_ice` are 0, 
the saturation specific humidity is that over liquid above freezing and over ice below 
freezing.
 
Otherwise, the fraction of liquid is given by the temperature dependent `liquid_fraction(param_set, T)`.
"""
@inline function q_vap_saturation(
    param_set::APS,
    T,
    ρ,
    q_liq,
    q_ice,
)
    p_v_sat = saturation_vapor_pressure(param_set, T, q_liq, q_ice)
    return q_vap_from_p_vap(param_set, T, ρ, p_v_sat)
end

@inline function q_vap_saturation(
    param_set::APS,
    T,
    ρ,
)
    p_v_sat = saturation_vapor_pressure(param_set, T)
    return q_vap_from_p_vap(param_set, T, ρ, p_v_sat)
end

"""
    q_vap_saturation(param_set, T, ρ, phase::Phase)

The saturation specific humidity, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `phase` the phase to compute saturation over (either `Liquid()` or `Ice()`)
"""
@inline function q_vap_saturation(
    param_set::APS,
    T,
    ρ,
    phase::Phase,
)
    p_v_sat = saturation_vapor_pressure(param_set, T, phase)
    return q_vap_from_p_vap(param_set, T, ρ, p_v_sat)
end

"""
    q_vap_saturation_generic(param_set, T, ρ[, phase=Liquid()])

The saturation specific humidity over a plane surface of condensate, given
    - `param_set`: an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
    - `T`: temperature
    - `ρ`: air density
    - (optional) `Liquid()`: indicating condensate is liquid (default)
    - (optional) `Ice()`: indicating condensate is ice

The saturation specific humidity is computed as `q_v^* = p_v^*(T) / (ρ * R_v * T)`,
where `p_v^*` is the saturation vapor pressure.
"""
@inline function q_vap_saturation_generic(
    param_set::APS,
    T,
    ρ,
    phase::Phase = Liquid(),
)
    p_v_sat = saturation_vapor_pressure(param_set, T, phase)
    return q_vap_from_p_vap(param_set, T, ρ, p_v_sat)
end

"""
    q_vap_saturation_from_pressure(param_set, q_tot, p, T)
    q_vap_saturation_from_pressure(param_set, q_tot, p, T, q_liq, q_ice)

The saturation specific humidity, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q_tot` total water specific humidity,
 - `p` air pressure,
 - `T` air temperature
 - (optional) `q_liq` liquid specific humidity
 - (optional) `q_ice` ice specific humidity

If `q_liq` and `q_ice` are provided, the saturation vapor pressure is computed
from a weighted mean of the latent heats of vaporization and sublimation, with
the weights given by the liquid fraction `liquid_fraction(param_set, T, q_liq, q_ice)`.
If `q_liq` and `q_ice` are 0, the saturation vapor pressure is that over liquid
above freezing and over ice below freezing.

Otherwise, the liquid fraction is computed from a temperature dependent parameterization
`liquid_fraction(param_set, T)`.

The saturation specific humidity is computed as:
``q_v^* = (R_d / R_v) * (1 - q_{tot}) * p_v^* / (p - p_v^*)``
where `p_v^*` is the saturation vapor pressure.
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

@inline function q_vap_saturation_from_pressure(
    param_set::APS,
    q_tot,
    p,
    T,
)
    p_v_sat = saturation_vapor_pressure(param_set, T)
    return q_vap_saturation_from_pressure_calc(param_set, q_tot, p, p_v_sat)
end

"""
    q_vap_saturation_from_pressure_calc(param_set, q_tot, p, p_v_sat)

Compute the saturation specific humidity from the saturation vapor pressure `p_v_sat`.

The saturation specific humidity is computed as:
``q_v^* = (R_d / R_v) * (1 - q_{tot}) * p_v^* / (p - p_v^*)``
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
        p - p_v_sat ≥ eps(FT),
        R_d / R_v * (1 - q_tot) * p_v_sat / (p - p_v_sat),
        FT(1),
    )
    return q_v_sat
end

"""
    supersaturation(param_set, q_vap, ρ, T[, phase::Phase = Liquid()])

The supersaturation (pv/pv_sat -1) over water or ice, given
 - `param_set` - abstract set with earth parameters
 - `q_vap` - vapor specific humidity
 - `ρ` - air density,
 - `T` - air temperature
 - `phase` - (optional) liquid or ice phase to dispatch over (default: `Liquid()`).
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

The supersaturation (pv/pv_sat - 1) given the saturation vapor pressure `p_v_sat`:

- `param_set`: Thermodynamic parameter set
- `q_vap`: Vapor specific humidity
- `ρ`: Air density
- `T`: Temperature
- `p_v_sat`: Saturation vapor pressure

"""
@inline function supersaturation(
    param_set::APS,
    q_vap,
    ρ,
    T,
    p_v_sat,
)

    p_v = q_vap * (ρ * TP.R_v(param_set) * T)

    return p_v / p_v_sat - 1
end

"""
    saturation_excess(param_set, T, ρ, q_tot)
    saturation_excess(param_set, T, ρ, q_tot, q_liq, q_ice)

The saturation excess in equilibrium, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `q_tot` total specific humidity
 - (optional) `q_liq` liquid specific humidity
 - (optional) `q_ice` ice specific humidity

The saturation excess is the difference between the total specific humidity `q_tot`
and the saturation specific humidity in equilibrium, and it is defined to be
nonzero only if this difference is positive.
"""
@inline function saturation_excess(
    param_set::APS,
    T,
    ρ,
    q_tot,
    q_liq,
    q_ice,
)
    p_vap_sat = saturation_vapor_pressure(param_set, T, q_liq, q_ice)
    return saturation_excess(param_set, T, ρ, q_tot, p_vap_sat)
end

@inline function saturation_excess(
    param_set::APS,
    T,
    ρ,
    q_tot,
)
    p_vap_sat = saturation_vapor_pressure(param_set, T)
    return saturation_excess(param_set, T, ρ, q_tot, p_vap_sat)
end

"""
    saturation_excess(param_set, T, ρ, q_tot, p_vap_sat)

The saturation excess given the saturation vapor pressure `p_vap_sat`:

- `param_set`: Thermodynamic parameter set
- `T`: Temperature
- `ρ`: Air density
- `q_tot`: Total specific humidity
- `p_vap_sat`: Saturation vapor pressure

The saturation excess is the difference between the total specific humidity `q.tot`
and the saturation specific humidity, and it is defined to be nonzero only if
this difference is positive.
"""
@inline function saturation_excess(
    param_set::APS,
    T,
    ρ,
    q_tot,
    p_vap_sat,
)
    q_vap_sat = q_vap_from_p_vap(param_set, T, ρ, p_vap_sat)
    return max(0, q_tot - q_vap_sat)
end

"""
    vapor_pressure_deficit(param_set, T, p, q_tot=0, q_liq=0, q_ice=0)

The vapor pressure deficit (saturation vapor pressure minus actual 
vapor pressure, truncated to be non-negative) over liquid water for temperatures 
above freezing and over ice for temperatures below freezing, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` air temperature
 - `p` air pressure
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

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
    es = ifelse(
        above_freezing,
        saturation_vapor_pressure(param_set, T, Liquid()),
        saturation_vapor_pressure(param_set, T, Ice()),
    )

    ea = partial_pressure_vapor(param_set, p, q_tot, q_liq, q_ice)
    return ReLU(es - ea)
end
