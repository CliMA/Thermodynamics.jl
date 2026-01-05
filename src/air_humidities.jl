export vapor_specific_humidity
export condensate_specific_humidity
export vol_vapor_mixing_ratio
export partial_pressure_dry
export partial_pressure_vapor
export shum_to_mixing_ratio
export q_vap_from_p_vap
export q_vap_from_RH
export q_vap_from_RH_liquid
export relative_humidity

"""
    vapor_specific_humidity(q_tot=0, q_liq=0, q_ice=0)

The vapor specific humidity (clamped to be non-negative), given
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is zero.
"""
@inline function vapor_specific_humidity(q_tot = 0, q_liq = 0, q_ice = 0)
    return max(0, q_tot - q_liq - q_ice)
end

"""
    condensate_specific_humidity(q_liq=0, q_ice=0)

The condensate specific humidity, given
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is zero.
"""
@inline function condensate_specific_humidity(q_liq = 0, q_ice = 0)
    return q_liq + q_ice
end

"""
    vol_vapor_mixing_ratio(param_set, q_tot=0, q_liq=0, q_ice=0)

The molar (volume) mixing ratio of water vapor, given
 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is zero.
"""
@inline function vol_vapor_mixing_ratio(
    param_set::APS,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    Rv_over_Rd = TP.Rv_over_Rd(param_set)
    q_vap = vapor_specific_humidity(q_tot, q_liq, q_ice)
    return Rv_over_Rd * shum_to_mixing_ratio(q_vap, q_tot)
end

"""
    partial_pressure_dry(param_set, p, q_tot=0, q_liq=0, q_ice=0)

The partial pressure of dry air, given
 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the partial pressure is the total pressure.
"""
@inline function partial_pressure_dry(
    param_set::APS,
    p,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    Rv_over_Rd = TP.Rv_over_Rd(param_set)
    q_vap = vapor_specific_humidity(q_tot, q_liq, q_ice)
    return p * (1 - q_tot) / (1 - q_tot + q_vap * Rv_over_Rd)
end

"""
    partial_pressure_vapor(param_set, p, q_tot=0, q_liq=0, q_ice=0)

The partial pressure of water vapor, given
 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the partial pressure is zero.
"""
@inline function partial_pressure_vapor(
    param_set::APS,
    p,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    Rv_over_Rd = TP.Rv_over_Rd(param_set)
    q_vap = vapor_specific_humidity(q_tot, q_liq, q_ice)
    return p * q_vap * Rv_over_Rd / (1 - q_tot + q_vap * Rv_over_Rd)
end

"""
    shum_to_mixing_ratio(q, q_tot)

The mixing ratio, given
 - `q` specific humidity
 - `q_tot` total specific humidity
"""
@inline function shum_to_mixing_ratio(q, q_tot)
    return q / (1 - q_tot)
end

"""
    q_vap_from_p_vap(param_set, T, ρ, p_v)

The vapor specific humidity, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature,
 - `ρ` (moist-)air density
 - `p_v` partial pressure of vapor
"""
@inline function q_vap_from_p_vap(param_set::APS, T, ρ, p_v)
    R_v = TP.R_v(param_set)
    return p_v / (ρ * R_v * T)
end

"""
    q_vap_from_RH(param_set, p, T, RH, phase)

The water vapor specific humidity, given
 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `RH` relative humidity with respect to `phase`
 - `phase` the phase to compute saturation over (either `Liquid()` or `Ice()`)
"""
@inline function q_vap_from_RH(param_set::APS, p, T, RH, phase::Phase)
    p_vap_sat = saturation_vapor_pressure(param_set, T, phase)
    p_vap = RH * p_vap_sat
    Rv_over_Rd = TP.Rv_over_Rd(param_set)
    return p_vap / Rv_over_Rd / (p - (1 - 1 / Rv_over_Rd) * p_vap)
end

"""
    q_vap_from_RH_liquid(param_set, p, T, RH)

The water vapor specific humidity, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `RH` relative humidity with respect to liquid water

 This function is deprecated. Use `q_vap_from_RH` with `Liquid()` instead.
"""
@inline function q_vap_from_RH_liquid(param_set::APS, p, T, RH)
    return q_vap_from_RH(param_set, p, T, RH, Liquid())
end

"""
    relative_humidity(param_set, T, p, q_tot=0, q_liq=0, q_ice=0)

The relative humidity (clipped between 0 and 1), given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `p` pressure
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

 If `q_liq` and `q_ice` are zero (or not given), the relative humidity is computed relative
 to saturation over ice below freezing and over liquid above freezing. If condensate is
 present, the relative humidity is computed relative to saturation over a mixture of liquid
 and ice, with the liquid fraction given by the ratio `q_liq / (q_liq + q_ice)`.

Note: `relative_humidity` uses `saturation_vapor_pressure(param_set, T, q_liq, q_ice)`. In
particular, for `q_liq == q_ice == 0` it includes a small smooth transition around freezing
(via `liquid_fraction(param_set, T, q_liq, q_ice)`).
"""
@inline function relative_humidity(
    param_set::APS,
    T,
    p,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    FT = eltype(param_set)
    p_vap = partial_pressure_vapor(param_set, p, q_tot, q_liq, q_ice)
    p_vap_sat = saturation_vapor_pressure(param_set, T, q_liq, q_ice)
    return max(0, min(1, p_vap / (p_vap_sat + eps(FT(0)))))
end
