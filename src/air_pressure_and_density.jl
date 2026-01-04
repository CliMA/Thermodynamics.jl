# Air pressure and density from the ideal gas law
export air_pressure
export air_density
export exner
export exner_given_pressure

"""
    air_pressure(param_set, T, ρ, q_tot=0, q_liq=0, q_ice=0)

The air pressure from the equation of state (ideal gas law), given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` air temperature
 - `ρ` (moist-)air density
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is for dry air.
"""
@inline function air_pressure(
    param_set::APS,
    T,
    ρ,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    return gas_constant_air(param_set, q_tot, q_liq, q_ice) * ρ * T
end

"""
    air_density(param_set, T, p, q_tot=0, q_liq=0, q_ice=0)

The (moist-)air density from the equation of state (ideal gas law), given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` air temperature
 - `p` pressure
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is for dry air.
"""
@inline function air_density(
    param_set::APS,
    T,
    p,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    return p / (gas_constant_air(param_set, q_tot, q_liq, q_ice) * T)
end

"""
    exner(param_set, T, ρ, q_tot=0, q_liq=0, q_ice=0)

The Exner function `Π = (p/p₀)^(R_m/cp_m)`, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is for dry air.
The pressure is computed internally from the equation of state.
"""
@inline function exner(param_set::APS, T, ρ, q_tot = 0, q_liq = 0, q_ice = 0)
    p = air_pressure(param_set, T, ρ, q_tot, q_liq, q_ice)
    return exner_given_pressure(param_set, p, q_tot, q_liq, q_ice)
end

"""
    exner_given_pressure(param_set, p, q_tot=0, q_liq=0, q_ice=0)

The Exner function `Π = (p/p₀)^(R_m/cp_m)`, where `p₀` is the reference pressure,
given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is for dry air.
"""
@inline function exner_given_pressure(
    param_set::APS,
    p,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    p0 = TP.p_ref_theta(param_set)
    R_m = gas_constant_air(param_set, q_tot, q_liq, q_ice)
    cpm = cp_m(param_set, q_tot, q_liq, q_ice)

    return fast_power(p / p0, R_m / cpm)
end
