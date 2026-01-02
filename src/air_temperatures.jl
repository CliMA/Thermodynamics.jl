# Temperature functions
export virtual_temperature
export air_temperature
export air_temperature_from_enthalpy
export air_temperature_given_ρp
export dry_pottemp
export dry_pottemp_given_pressure
export liquid_ice_pottemp
export liquid_ice_pottemp_given_pressure
export air_temperature_given_pθq
export air_temperature_given_ρθq
export virtual_pottemp

"""
    air_temperature(param_set, e_int, q_tot=0, q_liq=0, q_ice=0)

The air temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `e_int` specific internal energy of moist air
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is for dry air.
"""
@inline function air_temperature(
    param_set::APS,
    e_int,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    T_0 = TP.T_0(param_set)
    R_d = TP.R_d(param_set)
    e_int_v0 = TP.e_int_v0(param_set)
    e_int_i0 = TP.e_int_i0(param_set)
    cvm = cv_m(param_set, q_tot, q_liq, q_ice)
    q_vap = vapor_specific_humidity(q_tot, q_liq, q_ice)
    return T_0 +
           (
        e_int - q_vap * e_int_v0 +
        q_ice * e_int_i0 +
        (1 - q_tot) * R_d * T_0
    ) / cvm
end

"""
    air_temperature_from_enthalpy(param_set, h, q_tot=0, q_liq=0, q_ice=0)

The air temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `h` specific enthalpy of moist air
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is for dry air.
"""
@inline function air_temperature_from_enthalpy(
    param_set::APS,
    h,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    cpm = cp_m(param_set, q_tot, q_liq, q_ice)
    T_0 = TP.T_0(param_set)
    LH_v0 = TP.LH_v0(param_set)
    LH_f0 = TP.LH_f0(param_set)
    q_vap = vapor_specific_humidity(q_tot, q_liq, q_ice)
    return T_0 + (h - q_vap * LH_v0 + q_ice * LH_f0) / cpm
end

"""
    air_temperature_given_ρp(param_set, p, ρ, q_tot=0, q_liq=0, q_ice=0)

The air temperature, where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `ρ` air density
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is for dry air.
"""
@inline function air_temperature_given_ρp(
    param_set::APS,
    p,
    ρ,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    R_m = gas_constant_air(param_set, q_tot, q_liq, q_ice)
    return p / (R_m * ρ)
end

"""
    dry_pottemp(param_set, T, ρ, q_tot=0, q_liq=0, q_ice=0)

The dry potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is for dry air.
"""
@inline function dry_pottemp(
    param_set::APS,
    T,
    ρ,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    return T / exner(param_set, T, ρ, q_tot, q_liq, q_ice)
end

"""
    dry_pottemp_given_pressure(param_set, T, p, q_tot=0, q_liq=0, q_ice=0)

The dry potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `p` pressure
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is for dry air.
"""
@inline function dry_pottemp_given_pressure(
    param_set::APS,
    T,
    p,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    return T / exner_given_pressure(param_set, p, q_tot, q_liq, q_ice)
end

"""
    virtual_temperature(param_set, T, q_tot=0, q_liq=0, q_ice=0)

The virtual temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is for dry air.
"""
@inline function virtual_temperature(
    param_set::APS,
    T,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    R_d = TP.R_d(param_set)
    R_m = gas_constant_air(param_set, q_tot, q_liq, q_ice)
    return T * R_m / R_d
end

"""
    virtual_pottemp(param_set, T, ρ, q_tot=0, q_liq=0, q_ice=0)

The virtual potential temperature, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is the dry-air potential temperature.
"""
@inline function virtual_pottemp(
    param_set::APS,
    T,
    ρ,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    R_d = TP.R_d(param_set)
    R_m = gas_constant_air(param_set, q_tot, q_liq, q_ice)
    return dry_pottemp(param_set, T, ρ, q_tot, q_liq, q_ice) * R_m / R_d
end

"""
    liquid_ice_pottemp(param_set, T, ρ, q_tot=0, q_liq=0, q_ice=0)

The liquid-ice potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is the dry-air potential temperature.
"""
@inline function liquid_ice_pottemp(
    param_set::APS,
    T,
    ρ,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    return liquid_ice_pottemp_given_pressure(
        param_set,
        T,
        air_pressure(param_set, T, ρ, q_tot, q_liq, q_ice),
        q_tot,
        q_liq,
        q_ice,
    )
end

"""
    liquid_ice_pottemp_given_pressure(param_set, T, p, q_tot=0, q_liq=0, q_ice=0)

The liquid-ice potential temperature, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `p` pressure
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is the dry-air potential temperature.
"""
@inline function liquid_ice_pottemp_given_pressure(
    param_set::APS,
    T,
    p,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    # liquid-ice potential temperature, approximating latent heats
    # of phase transitions as constants
    cpm = cp_m(param_set, q_tot, q_liq, q_ice)
    return dry_pottemp_given_pressure(param_set, T, p, q_tot, q_liq, q_ice) *
           (1 - humidity_weighted_latent_heat(param_set, q_liq, q_ice) / (cpm * T))
end

"""
    air_temperature_given_pθq(
        param_set,
        p,
        θ_liq_ice,
        q_tot=0,
        q_liq=0,
        q_ice=0
    )

The air temperature obtained by inverting the liquid-ice potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the `θ_liq_ice` is assumed to be the dry-air potential temperature.
"""
@inline function air_temperature_given_pθq(
    param_set::APS,
    p,
    θ_liq_ice,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    cpm = cp_m(param_set, q_tot, q_liq, q_ice)
    return θ_liq_ice *
           exner_given_pressure(param_set, p, q_tot, q_liq, q_ice) +
           humidity_weighted_latent_heat(param_set, q_liq, q_ice) / cpm
end

"""
    air_temperature_given_ρθq(param_set, ρ, θ_liq_ice, q_tot=0, q_liq=0, q_ice=0)

The air temperature obtained by inverting the liquid-ice potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the results are for dry air.
"""
@inline function air_temperature_given_ρθq(
    param_set::APS,
    ρ,
    θ_liq_ice,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)

    p0 = TP.p_ref_theta(param_set)
    cvm = cv_m(param_set, q_tot, q_liq, q_ice)
    cpm = cp_m(param_set, q_tot, q_liq, q_ice)
    R_m = gas_constant_air(param_set, q_tot, q_liq, q_ice)
    κ = 1 - cvm / cpm
    T_u = (ρ * R_m * θ_liq_ice / p0)^(R_m / cvm) * θ_liq_ice
    L_q = humidity_weighted_latent_heat(param_set, q_liq, q_ice)
    T_1 = L_q / cvm
    T_2 = -κ / (2 * T_u) * (L_q / cvm)^2
    return T_u + T_1 + T_2
end
