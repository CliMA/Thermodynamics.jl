# Entropies

export entropy
export entropy_dry
export entropy_vapor

"""
    entropy(param_set, p, T, q_tot=0, q_liq=0, q_ice=0)

The specific entropy in thermodynamic equilibrium, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is for dry air.
The specific entropy is computed from equations (29)-(33) of [Pressel2015](@cite).
"""
@inline function entropy(
    param_set::APS,
    p,
    T,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    L_v = latent_heat_vapor(param_set, T)
    L_s = latent_heat_sublim(param_set, T)
    s_d = entropy_dry(param_set, p, T, q_tot, q_liq, q_ice)
    s_v = entropy_vapor(param_set, p, T, q_tot, q_liq, q_ice)
    return (1 - q_tot) * s_d + q_tot * s_v - (q_liq * L_v + q_ice * L_s) / T
end

"""
    entropy_dry(param_set, p, T, q_tot=0, q_liq=0, q_ice=0)

The dry air specific entropy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is for dry air.
"""
@inline function entropy_dry(
    param_set::APS,
    p,
    T,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    FT = eltype(param_set)
    T_ref = TP.entropy_reference_temperature(param_set)
    p_ref = TP.MSLP(param_set)
    s_d_ref = TP.entropy_dry_air(param_set)
    R_d = TP.R_d(param_set)
    cp_d = TP.cp_d(param_set)
    p_d = partial_pressure_dry(param_set, p, q_tot, q_liq, q_ice)
    return s_d_ref + cp_d * log(T / T_ref) - R_d * log((p_d + eps(FT)) / p_ref)
end

"""
    entropy_vapor(param_set, p, T, q_tot=0, q_liq=0, q_ice=0)

The specific entropy of water vapor, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

Note: the entropy of water vapor diverges logarithmically as `q_tot → 0` (since `p_v → 0`).
"""
@inline function entropy_vapor(
    param_set::APS,
    p,
    T,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    FT = eltype(param_set)
    T_ref = TP.entropy_reference_temperature(param_set)
    p_ref = TP.MSLP(param_set)
    s_v_ref = TP.entropy_water_vapor(param_set)
    R_v = TP.R_v(param_set)
    cp_v = TP.cp_v(param_set)
    p_v = partial_pressure_vapor(param_set, p, q_tot, q_liq, q_ice)
    return s_v_ref + cp_v * log(T / T_ref) - R_v * log((p_v + eps(FT)) / p_ref)
end
