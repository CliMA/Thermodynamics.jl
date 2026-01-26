# Material properties of moist air 

# Specific heats and gas constant of moist air
export gas_constant_air
export cp_m, cv_m

# Speed of sound
export soundspeed_air

"""
    gas_constant_air(param_set, q_tot=0, q_liq=0, q_ice=0)

The specific gas constant of moist air `R_m`.

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `R_m`: specific gas constant of moist air [J/(kg·K)]

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), this reduces to the dry-air expression.
"""
@inline function gas_constant_air(
    param_set::APS,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    R_d = TP.R_d(param_set)
    R_v = TP.R_v(param_set)
    q_vap = vapor_specific_humidity(q_tot, q_liq, q_ice)
    return R_d * (1 - q_tot) + R_v * q_vap
end

"""
    cp_m(param_set, q_tot=0, q_liq=0, q_ice=0)

The isobaric specific heat capacity of moist air `cp_m`.

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `q_tot`: total specific humidity of water [kg/kg]
 - `q_liq`: specific humidity of liquid [kg/kg]
 - `q_ice`: specific humidity of ice [kg/kg]

# Returns
 - `cp_m`: isobaric specific heat capacity [J/(kg·K)]

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), this reduces to the dry-air expression.
The specific heat capacities are assumed to be constant (calorically perfect air).
"""
@inline function cp_m(param_set::APS, q_tot = 0, q_liq = 0, q_ice = 0)
    cp_d = TP.cp_d(param_set)
    cp_v = TP.cp_v(param_set)
    cp_l = TP.cp_l(param_set)
    cp_i = TP.cp_i(param_set)
    # rearranged formula for cp_m to avoid computation of vapor specific humidity
    return cp_d +
           (cp_v - cp_d) * q_tot +
           (cp_l - cp_v) * q_liq +
           (cp_i - cp_v) * q_ice
end

"""
    cv_m(param_set, q_tot=0, q_liq=0, q_ice=0)

The isochoric specific heat capacity of moist air `cv_m`.

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `cv_m`: isochoric specific heat capacity [J/(kg·K)]

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), this reduces to the dry-air expression.
The specific heat capacities are assumed to be constant (calorically perfect air).
"""
@inline function cv_m(param_set::APS, q_tot = 0, q_liq = 0, q_ice = 0)
    cv_d = TP.cv_d(param_set)
    cv_v = TP.cv_v(param_set)
    cv_l = TP.cv_l(param_set)
    cv_i = TP.cv_i(param_set)
    return cv_d +
           (cv_v - cv_d) * q_tot +
           (cv_l - cv_v) * q_liq +
           (cv_i - cv_v) * q_ice
end



"""
    soundspeed_air(param_set, T, q_tot=0, q_liq=0, q_ice=0)

The speed of sound in unstratified air.

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T`: temperature [K]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `c`: speed of sound [m/s]

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), this reduces to the dry-air expression.
The formula is `c = √(γ R_m T)` where `γ = cp_m/cv_m`.
"""
@inline function soundspeed_air(
    param_set::APS,
    T,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    γ =
        cp_m(param_set, q_tot, q_liq, q_ice) /
        cv_m(param_set, q_tot, q_liq, q_ice)
    R_m = gas_constant_air(param_set, q_tot, q_liq, q_ice)
    return sqrt(γ * R_m * T)
end
