# Material properties of moist air 

# Specific heats and gas constant of moist air
export gas_constant_air
export cp_m, cv_m

# Specific latent heats
export latent_heat_vapor
export latent_heat_sublim
export latent_heat_fusion
export latent_heat_mixed
export humidity_weighted_latent_heat

# Speed of sound
export soundspeed_air

"""
    gas_constant_air(param_set, q_tot=0, q_liq=0, q_ice=0)

The specific gas constant of moist air `R_m`, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

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

The isobaric specific heat capacity of moist air `cp_m`, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `q_tot` total specific humidity of water
 - `q_liq` specific humidity of liquid
 - `q_ice` specific humidity of ice

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

The isochoric specific heat capacity of moist air `cv_m`, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

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
    latent_heat_vapor(param_set, T)

The specific latent heat of vaporization `L_v`, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature

 Because the specific heats are assumed constant, the latent heats are linear functions of
 temperature by Kirchhoff's law.
"""
@inline function latent_heat_vapor(param_set::APS, T)
    cp_l = TP.cp_l(param_set)
    cp_v = TP.cp_v(param_set)
    LH_v0 = TP.LH_v0(param_set)
    return latent_heat_generic(param_set, T, LH_v0, cp_v - cp_l)
end

"""
    latent_heat_sublim(param_set, T)

The specific latent heat of sublimation `L_s`, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature

 Because the specific heats are assumed constant, the latent heats are linear functions of
 temperature by Kirchhoff's law.
"""
@inline function latent_heat_sublim(param_set::APS, T)
    LH_s0 = TP.LH_s0(param_set)
    cp_v = TP.cp_v(param_set)
    cp_i = TP.cp_i(param_set)
    return latent_heat_generic(param_set, T, LH_s0, cp_v - cp_i)
end

"""
    latent_heat_fusion(param_set, T)

The specific latent heat of fusion `L_f`, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature

 Because the specific heats are assumed constant, the latent heats are linear functions of
 temperature by Kirchhoff's law.
"""
@inline function latent_heat_fusion(param_set::APS, T)
    LH_f0 = TP.LH_f0(param_set)
    cp_l = TP.cp_l(param_set)
    cp_i = TP.cp_i(param_set)
    return latent_heat_generic(param_set, T, LH_f0, cp_l - cp_i)
end

"""
    latent_heat_generic(param_set, T, LH_0, Δcp)

Internal function. The specific latent heat of a generic phase transition between
two phases, computed using Kirchhoff's law, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `LH_0` latent heat at the reference temperature `T_0`
 - `Δcp` difference in isobaric specific heat capacities between the two phases

 Because the specific heats are assumed constant, the latent heats are linear functions of
 temperature by Kirchhoff's law.
"""
@inline function latent_heat_generic(param_set::APS, T, LH_0, Δcp)
    T_0 = TP.T_0(param_set)
    return LH_0 + Δcp * (T - T_0)
end

latent_heat_generic(param_set, T, LH_0, Δcp) =
    latent_heat_generic(param_set, promote(T, LH_0, Δcp)...)

"""
    latent_heat_mixed(param_set, T, λ)

The specific latent heat of a mixed phase, weighted by the liquid fraction `λ`, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` air temperature
 - `λ` liquid fraction

 Because the specific heats are assumed constant, the latent heats are linear functions of
 temperature by Kirchhoff's law.
"""
@inline function latent_heat_mixed(param_set::APS, T, λ)
    L_v = latent_heat_vapor(param_set, T)
    L_s = latent_heat_sublim(param_set, T)
    return λ * L_v + (1 - λ) * L_s
end

"""
    humidity_weighted_latent_heat(param_set, q_liq=0, q_ice=0)

Specific-humidity weighted sum of latent heats of liquid and ice evaluated at reference
temperature `T_0`, given
 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If `q_liq` and `q_ice` are not provided, `humidity_weighted_latent_heat` is zero.
"""
@inline function humidity_weighted_latent_heat(
    param_set::APS,
    q_liq = 0,
    q_ice = 0,
)
    LH_v0 = TP.LH_v0(param_set)
    LH_s0 = TP.LH_s0(param_set)
    return LH_v0 * q_liq + LH_s0 * q_ice
end

"""
    soundspeed_air(param_set, T, q_tot=0, q_liq=0, q_ice=0)

The speed of sound in unstratified air, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), this reduces to the dry-air expression.
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
