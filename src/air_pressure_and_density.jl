# Air pressure and density from the ideal gas law
export air_pressure
export air_density
export exner
export exner_given_pressure

"""
    air_pressure(param_set, T, ρ, q_tot=0, q_liq=0, q_ice=0)

The air pressure from the equation of state (ideal gas law), given

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T`: air temperature [K]
 - `ρ`: (moist-)air density [kg/m³]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `p`: air pressure [Pa]

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), this reduces to the dry-air expression.
The formula is `p = ρ R_m T` where `R_m` is the gas constant of moist air (see [`gas_constant_air`](@ref)).
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

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T`: air temperature [K]
 - `p`: pressure [Pa]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `ρ`: (moist-)air density [kg/m³]

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), this reduces to the dry-air expression.
The formula is `ρ = p / (R_m T)` where `R_m` is the gas constant of moist air (see [`gas_constant_air`](@ref)).
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

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T`: temperature [K]
 - `ρ`: (moist-)air density [kg/m³]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `Π`: Exner function [dimensionless]

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), this reduces to the dry-air expression.
The pressure is computed internally from the equation of state (see [`air_pressure`](@ref)).
"""
@inline function exner(param_set::APS, T, ρ, q_tot = 0, q_liq = 0, q_ice = 0)
    p = air_pressure(param_set, T, ρ, q_tot, q_liq, q_ice)
    return exner_given_pressure(param_set, p, q_tot, q_liq, q_ice)
end

"""
    exner_given_pressure(param_set, p, q_tot=0, q_liq=0, q_ice=0)

The Exner function `Π = (p/p₀)^(R_m/cp_m)`, where `p₀` is the reference pressure,
given

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `p`: pressure [Pa]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `Π`: Exner function [dimensionless]

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), this reduces to the dry-air expression.
The reference pressure `p₀` is `p_ref_theta` from the parameter set.
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
