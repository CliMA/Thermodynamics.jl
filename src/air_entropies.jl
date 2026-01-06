# Entropies

export entropy
export entropy_dry
export entropy_vapor

"""
    entropy(param_set, p, T, q_tot=0, q_liq=0, q_ice=0)

The specific entropy in thermodynamic equilibrium.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `p`: pressure [Pa]
 - `T`: temperature [K]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `s`: specific entropy [J/(kg·K)]

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), this reduces to the dry-air expression.
The entropy is computed as a mass-weighted sum of the entropies of each component (dry air, vapor, liquid, ice).

# Reference
Pressel et al. (2015), "Numerics and subgrid-scale modeling in large eddy simulations of
stratocumulus clouds," *Journal of Advances in Modeling Earth Systems*, **7**(3), 1199-1220,
doi:[10.1002/2015MS000496](https://doi.org/10.1002/2015MS000496). (Their Eqs. (29)-(33))
"""
@inline function entropy(param_set::APS, p, T, q_tot = 0, q_liq = 0, q_ice = 0)
    L_v = latent_heat_vapor(param_set, T)
    L_s = latent_heat_sublim(param_set, T)
    s_d = entropy_dry(param_set, p, T, q_tot, q_liq, q_ice)
    s_v = entropy_vapor(param_set, p, T, q_tot, q_liq, q_ice)
    return (1 - q_tot) * s_d + q_tot * s_v - (q_liq * L_v + q_ice * L_s) / T
end

"""
    entropy_dry(param_set, p, T, q_tot=0, q_liq=0, q_ice=0)

The specific entropy of dry air at its partial pressure.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `p`: total air pressure [Pa]
 - `T`: temperature [K]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `s_d`: specific entropy of dry air [J/(kg·K)]

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), the dry-air partial pressure equals the total pressure.
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

The specific entropy of water vapor at its partial pressure.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `p`: total air pressure [Pa]
 - `T`: temperature [K]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `s_v`: specific entropy of water vapor [J/(kg·K)]

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
