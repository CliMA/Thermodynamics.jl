export internal_energy
export internal_energy_dry
export internal_energy_vapor
export internal_energy_liquid
export internal_energy_ice
export enthalpy
export enthalpy_dry
export enthalpy_vapor
export enthalpy_liquid
export enthalpy_ice
export total_energy
export total_enthalpy
export dry_static_energy
export vapor_static_energy
export moist_static_energy
export virtual_dry_static_energy


"""
    internal_energy(param_set, T, q_tot=0, q_liq=0, q_ice=0)

The internal energy per unit mass.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `e_int`: specific internal energy [J/kg]

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), this reduces to the dry-air expression.
The internal energy is computed as a mass-weighted sum of the internal energies of each component
(dry air, vapor, liquid, ice), referenced to `T_0`.
"""
@inline function internal_energy(
    param_set::APS,
    T,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    q_vap = vapor_specific_humidity(q_tot, q_liq, q_ice)
    q_dry = 1 - q_tot

    return q_dry * internal_energy_dry(param_set, T) +
           q_vap * internal_energy_vapor(param_set, T) +
           q_liq * internal_energy_liquid(param_set, T) +
           q_ice * internal_energy_ice(param_set, T)
end

"""
    internal_energy_dry(param_set, T)

The dry air internal energy.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]

# Returns
 - `e_d`: specific internal energy of dry air [J/kg]
"""
@inline function internal_energy_dry(param_set::APS, T)
    T_0 = TP.T_0(param_set)
    cv_d = TP.cv_d(param_set)
    R_d = TP.R_d(param_set)
    return cv_d * (T - T_0) - R_d * T_0
end

"""
    internal_energy_vapor(param_set, T)

The water vapor internal energy.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]

# Returns
 - `e_v`: specific internal energy of water vapor [J/kg]
"""
@inline function internal_energy_vapor(param_set::APS, T)
    T_0 = TP.T_0(param_set)
    cv_v = TP.cv_v(param_set)
    e_int_v0 = TP.e_int_v0(param_set)

    return cv_v * (T - T_0) + e_int_v0
end

"""
    internal_energy_liquid(param_set, T)

The liquid water internal energy.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]

# Returns
 - `e_l`: specific internal energy of liquid water [J/kg]
"""
@inline function internal_energy_liquid(param_set::APS, T)
    T_0 = TP.T_0(param_set)
    cv_l = TP.cv_l(param_set)

    return cv_l * (T - T_0)
end

"""
    internal_energy_ice(param_set, T)

The ice internal energy.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]

# Returns
 - `e_i`: specific internal energy of ice [J/kg]
"""
@inline function internal_energy_ice(param_set::APS, T)
    T_0 = TP.T_0(param_set)
    cv_i = TP.cv_i(param_set)
    e_int_i0 = TP.e_int_i0(param_set)

    return cv_i * (T - T_0) - e_int_i0
end

"""
    enthalpy(e_int, R_m, T)

The specific enthalpy, given

# Arguments
 - `e_int`: internal specific energy [J/kg]
 - `R_m`: gas constant of moist air [J/(kg·K)], see [`gas_constant_air`](@ref)
 - `T`: air temperature [K]

# Returns
 - `h`: specific enthalpy [J/kg]

The enthalpy is computed as `h = e_int + R_m T`, which follows from `h = e_int + p v`
with the ideal gas law `p v = R_m T` (specific volume v = 1/ρ).
"""
@inline function enthalpy(e_int, R_m, T)
    return e_int + R_m * T
end

"""
    enthalpy(param_set, T, q_tot=0, q_liq=0, q_ice=0)

The specific enthalpy.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `h`: specific enthalpy [J/kg]

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), this reduces to the dry-air expression.
The enthalpy is computed as a mass-weighted sum of the enthalpies of each component (dry air, vapor, liquid, ice).
"""
@inline function enthalpy(param_set::APS, T, q_tot = 0, q_liq = 0, q_ice = 0)
    q_vap = vapor_specific_humidity(q_tot, q_liq, q_ice)
    q_dry = 1 - q_tot

    return q_dry * enthalpy_dry(param_set, T) +
           q_vap * enthalpy_vapor(param_set, T) +
           q_liq * enthalpy_liquid(param_set, T) +
           q_ice * enthalpy_ice(param_set, T)
end

"""
    enthalpy_dry(param_set, T)

The specific enthalpy of dry air.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]

# Returns
 - `h_d`: specific enthalpy of dry air [J/kg]
"""
@inline function enthalpy_dry(param_set::APS, T)
    cp_d = TP.cp_d(param_set)
    T_0 = TP.T_0(param_set)
    return cp_d * (T - T_0)
end

"""
    enthalpy_vapor(param_set, T)

The specific enthalpy of water vapor.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]

# Returns
 - `h_v`: specific enthalpy of water vapor [J/kg]
"""
@inline function enthalpy_vapor(param_set::APS, T)
    cp_v = TP.cp_v(param_set)
    T_0 = TP.T_0(param_set)
    LH_v0 = TP.LH_v0(param_set)
    return cp_v * (T - T_0) + LH_v0
end

"""
    enthalpy_liquid(param_set, T)

The specific enthalpy of liquid, given

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T`: temperature [K]

# Returns
 - `h_l`: specific enthalpy of liquid [J/kg]

The specific enthalpy of liquid is equal to the internal energy of liquid because the
specific volume of condensed water is neglected (i.e., `p v_l ≈ 0`).
"""
@inline enthalpy_liquid(param_set::APS, T) =
    internal_energy_liquid(param_set, T)

"""
    enthalpy_ice(param_set, T)

The specific enthalpy of ice, given

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T`: temperature [K]

# Returns
 - `h_i`: specific enthalpy of ice [J/kg]

The specific enthalpy of ice is equal to the internal energy of ice because the
specific volume of condensed water is neglected (i.e., `p v_i ≈ 0`).
"""
@inline enthalpy_ice(param_set::APS, T) = internal_energy_ice(param_set, T)

"""
    total_energy(param_set, e_kin, e_pot, T, q_tot=0, q_liq=0, q_ice=0)

The total energy per unit mass.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `e_kin`: kinetic energy per unit mass [J/kg]
 - `e_pot`: gravitational potential energy per unit mass [J/kg] (geopotential)
 - `T`: temperature [K]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `e_tot`: total specific energy [J/kg]

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), this reduces to the dry-air expression.
The total energy is `e_tot = e_int + e_kin + e_pot`.
"""
@inline function total_energy(
    param_set::APS,
    e_kin,
    e_pot,
    T,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    return internal_energy(param_set, T, q_tot, q_liq, q_ice) + e_pot + e_kin
end

"""
    total_enthalpy(e_tot, R_m, T)

The total specific enthalpy.

# Arguments
 - `e_tot`: total specific energy [J/kg]
 - `R_m`: gas constant of moist air [J/(kg·K)], see [`gas_constant_air`](@ref)
 - `T`: air temperature [K]

# Returns
 - `h_tot`: total specific enthalpy [J/kg]

The total enthalpy is computed as `h_tot = e_tot + R_m T`.
"""
@inline function total_enthalpy(e_tot, R_m, T)
    return e_tot + R_m * T
end

"""
    dry_static_energy(param_set, T, e_pot)

The dry static energy.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `e_pot`: gravitational potential energy per unit mass (geopotential) [J/kg]

# Returns
 - `s_d`: dry static energy [J/kg]

The dry static energy is the sum of the dry enthalpy and geopotential: `s_d = h_d + e_pot`.
"""
@inline function dry_static_energy(param_set::APS, T, e_pot)
    return enthalpy_dry(param_set, T) + e_pot
end

"""
    vapor_static_energy(param_set, T, e_pot)

The static energy (sensible heat only) of water vapor.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `e_pot`: gravitational potential energy per unit mass (geopotential) [J/kg]

# Returns
 - `s_v`: vapor static energy [J/kg]

The formula is `s_v = cp_v * (T - T_0) + e_pot`, where `T_0` is the reference temperature.
"""
@inline function vapor_static_energy(param_set::APS, T, e_pot)
    cp_v = TP.cp_v(param_set)
    T_0 = TP.T_0(param_set)
    return cp_v * (T - T_0) + e_pot
end

"""
    moist_static_energy(param_set, T, e_pot, q_tot=0, q_liq=0, q_ice=0)

The moist static energy.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `e_pot`: gravitational potential energy per unit mass (geopotential) [J/kg]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `s_m`: moist static energy [J/kg]

If the specific humidities are not given, the result is for dry air.
The moist static energy is the sum of the moist enthalpy and geopotential: `s_m = h + e_pot`,
where `h` is computed from [`enthalpy`](@ref).
"""
@inline function moist_static_energy(
    param_set::APS,
    T,
    e_pot,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    h = enthalpy(param_set, T, q_tot, q_liq, q_ice)
    return h + e_pot
end

"""
    virtual_dry_static_energy(param_set, T, e_pot, q_tot=0, q_liq=0, q_ice=0)

The virtual dry static energy.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `e_pot`: gravitational potential energy per unit mass (geopotential) [J/kg]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `s_vd`: virtual dry static energy [J/kg]

If the specific humidities are not given, the result is for dry air.
The virtual dry static energy is `s_vd = cp_d * T_virt + e_pot`, where `T_virt` is
the virtual temperature (see [`virtual_temperature`](@ref)).
"""
@inline function virtual_dry_static_energy(
    param_set::APS,
    T,
    e_pot,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    cp_d = TP.cp_d(param_set)
    T_virt = virtual_temperature(param_set, T, q_tot, q_liq, q_ice)
    return cp_d * T_virt + e_pot
end
