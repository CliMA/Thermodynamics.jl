export internal_energy
export internal_energy_dry
export internal_energy_vapor
export internal_energy_liquid
export internal_energy_ice
export internal_energy_sat
export enthalpy
export enthalpy_dry
export enthalpy_vapor
export enthalpy_liquid
export enthalpy_ice
export enthalpy_sat
export total_energy
export total_enthalpy
export dry_static_energy
export vapor_static_energy
export moist_static_energy
export virtual_dry_static_energy


"""
    internal_energy(param_set, T, q_tot=0, q_liq=0, q_ice=0)

The internal energy per unit mass, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), this reduces to the dry-air expression.
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

The dry air internal energy, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline function internal_energy_dry(param_set::APS, T)
    T_0 = TP.T_0(param_set)
    cv_d = TP.cv_d(param_set)
    R_d = TP.R_d(param_set)
    return cv_d * (T - T_0) - R_d * T_0
end

"""
    internal_energy_vapor(param_set, T)

The water vapor internal energy, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline function internal_energy_vapor(param_set::APS, T)
    T_0 = TP.T_0(param_set)
    cv_v = TP.cv_v(param_set)
    e_int_v0 = TP.e_int_v0(param_set)

    return cv_v * (T - T_0) + e_int_v0
end

"""
    internal_energy_liquid(param_set, T)

The liquid water internal energy, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline function internal_energy_liquid(param_set::APS, T)
    T_0 = TP.T_0(param_set)
    cv_l = TP.cv_l(param_set)

    return cv_l * (T - T_0)
end

"""
    internal_energy_ice(param_set, T)

The ice internal energy, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline function internal_energy_ice(param_set::APS, T)
    T_0 = TP.T_0(param_set)
    cv_i = TP.cv_i(param_set)
    e_int_i0 = TP.e_int_i0(param_set)

    return cv_i * (T - T_0) - e_int_i0
end

"""
    internal_energy_sat(param_set, T, ρ, q_tot)

The internal energy per unit mass in thermodynamic equilibrium at saturation, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `q_tot` total specific humidity

The phase partition into liquid and ice is computed internally from `q_tot` using the 
temperature-dependent liquid fraction and saturation excess.
"""
@inline function internal_energy_sat(param_set::APS, T, ρ, q_tot)
    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)
    return internal_energy(param_set, T, q_tot, q_liq, q_ice)
end

"""
    enthalpy(e_int, R_m, T)

The specific enthalpy, given

 - `e_int` internal specific energy
 - `R_m` [`gas_constant_air`](@ref)
 - `T` air temperature
"""
@inline function enthalpy(e_int, R_m, T)
    return e_int + R_m * T
end

"""
    enthalpy(param_set, T, q_tot=0, q_liq=0, q_ice=0)

The specific enthalpy, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), this reduces to the dry-air expression.
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

The specific enthalpy of dry air, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline function enthalpy_dry(param_set::APS, T)
    cp_d = TP.cp_d(param_set)
    T_0 = TP.T_0(param_set)
    return cp_d * (T - T_0)
end

"""
    enthalpy_vapor(param_set, T)

The specific enthalpy of vapor, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
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

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature

 The specific enthalpy of liquid is equal to the internal energy of liquid because the
 specific volume of condensed water is neglected.
"""
@inline enthalpy_liquid(param_set::APS, T) =
    internal_energy_liquid(param_set, T)

"""
    enthalpy_ice(param_set, T)

The specific enthalpy of ice, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature

 The specific enthalpy of ice is equal to the internal energy of ice because the
 specific volume of condensed water is neglected.
"""
@inline enthalpy_ice(param_set::APS, T) = internal_energy_ice(param_set, T)

"""
    enthalpy_sat(param_set, T, ρ, q_tot)

The specific enthalpy in thermodynamic equilibrium at saturation, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `q_tot` total specific humidity

The phase partition into liquid and ice is computed internally from `q_tot` using the 
temperature-dependent liquid fraction and saturation excess.
"""
@inline function enthalpy_sat(param_set::APS, T, ρ, q_tot)
    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)
    return enthalpy(param_set, T, q_tot, q_liq, q_ice)
end

"""
    total_energy(param_set, e_kin, e_pot, T, q_tot=0, q_liq=0, q_ice=0)

The total energy per unit mass, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `e_kin` kinetic energy per unit mass
 - `e_pot` gravitational potential energy per unit mass
 - `T` temperature
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

In the dry limit (`q_tot = q_liq = q_ice = 0`, the default), this reduces to the dry-air expression.
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

The total specific enthalpy, given

 - `e_tot` total specific energy
 - `R_m` [`gas_constant_air`](@ref)
 - `T` air temperature
"""
@inline function total_enthalpy(e_tot, R_m, T)
    return e_tot + R_m * T
end

"""
    dry_static_energy(param_set, T, e_pot)

The dry static energy, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `e_pot` gravitational potential energy per unit mass (geopotential)
"""
@inline function dry_static_energy(param_set::APS, T, e_pot)
    return enthalpy_dry(param_set, T) + e_pot
end

"""
    vapor_static_energy(param_set, T, e_pot)

The static energy (sensible heat only) of water vapor, `cp_v * (T - T_0) + e_pot`, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `e_pot` gravitational potential energy per unit mass (geopotential)
"""
@inline function vapor_static_energy(param_set::APS, T, e_pot)
    cp_v = TP.cp_v(param_set)
    T_0 = TP.T_0(param_set)
    return cp_v * (T - T_0) + e_pot
end

"""
    moist_static_energy(param_set, T, e_pot, q_tot=0, q_liq=0, q_ice=0)

The moist static energy, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `e_pot` gravitational potential energy per unit mass (geopotential)
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is for dry air.
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

The virtual dry static energy, given

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `e_pot` gravitational potential energy per unit mass (geopotential)
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is for dry air.
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
