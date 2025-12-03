# Energies
export internal_energy
export internal_energy_dry
export internal_energy_vapor
export internal_energy_liquid
export internal_energy_ice
export internal_energy_sat
export total_energy
export moist_static_energy
export specific_enthalpy
export specific_enthalpy_dry
export specific_enthalpy_vapor
export specific_enthalpy_liquid
export specific_enthalpy_ice
export total_specific_enthalpy
export dry_static_energy
export virtual_dry_static_energy

"""
    internal_energy(param_set, T[, q::PhasePartition])

The internal energy per unit mass, given 

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air.
"""
@inline function internal_energy(
    param_set::APS,
    T,
    q::PhasePartition = q_pt_0(param_set),
)
    q_vap = vapor_specific_humidity(q)
    q_dry = 1 - q.tot

    return q_dry * internal_energy_dry(param_set, T) +
           q_vap * internal_energy_vapor(param_set, T) +
           q.liq * internal_energy_liquid(param_set, T) +
           q.ice * internal_energy_ice(param_set, T)
end
@inline internal_energy(param_set, T, q) =
    internal_energy(param_set, promote_phase_partition(T, q)...)

"""
    internal_energy_dry(param_set, T)

The dry air internal energy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
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

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
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

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
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

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline function internal_energy_ice(param_set::APS, T)
    T_0 = TP.T_0(param_set)
    cv_i = TP.cv_i(param_set)
    e_int_i0 = TP.e_int_i0(param_set)

    return cv_i * (T - T_0) - e_int_i0
end

"""
    internal_energy_sat(param_set, T, ρ, q_tot, phase_type)

The internal energy per unit mass in thermodynamic equilibrium 
at saturation with a fixed temperature and total specific humidity, 
given 

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
"""
@inline function internal_energy_sat(
    param_set::APS,
    T,
    ρ,
    q_tot,
    ::Type{phase_type},
    q_pt::PhasePartition = PhasePartition_equil(
        param_set,
        T,
        ρ,
        q_tot,
        phase_type,
    ),
) where {phase_type <: ThermodynamicState}
    return internal_energy(param_set, T, q_pt)
end
@inline internal_energy_sat(param_set, T, ρ, q_tot, phase_type) =
    internal_energy_sat(param_set, promote(T, ρ, q_tot)..., phase_type)

"""
    total_energy(param_set, e_kin, e_pot, T[, q::PhasePartition])

The total energy per unit mass, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `e_kin` kinetic energy per unit mass
 - `e_pot` gravitational potential energy per unit mass
 - `T` temperature

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air.
"""
@inline function total_energy(
    param_set::APS,
    e_kin,
    e_pot,
    T,
    q::PhasePartition = q_pt_0(param_set),
)
    return internal_energy(param_set, T, q) + e_pot + e_kin
end

# TODO Remove the method specific_enthalpy(e_int, R_m, T) in a future release (after ClimaAtmos update to not using this)
"""
    specific_enthalpy(e_int, R_m, T)

The specific enthalpy, given

 - `e_int` internal specific energy
 - `R_m` [`gas_constant_air`](@ref)
 - `T` air temperature

This method is deprecated and will be removed in a future release.
"""
@inline function specific_enthalpy(e_int, R_m, T)
    return e_int + R_m * T
end

"""
    specific_enthalpy(param_set, T[, q::PhasePartition])

The specific enthalpy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air.
"""
@inline function specific_enthalpy(
    param_set::APS,
    T,
    q::PhasePartition = q_pt_0(param_set),
)
    R_m = gas_constant_air(param_set, q)
    e_int = internal_energy(param_set, T, q)
    return e_int + R_m * T
end


"""
    specific_enthalpy_sat(param_set, T, ρ, q_tot, phase_type)

The specific enthalpy in thermodynamic equilibrium at saturation with a fixed temperature 
and total specific humidity, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
"""
@inline function specific_enthalpy_sat(
    param_set::APS,
    T,
    ρ,
    q_tot,
    ::Type{phase_type},
) where {phase_type <: ThermodynamicState}
    return specific_enthalpy(
        param_set,
        T,
        PhasePartition_equil(param_set, T, ρ, q_tot, phase_type),
    )
end

"""
    specific_enthalpy_dry(param_set, T)

The specific enthalpy of dry air, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline function specific_enthalpy_dry(param_set::APS, T)
    cp_d = TP.cp_d(param_set)
    T_0 = TP.T_0(param_set)
    return cp_d * (T - T_0)
end

"""
    specific_enthalpy_vapor(param_set, T)

The specific enthalpy of vapor, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline function specific_enthalpy_vapor(param_set::APS, T)
    cp_v = TP.cp_v(param_set)
    T_0 = TP.T_0(param_set)
    LH_v0 = TP.LH_v0(param_set)
    return cp_v * (T - T_0) + LH_v0
end

"""
    specific_enthalpy_liquid(param_set, T)

The specific enthalpy of liquid, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline specific_enthalpy_liquid(param_set::APS, T) =
    internal_energy_liquid(param_set, T)

"""
    specific_enthalpy_ice(param_set, T)

The specific enthalpy of ice, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline specific_enthalpy_ice(param_set::APS, T) =
    internal_energy_ice(param_set, T)

"""
    total_specific_enthalpy(param_set, e_tot, T[, q::PhasePartition])

The total specific enthalpy, defined as `e_tot + R_m * T`, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `e_tot` total specific energy
 - `T` temperature

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air.
"""
@inline function total_specific_enthalpy(
    param_set::APS,
    e_tot,
    T,
    q::PhasePartition = q_pt_0(param_set),
)
    R_m = gas_constant_air(param_set, q)
    return e_tot + R_m * T
end

# TODO: Remove the following method in a future release (after ClimaAtmos update to not using this)
"""
    total_specific_enthalpy(e_tot, R_m, T)

The total specific enthalpy, given

 - `e_tot` total specific energy
 - `R_m` [`gas_constant_air`](@ref)
 - `T` air temperature
"""
@inline function total_specific_enthalpy(e_tot, R_m, T)
    return e_tot + R_m * T
end

"""
    dry_static_energy(param_set, T, e_pot)

The dry static energy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `e_pot` gravitational potential energy per unit mass (geopotential)

"""
@inline function dry_static_energy(param_set::APS, T, e_pot)
    return specific_enthalpy_dry(param_set, T) + e_pot
end
