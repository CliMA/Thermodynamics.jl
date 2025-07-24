# Energies
export total_energy
export internal_energy
export internal_energy_sat
export internal_energy_dry
export internal_energy_vapor
export internal_energy_liquid
export internal_energy_ice
export moist_static_energy
export specific_enthalpy
export total_specific_enthalpy
export virtual_dry_static_energy


"""
    internal_energy(param_set, T[, q::PhasePartition])

The internal energy per unit mass, given 

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature

and, optionally,

 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function internal_energy(
    param_set::APS,
    T::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT <: Real}
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
    internal_energy(param_set, ts::ThermodynamicState)

The internal energy per unit mass, given a thermodynamic state `ts`.
"""
@inline internal_energy(param_set::APS, ts::ThermodynamicState) = ts.e_int

"""
    internal_energy_dry(param_set, T)

The dry air internal energy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline function internal_energy_dry(param_set::APS, T::FT) where {FT <: Real}
    T_0 = TP.T_0(param_set)
    cv_d = TP.cv_d(param_set)
    R_d = TP.R_d(param_set)
    return cv_d * (T - T_0) - R_d * T_0
end

"""
    internal_energy_dry(param_set, ts::ThermodynamicState)

The dry air internal energy, given a thermodynamic state `ts`.
"""
@inline internal_energy_dry(param_set::APS, ts::ThermodynamicState) =
    internal_energy_dry(param_set, air_temperature(param_set, ts))

"""
    internal_energy_vapor(param_set, T)

The water vapor internal energy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline function internal_energy_vapor(param_set::APS, T::FT) where {FT <: Real}
    T_0 = TP.T_0(param_set)
    cv_v = TP.cv_v(param_set)
    e_int_v0 = TP.e_int_v0(param_set)

    return cv_v * (T - T_0) + e_int_v0
end

"""
    internal_energy_vapor(param_set, ts::ThermodynamicState)

The water vapor internal energy, given a thermodynamic state `ts`.
"""
@inline internal_energy_vapor(param_set::APS, ts::ThermodynamicState) =
    internal_energy_vapor(param_set, air_temperature(param_set, ts))

"""
    internal_energy_liquid(param_set, T)

The liquid water internal energy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline function internal_energy_liquid(
    param_set::APS,
    T::FT,
) where {FT <: Real}
    T_0 = TP.T_0(param_set)
    cv_l = TP.cv_l(param_set)

    return cv_l * (T - T_0)
end

"""
    internal_energy_liquid(param_set, ts::ThermodynamicState)

The liquid water internal energy, given a thermodynamic state `ts`.
"""
@inline internal_energy_liquid(param_set::APS, ts::ThermodynamicState) =
    internal_energy_liquid(param_set, air_temperature(param_set, ts))

"""
    internal_energy_ice(param_set, T)

The ice internal energy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline function internal_energy_ice(param_set::APS, T::FT) where {FT <: Real}
    T_0 = TP.T_0(param_set)
    cv_i = TP.cv_i(param_set)
    e_int_i0 = TP.e_int_i0(param_set)

    return cv_i * (T - T_0) - e_int_i0
end

"""
    internal_energy_ice(param_set, ts::ThermodynamicState)

The ice internal energy, given a thermodynamic state `ts`.
"""
@inline internal_energy_ice(param_set::APS, ts::ThermodynamicState) =
    internal_energy_ice(param_set, air_temperature(param_set, ts))

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
    T::FT,
    ρ::FT,
    q_tot::FT,
    ::Type{phase_type},
    q_pt::PhasePartition{FT} = PhasePartition_equil(
        param_set,
        T,
        ρ,
        q_tot,
        phase_type,
    ),
) where {FT <: Real, phase_type <: ThermodynamicState}
    return internal_energy(param_set, T, q_pt)
end
@inline internal_energy_sat(param_set, T, ρ, q_tot, phase_type) =
    internal_energy_sat(param_set, promote(T, ρ, q_tot)..., phase_type)

"""
    internal_energy_sat(param_set, ts::ThermodynamicState)

The internal energy per unit mass in thermodynamic equilibrium 
at saturation with a fixed temperature and total specific humidity, 
given a thermodynamic state `ts`.
"""
@inline internal_energy_sat(param_set::APS, ts::ThermodynamicState) =
    internal_energy_sat(
        param_set,
        air_temperature(param_set, ts),
        air_density(param_set, ts),
        total_specific_humidity(param_set, ts),
        typeof(ts),
    )

"""
    total_energy(param_set, e_kin, e_pot, T[, q::PhasePartition])

The total energy per unit mass, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `e_kin` kinetic energy per unit mass
 - `e_pot` gravitational potential energy per unit mass
 - `T` temperature

and, optionally,

 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.

"""
@inline function total_energy(
    param_set::APS,
    e_kin::FT,
    e_pot::FT,
    T::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT <: Real}
    return internal_energy(param_set, T, q) + e_pot + e_kin
end

"""
    total_energy(param_set, ts::ThermodynamicState, e_kin, e_pot)

The total energy per unit mass, given a thermodynamic state `ts`.
"""
@inline function total_energy(
    param_set::APS,
    ts::ThermodynamicState{FT},
    e_kin::FT,
    e_pot::FT,
) where {FT <: Real}
    return internal_energy(param_set, ts) + e_pot + e_kin
end

"""
    total_energy_given_ρp(param_set, e_kin, e_pot, ρ, p[, q::PhasePartition])

The total energy per unit mass, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `p` pressure
 - `e_kin` kinetic energy per unit mass
 - `e_pot` gravitational potential energy per unit mass

and, optionally,

 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function total_energy_given_ρp(
    param_set::APS,
    ρ::FT,
    p::FT,
    e_kin::FT,
    e_pot::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT <: Real}
    T = air_temperature_from_ideal_gas_law(param_set, p, ρ, q)
    return total_energy(param_set, e_kin, e_pot, T, q)
end

"""
    specific_enthalpy(e_int, R_m, T)

The specific enthalpy, given

 - `e_int` internal specific energy
 - `R_m` [`gas_constant_air`](@ref)
 - `T` air temperature
"""
@inline function specific_enthalpy(e_int::FT, R_m::FT, T::FT) where {FT <: Real}
    return e_int + R_m * T
end

"""
    specific_enthalpy(param_set, ts)

The specific enthalpy, given a thermodynamic state `ts`.
"""
@inline function specific_enthalpy(
    param_set::APS,
    ts::ThermodynamicState{FT},
) where {FT <: Real}
    e_int = internal_energy(param_set, ts)
    R_m = gas_constant_air(param_set, ts)
    T = air_temperature(param_set, ts)
    return specific_enthalpy(e_int, R_m, T)
end

"""
    specific_enthalpy(param_set, T[, q::PhasePartition])

The specific enthalpy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature

and, optionally,

 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function specific_enthalpy(
    param_set::APS,
    T::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT <: Real}
    R_m = gas_constant_air(param_set, q)
    e_int = internal_energy(param_set, T, q)
    return specific_enthalpy(e_int, R_m, T)
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
    T::FT,
    ρ::FT,
    q_tot::FT,
    ::Type{phase_type},
) where {FT <: Real, phase_type <: ThermodynamicState}
    return specific_enthalpy(
        param_set,
        T,
        PhasePartition_equil(param_set, T, ρ, q_tot, phase_type),
    )
end


"""
    moist_static_energy(param_set, ts, e_pot)

The moist static energy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ts` a thermodynamic state
 - `e_pot` gravitational potential energy per unit mass
"""
@inline function moist_static_energy(
    param_set::APS,
    ts::ThermodynamicState,
    e_pot,
)
    return specific_enthalpy(param_set, ts) + e_pot
end

"""
    total_specific_enthalpy(e_tot, R_m, T)

The total specific enthalpy, given

 - `e_tot` total specific energy
 - `R_m` [`gas_constant_air`](@ref)
 - `T` air temperature
"""
@inline function total_specific_enthalpy(
    e_tot::FT,
    R_m::FT,
    T::FT,
) where {FT <: Real}
    return e_tot + R_m * T
end

"""
    total_specific_enthalpy(param_set, ts, e_tot::Real)

The total specific enthalpy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ts` a thermodynamic state
 - `e_tot` total specific energy
"""
@inline function total_specific_enthalpy(
    param_set::APS,
    ts::ThermodynamicState{FT},
    e_tot::FT,
) where {FT <: Real}
    R_m = gas_constant_air(param_set, ts)
    T = air_temperature(param_set, ts)
    return total_specific_enthalpy(e_tot, R_m, T)
end

"""
    virtual_dry_static_energy(param_set, ts, e_pot)

The virtual dry static energy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ts` a thermodynamic state
 - `e_pot` gravitational potential energy per unit mass

 Note that this static energy does not include the constant offset ``cp_d * T_0`` which is 
 present in the moist static energy.
"""
@inline function virtual_dry_static_energy(
    param_set::APS,
    ts::ThermodynamicState,
    e_pot,
)
    T_0 = TP.T_0(param_set)
    cp_d = TP.cp_d(param_set)
    T_virt = virtual_temperature(param_set, ts)
    return cp_d * T_virt + e_pot
end
