export air_pressure
export exner

export air_density
export specific_volume

"""
    air_pressure(param_set, T, ρ[, q::PhasePartition])

The air pressure from the equation of state (ideal gas law), given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` air temperature
 - `ρ` (moist-)air density

and, optionally,

 - `q` [`PhasePartition`](@ref). 

When `q` is not provided, the results are for dry air.
"""
@inline function air_pressure(
    param_set::APS,
    T,
    ρ,
    q::PhasePartition = q_pt_0(param_set),
)
    return gas_constant_air(param_set, q) * ρ * T
end
@inline air_pressure(param_set, T, ρ, q) =
    air_pressure(param_set, promote_phase_partition(T, ρ, q)...)

"""
    air_density(param_set, T, p[, q::PhasePartition])

The (moist-)air density from the equation of state (ideal gas law), given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` air temperature
 - `p` pressure

and, optionally,

 - `q` [`PhasePartition`](@ref). 

When `q` is not provided, the results are for dry air.
"""
@inline function air_density(
    param_set::APS,
    T,
    p,
    q::PhasePartition = q_pt_0(param_set),
)
    return p / (gas_constant_air(param_set, q) * T)
end
@inline air_density(param_set, T, p, q) =
    air_density(param_set, promote_phase_partition(T, p, q)...)

"""
    exner(param_set, T, ρ[, q::PhasePartition)])

The Exner function, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density

and, optionally,

 - `q` [`PhasePartition`](@ref). 

When `q` is not provided, the results are for dry air.
"""
@inline function exner(
    param_set::APS,
    T,
    ρ,
    q::PhasePartition = q_pt_0(param_set),
    cpm = cp_m(param_set, q),
)
    p = air_pressure(param_set, T, ρ, q)
    return exner_given_pressure(param_set, p, q, cpm)
end

"""
    exner_given_pressure(param_set, p[, q::PhasePartition, cpm])

The Exner function, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure

and, optionally,

 - `q` [`PhasePartition`](@ref). 

When `q` is not provided, the results are for dry air.
"""
@inline function exner_given_pressure(
    param_set::APS,
    p,
    q::PhasePartition = q_pt_0(param_set),
    cpm = cp_m(param_set, q),
)
    p0 = TP.p_ref_theta(param_set)
    # gas constant and isobaric specific heat of moist air
    _R_m = gas_constant_air(param_set, q)

    # return (p / p0)^(_R_m / cpm)
    return fast_power(p / p0, _R_m / cpm)
end
