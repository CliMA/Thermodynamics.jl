# Phase partition functions

export liquid_fraction, PhasePartition_equil
export has_condensate

# From air_humidities.jl
export liquid_specific_humidity
export ice_specific_humidity
export vapor_specific_humidity
export condensate_specific_humidity
export mixing_ratios
export vol_vapor_mixing_ratio
export partial_pressure_dry
export partial_pressure_vapor
export vapor_pressure_deficit

# From air_pressure_and_density.jl
export air_pressure, air_density, exner, exner_given_pressure

# From air_energies.jl
export internal_energy, total_energy, enthalpy, total_enthalpy
export moist_static_energy, virtual_dry_static_energy
export entropy, entropy_dry, entropy_vapor

# From air_saturation_functions.jl
export saturation_vapor_pressure

# From air_temperatures.jl
export air_temperature
export air_temperature_given_ρp
export air_temperature_from_enthalpy
export dry_pottemp
export dry_pottemp_given_pressure
export latent_heat_liq_ice
export liquid_ice_pottemp
export liquid_ice_pottemp_given_pressure
export liquid_ice_pottemp_sat
export air_temperature_given_pθq
export air_temperature_given_ρθq
export air_temperature_given_ρθq_nonlinear
export virtual_pottemp
export virtual_temperature




"""
    has_condensate(q::PhasePartition{FT})

Bool indicating if condensate exists in the phase partition
"""
@inline has_condensate(q::PhasePartition) =
    has_condensate(condensate_specific_humidity(q))

@inline function liquid_fraction(
    param_set::APS,
    T,
    ::Type{phase_type},
    q::PhasePartition = q_pt_0(param_set),
) where {phase_type <: ThermodynamicState}
    return phase_type <: PhaseNonEquil ?
           liquid_fraction(param_set, T, q.liq, q.ice) :
           liquid_fraction(param_set, T)
end


"""
    PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    PhasePartition_equil(param_set, T, ρ, q_tot, p_vap_sat, λ)

Partition the phases in equilibrium, returning a [`PhasePartition`](@ref) object, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
 - `p_vap_sat` saturation vapor pressure
 - `λ` liquid fraction

The residual `q.tot - q.liq - q.ice` is the vapor specific humidity.
"""
@inline function PhasePartition_equil(param_set::APS, T, ρ, q_tot, p_vap_sat, λ)
    q_c = saturation_excess(param_set, T, ρ, p_vap_sat, PhasePartition(q_tot)) # condensate specific humidity
    q_liq = λ * q_c                                                            # liquid specific humidity
    q_ice = (1 - λ) * q_c                                                      # ice specific humidity
    return PhasePartition(q_tot, q_liq, q_ice)
end


"""
    PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)

Partition the phases in equilibrium, returning a [`PhasePartition`](@ref) object, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` air density
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type

If the specific humidity `q_tot` exceeds the saturation specific humidity `q_vap_sat`
(computed using `ρ`, `T` and the saturation vapor pressure), the condensate
is partitioned into liquid and ice. The fraction of liquid is given by the
temperature dependent `liquid_fraction(param_set, T)`.
"""
@inline function PhasePartition_equil(
    param_set::APS,
    T,
    ρ,
    q_tot,
    ::Type{phase_type},
) where {phase_type <: ThermodynamicState}
    # q_vap_sat and λ are needed in PhasePartition_equil, so we pass
    # the pre-computed values into the function to minimize the
    # number of times they are computed.
    λ = liquid_fraction(param_set, T)
    FT = eltype(param_set)
    p_vap_sat = saturation_vapor_pressure(param_set, T)
    return PhasePartition_equil(param_set, T, ρ, q_tot, p_vap_sat, λ)
end

@inline PhasePartition_equil(param_set, T, ρ, q_tot, phase_type) =
    PhasePartition_equil(param_set, promote(T, ρ, q_tot)..., phase_type)

"""
    PhasePartition_equil_given_p(param_set, T, p, q_tot, phase_type)

Partition the phases in equilibrium, returning a [`PhasePartition`](@ref) object using the
[`liquid_fraction`](@ref) function, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `p` air pressure
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
The residual `q.tot - q.liq - q.ice` is the vapor specific humidity.
"""
@inline function PhasePartition_equil_given_p(
    param_set::APS,
    T,
    p,
    q_tot,
    ::Type{phase_type},
    λ = liquid_fraction(param_set, T, phase_type),
) where {phase_type <: ThermodynamicState}

    q_v_sat = q_vap_saturation_from_pressure(param_set, q_tot, p, T)
    q_c = q_tot - q_v_sat
    q_liq = λ * q_c
    q_ice = (1 - λ) * q_c
    return PhasePartition(q_tot, q_liq, q_ice)
end
@inline PhasePartition_equil_given_p(param_set, T, p, q_tot, phase_type) =
    PhasePartition_equil_given_p(param_set, promote(T, p, q_tot)..., phase_type)

# -------------------------------------------------------------------------
# From air_humidities.jl
# -------------------------------------------------------------------------

"""
    liquid_specific_humidity(q::PhasePartition)

The liquid specific humidity, given

 - `q` a `PhasePartition`
"""
@inline liquid_specific_humidity(q::PhasePartition) = q.liq

"""
    ice_specific_humidity(q::PhasePartition)

The ice specific humidity, given

 - `q` a `PhasePartition`
"""
@inline ice_specific_humidity(q::PhasePartition) = q.ice

"""
    vapor_specific_humidity(q::PhasePartition)

The vapor specific humidity, given a 

- `q` a `PhasePartition` 
"""
@inline vapor_specific_humidity(q::PhasePartition) =
    vapor_specific_humidity(q.tot, q.liq, q.ice)

"""
    condensate_specific_humidity(q::PhasePartition{FT})

The condensate specific humidity (liquid + ice) of the phase 
partition `q`.
"""
@inline condensate_specific_humidity(q::PhasePartition) =
    condensate_specific_humidity(q.liq, q.ice)

"""
    mixing_ratios(q::PhasePartition)

The mixing ratios, given a specific humidity phase partition, `q`, returned in a 
`PhasePartition` with the fields
 - `r.tot` total mixing ratio
 - `r.liq` liquid mixing ratio
 - `r.ice` ice mixing ratio
"""
@inline function mixing_ratios(q::PhasePartition)
    return PhasePartition(
        shum_to_mixing_ratio(q.tot, q.tot),
        shum_to_mixing_ratio(q.liq, q.tot),
        shum_to_mixing_ratio(q.ice, q.tot),
    )
end

"""
    vol_vapor_mixing_ratio(param_set, q::PhasePartition)

The volume mixing ratio of water vapor, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q` [`PhasePartition`](@ref)
"""
@inline function vol_vapor_mixing_ratio(param_set::APS, q::PhasePartition)
    return vol_vapor_mixing_ratio(param_set, q.tot, q.liq, q.ice)
end

"""
    partial_pressure_dry(param_set, p, q::PhasePartition)

The partial pressure of dry air, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `q` phase partition

"""
@inline function partial_pressure_dry(
    param_set::APS,
    p,
    q::PhasePartition,
)
    return partial_pressure_dry(param_set, p, q.tot, q.liq, q.ice)
end

"""
    partial_pressure_vapor(param_set, p, q::PhasePartition)

The partial pressure of water vapor, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `q` phase partition

"""
@inline function partial_pressure_vapor(
    param_set::APS,
    p,
    q::PhasePartition,
)
    return partial_pressure_vapor(param_set, p, q.tot, q.liq, q.ice)
end

"""
    vapor_pressure_deficit(param_set, T, p[, q::PhasePartition])

The vapor pressure deficit (saturation vapor pressure minus actual 
vapor pressure, truncated to be non-negative) over liquid water for temperatures 
above freezing and over ice for temperatures below freezing, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` air temperature
 - `p` air pressure
 - `q` [`PhasePartition`](@ref)

When `q` is not provided, the vapor pressure deficit is the saturation vapor pressure.
"""
@inline function vapor_pressure_deficit(
    param_set::APS,
    T,
    p,
    q::PhasePartition = q_pt_0(param_set),
    Tᶠ = TP.T_freeze(param_set),
)
    above_freezing = T > Tᶠ
    es = ifelse(
        above_freezing,
        saturation_vapor_pressure(param_set, T, Liquid()),
        saturation_vapor_pressure(param_set, T, Ice()),
    )

    ea = partial_pressure_vapor(param_set, p, q)
    return ReLU(es - ea)
end

"""
    relative_humidity(param_set, T, p, phase_type, q::PhasePartition)

The relative humidity (clipped between 0 and 1), given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `p` pressure
 - `phase_type` a thermodynamic state type
 - `q` [`PhasePartition`](@ref).

"""
@inline function relative_humidity(
    param_set::APS,
    T,
    p,
    ::Type{phase_type},
    q::PhasePartition,
) where {phase_type <: ThermodynamicState}
    return relative_humidity(param_set, T, p, phase_type, q.tot, q.liq, q.ice)
end

# -------------------------------------------------------------------------
# From air_properties.jl
# -------------------------------------------------------------------------

"""
    gas_constant_air(param_set[, q::PhasePartition])

The specific gas constant of moist air, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q` [`PhasePartition`](@ref).

When `q` is not provided, the results are for dry air.
"""
@inline function gas_constant_air(param_set::APS, q::PhasePartition)
    return gas_constant_air(param_set, q.tot, q.liq, q.ice)
end

"""
    cp_m(param_set[, q::PhasePartition])

The isobaric specific heat capacity of moist air given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q` [`PhasePartition`](@ref).

When `q` is not provided, the results are for dry air.
"""
@inline function cp_m(param_set::APS, q::PhasePartition)
    return cp_m(param_set, q.tot, q.liq, q.ice)
end

"""
    cv_m(param_set[, q::PhasePartition])

The isochoric specific heat capacity of moist air, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q` [`PhasePartition`](@ref). Without humidity argument, the results are for dry air.
"""
@inline function cv_m(param_set::APS, q::PhasePartition)
    return cv_m(param_set, q.tot, q.liq, q.ice)
end

"""
    soundspeed_air(param_set, T[, q::PhasePartition])

The speed of sound in unstratified air, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature

and, optionally,

 - `q` [`PhasePartition`](@ref). 

When `q` is not provided, the results are for dry air.
"""
@inline function soundspeed_air(
    param_set::APS,
    T,
    q::PhasePartition,
)
    return soundspeed_air(param_set, T, q.tot, q.liq, q.ice)
end

# -------------------------------------------------------------------------
# From air_pressure_and_density.jl
# -------------------------------------------------------------------------

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

# -------------------------------------------------------------------------
# From air_energies.jl
# -------------------------------------------------------------------------

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
    q::PhasePartition,
)
    return internal_energy(param_set, T, q.tot, q.liq, q.ice)
end
@inline internal_energy(param_set, T, q) =
    internal_energy(param_set, promote_phase_partition(T, q)...)


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
    q::PhasePartition,
)
    return total_energy(param_set, e_kin, e_pot, T, q.tot, q.liq, q.ice)
end

"""
    enthalpy(param_set, T, q::PhasePartition)

The specific enthalpy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `q` [`PhasePartition`](@ref). 

"""
@inline function enthalpy(
    param_set::APS,
    T,
    q::PhasePartition,
)
    return enthalpy(param_set, T, q.tot, q.liq, q.ice)
end


"""
    moist_static_energy(param_set, T, e_pot, q::PhasePartition)

The moist static energy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `e_pot` gravitational potential energy per unit mass (geopotential)
 - `q` [`PhasePartition`](@ref). 

"""
@inline function moist_static_energy(
    param_set::APS,
    T,
    e_pot,
    q::PhasePartition,
)
    return moist_static_energy(param_set, T, e_pot, q.tot, q.liq, q.ice)
end

"""
    virtual_dry_static_energy(param_set, T, e_pot, q::PhasePartition)

The virtual dry static energy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `e_pot` gravitational potential energy per unit mass (geopotential)
 - `q` [`PhasePartition`](@ref). 

"""
@inline function virtual_dry_static_energy(
    param_set::APS,
    T,
    e_pot,
    q::PhasePartition,
)
    return virtual_dry_static_energy(param_set, T, e_pot, q.tot, q.liq, q.ice)
end



"""
    enthalpy_sat(param_set, T, ρ, q_tot, phase_type)

The specific enthalpy in thermodynamic equilibrium at saturation with a fixed temperature 
and total specific humidity, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
"""
@inline function enthalpy_sat(
    param_set::APS,
    T,
    ρ,
    q_tot,
    ::Type{phase_type},
) where {phase_type <: ThermodynamicState}
    return enthalpy(
        param_set,
        T,
        PhasePartition_equil(param_set, T, ρ, q_tot, phase_type),
    )
end

"""
    total_enthalpy(param_set, e_tot, T[, q::PhasePartition])

The total specific enthalpy, defined as `e_tot + R_m * T`, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `e_tot` total specific energy
 - `T` temperature

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air.
"""
@inline function total_enthalpy(
    param_set::APS,
    e_tot,
    T,
    q::PhasePartition = q_pt_0(param_set),
)
    R_m = gas_constant_air(param_set, q)
    return e_tot + R_m * T
end

# -------------------------------------------------------------------------
# From air_entropies.jl
# -------------------------------------------------------------------------

"""
    entropy(param_set, p, T, q::PhasePartition)

The specific entropy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `q` [`PhasePartition`](@ref).

The specific entropy is computed from equations (29)-(33) of [Pressel2015](@cite).
"""
@inline function entropy(param_set::APS, p, T, q::PhasePartition)
    return entropy(param_set, p, T, q.tot, q.liq, q.ice)
end

"""
    entropy_dry(param_set, p, T, q::PhasePartition)

The dry air specific entropy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `q` [`PhasePartition`](@ref).
"""
@inline function entropy_dry(
    param_set::APS,
    p,
    T,
    q::PhasePartition,
)
    return entropy_dry(param_set, p, T, q.tot, q.liq, q.ice)
end

"""
    entropy_vapor(param_set, p, T, q::PhasePartition)

The specific entropy of water vapor, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `q` [`PhasePartition`](@ref).
"""
@inline function entropy_vapor(
    param_set::APS,
    p,
    T,
    q::PhasePartition,
)
    return entropy_vapor(param_set, p, T, q.tot, q.liq, q.ice)
end


# -------------------------------------------------------------------------
# From air_saturation_functions.jl
# -------------------------------------------------------------------------

"""
    saturation_vapor_pressure(param_set, T, q::PhasePartition)

The saturation vapor pressure, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `q` [`PhasePartition`](@ref).
"""
@inline function saturation_vapor_pressure(
    param_set::APS,
    T::Number,
    q::PhasePartition,
)
    return saturation_vapor_pressure(param_set, T, q.liq, q.ice)
end


"""

    q_vap_saturation(param_set, T, ρ, q::PhasePartition)

The saturation specific humidity, given

- `param_set`: an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
- `T`: temperature
- `ρ`: air density
- `q`: [`PhasePartition`](@ref)

If the `PhasePartition` `q` is given, the saturation specific humidity is that over a
mixture of liquid and ice, computed in a thermodynamically consistent way from the
weighted sum of the latent heats of the respective phase transitions (Pressel et al.,
JAMES, 2015).
"""
@inline function q_vap_saturation(
    param_set::APS,
    T,
    ρ,
    q::PhasePartition,
)
    return q_vap_saturation(param_set, T, ρ, q.liq, q.ice)
end

# -------------------------------------------------------------------------
# From air_temperatures.jl
# -------------------------------------------------------------------------

"""
    air_temperature(param_set, e_int[, q::PhasePartition])

The air temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `e_int` specific internal energy of moist air

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air.
"""
@inline function air_temperature(
    param_set::APS,
    e_int,
    q::PhasePartition = q_pt_0(param_set),
    cvm::Number = cv_m(param_set, q),
)
    T_0 = TP.T_0(param_set)
    R_d = TP.R_d(param_set)
    e_int_v0 = TP.e_int_v0(param_set)
    e_int_i0 = TP.e_int_i0(param_set)
    return T_0 +
           (
        e_int - (q.tot - q.liq - q.ice) * e_int_v0 +
        q.ice * e_int_i0 +
        (1 - q.tot) * R_d * T_0
    ) / cvm
end

"""
    air_temperature_from_enthalpy(param_set, h[, q::PhasePartition]s)

The air temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `h` specific enthalpy of moist air

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air.
"""
@inline function air_temperature_from_enthalpy(
    param_set::APS,
    h,
    q::PhasePartition = q_pt_0(param_set),
)
    cp_m_ = cp_m(param_set, q)
    T_0 = TP.T_0(param_set)
    LH_v0 = TP.LH_v0(param_set)
    LH_f0 = TP.LH_f0(param_set)
    return T_0 + (h - (q.tot - q.liq - q.ice) * LH_v0 + q.ice * LH_f0) / cp_m_
end

"""
    air_temperature_given_ρp(param_set, p, ρ[, q::PhasePartition])

The air temperature, where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `ρ` air density

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air.
"""
@inline function air_temperature_given_ρp(
    param_set::APS,
    p,
    ρ,
    q::PhasePartition = q_pt_0(param_set),
)
    R_m = gas_constant_air(param_set, q)
    return p / (R_m * ρ)
end

"""
    dry_pottemp(param_set, T, ρ[, q::PhasePartition])

The dry potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air.
"""
@inline function dry_pottemp(
    param_set::APS,
    T,
    ρ,
    q::PhasePartition = q_pt_0(param_set),
    cpm = cp_m(param_set, q),
)
    return T / exner(param_set, T, ρ, q, cpm)
end

"""
    dry_pottemp_given_pressure(param_set, T, p[, q::PhasePartition])

The dry potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `p` pressure

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air, i.e., using the adiabatic exponent 
 for dry air.
"""
@inline function dry_pottemp_given_pressure(
    param_set::APS,
    T,
    p,
    q::PhasePartition = q_pt_0(param_set),
    cpm = cp_m(param_set, q),
)
    return T / exner_given_pressure(param_set, p, q, cpm)
end

"""
    latent_heat_liq_ice(param_set[, q::PhasePartition])

Specific-humidity weighted sum of latent heats of liquid and ice evaluated at reference temperature 
`T_0`, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q` [`PhasePartition`](@ref). 

When `q` is not provided, `latent_heat_liq_ice` is zero.

 This is used in the evaluation of the liquid-ice potential temperature.
"""
@inline function latent_heat_liq_ice(
    param_set::APS,
    q::PhasePartition = q_pt_0(param_set),
)
    LH_v0 = TP.LH_v0(param_set)
    LH_s0 = TP.LH_s0(param_set)
    return LH_v0 * q.liq + LH_s0 * q.ice
end

"""
    liquid_ice_pottemp_given_pressure(param_set, T, p[, q::PhasePartition])

The liquid-ice potential temperature, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `p` pressure

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the result is the dry-air potential temperature.
"""
@inline function liquid_ice_pottemp_given_pressure(
    param_set::APS,
    T,
    p,
    q::PhasePartition = q_pt_0(param_set),
    cpm = cp_m(param_set, q),
)
    # liquid-ice potential temperature, approximating latent heats
    # of phase transitions as constants
    return dry_pottemp_given_pressure(param_set, T, p, q, cpm) *
           (1 - latent_heat_liq_ice(param_set, q) / (cpm * T))
end

"""
    liquid_ice_pottemp(param_set, T, ρ[, q::PhasePartition])

The liquid-ice potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the result is the dry-air potential temperature.
"""
@inline function liquid_ice_pottemp(
    param_set::APS,
    T,
    ρ,
    q::PhasePartition = q_pt_0(param_set),
    cpm = cp_m(param_set, q),
)
    return liquid_ice_pottemp_given_pressure(
        param_set,
        T,
        air_pressure(param_set, T, ρ, q),
        q,
        cpm,
    )
end

"""
    air_temperature_given_pθq(
        param_set,
        p,
        θ_liq_ice,
        [q::PhasePartition]
    )

The air temperature obtained by inverting the liquid-ice potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `θ_liq_ice` liquid-ice potential temperature
 
and, optionally,
 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the `θ_liq_ice` is assumed to be the dry-air potential temperature.
"""
@inline function air_temperature_given_pθq(
    param_set::APS,
    p,
    θ_liq_ice,
    q::PhasePartition = q_pt_0(param_set),
    cpm = cp_m(param_set, q),
)
    return θ_liq_ice * exner_given_pressure(param_set, p, q, cpm) +
           latent_heat_liq_ice(param_set, q) / cpm
end

"""
    air_temperature_given_ρθq(param_set, ρ, θ_liq_ice[, q::PhasePartition])

The air temperature obtained by inverting the liquid-ice potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `θ_liq_ice` liquid-ice potential temperature

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air.
"""
@inline function air_temperature_given_ρθq(
    param_set::APS,
    ρ,
    θ_liq_ice,
    q::PhasePartition = q_pt_0(param_set),
)

    p0 = TP.p_ref_theta(param_set)
    cvm = cv_m(param_set, q)
    cpm = cp_m(param_set, q)
    R_m = gas_constant_air(param_set, q)
    κ = 1 - cvm / cpm
    T_u = (ρ * R_m * θ_liq_ice / p0)^(R_m / cvm) * θ_liq_ice
    T_1 = latent_heat_liq_ice(param_set, q) / cvm
    T_2 = -κ / (2 * T_u) * (latent_heat_liq_ice(param_set, q) / cvm)^2
    return T_u + T_1 + T_2
end


"""
    liquid_ice_pottemp_sat(param_set, T, ρ, phase_type[, q::PhasePartition, cpm])

The saturated liquid ice potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `phase_type` a thermodynamic state type

and, optionally,

 - `q` [`PhasePartition`](@ref). 

When `q` is not provided, the air assumed to be dry.
"""
@inline function liquid_ice_pottemp_sat(
    param_set::APS,
    T,
    ρ,
    ::Type{phase_type},
    q::PhasePartition = q_pt_0(param_set),
    cpm = cp_m(param_set, q),
) where {phase_type <: ThermodynamicState}
    q_v_sat = q_vap_saturation(param_set, T, ρ, phase_type, q)
    return liquid_ice_pottemp(param_set, T, ρ, PhasePartition(q_v_sat), cpm)
end

"""
    liquid_ice_pottemp_sat(param_set, T, ρ, phase_type, q_tot)

The saturated liquid ice potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `phase_type` a thermodynamic state type
 - `q_tot` total specific humidity
"""
@inline function liquid_ice_pottemp_sat(
    param_set::APS,
    T,
    ρ,
    ::Type{phase_type},
    q_tot,
) where {phase_type <: ThermodynamicState}
    q = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    cpm = cp_m(param_set, q)
    return liquid_ice_pottemp(param_set, T, ρ, q, cpm)
end

"""
    virtual_pottemp(param_set, T, ρ[, q::PhasePartition])

The virtual potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
The virtual potential temperature is defined as `θ_v = (R_m/R_d) * θ`, where `θ` is the
potential temperature. It is the potential temperature a dry air parcel would need to have to
have the same density as the moist air parcel at the same pressure.
"""
@inline function virtual_pottemp(
    param_set::APS,
    T,
    ρ,
    q::PhasePartition,
    cpm = cp_m(param_set, q),
)
    return virtual_pottemp(param_set, T, ρ, q.tot, q.liq, q.ice, cpm)
end

"""
    virtual_temperature(param_set, T, q::PhasePartition)

The virtual temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `q` [`PhasePartition`](@ref). 
 
The virtual temperature is defined as `T_v = (R_m/R_d) * T`. It is the temperature
a dry air parcel would need to have to have the same density as the moist air parcel
at the same pressure.
"""
@inline function virtual_temperature(param_set::APS, T, q::PhasePartition)
    return virtual_temperature(param_set, T, q.tot, q.liq, q.ice)
end


"""
    air_temperature_given_ρθq_nonlinear(param_set, ρ, θ_liq_ice, maxiter, tol, q::PhasePartition)

Computes temperature `T`, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `θ_liq_ice` liquid-ice potential temperature

and, optionally,
 - `maxiter` maximum iterations for non-linear equation solve
- `tol` absolute tolerance for non-linear equation iterations. Can be one of:
    - `RelativeSolutionTolerance()` to stop when `|x_2 - x_1|/x_1 < tol`
    - `ResidualTolerance()` to stop when `|f(x)| < tol`
 - `q` [`PhasePartition`](@ref).When `q` is not provided, the results are for dry air,

The temperature `T` is found by finding the root of
`T - air_temperature_given_pθq(param_set,
                               air_pressure(param_set, T, ρ, q),
                               θ_liq_ice,
                               q) = 0`
"""
@inline function air_temperature_given_ρθq_nonlinear(
    param_set::APS,
    ρ,
    θ_liq_ice,
    maxiter::Int,
    tol::RS.AbstractTolerance,
    q::PhasePartition = q_pt_0(param_set),
)
    T_init_min = TP.T_init_min(param_set)
    _T_max = TP.T_max(param_set)
    @inline roots(T) =
        T - air_temperature_given_pθq(
            param_set,
            air_pressure(param_set, ReLU(T), ρ, q),
            θ_liq_ice,
            q,
        )
    sol = RS.find_zero(
        roots,
        RS.SecantMethod(T_init_min, _T_max),
        RS.CompactSolution(),
        tol,
        maxiter,
    )
    if !sol.converged
        if print_warning()
            print_warning_ρθq_nonlinear(
                RS.SecantMethod,
                θ_liq_ice,
                ρ,
                q.tot,
                q.liq,
                q.ice,
                sol.root,
                maxiter,
                tol.tol,
            )
        end
        if error_on_non_convergence()
            error("Failed to converge with printed set of inputs.")
        end
    end
    return sol.root
end
