# Atmospheric equation of state
export air_pressure
export air_temperature
export air_density
export specific_volume
export soundspeed_air
export total_specific_humidity
export liquid_specific_humidity
export ice_specific_humidity
export vapor_specific_humidity

# Energies
export total_energy
export internal_energy
export internal_energy_sat
export internal_energy_dry
export internal_energy_vapor
export internal_energy_liquid
export internal_energy_ice

# Specific heats and gas constants of moist air
export cp_m, cv_m, gas_constant_air, gas_constants

# Latent heats
export latent_heat_vapor
export latent_heat_sublim
export latent_heat_fusion
export latent_heat_liq_ice

# Saturation vapor pressures and specific humidities over liquid and ice
export Liquid, Ice
export saturation_vapor_pressure
export q_vap_saturation_generic
export q_vap_saturation
export q_vap_saturation_liquid
export q_vap_saturation_ice
export saturation_excess
export supersaturation

# Functions used in thermodynamic equilibrium among phases (liquid and ice
# determined diagnostically from total water specific humidity)
export liquid_fraction, PhasePartition_equil

# Auxiliary functions, e.g., for diagnostic purposes
export dry_pottemp
export virtual_pottemp
export virtual_dry_static_energy
export exner
export shum_to_mixing_ratio
export mixing_ratios
export vol_vapor_mixing_ratio
export liquid_ice_pottemp
export liquid_ice_pottemp_sat
export relative_humidity
export virtual_temperature
export condensate
export has_condensate
export specific_enthalpy
export total_specific_enthalpy
export moist_static_energy
export specific_entropy
export saturated
export q_vap_from_RH_liquid

@inline ReLU(x) = (x > 0) * x

"""
    gas_constant_air(param_set, [q::PhasePartition])
    gas_constant_air(param_set, q_tot, q_liq, q_ice)

The specific gas constant of moist air given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
 - `q_tot`, `q_liq` `q_ice` - specific liquid water contents for total water, liquid water and ice
"""
@inline function gas_constant_air(param_set::APS, q_tot, q_liq, q_ice)
    R_d = TP.R_d(param_set)
    R_v = TP.R_v(param_set)
    q_vap = q_tot - q_liq - q_ice
    return R_d * (1 - q_tot) + R_v * q_vap
end

@inline function gas_constant_air(
    param_set::APS,
    q::PhasePartition{FT},
) where {FT}
    return gas_constant_air(param_set, q.tot, q.liq, q.ice)
end

@inline gas_constant_air(param_set::APS, ::Type{FT}) where {FT} =
    gas_constant_air(param_set, q_pt_0(FT))

"""
    gas_constant_air(param_set::APS, ts::ThermodynamicState)

The specific gas constant of moist air given
a thermodynamic state `ts`.
"""
@inline gas_constant_air(param_set::APS, ts::ThermodynamicState) =
    gas_constant_air(param_set, PhasePartition(param_set, ts))
@inline gas_constant_air(
    param_set::APS,
    ts::AbstractPhaseDry{FT},
) where {FT <: Real} = FT(TP.R_d(param_set))


"""
    air_pressure(param_set, T, ρ[, q::PhasePartition])

The air pressure from the equation of state
(ideal gas law) where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` air temperature
 - `ρ` (moist-)air density
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function air_pressure(
    param_set::APS,
    T::FT,
    ρ::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT <: Real}
    return gas_constant_air(param_set, q) * ρ * T
end
@inline air_pressure(param_set, T, ρ, q) =
    air_pressure(param_set, promote_phase_partition(T, ρ, q)...)

"""
    air_pressure(param_set::APS, ts::ThermodynamicState)

The air pressure from the equation of state
(ideal gas law), given a thermodynamic state `ts`.
"""
@inline air_pressure(param_set::APS, ts::ThermodynamicState) = air_pressure(
    param_set,
    air_temperature(param_set, ts),
    air_density(param_set, ts),
    PhasePartition(param_set, ts),
)

@inline air_pressure(param_set::APS, ts::PhaseEquil) = ts.p

"""
    air_density(param_set, T, p[, q::PhasePartition])

The (moist-)air density from the equation of state
(ideal gas law) where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` air temperature
 - `p` pressure
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function air_density(
    param_set::APS,
    T::FT,
    p::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT <: Real}
    return p / (gas_constant_air(param_set, q) * T)
end
@inline air_density(param_set, T, p, q) =
    air_density(param_set, promote_phase_partition(T, p, q)...)

"""
    air_density(param_set::APS, ts::ThermodynamicState)

The (moist-)air density, given a thermodynamic state `ts`.
"""
@inline air_density(param_set::APS, ts::ThermodynamicState) = ts.ρ

"""
    specific_volume(param_set::APS, ts::ThermodynamicState)

The (moist-)air specific volume, given a thermodynamic state `ts`.
"""
@inline specific_volume(param_set::APS, ts::ThermodynamicState) =
    1 / air_density(param_set, ts)

"""
    total_specific_humidity(param_set::APS, ts::ThermodynamicState)

Total specific humidity given
 - `ts` a thermodynamic state
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for
   more details
"""
@inline total_specific_humidity(param_set::APS, ts::ThermodynamicState) =
    ts.q_tot
@inline total_specific_humidity(
    param_set::APS,
    ts::AbstractPhaseDry{FT},
) where {FT} = FT(0)
@inline total_specific_humidity(param_set::APS, ts::AbstractPhaseNonEquil) =
    ts.q.tot

"""
    liquid_specific_humidity(param_set::APS, ts::ThermodynamicState)
    liquid_specific_humidity(q::PhasePartition)

Liquid specific humidity given
 - `ts` a thermodynamic state
or
 - `q` a `PhasePartition`
"""
@inline liquid_specific_humidity(q::PhasePartition) = q.liq
@inline liquid_specific_humidity(param_set::APS, ts::ThermodynamicState) =
    PhasePartition(param_set, ts).liq
@inline liquid_specific_humidity(
    param_set::APS,
    ts::AbstractPhaseDry{FT},
) where {FT} = FT(0)
@inline liquid_specific_humidity(param_set::APS, ts::AbstractPhaseNonEquil) =
    ts.q.liq

"""
    ice_specific_humidity(param_set::APS, ts::ThermodynamicState)
    ice_specific_humidity(q::PhasePartition)

Ice specific humidity given
 - `ts` a thermodynamic state
or
 - `q` a `PhasePartition`
"""
@inline ice_specific_humidity(q::PhasePartition) = q.ice
@inline ice_specific_humidity(param_set::APS, ts::ThermodynamicState) =
    PhasePartition(param_set, ts).ice
@inline ice_specific_humidity(
    param_set::APS,
    ts::AbstractPhaseDry{FT},
) where {FT} = FT(0)
@inline ice_specific_humidity(param_set::APS, ts::AbstractPhaseNonEquil) =
    ts.q.ice

"""
    vapor_specific_humidity(q::PhasePartition{FT})

The vapor specific humidity, given a `PhasePartition` `q`.
"""
@inline vapor_specific_humidity(q::PhasePartition) =
    max(0, q.tot - q.liq - q.ice)
@inline vapor_specific_humidity(param_set::APS, ts::ThermodynamicState) =
    vapor_specific_humidity(PhasePartition(param_set, ts))

"""
    ∂q_vap_sat_∂T(param_set, ts)

The change in saturation vapor specific humidity with temperature given by the Clausius Clapeyron relation and the definition of specific humidity `q = p_v/(ρ R_v T)` in terms of vapor pressure `p_v`.
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ts` ThermodynamicState
"""
@inline function ∂q_vap_sat_∂T(
    param_set::APS,
    λ::FT,
    T::FT,
    q_vap_sat::FT,
    L = weighted_latent_heat(param_set, T, λ),
) where {FT <: Real}
    R_v = TP.R_v(param_set)
    return q_vap_sat * (L / (R_v * T^2) - 1 / T)
end

@inline function ∂q_vap_sat_∂T(param_set::APS, ts::ThermodynamicState)
    λ = liquid_fraction(param_set, ts)
    T = air_temperature(param_set, ts)
    q_vap_sat = vapor_specific_humidity(param_set, ts)
    return ∂q_vap_sat_∂T(param_set, λ, T, q_vap_sat)
end

"""
    cp_m(param_set, q::PhasePartition)

The isobaric specific heat capacity of moist air given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `qₜ` total specific humidity of water
 - `qₗ` specific humidity of liquid
 - `qᵢ` specific humidity of ice
"""
@inline function cp_m(
    param_set::APS,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
) where {FT <: Real}

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
    cp_m(param_set, q::PhasePartition)

The isobaric specific heat capacity of moist air given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q` [`PhasePartition`](@ref).
"""
@inline function cp_m(param_set::APS, q::PhasePartition{FT}) where {FT <: Real}
    return cp_m(param_set, q.tot, q.liq, q.ice)
end

@inline cp_m(param_set::APS, ::Type{FT}) where {FT <: Real} =
    cp_m(param_set, q_pt_0(FT))

"""
    cp_m(param_set::APS, ts::ThermodynamicState)

The isobaric specific heat capacity of moist air, given a thermodynamic state `ts`.
"""
@inline cp_m(param_set::APS, ts::ThermodynamicState) =
    cp_m(param_set, PhasePartition(param_set, ts))
@inline cp_m(param_set::APS, ts::AbstractPhaseDry{FT}) where {FT <: Real} =
    FT(TP.cp_d(param_set))

"""
    cv_m(param_set, q::PhasePartition)

The isochoric specific heat capacity of moist air given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q` [`PhasePartition`](@ref).
"""
@inline function cv_m(param_set::APS, q::PhasePartition{FT}) where {FT <: Real}
    cv_d = TP.cv_d(param_set)
    cv_v = TP.cv_v(param_set)
    cv_l = TP.cv_l(param_set)
    cv_i = TP.cv_i(param_set)
    return cv_d +
           (cv_v - cv_d) * q.tot +
           (cv_l - cv_v) * q.liq +
           (cv_i - cv_v) * q.ice
end

@inline cv_m(param_set::APS, ::Type{FT}) where {FT <: Real} =
    cv_m(param_set, q_pt_0(FT))

"""
    cv_m(param_set::APS, ts::ThermodynamicState)

The isochoric specific heat capacity of moist air, given a thermodynamic state `ts`.
"""
@inline cv_m(param_set::APS, ts::ThermodynamicState) =
    cv_m(param_set, PhasePartition(param_set, ts))
@inline cv_m(param_set::APS, ts::AbstractPhaseDry{FT}) where {FT <: Real} =
    FT(TP.cv_d(param_set))


"""
    (R_m, cp_m, cv_m, γ_m) = gas_constants(param_set, q::PhasePartition)

Wrapper to compute all gas constants at once, where optionally,
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q` [`PhasePartition`](@ref).

The function returns a tuple of
 - `R_m` [`gas_constant_air`](@ref)
 - `cp_m` [`cp_m`](@ref)
 - `cv_m` [`cv_m`](@ref)
 - `γ_m = cp_m/cv_m`
"""
@inline function gas_constants(
    param_set::APS,
    q::PhasePartition{FT},
) where {FT <: Real}
    R_gas = gas_constant_air(param_set, q)
    cp = cp_m(param_set, q)
    cv = cv_m(param_set, q)
    γ = cp / cv
    return (R_gas, cp, cv, γ)
end

"""
    (R_m, cp_m, cv_m, γ_m) = gas_constants(param_set::APS, ts::ThermodynamicState)

Wrapper to compute all gas constants at once, given a thermodynamic state `ts`.

The function returns a tuple of
 - `R_m` [`gas_constant_air`](@ref)
 - `cp_m` [`cp_m`](@ref)
 - `cv_m` [`cv_m`](@ref)
 - `γ_m = cp_m/cv_m`

"""
@inline gas_constants(param_set::APS, ts::ThermodynamicState) =
    gas_constants(param_set, PhasePartition(param_set, ts))

"""
    air_temperature(param_set, e_int, q::PhasePartition)

The air temperature, where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `e_int` internal energy per unit mass
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function air_temperature(
    param_set::APS,
    e_int::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
    cvm = cv_m(param_set, q),
) where {FT <: Real}
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
    air_temperature_from_enthalpy(param_set, h, q::PhasePartition)

The air temperature, where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `h` internal energy per unit mass
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function air_temperature_from_enthalpy(
    param_set::APS,
    h::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT <: Real}
    cp_m_ = cp_m(param_set, q)
    T_0 = TP.T_0(param_set)
    LH_v0 = TP.LH_v0(param_set)
    LH_f0 = TP.LH_f0(param_set)
    return T_0 + (h - (q.tot - q.liq - q.ice) * LH_v0 + q.ice * LH_f0) / cp_m_
end

"""
    air_temperature(param_set::APS, ts::ThermodynamicState)

The air temperature, given a thermodynamic state `ts`.
"""
@inline air_temperature(param_set::APS, ts::ThermodynamicState) =
    air_temperature(
        param_set,
        internal_energy(param_set, ts),
        PhasePartition(param_set, ts),
    )
@inline air_temperature(param_set::APS, ts::AbstractPhaseEquil) = ts.T

"""
    air_temperature_from_ideal_gas_law(param_set, p, ρ, q::PhasePartition)

The air temperature, where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `ρ` air density
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function air_temperature_from_ideal_gas_law(
    param_set::APS,
    p::FT,
    ρ::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT <: Real}
    R_m = gas_constant_air(param_set, q)
    return p / (R_m * ρ)
end

"""
    internal_energy(param_set, T[, q::PhasePartition])

The internal energy per unit mass, given a thermodynamic state `ts` or

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
    internal_energy(param_set::APS, ts::ThermodynamicState)

The internal energy per unit mass, given a thermodynamic state `ts`.
"""
@inline internal_energy(param_set::APS, ts::ThermodynamicState) = ts.e_int

"""
    internal_energy(ρ::FT, ρe::FT, ρu::AbstractVector{FT}, e_pot::FT)

The internal energy per unit mass, given
 - `ρ` (moist-)air density
 - `ρe` total energy per unit volume
 - `ρu` momentum vector
 - `e_pot` potential energy (e.g., gravitational) per unit mass
"""
@inline function internal_energy(
    ρ::FT,
    ρe::FT,
    ρu::AbstractVector{FT},
    e_pot::FT,
) where {FT <: Real}
    ρinv = 1 / ρ
    ρe_kin = ρinv * sum(abs2, ρu) / 2
    ρe_pot = ρ * e_pot
    ρe_int = ρe - ρe_kin - ρe_pot
    e_int = ρinv * ρe_int
    return e_int
end

"""
    internal_energy_dry(param_set, T)

The dry air internal energy

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
    internal_energy_dry(param_set::APS, ts::ThermodynamicState)

The the dry air internal energy, given a thermodynamic state `ts`.
"""
@inline internal_energy_dry(param_set::APS, ts::ThermodynamicState) =
    internal_energy_dry(param_set, air_temperature(param_set, ts))

"""
    internal_energy_vapor(param_set, T)

The water vapor internal energy

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
    internal_energy_vapor(param_set::APS, ts::ThermodynamicState)

The the water vapor internal energy, given a thermodynamic state `ts`.
"""
@inline internal_energy_vapor(param_set::APS, ts::ThermodynamicState) =
    internal_energy_vapor(param_set, air_temperature(param_set, ts))

"""
    internal_energy_liquid(param_set, T)

The liquid water internal energy

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
    internal_energy_liquid(param_set::APS, ts::ThermodynamicState)

The the liquid water internal energy, given a thermodynamic state `ts`.
"""
@inline internal_energy_liquid(param_set::APS, ts::ThermodynamicState) =
    internal_energy_liquid(param_set, air_temperature(param_set, ts))

"""
    internal_energy_ice(param_set, T)

The ice internal energy

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
    internal_energy_ice(param_set::APS, ts::ThermodynamicState)

The the ice internal energy, given a thermodynamic state `ts`.
"""
@inline internal_energy_ice(param_set::APS, ts::ThermodynamicState) =
    internal_energy_ice(param_set, air_temperature(param_set, ts))

"""
    internal_energy_sat(param_set, T, ρ, q_tot, phase_type)

The internal energy per unit mass in thermodynamic equilibrium at saturation where

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
    cvm = cv_m(param_set, q_pt),
) where {FT <: Real, phase_type <: ThermodynamicState}
    return internal_energy(param_set, T, q_pt)
end
@inline internal_energy_sat(param_set, T, ρ, q_tot, phase_type) =
    internal_energy_sat(param_set, promote(T, ρ, q_tot)..., phase_type)

"""
    internal_energy_sat(param_set::APS, ts::ThermodynamicState)

The internal energy per unit mass in
thermodynamic equilibrium at saturation,
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
 - `e_pot` potential energy per unit mass
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
    return e_kin + e_pot + internal_energy(param_set, T, q)
end

"""
    total_energy(param_set::APS, ts::ThermodynamicState, e_kin, e_pot)

The total energy per unit mass
given a thermodynamic state `ts`.
"""
@inline function total_energy(
    param_set::APS,
    ts::ThermodynamicState{FT},
    e_kin::FT,
    e_pot::FT,
) where {FT <: Real}
    return internal_energy(param_set, ts) + e_kin + e_pot
end

"""
    total_energy_given_ρp(param_set, ρ, p, e_kin, e_pot[, q::PhasePartition])

The total energy per unit mass, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `e_kin` kinetic energy per unit mass
 - `e_pot` potential energy per unit mass
 - `ρ` (moist-)air density
 - `p` pressure
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
    soundspeed_air(param_set, T[, q::PhasePartition])

The speed of sound in unstratified air, where
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function soundspeed_air(
    param_set::APS,
    T::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT <: Real}
    γ = cp_m(param_set, q) / cv_m(param_set, q)
    R_m = gas_constant_air(param_set, q)
    return sqrt(γ * R_m * T)
end

"""
    soundspeed_air(param_set::APS, ts::ThermodynamicState)

The speed of sound in unstratified air given a thermodynamic state `ts`.
"""
@inline soundspeed_air(param_set::APS, ts::ThermodynamicState) = soundspeed_air(
    param_set,
    air_temperature(param_set, ts),
    PhasePartition(param_set, ts),
)


"""
    latent_heat_vapor(param_set, T::FT) where {FT<:Real}

The specific latent heat of vaporization where
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline function latent_heat_vapor(param_set::APS, T::FT) where {FT <: Real}
    cp_l = TP.cp_l(param_set)
    cp_v = TP.cp_v(param_set)
    LH_v0 = TP.LH_v0(param_set)
    return latent_heat_generic(param_set, T, LH_v0, cp_v - cp_l)
end

"""
    latent_heat_vapor(param_set::APS, ts::ThermodynamicState)

The specific latent heat of vaporization
given a thermodynamic state `ts`.
"""
@inline latent_heat_vapor(param_set::APS, ts::ThermodynamicState) =
    latent_heat_vapor(param_set, air_temperature(param_set, ts))

"""
    latent_heat_sublim(param_set, T::FT) where {FT<:Real}

The specific latent heat of sublimation where
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline function latent_heat_sublim(param_set::APS, T::FT) where {FT <: Real}
    LH_s0 = TP.LH_s0(param_set)
    cp_v = TP.cp_v(param_set)
    cp_i = TP.cp_i(param_set)
    return latent_heat_generic(param_set, T, LH_s0, cp_v - cp_i)
end

"""
    latent_heat_sublim(param_set::APS, ts::ThermodynamicState)

The specific latent heat of sublimation
given a thermodynamic state `ts`.
"""
@inline latent_heat_sublim(param_set::APS, ts::ThermodynamicState) =
    latent_heat_sublim(param_set, air_temperature(param_set, ts))

"""
    latent_heat_fusion(param_set, T::FT) where {FT<:Real}

The specific latent heat of fusion where
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
"""
@inline function latent_heat_fusion(param_set::APS, T::FT) where {FT <: Real}
    LH_f0 = TP.LH_f0(param_set)
    cp_l = TP.cp_l(param_set)
    cp_i = TP.cp_i(param_set)
    return latent_heat_generic(param_set, T, LH_f0, cp_l - cp_i)
end

"""
    latent_heat_fusion(param_set::APS, ts::ThermodynamicState)

The specific latent heat of fusion
given a thermodynamic state `ts`.
"""
@inline latent_heat_fusion(param_set::APS, ts::ThermodynamicState) =
    latent_heat_fusion(param_set, air_temperature(param_set, ts))

"""
    latent_heat_generic(param_set, T::FT, LH_0::FT, Δcp::FT) where {FT<:Real}

The specific latent heat of a generic phase transition between
two phases, computed using Kirchhoff's relation with constant
isobaric specific heat capacities of the two phases, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `LH_0` latent heat of the phase transition at `T_0`
 - `Δcp` difference between the isobaric specific heat capacities
         (heat capacity in the higher-temperature phase minus that
         in the lower-temperature phase).
"""
@inline function latent_heat_generic(
    param_set::APS,
    T::FT,
    LH_0::FT,
    Δcp::FT,
) where {FT <: Real}
    T_0 = TP.T_0(param_set)
    return LH_0 + Δcp * (T - T_0)
end
latent_heat_generic(param_set, T, LH_0, Δcp) =
    latent_heat_generic(param_set, promote(T, LH_0, Δcp)...)

"""
    weighted_latent_heat(param_set, T, λ)

Weighted latent heats, computed from
 - `param_set` - a `ThermodynamicsParameters` struct
 - `T` air temperature
 - `λ` liquid fraction
"""
@inline function weighted_latent_heat(
    param_set::APS,
    T::FT,
    λ::FT,
) where {FT <: Real}
    L_v = latent_heat_vapor(param_set, T)
    L_s = latent_heat_sublim(param_set, T)
    return λ * L_v + (1 - λ) * L_s
end

"""
    Phase

A condensed phase, to dispatch over
[`saturation_vapor_pressure`](@ref) and
[`q_vap_saturation_generic`](@ref).
"""
abstract type Phase end

"""
    Liquid <: Phase

A liquid phase, to dispatch over
[`saturation_vapor_pressure`](@ref) and
[`q_vap_saturation_generic`](@ref).
"""
struct Liquid <: Phase end

"""
    Ice <: Phase

An ice phase, to dispatch over
[`saturation_vapor_pressure`](@ref) and
[`q_vap_saturation_generic`](@ref).
"""
struct Ice <: Phase end

"""
    saturation_vapor_pressure(param_set, T, ::Phase)
    saturation_vapor_pressure(param_set, T, LH_0, Δcp)

Return the saturation vapor pressure over a plane surface of condensate, given
 - `T` temperature
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `Phase` either `Liquid()` or `Ice()` to dispatch over the condensate type.
 - `LH_0` latent heat at reference temperature `T_0`
 - `Δcp` difference in isobaric specific heat capacity between the two phases

The saturation vapor pressure is computed by integration of the Clausius-Clapeyron
relation, assuming constant specific heat capacities. The closed-form solution is:
`p_v^*(T) = p_tr * (T/T_tr)^(Δc_p/R_v) * exp((L_0 - Δc_p*T_0)/R_v * (1/T_tr - 1/T))`,
where `p_tr` is the triple point pressure, `T_tr` is the triple point temperature,
`L_0` is the latent heat at the reference temperature `T_0`, and `Δc_p` is the
difference in isobaric specific heat capacities between the phases.
"""
@inline function saturation_vapor_pressure(
    param_set::APS,
    T::FT,
    ::Liquid,
) where {FT <: Real}
    LH_v0 = TP.LH_v0(param_set)
    cp_v = TP.cp_v(param_set)
    cp_l = TP.cp_l(param_set)
    return saturation_vapor_pressure(param_set, T, LH_v0, cp_v - cp_l)
end

@inline function saturation_vapor_pressure(
    param_set::APS,
    ts::ThermodynamicState{FT},
    ::Liquid,
) where {FT <: Real}
    LH_v0 = TP.LH_v0(param_set)
    cp_v = TP.cp_v(param_set)
    cp_l = TP.cp_l(param_set)
    return saturation_vapor_pressure(
        param_set,
        air_temperature(param_set, ts),
        LH_v0,
        cp_v - cp_l,
    )
end

@inline function saturation_vapor_pressure(
    param_set::APS,
    T::FT,
    ::Ice,
) where {FT <: Real}
    LH_s0 = TP.LH_s0(param_set)
    cp_v = TP.cp_v(param_set)
    cp_i = TP.cp_i(param_set)
    return saturation_vapor_pressure(param_set, T, LH_s0, cp_v - cp_i)
end

@inline function saturation_vapor_pressure(
    param_set::APS,
    ts::ThermodynamicState{FT},
    ::Ice,
) where {FT <: Real}
    LH_s0 = TP.LH_s0(param_set)
    cp_v = TP.cp_v(param_set)
    cp_i = TP.cp_i(param_set)
    return saturation_vapor_pressure(
        param_set,
        air_temperature(param_set, ts),
        LH_s0,
        cp_v - cp_i,
    )
end

@inline function saturation_vapor_pressure(
    param_set::APS,
    ::Type{phase_type},
    T::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
    λ = liquid_fraction(param_set, T, phase_type, q),
) where {FT <: Real, phase_type <: ThermodynamicState}

    LH_v0 = TP.LH_v0(param_set)
    LH_s0 = TP.LH_s0(param_set)
    cp_v = TP.cp_v(param_set)
    cp_l = TP.cp_l(param_set)
    cp_i = TP.cp_i(param_set)
    # get phase partitioning

    # effective latent heat at T_0 and effective difference in isobaric specific
    # heats of the mixture
    LH_0 = λ * LH_v0 + (1 - λ) * LH_s0
    Δcp = λ * (cp_v - cp_l) + (1 - λ) * (cp_v - cp_i)

    # saturation vapor pressure over possible mixture of liquid and ice
    return saturation_vapor_pressure(param_set, T, LH_0, Δcp)
end
saturation_vapor_pressure(param_set, T, LH_0, Δcp) =
    saturation_vapor_pressure(param_set, promote(T, LH_0, Δcp)...)

# we may be hitting a slow path:
# https://stackoverflow.com/questions/14687665/very-slow-stdpow-for-bases-very-close-to-1
@inline pow_hack(x, y) = exp(y * log(x))

@inline function saturation_vapor_pressure(
    param_set::APS,
    T::FT,
    LH_0::FT,
    Δcp::FT,
) where {FT <: Real}
    press_triple = TP.press_triple(param_set)
    R_v = TP.R_v(param_set)
    T_triple = TP.T_triple(param_set)
    T_0 = TP.T_0(param_set)

    return press_triple *
           # (T / T_triple)^(Δcp / R_v) *
           pow_hack(T / T_triple, Δcp / R_v) *
           exp((LH_0 - Δcp * T_0) / R_v * (1 / T_triple - 1 / T))

end

"""
    q_vap_saturation_generic(param_set, T, ρ[, phase=Liquid()])

Compute the saturation specific humidity over a plane surface of condensate, given
    - `param_set`: an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
    - `T`: temperature
    - `ρ`: air density
    - (optional) `Liquid()`: indicating condensate is liquid
    - (optional) `Ice()`: indicating condensate is ice

The saturation specific humidity is computed as `q_v^* = p_v^*(T) / (ρ * R_v * T)`,
where `p_v^*` is the saturation vapor pressure.
"""
@inline function q_vap_saturation_generic(
    param_set::APS,
    T::FT,
    ρ::FT,
    phase::Phase,
) where {FT <: Real}
    p_v_sat = saturation_vapor_pressure(param_set, T, phase)
    return q_vap_saturation_from_density(param_set, T, ρ, p_v_sat)
end
q_vap_saturation_generic(param_set::APS, T, ρ, phase::Phase) =
    q_vap_saturation_generic(param_set, promote(T, ρ)..., phase)

"""
    q_vap_saturation(param_set, T, ρ, phase_type[, q::PhasePartition])

Compute the saturation specific humidity, given
- `param_set`: an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
- `T`: temperature
- `ρ`: air density
- `phase_type`: a thermodynamic state type
- (optional) `q`: [`PhasePartition`](@ref)

If the `PhasePartition` `q` is given, the saturation specific humidity is that of a
mixture of liquid and ice, computed in a thermodynamically consistent way from the
weighted sum of the latent heats of the respective phase transitions (Pressel et al.,
JAMES, 2015). That is, the saturation vapor pressure and from it the saturation specific
humidity are computed from a weighted mean of the latent heats of vaporization and
sublimation, with the weights given by the liquid fraction.

If the `PhasePartition` `q` is not given, or has zero liquid and ice specific humidities,
the saturation specific humidity is that over a mixture of liquid and ice, with the
fraction of liquid given by temperature dependent `liquid_fraction(param_set, T, phase_type)`
and the fraction of ice by the complement `1 - liquid_fraction(param_set, T, phase_type)`.
"""
@inline function q_vap_saturation(
    param_set::APS,
    T::FT,
    ρ::FT,
    ::Type{phase_type},
    q::PhasePartition{FT} = q_pt_0(FT),
    λ = liquid_fraction(param_set, T, phase_type, q),
) where {FT <: Real, phase_type <: ThermodynamicState}
    p_v_sat = saturation_vapor_pressure(param_set, phase_type, T, q, λ)
    return q_vap_saturation_from_density(param_set, T, ρ, p_v_sat)
end

"""
    q_vap_saturation(param_set::APS, ts::ThermodynamicState)

Compute the saturation specific humidity, given a thermodynamic state `ts`.
"""
@inline function q_vap_saturation(param_set::APS, ts::ThermodynamicState)
    T = air_temperature(param_set, ts)
    ρ = air_density(param_set, ts)
    q = PhasePartition(param_set, ts)
    p_v_sat = saturation_vapor_pressure(param_set, typeof(ts), T, q)
    return q_vap_saturation_from_density(param_set, T, ρ, p_v_sat)
end

"""
    q_vap_saturation_liquid(param_set::APS, ts::ThermodynamicState)

Compute the saturation specific humidity over liquid,
given a thermodynamic state `ts`.
"""
@inline function q_vap_saturation_liquid(param_set::APS, ts::ThermodynamicState)
    T = air_temperature(param_set, ts)
    ρ = air_density(param_set, ts)
    p_v_sat = saturation_vapor_pressure(param_set, T, Liquid())
    return q_vap_saturation_from_density(param_set, T, ρ, p_v_sat)
end

"""
    q_vap_saturation_ice(param_set::APS, ts::ThermodynamicState)

Compute the saturation specific humidity over ice,
given a thermodynamic state `ts`.
"""
@inline function q_vap_saturation_ice(param_set::APS, ts::ThermodynamicState)
    T = air_temperature(param_set, ts)
    ρ = air_density(param_set, ts)
    p_v_sat = saturation_vapor_pressure(param_set, T, Ice())
    return q_vap_saturation_from_density(param_set, T, ρ, p_v_sat)
end

"""
    q_vap_saturation_from_density(param_set, T, ρ, p_v_sat)

Compute the saturation specific humidity, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature,
 - `ρ` (moist-)air density
 - `p_v_sat` saturation vapor pressure
"""
@inline function q_vap_saturation_from_density(
    param_set::APS,
    T::FT,
    ρ::FT,
    p_v_sat::FT,
) where {FT <: Real}
    R_v = TP.R_v(param_set)
    return p_v_sat / (ρ * R_v * T)
end

"""
    q_vap_saturation_from_pressure(param_set, q_tot, p, T, phase_type)
Compute the saturation specific humidity, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q_tot` total water specific humidity,
 - `p` air pressure,
 - `T` air tempearture
 - `phase_type` a thermodynamic state type

"""
@inline function q_vap_saturation_from_pressure(
    param_set::APS,
    q_tot::FT,
    p::FT,
    T::FT,
    ::Type{phase_type},
    λ = liquid_fraction(param_set, T, phase_type),
) where {FT <: Real, phase_type <: ThermodynamicState}
    R_v = TP.R_v(param_set)
    R_d = TP.R_d(param_set)
    p_v_sat = saturation_vapor_pressure(
        param_set,
        phase_type,
        T,
        PhasePartition(FT(0)),
        λ,
    )
    q_v_sat = if p - p_v_sat ≥ eps(FT)
        R_d / R_v * (1 - q_tot) * p_v_sat / (p - p_v_sat)
    else
        FT(1)
    end
    return q_v_sat
end

"""
    supersaturation(param_set, q, ρ, T, Liquid())
    supersaturation(param_set, q, ρ, T, Ice())
    supersaturation(ts, Ice())
    supersaturation(ts, Liquid())

 - `param_set` - abstract set with earth parameters
 - `q` - phase partition
 - `ρ` - air density,
 - `T` - air temperature
 - `Liquid()`, `Ice()` - liquid or ice phase to dispatch over.
 - `ts` thermodynamic state

Returns supersaturation (pv/pv_sat -1) over water or ice.
"""
@inline function supersaturation(
    param_set::APS,
    q::PhasePartition{FT},
    ρ::FT,
    T::FT,
    phase::Phase,
) where {FT <: Real}

    p_v_sat = saturation_vapor_pressure(param_set, T, phase)

    return supersaturation(param_set, q, ρ, T, p_v_sat)
end
@inline function supersaturation(
    param_set::APS,
    q::PhasePartition{FT},
    ρ::FT,
    T::FT,
    p_v_sat::FT,
) where {FT <: Real}

    p_v = vapor_specific_humidity(q) * (ρ * TP.R_v(param_set) * T)

    return p_v / p_v_sat - FT(1)
end

@inline supersaturation(param_set::APS, ts::ThermodynamicState, phase::Phase) =
    supersaturation(
        param_set,
        PhasePartition(param_set, ts),
        air_density(param_set, ts),
        air_temperature(param_set, ts),
        phase,
    )

"""
    saturation_excess(param_set, T, ρ, phase_type, q::PhasePartition)

The saturation excess in equilibrium where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `phase_type` a thermodynamic state type
 - `q` [`PhasePartition`](@ref)

The saturation excess is the difference between the total specific humidity `q.tot`
and the saturation specific humidity in equilibrium, and it is defined to be
nonzero only if this difference is positive.
"""
@inline function saturation_excess(
    param_set::APS,
    T::FT,
    ρ::FT,
    p_vap_sat::FT,
    q::PhasePartition{FT},
) where {FT <: Real}
    q_vap_sat = q_vap_saturation_from_density(param_set, T, ρ, p_vap_sat)
    return max(0, q.tot - q_vap_sat)
end

@inline function saturation_excess(
    param_set::APS,
    T::FT,
    ρ::FT,
    ::Type{phase_type},
    q::PhasePartition{FT},
    λ = liquid_fraction(param_set, T, phase_type, q),
) where {FT <: Real, phase_type <: ThermodynamicState}
    p_vap_sat = saturation_vapor_pressure(
        param_set,
        phase_type,
        T,
        PhasePartition(FT(0)),
        λ,
    )
    return saturation_excess(param_set, T, ρ, p_vap_sat, q)
end

"""
    saturation_excess(param_set::APS, ts::ThermodynamicState)

Compute the saturation excess in equilibrium,
given a thermodynamic state `ts`.
"""
@inline saturation_excess(param_set::APS, ts::ThermodynamicState) =
    saturation_excess(
        param_set,
        air_temperature(param_set, ts),
        air_density(param_set, ts),
        typeof(ts),
        PhasePartition(param_set, ts),
    )

"""
    condensate(q::PhasePartition{FT})
    condensate(param_set::APS, ts::ThermodynamicState)

Condensate of the phase partition.
"""
@inline condensate(q::PhasePartition) = q.liq + q.ice
@inline condensate(param_set::APS, ts::ThermodynamicState) =
    condensate(PhasePartition(param_set, ts))

"""
    has_condensate(q::PhasePartition{FT})
    has_condensate(param_set::APS, ts::ThermodynamicState)

Bool indicating if condensate exists in the phase
partition
"""
@inline has_condensate(q_c::FT) where {FT <: Real} = q_c > eps(FT)
@inline has_condensate(q::PhasePartition) = has_condensate(condensate(q))
@inline has_condensate(param_set::APS, ts::ThermodynamicState) =
    has_condensate(PhasePartition(param_set, ts))


"""
    liquid_fraction(param_set, T, phase_type[, q])

The fraction of condensate that is liquid, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `phase_type` a thermodynamic state type

# PhaseNonEquil behavior
If `q.liq` or `q.ice` are nonzero, the liquid fraction is computed from
them.

# ThermodynamicState
Otherwise, phase equilibrium is assumed so that the fraction of liquid
is a function that is 1 above `T_freeze` and goes to zero below `T_icenuc`.
"""
@inline function liquid_fraction(
    param_set::APS,
    T::FT,
    ::Type{phase_type},
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT, phase_type <: ThermodynamicState}

    Tᶠ = TP.T_freeze(param_set)   # freezing temperature
    Tⁿ = TP.T_icenuc(param_set)   # temperature of homogeneous ice nucleation
    α = TP.pow_icenuc(param_set) # power law partial ice nucleation parameter

    # Model for cloud water condensate probabilty when temperature is below
    # freezing (all liquid condensate), but above the temperature of homogeneous
    # ice nucleation (all ice condensate). In it's simplest form, this is simply
    # a number that varies between 0 and 1. For example, see figure 6 in
    # Hu et al., https://doi.org/10.1029/2009JD012384, JGR (2010).
    λᵖ = ((T - Tⁿ) / (Tᶠ - Tⁿ))^α

    above_freezing = T > Tᶠ
    supercooled_liquid = (T ≤ Tᶠ) & (T > Tⁿ)

    return ifelse(
        above_freezing,
        one(FT),
        ifelse(supercooled_liquid, λᵖ, zero(FT)),
    )
end

@inline function liquid_fraction(
    param_set::APS,
    T::FT,
    ::Type{phase_type},
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT <: Real, phase_type <: PhaseNonEquil}
    q_c = condensate(q)     # condensate specific humidity
    if has_condensate(q_c)
        return q.liq / q_c
    else
        return liquid_fraction(param_set, T, PhaseEquil, q)
    end
end

"""
    liquid_fraction(param_set::APS, ts::ThermodynamicState)

The fraction of condensate that is liquid given a thermodynamic state `ts`.
"""
@inline liquid_fraction(param_set::APS, ts::ThermodynamicState) =
    liquid_fraction(
        param_set,
        air_temperature(param_set, ts),
        typeof(ts),
        PhasePartition(param_set, ts),
    )

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
@inline function PhasePartition_equil(
    param_set::APS,
    T::FT,
    ρ::FT,
    q_tot::FT,
    p_vap_sat::FT,
    λ::FT,
) where {FT <: Real}
    q_c = saturation_excess(param_set, T, ρ, p_vap_sat, PhasePartition(q_tot)) # condensate specific humidity
    q_liq = λ * q_c                                                            # liquid specific humidity
    q_ice = (1 - λ) * q_c                                                      # ice specific humidity
    return PhasePartition(q_tot, q_liq, q_ice)
end

@inline function PhasePartition_equil(
    param_set::APS,
    T::FT,
    ρ::FT,
    q_tot::FT,
    ::Type{phase_type},
    λ = liquid_fraction(param_set, T, phase_type), # fraction of condensate that is liquid
) where {FT <: Real, phase_type <: ThermodynamicState}
    p_vap_sat = saturation_vapor_pressure(
        param_set,
        phase_type,
        T,
        PhasePartition(FT(0)),
        λ,
    )
    return PhasePartition_equil(param_set, T, ρ, q_tot, p_vap_sat, λ)
end
@inline PhasePartition_equil(param_set, T, ρ, q_tot, phase_type) =
    PhasePartition_equil(param_set, promote(T, ρ, q_tot)..., phase_type)

@inline PhasePartition_equil(param_set::APS, ts::ThermodynamicState) =
    PhasePartition_equil(
        param_set,
        air_temperature(param_set, ts),
        air_density(param_set, ts),
        total_specific_humidity(param_set, ts),
        typeof(ts),
    )


"""
    PhasePartition_equil_given_p(param_set, T, p, q_tot, phase_type)
Partition the phases in equilibrium, returning a [`PhasePartition`](@ref) object using the
[`liquid_fraction`](@ref) function where
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `p` air pressure
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
The residual `q.tot - q.liq - q.ice` is the vapor specific humidity.
"""
@inline function PhasePartition_equil_given_p(
    param_set::APS,
    T::FT,
    p::FT,
    q_tot::FT,
    ::Type{phase_type},
    λ = liquid_fraction(param_set, T, phase_type),
) where {FT <: Real, phase_type <: ThermodynamicState}

    q_v_sat =
        q_vap_saturation_from_pressure(param_set, q_tot, p, T, phase_type, λ)
    q_c = q_tot - q_v_sat
    q_liq = λ * q_c
    q_ice = (1 - λ) * q_c
    return PhasePartition(q_tot, q_liq, q_ice)
end
@inline PhasePartition_equil_given_p(param_set, T, p, q_tot, phase_type) =
    PhasePartition_equil_given_p(param_set, promote(T, p, q_tot)..., phase_type)

@inline PhasePartition(
    param_set::APS,
    ts::AbstractPhaseDry{FT},
) where {FT <: Real} = q_pt_0(FT)
@inline function PhasePartition(param_set::APS, ts::AbstractPhaseEquil)
    T = air_temperature(param_set, ts)
    ρ = air_density(param_set, ts)
    q_tot = total_specific_humidity(param_set, ts)
    phase_type = typeof(ts)
    λ = liquid_fraction(param_set, T, phase_type) # fraction of condensate that is liquid
    qpt0 = PhasePartition(typeof(λ)(0))
    p_vap_sat = saturation_vapor_pressure(param_set, phase_type, T, qpt0, λ)

    return PhasePartition_equil(param_set, T, ρ, q_tot, p_vap_sat, λ)
end
@inline PhasePartition(param_set::APS, ts::AbstractPhaseNonEquil) = ts.q

@inline function ∂e_int_∂T(
    param_set::APS,
    T::FT,
    e_int::FT,
    ρ::FT,
    q_tot::FT,
    ::Type{phase_type},
    λ = liquid_fraction(param_set, T, phase_type),
    p_vap_sat = saturation_vapor_pressure(
        param_set,
        phase_type,
        T,
        PhasePartition(FT(0)),
        λ,
    ),
    q = PhasePartition_equil(param_set, T, ρ, q_tot, p_vap_sat, λ),
    cvm = cv_m(param_set, q),
) where {FT <: Real, phase_type <: PhaseEquil}
    T_0 = TP.T_0(param_set)
    cv_v = TP.cv_v(param_set)
    cv_l = TP.cv_l(param_set)
    cv_i = TP.cv_i(param_set)
    e_int_v0 = TP.e_int_v0(param_set)
    e_int_i0 = TP.e_int_i0(param_set)
    T_f = TP.T_freeze(param_set)
    T_i = TP.T_icenuc(param_set)
    n_i = TP.pow_icenuc(param_set)

    q_c = condensate(q)
    q_vap_sat = q_vap_saturation_from_density(param_set, T, ρ, p_vap_sat)
    L = weighted_latent_heat(param_set, T, λ)

    ∂λ_∂T = if (T_i < T < T_f)
        if n_i == 1
            (1 / (T_f - T_i))
        else
            (1 / (T_f - T_i))^n_i * n_i * T^(n_i - 1)
        end
    else
        FT(0)
    end

    _∂q_vap_sat_∂T = ∂q_vap_sat_∂T(param_set, λ, T, q_vap_sat, L)
    dcvm_dq_vap = cv_v - λ * cv_l - (1 - λ) * cv_i
    return cvm +
           (e_int_v0 + (1 - λ) * e_int_i0 + (T - T_0) * dcvm_dq_vap) *
           _∂q_vap_sat_∂T +
           q_c * e_int_i0 * ∂λ_∂T
end

"""
    saturation_adjustment(
        sat_adjust_method,
        param_set,
        e_int,
        ρ,
        q_tot,
        phase_type,
        maxiter,
        relative_temperature_tol,
        T_guess,
    )

Compute the temperature that is consistent with

 - `sat_adjust_method` the numerical method to use.
    See the [`Thermodynamics`](@ref) for options.
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `e_int` internal energy
 - `ρ` (moist-)air density
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
 - `maxiter` maximum iterations for non-linear equation solve
 - `relative_temperature_tol` relative temperature tolerance
 - `T_guess` initial temperature guess

by finding the root of

`e_int - internal_energy_sat(param_set, T, ρ, q_tot, phase_type) = 0`

using the given numerical method `sat_adjust_method`.

See also [`saturation_adjustment`](@ref).
"""
@inline function saturation_adjustment(
    ::Type{sat_adjust_method},
    param_set::APS,
    e_int::FT,
    ρ::FT,
    q_tot::FT,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::Real,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, sat_adjust_method, phase_type <: PhaseEquil}
    _T_min = TP.T_min(param_set)
    tol = RS.RelativeSolutionTolerance(relative_temperature_tol)

    T_1 = max(_T_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    if T_1 ≥ _T_min
        q_v_sat = q_vap_saturation(param_set, T_1, ρ, phase_type)
        unsaturated = q_tot <= q_v_sat
        if unsaturated
            return T_1
        end
    end
    _T_freeze = TP.T_freeze(param_set)
    @inline e_int_sat(T) =
        internal_energy_sat(param_set, ReLU(T), ρ, q_tot, phase_type)
    temperature_tol = _T_freeze * relative_temperature_tol
    e_int_upper = e_int_sat(_T_freeze + temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    if e_int < e_int_upper
        e_int_lower = e_int_sat(_T_freeze - temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
        if e_int_lower < e_int # < e_int_upper
            return _T_freeze
        end
    end
    @inline function roots(_T) # ff′
        T = ReLU(_T)
        if sat_adjust_method <: RS.NewtonsMethod
            λ = liquid_fraction(param_set, T, phase_type)
            p_vap_sat = saturation_vapor_pressure(
                param_set,
                phase_type,
                T,
                PhasePartition(FT(0)),
                λ,
            )
            q = PhasePartition_equil(param_set, T, ρ, q_tot, p_vap_sat, λ)

            cvm = cv_m(param_set, q)
            f′ = ∂e_int_∂T(
                param_set,
                T,
                e_int,
                ρ,
                q_tot,
                phase_type,
                λ,
                p_vap_sat,
                q,
                cvm,
            )
            _e_int_sat =
                internal_energy_sat(param_set, T, ρ, q_tot, phase_type, q, cvm)
            f = _e_int_sat - e_int
            return (f, f′)
        else
            return e_int_sat(T) - e_int
        end
    end
    sol = RS.find_zero(
        roots,
        sa_numerical_method(
            sat_adjust_method,
            param_set,
            ρ,
            e_int,
            q_tot,
            phase_type,
            T_guess,
        ),
        solution_type(),
        tol,
        maxiter,
    )

    DataCollection.log_meta(sol)
    if !sol.converged
        if print_warning()
            print_warning_ρeq(
                sat_adjust_method,
                ρ,
                e_int,
                q_tot,
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

"""
    saturation_adjustment_given_peq(
        sat_adjust_method,
        param_set,
        e_int,
        p,
        q_tot,
        phase_type,
        maxiter,
        relative_temperature_tol,
        T_guess,
    )

Compute the temperature that is consistent with

 - `sat_adjust_method` the numerical method to use.
    See the [`Thermodynamics`](@ref) for options.
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `e_int` internal energy
 - `p` air pressure
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
 - `maxiter` maximum iterations for non-linear equation solve
 - `relative_temperature_tol` relative temperature tolerance
 - `T_guess` initial temperature guess

by finding the root of

`e_int - internal_energy_sat(param_set, T, ρ(T), q_tot, phase_type) = 0`

where `ρ(T) = air_density(param_set, T, p, PhasePartition(q_tot))`

using the given numerical method `sat_adjust_method`.

See also [`saturation_adjustment`](@ref).
"""
@inline function saturation_adjustment_given_peq(
    ::Type{sat_adjust_method},
    param_set::APS,
    p::FT,
    e_int::FT,
    q_tot::FT,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::Real,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, sat_adjust_method, phase_type <: PhaseEquil}
    _T_min = TP.T_min(param_set)
    tol = RS.RelativeSolutionTolerance(relative_temperature_tol)

    T_1 = max(_T_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    @inline ρ_T(T) = air_density(param_set, T, p, PhasePartition(q_tot))
    ρ_1 = ρ_T(T_1)
    q_v_sat = q_vap_saturation(param_set, T_1, ρ_1, phase_type)
    unsaturated = q_tot <= q_v_sat
    if unsaturated && T_1 ≥ _T_min
        return T_1
    end
    _T_freeze = TP.T_freeze(param_set)
    @inline e_int_sat(T) =
        internal_energy_sat(param_set, ReLU(T), ρ_T(T), q_tot, phase_type)

    temperature_tol = _T_freeze * relative_temperature_tol
    e_int_upper = e_int_sat(_T_freeze + temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    e_int_lower = e_int_sat(_T_freeze - temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    if e_int_lower < e_int < e_int_upper
        return _T_freeze
    end
    @inline roots(T) = e_int_sat(T) - e_int
    sol = RS.find_zero(
        roots,
        sa_numerical_method_peq(
            sat_adjust_method,
            param_set,
            p,
            e_int,
            q_tot,
            phase_type,
            T_guess,
        ),
        RS.CompactSolution(),
        tol,
        maxiter,
    )
    if !sol.converged
        if print_warning()
            print_warning_peq(
                sat_adjust_method,
                p,
                e_int,
                q_tot,
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


"""
    saturation_adjustment_given_phq(
        sat_adjust_method,
        param_set,
        h,
        p,
        q_tot,
        phase_type,
        maxiter,
        relative_temperature_tol
        T_guess,
    )

Compute the temperature that is consistent with

 - `sat_adjust_method` the numerical method to use.
    See the [`Thermodynamics`](@ref) for options.
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `h` specific_enthalpy
 - `p` air pressure
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
 - `maxiter` maximum iterations for non-linear equation solve
 - `relative_temperature_tol` relative temperature tolerance
 - `T_guess` initial temperature guess

by finding the root of

`h - enthalpy_sat(param_set, T, ρ(T), q_tot, phase_type) = 0`

where `ρ(T) = air_density(param_set, T, p, PhasePartition(q_tot))`

using the given numerical method `sat_adjust_method`.

See also [`saturation_adjustment`](@ref).
"""
@inline function saturation_adjustment_given_phq(
    ::Type{sat_adjust_method},
    param_set::APS,
    p::FT,
    h::FT,
    q_tot::FT,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::Real,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, sat_adjust_method, phase_type <: PhaseEquil}
    _T_min = TP.T_min(param_set)
    tol = RS.RelativeSolutionTolerance(relative_temperature_tol)

    T_1 = max(
        _T_min,
        air_temperature_from_enthalpy(param_set, h, PhasePartition(q_tot)),
    ) # Assume all vapor
    @inline ρ_T(T) = air_density(param_set, T, p, PhasePartition(q_tot))
    ρ_1 = ρ_T(T_1)
    q_v_sat = q_vap_saturation(param_set, T_1, ρ_1, phase_type)
    unsaturated = q_tot <= q_v_sat
    if unsaturated && T_1 ≥ _T_min
        return T_1
    end
    _T_freeze = TP.T_freeze(param_set)
    @inline h_sat(T) =
        specific_enthalpy_sat(param_set, ReLU(T), ρ_T(T), q_tot, phase_type)

    temperature_tol = _T_freeze * relative_temperature_tol
    h_upper = h_sat(_T_freeze + temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    h_lower = h_sat(_T_freeze - temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    if h_lower < h < h_upper
        return _T_freeze
    end
    @inline roots(T) = h_sat(T) - h
    sol = RS.find_zero(
        roots,
        sa_numerical_method_phq(
            sat_adjust_method,
            param_set,
            p,
            h,
            q_tot,
            phase_type,
            T_guess,
        ),
        RS.CompactSolution(),
        tol,
        maxiter,
    )
    if !sol.converged
        if print_warning()
            # Must match argument order!
            print_warning_hpq(
                sat_adjust_method,
                h,
                p,
                q_tot,
                T_guess,
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

"""
    saturation_adjustment_ρpq(
        sat_adjust_method,
        param_set,
        ρ,
        p,
        q_tot,
        phase_type,
        maxiter,
        relative_temperature_tol,
        T_guess,
    )
Compute the temperature that is consistent with
 - `sat_adjust_method` the numerical method to use.
    See the [`Thermodynamics`](@ref) for options.
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `p` pressure
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
 - `maxiter` maximum iterations for non-linear equation solve
 - `relative_temperature_tol` relative temperature tolerance
 - `T_guess` initial temperature guess
by finding the root of

```
T - air_temperature_from_ideal_gas_law(
        param_set,
        p,
        ρ,
        PhasePartition_equil(param_set, T, ρ, q_tot, phase_type),
    )
```
using Newtons method using ForwardDiff.
See also [`saturation_adjustment`](@ref).
"""
@inline function saturation_adjustment_ρpq(
    ::Type{sat_adjust_method},
    param_set::APS,
    ρ::FT,
    p::FT,
    q_tot::FT,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::Real,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, sat_adjust_method, phase_type <: PhaseEquil}
    tol = RS.RelativeSolutionTolerance(relative_temperature_tol)
    # Use `oftype` to preserve diagonalized type signatures:
    @inline roots(T) =
        T - air_temperature_from_ideal_gas_law(
            param_set,
            oftype(T, p),
            oftype(T, ρ),
            PhasePartition_equil(
                param_set,
                T,
                oftype(T, ρ),
                oftype(T, q_tot),
                phase_type,
            ),
        )
    sol = RS.find_zero(
        roots,
        sa_numerical_method_ρpq(
            sat_adjust_method,
            param_set,
            ρ,
            p,
            q_tot,
            phase_type,
            T_guess,
        ),
        RS.CompactSolution(),
        tol,
        maxiter,
    )
    if !sol.converged
        if print_warning()
            print_warning_ρpq(
                sat_adjust_method,
                ρ,
                p,
                q_tot,
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

"""
    ΔT_min(::Type{FT})

Minimum interval for saturation adjustment using Secant method
"""
@inline ΔT_min(::Type{FT}) where {FT} = FT(3)

"""
    ΔT_max(::Type{FT})

Maximum interval for saturation adjustment using Secant method
"""
@inline ΔT_max(::Type{FT}) where {FT} = FT(10)

"""
    bound_upper_temperature(T_1::FT, T_2::FT) where {FT<:Real}

Bounds the upper temperature, `T_2`, for
saturation adjustment using Secant method
"""
@inline function bound_upper_temperature(T_1::FT, T_2::FT) where {FT <: Real}
    T_2 = max(T_1 + ΔT_min(FT), T_2)
    return min(T_1 + ΔT_max(FT), T_2)
end

"""
    saturation_adjustment_given_ρθq(
        param_set,
        ρ,
        θ_liq_ice,
        q_tot,
        phase_type,
        maxiter,
        tol,
        T_guess,
    )

Compute the temperature `T` that is consistent with

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
 - `maxiter` maximum iterations for non-linear equation solve
 - `tol` absolute tolerance for saturation adjustment iterations. Can be one of:
    - `RelativeSolutionTolerance()` to stop when `|x_2 - x_1| < tol`
    - `ResidualTolerance()` to stop when `|f(x)| < tol`
    - `RelativeRelativeSolutionTolerance()` to stop when `|x_2 - x_1|/x_1 < tol`
 - `T_guess` initial temperature guess

by finding the root of

`θ_{liq_ice} - liquid_ice_pottemp_sat(param_set, T, ρ, phase_type, q_tot) = 0`

See also [`saturation_adjustment`](@ref).
"""
@inline function saturation_adjustment_given_ρθq(
    param_set::APS,
    ρ::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    ::Type{phase_type},
    maxiter::Int,
    tol::RS.AbstractTolerance,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, phase_type <: PhaseEquil}
    _T_min = TP.T_min(param_set)
    T_init_min = TP.T_init_min(param_set)
    @inline air_temp(q) = air_temperature_given_ρθq(param_set, ρ, θ_liq_ice, q)
    T_1 = max(_T_min, air_temp(PhasePartition(q_tot))) # Assume all vapor
    q_v_sat = q_vap_saturation(param_set, T_1, ρ, phase_type)
    unsaturated = q_tot <= q_v_sat
    if unsaturated && T_1 ≥ _T_min
        return T_1
    end
    T_2 = air_temp(PhasePartition(q_tot, FT(0), q_tot)) # Assume all ice
    T_2 = bound_upper_temperature(T_1, T_2)
    @inline roots(T) =
        liquid_ice_pottemp_sat(param_set, ReLU(T), ρ, phase_type, q_tot) -
        θ_liq_ice
    sol = RS.find_zero(
        roots,
        RS.SecantMethod(T_init_min, T_2),
        RS.CompactSolution(),
        tol,
        maxiter,
    )
    if !sol.converged
        if print_warning()
            print_warning_ρθq(
                RS.SecantMethod,
                ρ,
                θ_liq_ice,
                q_tot,
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

"""
    saturation_adjustment_given_pθq(
        sat_adjust_method,
        param_set,
        p,
        θ_liq_ice,
        q_tot,
        phase_type,
        maxiter,
        relative_temperature_tol,
        T_guess
    )

Compute the temperature `T` that is consistent with

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
 - `relative_temperature_tol` relative temperature tolerance
 - `maxiter` maximum iterations for non-linear equation solve
 - `sat_adjust_method` the numerical method to use.
 - `T_guess` initial temperature guess

by finding the root of

`θ_{liq_ice} - liquid_ice_pottemp_given_pressure(param_set, T, p, phase_type, q_tot) = 0`

See also [`saturation_adjustment`](@ref).
"""
@inline function saturation_adjustment_given_pθq(
    ::Type{sat_adjust_method},
    param_set::APS,
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::Real,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, sat_adjust_method, phase_type <: PhaseEquil}
    tol = RS.RelativeSolutionTolerance(relative_temperature_tol)
    T_min = TP.T_min(param_set)
    T_freeze = TP.T_freeze(param_set)
    cp_d = TP.cp_d(param_set)
    cp_v = TP.cp_v(param_set)
    @inline air_temp(q) = air_temperature_given_pθq(param_set, p, θ_liq_ice, q)
    @inline function θ_liq_ice_closure(T)
        q = PhasePartition(oftype(T, 0))
        λ = liquid_fraction(param_set, T, phase_type, q)
        q_pt = PhasePartition_equil_given_p(
            param_set,
            T,
            oftype(T, p),
            oftype(T, q_tot),
            phase_type,
            λ,
        )
        return liquid_ice_pottemp_given_pressure(
            param_set,
            T,
            oftype(T, p),
            q_pt,
        )
    end
    @inline q_vap_sat(T) =
        q_vap_saturation_from_pressure(param_set, q_tot, p, T, phase_type)
    T_1 = max(T_min, air_temp(PhasePartition(q_tot))) # Assume all vapor
    q_v_sat_1 = q_vap_sat(T_1)
    unsaturated = q_tot <= q_v_sat_1
    if unsaturated && T_1 ≥ T_min
        return T_1
    end
    temperature_tol = T_freeze * relative_temperature_tol
    θ_liq_ice_upper = θ_liq_ice_closure(T_freeze + temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    θ_liq_ice_lower = θ_liq_ice_closure(T_freeze - temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    if θ_liq_ice_lower < θ_liq_ice < θ_liq_ice_upper
        return T_freeze
    end
    roots(T) = oftype(T, θ_liq_ice) - θ_liq_ice_closure(T)
    sol = RS.find_zero(
        roots,
        sa_numerical_method_pθq(
            sat_adjust_method,
            param_set,
            p,
            θ_liq_ice,
            q_tot,
            phase_type,
            T_guess,
        ),
        RS.CompactSolution(),
        tol,
        maxiter,
    )
    if !sol.converged
        if print_warning()
            print_warning_pθq(
                sat_adjust_method,
                p,
                θ_liq_ice,
                q_tot,
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

"""
    latent_heat_liq_ice(param_set, q::PhasePartition{FT})

Effective latent heat of condensate (weighted sum of liquid and ice),
with specific latent heat evaluated at reference temperature `T_0` given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function latent_heat_liq_ice(
    param_set::APS,
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT <: Real}
    LH_v0 = TP.LH_v0(param_set)
    LH_s0 = TP.LH_s0(param_set)
    return LH_v0 * q.liq + LH_s0 * q.ice
end
latent_heat_liq_ice(param_set::APS, ts::ThermodynamicState) =
    latent_heat_liq_ice(param_set, PhasePartition(param_set, ts))


"""
    liquid_ice_pottemp_given_pressure(param_set, T, p[, q::PhasePartition, cpm])

The liquid-ice potential temperature where
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `p` pressure
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function liquid_ice_pottemp_given_pressure(
    param_set::APS,
    T::FT,
    p::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
    cpm = cp_m(param_set, q),
) where {FT <: Real}
    # liquid-ice potential temperature, approximating latent heats
    # of phase transitions as constants
    return dry_pottemp_given_pressure(param_set, T, p, q, cpm) *
           (1 - latent_heat_liq_ice(param_set, q) / (cpm * T))
end


"""
    liquid_ice_pottemp(param_set, T, ρ[, q::PhasePartition, cpm])

The liquid-ice potential temperature, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function liquid_ice_pottemp(
    param_set::APS,
    T::FT,
    ρ::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
    cpm = cp_m(param_set, q),
) where {FT <: Real}
    return liquid_ice_pottemp_given_pressure(
        param_set,
        T,
        air_pressure(param_set, T, ρ, q),
        q,
        cpm,
    )
end

"""
    liquid_ice_pottemp(param_set::APS, ts::ThermodynamicState)

The liquid-ice potential temperature,
given a thermodynamic state `ts`.
"""
@inline liquid_ice_pottemp(param_set::APS, ts::ThermodynamicState) =
    liquid_ice_pottemp(
        param_set,
        air_temperature(param_set, ts),
        air_density(param_set, ts),
        PhasePartition(param_set, ts),
    )

"""
    dry_pottemp(param_set, T, ρ[, q::PhasePartition, cpm])

The dry potential temperature where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
 """
@inline function dry_pottemp(
    param_set::APS,
    T::FT,
    ρ::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
    cpm = cp_m(param_set, q),
) where {FT <: Real}
    return T / exner(param_set, T, ρ, q, cpm)
end

"""
    dry_pottemp_given_pressure(param_set, T, p[, q::PhasePartition])

The dry potential temperature where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `p` pressure
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
 """
@inline function dry_pottemp_given_pressure(
    param_set::APS,
    T::FT,
    p::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
    cpm = cp_m(param_set, q),
) where {FT <: Real}
    return T / exner_given_pressure(param_set, p, q, cpm)
end

"""
    dry_pottemp(param_set::APS, ts::ThermodynamicState)

The dry potential temperature, given a thermodynamic state `ts`.
"""
@inline dry_pottemp(param_set::APS, ts::ThermodynamicState) = dry_pottemp(
    param_set,
    air_temperature(param_set, ts),
    air_density(param_set, ts),
    PhasePartition(param_set, ts),
)

@inline function virt_temp_from_RH(
    param_set::APS,
    T::FT,
    ρ::FT,
    RH::FT,
    ::Type{phase_type},
) where {FT <: Real, phase_type <: ThermodynamicState}
    q_tot = RH * q_vap_saturation(param_set, T, ρ, phase_type)
    q_pt = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    return virtual_temperature(param_set, T, q_pt)
end
"""
    temperature_and_humidity_given_TᵥρRH(param_set, T_virt, ρ, RH)

The air temperature and `q_tot` where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T_virt` virtual temperature
 - `ρ` air density
 - `RH` relative humidity
 - `phase_type` a thermodynamic state type
"""
@inline function temperature_and_humidity_given_TᵥρRH(
    param_set::APS,
    T_virt::FT,
    ρ::FT,
    RH::FT,
    ::Type{phase_type},
    maxiter::Int = 100,
    tol::RS.AbstractTolerance = RS.ResidualTolerance{FT}(sqrt(eps(FT))),
) where {FT <: Real, phase_type <: ThermodynamicState}

    T_init_min = TP.T_init_min(param_set)
    _T_max = T_virt
    @inline roots(T) =
        T_virt - virt_temp_from_RH(param_set, ReLU(T), ρ, RH, phase_type)
    sol = RS.find_zero(
        roots,
        RS.SecantMethod(T_init_min, _T_max),
        RS.CompactSolution(),
        tol,
        maxiter,
    )
    if !sol.converged
        if print_warning()
            print_warning_TᵥρRH(
                RS.SecantMethod,
                T_virt,
                RH,
                ρ,
                sol.root,
                maxiter,
                tol.tol,
            )
        end
        if error_on_non_convergence()
            error("Failed to converge with printed set of inputs.")
        end
    end
    T = sol.root

    # Re-compute specific humidity and phase partitioning
    # given the temperature
    q_tot = RH * q_vap_saturation(param_set, T, ρ, phase_type)
    q_pt = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    return (T, q_pt)

end

"""
    air_temperature_given_ρθq(param_set, ρ, θ_liq_ice, q::PhasePartition)

The temperature given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `θ_liq_ice` liquid-ice potential temperature
 - `ρ` (moist-)air density
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function air_temperature_given_ρθq(
    param_set::APS,
    ρ::FT,
    θ_liq_ice::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT <: Real}

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
    air_temperature_given_ρθq_nonlinear(param_set, ρ, θ_liq_ice, q::PhasePartition)

Computes temperature `T` given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `θ_liq_ice` liquid-ice potential temperature
 - `ρ` (moist-)air density
 - `tol` absolute tolerance for non-linear equation iterations. Can be one of:
    - `RelativeSolutionTolerance()` to stop when `|x_2 - x_1|/x_1 < tol`
    - `ResidualTolerance()` to stop when `|f(x)| < tol`
 - `maxiter` maximum iterations for non-linear equation solve
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air,

by finding the root of
`T - air_temperature_given_pθq(param_set,
                               air_pressure(param_set, T, ρ, q),
                               θ_liq_ice,
                               q) = 0`
"""
@inline function air_temperature_given_ρθq_nonlinear(
    param_set::APS,
    ρ::FT,
    θ_liq_ice::FT,
    maxiter::Int,
    tol::RS.AbstractTolerance,
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT <: Real}
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

"""
    air_temperature_given_pθq(
        param_set,
        p,
        θ_liq_ice,
        [q::PhasePartition]
    )

The air temperature where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `θ_liq_ice` liquid-ice potential temperature
 - `p` pressure
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function air_temperature_given_pθq(
    param_set::APS,
    p::FT,
    θ_liq_ice::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
    cpm = cp_m(param_set, q),
) where {FT <: Real}
    return θ_liq_ice * exner_given_pressure(param_set, p, q, cpm) +
           latent_heat_liq_ice(param_set, q) / cpm
end

"""
    virtual_pottemp(param_set, T, ρ[, q::PhasePartition])

The virtual potential temperature, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.

The virtual potential temperature is defined as `θ_v = (R_m/R_d) * θ`, where `θ` is the
potential temperature. It is the temperature a dry air parcel would need to have to
have the same density as the moist air parcel at the same pressure.
"""
@inline function virtual_pottemp(
    param_set::APS,
    T::FT,
    ρ::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
    cpm = cp_m(param_set, q),
) where {FT <: Real}
    R_d = TP.R_d(param_set)
    return gas_constant_air(param_set, q) / R_d *
           dry_pottemp(param_set, T, ρ, q, cpm)
end

"""
    virtual_pottemp(param_set::APS, ts::ThermodynamicState)

The virtual potential temperature,
given a thermodynamic state `ts`.
"""
@inline virtual_pottemp(param_set::APS, ts::ThermodynamicState) =
    virtual_pottemp(
        param_set,
        air_temperature(param_set, ts),
        air_density(param_set, ts),
        PhasePartition(param_set, ts),
    )

"""
    virtual_temperature(param_set, T[, q::PhasePartition])

The virtual temperature, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.

The virtual temperature is defined as `T_v = (R_m/R_d) * T`. It is the temperature
a dry air parcel would need to have to have the same density as the moist air parcel
at the same pressure.
"""
@inline function virtual_temperature(
    param_set::APS,
    T::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT <: Real}
    R_d = TP.R_d(param_set)
    return gas_constant_air(param_set, q) / R_d * T
end

"""
    virtual_temperature(param_set::APS, ts::ThermodynamicState)

The virtual temperature,
given a thermodynamic state `ts`.
"""
@inline virtual_temperature(param_set::APS, ts::ThermodynamicState) =
    virtual_temperature(
        param_set,
        air_temperature(param_set, ts),
        PhasePartition(param_set, ts),
    )


"""
    liquid_ice_pottemp_sat(param_set, T, ρ, phase_type[, q::PhasePartition, cpm])

The saturated liquid ice potential temperature where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `phase_type` a thermodynamic state type
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function liquid_ice_pottemp_sat(
    param_set::APS,
    T::FT,
    ρ::FT,
    ::Type{phase_type},
    q::PhasePartition{FT} = q_pt_0(FT),
    cpm = cp_m(param_set, q),
) where {FT <: Real, phase_type <: ThermodynamicState}
    q_v_sat = q_vap_saturation(param_set, T, ρ, phase_type, q)
    return liquid_ice_pottemp(param_set, T, ρ, PhasePartition(q_v_sat), cpm)
end

"""
    liquid_ice_pottemp_sat(param_set, T, ρ, phase_type, q_tot)

The saturated liquid ice potential temperature where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `phase_type` a thermodynamic state type
 - `q_tot` total specific humidity
"""
@inline function liquid_ice_pottemp_sat(
    param_set::APS,
    T::FT,
    ρ::FT,
    ::Type{phase_type},
    q_tot::FT,
) where {FT <: Real, phase_type <: ThermodynamicState}
    q = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    cpm = cp_m(param_set, q)
    return liquid_ice_pottemp(param_set, T, ρ, q, cpm)
end

"""
    liquid_ice_pottemp_sat(param_set::APS, ts::ThermodynamicState)

The liquid potential temperature given a thermodynamic state `ts`.
"""
@inline liquid_ice_pottemp_sat(param_set::APS, ts::ThermodynamicState) =
    liquid_ice_pottemp_sat(
        param_set,
        air_temperature(param_set, ts),
        air_density(param_set, ts),
        typeof(ts),
        PhasePartition(param_set, ts),
    )

"""
    exner_given_pressure(param_set, p[, q::PhasePartition, cpm])

The Exner function, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function exner_given_pressure(
    param_set::APS,
    p::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
    cpm = cp_m(param_set, q),
) where {FT <: Real}
    p0 = TP.p_ref_theta(param_set)
    # gas constant and isobaric specific heat of moist air
    _R_m = gas_constant_air(param_set, q)

    # return (p / p0)^(_R_m / cpm)
    return pow_hack(p / p0, _R_m / cpm)
end

"""
    exner(param_set, T, ρ[, q::PhasePartition)])

The Exner function, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
"""
@inline function exner(
    param_set::APS,
    T::FT,
    ρ::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
    cpm = cp_m(param_set, q),
) where {FT <: Real}
    p = air_pressure(param_set, T, ρ, q)
    return exner_given_pressure(param_set, p, q, cpm)
end

"""
    exner(param_set::APS, ts::ThermodynamicState)

The Exner function, given a thermodynamic state `ts`.
"""
@inline exner(param_set::APS, ts::ThermodynamicState) = exner(
    param_set,
    air_temperature(param_set, ts),
    air_density(param_set, ts),
    PhasePartition(param_set, ts),
)

"""
    shum_to_mixing_ratio(q, q_tot)

Mixing ratio, from specific humidity
 - `q` specific humidity
 - `q_tot` total specific humidity
"""
@inline function shum_to_mixing_ratio(q::FT, q_tot::FT) where {FT <: Real}
    return q / (1 - q_tot)
end

"""
    mixing_ratios(q::PhasePartition)

Mixing ratios
 - `r.tot` total mixing ratio
 - `r.liq` liquid mixing ratio
 - `r.ice` ice mixing ratio
given a phase partition, `q`.
"""
@inline function mixing_ratios(q::PhasePartition{FT}) where {FT <: Real}
    return PhasePartition(
        shum_to_mixing_ratio(q.tot, q.tot),
        shum_to_mixing_ratio(q.liq, q.tot),
        shum_to_mixing_ratio(q.ice, q.tot),
    )
end

"""
    mixing_ratios(param_set::APS, ts::ThermodynamicState)

Mixing ratios stored, in a phase partition, for
 - total specific humidity
 - liquid specific humidity
 - ice specific humidity
"""
@inline mixing_ratios(param_set::APS, ts::ThermodynamicState) =
    mixing_ratios(PhasePartition(param_set, ts))

"""
    vol_vapor_mixing_ratio(param_set, q::PhasePartition)

Volume mixing ratio of water vapor
given a parameter set `param_set`
and a phase partition, `q`.
"""
@inline function vol_vapor_mixing_ratio(
    param_set::APS,
    q::PhasePartition{FT},
) where {FT <: Real}
    molmass_ratio = TP.molmass_ratio(param_set)
    q_vap = vapor_specific_humidity(q)
    return molmass_ratio * shum_to_mixing_ratio(q_vap, q.tot)
end
vol_vapor_mixing_ratio(param_set, ts::ThermodynamicState) =
    vol_vapor_mixing_ratio(param_set, PhasePartition(param_set, ts))

"""
    relative_humidity(param_set, T, p, phase_type, q::PhasePartition)

The relative humidity (clipped between 0 and 1), given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `phase_type` a thermodynamic state type
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument, the results are for dry air.
 """
@inline function relative_humidity(
    param_set::APS,
    T::FT,
    p::FT,
    ::Type{phase_type},
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT <: Real, phase_type <: ThermodynamicState}
    R_v = TP.R_v(param_set)
    q_vap = vapor_specific_humidity(q)
    p_vap = q_vap * air_density(param_set, T, p, q) * R_v * T
    p_vap_sat = saturation_vapor_pressure(param_set, phase_type, T)
    return max(FT(0), min(FT(1), p_vap / (p_vap_sat + eps(FT(0)))))
end

"""
    relative_humidity(param_set::APS, ts::ThermodynamicState)

The relative humidity, given a thermodynamic state `ts`.
"""
@inline relative_humidity(
    param_set::APS,
    ts::ThermodynamicState{FT},
) where {FT <: Real} = relative_humidity(
    param_set,
    air_temperature(param_set, ts),
    air_pressure(param_set, ts),
    typeof(ts),
    PhasePartition(param_set, ts),
)

@inline relative_humidity(
    param_set::APS,
    ts::AbstractPhaseDry{FT},
) where {FT <: Real} = FT(0)

"""
    total_specific_enthalpy(e_tot, R_m, T)

Total specific enthalpy, given
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

Total specific enthalpy, given
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
    specific_enthalpy(e_int, R_m, T)

Specific enthalpy, given
 - `e_int` internal specific energy
 - `R_m` [`gas_constant_air`](@ref)
 - `T` air temperature
"""
@inline function specific_enthalpy(e_int::FT, R_m::FT, T::FT) where {FT <: Real}
    return e_int + R_m * T
end

"""
    specific_enthalpy(param_set, ts)

Specific enthalpy, given a thermodynamic state `ts`.
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

The specific_enthalpy per unit mass, given a thermodynamic state `ts` or

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

The enthalpy per unit mass in thermodynamic equilibrium at saturation where

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

Moist static energy, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ts` a thermodynamic state
 - `e_pot` potential energy (e.g., gravitational) per unit mass
"""
@inline function moist_static_energy(
    param_set::APS,
    ts::ThermodynamicState,
    e_pot,
)
    return specific_enthalpy(param_set, ts) + e_pot
end

"""
    virtual_dry_static_energy(param_set, ts, e_pot)

Virtual dry static energy, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ts` a thermodynamic state
 - `e_pot` gravitational potential energy per unit mass
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

"""
    specific_entropy(param_set, p, T, q)
    specific_entropy(param_set, ts)

Specific entropy, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `q` phase partition

following equations (29)-(33) of [Pressel2015](@cite).
"""
@inline function specific_entropy(
    param_set::APS,
    p::FT,
    T::FT,
    q::PhasePartition{FT},
) where {FT <: Real}
    L_v = latent_heat_vapor(param_set, T)
    L_s = latent_heat_sublim(param_set, T)
    s_d = specific_entropy_dry(param_set, p, T, q)
    s_v = specific_entropy_vapor(param_set, p, T, q)
    return (1 - q.tot) * s_d + q.tot * s_v - (q.liq * L_v + q.ice * L_s) / T
end

@inline specific_entropy(param_set::APS, ts::ThermodynamicState) =
    specific_entropy(
        param_set,
        air_pressure(param_set, ts),
        air_temperature(param_set, ts),
        PhasePartition(param_set, ts),
    )

"""
    specific_entropy_dry(param_set, p, T, q)

The dry air specific entropy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `q` phase partition
"""
@inline function specific_entropy_dry(
    param_set::APS,
    p::FT,
    T::FT,
    q::PhasePartition{FT},
) where {FT <: Real}
    T_ref = TP.entropy_reference_temperature(param_set)
    p_ref = TP.MSLP(param_set)
    s_d_ref = TP.entropy_dry_air(param_set)
    R_d = TP.R_d(param_set)
    cp_d = TP.cp_d(param_set)
    p_d = partial_pressure_dry(param_set, p, q)
    return s_d_ref + cp_d * log(T / T_ref) - R_d * log((p_d + eps(FT)) / p_ref)
end

"""
    specific_entropy_vapor(param_set, p, T, q)

The specific entropy of water vapor, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `q` phase partition
"""
@inline function specific_entropy_vapor(
    param_set::APS,
    p::FT,
    T::FT,
    q::PhasePartition{FT},
) where {FT <: Real}
    T_ref = TP.entropy_reference_temperature(param_set)
    p_ref = TP.MSLP(param_set)
    s_v_ref = TP.entropy_water_vapor(param_set)
    R_v = TP.R_v(param_set)
    cp_v = TP.cp_v(param_set)
    p_v = partial_pressure_vapor(param_set, p, q)
    return s_v_ref + cp_v * log(T / T_ref) - R_v * log((p_v + eps(FT)) / p_ref)
end

"""
    partial_pressure_dry(param_set, p, q)

The partial pressure of water vapor, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `q` phase partition
"""
@inline function partial_pressure_dry(
    param_set::APS,
    p::FT,
    q::PhasePartition{FT},
) where {FT <: Real}
    molmass_ratio = TP.molmass_ratio(param_set)
    return p * (1 - q.tot) /
           (1 - q.tot + vapor_specific_humidity(q) * molmass_ratio)
end

"""
    partial_pressure_vapor(param_set, p, q)

The partial pressure of water vapor, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `q` phase partition
"""
@inline function partial_pressure_vapor(
    param_set::APS,
    p::FT,
    q::PhasePartition{FT},
) where {FT <: Real}
    molmass_ratio = TP.molmass_ratio(param_set)
    return p * vapor_specific_humidity(q) * molmass_ratio /
           (1 - q.tot + vapor_specific_humidity(q) * molmass_ratio)
end

"""
    saturated(param_set::APS, ts::ThermodynamicState)

Boolean indicating if thermodynamic
state is saturated.
"""
@inline function saturated(param_set::APS, ts::ThermodynamicState)
    RH = relative_humidity(param_set, ts)
    return RH ≈ 1 || RH > 1
end


"""
    q_vap_from_RH_liquid(param_set, p, T, RH)

Computes the water vapor specific humidity from relative humidity with
  regard to water (i.e. saturation vapor pressure over ice is not considered
  even in cold temperatures).

Inputs:
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `RH` relative humidity
"""
@inline function q_vap_from_RH_liquid(
    param_set::APS,
    p::FT,
    T::FT,
    RH::FT,
) where {FT <: Real}
    @assert RH <= FT(1)
    p_vap_sat = saturation_vapor_pressure(param_set, T, Liquid())
    p_vap = RH * p_vap_sat
    mmr = TP.molmass_ratio(param_set)
    return p_vap / mmr / (p - (1 - 1 / mmr) * p_vap)
end
