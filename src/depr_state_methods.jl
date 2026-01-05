# This file contains all methods that take a thermodynamic state (ts) as input

export specific_volume
export total_specific_humidity
export liquid_specific_humidity
export ice_specific_humidity
export mixing_ratios
export saturated

"""
    gas_constant_air(param_set, ts::ThermodynamicState)

The specific gas constant of moist air given a thermodynamic state `ts`.
"""
@inline gas_constant_air(param_set::APS, ts::ThermodynamicState) =
    gas_constant_air(param_set, PhasePartition(param_set, ts))
@inline gas_constant_air(
    param_set::APS,
    ts::AbstractPhaseDry{FT},
) where {FT <: Real} = FT(TP.R_d(param_set))

"""
    cp_m(param_set, ts::ThermodynamicState)

The isobaric specific heat capacity of moist air, given a thermodynamic state `ts`.
"""
@inline cp_m(param_set::APS, ts::ThermodynamicState) =
    cp_m(param_set, PhasePartition(param_set, ts))
@inline cp_m(param_set::APS, ts::AbstractPhaseDry{FT}) where {FT <: Real} =
    FT(TP.cp_d(param_set))

"""
    cv_m(param_set, ts::ThermodynamicState)

The isochoric specific heat capacity of moist air, given a thermodynamic state `ts`.
"""
@inline cv_m(param_set::APS, ts::ThermodynamicState) =
    cv_m(param_set, PhasePartition(param_set, ts))
@inline cv_m(param_set::APS, ts::AbstractPhaseDry{FT}) where {FT <: Real} =
    FT(TP.cv_d(param_set))

"""
    latent_heat_vapor(param_set, ts::ThermodynamicState)

The specific latent heat of vaporization, given a thermodynamic state `ts`.
"""
@inline latent_heat_vapor(param_set::APS, ts::ThermodynamicState) =
    latent_heat_vapor(param_set, air_temperature(param_set, ts))

"""
    latent_heat_sublim(param_set, ts::ThermodynamicState)

The specific latent heat of sublimation, given a thermodynamic state `ts`.
"""
@inline latent_heat_sublim(param_set::APS, ts::ThermodynamicState) =
    latent_heat_sublim(param_set, air_temperature(param_set, ts))

"""
    latent_heat_fusion(param_set, ts::ThermodynamicState)

The specific latent heat of fusion, given a thermodynamic state `ts`.
"""
@inline latent_heat_fusion(param_set::APS, ts::ThermodynamicState) =
    latent_heat_fusion(param_set, air_temperature(param_set, ts))

"""
    soundspeed_air(param_set, ts::ThermodynamicState)

The speed of sound in unstratified air, given a thermodynamic state `ts`.
"""
@inline soundspeed_air(param_set::APS, ts::ThermodynamicState) = soundspeed_air(
    param_set,
    air_temperature(param_set, ts),
    PhasePartition(param_set, ts),
)

"""
    air_temperature(param_set, ts::ThermodynamicState)

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
    air_pressure(param_set, ts::ThermodynamicState)

The air pressure from the equation of state (ideal gas law), given a thermodynamic state `ts`.
"""
@inline air_pressure(param_set::APS, ts::ThermodynamicState) = air_pressure(
    param_set,
    air_temperature(param_set, ts),
    air_density(param_set, ts),
    PhasePartition(param_set, ts),
)

@inline air_pressure(param_set::APS, ts::PhaseEquil) = ts.p

"""
    air_density(param_set, ts::ThermodynamicState)

The (moist-)air density, given a thermodynamic state `ts`.
"""
@inline air_density(param_set::APS, ts::ThermodynamicState) = ts.ρ

"""
    specific_volume(param_set, ts::ThermodynamicState)

The (moist-)air specific volume, given a thermodynamic state `ts`.
"""
@inline specific_volume(param_set::APS, ts::ThermodynamicState) =
    1 / air_density(param_set, ts)

"""
    exner(param_set, ts::ThermodynamicState)

The Exner function, given a thermodynamic state `ts`.
"""
@inline exner(param_set::APS, ts::ThermodynamicState) = exner(
    param_set,
    air_temperature(param_set, ts),
    air_density(param_set, ts),
    PhasePartition(param_set, ts),
)

"""
    dry_pottemp(param_set, ts::ThermodynamicState)

The dry potential temperature, given a thermodynamic state `ts`.
"""
@inline dry_pottemp(param_set::APS, ts::ThermodynamicState) = dry_pottemp(
    param_set,
    air_temperature(param_set, ts),
    air_density(param_set, ts),
    PhasePartition(param_set, ts),
)

"""
    entropy(param_set, ts)

The specific entropy, given a thermodynamic state `ts`.
"""
@inline entropy(param_set::APS, ts::ThermodynamicState) =
    entropy(
        param_set,
        air_pressure(param_set, ts),
        air_temperature(param_set, ts),
        PhasePartition(param_set, ts),
    )

"""
    humidity_weighted_latent_heat(param_set::APS, ts::ThermodynamicState)

Specific-humidity weighted sum of latent heats of liquid and ice evaluated at reference temperature `T_0`, given a thermodynamic state `ts`.
"""
humidity_weighted_latent_heat(param_set::APS, ts::ThermodynamicState) =
    humidity_weighted_latent_heat(param_set, PhasePartition(param_set, ts))

"""
    liquid_ice_pottemp(param_set, ts::ThermodynamicState)

The liquid-ice potential temperature, given a thermodynamic state `ts`.
"""
@inline liquid_ice_pottemp(param_set::APS, ts::ThermodynamicState) =
    liquid_ice_pottemp(
        param_set,
        air_temperature(param_set, ts),
        air_density(param_set, ts),
        PhasePartition(param_set, ts),
    )

"""
    liquid_ice_pottemp_sat(param_set, ts::ThermodynamicState)

The liquid potential temperature, given a thermodynamic state `ts`.
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
    virtual_temperature(param_set, ts::ThermodynamicState)

The virtual temperature, given a thermodynamic state `ts`.
"""
@inline virtual_temperature(param_set::APS, ts::ThermodynamicState) =
    virtual_temperature(
        param_set,
        air_temperature(param_set, ts),
        PhasePartition(param_set, ts),
    )

"""
    virtual_pottemp(param_set, ts::ThermodynamicState)

The virtual potential temperature, given a thermodynamic state `ts`.
"""
@inline virtual_pottemp(param_set::APS, ts::ThermodynamicState) =
    virtual_pottemp(
        param_set,
        air_temperature(param_set, ts),
        air_density(param_set, ts),
        PhasePartition(param_set, ts),
    )

"""
    internal_energy(param_set, ts::ThermodynamicState)

The internal energy per unit mass, given a thermodynamic state `ts`.
"""
@inline internal_energy(param_set::APS, ts::ThermodynamicState) = ts.e_int

"""
    internal_energy_dry(param_set, ts::ThermodynamicState)

The dry air internal energy, given a thermodynamic state `ts`.
"""
@inline internal_energy_dry(param_set::APS, ts::ThermodynamicState) =
    internal_energy_dry(param_set, air_temperature(param_set, ts))

"""
    internal_energy_vapor(param_set, ts::ThermodynamicState)

The water vapor internal energy, given a thermodynamic state `ts`.
"""
@inline internal_energy_vapor(param_set::APS, ts::ThermodynamicState) =
    internal_energy_vapor(param_set, air_temperature(param_set, ts))

"""
    internal_energy_liquid(param_set, ts::ThermodynamicState)

The liquid water internal energy, given a thermodynamic state `ts`.
"""
@inline internal_energy_liquid(param_set::APS, ts::ThermodynamicState) =
    internal_energy_liquid(param_set, air_temperature(param_set, ts))

"""
    internal_energy_ice(param_set, ts::ThermodynamicState)

The ice internal energy, given a thermodynamic state `ts`.
"""
@inline internal_energy_ice(param_set::APS, ts::ThermodynamicState) =
    internal_energy_ice(param_set, air_temperature(param_set, ts))

"""
    internal_energy_sat(param_set, ts::ThermodynamicState)

The internal energy per unit mass in thermodynamic equilibrium at saturation with a fixed temperature and total specific humidity, given a thermodynamic state `ts`.
"""
@inline internal_energy_sat(param_set::APS, ts::ThermodynamicState) =
    internal_energy_sat(
        param_set,
        air_temperature(param_set, ts),
        air_density(param_set, ts),
        total_specific_humidity(param_set, ts),
    )

"""
    total_energy(param_set, ts::ThermodynamicState, e_kin, e_pot)

The total energy per unit mass, given a thermodynamic state `ts`.
"""
@inline function total_energy(
    param_set::APS,
    ts::ThermodynamicState,
    e_kin,
    e_pot,
)
    return internal_energy(param_set, ts) + e_pot + e_kin
end

"""
    enthalpy(param_set, ts)

The specific enthalpy, given a thermodynamic state `ts`.
"""
@inline function enthalpy(
    param_set::APS,
    ts::ThermodynamicState{FT},
) where {FT <: Real}
    e_int = internal_energy(param_set, ts)
    R_m = gas_constant_air(param_set, ts)
    T = air_temperature(param_set, ts)
    return enthalpy(e_int, R_m, T)
end

"""
    moist_static_energy(param_set, ts, e_pot)

The moist static energy, given a thermodynamic state `ts`.
"""
@inline function moist_static_energy(
    param_set::APS,
    ts::ThermodynamicState,
    e_pot,
)
    return enthalpy(param_set, ts) + e_pot
end

"""
    total_enthalpy(param_set, ts, e_tot::Real)

The total specific enthalpy, given a thermodynamic state `ts`.
"""
@inline function total_enthalpy(
    param_set::APS,
    ts::ThermodynamicState,
    e_tot,
)
    R_m = gas_constant_air(param_set, ts)
    T = air_temperature(param_set, ts)
    return total_enthalpy(e_tot, R_m, T)
end

"""
    virtual_dry_static_energy(param_set, ts, e_pot)

The virtual dry static energy, given a thermodynamic state `ts`.
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
    total_specific_humidity(param_set, ts::ThermodynamicState)

The total specific humidity, given a thermodynamic state `ts`.
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
    liquid_specific_humidity(param_set, ts::ThermodynamicState)

The liquid specific humidity, given a thermodynamic state `ts`.
"""
@inline liquid_specific_humidity(param_set::APS, ts::ThermodynamicState) =
    PhasePartition(param_set, ts).liq
@inline liquid_specific_humidity(
    param_set::APS,
    ts::AbstractPhaseDry{FT},
) where {FT} = FT(0)
@inline liquid_specific_humidity(param_set::APS, ts::AbstractPhaseNonEquil) =
    ts.q.liq

"""
    ice_specific_humidity(param_set, ts::ThermodynamicState)

The ice specific humidity, given a thermodynamic state `ts`.
"""
@inline ice_specific_humidity(param_set::APS, ts::ThermodynamicState) =
    PhasePartition(param_set, ts).ice
@inline ice_specific_humidity(
    param_set::APS,
    ts::AbstractPhaseDry{FT},
) where {FT} = FT(0)
@inline ice_specific_humidity(param_set::APS, ts::AbstractPhaseNonEquil) =
    ts.q.ice

"""
    vapor_specific_humidity(param_set, ts::ThermodynamicState)

The vapor specific humidity, given a thermodynamic state `ts`.
"""
@inline vapor_specific_humidity(param_set::APS, ts::ThermodynamicState) =
    vapor_specific_humidity(PhasePartition(param_set, ts))

"""
    mixing_ratios(param_set, ts::ThermodynamicState)

The mixing ratios, stored in a `PhasePartition`, given a thermodynamic state `ts`.
"""
@inline mixing_ratios(param_set::APS, ts::ThermodynamicState) =
    mixing_ratios(PhasePartition(param_set, ts))

"""
    vol_vapor_mixing_ratio(param_set, ts::ThermodynamicState)

The volume mixing ratio of water vapor, given a thermodynamic state `ts`.
"""
vol_vapor_mixing_ratio(param_set::APS, ts::ThermodynamicState) =
    vol_vapor_mixing_ratio(param_set, PhasePartition(param_set, ts))

"""
    partial_pressure_dry(param_set, ts::ThermodynamicState)

The partial pressure of dry air, given a thermodynamic state `ts`.
"""
@inline partial_pressure_dry(
    param_set::APS,
    ts::ThermodynamicState{FT},
) where {FT <: Real} = partial_pressure_dry(
    param_set,
    air_pressure(param_set, ts),
    PhasePartition(param_set, ts),
)

@inline partial_pressure_dry(
    param_set::APS,
    ts::AbstractPhaseDry{FT},
) where {FT <: Real} = air_pressure(param_set, ts)

"""
    partial_pressure_vapor(param_set, ts::ThermodynamicState)

The partial pressure of water vapor, given a thermodynamic state `ts`.
"""
@inline partial_pressure_vapor(
    param_set::APS,
    ts::ThermodynamicState{FT},
) where {FT <: Real} = partial_pressure_vapor(
    param_set,
    air_pressure(param_set, ts),
    PhasePartition(param_set, ts),
)




@inline partial_pressure_vapor(
    param_set::APS,
    ts::AbstractPhaseDry{FT},
) where {FT <: Real} = FT(0)

"""
    relative_humidity(param_set, ts::ThermodynamicState)

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
    saturation_vapor_pressure(param_set, ts::ThermodynamicState, ::Liquid)

The saturation vapor pressure over a plane surface of condensate, given a thermodynamic state `ts` and phase `Liquid`.
"""
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

"""
    saturation_vapor_pressure(param_set, ts::ThermodynamicState, ::Ice)

The saturation vapor pressure over a plane surface of condensate, given a thermodynamic state `ts` and phase `Ice`.
"""
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

"""
    q_vap_saturation(param_set, ts::ThermodynamicState)

The saturation specific humidity, given a thermodynamic state `ts`.
"""
@inline function q_vap_saturation(param_set::APS, ts::ThermodynamicState)
    T = air_temperature(param_set, ts)
    ρ = air_density(param_set, ts)
    q = PhasePartition(param_set, ts)
    return q_vap_saturation(param_set, T, ρ, q.liq, q.ice)
end

@inline function q_vap_saturation(
    param_set::APS,
    T,
    ρ,
    ::Type{phase_type},
) where {phase_type <: ThermodynamicState}
    return q_vap_saturation(param_set, T, ρ)
end

@inline function q_vap_saturation_from_pressure(
    param_set::APS,
    q_tot,
    p,
    T,
    ::Type{phase_type},
) where {phase_type <: ThermodynamicState}
    return q_vap_saturation_from_pressure(param_set, q_tot, p, T)
end



"""
    saturation_excess(param_set, ts::ThermodynamicState)

The saturation excess in equilibrium, given a thermodynamic state `ts`.
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
    supersaturation(param_set::APS, ts::ThermodynamicState, phase::Phase)

The supersaturation (pv/pv_sat -1) over water or ice, given a thermodynamic state `ts`.
"""
@inline supersaturation(param_set::APS, ts::ThermodynamicState, phase::Phase) =
    supersaturation(
        param_set,
        PhasePartition(param_set, ts),
        air_density(param_set, ts),
        air_temperature(param_set, ts),
        phase,
    )

"""
    saturated(param_set, ts::ThermodynamicState)

Checks if the thermodynamic state `ts` is saturated or supersaturated.
"""
@inline function saturated(param_set::APS, ts::ThermodynamicState)
    FT = eltype(param_set)
    return supersaturation(param_set, ts, Liquid()) >= -sqrt(eps(FT)) ||
           supersaturation(param_set, ts, Ice()) >= -sqrt(eps(FT))
end

"""
    condensate_specific_humidity(param_set, ts::ThermodynamicState)

The condensate specific humidity (liquid + ice) given a thermodynamic state `ts`.
"""
@inline condensate_specific_humidity(param_set::APS, ts::ThermodynamicState) =
    condensate_specific_humidity(PhasePartition(param_set, ts))

"""
    has_condensate(param_set, ts::ThermodynamicState)

Bool indicating if condensate exists given a thermodynamic state `ts`.
"""
@inline has_condensate(param_set::APS, ts::ThermodynamicState) =
    has_condensate(PhasePartition(param_set, ts))

"""
    liquid_fraction(param_set, ts::ThermodynamicState)

The fraction of condensate that is liquid, given a thermodynamic state `ts`.
"""
@inline liquid_fraction(param_set::APS, ts::ThermodynamicState) =
    liquid_fraction(
        param_set,
        air_temperature(param_set, ts),
        typeof(ts),
        PhasePartition(param_set, ts),
    )

"""
    PhasePartition_equil(param_set, ts::ThermodynamicState)

Partition the phases in equilibrium, returning a [`PhasePartition`](@ref) object, given a thermodynamic state `ts`.
"""
@inline PhasePartition_equil(param_set::APS, ts::ThermodynamicState) =
    PhasePartition_equil(
        param_set,
        air_temperature(param_set, ts),
        air_density(param_set, ts),
        total_specific_humidity(param_set, ts),
        typeof(ts),
    )

"""
    PhasePartition(param_set, ts::AbstractPhaseDry)

Create a PhasePartition for dry air thermodynamic state.
"""
@inline PhasePartition(param_set::APS, ts::AbstractPhaseDry) = q_pt_0(param_set)

"""
    PhasePartition(param_set, ts::AbstractPhaseEquil)

Create a PhasePartition for equilibrium thermodynamic state.
"""
@inline function PhasePartition(param_set::APS, ts::AbstractPhaseEquil)
    FT = eltype(param_set)
    T = air_temperature(param_set, ts)
    ρ = air_density(param_set, ts)
    q_tot = total_specific_humidity(param_set, ts)
    phase_type = typeof(ts)
    λ = liquid_fraction(param_set, T) # fraction of condensate that is liquid
    p_vap_sat = saturation_vapor_pressure(param_set, T)

    return PhasePartition_equil(param_set, T, ρ, q_tot, p_vap_sat, λ)
end

"""
    PhasePartition(param_set, ts::AbstractPhaseNonEquil)

Create a PhasePartition for non-equilibrium thermodynamic state.
"""
@inline PhasePartition(param_set::APS, ts::AbstractPhaseNonEquil) = ts.q

@inline function ∂q_vap_sat_∂T(param_set::APS, ts::ThermodynamicState)
    λ = liquid_fraction(param_set, ts)
    T = air_temperature(param_set, ts)
    q_vap_sat = vapor_specific_humidity(param_set, ts)
    return ∂q_vap_sat_∂T(param_set, λ, T, q_vap_sat)
end
