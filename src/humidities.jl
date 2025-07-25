export total_specific_humidity
export liquid_specific_humidity
export ice_specific_humidity
export vapor_specific_humidity
export partial_pressure_vapor
export partial_pressure_dry
export vapor_pressure_deficit
export shum_to_mixing_ratio
export mixing_ratios
export vol_vapor_mixing_ratio
export relative_humidity
export q_vap_from_p_vap
export q_vap_from_RH_liquid
export q_vap_saturation_from_density  # TODO Remove after ClimaAtmos and ClimaLand are updated to use q_vap_from_p_vap


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
    vapor_specific_humidity(q::PhasePartition{FT})

The vapor specific humidity, given a 

- `q` a `PhasePartition` 
"""
@inline vapor_specific_humidity(q::PhasePartition) =
    max(0, q.tot - q.liq - q.ice)

"""
    shum_to_mixing_ratio(q, q_tot)

The mixing ratio, given
 - `q` specific humidity
 - `q_tot` total specific humidity
"""
@inline function shum_to_mixing_ratio(q::FT, q_tot::FT) where {FT <: Real}
    return q / (1 - q_tot)
end

"""
    mixing_ratios(q::PhasePartition)

The mixing ratios, given a specific humidity phase partition, `q`, returned in a 
`PhasePartition` with the fields
 - `r.tot` total mixing ratio
 - `r.liq` liquid mixing ratio
 - `r.ice` ice mixing ratio
"""
@inline function mixing_ratios(q::PhasePartition{FT}) where {FT <: Real}
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
@inline function vol_vapor_mixing_ratio(
    param_set::APS,
    q::PhasePartition{FT},
) where {FT <: Real}
    Rv_over_Rd = TP.Rv_over_Rd(param_set)
    q_vap = vapor_specific_humidity(q)
    return Rv_over_Rd * shum_to_mixing_ratio(q_vap, q.tot)
end

"""
    partial_pressure_dry(param_set, p, q)

The partial pressure of dry air, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `q` phase partition
"""
@inline function partial_pressure_dry(
    param_set::APS,
    p::FT,
    q::PhasePartition{FT},
) where {FT <: Real}
    Rv_over_Rd = TP.Rv_over_Rd(param_set)
    return p * (1 - q.tot) /
           (1 - q.tot + vapor_specific_humidity(q) * Rv_over_Rd)
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
    Rv_over_Rd = TP.Rv_over_Rd(param_set)
    return p * vapor_specific_humidity(q) * Rv_over_Rd /
           (1 - q.tot + vapor_specific_humidity(q) * Rv_over_Rd)
end

"""
    vapor_pressure_deficit(param_set, T, p, q::PhasePartition)

The vapor pressure deficit (saturation vapor pressure minus actual 
vapor pressure, truncated to be non-negative) over liquid water for temperatures 
above freezing and over ice for temperatures below freezing, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` air temperature
 - `p` air pressure
 - `q` [`PhasePartition`](@ref)
"""
@inline function vapor_pressure_deficit(
    param_set::APS,
    T::FT,
    p::FT,
    q::PhasePartition{FT},
    Tᶠ = TP.T_freeze(param_set),
) where {FT <: Real}
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
    vapor_pressure_deficit(param_set, T, p, q_vap)

The vapor pressure deficit over liquid water (saturation vapor pressure minus actual 
vapor pressure, truncated to be non-negative), given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` air temperature
 - `p` air pressure
 - `q_vap` vapor specific humidity
"""
@inline function vapor_pressure_deficit(
    param_set::APS,
    T::FT,
    p::FT,
    q_vap::FT,
) where {FT <: Real}
    # Create a PhasePartition with only vapor and call the existing method
    q = PhasePartition(q_vap)
    return vapor_pressure_deficit(param_set, T, p, q)
end

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
    q_vap_from_p_vap(param_set, T, ρ, p_v)

The vapor specific humidity, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature,
 - `ρ` (moist-)air density
 - `p_v` partial pressure of vapor
"""
@inline function q_vap_from_p_vap(
    param_set::APS,
    T::FT,
    ρ::FT,
    p_v::FT,
) where {FT <: Real}
    R_v = TP.R_v(param_set)
    return p_v / (ρ * R_v * T)
end

"""
    q_vap_saturation_from_density(param_set, T, ρ, p_v)

This function is identical to `q_vap_from_p_vap` and is provided for backward compatibility. 
It will be removed in a future release.
"""
const q_vap_saturation_from_density = q_vap_from_p_vap  # TODO Remove after ClimaAtmos and ClimaLand are updated to use q_vap_from_p_vap

"""
    q_vap_from_RH_liquid(param_set, p, T, RH)

The water vapor specific humidity, given 

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `RH` relative humidity with respect to liquid water
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
    _Rv_over_Rd = TP.Rv_over_Rd(param_set)
    return p_vap / _Rv_over_Rd / (p - (1 - 1 / _Rv_over_Rd) * p_vap)
end
