# Phase partition functions

export Liquid, Ice
export liquid_fraction, PhasePartition_equil
export has_condensate

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
    has_condensate(q::PhasePartition{FT})

Bool indicating if condensate exists in the phase partition
"""
@inline has_condensate(q_c::FT) where {FT} = q_c > eps(FT)
@inline has_condensate(q::PhasePartition) =
    has_condensate(condensate_specific_humidity(q))

"""
    liquid_fraction(param_set, T, phase_type[, q::PhasePartition])

The fraction of condensate that is liquid, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `phase_type` a thermodynamic state type

# PhaseNonEquil behavior
If `q.liq` or `q.ice` are nonzero, the liquid fraction is computed from
them.

# PhaseEquil, PhaseDry behavior
Otherwise, the liquid fraction goes from 0 below `T_icenuc` to 1 above `T_freeze`,
with a power law interpolation between the two temperatures based on Kaul et al., Monthly 
Weather Rev., 2015, https://doi.org/10.1029/2009JD012384
"""
@inline function liquid_fraction(
    param_set::APS{FT},
    T,
    ::Type{phase_type},
    q::PhasePartition = q_pt_0(param_set),
) where {FT, phase_type <: ThermodynamicState}

    # Interpolation between homogeneous nucleation and freezing temperatures
    Tᶠ = TP.T_freeze(param_set)   # freezing temperature
    Tⁱ = TP.T_icenuc(param_set)   # temperature of homogeneous ice nucleation
    n = TP.pow_icenuc(param_set)  # power law partial ice nucleation parameter
    λᵖ = ((T - Tⁱ) / (Tᶠ - Tⁱ))^n

    above_freezing = T > Tᶠ
    supercooled_liquid = (T ≤ Tᶠ) & (T > Tⁱ)

    return ifelse(
        above_freezing,
        one(FT),
        ifelse(supercooled_liquid, λᵖ, zero(FT)),
    )
end

@inline function liquid_fraction(
    param_set::APS,
    T,
    ::Type{phase_type},
    q::PhasePartition = q_pt_0(param_set),
) where {phase_type <: PhaseNonEquil}
    q_c = condensate_specific_humidity(q)     # condensate specific humidity
    return ifelse(
        has_condensate(q_c),
        q.liq / q_c,
        liquid_fraction(param_set, T, PhaseEquil, q),
    )
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

@inline function PhasePartition_equil(
    param_set::APS{FT},
    T,
    ρ,
    q_tot,
    ::Type{phase_type},
    λ = liquid_fraction(param_set, T, phase_type), # fraction of condensate that is liquid
) where {FT, phase_type <: ThermodynamicState}
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

    q_v_sat =
        q_vap_saturation_from_pressure(param_set, q_tot, p, T, phase_type, λ)
    q_c = q_tot - q_v_sat
    q_liq = λ * q_c
    q_ice = (1 - λ) * q_c
    return PhasePartition(q_tot, q_liq, q_ice)
end
@inline PhasePartition_equil_given_p(param_set, T, p, q_tot, phase_type) =
    PhasePartition_equil_given_p(param_set, promote(T, p, q_tot)..., phase_type)
