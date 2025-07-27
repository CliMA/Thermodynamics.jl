export saturation_vapor_pressure
export q_vap_saturation_generic
export q_vap_saturation
export q_vap_saturation_liquid
export q_vap_saturation_ice
export saturation_excess
export supersaturation
export saturated

"""
    saturation_vapor_pressure(param_set, T, ::Phase)
    saturation_vapor_pressure(param_set, T, LH_0, Δcp)

The saturation vapor pressure over a plane surface of condensate, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature or `ts` a thermodynamic state
 - `Phase` either `Liquid()` or `Ice()` to dispatch over the condensate type

 or, given

 - `param_set` an `AbstractParameterSet`
 - `T` temperature  
 - `LH_0` latent heat at reference temperature `T_0`
 - `Δcp` difference in isobaric specific heat capacity between the two phases

The saturation vapor pressure is computed by integration of the Clausius-Clapeyron
relation, assuming constant specific heat capacities. The closed-form solution is:
`p_v^*(T) = p_tr * (T/T_tr)^(Δc_p/R_v) * exp((L_0 - Δc_p*T_0)/R_v * (1/T_tr - 1/T))`,
where `p_tr` is the triple point pressure, `T_tr` is the triple point temperature,
`L_0` is the latent heat at the reference temperature `T_0`, and `Δc_p` is the
difference in isobaric specific heat capacities between the phases.
"""
@inline function saturation_vapor_pressure(param_set::APS, T, ::Liquid)
    LH_v0 = TP.LH_v0(param_set)
    cp_v = TP.cp_v(param_set)
    cp_l = TP.cp_l(param_set)
    return saturation_vapor_pressure(param_set, T, LH_v0, cp_v - cp_l)
end

@inline function saturation_vapor_pressure(param_set::APS, T, ::Ice)
    LH_s0 = TP.LH_s0(param_set)
    cp_v = TP.cp_v(param_set)
    cp_i = TP.cp_i(param_set)
    return saturation_vapor_pressure(param_set, T, LH_s0, cp_v - cp_i)
end

@inline function saturation_vapor_pressure(param_set::APS, T, LH_0, Δcp)
    press_triple = TP.press_triple(param_set)
    R_v = TP.R_v(param_set)
    T_triple = TP.T_triple(param_set)
    T_0 = TP.T_0(param_set)

    return press_triple *
           # (T / T_triple)^(Δcp / R_v) *
           fast_power(T / T_triple, Δcp / R_v) *
           exp((LH_0 - Δcp * T_0) / R_v * (1 / T_triple - 1 / T))
end

@inline function saturation_vapor_pressure(
    param_set::APS,
    ::Type{phase_type},
    T::Number,
    q::PhasePartition = q_pt_0(param_set),
    λ = liquid_fraction(param_set, T, phase_type, q),
) where {phase_type <: ThermodynamicState}

    LH_v0 = TP.LH_v0(param_set)
    LH_s0 = TP.LH_s0(param_set)
    cp_v = TP.cp_v(param_set)
    cp_l = TP.cp_l(param_set)
    cp_i = TP.cp_i(param_set)

    # effective latent heat at T_0 and effective difference in isobaric specific
    # heats of the mixture
    LH_0 = λ * LH_v0 + (1 - λ) * LH_s0
    Δcp = λ * (cp_v - cp_l) + (1 - λ) * (cp_v - cp_i)

    # saturation vapor pressure over possible mixture of liquid and ice
    return saturation_vapor_pressure(param_set, T, LH_0, Δcp)
end
saturation_vapor_pressure(param_set, T, LH_0, Δcp) =
    saturation_vapor_pressure(param_set, promote(T, LH_0, Δcp)...)

"""
    q_vap_saturation_generic(param_set, T, ρ[, phase=Liquid()])

The saturation specific humidity over a plane surface of condensate, given
    - `param_set`: an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
    - `T`: temperature
    - `ρ`: air density
    - (optional) `Liquid()`: indicating condensate is liquid (default)
    - (optional) `Ice()`: indicating condensate is ice

The saturation specific humidity is computed as `q_v^* = p_v^*(T) / (ρ * R_v * T)`,
where `p_v^*` is the saturation vapor pressure.
"""
@inline function q_vap_saturation_generic(
    param_set::APS,
    T,
    ρ,
    phase::Phase = Liquid(),
)
    p_v_sat = saturation_vapor_pressure(param_set, T, phase)
    return q_vap_from_p_vap(param_set, T, ρ, p_v_sat)
end

"""
    q_vap_saturation(param_set, T, ρ, phase_type[, q::PhasePartition])

The saturation specific humidity, given

- `param_set`: an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
- `T`: temperature
- `ρ`: air density
- `phase_type`: a thermodynamic state type
- (optional) `q`: [`PhasePartition`](@ref)

If the `PhasePartition` `q` is given, the saturation specific humidity is that over a
mixture of liquid and ice, computed in a thermodynamically consistent way from the
weighted sum of the latent heats of the respective phase transitions (Pressel et al.,
JAMES, 2015). That is, the saturation vapor pressure and from it the saturation specific
humidity are computed from a weighted mean of the latent heats of vaporization and
sublimation, with the weights given by the liquid fraction.

If the `PhasePartition` `q` is not given, or has zero liquid and ice specific humidities,
the saturation specific humidity is that over a mixture of liquid and ice, with the
fraction of liquid given by the temperature dependent `liquid_fraction(param_set, T, phase_type)`
and the fraction of ice by the complement `1 - liquid_fraction(param_set, T, phase_type)`.
"""
@inline function q_vap_saturation(
    param_set::APS,
    T,
    ρ,
    ::Type{phase_type},
    q::PhasePartition = q_pt_0(param_set),
    λ = liquid_fraction(param_set, T, phase_type, q),
) where {phase_type <: ThermodynamicState}
    p_v_sat = saturation_vapor_pressure(param_set, phase_type, T, q, λ)
    return q_vap_from_p_vap(param_set, T, ρ, p_v_sat)
end

"""
    q_vap_saturation_from_pressure(param_set, q_tot, p, T, phase_type)

The saturation specific humidity, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q_tot` total water specific humidity,
 - `p` air pressure,
 - `T` air tempearture
 - `phase_type` a thermodynamic state type

"""
@inline function q_vap_saturation_from_pressure(
    param_set::APS{FT},
    q_tot,
    p,
    T,
    ::Type{phase_type},
    λ = liquid_fraction(param_set, T, phase_type),
) where {FT, phase_type <: ThermodynamicState}
    R_v = TP.R_v(param_set)
    R_d = TP.R_d(param_set)
    p_v_sat = saturation_vapor_pressure(
        param_set,
        phase_type,
        T,
        PhasePartition(FT(0)),
        λ,
    )
    q_v_sat = ifelse(
        p - p_v_sat ≥ eps(FT),
        R_d / R_v * (1 - q_tot) * p_v_sat / (p - p_v_sat),
        FT(1),
    )
    return q_v_sat
end

"""
    supersaturation(param_set, q, ρ, T, Liquid())
    supersaturation(param_set, q, ρ, T, Ice())

The supersaturation (pv/pv_sat -1) over water or ice, given
 - `param_set` - abstract set with earth parameters
 - `q` - phase partition
 - `ρ` - air density,
 - `T` - air temperature
 - `Liquid()`, `Ice()` - liquid or ice phase to dispatch over.
"""
@inline function supersaturation(
    param_set::APS,
    q::PhasePartition,
    ρ,
    T,
    phase::Phase,
)

    p_v_sat = saturation_vapor_pressure(param_set, T, phase)

    return supersaturation(param_set, q, ρ, T, p_v_sat)
end

"""
    supersaturation(param_set, q::PhasePartition, ρ, T, p_v_sat)

The supersaturation (pv/pv_sat - 1) given the saturation vapor pressure `p_v_sat`:

- `param_set`: Thermodynamic parameter set
- `q`: Phase partition
- `ρ`: Air density
- `T`: Temperature
- `p_v_sat`: Saturation vapor pressure

"""
@inline function supersaturation(
    param_set::APS,
    q::PhasePartition,
    ρ,
    T,
    p_v_sat,
)

    p_v = vapor_specific_humidity(q) * (ρ * TP.R_v(param_set) * T)

    return p_v / p_v_sat - 1
end

"""
    saturation_excess(param_set, T, ρ, phase_type, q::PhasePartition)

The saturation excess in equilibrium, given

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
    param_set::APS{FT},
    T,
    ρ,
    ::Type{phase_type},
    q::PhasePartition,
    λ = liquid_fraction(param_set, T, phase_type, q),
) where {FT, phase_type <: ThermodynamicState}
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
    saturation_excess(param_set, T, ρ, p_vap_sat, q::PhasePartition)

The saturation excess given the saturation vapor pressure `p_vap_sat`:

- `param_set`: Thermodynamic parameter set
- `T`: Temperature
- `ρ`: Air density
- `p_vap_sat`: Saturation vapor pressure
- `q`: Phase partition

The saturation excess is the difference between the total specific humidity `q.tot`
and the saturation specific humidity, and it is defined to be nonzero only if
this difference is positive.
"""
@inline function saturation_excess(
    param_set::APS,
    T,
    ρ,
    p_vap_sat,
    q::PhasePartition,
)
    q_vap_sat = q_vap_from_p_vap(param_set, T, ρ, p_vap_sat)
    return max(0, q.tot - q_vap_sat)
end
