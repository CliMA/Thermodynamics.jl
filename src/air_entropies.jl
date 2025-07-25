# Entropy functions (for thermodynamic equilibrium)

export specific_entropy

"""
    specific_entropy(param_set, p, T, q)
    specific_entropy(param_set, ts)

The specific entropy, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `q` phase partition

 or a thermodynamic state `ts`.

The specific entropy is computed from equations (29)-(33) of [Pressel2015](@cite).
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
