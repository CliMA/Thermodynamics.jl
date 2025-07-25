# Material properties of moist air 

# Specific heats and gas constants of moist air
export cp_m, cv_m, gas_constant_air, gas_constants

# Latent heats
export latent_heat_vapor
export latent_heat_sublim
export latent_heat_fusion
export latent_heat_liq_ice

# Speed of sound
export soundspeed_air

"""
    gas_constant_air(param_set, [q::PhasePartition])
    gas_constant_air(param_set, q_tot, q_liq, q_ice)

The specific gas constant of moist air, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q` [`PhasePartition`](@ref). Without humidity argument, the results are for dry air.
 - `q_tot`, `q_liq`, `q_ice` - specific humidities for total water, liquid water and ice
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
    cp_m(param_set, q_tot, q_liq, q_ice)
    cp_m(param_set[, q::PhasePartition])

The isobaric specific heat capacity of moist air, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q` [`PhasePartition`](@ref). Without humidity argument, the results are for dry air.
 - `q_tot` total specific humidity of water
 - `q_liq` specific humidity of liquid
 - `q_ice` specific humidity of ice
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
    cv_m(param_set[, q::PhasePartition])

The isochoric specific heat capacity of moist air, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q` [`PhasePartition`](@ref). Without humidity argument, the results are for dry air.
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
    (R_m, cp_m, cv_m, γ_m) = gas_constants(param_set, q::PhasePartition)

Wrapper to compute the gas constant, specific heat capacities, and their 
ratio at once, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q` [`PhasePartition`](@ref)

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
    latent_heat_vapor(param_set, T)

The specific latent heat of vaporization, given

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
    latent_heat_sublim(param_set, T) 

The specific latent heat of sublimation, given

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
    latent_heat_fusion(param_set, T) 

The specific latent heat of fusion, given

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
    latent_heat_generic(param_set, T, LH_0, Δcp) 

The specific latent heat of a generic phase transition between two phases, 
computed using Kirchhoff's relation with constant isobaric specific heat 
capacities of the two phases, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `LH_0` latent heat at the reference temperature `T_0`
 - `Δcp` difference in isobaric specific heat capacities between the two phases
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
    latent_heat_mixed(param_set, T, λ)

Latent heat for a mixture of liquid and ice weighted by the 
liquid fraction, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` air temperature
 - `λ` liquid fraction
"""
@inline function latent_heat_mixed(
    param_set::APS,
    T::FT,
    λ::FT,
) where {FT <: Real}
    L_v = latent_heat_vapor(param_set, T)
    L_s = latent_heat_sublim(param_set, T)
    return λ * L_v + (1 - λ) * L_s
end

"""
    soundspeed_air(param_set, T[, q::PhasePartition])

The speed of sound in unstratified air, given

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
