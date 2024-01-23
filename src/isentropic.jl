#####
##### Formulas assuming dry (unsaturated) adiabatic (i.e. isentropic) processes
#####

export DryAdiabaticProcess

export air_pressure, air_temperature

"""
    DryAdiabaticProcess

For dispatching to isentropic formulas
"""
struct DryAdiabaticProcess end

"""
    air_pressure_given_θ(param_set, θ::FT, Φ::FT, ::DryAdiabaticProcess)

The air pressure for an isentropic process, where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `θ` potential temperature
 - `Φ` gravitational potential
"""
@inline function air_pressure_given_θ(
    param_set::APS,
    θ::FT,
    Φ::FT,
    ::DryAdiabaticProcess,
) where {FT <: AbstractFloat}
    p0::FT = TP.p_ref_theta(param_set)
    _R_d::FT = TP.R_d(param_set)
    _cp_d::FT = TP.cp_d(param_set)
    return p0 * (1 - Φ / (θ * _cp_d))^(_cp_d / _R_d)
end

"""
    air_pressure(param_set, T::FT, T∞::FT, p∞::FT, ::DryAdiabaticProcess)

The air pressure for an isentropic process, where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `T∞` ambient temperature
 - `p∞` ambient pressure
"""
@inline function air_pressure(
    param_set::APS,
    T::FT,
    T∞::FT,
    p∞::FT,
    ::DryAdiabaticProcess,
) where {FT <: AbstractFloat}
    _kappa_d::FT = TP.kappa_d(param_set)
    return p∞ * (T / T∞)^(FT(1) / _kappa_d)
end

"""
    air_temperature(param_set, p::FT, θ::FT, Φ::FT, ::DryAdiabaticProcess)

The air temperature for an isentropic process, where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `θ` potential temperature
"""
@inline function air_temperature(
    param_set::APS,
    p::FT,
    θ::FT,
    ::DryAdiabaticProcess,
) where {FT <: AbstractFloat}
    _R_d::FT = TP.R_d(param_set)
    _cp_d::FT = TP.cp_d(param_set)
    p0::FT = TP.p_ref_theta(param_set)
    return (p / p0)^(_R_d / _cp_d) * θ
end
