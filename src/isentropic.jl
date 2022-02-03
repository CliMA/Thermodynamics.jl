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
function air_pressure_given_θ(
    param_set::ThermodynamicsParameters,
    θ::FT,
    Φ::FT,
    ::DryAdiabaticProcess,
) where {FT <: AbstractFloat}
    MSLP::FT = param_set.MSLP
    _R_d::FT = param_set.R_d
    _cp_d::FT = param_set.cp_d
    return MSLP * (1 - Φ / (θ * _cp_d))^(_cp_d / _R_d)
end

"""
    air_pressure(param_set, T::FT, T∞::FT, p∞::FT, ::DryAdiabaticProcess)

The air pressure for an isentropic process, where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `T∞` ambient temperature
 - `p∞` ambient pressure
"""
function air_pressure(
    param_set::ThermodynamicsParameters,
    T::FT,
    T∞::FT,
    p∞::FT,
    ::DryAdiabaticProcess,
) where {FT <: AbstractFloat}
    _kappa_d::FT = param_set.kappa_d
    return p∞ * (T / T∞)^(FT(1) / _kappa_d)
end

"""
    air_temperature(param_set, p::FT, θ::FT, Φ::FT, ::DryAdiabaticProcess)

The air temperature for an isentropic process, where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `θ` potential temperature
"""
function air_temperature(
    param_set::ThermodynamicsParameters,
    p::FT,
    θ::FT,
    ::DryAdiabaticProcess,
) where {FT <: AbstractFloat}
    _R_d::FT = param_set.R_d
    _cp_d::FT = param_set.cp_d
    MSLP::FT = param_set.MSLP
    return (p / MSLP)^(_R_d / _cp_d) * θ
end
