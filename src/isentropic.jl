#####
##### Formulas assuming dry (unsaturated) adiabatic (i.e. isentropic) processes
#####

export DryAdiabaticProcess

export air_pressure, air_temperature

"""
    DryAdiabaticProcess

A type for dispatching to isentropic (dry adiabatic) formulas.

"""
struct DryAdiabaticProcess end

"""
    air_pressure_given_θ(param_set, θ, Φ, ::DryAdiabaticProcess)

The air pressure for an isentropic process, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `θ` potential temperature
 - `Φ` gravitational potential

The pressure is computed using the hydrostatic balance and the definition of potential temperature
for an isentropic process: `p = p₀ * (1 - Φ/(θ * cₚ))^(cₚ/R)`, where `p₀` is the reference pressure,
`cₚ` is the isobaric specific heat capacity of dry air, and `R` is the gas constant of dry air.
"""
@inline function air_pressure_given_θ(
    param_set::APS,
    θ::FT,
    Φ::FT,
    ::DryAdiabaticProcess,
) where {FT <: Real}
    p0 = TP.p_ref_theta(param_set)
    _R_d = TP.R_d(param_set)
    _cp_d = TP.cp_d(param_set)
    return p0 * (1 - Φ / (θ * _cp_d))^(_cp_d / _R_d)
end

"""
    air_pressure(param_set, T, T∞, p∞, ::DryAdiabaticProcess)

The air pressure for an isentropic process, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `T∞` ambient temperature
 - `p∞` ambient pressure

The pressure is computed using the isentropic relation: `p = p∞ * (T/T∞)^(1/κ)`,
where `κ = R/cₚ` is the ratio of the gas constant to the isobaric specific heat capacity
of dry air.
"""
@inline function air_pressure(
    param_set::APS,
    T::FT,
    T∞::FT,
    p∞::FT,
    ::DryAdiabaticProcess,
) where {FT <: Real}
    _kappa_d = TP.kappa_d(param_set)
    return p∞ * (T / T∞)^(FT(1) / _kappa_d)
end

"""
    air_temperature(param_set, p, θ, ::DryAdiabaticProcess)

The air temperature, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `θ` potential temperature

The temperature is computed using the definition of the dry potential temperature:
`T = θ * (p/p₀)^(R/cₚ)`, where `p₀` is the reference pressure, `R` is the gas constant of dry air,
and `cₚ` is the isobaric specific heat capacity of dry air.
"""
@inline function air_temperature(
    param_set::APS,
    p::FT,
    θ::FT,
    ::DryAdiabaticProcess,
) where {FT <: Real}
    _R_d = TP.R_d(param_set)
    _cp_d = TP.cp_d(param_set)
    p0 = TP.p_ref_theta(param_set)
    return (p / p0)^(_R_d / _cp_d) * θ
end
