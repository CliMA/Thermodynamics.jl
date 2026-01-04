# Formulas assuming dry (unsaturated) adiabatic (i.e. isentropic) processes

export DryAdiabaticProcess
export air_pressure_given_θ

"""
    DryAdiabaticProcess

A type for dispatching to isentropic (dry adiabatic) formulas.

"""
struct DryAdiabaticProcess end

"""
    air_pressure_given_θ(param_set, ::DryAdiabaticProcess, θ, Φ)

The air pressure for an isentropic process, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `θ` potential temperature
 - `Φ` gravitational potential

The pressure is computed using the hydrostatic balance and the definition of potential
temperature for an isentropic process: `p = p₀ * (1 - Φ/(θ * cp_d))^(1/κ_d)`, where `p₀`
is the reference pressure, `cp_d` is the isobaric specific heat capacity of dry air,
and `κ_d = R_d/cp_d`.
"""
@inline function air_pressure_given_θ(
    param_set::APS,
    ::DryAdiabaticProcess,
    θ,
    Φ,
)
    p0 = TP.p_ref_theta(param_set)
    R_d = TP.R_d(param_set)
    cp_d = TP.cp_d(param_set)
    return p0 * (1 - Φ / (θ * cp_d))^(cp_d / R_d)
end

"""
    air_pressure(param_set, ::DryAdiabaticProcess, T, T∞, p∞)

The air pressure for an isentropic process, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` current temperature
 - `T∞` reference temperature
 - `p∞` reference pressure

The pressure is computed using the isentropic relation: `p = p∞ * (T/T∞)^(1/κ_d)`,
where `κ_d = R_d/cp_d` is the ratio of the gas constant to the isobaric specific heat constant
of dry air.
"""
@inline function air_pressure(param_set::APS, ::DryAdiabaticProcess, T, T∞, p∞)
    kappa_d = TP.kappa_d(param_set)
    return p∞ * (T / T∞)^(1 / kappa_d)
end

"""
    air_temperature(param_set, ::DryAdiabaticProcess, p, θ)

The air temperature, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `θ` potential temperature

The temperature is computed using the definition of the dry potential temperature:
`T = θ * (p/p₀)^κ_d`, where `p₀` is the reference pressure and `κ_d = R_d/cp_d`.
"""
@inline function air_temperature(param_set::APS, ::DryAdiabaticProcess, p, θ)
    R_d = TP.R_d(param_set)
    cp_d = TP.cp_d(param_set)
    p0 = TP.p_ref_theta(param_set)
    return (p / p0)^(R_d / cp_d) * θ
end
