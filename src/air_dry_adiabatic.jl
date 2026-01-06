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

The air pressure for an isentropic process.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `θ`: potential temperature [K]
 - `Φ`: gravitational potential [J/kg]

# Returns
 - `p`: air pressure [Pa]

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

The air pressure for an isentropic process.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: current temperature [K]
 - `T∞`: reference temperature [K]
 - `p∞`: reference pressure [Pa]

# Returns
 - `p`: air pressure [Pa]

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

The air temperature for an isentropic process.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `p`: pressure [Pa]
 - `θ`: potential temperature [K]

# Returns
 - `T`: air temperature [K]

The temperature is computed using the definition of the dry potential temperature:
`T = θ * (p/p₀)^κ_d`, where `p₀` is the reference pressure and `κ_d = R_d/cp_d`.
"""
@inline function air_temperature(param_set::APS, ::DryAdiabaticProcess, p, θ)
    R_d = TP.R_d(param_set)
    cp_d = TP.cp_d(param_set)
    p0 = TP.p_ref_theta(param_set)
    return (p / p0)^(R_d / cp_d) * θ
end
