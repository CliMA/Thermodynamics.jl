#####
##### Formulas assuming dry (unsaturated) adiabatic (i.e. isentropic) processes
#####

export DryAdiabaticProcess

export air_pressure_given_θ, air_pressure, air_temperature

"""
    DryAdiabaticProcess

For dispatching to isentropic formulas
"""
struct DryAdiabaticProcess end

"""
    air_pressure_given_θ(θ::FT, Φ::FT, ::DryAdiabaticProcess)

The air pressure for an isentropic process, where

 - `θ` potential temperature
 - `Φ` gravitational potential
"""
function air_pressure_given_θ(
    param_set::PS,
    θ::FT,
    Φ::FT,
    ::DryAdiabaticProcess,
) where {FT, PS}
    return FT(MSLP(param_set)) * (1 - Φ / (θ * FT(cp_d(param_set))))^(FT(cp_d(param_set)) / FT(R_d(param_set)))
end

"""
    air_pressure(T::FT, T∞::FT, p∞::FT, ::DryAdiabaticProcess)

The air pressure for an isentropic process, where

 - `T` temperature
 - `T∞` ambient temperature
 - `p∞` ambient pressure
"""
function air_pressure(
    param_set::PS,
    T::FT,
    T∞::FT,
    p∞::FT,
    ::DryAdiabaticProcess,
) where {FT, PS}
    return p∞ * (T / T∞)^(FT(1) / FT(kappa_d(param_set)))
end

"""
    air_temperature(p::FT, θ::FT, Φ::FT, ::DryAdiabaticProcess)

The air temperature for an isentropic process, where

 - `p` pressure
 - `θ` potential temperature
"""
function air_temperature(
    param_set::PS,
    p::FT,
    θ::FT,
    ::DryAdiabaticProcess,
) where {FT, PS}
    return (p / FT(MSLP(param_set)))^(FT(R_d(param_set)) / FT(cp_d(param_set))) * θ
end
