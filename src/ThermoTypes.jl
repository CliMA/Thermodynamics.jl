export IndepVars
export ρeq
export peq
export phq
export pρq
export pθ_liq_ice_q
export ρθ_liq_ice_q

export Phase
export Liquid
export Ice

"""
    IndepVars

Abstract type to dispatch over the independent variables used to compute
a thermodynamic quantity (e.g., in [`air_temperature](@ref)).
"""
abstract type IndepVars end

"""
    ρeq <: IndepVars

Density, internal energy, and total specific humidity.
"""
struct ρeq <: IndepVars end

"""
    peq <: IndepVars

Pressure, internal energy, and total specific humidity.
"""
struct peq <: IndepVars end

"""
    phq <: IndepVars

Pressure, enthalpy, and total specific humidity.
"""
struct phq <: IndepVars end

"""
    pρq <: IndepVars

Pressure, density, and total specific humidity.
"""
struct pρq <: IndepVars end

"""
    pθ_liq_ice_q <: IndepVars

Pressure, liquid-ice potential temperature, and total specific humidity.
"""
struct pθ_liq_ice_q <: IndepVars end

"""
    ρθ_liq_ice_q <: IndepVars

Density, liquid-ice potential temperature, and total specific humidity.
"""
struct ρθ_liq_ice_q <: IndepVars end

"""
    Phase

A condensed phase, to dispatch over (e.g., in
[`saturation_vapor_pressure`](@ref)).
"""
abstract type Phase end

"""
    Liquid <: Phase

A liquid phase, to dispatch over (e.g., in
[`saturation_vapor_pressure`](@ref)).
"""
struct Liquid <: Phase end

"""
    Ice <: Phase

An ice phase, to dispatch over (e.g., in
[`saturation_vapor_pressure`](@ref)).
"""
struct Ice <: Phase end
