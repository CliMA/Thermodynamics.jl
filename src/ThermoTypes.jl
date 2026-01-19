export IndepVars
export ρe
export pe
export ph
export pρ
export pθ_li
export ρθ_li

export Phase
export Liquid
export Ice

"""
    IndepVars

Abstract type to dispatch over the independent variables used to compute
a thermodynamic quantity (e.g., in [`air_temperature`](@ref)).
"""
abstract type IndepVars end

"""
    ρe <: IndepVars

Density, internal energy, and specific humidity
"""
struct ρe <: IndepVars end

"""
    pe <: IndepVars

Pressure, internal energy, and specific humidity
"""
struct pe <: IndepVars end

"""
    ph <: IndepVars

Pressure, enthalpy, and specific humidity
"""
struct ph <: IndepVars end

"""
    pρ <: IndepVars

Pressure, density, and specific humidity
"""
struct pρ <: IndepVars end

"""
    pθ_li <: IndepVars

Pressure, liquid-ice potential temperature, and specific humidity
"""
struct pθ_li <: IndepVars end

"""
    ρθ_li <: IndepVars

Density, liquid-ice potential temperature, and specific humidity
"""
struct ρθ_li <: IndepVars end

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
