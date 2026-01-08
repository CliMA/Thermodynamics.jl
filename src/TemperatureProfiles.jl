module TemperatureProfiles

export TemperatureProfile,
    IsothermalProfile, DecayingTemperatureProfile, DryAdiabaticProfile

import ..Parameters
const TP = Parameters
const APS = TP.AbstractThermodynamicsParameters

"""
    TemperatureProfile

Abstract type for temperature or virtual temperature reference profiles 
that can be used in atmosphere models.

Instances of this type are required to be callable objects with the 
following signature

    T,p = (::TemperatureProfile)(param_set::APS, z::FT) where {FT}

where `T` is the temperature or virtual temperature (K), and `p` is 
the pressure (Pa).
"""
abstract type TemperatureProfile{FT} end

"""
    IsothermalProfile(param_set, T_virt)
    IsothermalProfile(param_set, ::Type{FT<:Real})

Uniform virtual temperature profile implemented as a special case 
of [`DecayingTemperatureProfile`](@ref).
"""
IsothermalProfile(param_set::APS, T_virt::FT) where {FT} =
    DecayingTemperatureProfile{FT}(param_set, T_virt, T_virt)

function IsothermalProfile(param_set::APS, ::Type{FT}) where {FT}
    T_virt = FT(TP.T_surf_ref(param_set))
    return DecayingTemperatureProfile{FT}(param_set, T_virt, T_virt)
end

"""
    DryAdiabaticProfile{FT} <: TemperatureProfile{FT}

Temperature profile with uniform potential temperature `θ` up to 
the height where a minimum temperature is reached.

# Fields

 - `T_surface`: Surface temperature (K)
 - `T_min_ref`: Minimum temperature (K)
"""
struct DryAdiabaticProfile{FT} <: TemperatureProfile{FT}
    "Surface temperature (K)"
    T_surface::FT
    "Minimum temperature (K)"
    T_min_ref::FT
    function DryAdiabaticProfile{FT}(
        param_set::APS,
        T_surface::FT = FT(TP.T_surf_ref(param_set)),
        T_min_ref::FT = FT(TP.T_min_ref(param_set)),
    ) where {FT}
        return new{FT}(T_surface, T_min_ref)
    end
end

"""
    (profile::DryAdiabaticProfile)(param_set::APS, z::FT) where {FT}

Returns dry adiabatic temperature and pressure profiles with zero relative humidity.
Temperature is truncated to be ≥ `profile.T_min_ref`.
"""
function (profile::DryAdiabaticProfile)(param_set::APS, z::FT) where {FT}
    R_d = TP.R_d(param_set)
    cp_d = TP.cp_d(param_set)
    grav = TP.grav(param_set)
    MSLP = TP.MSLP(param_set)

    # Temperature
    Γ = grav / cp_d
    T = max(profile.T_surface - Γ * z, profile.T_min_ref)

    # Pressure
    p = MSLP * (T / profile.T_surface)^(grav / (R_d * Γ))
    if T == profile.T_min_ref
        z_top = (profile.T_surface - profile.T_min_ref) / Γ
        H_min = R_d * profile.T_min_ref / grav
        p *= exp(-(z - z_top) / H_min)
    end
    return (T, p)
end

"""
    DecayingTemperatureProfile{FT} <: TemperatureProfile{FT}

Virtual temperature profile that decays smoothly with height `z` from `T_virt_surf` to `T_min_ref`
over height scale `H_t` (default: density scale height with `T_virt_surf`).

```math
T_{\\text{v}}(z) = \\max(T_{\\text{v, sfc}} − (T_{\\text{v, sfc}} - T_{\\text{v, min}}) \\tanh(z/H_{\\text{t}}), T_{\\text{v, min}})
```

# Fields

 - `T_virt_surf`: Virtual temperature at surface (K)
 - `T_min_ref`: Minimum virtual temperature at the top of the atmosphere (K)
 - `H_t`: Height scale over which virtual temperature drops (m)
"""
struct DecayingTemperatureProfile{FT} <: TemperatureProfile{FT}
    "Virtual temperature at surface (K)"
    T_virt_surf::FT
    "Minimum virtual temperature at the top of the atmosphere (K)"
    T_min_ref::FT
    "Height scale over which virtual temperature drops (m)"
    H_t::FT
    function DecayingTemperatureProfile{FT}(
        param_set::APS,
        T_virt_surf::FT = FT(TP.T_surf_ref(param_set)),
        T_min_ref::FT = FT(TP.T_min_ref(param_set)),
        H_t::FT = FT(TP.R_d(param_set)) * T_virt_surf / FT(TP.grav(param_set)),
    ) where {FT}
        return new{FT}(T_virt_surf, T_min_ref, H_t)
    end
end

"""
    (profile::DecayingTemperatureProfile)(param_set::APS, z::FT) where {FT}

Returns decaying temperature and pressure profiles using hyperbolic tangent decay.
"""
function (profile::DecayingTemperatureProfile)(param_set::APS, z::FT) where {FT}
    R_d = TP.R_d(param_set)
    grav = TP.grav(param_set)
    MSLP = TP.MSLP(param_set)

    # Scale height for surface temperature
    H_sfc = R_d * profile.T_virt_surf / grav
    H_t = profile.H_t
    z′ = z / H_t
    tanh_z′ = tanh(z′)

    ΔTv = profile.T_virt_surf - profile.T_min_ref
    Tv = profile.T_virt_surf - ΔTv * tanh_z′

    ΔTv′ = ΔTv / profile.T_virt_surf
    p = -H_t * (z′ + ΔTv′ * (log(1 - ΔTv′ * tanh_z′) - log(1 + tanh_z′) + z′))
    p /= H_sfc * (1 - ΔTv′^2)
    p = MSLP * exp(p)
    return (Tv, p)
end

end
