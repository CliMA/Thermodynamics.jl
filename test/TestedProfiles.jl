"""
    TestedProfiles

This module contains functions to compute thermodynamic test profiles
for various atmospheric conditions.

Notes:
- These profiles are intended for **testing**; some are in thermodynamic equilibrium
  (`EquilMoistProfiles`) while others are intentionally non-equilibrium (`NonEquilMoistProfiles`).
- `TemperatureProfiles` return a temperature-like quantity which may be interpreted as
  virtual temperature, depending on the profile implementation.
"""
module TestedProfiles

import Thermodynamics as TD
import Thermodynamics.Parameters as TP
const APS = TP.ThermodynamicsParameters
const TDTP = TD.TemperatureProfiles

import Random

"""
    ProfileSet

A set of profiles used to test Thermodynamics.
"""
struct ProfileSet{AFT}
    z::AFT          # Altitude
    T::AFT          # Temperature-like quantity (T or T_virt)
    p::AFT          # Pressure
    RS::AFT         # Relative saturation factor used to scale q_vap_saturation
    e_int::AFT      # Internal energy
    h::AFT          # Specific enthalpy
    ρ::AFT          # Density
    θ_li::AFT  # Liquid ice potential temperature
    q_tot::AFT      # Total specific humidity
    q_liq::AFT      # Liquid specific humidity
    q_ice::AFT      # Ice specific humidity
    RH::AFT         # Relative humidity
    e_pot::AFT      # Gravitational potential
    u::AFT          # Velocity (x component)
    v::AFT          # Velocity (y component)
    w::AFT          # Velocity (z component)
    e_kin::AFT      # Kinetic energy
end

function Base.iterate(ps::ProfileSet, state = 1)
    state > length(ps.z) && return nothing
    nt = (
        z = ps.z[state],
        T = ps.T[state],
        p = ps.p[state],
        RS = ps.RS[state],
        e_int = ps.e_int[state],
        h = ps.h[state],
        ρ = ps.ρ[state],
        θ_li = ps.θ_li[state],
        q_tot = ps.q_tot[state],
        q_liq = ps.q_liq[state],
        q_ice = ps.q_ice[state],
        RH = ps.RH[state],
        e_pot = ps.e_pot[state],
        u = ps.u[state],
        v = ps.v[state],
        w = ps.w[state],
        e_kin = ps.e_kin[state],
    )
    return (nt, state + 1)
end
Base.IteratorSize(::ProfileSet) = Base.HasLength()
Base.length(ps::ProfileSet) = length(ps.z)

"""
    input_config(
        ArrayType;
        n=50,
        n_RS1=10,
        n_RS2=20,
        T_min=150,
        T_surface=340
    )

Return input arguments to construct profiles.

`relative_sat` is a dimensionless scaling applied to the saturation vapor specific humidity
(`q_vap_saturation`) to generate both unsaturated and slightly supersaturated cases (values > 1).
"""
function input_config(
    ArrayType;
    n = 50,
    n_RS1 = 10,
    n_RS2 = 20,
    T_surface = 340,
    T_min = 150,
)
    n_RS = n_RS1 + n_RS2
    z_range = ArrayType(range(0, stop = 2.5e4, length = n))
    relative_sat1 = ArrayType(range(0, stop = 1, length = n_RS1))
    relative_sat2 = ArrayType(range(1, stop = 1.02, length = n_RS2))
    relative_sat = vcat(relative_sat1, relative_sat2)
    return z_range, relative_sat, T_surface, T_min
end

"""
    shared_profiles(
        param_set::APS,
        z_range::AbstractArray,
        relative_sat::AbstractArray,
        T_surface,
        T_min,
    )

Compute profiles shared across dry and moist thermodynamic states:
 - `z` altitude
 - `T` temperature-like quantity returned by `TemperatureProfiles` (often interpreted as virtual temperature)
 - `p` pressure
 - `RS` relative saturation factor
"""
function shared_profiles(
    param_set::APS,
    z_range::AbstractArray,
    relative_sat::AbstractArray,
    T_surface,
    T_min,
)
    FT = eltype(z_range)
    n_RS = length(relative_sat)
    n = length(z_range)
    T = similar(z_range, n * n_RS)
    p = similar(z_range, n * n_RS)
    RS = similar(z_range, n * n_RS)
    z = similar(z_range, n * n_RS)
    linear_indices = LinearIndices((1:n, 1:n_RS))
    profile =
        TDTP.DecayingTemperatureProfile{FT}(param_set, FT(T_surface), FT(T_min))
    for i in linear_indices.indices[1]
        for j in linear_indices.indices[2]
            k = linear_indices[i, j]
            z[k] = z_range[i]
            T[k], p[k] = profile(param_set, z[k])
            RS[k] = relative_sat[j]
        end
    end
    return z, T, p, RS
end

"""
    DryProfiles(param_set, ::Type{ArrayType})

Returns a `ProfileSet` for dry atmospheric conditions (no moisture).
"""
function DryProfiles(param_set::APS, ::Type{ArrayType}) where {ArrayType}
    z_range, relative_sat, T_surface, T_min = input_config(ArrayType)
    z, T, p, RS =
        shared_profiles(param_set, z_range, relative_sat, T_surface, T_min)

    FT = eltype(T)
    R_d = TP.R_d(param_set)
    grav = TP.grav(param_set)
    ρ = p ./ (R_d .* T)

    # Dry state: no moisture
    q_tot = zero(T)
    q_liq = zero(T)
    q_ice = zero(T)

    e_int = TD.internal_energy.(Ref(param_set), T, q_tot, q_liq, q_ice)
    h = TD.enthalpy.(Ref(param_set), T, q_tot, q_liq, q_ice)
    θ_li =
        TD.liquid_ice_pottemp.(Ref(param_set), T, ρ, q_tot, q_liq, q_ice)
    RH = TD.relative_humidity.(Ref(param_set), T, p, q_tot, q_liq, q_ice)

    e_pot = grav * z
    # Use a local RNG for reproducibility without mutating global RNG state.
    rng = Random.MersenneTwister(15)
    u = rand(rng, FT, size(T)) * 50
    v = rand(rng, FT, size(T)) * 50
    w = rand(rng, FT, size(T)) * 50
    e_kin = (u .^ 2 .+ v .^ 2 .+ w .^ 2) / 2

    return ProfileSet{typeof(T)}(
        z,
        T,
        p,
        RS,
        e_int,
        h,
        ρ,
        θ_li,
        q_tot,
        q_liq,
        q_ice,
        RH,
        e_pot,
        u,
        v,
        w,
        e_kin,
    )
end

"""
    EquilMoistProfiles(param_set, ::Type{ArrayType})

Returns a `ProfileSet` for moist atmospheric conditions in thermodynamic equilibrium.

`RS` is interpreted as a scaling factor for the equilibrium saturation vapor specific humidity,
i.e. `q_tot = RS * q_vap_saturation(param_set, T, ρ)` (and any excess condenses).
"""
function EquilMoistProfiles(param_set::APS, ::Type{ArrayType}) where {ArrayType}
    z_range, relative_sat, T_surface, T_min = input_config(ArrayType)
    z, T_virt, p, RS =
        shared_profiles(param_set, z_range, relative_sat, T_surface, T_min)

    FT = eltype(T_virt)
    R_d = TP.R_d(param_set)
    grav = TP.grav(param_set)
    ρ = p ./ (R_d .* T_virt)

    # Compute equilibrium moisture.
    # Here we treat the profile's returned temperature-like quantity as `T`.
    T = T_virt
    q_tot = RS .* TD.q_vap_saturation.(Ref(param_set), T, ρ)

    # Compute phase partitioning for each element
    q_liq = similar(T)
    q_ice = similar(T)
    for i in eachindex(T)
        q_liq[i], q_ice[i] =
            TD.condensate_partition(param_set, T[i], ρ[i], q_tot[i])
    end

    # Update pressure to be thermodynamically consistent with moist air equation of state at fixed ρ.
    p = TD.air_pressure.(Ref(param_set), T, ρ, q_tot, q_liq, q_ice)

    e_int = TD.internal_energy.(Ref(param_set), T, q_tot, q_liq, q_ice)
    h = TD.enthalpy.(Ref(param_set), T, q_tot, q_liq, q_ice)
    θ_li =
        TD.liquid_ice_pottemp.(Ref(param_set), T, ρ, q_tot, q_liq, q_ice)
    RH = TD.relative_humidity.(Ref(param_set), T, p, q_tot, q_liq, q_ice)

    e_pot = grav * z
    # Use a local RNG for reproducibility without mutating global RNG state.
    rng = Random.MersenneTwister(15)
    u = rand(rng, FT, size(T)) * 50
    v = rand(rng, FT, size(T)) * 50
    w = rand(rng, FT, size(T)) * 50
    e_kin = (u .^ 2 .+ v .^ 2 .+ w .^ 2) / 2

    return ProfileSet{typeof(T)}(
        z,
        T,
        p,
        RS,
        e_int,
        h,
        ρ,
        θ_li,
        q_tot,
        q_liq,
        q_ice,
        RH,
        e_pot,
        u,
        v,
        w,
        e_kin,
    )
end

"""
    NonEquilMoistProfiles(param_set, ::Type{ArrayType})

Returns a `ProfileSet` for moist atmospheric conditions 
in thermodynamic non-equilibrium (supersaturated).

This constructor intentionally does **not** enforce thermodynamic equilibrium:
it prescribes `q_tot` and then assigns `q_liq`/`q_ice` via a fixed split for simplicity.
"""
function NonEquilMoistProfiles(
    param_set::APS,
    ::Type{ArrayType},
) where {ArrayType}
    z_range, relative_sat, T_surface, T_min = input_config(ArrayType)
    z, T, p, RS =
        shared_profiles(param_set, z_range, relative_sat, T_surface, T_min)

    FT = eltype(T)
    R_d = TP.R_d(param_set)
    grav = TP.grav(param_set)

    # Non-equilibrium: prescribe q_tot > q_sat for some profiles
    q_vap_sat = TD.q_vap_saturation.(Ref(param_set), T, p ./ (R_d .* T))
    q_tot = RS .* q_vap_sat

    # Prescribe non-equilibrium phase partitioning
    # For simplicity, prescribe q_liq and q_ice separately
    q_liq = max.(zero(FT), (q_tot .- q_vap_sat) .* FT(0.6))  # 60% liquid
    q_ice = max.(zero(FT), (q_tot .- q_vap_sat) .* FT(0.4))  # 40% ice

    ρ = TD.air_density.(Ref(param_set), T, p, q_tot, q_liq, q_ice)

    e_int = TD.internal_energy.(Ref(param_set), T, q_tot, q_liq, q_ice)
    h = TD.enthalpy.(Ref(param_set), T, q_tot, q_liq, q_ice)
    θ_li =
        TD.liquid_ice_pottemp.(Ref(param_set), T, ρ, q_tot, q_liq, q_ice)
    RH = TD.relative_humidity.(Ref(param_set), T, p, q_tot, q_liq, q_ice)

    e_pot = grav * z
    # Use a local RNG for reproducibility without mutating global RNG state.
    rng = Random.MersenneTwister(15)
    u = rand(rng, FT, size(T)) * 50
    v = rand(rng, FT, size(T)) * 50
    w = rand(rng, FT, size(T)) * 50
    e_kin = (u .^ 2 .+ v .^ 2 .+ w .^ 2) / 2

    return ProfileSet{typeof(T)}(
        z,
        T,
        p,
        RS,
        e_int,
        h,
        ρ,
        θ_li,
        q_tot,
        q_liq,
        q_ice,
        RH,
        e_pot,
        u,
        v,
        w,
        e_kin,
    )
end

# Backwards compatibility aliases (deprecatedto be removed)
const PhaseDryProfiles = DryProfiles
const PhaseEquilProfiles = EquilMoistProfiles

end # module
