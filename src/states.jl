export PhasePartition
# Thermodynamic states
export ThermodynamicState,
    PhaseDry,
    PhaseDry_ρT,
    PhaseDry_pT,
    PhaseDry_pe,
    PhaseDry_ρθ,
    PhaseDry_pθ,
    PhaseDry_ρp,
    PhaseEquil,
    PhaseEquil_ρeq,
    PhaseEquil_ρTq,
    PhaseEquil_pTq,
    PhaseEquil_peq,
    PhaseEquil_ρθq,
    PhaseEquil_pθq,
    PhaseEquil_ρpq,
    PhaseNonEquil,
    PhaseNonEquil_ρTq,
    PhaseNonEquil_pTq,
    PhaseNonEquil_ρθq,
    PhaseNonEquil_pθq,
    PhaseNonEquil_peq,
    PhaseNonEquil_ρpq

"""
    ThermodynamicState{FT}

A thermodynamic state, which can be initialized for
various thermodynamic formulations (via its sub-types).
All `ThermodynamicState`'s have access to functions to
compute all other thermodynamic properties.
"""
abstract type ThermodynamicState{FT} end

abstract type AbstractPhaseDry{FT} <: ThermodynamicState{FT} end
abstract type AbstractPhaseEquil{FT} <: ThermodynamicState{FT} end
abstract type AbstractPhaseNonEquil{FT} <: ThermodynamicState{FT} end

Base.eltype(::ThermodynamicState{FT}) where {FT} = FT

"""
    PhasePartition

Represents the mass fractions of the moist air mixture.

# Constructors

    PhasePartition(q_tot::Real[, q_liq::Real[, q_ice::Real]])
    PhasePartition(param_set::APS, ts::ThermodynamicState)

See also [`PhasePartition_equil`](@ref)

# Fields

$(DocStringExtensions.FIELDS)
"""
struct PhasePartition{FT <: Real}
    "total specific humidity"
    tot::FT
    "liquid water specific humidity (default: `0`)"
    liq::FT
    "ice specific humidity (default: `0`)"
    ice::FT
    function PhasePartition(tot::FT, liq::FT, ice::FT) where {FT}
        q_tot_safe = max(tot, 0)
        q_liq_safe = max(liq, 0)
        q_ice_safe = max(ice, 0)
        return new{FT}(q_tot_safe, q_liq_safe, q_ice_safe)
    end
end

PhasePartition(q_tot::FT, q_liq::FT) where {FT <: Real} =
    PhasePartition(q_tot, q_liq, zero(FT))
PhasePartition(q_tot::FT) where {FT <: Real} =
    PhasePartition(q_tot, zero(FT), zero(FT))

const ITERTYPE = Union{Int, Nothing}
TOLTYPE(FT) = Union{FT, Nothing}

#####
##### Dry states
#####

"""
    PhaseDry{FT} <: AbstractPhaseDry

A dry thermodynamic state (`q_tot = 0`).

# Constructors

    PhaseDry(param_set, e_int, ρ)

# Fields

$(DocStringExtensions.FIELDS)
"""
struct PhaseDry{FT} <: AbstractPhaseDry{FT}
    "internal energy"
    e_int::FT
    "density of dry air"
    ρ::FT
end
PhaseDry(param_set::APS, e_int::FT, ρ::FT) where {FT} = PhaseDry{FT}(e_int, ρ)

"""
    PhaseDry_pT(param_set, p, T)

Constructs a [`PhaseDry`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
"""
function PhaseDry_pT(param_set::APS, p::FT, T::FT) where {FT <: Real}
    e_int = internal_energy(param_set, T)
    ρ = air_density(param_set, T, p)
    return PhaseDry{FT}(e_int, ρ)
end

"""
    PhaseDry_pe(param_set, p, e_int)

Constructs a [`PhaseDry`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `e_int` internal energy
"""
function PhaseDry_pe(param_set::APS, p::FT, e_int::FT) where {FT <: Real}
    T = air_temperature(param_set, e_int)
    ρ = air_density(param_set, T, p)
    return PhaseDry{FT}(e_int, ρ)
end

"""
    PhaseDry_ρθ(param_set, ρ, θ_dry)

Constructs a [`PhaseDry`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `θ_dry` dry potential temperature
"""
function PhaseDry_ρθ(param_set::APS, ρ::FT, θ_dry::FT) where {FT <: Real}
    T = air_temperature_given_ρθq(param_set, ρ, θ_dry)
    e_int = internal_energy(param_set, T)
    return PhaseDry{FT}(e_int, ρ)
end

"""
    PhaseDry_pθ(param_set, p, θ_dry)

Constructs a [`PhaseDry`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `θ_dry` dry potential temperature
"""
function PhaseDry_pθ(param_set::APS, p::FT, θ_dry::FT) where {FT <: Real}
    T = exner_given_pressure(param_set, p) * θ_dry
    e_int = internal_energy(param_set, T)
    ρ = air_density(param_set, T, p)
    return PhaseDry{FT}(e_int, ρ)
end

"""
    PhaseDry_ρT(param_set, ρ, T)

Constructs a [`PhaseDry`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `T` temperature
"""
function PhaseDry_ρT(param_set::APS, ρ::FT, T::FT) where {FT <: Real}
    e_int = internal_energy(param_set, T)
    return PhaseDry{FT}(e_int, ρ)
end

"""
    PhaseDry_ρp(param_set, ρ, p)

Constructs a [`PhaseDry`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `p` pressure
"""
function PhaseDry_ρp(param_set::APS, ρ::FT, p::FT) where {FT <: Real}
    T = air_temperature_from_ideal_gas_law(param_set, p, ρ)
    e_int = internal_energy(param_set, T)
    return PhaseDry{FT}(e_int, ρ)
end

#####
##### Equilibrium states
#####

"""
    PhaseEquil{FT} <: AbstractPhaseEquil

A thermodynamic state assuming thermodynamic equilibrium (therefore, saturation adjustment
may be needed).

# Constructors

    PhaseEquil(param_set, ρ, e_int, q_tot)

# Fields

$(DocStringExtensions.FIELDS)
"""
struct PhaseEquil{FT} <: AbstractPhaseEquil{FT}
    "density of air (potentially moist)"
    ρ::FT
    "air pressure"
    p::FT
    "internal energy"
    e_int::FT
    "total specific humidity"
    q_tot::FT
    "temperature: computed via [`saturation_adjustment`](@ref)"
    T::FT
end

"""
    PhaseEquil_ρeq

Moist thermodynamic phase, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `e_int` internal energy
 - `q_tot` total specific humidity
and, optionally
 - `maxiter` maximum iterations for saturation adjustment
 - `temperature_tol` temperature tolerance for saturation adjustment
 - `sat_adjust_method` the numerical method to use.
    See the [`Thermodynamics`](@ref) for options.
"""
function PhaseEquil_ρeq(
    param_set::APS,
    ρ::FT,
    e_int::FT,
    q_tot::FT,
    maxiter::IT = nothing,
    temperature_tol::FTT = nothing,
    ::Type{sat_adjust_method} = RS.NewtonsMethod,
) where {FT <: Real, sat_adjust_method, IT <: ITERTYPE, FTT <: TOLTYPE(FT)}
    maxiter === nothing && (maxiter = 8)
    temperature_tol === nothing && (temperature_tol = FT(1e-1))
    phase_type = PhaseEquil{FT}
    q_tot_safe = clamp(q_tot, FT(0), FT(1))
    T = saturation_adjustment(
        sat_adjust_method,
        param_set,
        e_int,
        ρ,
        q_tot_safe,
        phase_type,
        maxiter,
        temperature_tol,
    )
    q_pt = PhasePartition_equil(param_set, T, ρ, q_tot_safe, phase_type)
    p = air_pressure(param_set, T, ρ, q_pt)
    return PhaseEquil{FT}(ρ, p, e_int, q_tot_safe, T)
end

# Convenience method for comparing Numerical
# methods without having to specify maxiter
# and temperature_tol. maxiter and temperature_tol
# should be in sync with the PhaseEquil(...) constructor
function PhaseEquil_dev_only(
    param_set::APS,
    ρ::FT,
    e_int::FT,
    q_tot::FT;
    maxiter::Int = 8,
    temperature_tol::FT = FT(1e-1),
    sat_adjust_method::Type{NM} = RS.NewtonsMethod,
) where {FT <: Real, NM}
    return PhaseEquil_ρeq(
        param_set,
        ρ,
        e_int,
        q_tot,
        maxiter,
        temperature_tol,
        sat_adjust_method,
    )
end

"""
    PhaseEquil_ρθq(param_set, ρ, θ_liq_ice, q_tot)

Constructs a [`PhaseEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_tot` total specific humidity
 - `temperature_tol` temperature tolerance for saturation adjustment
 - `maxiter` maximum iterations for saturation adjustment
"""
function PhaseEquil_ρθq(
    param_set::APS,
    ρ::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    maxiter::IT = nothing,
    temperature_tol::FTT = nothing,
) where {FT <: Real, IT <: ITERTYPE, FTT <: TOLTYPE(FT)}
    maxiter === nothing && (maxiter = 36)
    temperature_tol === nothing && (temperature_tol = FT(1e-1))
    phase_type = PhaseEquil{FT}
    tol = RS.ResidualTolerance(temperature_tol)
    T = saturation_adjustment_given_ρθq(
        param_set,
        ρ,
        θ_liq_ice,
        q_tot,
        phase_type,
        maxiter,
        tol,
    )
    q_pt = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    p = air_pressure(param_set, T, ρ, q_pt)
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseEquil{FT}(ρ, p, e_int, q_tot, T)
end

"""
    PhaseEquil_ρTq(param_set, ρ, T, q_tot)

Constructs a [`PhaseEquil`](@ref) thermodynamic state from temperature.

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `T` temperature
 - `q_tot` total specific humidity
"""
function PhaseEquil_ρTq(
    param_set::APS,
    ρ::FT,
    T::FT,
    q_tot::FT,
) where {FT <: Real}
    phase_type = PhaseEquil{FT}
    q_pt = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    p = air_pressure(param_set, T, ρ, q_pt)
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseEquil{FT}(ρ, p, e_int, q_tot, T)
end

"""
    PhaseEquil_pTq(param_set, p, T, q_tot)

Constructs a [`PhaseEquil`](@ref) thermodynamic state from temperature.

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `q_tot` total specific humidity
"""
function PhaseEquil_pTq(
    param_set::APS,
    p::FT,
    T::FT,
    q_tot::FT,
) where {FT <: Real}
    phase_type = PhaseEquil{FT}
    q_tot_safe = clamp(q_tot, FT(0), FT(1))
    q_pt = PhasePartition_equil_given_p(param_set, T, p, q_tot_safe, phase_type)
    ρ = air_density(param_set, T, p, q_pt)
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseEquil{FT}(ρ, p, e_int, q_tot_safe, T)
end

"""
    PhaseEquil_peq(param_set, p, e_int, q_tot)

Constructs a [`PhaseEquil`](@ref) thermodynamic state from temperature.

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `e_int` temperature
 - `q_tot` total specific humidity
"""
function PhaseEquil_peq(
    param_set::APS,
    p::FT,
    e_int::FT,
    q_tot::FT,
    maxiter::IT = nothing,
    temperature_tol::FTT = nothing,
    ::Type{sat_adjust_method} = RS.SecantMethod,
) where {FT <: Real, sat_adjust_method, IT <: ITERTYPE, FTT <: TOLTYPE(FT)}
    maxiter === nothing && (maxiter = 40)
    temperature_tol === nothing && (temperature_tol = FT(1e-2))
    phase_type = PhaseEquil{FT}
    q_tot_safe = clamp(q_tot, FT(0), FT(1))
    T = saturation_adjustment_given_peq(
        sat_adjust_method,
        param_set,
        p,
        e_int,
        q_tot_safe,
        phase_type,
        maxiter,
        temperature_tol,
    )
    q_pt = PhasePartition_equil_given_p(param_set, T, p, q_tot_safe, phase_type)
    ρ = air_density(param_set, T, p, q_pt)
    return PhaseEquil{FT}(ρ, p, e_int, q_tot_safe, T)
end

"""
    PhaseEquil_ρpq(param_set, ρ, p, q_tot, perform_sat_adjust=true)

Constructs a [`PhaseEquil`](@ref) thermodynamic state from temperature.

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `p` pressure
 - `q_tot` total specific humidity
 - `perform_sat_adjust` Bool indicating to perform saturation adjustment
 - `maxiter` maximum number of iterations to perform in saturation adjustment
 - `sat_adjust_method` the numerical method to use.
    See the [`Thermodynamics`](@ref) for options.
"""
function PhaseEquil_ρpq(
    param_set::APS,
    ρ::FT,
    p::FT,
    q_tot::FT,
    perform_sat_adjust = false,
    maxiter::Int = 5,
    ::Type{sat_adjust_method} = RS.NewtonsMethodAD,
) where {FT <: Real, sat_adjust_method}

    phase_type = PhaseEquil{FT}
    if perform_sat_adjust
        T = saturation_adjustment_ρpq(
            sat_adjust_method,
            param_set,
            ρ,
            p,
            q_tot,
            phase_type,
            maxiter,
        )
        q_pt = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
        e_int = internal_energy(param_set, T, q_pt)
    else
        q_pt = PhasePartition(q_tot)
        T = air_temperature_from_ideal_gas_law(param_set, p, ρ, q_pt)
        e_int = internal_energy(param_set, T, q_pt)
    end
    return PhaseEquil{FT}(ρ, p, e_int, q_tot, T)
end


"""
    PhaseEquil_pθq(param_set, θ_liq_ice, q_tot[, maxiter, temperature_tol, sat_adjust_method])

Constructs a [`PhaseEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_tot` total specific humidity
 - `temperature_tol` temperature tolerance for saturation adjustment
 - `maxiter` maximum iterations for saturation adjustment
 - `sat_adjust_method` the numerical method to use.
"""
function PhaseEquil_pθq(
    param_set::APS,
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    maxiter::IT = nothing,
    temperature_tol::FTT = nothing,
    ::Type{sat_adjust_method} = RS.SecantMethod,
) where {FT <: Real, IT <: ITERTYPE, FTT <: TOLTYPE(FT), sat_adjust_method}
    maxiter === nothing && (maxiter = 50)
    temperature_tol === nothing && (temperature_tol = FT(1e-3))
    phase_type = PhaseEquil{FT}
    q_tot_safe = clamp(q_tot, FT(0), FT(1))
    T = saturation_adjustment_given_pθq(
        sat_adjust_method,
        param_set,
        p,
        θ_liq_ice,
        q_tot_safe,
        phase_type,
        maxiter,
        temperature_tol,
    )
    q_pt = PhasePartition_equil_given_p(param_set, T, p, q_tot_safe, phase_type)
    ρ = air_density(param_set, T, p, q_pt)
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseEquil{FT}(ρ, p, e_int, q_tot_safe, T)
end


#####
##### Non-equilibrium states
#####

"""
     PhaseNonEquil{FT} <: ThermodynamicState

A thermodynamic state assuming thermodynamic non-equilibrium (therefore, temperature can
be computed directly).

# Constructors

    PhaseNonEquil(param_set, e_int, q::PhasePartition, ρ)

# Fields

$(DocStringExtensions.FIELDS)

"""
struct PhaseNonEquil{FT} <: AbstractPhaseNonEquil{FT}
    "internal energy"
    e_int::FT
    "density of air (potentially moist)"
    ρ::FT
    "phase partition"
    q::PhasePartition{FT}
end
function PhaseNonEquil(
    param_set::APS,
    e_int::FT,
    ρ::FT,
    q::PhasePartition{FT} = q_pt_0(FT),
) where {FT}
    return PhaseNonEquil{FT}(e_int, ρ, q)
end

"""
    PhaseNonEquil_ρTq(param_set, ρ, T, q_pt)

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `T` temperature
 - `q_pt` phase partition
"""
function PhaseNonEquil_ρTq(
    param_set::APS,
    ρ::FT,
    T::FT,
    q_pt::PhasePartition{FT},
) where {FT <: Real}
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseNonEquil{FT}(e_int, ρ, q_pt)
end

"""
    PhaseNonEquil_ρθq(param_set, ρ, θ_liq_ice, q_pt)

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_pt` phase partition
and, optionally
 - `potential_temperature_tol` potential temperature for non-linear equation solve
 - `maxiter` maximum iterations for non-linear equation solve
"""
function PhaseNonEquil_ρθq(
    param_set::APS,
    ρ::FT,
    θ_liq_ice::FT,
    q_pt::PhasePartition{FT},
    maxiter::Int = 10,
    potential_temperature_tol::FT = FT(1e-2),
) where {FT <: Real}
    phase_type = PhaseNonEquil{FT}
    tol = RS.ResidualTolerance(potential_temperature_tol)
    T = air_temperature_given_ρθq_nonlinear(
        param_set,
        ρ,
        θ_liq_ice,
        maxiter,
        tol,
        q_pt,
    )
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseNonEquil{FT}(e_int, ρ, q_pt)
end

"""
    PhaseNonEquil_pθq(param_set, p, θ_liq_ice, q_pt)

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_pt` phase partition
"""
function PhaseNonEquil_pθq(
    param_set::APS,
    p::FT,
    θ_liq_ice::FT,
    q_pt::PhasePartition{FT},
) where {FT <: Real}
    T = air_temperature_given_pθq(param_set, p, θ_liq_ice, q_pt)
    ρ = air_density(param_set, T, p, q_pt)
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseNonEquil{FT}(e_int, ρ, q_pt)
end

"""
    PhaseNonEquil_pTq(param_set, p, T, q_pt)

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` air temperature
 - `q_pt` phase partition
"""
function PhaseNonEquil_pTq(
    param_set::APS,
    p::FT,
    T::FT,
    q_pt::PhasePartition{FT},
) where {FT <: Real}
    ρ = air_density(param_set, T, p, q_pt)
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseNonEquil{FT}(e_int, ρ, q_pt)
end

"""
    PhaseNonEquil_peq(param_set, p, e_int, q_pt)

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `e_int` internal energy
 - `q_pt` phase partition
"""
function PhaseNonEquil_peq(
    param_set::APS,
    p::FT,
    e_int::FT,
    q_pt::PhasePartition{FT},
) where {FT <: Real}
    T = air_temperature(param_set, e_int, q_pt)
    ρ = air_density(param_set, T, p, q_pt)
    return PhaseNonEquil{FT}(e_int, ρ, q_pt)
end

"""
    PhaseNonEquil_ρpq(param_set, ρ, p, q_pt)

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `p` pressure
 - `q_pt` phase partition
"""
function PhaseNonEquil_ρpq(
    param_set::APS,
    ρ::FT,
    p::FT,
    q_pt::PhasePartition{FT},
) where {FT <: Real}
    T = air_temperature_from_ideal_gas_law(param_set, p, ρ, q_pt)
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseNonEquil{FT}(e_int, ρ, q_pt)
end
