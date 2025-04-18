export PhasePartition
# Thermodynamic states
export ThermodynamicState,
    PhaseDry,
    PhaseDry_ρe,
    PhaseDry_ρT,
    PhaseDry_pT,
    PhaseDry_pe,
    PhaseDry_ph,
    PhaseDry_ρθ,
    PhaseDry_pθ,
    PhaseDry_ρp,
    PhaseEquil,
    PhaseEquil_ρeq,
    PhaseEquil_ρTq,
    PhaseEquil_pTq,
    PhaseEquil_pTRH,
    PhaseEquil_peq,
    PhaseEquil_phq,
    PhaseEquil_ρθq,
    PhaseEquil_pθq,
    PhaseEquil_ρpq,
    PhaseNonEquil,
    PhaseNonEquil_ρTq,
    PhaseNonEquil_pTq,
    PhaseNonEquil_ρθq,
    PhaseNonEquil_pθq,
    PhaseNonEquil_peq,
    PhaseNonEquil_phq,
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
    @inline function PhasePartition(tot::FT, liq::FT, ice::FT) where {FT}
        q_tot_safe = max(tot, 0)
        q_liq_safe = max(liq, 0)
        q_ice_safe = max(ice, 0)
        return new{FT}(q_tot_safe, q_liq_safe, q_ice_safe)
    end
end

@inline Base.zero(::Type{PhasePartition{FT}}) where {FT} =
    PhasePartition(FT(0), FT(0), FT(0))

@inline PhasePartition(q_tot::FT, q_liq::FT) where {FT <: Real} =
    PhasePartition(q_tot, q_liq, zero(FT))
@inline PhasePartition(q_tot::FT) where {FT <: Real} =
    PhasePartition(q_tot, zero(FT), zero(FT))

Base.convert(::Type{PhasePartition{FT}}, q_pt::PhasePartition) where {FT} =
    PhasePartition(FT(q_pt.tot), FT(q_pt.liq), FT(q_pt.ice))

function promote_phase_partition(x, q_pt::PhasePartition)
    (x′, tot, liq, ice) = promote(x, q_pt.tot, q_pt.liq, q_pt.ice)
    return (x′, PhasePartition(tot, liq, ice))
end
function promote_phase_partition(x1, x2, q_pt::PhasePartition)
    (x1′, x2′, tot, liq, ice) = promote(x1, x2, q_pt.tot, q_pt.liq, q_pt.ice)
    return (x1′, x2′, PhasePartition(tot, liq, ice))
end

const ITERTYPE = Union{Int, Nothing}
const TOLTYPE = Union{Real, Nothing}

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
    # TODO: swap order of variables (breaking change)
    "internal energy"
    e_int::FT
    "density of dry air"
    ρ::FT
end
@inline PhaseDry(param_set::APS, e_int::FT, ρ::FT) where {FT} =
    PhaseDry{FT}(e_int, ρ)

Base.zero(::Type{PhaseDry{FT}}) where {FT} = PhaseDry{FT}(0, 0)

Base.convert(::Type{PhaseDry{FT}}, ts::PhaseDry) where {FT} =
    PhaseDry{FT}(ts.e_int, ts.ρ)

"""
    PhaseDry_ρe(param_set, ρ, e_int)

A dry thermodynamic state (`q_tot = 0`) from

 - `param_set` an [`ThermodynamicsParameters`](@ref Parameters.ThermodynamicsParameters) for more details
 - `ρ`
 - `e_int` internal energy
"""
PhaseDry_ρe(param_set::APS, ρ::FT, e_int::FT) where {FT} =
    PhaseDry{FT}(e_int, ρ)
PhaseDry_ρe(param_set::APS, ρ, e_int) =
    PhaseDry_ρe(param_set, promote(ρ, e_int)...)

"""
    PhaseDry_pT(param_set, p, T)

Constructs a [`PhaseDry`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
"""
@inline function PhaseDry_pT(param_set::APS, p::FT, T::FT) where {FT <: Real}
    e_int = internal_energy(param_set, T)
    ρ = air_density(param_set, T, p)
    return PhaseDry{FT}(e_int, ρ)
end
PhaseDry_pT(param_set::APS, p, T) = PhaseDry_pT(param_set, promote(p, T)...)

"""
    PhaseDry_pe(param_set, p, e_int)

Constructs a [`PhaseDry`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `e_int` internal energy
"""
@inline function PhaseDry_pe(
    param_set::APS,
    p::FT,
    e_int::FT,
) where {FT <: Real}
    T = air_temperature(param_set, e_int)
    ρ = air_density(param_set, T, p)
    return PhaseDry{FT}(e_int, ρ)
end
PhaseDry_pe(param_set::APS, p, e_int) =
    PhaseDry_pe(param_set, promote(p, e_int)...)

"""
     PhaseDry_ph(param_set, p, h)

 Constructs a [`PhaseDry`](@ref) thermodynamic state from:

  - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
  - `p` pressure
  - `h` specific enthalpy
 """
@inline function PhaseDry_ph(param_set::APS, p::FT, h::FT) where {FT <: Real}
    T = air_temperature_from_enthalpy(param_set, h)
    ρ = air_density(param_set, T, p)
    e_int = internal_energy(param_set, T)
    return PhaseDry{FT}(e_int, ρ)
end
PhaseDry_ph(param_set::APS, p, h) = PhaseDry_ph(param_set, promote(p, h)...)

"""
    PhaseDry_ρθ(param_set, ρ, θ_dry)

Constructs a [`PhaseDry`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `θ_dry` dry potential temperature
"""
@inline function PhaseDry_ρθ(
    param_set::APS,
    ρ::FT,
    θ_dry::FT,
) where {FT <: Real}
    T = air_temperature_given_ρθq(param_set, ρ, θ_dry)
    e_int = internal_energy(param_set, T)
    return PhaseDry{FT}(e_int, ρ)
end
PhaseDry_ρθ(param_set::APS, ρ, θ_dry) =
    PhaseDry_ρθ(param_set, promote(ρ, θ_dry)...)

"""
    PhaseDry_pθ(param_set, p, θ_dry)

Constructs a [`PhaseDry`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `θ_dry` dry potential temperature
"""
@inline function PhaseDry_pθ(
    param_set::APS,
    p::FT,
    θ_dry::FT,
) where {FT <: Real}
    T = exner_given_pressure(param_set, p) * θ_dry
    e_int = internal_energy(param_set, T)
    ρ = air_density(param_set, T, p)
    return PhaseDry{FT}(e_int, ρ)
end
PhaseDry_pθ(param_set::APS, p, θ_dry) =
    PhaseDry_pθ(param_set, promote(p, θ_dry)...)

"""
    PhaseDry_ρT(param_set, ρ, T)

Constructs a [`PhaseDry`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `T` temperature
"""
@inline function PhaseDry_ρT(param_set::APS, ρ::FT, T::FT) where {FT <: Real}
    e_int = internal_energy(param_set, T)
    return PhaseDry{FT}(e_int, ρ)
end
PhaseDry_ρT(param_set::APS, ρ, T) = PhaseDry_ρT(param_set, promote(ρ, T)...)

"""
    PhaseDry_ρp(param_set, ρ, p)

Constructs a [`PhaseDry`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `p` pressure
"""
@inline function PhaseDry_ρp(param_set::APS, ρ::FT, p::FT) where {FT <: Real}
    T = air_temperature_from_ideal_gas_law(param_set, p, ρ)
    e_int = internal_energy(param_set, T)
    return PhaseDry{FT}(e_int, ρ)
end
PhaseDry_ρp(param_set::APS, ρ, p) = PhaseDry_ρp(param_set, promote(ρ, p)...)

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

@inline Base.zero(::Type{PhaseEquil{FT}}) where {FT} =
    PhaseEquil{FT}(0, 0, 0, 0, 0)

Base.convert(::Type{PhaseEquil{FT}}, ts::PhaseEquil) where {FT} =
    PhaseEquil{FT}(ts.ρ, ts.p, ts.e_int, ts.q_tot, ts.T)

"""
    PhaseEquil_ρeq(param_set, ρ, e_int, q_tot[, maxiter, relative_temperature_tol, sat_adjust_method, T_guess])

Moist thermodynamic phase, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `e_int` internal energy
 - `q_tot` total specific humidity
and, optionally
 - `maxiter` maximum iterations for saturation adjustment
 - `relative_temperature_tol` relative temperature tolerance for saturation adjustment
 - `sat_adjust_method` the numerical method to use.
    See the [`Thermodynamics`](@ref) for options.
 - `T_guess` initial guess for temperature in saturation adjustment
"""
@inline function PhaseEquil_ρeq(
    param_set::APS,
    ρ::FT,
    e_int::FT,
    q_tot::FT,
    maxiter::IT = nothing,
    relative_temperature_tol::TT = nothing,
    ::Type{sat_adjust_method} = RS.NewtonsMethod,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, sat_adjust_method, IT <: ITERTYPE, TT <: TOLTYPE}
    maxiter === nothing && (maxiter = 8)
    relative_temperature_tol === nothing &&
        (relative_temperature_tol = FT(1e-4))
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
        relative_temperature_tol,
        T_guess,
    )
    q_pt = PhasePartition_equil(param_set, T, ρ, q_tot_safe, phase_type)
    p = air_pressure(param_set, T, ρ, q_pt)
    return PhaseEquil{FT}(ρ, p, e_int, q_tot_safe, T)
end
PhaseEquil_ρeq(param_set::APS, ρ, e_int, q_tot, args...) =
    PhaseEquil_ρeq(param_set, promote(ρ, e_int, q_tot)..., args...)

# Convenience method for comparing Numerical
# methods without having to specify maxiter
# and relative_temperature_tol. maxiter and relative_temperature_tol
# should be in sync with the PhaseEquil(...) constructor
@inline function PhaseEquil_dev_only(
    param_set::APS,
    ρ::FT,
    e_int::FT,
    q_tot::FT;
    maxiter::Int = 8,
    relative_temperature_tol::FT = FT(1e-4),
    sat_adjust_method::Type{NM} = RS.NewtonsMethod,
) where {FT <: Real, NM}
    return PhaseEquil_ρeq(
        param_set,
        ρ,
        e_int,
        q_tot,
        maxiter,
        relative_temperature_tol,
        sat_adjust_method,
    )
end

"""
    PhaseEquil_ρθq(param_set, ρ, θ_liq_ice, q_tot[, maxiter, relative_temperature_tol, sat_adjust_method, T_guess])

Constructs a [`PhaseEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_tot` total specific humidity
 - `relative_temperature_tol` relative temperature tolerance for saturation adjustment
 - `maxiter` maximum iterations for saturation adjustment
 - `T_guess` initial guess for temperature in saturation adjustment
"""
@inline function PhaseEquil_ρθq(
    param_set::APS,
    ρ::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    maxiter::IT = nothing,
    relative_temperature_tol::TT = nothing,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, IT <: ITERTYPE, TT <: TOLTYPE}
    maxiter === nothing && (maxiter = 36)
    relative_temperature_tol === nothing &&
        (relative_temperature_tol = FT(1e-4))
    phase_type = PhaseEquil{FT}
    tol = RS.RelativeSolutionTolerance(relative_temperature_tol)
    T = saturation_adjustment_given_ρθq(
        param_set,
        ρ,
        θ_liq_ice,
        q_tot,
        phase_type,
        maxiter,
        tol,
        T_guess,
    )
    q_pt = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    p = air_pressure(param_set, T, ρ, q_pt)
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseEquil{FT}(ρ, p, e_int, q_tot, T)
end
PhaseEquil_ρθq(param_set::APS, ρ, θ_liq_ice, q_tot, args...) =
    PhaseEquil_ρθq(param_set, promote(ρ, θ_liq_ice, q_tot)..., args...)

"""
    PhaseEquil_ρTq(param_set, ρ, T, q_tot)

Constructs a [`PhaseEquil`](@ref) thermodynamic state from temperature.

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `T` temperature
 - `q_tot` total specific humidity
"""
@inline function PhaseEquil_ρTq(
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
PhaseEquil_ρTq(param_set::APS, ρ, T, q_tot) =
    PhaseEquil_ρTq(param_set, promote(ρ, T, q_tot)...)

"""
    PhaseEquil_pTq(param_set, p, T, q_tot)

Constructs a [`PhaseEquil`](@ref) thermodynamic state from temperature.

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `q_tot` total specific humidity
"""
@inline function PhaseEquil_pTq(
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
PhaseEquil_pTq(param_set::APS, p, T, q_tot) =
    PhaseEquil_pTq(param_set, promote(p, T, q_tot)...)

"""
    PhaseEquil_peq(param_set, p, e_int, q_tot[, maxiter, relative_temperature_tol, sat_adjust_method, T_guess])

Constructs a [`PhaseEquil`](@ref) thermodynamic state from temperature.

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `e_int` temperature
 - `q_tot` total specific humidity
 - `T_guess` initial guess for temperature in saturation adjustment
"""
@inline function PhaseEquil_peq(
    param_set::APS,
    p::FT,
    e_int::FT,
    q_tot::FT,
    maxiter::IT = nothing,
    relative_temperature_tol::TT = nothing,
    ::Type{sat_adjust_method} = RS.SecantMethod,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, sat_adjust_method, IT <: ITERTYPE, TT <: TOLTYPE}
    maxiter === nothing && (maxiter = 40)
    relative_temperature_tol === nothing &&
        (relative_temperature_tol = FT(1e-4))
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
        relative_temperature_tol,
        T_guess,
    )
    q_pt = PhasePartition_equil_given_p(param_set, T, p, q_tot_safe, phase_type)
    ρ = air_density(param_set, T, p, q_pt)
    return PhaseEquil{FT}(ρ, p, e_int, q_tot_safe, T)
end
PhaseEquil_peq(param_set::APS, p, e_int, q_tot, args...) =
    PhaseEquil_peq(param_set, promote(p, e_int, q_tot)..., args...)


"""
    PhaseEquil_phq(param_set, p, h, q_tot[, maxiter, relative_temperature_tol, sat_adjust_method, T_guess])

Constructs a [`PhaseEquil`](@ref) thermodynamic state from temperature.

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `h` specific enthalpy
 - `q_tot` total specific humidity
 - `T_guess` initial guess for temperature in saturation adjustment
"""
@inline function PhaseEquil_phq(
    param_set::APS,
    p::FT,
    h::FT,
    q_tot::FT,
    maxiter::IT = nothing,
    relative_temperature_tol::TT = nothing,
    ::Type{sat_adjust_method} = RS.SecantMethod,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, sat_adjust_method, IT <: ITERTYPE, TT <: TOLTYPE}
    maxiter === nothing && (maxiter = 40)
    relative_temperature_tol === nothing &&
        (relative_temperature_tol = FT(1e-4))
    phase_type = PhaseEquil{FT}
    q_tot_safe = clamp(q_tot, FT(0), FT(1))
    T = saturation_adjustment_given_phq(
        sat_adjust_method,
        param_set,
        p,
        h,
        q_tot_safe,
        phase_type,
        maxiter,
        relative_temperature_tol,
        T_guess,
    )
    q_pt = PhasePartition_equil_given_p(param_set, T, p, q_tot_safe, phase_type)
    ρ = air_density(param_set, T, p, q_pt)
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseEquil{FT}(ρ, p, e_int, q_tot_safe, T)
end
PhaseEquil_phq(param_set::APS, p, h, q_tot, args...) =
    PhaseEquil_phq(param_set, promote(p, h, q_tot)..., args...)

"""
    PhaseEquil_ρpq(param_set, ρ, p, q_tot[, perform_sat_adjust=true, maxiter, sat_adjust_method, T_guess])

Constructs a [`PhaseEquil`](@ref) thermodynamic state from temperature.

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `p` pressure
 - `q_tot` total specific humidity
 - `perform_sat_adjust` Bool indicating to perform saturation adjustment
 - `maxiter` maximum number of iterations to perform in saturation adjustment
 - `sat_adjust_method` the numerical method to use.
    See the [`Thermodynamics`](@ref) for options.
 - `T_guess` initial guess for temperature in saturation adjustment

TODO: change input argument order: perform_sat_adjust is
      unique to this constructor, so it should be last.
      (breaking change)
"""
@inline function PhaseEquil_ρpq(
    param_set::APS,
    ρ::FT,
    p::FT,
    q_tot::FT,
    perform_sat_adjust = false,
    maxiter::IT = nothing,
    relative_temperature_tol::TT = nothing,
    ::Type{sat_adjust_method} = RS.NewtonsMethodAD,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, sat_adjust_method, IT <: ITERTYPE, TT <: TOLTYPE}
    maxiter === nothing && (maxiter = 5)
    relative_temperature_tol === nothing &&
        (relative_temperature_tol = FT(1e-4))
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
            relative_temperature_tol,
            T_guess,
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
PhaseEquil_ρpq(param_set::APS, ρ, p, q_tot, args...) =
    PhaseEquil_ρpq(param_set, promote(ρ, p, q_tot)..., args...)


"""
    PhaseEquil_pθq(param_set, θ_liq_ice, q_tot[, maxiter, relative_temperature_tol, sat_adjust_method, T_guess])

Constructs a [`PhaseEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_tot` total specific humidity
 - `relative_temperature_tol` relative temperature tolerance for saturation adjustment
 - `maxiter` maximum iterations for saturation adjustment
 - `sat_adjust_method` the numerical method to use.
 - `T_guess` initial guess for temperature in saturation adjustment
"""
@inline function PhaseEquil_pθq(
    param_set::APS,
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    maxiter::IT = nothing,
    relative_temperature_tol::TT = nothing,
    ::Type{sat_adjust_method} = RS.SecantMethod,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, IT <: ITERTYPE, TT <: TOLTYPE, sat_adjust_method}
    maxiter === nothing && (maxiter = 50)
    relative_temperature_tol === nothing &&
        (relative_temperature_tol = FT(1e-4))
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
        relative_temperature_tol,
        T_guess,
    )
    q_pt = PhasePartition_equil_given_p(param_set, T, p, q_tot_safe, phase_type)
    ρ = air_density(param_set, T, p, q_pt)
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseEquil{FT}(ρ, p, e_int, q_tot_safe, T)
end
PhaseEquil_pθq(param_set::APS, p, θ_liq_ice, q_tot, args...) =
    PhaseEquil_pθq(param_set, promote(p, θ_liq_ice, q_tot)..., args...)

"""
    PhaseEquil_pTRH(param_set, p, T, RH)

Constructs a [`PhaseEquil`](@ref) thermodynamic state from temperature.

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `RH` relative humidity
"""
@inline function PhaseEquil_pTRH(
    param_set::APS,
    p::FT,
    T::FT,
    RH::FT,
) where {FT <: Real}
    phase_type = PhaseEquil{FT}
    p_vap_sat = saturation_vapor_pressure(param_set, T, Liquid())
    p_vap = RH * p_vap_sat
    mmr = TP.molmass_ratio(param_set)
    q_tot = p_vap * mmr / (p - (1 - mmr) * p_vap)
    q_tot_safe = clamp(q_tot, FT(0), FT(1))
    q_pt = PhasePartition_equil_given_p(param_set, T, p, q_tot_safe, phase_type)
    ρ = air_density(param_set, T, p, q_pt)
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseEquil{FT}(ρ, p, e_int, q_tot_safe, T)
end
PhaseEquil_pTRH(param_set::APS, p, T, q_tot) =
    PhaseEquil_pTRH(param_set, promote(p, T, q_tot)...)

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
@inline Base.zero(::Type{PhaseNonEquil{FT}}) where {FT} =
    PhaseNonEquil{FT}(0, 0, zero(PhasePartition{FT}))

Base.convert(::Type{PhaseNonEquil{FT}}, ts::PhaseNonEquil) where {FT} =
    PhaseNonEquil{FT}(ts.e_int, ts.ρ, ts.q)

@inline function PhaseNonEquil(
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
@inline function PhaseNonEquil_ρTq(
    param_set::APS,
    ρ::FT,
    T::FT,
    q_pt::PhasePartition{FT},
) where {FT <: Real}
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseNonEquil{FT}(e_int, ρ, q_pt)
end
PhaseNonEquil_ρTq(param_set::APS, ρ, T, q_pt) =
    PhaseNonEquil_ρTq(param_set, promote_phase_partition(ρ, T, q_pt)...)

"""
    PhaseNonEquil_ρθq(param_set, ρ, θ_liq_ice, q_pt)

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_pt` phase partition
and, optionally
 - `relative_temperature_tol` potential temperature for non-linear equation solve
 - `maxiter` maximum iterations for non-linear equation solve
"""
@inline function PhaseNonEquil_ρθq(
    param_set::APS,
    ρ::FT,
    θ_liq_ice::FT,
    q_pt::PhasePartition{FT},
    maxiter::Int = 10,
    relative_temperature_tol::FT = FT(1e-2),
) where {FT <: Real}
    phase_type = PhaseNonEquil{FT}
    tol = RS.RelativeSolutionTolerance(relative_temperature_tol)
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
PhaseNonEquil_ρθq(param_set::APS, ρ, θ_liq_ice, q_pt) =
    PhaseNonEquil_ρθq(param_set, promote_phase_partition(ρ, θ_liq_ice, q_pt)...)

"""
    PhaseNonEquil_pθq(param_set, p, θ_liq_ice, q_pt)

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_pt` phase partition
"""
@inline function PhaseNonEquil_pθq(
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
PhaseNonEquil_pθq(param_set::APS, p, θ_liq_ice, q_pt) =
    PhaseNonEquil_pθq(param_set, promote_phase_partition(p, θ_liq_ice, q_pt)...)

"""
    PhaseNonEquil_pTq(param_set, p, T, q_pt)

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` air temperature
 - `q_pt` phase partition
"""
@inline function PhaseNonEquil_pTq(
    param_set::APS,
    p::FT,
    T::FT,
    q_pt::PhasePartition{FT},
) where {FT <: Real}
    ρ = air_density(param_set, T, p, q_pt)
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseNonEquil{FT}(e_int, ρ, q_pt)
end
PhaseNonEquil_pTq(param_set::APS, p, T, q_pt) =
    PhaseNonEquil_pTq(param_set, promote_phase_partition(p, T, q_pt)...)

"""
    PhaseNonEquil_peq(param_set, p, e_int, q_pt)

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `e_int` internal energy
 - `q_pt` phase partition
"""
@inline function PhaseNonEquil_peq(
    param_set::APS,
    p::FT,
    e_int::FT,
    q_pt::PhasePartition{FT},
) where {FT <: Real}
    T = air_temperature(param_set, e_int, q_pt)
    ρ = air_density(param_set, T, p, q_pt)
    return PhaseNonEquil{FT}(e_int, ρ, q_pt)
end
PhaseNonEquil_peq(param_set::APS, p, e_int, q_pt) =
    PhaseNonEquil_peq(param_set, promote_phase_partition(p, e_int, q_pt)...)

"""
    PhaseNonEquil_phq(param_set, p, e_int, q_pt)

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `h` specific enthalpy
 - `q_pt` phase partition
"""
@inline function PhaseNonEquil_phq(
    param_set::APS,
    p::FT,
    h::FT,
    q_pt::PhasePartition{FT},
) where {FT <: Real}
    T = air_temperature_from_enthalpy(param_set, h, q_pt)
    ρ = air_density(param_set, T, p, q_pt)
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseNonEquil{FT}(e_int, ρ, q_pt)
end
PhaseNonEquil_phq(param_set::APS, p, h, q_pt) =
    PhaseNonEquil_phq(param_set, promote_phase_partition(p, h, q_pt)...)

"""
    PhaseNonEquil_ρpq(param_set, ρ, p, q_pt)

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `p` pressure
 - `q_pt` phase partition
"""
@inline function PhaseNonEquil_ρpq(
    param_set::APS,
    ρ::FT,
    p::FT,
    q_pt::PhasePartition{FT},
) where {FT <: Real}
    T = air_temperature_from_ideal_gas_law(param_set, p, ρ, q_pt)
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseNonEquil{FT}(e_int, ρ, q_pt)
end
PhaseNonEquil_ρpq(param_set::APS, ρ, p, q_pt) =
    PhaseNonEquil_ρpq(param_set, promote_phase_partition(ρ, p, q_pt)...)
