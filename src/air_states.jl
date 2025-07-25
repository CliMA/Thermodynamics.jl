# Constructors for thermodynamic states given various input variables

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

A thermodynamic state representing the complete thermodynamic properties of a moist air parcel.
All `ThermodynamicState` subtypes provide access to functions to compute all other thermodynamic
properties through the equation of state and thermodynamic relations.

The state can be initialized using various thermodynamic formulations (via its subtypes),
each representing different assumptions about phase equilibrium and the specific variables
used to define the state.
"""
abstract type ThermodynamicState{FT} end

abstract type AbstractPhaseDry{FT} <: ThermodynamicState{FT} end
abstract type AbstractPhaseEquil{FT} <: ThermodynamicState{FT} end
abstract type AbstractPhaseNonEquil{FT} <: ThermodynamicState{FT} end

Base.eltype(::ThermodynamicState{FT}) where {FT} = FT

"""
    PhasePartition

Represents the mass fractions of the moist air mixture (the partitioning of water substance 
between vapor, liquid, and ice phases).

The total specific humidity `q_tot` represents the total water content, while `q_liq` and
`q_ice` represent the liquid and ice specific humidities, respectively. The vapor specific
humidity is computed as `q_vap = q_tot - q_liq - q_ice`.

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

A dry thermodynamic state representing air with no water vapor (`q_tot = 0`).
This state assumes the air parcel contains only dry air components.

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

Constructs a [`PhaseDry`](@ref) thermodynamic state from density and internal energy, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` air density
 - `e_int` specific internal energy

This constructor directly stores the provided density and internal energy without any
additional computations, assuming the air is completely dry.
"""
PhaseDry_ρe(param_set::APS, ρ::FT, e_int::FT) where {FT} =
    PhaseDry{FT}(e_int, ρ)
PhaseDry_ρe(param_set::APS, ρ, e_int) =
    PhaseDry_ρe(param_set, promote(ρ, e_int)...)

"""
    PhaseDry_pT(param_set, p, T)

Constructs a [`PhaseDry`](@ref) thermodynamic state from pressure and temperature, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature

The internal energy is computed from the temperature using the dry air equation of state,
and the density is computed from the ideal gas law using the pressure and temperature.
"""
@inline function PhaseDry_pT(param_set::APS, p::FT, T::FT) where {FT <: Real}
    e_int = internal_energy(param_set, T)
    ρ = air_density(param_set, T, p)
    return PhaseDry{FT}(e_int, ρ)
end
PhaseDry_pT(param_set::APS, p, T) = PhaseDry_pT(param_set, promote(p, T)...)

"""
    PhaseDry_pe(param_set, p, e_int)

Constructs a [`PhaseDry`](@ref) thermodynamic state from pressure and internal energy, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `e_int` specific internal energy

The temperature is computed from the internal energy using the dry air equation of state,
and the density is computed from the ideal gas law using the pressure and temperature.
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

Constructs a [`PhaseDry`](@ref) thermodynamic state from pressure and specific enthalpy, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `h` specific enthalpy

The temperature is computed from the specific enthalpy using the dry air equation of state,
and the density is computed from the ideal gas law using the pressure and temperature.
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

Constructs a [`PhaseDry`](@ref) thermodynamic state from density and dry potential temperature, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `θ_dry` dry potential temperature

The temperature is computed from the density and potential temperature using the dry air
equation of state, and the internal energy is computed from the temperature.
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

Constructs a [`PhaseDry`](@ref) thermodynamic state from pressure and dry potential temperature, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `θ_dry` dry potential temperature

The temperature is computed from the pressure and potential temperature using the Exner function,
and the density is computed from the ideal gas law using the pressure and temperature.
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

Constructs a [`PhaseDry`](@ref) thermodynamic state from density and temperature, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `T` temperature

The internal energy is computed directly from the temperature using the dry air equation of state.
"""
@inline function PhaseDry_ρT(param_set::APS, ρ::FT, T::FT) where {FT <: Real}
    e_int = internal_energy(param_set, T)
    return PhaseDry{FT}(e_int, ρ)
end
PhaseDry_ρT(param_set::APS, ρ, T) = PhaseDry_ρT(param_set, promote(ρ, T)...)

"""
    PhaseDry_ρp(param_set, ρ, p)

Constructs a [`PhaseDry`](@ref) thermodynamic state from density and pressure, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `p` pressure

The temperature is computed from the ideal gas law using the pressure and density,
and the internal energy is computed from the temperature.
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

A thermodynamic state assuming thermodynamic equilibrium between water phases.
This state assumes that the water vapor is in equilibrium with liquid and/or ice,
requiring saturation adjustment to compute the temperature and phase partitioning.

The state stores the density, pressure, internal energy, total specific humidity,
and the computed temperature from saturation adjustment.

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

Constructs a [`PhaseEquil`](@ref) thermodynamic state from density, internal energy, and total specific humidity, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `e_int` specific internal energy
 - `q_tot` total specific humidity
and, optionally
 - `maxiter` maximum iterations for saturation adjustment (default: 8)
 - `relative_temperature_tol` relative temperature tolerance for saturation adjustment (default: 1e-4)
 - `sat_adjust_method` the numerical method to use for saturation adjustment (default: NewtonsMethod)
    See the [`Thermodynamics`](@ref) for options.
 - `T_guess` initial guess for temperature in saturation adjustment

The temperature is computed using saturation adjustment to ensure thermodynamic equilibrium,
and the pressure is computed from the equation of state using the temperature and density.
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
    PhaseEquil_ρθq(param_set, ρ, θ_liq_ice, q_tot[, maxiter, relative_temperature_tol, T_guess])

Constructs a [`PhaseEquil`](@ref) thermodynamic state from density, liquid-ice potential temperature, and total specific humidity, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_tot` total specific humidity
and, optionally
 - `maxiter` maximum iterations for saturation adjustment (default: 36)
 - `relative_temperature_tol` relative temperature tolerance for saturation adjustment (default: 1e-4)
 - `T_guess` initial guess for temperature in saturation adjustment

The temperature is computed using saturation adjustment with respect to the liquid-ice potential temperature,
and the pressure and internal energy are computed from the equation of state.
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

Constructs a [`PhaseEquil`](@ref) thermodynamic state from density, temperature, and total specific humidity, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `T` temperature
 - `q_tot` total specific humidity

The phase partitioning is computed assuming thermodynamic equilibrium at the given temperature,
and the pressure and internal energy are computed from the equation of state.
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

Constructs a [`PhaseEquil`](@ref) thermodynamic state from pressure, temperature, and total specific humidity, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `q_tot` total specific humidity

The phase partitioning is computed assuming thermodynamic equilibrium at the given temperature,
and the density and internal energy are computed from the equation of state.
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

Constructs a [`PhaseEquil`](@ref) thermodynamic state from pressure, internal energy, and total specific humidity, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `e_int` specific internal energy
 - `q_tot` total specific humidity
and, optionally
 - `maxiter` maximum iterations for saturation adjustment (default: 40)
 - `relative_temperature_tol` relative temperature tolerance for saturation adjustment (default: 1e-4)
 - `sat_adjust_method` the numerical method to use for saturation adjustment (default: SecantMethod)
    See the [`Thermodynamics`](@ref) for options.
 - `T_guess` initial guess for temperature in saturation adjustment

The temperature is computed using saturation adjustment given pressure and internal energy,
and the density is computed from the equation of state using the pressure and temperature.
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

Constructs a [`PhaseEquil`](@ref) thermodynamic state from pressure, specific enthalpy, and total specific humidity, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `h` specific enthalpy
 - `q_tot` total specific humidity
and, optionally
 - `maxiter` maximum iterations for saturation adjustment (default: 40)
 - `relative_temperature_tol` relative temperature tolerance for saturation adjustment (default: 1e-4)
 - `sat_adjust_method` the numerical method to use for saturation adjustment (default: SecantMethod)
    See the [`Thermodynamics`](@ref) for options.
 - `T_guess` initial guess for temperature in saturation adjustment

The temperature is computed using saturation adjustment given pressure and specific enthalpy,
and the density and internal energy are computed from the equation of state.
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

Constructs a [`PhaseEquil`](@ref) thermodynamic state from density, pressure, and total specific humidity, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `p` pressure
 - `q_tot` total specific humidity
and, optionally
 - `perform_sat_adjust` Boolean indicating whether to perform saturation adjustment (default: false)
 - `maxiter` maximum number of iterations to perform in saturation adjustment (default: 5)
 - `sat_adjust_method` the numerical method to use for saturation adjustment (default: NewtonsMethodAD)
    See the [`Thermodynamics`](@ref) for options.
 - `T_guess` initial guess for temperature in saturation adjustment

If `perform_sat_adjust` is true, the temperature is computed using saturation adjustment.
Otherwise, the temperature is computed directly from the ideal gas law.
The internal energy is computed from the temperature and phase partitioning.

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
    PhaseEquil_pθq(param_set, p, θ_liq_ice, q_tot[, maxiter, relative_temperature_tol, sat_adjust_method, T_guess])

Constructs a [`PhaseEquil`](@ref) thermodynamic state from pressure, liquid-ice potential temperature, and total specific humidity, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_tot` total specific humidity
and, optionally
 - `maxiter` maximum iterations for saturation adjustment (default: 50)
 - `relative_temperature_tol` relative temperature tolerance for saturation adjustment (default: 1e-4)
 - `sat_adjust_method` the numerical method to use for saturation adjustment (default: SecantMethod)
    See the [`Thermodynamics`](@ref) for options.
 - `T_guess` initial guess for temperature in saturation adjustment

The temperature is computed using saturation adjustment with respect to the liquid-ice potential temperature,
and the density and internal energy are computed from the equation of state.
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

#####
##### Non-equilibrium states
#####

"""
     PhaseNonEquil{FT} <: ThermodynamicState

A thermodynamic state assuming thermodynamic non-equilibrium between water phases.
This state allows for arbitrary phase partitioning without requiring saturation adjustment,
enabling direct computation of temperature from the given thermodynamic variables.

The state stores the internal energy, density, and a complete phase partition
specifying the distribution of water substance between vapor, liquid, and ice phases.

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

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from density, temperature, and phase partition, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `T` temperature
 - `q_pt` phase partition

The internal energy is computed from the temperature and phase partition using the equation of state.
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

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from density, liquid-ice potential temperature, and phase partition, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_pt` phase partition
and, optionally
 - `maxiter` maximum iterations for non-linear equation solve (default: 10)
 - `relative_temperature_tol` relative temperature tolerance for non-linear equation solve (default: 1e-2)

The temperature is computed from the density and liquid-ice potential temperature using a non-linear solver,
and the internal energy is computed from the temperature and phase partition.
"""
@inline function PhaseNonEquil_ρθq(
    param_set::APS,
    ρ::FT,
    θ_liq_ice::FT,
    q_pt::PhasePartition{FT},
    maxiter::Int = 10,
    relative_temperature_tol::FT = FT(1e-4),
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

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from pressure, liquid-ice potential temperature, and phase partition, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_pt` phase partition

The temperature is computed from the pressure and liquid-ice potential temperature,
and the density and internal energy are computed from the equation of state.
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

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from pressure, temperature, and phase partition, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` air temperature
 - `q_pt` phase partition

The density is computed from the ideal gas law using the pressure and temperature,
and the internal energy is computed from the temperature and phase partition.
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

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from pressure, internal energy, and phase partition, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `e_int` specific internal energy
 - `q_pt` phase partition

The temperature is computed from the internal energy and phase partition using the equation of state,
and the density is computed from the ideal gas law using the pressure and temperature.
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
    PhaseNonEquil_phq(param_set, p, h, q_pt)

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from pressure, specific enthalpy, and phase partition, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `h` specific enthalpy
 - `q_pt` phase partition

The temperature is computed from the specific enthalpy and phase partition using the equation of state,
and the density and internal energy are computed from the equation of state.
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

Constructs a [`PhaseNonEquil`](@ref) thermodynamic state from density, pressure, and phase partition, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` density
 - `p` pressure
 - `q_pt` phase partition

The temperature is computed from the ideal gas law using the pressure and density,
and the internal energy is computed from the temperature and phase partition.
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
