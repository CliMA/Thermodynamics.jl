"""
    Parameters

Submodule defining thermodynamic parameter types and accessor functions.
"""
module Parameters

export ThermodynamicsParameters

"""
    AbstractThermodynamicsParameters{FT}

Abstract supertype for thermodynamic parameter sets.

Subtypes:
- [`ThermodynamicsParameters`](@ref): default concrete parameter set.

Subtypes must implement accessor functions for all thermodynamic constants
(e.g., `T_0`, `R_d`, `cp_d`, …). See [`ThermodynamicsParameters`](@ref) for
the full list of required fields.
"""
abstract type AbstractThermodynamicsParameters{FT} end
const ATP = AbstractThermodynamicsParameters

"""
    ThermodynamicsParameters{FT} <: AbstractThermodynamicsParameters{FT}

Concrete parameter set for Thermodynamics.jl.

# Fields
- `T_0`: reference temperature [K].
- `T_triple`: triple-point temperature of water [K].
- `T_freeze`: freezing temperature of water [K].
- `T_icenuc`: temperature of homogeneous ice nucleation [K].
- `T_min`: minimum allowable temperature for saturation adjustment [K].
- `T_max`: maximum allowable temperature for saturation adjustment [K].
- `T_init_min`: minimum initial temperature guess for saturation adjustment [K].
- `T_surf_ref`: reference surface temperature [K].
- `T_min_ref`: minimum reference temperature (for temperature profiles) [K].
- `entropy_reference_temperature`: reference temperature for entropy [K].
- `MSLP`: mean sea-level pressure [Pa].
- `p_ref_theta`: reference pressure for potential temperature [Pa].
- `press_triple`: triple-point pressure of water [Pa].
- `R_d`: specific gas constant of dry air [J/(kg·K)].
- `R_v`: specific gas constant of water vapor [J/(kg·K)].
- `cp_d`: isobaric specific heat capacity of dry air [J/(kg·K)].
- `cp_v`: isobaric specific heat capacity of water vapor [J/(kg·K)].
- `cp_l`: isobaric specific heat capacity of liquid water [J/(kg·K)].
- `cp_i`: isobaric specific heat capacity of ice [J/(kg·K)].
- `LH_v0`: latent heat of vaporization at `T_0` [J/kg].
- `LH_s0`: latent heat of sublimation at `T_0` [J/kg].
- `entropy_dry_air`: reference specific entropy of dry air [J/(kg·K)].
- `entropy_water_vapor`: reference specific entropy of water vapor [J/(kg·K)].
- `grav`: gravitational acceleration [m/s²].
- `pow_icenuc`: exponent for power-law ice nucleation fraction [-].
- `q_min`: minimum condensate threshold for `has_condensate` [kg/kg].

# Examples
```julia
import ClimaParams as CP
import Thermodynamics.Parameters as TP

FT = Float64;
param_set = TP.ThermodynamicsParameters(FT)

# Alternatively:

toml_dict = CP.create_toml_dict(FT)
param_set = TP.ThermodynamicsParameters(toml_dict)
```
"""
Base.@kwdef struct ThermodynamicsParameters{FT} <:
                   AbstractThermodynamicsParameters{FT}
    # Reference temperatures
    T_0::FT
    T_triple::FT
    T_freeze::FT
    T_icenuc::FT
    T_min::FT
    T_max::FT
    T_init_min::FT
    T_surf_ref::FT
    T_min_ref::FT
    entropy_reference_temperature::FT
    # Reference pressures
    MSLP::FT
    p_ref_theta::FT
    press_triple::FT
    # Gas constants
    R_d::FT
    R_v::FT
    # Isobaric specific heat capacities
    cp_d::FT
    cp_v::FT
    cp_l::FT
    cp_i::FT
    # Latent heats at reference temperature
    LH_v0::FT
    LH_s0::FT
    # Entropy reference values
    entropy_dry_air::FT
    entropy_water_vapor::FT
    # Other
    grav::FT
    pow_icenuc::FT
    q_min::FT
end

Base.broadcastable(ps::ATP) = tuple(ps)
Base.eltype(::ThermodynamicsParameters{FT}) where {FT} = FT

# wrappers
for fn in fieldnames(ThermodynamicsParameters)
    @eval $(fn)(ps::ThermodynamicsParameters) = ps.$(fn)
end

# Derived parameters

"Ratio of gas constants `R_v / R_d`."
@inline Rv_over_Rd(ps::ATP) = R_v(ps) / R_d(ps)

"Latent heat of fusion at reference temperature `T_0`."
@inline LH_f0(ps::ATP) = LH_s0(ps) - LH_v0(ps)

"Internal energy of vaporization at reference temperature `T_0`."
@inline e_int_v0(ps::ATP) = LH_v0(ps) - R_v(ps) * T_0(ps)

"Internal energy of fusion at reference temperature `T_0`. Equal to the latent heat of fusion
because for incompressible condensed phases the `p Δv` term is negligible (`p Δv ≪ L_f`)."
@inline e_int_i0(ps::ATP) = LH_f0(ps)

"Ratio of dry air gas constant to isobaric specific heat `R_d / cp_d`."
@inline kappa_d(ps::ATP) = R_d(ps) / cp_d(ps)

"Isochoric specific heat capacity of dry air."
@inline cv_d(ps::ATP) = cp_d(ps) - R_d(ps)

"Isochoric specific heat capacity of water vapor."
@inline cv_v(ps::ATP) = cp_v(ps) - R_v(ps)

"Isochoric specific heat capacity of liquid water."
@inline cv_l(ps::ATP) = cp_l(ps)

"Isochoric specific heat capacity of ice."
@inline cv_i(ps::ATP) = cp_i(ps)

end
