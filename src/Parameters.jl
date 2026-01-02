"""
    Parameters

Submodule defining thermodynamic parameter types and accessor functions.
"""
module Parameters

export ThermodynamicsParameters

abstract type AbstractThermodynamicsParameters{FT} end
const ATP = AbstractThermodynamicsParameters

"""
    ThermodynamicsParameters

Parameters for Thermodynamics.jl.

# Example
```
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

"Internal energy of fusion at reference temperature `T_0`."
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
