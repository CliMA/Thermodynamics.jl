module Parameters

export ThermodynamicsParameters, Rv_over_Rd

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
    T_0::FT
    MSLP::FT
    p_ref_theta::FT
    cp_v::FT
    cp_l::FT
    cp_i::FT
    LH_v0::FT
    LH_s0::FT
    press_triple::FT
    T_triple::FT
    T_freeze::FT
    T_min::FT
    T_max::FT
    T_init_min::FT
    entropy_reference_temperature::FT
    entropy_dry_air::FT
    entropy_water_vapor::FT
    kappa_d::FT
    gas_constant::FT
    molmass_dryair::FT
    molmass_water::FT
    T_surf_ref::FT
    T_min_ref::FT
    grav::FT
    T_icenuc::FT
    pow_icenuc::FT
end

Base.broadcastable(ps::ATP) = tuple(ps)
Base.eltype(::ThermodynamicsParameters{FT}) where {FT} = FT

# wrappers
for fn in fieldnames(ThermodynamicsParameters)
    @eval $(fn)(ps::ThermodynamicsParameters) = ps.$(fn)
end

# Derived parameters
@inline R_d(ps::ATP) = ps.gas_constant / ps.molmass_dryair
@inline R_v(ps::ATP) = ps.gas_constant / ps.molmass_water
@inline Rv_over_Rd(ps::ATP) = ps.molmass_dryair / ps.molmass_water
@inline molmass_ratio(ps::ATP) = Rv_over_Rd(ps)  # For backward compatibility
@inline LH_f0(ps::ATP) = ps.LH_s0 - ps.LH_v0
@inline e_int_v0(ps::ATP) = ps.LH_v0 - R_v(ps) * ps.T_0
@inline e_int_i0(ps::ATP) = LH_f0(ps)
@inline cp_d(ps::ATP) = R_d(ps) / ps.kappa_d
@inline cv_d(ps::ATP) = cp_d(ps) - R_d(ps)
@inline cv_v(ps::ATP) = ps.cp_v - R_v(ps)
@inline cv_l(ps::ATP) = ps.cp_l
@inline cv_i(ps::ATP) = ps.cp_i

end
