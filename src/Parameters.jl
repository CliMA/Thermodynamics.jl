module Parameters

export ThermodynamicsParameters

"""
    ThermodynamicsParameters

Parameters for Thermodynamics.jl.

# Example
```
import CLIMAParameters as CP
import Thermodynamics.Parameters as TP
FT = Float64;
toml_dict = CP.create_toml_dict(FT; dict_type = "alias");
aliases = string.(fieldnames(TP.ThermodynamicsParameters));
param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics");
param_set = TP.ThermodynamicsParameters{FT}(; param_pairs...);
```
"""
Base.@kwdef struct ThermodynamicsParameters{FT}
    T_0::FT
    MSLP::FT
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

const ATP = ThermodynamicsParameters

Base.broadcastable(ps::ATP) = Ref(ps)
Base.eltype(::ThermodynamicsParameters{FT}) where {FT} = FT

# wrappers
for fn in fieldnames(ATP)
    @eval $(fn)(ps::ATP) = ps.$(fn)
end

# Derived parameters
R_d(ps::ATP) = ps.gas_constant / ps.molmass_dryair
R_v(ps::ATP) = ps.gas_constant / ps.molmass_water
molmass_ratio(ps::ATP) = ps.molmass_dryair / ps.molmass_water
LH_f0(ps::ATP) = ps.LH_s0 - ps.LH_v0
e_int_v0(ps::ATP) = ps.LH_v0 - R_v(ps) * ps.T_0
e_int_i0(ps::ATP) = LH_f0(ps)
cp_d(ps::ATP) = R_d(ps) / ps.kappa_d
cv_d(ps::ATP) = cp_d(ps) - R_d(ps)
cv_v(ps::ATP) = ps.cp_v - R_v(ps)
cv_l(ps::ATP) = ps.cp_l
cv_i(ps::ATP) = ps.cp_i

end
