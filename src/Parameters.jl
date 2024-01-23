module Parameters

export ThermodynamicsParameters

abstract type AbstractThermodynamicsParameters{FT} end
const ATP = AbstractThermodynamicsParameters

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
Base.@kwdef struct ThermodynamicsParameters{FT} <: AbstractThermodynamicsParameters{FT}
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
Base.eltype(::AbstractThermodynamicsParameters{FT}) where {FT} = FT

# wrappers
for fn in fieldnames(ThermodynamicsParameters)
    @eval $(fn)(ps::ThermodynamicsParameters) = ps.$(fn)
end

# Derived parameters
R_d(ps::ATP) = gas_constant(ps) / molmass_dryair(ps)
R_v(ps::ATP) = gas_constant(ps) / molmass_water(ps)
molmass_ratio(ps::ATP) = molmass_dryair(ps) / molmass_water(ps)
LH_f0(ps::ATP) = LH_s0(ps) - LH_v0(ps)
e_int_v0(ps::ATP) = LH_v0(ps) - R_v(ps) * T_0(ps)
e_int_i0(ps::ATP) = LH_f0(ps)
cp_d(ps::ATP) = R_d(ps) / kappa_d(ps)
cv_d(ps::ATP) = cp_d(ps) - R_d(ps)
cv_v(ps::ATP) = cp_v(ps) - R_v(ps)
cv_l(ps::ATP) = cp_l(ps)
cv_i(ps::ATP) = cp_i(ps)

#=
#####
##### Hierarchical parameter set
#####

struct ConstitutiveProperties{FT} <: AbstractThermodynamicsParameters{FT}
    gas_constant               :: FT
    dry_air_molar_mass         :: FT
    water_molar_mass           :: FT
end

"""
    ConstitutiveProperties(FT; gas_constant               = 8.3144598,
                               dry_air_molar_mass         = 0.02897,
                               water_molar_mass           = 0.018015)

Construct a set of parameters that define the density of moist air,

```math
ρ = p / Rᵐ(q) T,
```

where ``p`` is pressure, ``T`` is temperature, ``q`` defines the partition
of total mass into vapor, liqiud, and ice mass fractions, and
``Rᵐ`` is the effective specific gas constant for the mixture,

```math
Rᵐ(q) =
```

where

For more information see [reference docs].
"""
function ConstitutiveProperties(FT = Float64;
                                gas_constant       = 8.3144598,
                                dry_air_molar_mass = 0.02897,
                                water_molar_mass   = 0.018015)

    return ConstitutiveProperties(convert(FT, gas_constant),
                                  convert(FT, dry_air_molar_mass),
                                  convert(FT, water_molar_mass))
end

const CP = ConstitutiveProperties

gas_constant(p::CP)   = p.gas_constant
molmass_dryair(p::CP) = p.dry_air_molar_mass
molmass_water(p::CP)  = p.water_molar_mass

struct HeatCapacityProperties{FT} <: AbstractThermodynamicsParameters{FT}
    dry_air_adiabatic_exponent :: FT
    water_vapor_heat_capacity  :: FT
    liquid_water_heat_capacity :: FT
    ice_heat_capacity          :: FT
end

function HeatCapacityProperties(FT = Float64;
                           dry_air_adiabatic_exponent = 2/7,
                           water_vapor_heat_capacity = 1859,
                           liquid_water_heat_capacity = 4181,
                           ice_heat_capacity = 2100)
                                  
    return HeatCapacityProperties(convert(FT, dry_air_adiabatic_exponent),
                             convert(FT, water_vapor_heat_capacity),
                             convert(FT, liquid_water_heat_capacity),
                             convert(FT, ice_heat_capacity))
end

const HCP = HeatCapacityProperties
cp_v(p::HCP)    = p.water_vapor_heat_capacity
cp_l(p::HCP)    = p.liquid_water_heat_capacity
cp_i(p::HCP)    = p.ice_heat_capacity
kappa_d(p::HCP) = p.dry_air_adiabatic_exponent

struct PhaseTransitionParameters{FT} <: AbstractThermodynamicsParameters{FT}
    reference_vaporization_enthalpy :: FT
    reference_sublimation_enthalpy  :: FT
    reference_temperature           :: FT
    triple_point_temperature        :: FT
    triple_point_pressure           :: FT
    water_freezing_temperature      :: FT
    ice_nucleation_temperature      :: FT
end

function PhaseTransitionParameters(FT = Float64;
                                   reference_vaporization_enthalpy = 2500800,
                                   reference_sublimation_enthalpy = 2834400,
                                   reference_temperature = 273.16,
                                   triple_point_temperature = 273.16,
                                   triple_point_pressure = 611.657,
                                   water_freezing_temperature = 273.16,
                                   ice_nucleation_temperature = 233)

   return PhaseTransitionParameters(convert(FT, reference_vaporization_enthalpy),
                                    convert(FT, reference_sublimation_enthalpy),
                                    convert(FT, reference_temperature),
                                    convert(FT, triple_point_temperature),
                                    convert(FT, triple_point_pressure),
                                    convert(FT, water_freezing_temperature),
                                    convert(FT, ice_nucleation_temperature))
end

const PTP = PhaseTransitionParameters
LH_v0(p::PTP)        = p.reference_vaporization_enthalpy
LH_s0(p::PTP)        = p.reference_sublimation_enthalpy
T_freeze(p::PTP)     = p.water_freezing_temperature
T_triple(p::PTP)     = p.triple_point_temperature
T_icenuc(p::PTP)     = p.ice_nucleation_temperature
press_triple(p::PTP) = p.triple_point_pressure
T_0(p::PTP)          = p.reference_temperature

struct HierarchicalThermodynamicsParameters{FT} <: AbstractThermodynamicsParameters{FT}
    constitutive       :: ConstitutiveProperties{FT}
    phase_transitions  :: PhaseTransitionParameters{FT}
    heat_capacity      :: HeatCapacityProperties{FT}
end

function HierarchicalThermodynamicsParameters(FT = Float64;
                                              constitutive = ConstitutiveProperties(FT),
                                              phase_transitions = PhaseTransitionParameters(FT),
                                              heat_capacity = HeatCapacityProperties(FT))

    return HierarchicalThermodynamicsParameters(constitutive, phase_transitions, heat_capacity)
end

const HTP = HierarchicalThermodynamicsParameters

gas_constant(p::HTP)   = gas_constant(p.constitutive)
molmass_dryair(p::HTP) = molmass_dryair(p.constitutive)
molmass_water(p::HTP)  = molmass_water(p.constitutive)
kappa_d(p::HTP)        = kappa_d(p.heat_capacity)
LH_v0(p::HTP)          = LH_v0(p.phase_transitions)
LH_s0(p::HTP)          = LH_s0(p.phase_transitions)
cp_v(p::HTP)           = cp_v(p.heat_capacity)
cp_l(p::HTP)           = cp_l(p.heat_capacity)
cp_i(p::HTP)           = cp_i(p.heat_capacity)
T_freeze(p::HTP)       = T_freeze(p.phase_transitions)
T_triple(p::HTP)       = T_triple(p.phase_transitions)
T_icenuc(p::HTP)       = T_icenuc(p.phase_transitions)
press_triple(p::HTP)   = press_triple(p.phase_transitions)
T_0(p::HTP)            = T_0(p.phase_transitions)
=#

end # module

