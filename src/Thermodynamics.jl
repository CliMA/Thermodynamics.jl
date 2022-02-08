"""
    Thermodynamics

Moist thermodynamic functions, e.g., for air pressure (atmosphere equation
of state), latent heats of phase transitions, saturation vapor pressures, and
saturation specific humidities.


## AbstractParameterSet's

Many functions defined in this module rely on CLIMAParameters.jl.
CLIMAParameters.jl defines several functions (e.g., many planet
parameters). For example, to compute the mole-mass ratio:

```julia
using CLIMAParameters.Planet: molmass_ratio
using CLIMAParameters: AbstractEarthParameterSet
struct EarthParameterSet <: AbstractEarthParameterSet end
param_set = EarthParameterSet()
_molmass_ratio = molmass_ratio(param_set)
```

Because these parameters are widely used throughout this module,
`param_set` is an argument for many Thermodynamics functions.

## Numerical methods for saturation adjustment

Saturation adjustment function are designed to accept
 - `sat_adjust_method` a type used to dispatch which numerical method to use

and a function to return an instance of the numerical method. For example:

 - `sa_numerical_method_ρpq` returns an instance of the numerical
    method. One of these functions must be defined for the particular
    numerical method and the particular formulation (`ρ-p-q_tot` in this case).

The currently supported numerical methods, in RootSolvers.jl, are:
 - `NewtonsMethod` uses Newton method with analytic gradients
 - `NewtonsMethodAD` uses Newton method with autodiff
 - `SecantMethod` uses Secant method
 - `RegulaFalsiMethod` uses Regula-Falsi method
"""
module Thermodynamics

import DocStringExtensions
const DSE = DocStringExtensions

import RootSolvers
const RS = RootSolvers

import KernelAbstractions
const KA = KernelAbstractions

#import CLIMAParameters
#const CP = CLIMAParameters
#const CPP = CP.Planet
#const APS = CP.AbstractParameterSet

####
#### ThermodynamicsParameters object
####

#=
# When we get proper name conventions, otherwise we use aliases below
struct ThermodynamicsParameters{FT} 
    thermodynamics_temperature_reference::FT
    mean_sea_level_pressure::FT
    gas_constant_dry_air::FT
    isobaric_specific_heat_dry_air::FT
    isochoric_specific_heat_dry_air::FT
    isobaric_specific_heat_vapor::FT
    isochoric_specific_heat_vapor::FT
    isobaric_specific_heat_liquid::FT
    isochoric_specific_heat_liquid::FT
    isobaric_specific_heat_ice::FT
    isochoric_specific_heat_ice::FT
    latent_heat_vaporization_at_reference    ::FT
    latent_heat_sublimation_at_reference::FT
    latent_heat_fusion_at_reference::FT
    specific_internal_energy_vapor_at_reference::FT
    specific_internal_energy_ice_at_reference::FT
    pressure_triple_point::FT
    temperature_triple_point::FT
    gas_constant_vapor::FT
    temperature_water_freeze::FT
    temperature_saturation_adjustment_min::FT
    temperature_saturation_adjustment_max::FT
    entropy_reference_temperature::FT
    entropy_dry_air::FT
    entropy_water_vapor::FT
    adiabatic_exponent_dry_air::FT
    gas_constant::FT
    molar_mass_dry_air::FT
    molar_mass_water::FT
    temperature_mean_at_reference::FT
    temperature_min_at_reference::FT
    gravitational_acceleration::FT
    temperature_homogenous_nucleation::FT
end

function ThermodynamicsParameters(param_set::NamedTuple)

    # Used in thermodynamics, from parameter file
    thermodynamics_temperature_reference = param_set.thermodynamics_temperature_reference
    mean_sea_level_pressure = param_set.mean_sea_level_pressure
    isobaric_specific_heat_dry_air = param_set.isobaric_specific_heat_dry_air
    isobaric_specific_heat_vapor = param_set.isobaric_specific_heat_vapor
    isobaric_specific_heat_liquid = param_set.isobaric_specific_heat_liquid
    isobaric_specific_heat_ice = param_set.isobaric_specific_heat_ice
    latent_heat_vaporization_at_reference = param_set.latent_heat_vaporization_at_reference
    latent_heat_sublimation_at_reference = param_set.latent_heat_sublimation_at_reference
    pressure_triple_point = param_set.pressure_triple_point
    temperature_triple_point = param_set.temperature_triple_point
    temperature_water_freeze = param_set.temperature_water_freeze
    temperature_saturation_adjustment_min = param_set.temperature_saturation_adjustment_min
    temperature_saturation_adjustment_max = param_set.temperature_saturation_adjustment_max 
    entropy_reference_temperature = param_set.entropy_reference_temperature
    entropy_dry_air = param_set.entropy_dry_air
    entropy_water_vapor = param_set.entropy_water_vapor
    adiabatic_exponent_dry_air = param_set.adiabatic_exponent_dry_air
    gas_constant = param_set.gas_constant
    molar_mass_dry_air = param_set.molar_mass_dry_air
    molar_mass_water = param_set.molar_mass_water
    temperature_mean_at_reference = param_set.temperature_mean_at_reference
    temperature_min_at_reference = param_set.temperature_min_at_reference
    gravitational_acceleration = param_set.gravitational_acceleration
    temperature_homogenous_nucleation = param_set.temperature_homogenous_nucleation

    # derived parameters from parameter file
    gas_constant_dry_air = gas_constant/ molar_mass_dry_air
    molar_mass_ratio_dry_air_water = molar_mass_dry_air / molar_mass_water
    gas_constant_vapor = gas_constant / molar_mass_water
    latent_heat_fusion_at_reference = latent_heat_sublimation_at_reference - latent_heat_vaporization_at_reference
    specific_internal_energy_vapor_at_reference =  latent_heat_vaporization_at_reference - gas_constant_vapor * thermodynamics_temperature_reference
    specific_internal_energy_ice_at_reference =  latent_heat_fusion_at_reference
    isobaric_specific_heat_dry_air = gas_constant_dry_air / adiabatic_exponent_dry_air
    isochoric_specific_heat_dry_air = isobaric_specific_heat_dry_air -  gas_constant_dry_air 
    isochoric_specific_heat_vapor = isobaric_specific_heat_vapor - gas_constant_vapor
    isochoric_specific_heat_liquid = isobaric_specific_heat_liquid
    isochoric_specific_heat_ice = isobaric_specific_heat_ice
    
    return ThermoDynamicsParameters{typeof(param_set)[1]}(
        thermodynamics_temperature_reference,
        mean_sea_level_pressure,
        gas_constant_dry_air,
        isobaric_specific_heat_dry_air,
        isochoric_specific_heat_dry_air,
        isobaric_specific_heat_vapor,
        isochoric_specific_heat_vapor,
        isobaric_specific_heat_liquid,
        isochoric_specific_heat_liquid,
        isobaric_specific_heat_ice,
        isochoric_specific_heat_ice,
        latent_heat_vaporization_at_reference    ,
        latent_heat_sublimation_at_reference,
        latent_heat_fusion_at_reference,
        specific_internal_energy_vapor_at_reference,
        specific_internal_energy_ice_at_reference,
        pressure_triple_point,
        temperature_triple_point,
        gas_constant_vapor,
        temperature_water_freeze,
        temperature_saturation_adjustment_min,
        temperature_saturation_adjustment_max,
        entropy_reference_temperature,
        entropy_dry_air,
        entropy_water_vapor,
        adiabatic_exponent_dry_air,
        gas_constant,
        molar_mass_dry_air,    
        molar_mass_water,
        temperature_mean_at_reference,
        temperature_min_at_reference,
        gravitational_acceleration,
        temperature_homogenous_nucleation
    )
end
=#

struct ThermodynamicsParameters{FT}
    T_0::FT
    MSLP::FT
    R_d::FT
    cp_d::FT
    cv_d::FT
    cp_v::FT
    cv_v::FT
    cp_l::FT
    cv_l::FT
    cp_i::FT
    cv_i::FT
    LH_v0::FT
    LH_s0::FT
    LH_f0::FT
    e_int_v0::FT
    e_int_i0::FT
    press_triple::FT
    T_triple::FT
    R_v::FT
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
    molmass_ratio::FT
    T_icenuc::FT
end

function ThermodynamicsParameters(param_set::NamedTuple)

    # Used in thermodynamics, from parameter file
    T_0 = param_set.T_0
    MSLP = param_set.MSLP
    cp_v = param_set.cp_v
    cp_l = param_set.cp_l
    cp_i = param_set.cp_i
    LH_v0 = param_set.LH_v0
    LH_s0 = param_set.LH_s0
    press_triple = param_set.press_triple
    T_triple = param_set.T_triple
    T_freeze = param_set.T_freeze
    T_min = param_set.T_min
    T_max = param_set.T_max
    entropy_reference_temperature = param_set.entropy_reference_temperature
    entropy_dry_air = param_set.entropy_dry_air
    entropy_water_vapor = param_set.entropy_water_vapor
    kappa_d = param_set.kappa_d
    gas_constant = param_set.gas_constant
    molmass_dryair = param_set.molmass_dryair
    molmass_water = param_set.molmass_water
    T_surf_ref = param_set.T_surf_ref
    T_min_ref = param_set.T_min_ref
    grav = param_set.grav
    T_icenuc = param_set.T_icenuc

    # derived parameters from parameter file
    R_d = gas_constant / molmass_dryair
    molmass_ratio = molmass_dryair / molmass_water
    R_v = gas_constant / molmass_water
    LH_f0 = LH_s0 - LH_v0
    e_int_v0 = LH_v0 - R_v * T_0
    e_int_i0 = LH_f0
    cp_d = R_d / kappa_d
    cv_d = cp_d - R_d
    cv_v = cp_v - R_v
    cv_l = cp_l
    cv_i = cp_i

    return ThermodynamicsParameters{typeof(param_set[1])}(
        T_0,
        MSLP,
        R_d,
        cp_d,
        cv_d,
        cp_v,
        cv_v,
        cp_l,
        cv_l,
        cp_i,
        cv_i,
        LH_v0,
        LH_s0,
        LH_f0,
        e_int_v0,
        e_int_i0,
        press_triple,
        T_triple,
        R_v,
        T_freeze,
        T_min,
        T_max,
        entropy_reference_temperature,
        entropy_dry_air,
        entropy_water_vapor,
        kappa_d,
        gas_constant,
        molmass_dryair,
        molmass_water,
        T_surf_ref,
        T_min_ref,
        grav,
        molmass_ratio,
        T_icenuc,
    )
end

const TPS = ThermodynamicsParameters



Base.broadcastable(param_set::ThermodynamicsParameters) = Ref(param_set)

# Allow users to skip error on non-convergence
# by importing:
# ```julia
# import Thermodynamics
# Thermodynamics.error_on_non_convergence() = false
# ```
# Error on convergence must be the default
# behavior because this can result in printing
# very large logs resulting in CI to seemingly hang.
error_on_non_convergence() = true

# Allow users to skip printing warnings on non-convergence
print_warning() = true

@inline q_pt_0(::Type{FT}) where {FT} = PhasePartition(FT(0), FT(0), FT(0))

include("states.jl")
include("relations.jl")
include("isentropic.jl")
include("config_numerical_method.jl")
include("TemperatureProfiles.jl")
include("TestedProfiles.jl")

Base.broadcastable(dap::DryAdiabaticProcess) = Ref(dap)
Base.broadcastable(phase::Phase) = Ref(phase)

end #module Thermodynamics.jl
