# Tutorial on defining parameter sets and performing thermodynamic computations.
#
# This script gives a short tutorial on defining thermodynamic parameter sets,
# and performing some basic computations.

# # Define parameters for computing the density of air
#
# First, we build a parameter set suitable for computing the density of
# _moist_ air --- that is, a mixture of dry air, water vapor, liquid droplets,
# and ice crystals (the latter two are called "condensates").

using Thermodynamics
using Thermodynamics.Parameters: AbstractThermodynamicsParameters

struct ConstitutiveParameters{FT} <: AbstractThermodynamicsParameters{FT}
    gas_constant       :: FT
    dry_air_molar_mass :: FT
    water_molar_mass   :: FT
end

"""
    ConstitutiveParameters(FT; gas_constant       = 8.3144598,
                               dry_air_molar_mass = 0.02897,
                               water_molar_mass   = 0.018015)

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
function ConstitutiveParameters(FT = Float64;
                                gas_constant       = 8.3144598,
                                dry_air_molar_mass = 0.02897,
                                water_molar_mass   = 0.018015)

    return ConstitutiveParameters(convert(FT, gas_constant),
                                convert(FT, dry_air_molar_mass),
                                convert(FT, water_molar_mass))
end

# Next, we define functions that return:
#   1. The specific gas constant for dry air
#   2. The specific gas constant for water vapor
#   3. The ratio between the dry air molar mass and water molar mass

const VTP = ConstitutiveParameters
using Thermodynamics: Parameters

Parameters.R_d(p::VTP) = p.gas_constant / p.dry_air_molar_mass
Parameters.R_v(p::VTP) = p.gas_constant / p.water_molar_mass
Parameters.molmass_ratio(p::VTP) = p.dry_air_molar_mass / p.water_molar_mass

# # The density of dry air

import Thermodynamics as AtmosphericThermodynamics

# To compute the density of dry air, we first build a default
# parameter set,

parameters = ConstitutiveParameters()

# Next, we construct a phase partitioning representing dry air --- the trivial
# case where the all mass ratios are 0.
#
# Note that the syntax for `PhasePartition` is
#
# ```
# PhasePartition(q_total, q_liquid, q_ice)
# ```
# where the q's are mass fractions of total water components,
# liquid droplet condensate, and ice crystal condensate.

q_dry = AtmosphericThermodynamics.PhasePartition(0.0)

# Finally, we define the pressure and temperature, which constitute
# the "state" of our atmosphere,

p = 101325.0 # pressure in Pascals (here taken to be mean sea level pressure)
T = 273.15   # temperature in Kelvin

# Note that the above must be defined with the same float point precision as
# used for the parameters and PhasePartition.

# We're now ready to compute the density of dry air,

ρ = air_density(parameters, T, p, q_dry)

struct HeatCapacities{FT} <: AbstractThermodynamicsParameters{FT}
    dry :: FT
    vapor :: FT
    liquid :: FT
    ice :: FT
end

const HC = HeatCapacities
Parameters.cp_d(heat_capacities::HC) = heat_capacities.dry
Parameters.cp_v(heat_capacities::HC) = heat_capacities.vapor
Parameters.cp_l(heat_capacities::HC) = heat_capacities.liquid
Parameters.cp_i(heat_capacities::HC) = heat_capacities.ice

struct PhaseTransitionParameters{FT} <: AbstractThermodynamicsParameters{FT}
    reference_fusion_enthalpy :: FT
    reference_vaporization_enthalpy :: FT
    reference_temperature :: FT
    triple_point_temperature :: FT
    triple_point_pressure :: FT
end

const PTP = PhaseTransitionParameters
Parameters.LH_v0(p::PTP) = p.reference_vaporization_enthalpy

struct MyThermodynamicsParameters{FT}
    constitutive :: ConstitutiveParameters{FT}
    phase_transitions :: PhaseTransitionParameters{FT}
    heat_capacities :: HeatCapacities{FT}
end

const MTP = MyThermodynamicsParameters

Parameters.R_d(p::MTP)           = R_d(p.constitutive)
Parameters.R_v(p::MTP)           = R_v(p.constitutive)
Parameters.molmass_ratio(p::MTP) = molmass_ratio(p.constitutive)
Parameters.LH_v0(p.MTP)          = LH_v0(p.phase_transitions)
Parameters.cp_d(p::MTP)          = cp_d(p.heat_capacities)
Parameters.cp_v(p::MTP)          = cp_v(p.heat_capacities)
Parameters.cp_l(p::MTP)          = cp_l(p.heat_capacities)
Parameters.cp_i(p::MTP)          = cp_i(p.heat_capacities)

# and then virtual temperature:
Tᵛ = virtual_temperature(parameters, T, ρ, q)

@info """ With the parameters $p, and for the state
        
        - pressure:    $p Pa
        - temperature: $T K

    We computed

        - density:             $ρ
        - virtual temperature: $Tᵛ
"""
