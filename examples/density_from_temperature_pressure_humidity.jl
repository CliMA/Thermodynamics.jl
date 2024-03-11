#! format: off #src

# # Defining a simple parameter set and using it to compute density
#
# This script shows how to define a simple parameter set, and then using it to
# compute density as a function of pressure, temperature, and humidity.

# # Define parameters for computing the density of air
#
# First, we build a parameter set suitable for computing the density of
# _moist_ air --- that is, a mixture of dry air, water vapor, liquid droplets,
# and ice crystals (the latter two are called "condensates").

using Thermodynamics
using Thermodynamics.Parameters: AbstractThermodynamicsParameters

struct ConstitutiveParameters{FT} <: AbstractThermodynamicsParameters{FT}
    gas_constant::FT
    dry_air_molar_mass::FT
    water_molar_mass::FT
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
function ConstitutiveParameters(
    FT = Float64;
    gas_constant = 8.3144598,
    dry_air_molar_mass = 0.02897,
    water_molar_mass = 0.018015,
)

    return ConstitutiveParameters(
        convert(FT, gas_constant),
        convert(FT, dry_air_molar_mass),
        convert(FT, water_molar_mass),
    )
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

p₀ = 101325.0 # sea level pressure in Pascals (here taken to be mean sea level pressure)
T₀ = 273.15   # temperature in Kelvin

# Note that the above must be defined with the same float point precision as
# used for the parameters and PhasePartition.
# We're now ready to compute the density of dry air,

ρ = air_density(parameters, T₀, p₀, q_dry)

@show ρ

# Note that the above must be defined with the same float point precision as
# used for the parameters and PhasePartition.
# We're now ready to compute the density of dry air,

using JLD2

# Next, we load an atmospheric state correpsonding to atmospheric surface
# variables substampled from the JRA55 dataset, from the date Jan 1, 1991:

@load "JRA55_atmospheric_state_Jan_1_1991.jld2" q T p

# The variables `q`, `T`, and `p` correspond to the total specific humidity
# (a mass fraction), temperature (Kelvin), and sea level pressure (Pa)
#
# We use `q` to build a vector of PhasePartition,

qp = PhasePartition.(q)

# And then compute the density using the same parameters as before:

ρ = air_density.(parameters, T, p, qp)

# Finally, we plot the density as a function of temperature and specific humidity,

using CairoMakie

## Pressure range, centered around the mean sea level pressure defined above
pmax = maximum(abs, p)
dp = 3 / 4 * (pmax - p₀)
prange = (p₀ - dp, p₀ + dp)
pmap = :balance

## Compute temperature range
Tmin = minimum(T)
Tmax = maximum(T)
Trange = (Tmin, Tmax)
Tmap = :viridis

fig = Figure(size = (1000, 450))

axρ = Axis(fig[2, 1], xlabel = "Temperature (K) ", ylabel = "Density (kg m⁻³)")
axq = Axis(fig[2, 2], xlabel = "Specific humidity", ylabel = "Density (kg m⁻³)")

scatter!(axρ, T[:], ρ[:], color = p[:], colorrange = prange, colormap = pmap, alpha = 0.1)
scatter!(axq, q[:], ρ[:], color = T[:], colorrange = Trange, colormap = Tmap, alpha = 0.1)

Colorbar(fig[1, 1], label = "Pressure (Pa)", vertical = false, colorrange = prange, colormap = pmap)
Colorbar(fig[1, 2], label = "Temperature (K)", vertical = false, colorrange = Trange, colormap = Tmap)

current_figure()

#! format: on #src