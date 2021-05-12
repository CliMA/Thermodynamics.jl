# Tested Profiles

Thermodynamics.jl is tested using a set of profiles, or thermodynamic state regimes, specified in [`tested_profiles`](@ref Thermodynamics.tested_profiles).

## Dry Phase

```@example
using Thermodynamics
MT = Thermodynamics
using CLIMAParameters
using CLIMAParameters.Planet
using Plots

struct EarthParameterSet <: AbstractEarthParameterSet end;
const param_set = EarthParameterSet();

FT = Float64;
e_int, ρ, q_tot, q_pt, T, p, θ_liq_ice = MT.tested_profiles(param_set, 50, FT);

mask_dry = q_tot .≈ 0;
ρ_dry = ρ[mask_dry];
T_dry = T[mask_dry];

scatter(ρ_dry, T_dry, xlabel="density [kg/m^3]", ylabel="T [K]", title="Tested states for dry thermodynamic phase", legend=false);
savefig("tested_profiles_dry.svg");
```
![](tested_profiles_dry.svg)

## Moist Phase in thermodynamic equilibrium

Here is a 2D representation:
```@example
using Thermodynamics
MT = Thermodynamics
using CLIMAParameters
using CLIMAParameters.Planet
using Plots

struct EarthParameterSet <: AbstractEarthParameterSet end;
const param_set = EarthParameterSet();

FT = Float64;
e_int, ρ, q_tot, q_pt, T, p, θ_liq_ice = MT.tested_profiles(param_set, 50, FT);
scatter(ρ, T, xlabel="density [kg/m^3]", ylabel="T [K]", marker_z=q_tot, title="Tested states for moist thermodynamic phase", label="q_tot");
savefig("tested_profiles.svg")
```
![](tested_profiles.svg)

## Moist Phase in thermodynamic non-equilibrium

In progress...
