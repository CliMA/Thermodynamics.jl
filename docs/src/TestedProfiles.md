# Tested Profiles

MoistThermodynamics.jl is tested using a set of profiles, or thermodynamic state regimes, specified in [`tested_profiles`](@ref MoistThermodynamics.tested_profiles).

## Dry Phase

```@example
using MoistThermodynamics
MT = MoistThermodynamics
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
using MoistThermodynamics
MT = MoistThermodynamics
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

And here is a 3D representation:
```@example
using MoistThermodynamics
MT = MoistThermodynamics
using CLIMAParameters
using CLIMAParameters.Planet
using Plots

struct EarthParameterSet <: AbstractEarthParameterSet end;
const param_set = EarthParameterSet();

FT = Float64;
e_int, ρ, q_tot, q_pt, T, p, θ_liq_ice = MT.tested_profiles(param_set, 50, FT);

# initialize a 3D plot with 1 empty series
plt = plot3d(
    1,
    xlim = (min(ρ...), max(ρ...)),
    ylim = (min(T...), max(T...)),
    zlim = (min(q_tot...), max(q_tot...)),
    xlabel="density [kg/m^3]",
    ylabel="T [K]",
    zlabel="total specific humidity []",
    legend=false,
    title = "Tested states for moist thermodynamic phase",
    marker = 2,
)

# build an animated gif by pushing new points to the plot, saving every nth frame
@gif for i=1:length(ρ)
    push!(plt, ρ[i], T[i], q_tot[i])
end every 5
```

## Moist Phase in thermodynamic non-equilibrium

In progress...
