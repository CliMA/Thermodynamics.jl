# Tested Profiles

Thermodynamics.jl is tested using a set of profiles specified in `src/TestedProfiles.jl`.

## Dry Phase

```@example
import Thermodynamics as TD
import Plots
import CLIMAParameters as CP
import Thermodynamics.InternalClimaParams as ICP
FT = Float64;
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
aliases = string.(fieldnames(ICP.ThermodynamicsParameters))
param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
param_set = ICP.ThermodynamicsParameters{FT}(; param_pairs...)

profiles = TD.TestedProfiles.PhaseDryProfiles(param_set, Array{FT});
(;T, ρ, z) = profiles
p1 = Plots.scatter(ρ, z./10^3, xlabel="Density [kg/m^3]", ylabel="z [km]", title="Density");
p2 = Plots.scatter(T, z./10^3, xlabel="Temperature [K]", ylabel="z [km]", title="Temperature");
Plots.plot(p1, p2, layout=(1,2))
Plots.savefig("tested_profiles_dry.svg");
```
![](tested_profiles_dry.svg)

## Moist Phase in thermodynamic equilibrium

```@example
import Thermodynamics as TD
import Plots
import CLIMAParameters as CP
import Thermodynamics.InternalClimaParams as ICP
FT = Float64;
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
aliases = string.(fieldnames(ICP.ThermodynamicsParameters))
param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
param_set = ICP.ThermodynamicsParameters{FT}(; param_pairs...)

profiles = TD.TestedProfiles.PhaseEquilProfiles(param_set, Array{FT});
(;T, ρ, q_tot, z) = profiles
p1 = Plots.scatter(ρ, z./10^3, xlabel="Density [kg/m^3]", ylabel="z [km]", title="Density");
p2 = Plots.scatter(T, z./10^3, xlabel="Temperature [K]", ylabel="z [km]", title="Temperature");
p3 = Plots.scatter(q_tot*1000, z./10^3, xlabel="Total specific\nhumidity [g/kg]", ylabel="z [km]", title="Total specific\nhumidity");
Plots.plot(p1, p2, p3, layout=(1,3))
Plots.savefig("tested_profiles_virt_temp.svg")
```
![](tested_profiles_virt_temp.svg)
