# Thermodynamics.jl

Thermodynamics.jl provides flexible and performant functions for computing various thermodynamic variables for dry and moist atmospheres.

|||
|---------------------:|:----------------------------------------------|
| **Docs Build**       | [![docs build][docs-bld-img]][docs-bld-url]   |
| **Documentation**    | [![dev][docs-dev-img]][docs-dev-url]          |
| **GHA CI**           | [![gha ci][gha-ci-img]][gha-ci-url]           |
| **Code Coverage**    | [![codecov][codecov-img]][codecov-url]        |

[docs-bld-img]: https://github.com/CliMA/Thermodynamics.jl/actions/workflows/docs.yml/badge.svg
[docs-bld-url]: https://github.com/CliMA/Thermodynamics.jl/actions/workflows/docs.yml

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://CliMA.github.io/Thermodynamics.jl/dev/

[gha-ci-img]: https://github.com/CliMA/Thermodynamics.jl/actions/workflows/ci.yml/badge.svg
[gha-ci-url]: https://github.com/CliMA/Thermodynamics.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/CliMA/Thermodynamics.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/Thermodynamics.jl

# Principal design

This package uses an abstraction that leverages the idea of a thermodynamic state:

 - given two (or more) independent intrinsic thermodynamic properties, we can establish a thermodynamic state and
 - given a thermodynamic state, we can compute any thermodynamic property

## Example
For our example, we first load packages, and create a set of thermodynamic parameters, using a convenience constructor, offered through [CLIMAParameters.jl](https://github.com/CliMA/CLIMAParameters.jl). Then we create a thermodynamic state using density, liquid-ice potential temperature, and total specific humidity. Finally, we compute air temperature from the thermodynamic state.

```julia
import Thermodynamics as TD
using CLIMAParameters # load convenience parameter struct wrappers
params = TD.Parameters.ThermodynamicsParameters(Float64);

ts = TD.PhaseEquil_ρθq(params, 1.0, 374.0, 0.01);
T = TD.air_temperature(params, ts)
```

See a full list of different thermodynamic state constructors, in case you want to create a thermodynamic state with different variables, [here](https://clima.github.io/Thermodynamics.jl/dev/API/#Thermodynamic-State-Constructors).

See a full list of quantities that you can compute from a thermodynamic state, see our thermodynamic state-compatible methods [here](https://clima.github.io/Thermodynamics.jl/dev/API/#Thermodynamic-state-methods).
