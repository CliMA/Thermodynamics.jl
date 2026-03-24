Thermodynamics.jl Release Notes
========================

v1.1.0
--------

- Use a fixed `q_min` parameter as a lower bound for condensate specific humidity in the function `has_condensate`
  PR [315](https://github.com/CliMA/Thermodynamics.jl/pull/315)

v1.0.1
--------

- ![][[badge-🐛bugfix]] Fixed limiting behavior of saturation vapor pressure at very low temperatures.
  PR [313](https://github.com/CliMA/Thermodynamics.jl/pull/313)

v1.0.0
--------

- ![][[badge-💥breaking]] Removal of the deprecated object-oriented and state-based API (e.g., `ThermodynamicState`, `PhasePartition`, `PhaseEquil`). The package now relies exclusively on a stateless, completely functional API.
- ![][[badge-🚀performance]] Enforced zero-allocations across core routines (`air_temperature`, `air_pressure`, `saturation_adjustment` with fixed iterations) backed by explicit testing.
- ![][[badge-✨feature/enhancement]] Upgraded testing infrastructure: `Documenter.doctest` continuously checks all mathematical examples in the docstrings.
- ![][[badge-✨feature/enhancement]] Rewrite of documentation
  - Restructured documentation for better organization and clarity
  - Updated all function documentation to be consistent and comprehensive
  - Improved code examples and usage patterns throughout documentation
  - Enhanced cross-references and internal documentation links
  - Split monolithic src/relations.jl and test/relations.jl into multiple files

- ![][[badge-✨feature/enhancement]] Added new thermodynamic functions for export
  - `vapor_pressure_deficit` function for computing vapor pressure deficit
  - New methods for `partial_pressure_vapor` and `partial_pressure_dry` functions
  - Added comprehensive tests and physical consistency validation for these functions
  PR [259](https://github.com/CliMA/Thermodynamics.jl/pull/259)
  PR [263](https://github.com/CliMA/Thermodynamics.jl/pull/263)
- Renamed `specific_enthalpy*` to `enthalpy*`.
- Renamed `specific_entropy*` to `entropy*`.
- Renamed `latent_heat_liq_ice` to `humidity_weighted_latent_heat`.
- Renamed `air_temperature_given_ρp` to `air_temperature_given_pρq`.
- Renamed `air_temperature_from_enthalpy` to `air_temperature_given_hq`.
- Removed `phase_type` from `relative_humidity`.
- Removed `q_vap_saturation_generic`. Use `q_vap_saturation` instead.
- Removed `universal_gas_constant` (or `gas_constant`), `molar_mass_dryair` and `molar_mass_water` from thermo parameters.
- Fixed bug in `liquid_fraction` logic for nonequilibrium phases

v0.16.0
--------

- Added functional methods for cv_m and virtual (potential) temperature.

v0.15.0
--------

- Remove `q_vap_saturation_from_density` and `condensate`.
  PR [284](https://github.com/CliMA/Thermodynamics.jl/pull/284)
- Set error_on_non_covergence and print_warning to false by default.
  PR [283](https://github.com/CliMA/Thermodynamics.jl/pull/283)

v0.14.2
--------

- Add specific enthalpy functions. PR [282](https://github.com/CliMA/Thermodynamics.jl/pull/282)

v0.12.15
--------

- Fix inverse molmass bug. PR [258](https://github.com/CliMA/Thermodynamics.jl/pull/258)
- Thermodynamics.jl is no longer tested on Julia versions before 1.10.
  Please do not expect compatibility with those versions.
  PR [257](https://github.com/CliMA/Thermodynamics.jl/pull/257)

v0.12.14
--------

- Added an option to call `cp_m` without using `PhasePartition`
  PR [256](https://github.com/CliMA/Thermodynamics.jl/pull/256)

v0.12.13
-------

- ![][[badge-✨feature/enhancement]] Allow cp_m calculation without using PhasePartition.
  PR [254](https://github.com/CliMA/Thermodynamics.jl/pull/254)

v0.12.10
-------

- ![][[badge-🐛bugfix]] Asynchronous printing on the gpu has been fixed.
  PR [239](https://github.com/CliMA/Thermodynamics.jl/pull/239)

v0.12.9
-------

- ![][[badge-🐛bugfix]] Protest against zero division in relative humidity
  calculation and limit relative humidity between 0 and 1.
  PR [230](https://github.com/CliMA/Thermodynamics.jl/pull/230)

v0.12.8
-------

- ![][badge-🤖precisionΔ] Change the tolerance of PhaseEquil constructor to 1e-4
- ![][badge-🔥behavioralΔ] Change the definition of dry air internal energy and enthalpy

v0.12.7
-------

- ![][badge-🔥behavioralΔ] Change the upper limit of saturation specific humidity

v0.12.4
-------

- Upgraded to use ClimaParams.jl

v0.12.3
-------

- ![][badge-✨feature/enhancement] Additional Dual number support

v0.12.2
-------

- ![][badge-✨feature/enhancement] Additional Dual number support

v0.12.1
-------

- Started changelog
- ![][badge-✨feature/enhancement] Added support for Dual numbers

<!--
Contributors are welcome to begin the description of changelog items with badge(s) below. Here is a brief description of when to use badges for a particular pull request / set of changes:
 - 🔥behavioralΔ - behavioral changes. For example: a new model is used, yielding more accurate results.
 - 🤖precisionΔ - machine-precision changes. For example, swapping the order of summed arguments can result in machine-precision changes.
 - 💥breaking - breaking changes. For example: removing deprecated functions/types, removing support for functionality, API changes.
 - 🚀performance - performance improvements. For example: improving type inference, reducing allocations, or code hoisting.
 - ✨feature - new feature added. For example: adding support for a cubed-sphere grid
 - 🐛bugfix - bugfix. For example: fixing incorrect logic, resulting in incorrect results, or fixing code that otherwise might give a `MethodError`.
-->

[badge-🔥behavioralΔ]: https://img.shields.io/badge/🔥behavioralΔ-orange.svg
[badge-🤖precisionΔ]: https://img.shields.io/badge/🤖precisionΔ-black.svg
[badge-✨feature/enhancement]: https://img.shields.io/badge/feature/enhancement-blue.svg
[badge-🐛bugfix]: https://img.shields.io/badge/🐛bugfix-purple.svg
