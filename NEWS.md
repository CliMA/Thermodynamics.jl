Thermodynamics.jl Release Notes
========================

v1.xxx  (TODO: release once refactoring is done)
--------

- ![][[badge-✨feature/enhancement]] Complete rewrite of documentation 
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

main
----

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
[badge-💥breaking]: https://img.shields.io/badge/💥BREAKING-red.svg
[badge-🚀performance]: https://img.shields.io/badge/🚀performance-green.svg
[badge-✨feature/enhancement]: https://img.shields.io/badge/feature/enhancement-blue.svg
[badge-🐛bugfix]: https://img.shields.io/badge/🐛bugfix-purple.svg
