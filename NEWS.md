Thermodynamics.jl Release Notes
========================

main
----

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
