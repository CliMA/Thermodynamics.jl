Thermodynamics.jl Release Notes
========================

main
--------

### Saturation adjustment: correctness fixes (changes model output)

- ![][badge-🐛bugfix] The fixed-iteration saturation adjustment used by the GPU path could
  fail to converge and silently return a wrong temperature. The residual's derivative is discontinuous at the saturation boundary, so iterates could jump across the boundary and settle into a limit cycle. Over a sweep of 1335 saturated states, 12.6% exceeded the documented 0.1 K accuracy and 6% exceeded 10 K, with a worst case of 25.6 K; `converged` was hard-coded to `true` throughout. The solver now iterates on the analytic continuation of the saturated branch (a continuous residual *and* derivative), limits its step size, keeps the iterate positive, and selects the result against the exact unsaturated solution. Every one of those states now converges to the exact answer given enough iterations (worst error 0.0 K), and the residual error decreases monotonically with `maxiter` instead of stalling.

- ![][badge-🐛bugfix] The pressure-based formulations (`pe`, `ph`, `pθ_li`) computed their
  saturation humidity from a density that assumed no condensate, biasing the result. Fully converged solves of equilibrium states returned temperatures biased by −0.16 K at 300 K and 1000 hPa. These paths now obtain the equilibrium partition from the pressure directly, and round-trip to machine precision.

- ![][badge-🐛bugfix] `∂q_vap_sat_∂T` omitted the contribution of the temperature-dependent
  liquid fraction, making it wrong by up to 7.4% in the mixed-phase regime. The missing term has a closed form requiring no extra `exp` call, since the mixed-phase saturation vapor pressure is the geometric mean `(p_liq)^λ (p_ice)^(1-λ)`. The derivative now matches finite differences to round-off.

- The default `maxiter` of the convenience methods stays at `2`, which remains the right
  choice for time-stepping models: the iteration count is set by the gap between the
  unsaturated first guess and the solution, that gap stays under 4 K across the tested
  profiles, and two iterations hold the error below 2e-3 K there. What changed is the failure
  mode — remaining error is now a truncation error that decreases monotonically with
  `maxiter`, rather than a limit cycle that no iteration count could escape. Per-iteration
  cost also fell about 11%.

- The `converged` field returned by the full signature is now computed from the final Newton
  increment rather than hard-coded to `true`.

### Robustness

- ![][badge-🐛bugfix] `liquid_fraction_ramp` threw a `DomainError` below `T_icenuc` for any
  non-integer `pow_icenuc`, an unrecoverable failure inside a GPU kernel. The ramp argument
  is now clamped before exponentiation.
- ![][badge-🐛bugfix] `liquid_fraction` and `q_vap_saturation_from_pressure_calc` returned
  `Union{Float32, Float64}` when called with arguments wider than the parameter set, which
  de-optimizes GPU kernels. Both branches now share a type.
- Temperature guards are applied consistently across all root-finding closures; previously
  only the `ρe` non-Newton path was guarded.

### Documentation

- Corrected the How-To Guide's description of the GPU solver, which had the API backwards:
  the convenience methods take `maxiter` as a keyword, have no `forced_fixed_iters`
  argument, and *are* the fixed-iteration path.
- Fixed a dangling reference to a nonexistent `liquid_fraction(param_set, T)` method, the
  liquid-fraction ramp interval (`[T_freeze - 0.2 K, T_freeze]`, not ±0.1 K around
  freezing), and the `has_condensate` threshold (`q_min`, not `eps`).
- Replaced the heat-capacity table in the Mathematical Formulation with values computed from
  `ClimaParams` at build time; the hard-coded ones had drifted by up to 2%. Corrected the
  specific-enthalpy example, which had quoted `c_p T` rather than a value referenced to `T_0`.
- Documented all six saturation-adjustment formulations; the README, home page, and How-To
  Guide had listed four.
- `julia` compat raised to `1.10` in `Project.toml`, `docs/Project.toml`, and
  `test/Project.toml`, matching what CI actually tests.

v1.2.2
--------

- Bumped `RootSolvers` compatibility.
- Added developer guides and updated/corrected documentation.

v1.2.1
--------

- ![][badge-🐛bugfix] Fixed type stability in `liquid_fraction`.
  PR [320](https://github.com/CliMA/Thermodynamics.jl/pull/320)

v1.2.0
--------

- Update `has_condensate` to accept thermodynamic parameters as input.
  PR [316](https://github.com/CliMA/Thermodynamics.jl/pull/316)

v1.1.0
--------

- Use a fixed `q_min` parameter as a lower bound for condensate specific humidity in the function `has_condensate`.
  PR [315](https://github.com/CliMA/Thermodynamics.jl/pull/315)

v1.0.1
--------

- ![][badge-🐛bugfix] Fixed limiting behavior of saturation vapor pressure at very low temperatures.
  PR [313](https://github.com/CliMA/Thermodynamics.jl/pull/313)

v1.0.0
--------

- ![][badge-💥breaking] Removal of the deprecated object-oriented and state-based API (e.g., `ThermodynamicState`, `PhasePartition`, `PhaseEquil`). The package now relies exclusively on a stateless, completely functional API.
- ![][badge-🚀performance] Enforced zero-allocations across core routines (`air_temperature`, `air_pressure`, `saturation_adjustment` with fixed iterations) backed by explicit testing.
- ![][badge-✨feature/enhancement] Upgraded testing infrastructure: `Documenter.doctest` continuously checks all mathematical examples in the docstrings.
- ![][badge-✨feature/enhancement] Rewrite of documentation
  - Restructured documentation for better organization and clarity
  - Updated all function documentation to be consistent and comprehensive
  - Improved code examples and usage patterns throughout documentation
  - Enhanced cross-references and internal documentation links
  - Split monolithic src/relations.jl and test/relations.jl into multiple files

- ![][badge-✨feature/enhancement] Added new thermodynamic functions for export
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

- ![][badge-✨feature/enhancement] Allow cp_m calculation without using PhasePartition.
  PR [254](https://github.com/CliMA/Thermodynamics.jl/pull/254)

v0.12.10
-------

- ![][badge-🐛bugfix] Asynchronous printing on the gpu has been fixed.
  PR [239](https://github.com/CliMA/Thermodynamics.jl/pull/239)

v0.12.9
-------

- ![][badge-🐛bugfix] Protest against zero division in relative humidity
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
[badge-💥breaking]: https://img.shields.io/badge/💥breaking-red.svg
[badge-🚀performance]: https://img.shields.io/badge/🚀performance-green.svg
