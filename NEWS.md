Thermodynamics.jl Release Notes
========================

main
--------

v1.3.0
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
  `maxiter`, rather than a limit cycle that no iteration count could escape.

- The `converged` field returned by the full signature is now computed from the final Newton
  increment rather than hard-coded to `true`.

### Internals

- ![][badge-🚀performance] The saturation derivatives evaluated the saturation vapor pressure
  several times per call. Each now computes it once and passes it along, so
  `∂e_int_∂T_sat_ρ` costs one `exp`/two `log` instead of three/four, and `∂θ_li_∂T_sat_ρ` two
  and three instead of four and six. `vapor_pressure_deficit` selects the phase's parameters
  rather than evaluating both phases. At the default `maxiter = 2`, a `ρe` CPU solve is about 11% faster and a `ρθ_li` solve about 22% faster than before this release, with
  `∂e_int_∂T_sat_ρ` about 42% faster; the numerical results are unchanged.

- ![][badge-🐛bugfix] `saturation_vapor_pressure_calc` and `latent_heat_generic` each carried
  a fallback method that promoted its arguments to a common type "to allow AD with dual
  numbers". Those fallbacks were unreachable, so the promotion never happened. For
  `saturation_vapor_pressure_calc` this was a defect: differentiating with respect to
  `LH_0` or `Δcp` left the zero-temperature guard comparing a float against a dual, and the
  function inferred as `Any`. It now promotes for real and infers concretely. The
  `latent_heat_generic` fallback was removed: that function is unbranched arithmetic and
  promotes on its own.

- ![][badge-🐛bugfix] `solution_type()` had no call sites (`_saturation_adjustment_generic`
  hard-coded `RS.CompactSolution()`), so the `DataCollection` workflow documented in that
  module could never collect anything and always reported zeros. The solver now calls
  `solution_type()`.

- `Base.broadcastable` is now defined for `IndepVars`, which is the most frequently broadcast
  of the dispatch singletons, so it no longer needs a `Ref` wrapper at broadcast call sites.

- The six `saturation_adjustment` methods, and the six convenience methods, were near
  duplicates of each other. Each formulation is now described by four small dispatch methods
  (`_temperature_unsaturated`, `_temperature_all_ice`, `_q_vap_sat_at`,
  `_equilibrium_partition`) collected in one place so the six can be read side by side, and
  the solver itself is written once. Similarly, the internal-energy and enthalpy derivatives
  now share `_∂energy_∂T_sat`, and the two liquid-ice potential temperature derivatives share
  `_θ_li_derivative_state`. This removes roughly 700 lines of duplicated code.

### Robustness

- ![][badge-🐛bugfix] `liquid_fraction_ramp` threw a `DomainError` below `T_icenuc` for any
  non-integer `pow_icenuc`, an unrecoverable failure inside a GPU kernel. The ramp argument
  is now clamped before exponentiation.
- ![][badge-🐛bugfix] `liquid_fraction` and `q_vap_saturation_from_pressure_calc` returned
  `Union{Float32, Float64}` when called with arguments wider than the parameter set, which
  de-optimizes GPU kernels. Both branches now share a type.
- Temperature guards are applied consistently across all root-finding closures; previously
  only the `ρe` non-Newton path was guarded.

### Testing

- Added direct tests for exported functions missing before:
  `soundspeed_air`, both `supersaturation` methods, `vol_vapor_mixing_ratio`, `exner`,
  `potential_temperature`, `virtual_pottemp` and `saturation_vapor_pressure_mixture`.
- Added property tests that constrain the physics rather than restating the implementation:
  Clausius-Clapeyron for both pure phases at several temperatures and across the
  mixed-phase band, `p_ice^* < p_liq^*` below the triple point, monotonicity of the
  saturation curves and of `q_vap_saturation` in temperature and density, bounds on the
  liquid fraction and relative humidity, and the dry limit of the mixture properties.
- Added round-trip tests for `enthalpy`, `internal_energy`, `pρ` and `pθ_li`, and quantified
  the truncation error of `air_temperature(::ρθ_li)`, which is a second-order Taylor
  approximation rather than an exact inverse.
- Added regression tests pinning this release's fixes: exact recovery of equilibrium states
  by all six formulations, monotone convergence of the fixed-iteration solver with
  `maxiter`, the accuracy envelope of the default iteration count, the saturation-humidity
  derivative against automatic differentiation, safety of non-integer `pow_icenuc`, AD type
  stability of `saturation_vapor_pressure_calc`, and that `solution_type` reaches the solver.
- Extended type-stability and allocation coverage from the `ρe` formulation alone to all
  six, plus the phase partitioning and the six analytic derivatives the solvers call every
  iteration.
- Every `RootSolvers` method is now checked for accuracy, not merely for returning a finite
  number; `T_guess` (good, poor, and absent) and varying `maxiter` are exercised.
- `exceptions.jl` and `liquid_fraction_ramp_tests.jl` now run in both precisions, and the
  suite has its first `@test_throws` coverage.
- Removed the `Documenter.doctest` call from the test suite and the corresponding test
  dependency. There are no `jldoctest` blocks in the package, so it asserted nothing.
- ![][badge-🐛bugfix] `test/runtests_gpu.jl` fell back to `Array` whenever CUDA was
  unavailable, so a CI agent with a broken GPU produced a green run labelled "GPU tests"
  while exercising no kernel at all. Reaching a GPU is now required: `CuArray` fails if CUDA
  is unusable, an argument-less run fails unless `THERMODYNAMICS_ALLOW_CPU_FALLBACK=true`,
  and the Buildkite step passes `CuArray` explicitly. `Array` still selects the CPU
  deliberately.
- Broadened the device-broadcast tests from three functions to the full saturation
  adjustment surface: the saturation functions, phase partitioning, latent heats, entropy,
  relative humidity, all six analytic derivatives that the fixed-iteration solver evaluates
  each iteration, and `saturation_adjustment` itself for all six formulations. These now run
  in both `Float32` and `Float64`; previously the GPU suite was `Float32`-only, so promotion
  bugs on device were invisible.

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
- `relative_humidity` now warns that it uses a different liquid-fraction parameterization
  from `q_vap_saturation` and the saturation adjustment solvers, so that the two can differ
  appreciably in condensate-free air below freezing. The divergence is intentional; the
  docstring says which to use when consistency with the solver matters.
- `entropy_dry` and `entropy_vapor` now state that their reference pressure is `MSLP`, not
  the `p_ref_theta` used by `exner` and the potential temperatures.
- `bound_upper_temperature` documented the case where its two requirements conflict and the
  returned bound exceeds `T_max`, rather than claiming it never does.
- `internal_energy_sat` and `enthalpy_sat` added to the API reference; they are the physical
  core of saturation adjustment but were absent while their derivatives were listed.
- Added a "Citing" section to the README and a `CITATION.cff`, so the paper the package
  implements can be cited from GitHub's own interface.
- Documented the `Parameters` accessors in the API reference. The derived ones
  (`Rv_over_Rd`, `LH_f0`, `e_int_v0`, `e_int_i0`, `kappa_d`, `cv_d`, `cv_v`, `cv_l`, `cv_i`)
  are used by the package's own examples and by downstream models, but appeared nowhere in
  the docs; `checkdocs = :exports` does not reach into the submodule.
- Fixed the last Documenter warnings: two inline math spans in `Formulation.md` began with a
  bare identifier, which Julia's Markdown parser read as an interpolation. The docs now
  build warning-free apart from a size hint on the API page.
- Corrected the `TemperatureProfiles` usage description, which claimed the constructor takes
  an altitude. It takes a parameter set; the resulting object is callable and takes the
  altitude.

### Repository

- Populated `.github/pull_request_template.md`, which was an empty file, with a checklist
  covering formatting, docstrings, tests, `NEWS.md`, and the extra steps that apply when a
  change affects numerical results.
- Fixed a broken link in `AGENTS.md` (`software_design_patterns.md` is under `code-quality/`,
  not `architecture/`) and removed references to `perf/jet.jl` from `perf/README.md`; that
  file does not exist, and its JET checks live in `test/optimization_tests.jl`.
- Removed unused documentation dependencies (`CairoMakie`, `ExprTools`, `JLD2`,
  `KernelAbstractions`, `Literate`) and the orphaned `docs/plot_helpers.jl`, and added the
  missing `Plots` compat entry. `CairoMakie` had been pinned to `0.11`, constraining
  resolution for a package the docs never loaded.
- Deleted the stale `docs_output_api_check.txt` and `docs_output_final.txt` build artifacts
  from the repository root and gitignored the pattern.
- Aligned the copyright years in `LICENSE` (2026) with `NOTICE` (2022-2026).

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
