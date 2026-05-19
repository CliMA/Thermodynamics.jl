# Thermodynamics.jl Agent Guide

## Ecosystem Guidelines

Please refer to the shared CliMA agent index for ecosystem-wide rules regarding architecture, performance, code quality, infrastructure, and workflows:

- [docs/dev-guides/AGENTS.md](docs/dev-guides/AGENTS.md) — Shared CliMA agent guidelines.

> Shared guides live at `docs/dev-guides/` and are vendored from the canonical source:
> <https://github.com/CliMA/DeveloperGuides>. Edit shared guides there, not here.

## Before You Act: Agent Autonomy

Before making changes that are externally visible or scientifically consequential (`git push`, version bumps, reproducibility-test edits, CI config changes, public API renames), check [docs/dev-guides/workflow/agent_autonomy.md](docs/dev-guides/workflow/agent_autonomy.md). The boundaries listed there require explicit user approval.

## Repo-Specific Guidelines

Thermodynamics.jl provides a library of thermodynamic functions for the CliMA ecosystem. It is a relatively small, focused package.

### Architecture

- **Pure-function library**: All public API functions are pure (no mutation of global state). They accept a parameter set (`AbstractThermodynamicsParameters`) and thermodynamic state variables as arguments.
- **Saturation adjustment**: The core numerical solver lives in `src/saturation_adjustment.jl`. Multiple methods are available (Newton, Newton-AD, secant, Brent's). This is the most performance-critical code path.

### Source layout

| Path | Purpose |
|------|---------|
| `src/Thermodynamics.jl` | Module definition, exports |
| `src/Parameters.jl` | Parameter interface (`AbstractThermodynamicsParameters`) |
| `src/ThermoTypes.jl` | Independent-variable dispatch types (`IndepVars`) and phase types |
| `src/saturation_adjustment.jl` | Saturation adjustment solvers |
| `src/air_*.jl` | Thermodynamic property functions |
| `src/TemperatureProfiles.jl` | Reference temperature/pressure profiles |
| `src/DataCollection.jl` | Internal solver diagnostic statistics |
| `src/aux_functions.jl` | Internal helpers (`ReLU`, `fast_power`, `ϵ_numerics`) |
| `ext/CreateParametersExt.jl` | ClimaParams weak-dep extension for parameter construction |
| `test/` | Test suite |
| `perf/` | Performance benchmarks |
| `docs/` | Documentation source |

## Local norms

- For package tests, prefer `Pkg.test()` over manually `include`ing `test/runtests.jl` because test-only deps are loaded through the package test path.
- Match existing style: explicit names, narrow imports, comments that explain why.
- Follow the software design patterns in [docs/dev-guides/architecture/software_design_patterns.md](docs/dev-guides/architecture/software_design_patterns.md) for new code and refactor toward them when touching existing code.
- Run `julia -e 'using JuliaFormatter; format(".")'` before committing code.

## Self-correction

- If the source layout table above is discovered to be stale, update it.
- If the user gives a correction about how work should be done in this repo, add it to `Local norms` or another clearly labeled persistent section in this file so future sessions inherit it.
