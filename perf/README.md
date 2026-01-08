## perf/ (optional developer tooling)

This directory contains **optional** scripts for performance and optimization work. It is **not**
used by CI.

### What’s here

- `perf/microbenchmarks.jl`: microbenchmarks for the functional `saturation_adjustment` API
- `perf/jet.jl`: JET optimization checks for representative `saturation_adjustment` call paths
- `perf/common*.jl`: shared helpers (inputs based on `test/TestedProfiles.jl`)

### How to run

Use a dedicated environment so benchmark/JET dependencies don’t affect the main test env:

```bash
julia --project=perf -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=perf perf/microbenchmarks.jl
julia --project=perf perf/jet.jl
```


