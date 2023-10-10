# Saturation adjustment input space convergence maps

The saturation adjustment procedure requires solving a non-linear
equation.

In the [Tested Profiles](@ref) section, we plotted the tested thermodynamic
states. In this section, we explore the convergence of the
input space beyond what is tested. In particular, rather than
being interested in physically meaningful combinations of
constructor inputs (e.g., `ρ, e_int, q_tot`), we are interested
in all permutations of inputs within a given range of `ρ, e_int,
q_tot`. Some of these permutations may not be physically meaningful,
or likely to be observed in climate simulations, but showing
the convergence space helps illustrate the buffer between our
tested profiles and the nearest space where convergence fails.

This section is dedicated to monitoring the status and improvement
of the performance and robustness of various numerical methods
in solving the saturation adjustment equations for various thermodynamic
formulations.

!!! note

    `dims` in `docs/src/saturation_adjustment.jl` is currently set to
    ``dims = (6, 6, 6);`` to avoid heavy computations in the doc build, but you
    may want to increase it to, e.g., ``dims = (10, 10, 10);`` when running locally
    to see a higher resolution map.

```@example
include("saturation_adjustment.jl")
```

## 3D space
| Numerical method  | Converged  |  Non-converged |
:-----------------:|:-----------------:|:---------------------:
SecantMethod | ![](3DSpace_converged_SecantMethod.svg)       |  ![](3DSpace_non_converged_SecantMethod.svg)
NewtonsMethod | ![](3DSpace_converged_NewtonsMethod.svg)      |  ![](3DSpace_non_converged_NewtonsMethod.svg)

## 2D slices, binned by total specific humidity

| Numerical method  | Converged  |  Non-converged |
:-----------------:|:-----------------:|:---------------------:
SecantMethod | ![](2DSlice_converged_SecantMethod.svg)  |  ![](2DSlice_non_converged_SecantMethod.svg)
NewtonsMethod | ![](2DSlice_converged_NewtonsMethod.svg)  |  ![](2DSlice_non_converged_NewtonsMethod.svg)
