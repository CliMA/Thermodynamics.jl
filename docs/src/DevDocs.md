# Input space exploration

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

```@example
include(joinpath(@__DIR__, "..", "ThreeDimensionalInput.jl"))
```

## Converged cases (3D view)
![](Scatter3DConverged.svg)

## Non-converged cases (3D view)
![](Scatter3DNonConverged.svg)

## Converged cases (2D view), binned by total specific humidity
![](Slices2DConverged.svg)

## Non-converged cases (2D view), binned by total specific humidity
![](Slices2DNonConverged.svg)
