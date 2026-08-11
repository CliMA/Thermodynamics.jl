# Saturation adjustment functions for various combinations of input variables

import RootSolvers as RS

export saturation_adjustment
export ∂e_int_∂T_sat_ρ
export ∂e_int_∂T_sat_p
export ∂h_∂T_sat_p
export ∂θ_li_∂T_sat_ρ
export ∂θ_li_∂T_sat_p
export ∂p_∂T_sat_ρ

# ---------------------------------------------
# Public API: full solver methods
# ---------------------------------------------

# ---------------------------------------------
# Per-formulation definitions
# ---------------------------------------------

# Each `IndepVars` formulation is defined by four small methods, collected here so that the
# six formulations can be read side by side. Everything else about saturation adjustment is
# shared. `var₁` and `var₂` are the formulation's two independent variables, in the order
# the type name spells them: `ρe` is `(ρ, e_int)`, `pθ_li` is `(p, θ_li)`, and so on.

"""
    _temperature_unsaturated(param_set, indep_vars, var₁, var₂, q_tot)

Internal function. Temperature the state would have if all of its water were vapor.

This is the exact answer whenever the state turns out to be unsaturated, and the starting
point for the iteration when it is not.
"""
@inline _temperature_unsaturated(param_set::APS, ::ρe, ρ, e_int, q_tot) =
    air_temperature(param_set, e_int, q_tot)

@inline _temperature_unsaturated(param_set::APS, ::pe, p, e_int, q_tot) =
    air_temperature(param_set, e_int, q_tot)

@inline _temperature_unsaturated(param_set::APS, ::ph, p, h, q_tot) =
    air_temperature(param_set, ph(), h, q_tot, zero(q_tot), zero(q_tot))

@inline _temperature_unsaturated(param_set::APS, ::pρ, p, ρ, q_tot) =
    air_temperature(param_set, pρ(), p, ρ, q_tot)

@inline _temperature_unsaturated(param_set::APS, ::pθ_li, p, θ_li, q_tot) =
    air_temperature(param_set, pθ_li(), p, θ_li, q_tot)

@inline _temperature_unsaturated(param_set::APS, ::ρθ_li, ρ, θ_li, q_tot) =
    air_temperature(param_set, ρθ_li(), ρ, θ_li, q_tot)

"""
    _temperature_all_ice(param_set, indep_vars, var₁, var₂, q_tot)

Internal function. Temperature the state would have if all of its water were ice.

Freezing all water releases the most latent heat the state can supply, so this
bounds the solution from above and is used to bracket the root.
"""
@inline _temperature_all_ice(param_set::APS, ::ρe, ρ, e_int, q_tot) =
    air_temperature(param_set, ρe(), e_int, q_tot, zero(q_tot), q_tot)

@inline _temperature_all_ice(param_set::APS, ::pe, p, e_int, q_tot) =
    air_temperature(param_set, pe(), e_int, q_tot, zero(q_tot), q_tot)

@inline _temperature_all_ice(param_set::APS, ::ph, p, h, q_tot) =
    air_temperature(param_set, ph(), h, q_tot, zero(q_tot), q_tot)

@inline _temperature_all_ice(param_set::APS, ::pρ, p, ρ, q_tot) =
    air_temperature(param_set, pρ(), p, ρ, q_tot, zero(q_tot), q_tot)

@inline _temperature_all_ice(param_set::APS, ::pθ_li, p, θ_li, q_tot) =
    air_temperature(param_set, pθ_li(), p, θ_li, q_tot, zero(q_tot), q_tot)

@inline _temperature_all_ice(param_set::APS, ::ρθ_li, ρ, θ_li, q_tot) =
    air_temperature(param_set, ρθ_li(), ρ, θ_li, q_tot, zero(q_tot), q_tot)

"""
    _q_vap_sat_at(param_set, indep_vars, var₁, var₂, T, q_tot)

Internal function. Saturation specific humidity at temperature `T` for this formulation.

Formulations that carry a density use it directly. Those that carry only a pressure go
through [`q_vap_saturation_from_pressure`](@ref) rather than forming a density first, which
would require the phase partition that is still being solved for.
"""
@inline _q_vap_sat_at(param_set::APS, ::ρe, ρ, e_int, T, q_tot) =
    q_vap_saturation(param_set, T, ρ)

@inline _q_vap_sat_at(param_set::APS, ::pe, p, e_int, T, q_tot) =
    q_vap_saturation_from_pressure(param_set, q_tot, p, T)

@inline _q_vap_sat_at(param_set::APS, ::ph, p, h, T, q_tot) =
    q_vap_saturation_from_pressure(param_set, q_tot, p, T)

@inline _q_vap_sat_at(param_set::APS, ::pρ, p, ρ, T, q_tot) =
    q_vap_saturation(param_set, T, ρ)

@inline _q_vap_sat_at(param_set::APS, ::pθ_li, p, θ_li, T, q_tot) =
    q_vap_saturation_from_pressure(param_set, q_tot, p, T)

@inline _q_vap_sat_at(param_set::APS, ::ρθ_li, ρ, θ_li, T, q_tot) =
    q_vap_saturation(param_set, T, ρ)

"""
    _equilibrium_partition(param_set, indep_vars, var₁, var₂, T, q_tot)

Internal function. Equilibrium `(q_liq, q_ice)` at the solved temperature.

Mirrors [`_q_vap_sat_at`](@ref) in its choice of density or pressure as the second
independent variable.
"""
@inline _equilibrium_partition(param_set::APS, ::ρe, ρ, e_int, T, q_tot) =
    condensate_partition(param_set, T, ρ, q_tot)

@inline _equilibrium_partition(param_set::APS, ::pe, p, e_int, T, q_tot) =
    _condensate_partition_from_p(param_set, T, p, q_tot)

@inline _equilibrium_partition(param_set::APS, ::ph, p, h, T, q_tot) =
    _condensate_partition_from_p(param_set, T, p, q_tot)

@inline _equilibrium_partition(param_set::APS, ::pρ, p, ρ, T, q_tot) =
    condensate_partition(param_set, T, ρ, q_tot)

@inline _equilibrium_partition(param_set::APS, ::pθ_li, p, θ_li, T, q_tot) =
    _condensate_partition_from_p(param_set, T, p, q_tot)

@inline _equilibrium_partition(param_set::APS, ::ρθ_li, ρ, θ_li, T, q_tot) =
    condensate_partition(param_set, T, ρ, q_tot)

# ---------------------------------------------
# Public API: full solver methods
# ---------------------------------------------

"""
    saturation_adjustment(
        ::Type{M},  # RS.RootSolvingMethod type
        param_set,
        ::ρe,    ρ, e_int, q_tot, maxiter, tol, [T_guess], [forced_fixed_iters]
    )
    saturation_adjustment(
        ::Type{M}, param_set,
        ::pe,    p, e_int, q_tot, maxiter, tol, [T_guess], [forced_fixed_iters]
    )
    saturation_adjustment(
        ::Type{M}, param_set,
        ::ph,    p, h,     q_tot, maxiter, tol, [T_guess], [forced_fixed_iters]
    )
    saturation_adjustment(
        ::Type{M}, param_set,
        ::pρ,    p, ρ,     q_tot, maxiter, tol, [T_guess], [forced_fixed_iters]
    )
    saturation_adjustment(
        ::Type{M}, param_set,
        ::pθ_li, p, θ_li,  q_tot, maxiter, tol, [T_guess], [forced_fixed_iters]
    )
    saturation_adjustment(
        ::Type{M}, param_set,
        ::ρθ_li, ρ, θ_li,  q_tot, maxiter, tol, [T_guess], [forced_fixed_iters]
    )

Solve for the temperature at which the given state is in phase equilibrium, and return that
temperature together with the resulting condensate.

# Arguments
 - `M`: root-solving method from `RootSolvers`, e.g. `RS.NewtonsMethod` or `RS.SecantMethod`.
   `RS.NewtonsMethod` is recommended: analytic derivatives are available for every
   formulation.
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `indep_vars`: an [`IndepVars`](@ref) singleton naming the two independent variables
 - `var₁`, `var₂`: those two variables, in the order the type name spells them. Units follow
   the variable: `ρ` [kg/m³], `p` [Pa], `e_int` and `h` [J/kg], `θ_li` [K].
 - `q_tot`: total specific humidity [kg/kg]
 - `maxiter`: maximum number of solver iterations
 - `tol`: relative tolerance on the temperature, or an `RS.RelativeSolutionTolerance`
 - `T_guess`: optional initial temperature guess [K]. Defaults to `nothing`.
 - `forced_fixed_iters`: run `maxiter` iterations without testing convergence, avoiding
   branch divergence on GPUs. `T_guess` and `tol` are ignored when `true`. Defaults to
   `false`.

# Returns
 - `NamedTuple` `(; T, q_liq, q_ice, converged)`:
     - `T`: temperature [K]
     - `q_liq`: liquid specific humidity [kg/kg]
     - `q_ice`: ice specific humidity [kg/kg]
     - `converged`: whether the solver converged

# Notes
 - The state is first tested against the temperature it would have with no condensate; if
   that temperature is already subsaturated it is returned directly, since it solves the
   problem exactly. Otherwise the equilibrium condition is solved by root-finding and the
   condensate follows from [`condensate_partition`](@ref).
 - **GPU broadcasting**: pass `forced_fixed_iters` as a positional `Bool`.

# Examples
```julia
import RootSolvers as RS
import Thermodynamics as TD
using ClimaParams

param_set = TD.Parameters.ThermodynamicsParameters(Float64)

# A cloudy state near 290 K and 850 hPa. Internal energy is measured relative to the
# reference temperature T_0, so it is negative here.
ρ, e_int, q_tot = 1.0175, -30923.0, 0.019
sol = TD.saturation_adjustment(
    RS.NewtonsMethod, param_set, TD.ρe(), ρ, e_int, q_tot, 20, 1e-4,
)
sol.T, sol.q_liq, sol.q_ice, sol.converged   # ≈ (290.0, 0.0049, 0.0, true)
```

See also the convenience methods below, which pick GPU-friendly defaults for you.
"""
function saturation_adjustment(
    ::Type{M},  # RS.AbstractMethod type
    param_set::APS,
    indep_vars::IndepVars,
    var₁,
    var₂,
    q_tot,
    maxiter::Int,
    tol,
    T_guess = nothing,
    forced_fixed_iters::Bool = false,
) where {M}
    if forced_fixed_iters
        return saturation_adjustment_fixed_iters(
            param_set,
            indep_vars,
            var₁,
            var₂,
            q_tot,
            maxiter,
        )
    end

    (T, converged) = _saturation_adjustment_generic(
        M,
        param_set,
        indep_vars,
        var₁,
        var₂,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )

    (q_liq, q_ice) =
        _equilibrium_partition(param_set, indep_vars, var₁, var₂, T, q_tot)
    return (; T, q_liq, q_ice, converged)
end

"""
    saturation_adjustment(param_set, ::ρe,    ρ, e_int, q_tot; maxiter = 2)
    saturation_adjustment(param_set, ::pe,    p, e_int, q_tot; maxiter = 2)
    saturation_adjustment(param_set, ::ph,    p, h,     q_tot; maxiter = 2)
    saturation_adjustment(param_set, ::pρ,    p, ρ,     q_tot; maxiter = 2)
    saturation_adjustment(param_set, ::pθ_li, p, θ_li,  q_tot; maxiter = 2)
    saturation_adjustment(param_set, ::ρθ_li, ρ, θ_li,  q_tot; maxiter = 2)

Convenience methods with GPU-friendly defaults.

Runs a fixed number of safeguarded Newton iterations (see
`saturation_adjustment_fixed_iters`) rather than a convergence-tested solve, so every GPU
lane executes the same instructions.

The default `maxiter = 2` is chosen for how saturation adjustment is used in a time-stepping
model, where it is applied every step to a state that was near equilibrium at the previous
one. What sets the iteration count is the gap between the unsaturated first guess and the
solution, i.e. the warming from condensing the excess vapor. Across the tested profiles that
gap never exceeds 4 K, and two iterations hold the temperature error below 2e-3 K. The error
grows for states quenched from much further out of equilibrium (about 0.1 K at a 10 K gap,
rising to a few K beyond 20 K), which is where a larger `maxiter` is worth passing.

Overshooting the iteration budget is no longer dangerous: the error is a truncation error
that decreases monotonically with `maxiter`, not a stalled iteration. Earlier versions could
settle into a limit cycle and return an answer tens of K off no matter how many iterations
were allowed.

For more control over solver parameters, use the full signature with an explicit method type.

# Returns
 - `NamedTuple` `(; T, q_liq, q_ice)`. Use the full signature if you need the `converged`
   flag.

# Examples
```julia
import Thermodynamics as TD
using ClimaParams

param_set = TD.Parameters.ThermodynamicsParameters(Float64)
sol = TD.saturation_adjustment(param_set, TD.ρe(), 1.0175, -30923.0, 0.019)
sol.T, sol.q_liq, sol.q_ice   # ≈ (290.4, 0.0046, 0.0)
```
"""
function saturation_adjustment(
    param_set::APS,
    indep_vars::IndepVars,
    var₁,
    var₂,
    q_tot;
    maxiter::Int = 2,
)
    sa_result = saturation_adjustment_fixed_iters(
        param_set,
        indep_vars,
        var₁,
        var₂,
        q_tot,
        maxiter,
    )
    return (; sa_result.T, sa_result.q_liq, sa_result.q_ice)
end

# ---------------------------------------------
# GPU-optimized: fixed iteration methods
# ---------------------------------------------

"""
    saturation_adjustment_fixed_iters(param_set, ::ThermoType, args..., maxiter)

GPU-optimized saturation adjustment using a fixed number of Newton iterations.

Bypasses bracketing and convergence testing to avoid branch divergence on GPUs, and
dispatches on thermodynamic formulation type.

# Algorithm

The iteration is Newton's method applied to the analytic continuation of the saturated
branch, i.e. the saturation excess is allowed to go negative
(see `_clamp_excess`). Both the residual and its derivative are continuous there,
which matters: the physical residual has a kink at the saturation boundary, and taking
Newton steps across it with the one-sided derivative produces a limit cycle rather than
convergence. Because the continuation does not describe subsaturated air, the result is
selected against the exact unsaturated solution (see `_is_saturated`) by a
branchless `ifelse` after the loop. Steps are limited and the iterate is kept positive by
`_newton_update`.

Note that this is *not* `RS.NewtonsMethod`, which additionally applies a backtracking line
search and would introduce data-dependent iteration counts.

# Supported formulations
- `ρe`: `saturation_adjustment_fixed_iters(param_set, ρe(), ρ, e_int, q_tot, maxiter)`
- `pe`: `saturation_adjustment_fixed_iters(param_set, pe(), p, e_int, q_tot, maxiter)`
- `ph`: `saturation_adjustment_fixed_iters(param_set, ph(), p, h, q_tot, maxiter)`
- `pθ_li`: `saturation_adjustment_fixed_iters(param_set, pθ_li(), p, θ_li, q_tot, maxiter)`
- `ρθ_li`: `saturation_adjustment_fixed_iters(param_set, ρθ_li(), ρ, θ_li, q_tot, maxiter)`
- `pρ`: `saturation_adjustment_fixed_iters(param_set, pρ(), p, ρ, q_tot, maxiter)`

# Returns
- `NamedTuple` `(; T, q_liq, q_ice, converged)`

# Notes
- `converged` reports whether the final iteration left the temperature essentially
  unchanged (see `_fixed_iters_converged`). It is a statement about the iteration
  settling, not a residual test.
- Accuracy is governed by the gap between the unsaturated first guess and the solution —
  the warming from condensing the excess vapor — rather than by the absolute humidity. Two
  iterations hold the error below 2e-3 K for gaps up to the ~4 K seen across the tested
  profiles, reaching about 0.1 K at a 10 K gap and a few K beyond 20 K. The convenience
  methods of [`saturation_adjustment`](@ref) default to `maxiter = 2`; states quenched far
  from equilibrium warrant more.
- This is an internal helper function. For the public API, use [`saturation_adjustment`](@ref)
  with `forced_fixed_iters=true` as a positional argument.
"""
function saturation_adjustment_fixed_iters end

"""
    _newton_update(param_set, T, ΔT_raw)

Internal function. Apply one safeguarded Newton increment to the temperature.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: current temperature iterate [K], strictly positive
 - `ΔT_raw`: unsafeguarded Newton increment [K]

# Returns
 - `(T_new, ΔT)`: updated temperature and the increment Newton requested [K]

Two branchless guards constrain the *search*, never the answer. The increment is limited to
`ΔT_max`, because the residual's derivative changes sharply across the saturation boundary
and an unlimited step can traverse the whole physical temperature range in one iteration.
And the iterate can lose at most half its value per step, which keeps it strictly positive —
the saturation functions are evaluable at any positive temperature but not at zero — while
still allowing descent to an arbitrarily cold solution at a geometric rate. There is no
fixed lower bound: a saturated solution at, say, 130 K is reachable, which it was not when
iterates were floored at `T_init_min`.

The returned increment is the one Newton *asked* for, not the one that survived the guards:
a guarded step reading as "settled" is exactly the false convergence signal to avoid.
"""
@inline function _newton_update(param_set::APS, T, ΔT_raw)
    FT = eltype(param_set)
    # Wide enough never to bind near a solution, small enough to keep a diverging
    # iterate within the range where the saturation functions remain informative.
    ΔT_max = FT(50)
    ΔT = clamp(ΔT_raw, -ΔT_max, ΔT_max)
    T_new = max(T + ΔT, T / 2)
    return (T_new, ΔT_raw)
end

"""
    _is_saturated(param_set, T_unsat, ρ, q_tot)

Internal function. Decide whether a state is saturated, given the temperature it would have
if all of its water were vapor.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T_unsat`: temperature obtained by assuming no condensate [K]
 - `ρ`: (moist-)air density [kg/m³]
 - `q_tot`: total specific humidity [kg/kg]

# Returns
 - `saturated`: `true` when condensate must be present

The test is exact rather than heuristic: `T_unsat` solves the no-condensate problem, so if
no condensation is required there, it is already the answer.
"""
@inline function _is_saturated(param_set::APS, T_unsat, ρ, q_tot)
    return q_tot > q_vap_saturation(param_set, T_unsat, ρ)
end

"""
    _is_saturated_from_p(param_set, T_unsat, p, q_tot)

Internal function. Pressure-based counterpart of `_is_saturated`.
"""
@inline function _is_saturated_from_p(param_set::APS, T_unsat, p, q_tot)
    return q_tot > q_vap_saturation_from_pressure(param_set, q_tot, p, T_unsat)
end

"""
    _select_solution(saturated, T_sat, ΔT_sat, T_unsat)

Internal function. Choose between the iterated saturated solution and the exact unsaturated
one, returning the corresponding temperature and final increment.

`T_unsat` must be the *unclamped* no-condensate temperature. It solves the problem exactly
when the state is subsaturated, so passing the clamped starting guess instead would return
the clamp for any state colder than it.

Both are computed unconditionally so that every GPU lane executes the same instructions; the
selection is a branchless `ifelse`.
"""
@inline function _select_solution(saturated, T_sat, ΔT_sat, T_unsat)
    T = ifelse(saturated, T_sat, T_unsat)
    ΔT = ifelse(saturated, ΔT_sat, zero(ΔT_sat))
    return (T, ΔT)
end

"""
    _fixed_iters_converged(T, ΔT)

Internal function. Report whether the fixed-iteration solver reached a self-consistent
temperature.

# Arguments
 - `T`: final temperature [K]
 - `ΔT`: increment requested by the final Newton step [K], before step limiting

# Returns
 - `converged`: `true` when the remaining error is at the level of round-off

This is a test on the iteration itself, not on a residual: the fixed-iteration solvers take
a set number of steps and never evaluate a stopping criterion.

The tolerance is a fixed relative temperature, not one derived from `eps`. Because Newton
converges quadratically here, an increment this small leaves an error far below it, so the
test is conservative: measured against the fixed point of the iteration, states reported as
converged agree to better than 3e-3 K. Tying the tolerance to `eps` instead would make the
flag mean different things in different precisions — `sqrt(eps)` is so tight in `Float64`
that solutions accurate to 1e-9 K report failure, while `eps^(1/4)` is loose enough in
`Float32` to accept errors of several K.
"""
@inline function _fixed_iters_converged(T, ΔT)
    FT = typeof(T)
    # An increment below 1e-5 T is a few mK at atmospheric temperatures, and is reachable
    # in Float32 (~80x its eps) as well as Float64.
    rtol = FT(1e-5)
    return isfinite(T) & (abs(ΔT) <= rtol * abs(T))
end

@inline function saturation_adjustment_fixed_iters(
    param_set::APS,
    ::ρe,
    ρ,
    e_int,
    q_tot,
    maxiter,
)
    T_unsat = air_temperature(param_set, e_int, q_tot)
    saturated = _is_saturated(param_set, T_unsat, ρ, q_tot)

    # Start strictly positive: the saturation functions are evaluable at any positive
    # temperature but not at zero or below. The floor is numerics, not physics, and must
    # not reach the returned value — for an unsaturated state T_unsat is exact, and a
    # saturated solution may lie at any positive temperature.
    T = max(T_unsat, T_positive_floor(eltype(param_set)))
    ΔT = zero(T)
    for _ in 1:maxiter
        e_val = internal_energy_sat(param_set, T, ρ, q_tot, Val(false))
        de_int_dT = ∂e_int_∂T_sat_ρ(param_set, T, ρ, q_tot, Val(false))
        (T, ΔT) = _newton_update(param_set, T, (e_int - e_val) / de_int_dT)
    end
    (T, ΔT) = _select_solution(saturated, T, ΔT, T_unsat)

    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)
    return (; T, q_liq, q_ice, converged = _fixed_iters_converged(T, ΔT))
end

@inline function saturation_adjustment_fixed_iters(
    param_set::APS,
    ::pe,
    p,
    e_int,
    q_tot,
    maxiter,
)
    T_unsat = air_temperature(param_set, e_int, q_tot)
    saturated = _is_saturated_from_p(param_set, T_unsat, p, q_tot)

    # Start strictly positive: the saturation functions are evaluable at any positive
    # temperature but not at zero or below. The floor is numerics, not physics, and must
    # not reach the returned value — for an unsaturated state T_unsat is exact, and a
    # saturated solution may lie at any positive temperature.
    T = max(T_unsat, T_positive_floor(eltype(param_set)))
    ΔT = zero(T)
    for _ in 1:maxiter
        (q_liq, q_ice) = _condensate_partition_from_p(param_set, T, p, q_tot, Val(false))
        e_val = internal_energy(param_set, T, q_tot, q_liq, q_ice)
        de_int_dT = ∂e_int_∂T_sat_p(param_set, T, p, q_tot, Val(false))
        (T, ΔT) = _newton_update(param_set, T, (e_int - e_val) / de_int_dT)
    end
    (T, ΔT) = _select_solution(saturated, T, ΔT, T_unsat)

    (ρ, q_liq, q_ice) = _phase_partition_from_T_p(param_set, T, p, q_tot)
    return (; T, q_liq, q_ice, converged = _fixed_iters_converged(T, ΔT))
end

@inline function saturation_adjustment_fixed_iters(
    param_set::APS,
    ::ph,
    p,
    h,
    q_tot,
    maxiter,
)
    T_unsat = air_temperature(param_set, ph(), h, q_tot, zero(q_tot), zero(q_tot))
    saturated = _is_saturated_from_p(param_set, T_unsat, p, q_tot)

    # Start strictly positive: the saturation functions are evaluable at any positive
    # temperature but not at zero or below. The floor is numerics, not physics, and must
    # not reach the returned value — for an unsaturated state T_unsat is exact, and a
    # saturated solution may lie at any positive temperature.
    T = max(T_unsat, T_positive_floor(eltype(param_set)))
    ΔT = zero(T)
    for _ in 1:maxiter
        (q_liq, q_ice) = _condensate_partition_from_p(param_set, T, p, q_tot, Val(false))
        h_val = enthalpy(param_set, T, q_tot, q_liq, q_ice)
        dh_dT = ∂h_∂T_sat_p(param_set, T, p, q_tot, Val(false))
        (T, ΔT) = _newton_update(param_set, T, (h - h_val) / dh_dT)
    end
    (T, ΔT) = _select_solution(saturated, T, ΔT, T_unsat)

    (ρ, q_liq, q_ice) = _phase_partition_from_T_p(param_set, T, p, q_tot)
    return (; T, q_liq, q_ice, converged = _fixed_iters_converged(T, ΔT))
end

@inline function saturation_adjustment_fixed_iters(
    param_set::APS,
    ::pθ_li,
    p,
    θ_li,
    q_tot,
    maxiter,
)
    T_unsat = air_temperature(param_set, pθ_li(), p, θ_li, q_tot)
    saturated = _is_saturated_from_p(param_set, T_unsat, p, q_tot)

    # Start strictly positive: the saturation functions are evaluable at any positive
    # temperature but not at zero or below. The floor is numerics, not physics, and must
    # not reach the returned value — for an unsaturated state T_unsat is exact, and a
    # saturated solution may lie at any positive temperature.
    T = max(T_unsat, T_positive_floor(eltype(param_set)))
    ΔT = zero(T)
    for _ in 1:maxiter
        (q_liq, q_ice) = _condensate_partition_from_p(param_set, T, p, q_tot, Val(false))
        θ_li_val = liquid_ice_pottemp_given_pressure(param_set, T, p, q_tot, q_liq, q_ice)
        dθ_li_dT = ∂θ_li_∂T_sat_p(param_set, T, p, q_tot, Val(false))
        (T, ΔT) = _newton_update(param_set, T, (θ_li - θ_li_val) / dθ_li_dT)
    end
    (T, ΔT) = _select_solution(saturated, T, ΔT, T_unsat)

    (ρ, q_liq, q_ice) = _phase_partition_from_T_p(param_set, T, p, q_tot)
    return (; T, q_liq, q_ice, converged = _fixed_iters_converged(T, ΔT))
end

@inline function saturation_adjustment_fixed_iters(
    param_set::APS,
    ::ρθ_li,
    ρ,
    θ_li,
    q_tot,
    maxiter,
)
    T_unsat = air_temperature(param_set, ρθ_li(), ρ, θ_li, q_tot)
    saturated = _is_saturated(param_set, T_unsat, ρ, q_tot)

    # Start strictly positive: the saturation functions are evaluable at any positive
    # temperature but not at zero or below. The floor is numerics, not physics, and must
    # not reach the returned value — for an unsaturated state T_unsat is exact, and a
    # saturated solution may lie at any positive temperature.
    T = max(T_unsat, T_positive_floor(eltype(param_set)))
    ΔT = zero(T)
    for _ in 1:maxiter
        (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot, Val(false))
        θ_li_val = liquid_ice_pottemp(param_set, T, ρ, q_tot, q_liq, q_ice)
        dθ_li_dT = ∂θ_li_∂T_sat_ρ(param_set, T, ρ, q_tot, Val(false))
        (T, ΔT) = _newton_update(param_set, T, (θ_li - θ_li_val) / dθ_li_dT)
    end
    (T, ΔT) = _select_solution(saturated, T, ΔT, T_unsat)

    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)
    return (; T, q_liq, q_ice, converged = _fixed_iters_converged(T, ΔT))
end

@inline function saturation_adjustment_fixed_iters(
    param_set::APS,
    ::pρ,
    p,
    ρ,
    q_tot,
    maxiter,
)
    T_unsat = air_temperature(param_set, pρ(), p, ρ, q_tot)
    saturated = _is_saturated(param_set, T_unsat, ρ, q_tot)

    # Start strictly positive: the saturation functions are evaluable at any positive
    # temperature but not at zero or below. The floor is numerics, not physics, and must
    # not reach the returned value — for an unsaturated state T_unsat is exact, and a
    # saturated solution may lie at any positive temperature.
    T = max(T_unsat, T_positive_floor(eltype(param_set)))
    ΔT = zero(T)
    for _ in 1:maxiter
        (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot, Val(false))
        p_val = air_pressure(param_set, T, ρ, q_tot, q_liq, q_ice)
        dp_dT = ∂p_∂T_sat_ρ(param_set, T, ρ, q_tot, Val(false))
        (T, ΔT) = _newton_update(param_set, T, (p - p_val) / dp_dT)
    end
    (T, ΔT) = _select_solution(saturated, T, ΔT, T_unsat)

    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)
    return (; T, q_liq, q_ice, converged = _fixed_iters_converged(T, ΔT))
end

# ---------------------------------------------
# Internal helpers
# ---------------------------------------------

"""
    bound_upper_temperature(param_set, T_lo, T_hi)

Internal function. Bounds the upper temperature guess `T_hi` for bracket methods.

Returns `T_hi` capped at `T_max`, then raised if necessary to keep it above `T_lo` so that
the bracket is non-degenerate.

Note that the two requirements can conflict: when `T_lo` is itself at or above `T_max`, the
separation requirement wins and the result exceeds `T_max`. That is deliberate, since a
collapsed bracket would break the solver outright, whereas a slightly too-warm upper bound
only widens the search.
"""
@inline function bound_upper_temperature(param_set, T_lo, T_hi)
    FT = eltype(param_set)
    T_max = TP.T_max(param_set)
    # Ensure T_hi is physically valid (<= T_max)
    T_hi_phys = min(T_max, T_hi)
    # Ensure T_hi > T_lo for numerical initialization (use relative tolerance)
    return max(T_lo * (1 + FT(1e-3)), T_hi_phys)
end


"""
    _make_sa_solver(::Type{M}, param_set, T_unsat, T_ice, T_guess)

Internal helper to construct a root-solving method instance for saturation adjustment.

# Arguments
- `M`: Root-solving method type (e.g., `RS.NewtonsMethod`, `RS.SecantMethod`, `RS.BrentsMethod`).
- `param_set`: Thermodynamics parameter set.
- `T_unsat`: Unsaturated temperature estimate (used for initialization or lower bound) [K].
- `T_ice`: Temperature with all water as ice (used for upper bound) [K].
- `T_guess`: Optional user-provided initial guess [K] (or `nothing`).

# Returns
- Instantiated solver method ready for `RootSolvers.find_zero`.

# Notes
- For Newton-type methods (`NewtonsMethod`, `NewtonsMethodAD`): Uses `T_guess` if provided,
  otherwise `T_unsat` (floored at `T_positive_floor`, a numerics bound far below
  any physical temperature).
- For bracket methods (`SecantMethod`, `BrentsMethod`): Constructs bracket `[T_lo, T_hi]`
  where `T_lo` is `T_unsat` with the same floor and `T_hi` is bounded by `T_ice` and
  `T_max`. `T_unsat` is a true lower bound on the saturated solution, since condensation
  can only warm the state.
"""
@inline function _make_sa_solver(
    ::Type{RS.NewtonsMethod},
    param_set::APS,
    T_unsat,
    T_ice,
    T_guess,
)
    T_init =
        T_guess isa Nothing ?
        max(T_unsat, T_positive_floor(eltype(param_set))) : T_guess
    return RS.NewtonsMethod(T_init)
end

@inline function _make_sa_solver(
    ::Type{RS.NewtonsMethodAD},
    param_set::APS,
    T_unsat,
    T_ice,
    T_guess,
)
    T_init =
        T_guess isa Nothing ?
        max(T_unsat, T_positive_floor(eltype(param_set))) : T_guess
    return RS.NewtonsMethodAD(T_init)
end

@inline function _make_sa_solver(
    ::Type{RS.SecantMethod},
    param_set::APS,
    T_unsat,
    T_ice,
    T_guess,
)
    T_floor = T_positive_floor(eltype(param_set))
    T_lo = T_guess isa Nothing ? max(T_unsat, T_floor) : max(T_guess, T_floor)
    T_hi = bound_upper_temperature(param_set, T_lo, T_ice)
    return RS.SecantMethod(T_lo, T_hi)
end

@inline function _make_sa_solver(
    ::Type{RS.BrentsMethod},
    param_set::APS,
    T_unsat,
    T_ice,
    T_guess,
)
    # BrentsMethod requires strict bracketing - ignore T_guess
    T_lo = max(T_unsat, T_positive_floor(eltype(param_set)))
    T_hi = bound_upper_temperature(param_set, T_lo, T_ice)
    return RS.BrentsMethod(T_lo, T_hi)
end


"""
    internal_energy_sat(param_set, T, ρ, q_tot)

The internal energy per unit mass in thermodynamic equilibrium at saturation.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `ρ`: (moist-)air density [kg/m³]
 - `q_tot`: total specific humidity [kg/kg]

# Returns
 - `e_int`: specific internal energy [J/kg]

The phase partition into liquid and ice is computed internally from `q_tot` using the 
temperature-dependent liquid fraction (see [`liquid_fraction_ramp`](@ref)) and saturation 
excess (see [`saturation_excess`](@ref)).
"""
@inline function internal_energy_sat(
    param_set::APS,
    T,
    ρ,
    q_tot,
    clamped::Val = Val(true),
)
    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot, clamped)
    return internal_energy(param_set, T, q_tot, q_liq, q_ice)
end

"""
    enthalpy_sat(param_set, T, ρ, q_tot)

The specific enthalpy in thermodynamic equilibrium at saturation.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `ρ`: (moist-)air density [kg/m³]
 - `q_tot`: total specific humidity [kg/kg]

# Returns
 - `h`: specific enthalpy [J/kg]

The phase partition into liquid and ice is computed internally from `q_tot` using the 
temperature-dependent liquid fraction (see [`liquid_fraction_ramp`](@ref)) and saturation 
excess (see [`saturation_excess`](@ref)).
"""
@inline function enthalpy_sat(
    param_set::APS,
    T,
    ρ,
    q_tot,
    clamped::Val = Val(true),
)
    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot, clamped)
    return enthalpy(param_set, T, q_tot, q_liq, q_ice)
end

"""
    _make_roots_function(::Type{M}, param_set, ::ThermoType, args..., q_tot)

Helper function to create the root function for Newton's method (with derivative) or
other methods (without derivative), dispatching on the method type and thermo type.

Returns `(f, f')` for `NewtonsMethod`, or just `f` for other methods.
"""
function _make_roots_function end

# ρe formulation
@inline function _make_roots_function(
    ::Type{RS.NewtonsMethod},
    param_set::APS,
    ::ρe,
    ρ,
    e_int,
    q_tot,
)
    return _T -> begin
        T_val = ReLU(_T)
        f = e_int - internal_energy_sat(param_set, T_val, ρ, q_tot)
        (f, -∂e_int_∂T_sat_ρ(param_set, T_val, ρ, q_tot))
    end
end

@inline function _make_roots_function(
    ::Type{M},
    param_set::APS,
    ::ρe,
    ρ,
    e_int,
    q_tot,
) where {M}
    return _T -> begin
        T_val = ReLU(_T)
        e_int - internal_energy_sat(param_set, T_val, ρ, q_tot)
    end
end

# pe formulation
@inline function _make_roots_function(
    ::Type{RS.NewtonsMethod},
    param_set::APS,
    ::pe,
    p,
    e_int,
    q_tot,
)
    return _T -> begin
        T_val = ReLU(_T)
        f = _internal_energy_sat_from_p(param_set, T_val, p, q_tot) - e_int
        (f, ∂e_int_∂T_sat_p(param_set, T_val, p, q_tot))
    end
end

@inline function _make_roots_function(
    ::Type{M},
    param_set::APS,
    ::pe,
    p,
    e_int,
    q_tot,
) where {M}
    return _T -> begin
        T_val = ReLU(_T)
        _internal_energy_sat_from_p(param_set, T_val, p, q_tot) - e_int
    end
end

# ph formulation
@inline function _make_roots_function(
    ::Type{RS.NewtonsMethod},
    param_set::APS,
    ::ph,
    p,
    h,
    q_tot,
)
    return _T -> begin
        T_val = ReLU(_T)
        f = _enthalpy_sat_from_p(param_set, T_val, p, q_tot) - h
        (f, ∂h_∂T_sat_p(param_set, T_val, p, q_tot))
    end
end

@inline function _make_roots_function(
    ::Type{M},
    param_set::APS,
    ::ph,
    p,
    h,
    q_tot,
) where {M}
    return _T -> begin
        T_val = ReLU(_T)
        _enthalpy_sat_from_p(param_set, T_val, p, q_tot) - h
    end
end

# pθ_li formulation
@inline function _make_roots_function(
    ::Type{RS.NewtonsMethod},
    param_set::APS,
    ::pθ_li,
    p,
    θ_li,
    q_tot,
)
    return _T -> begin
        T_val = ReLU(_T)
        (_q_liq, _q_ice) = _condensate_partition_from_p(param_set, T_val, p, q_tot)
        f =
            liquid_ice_pottemp_given_pressure(param_set, T_val, p, q_tot, _q_liq, _q_ice) - θ_li
        (f, ∂θ_li_∂T_sat_p(param_set, T_val, p, q_tot))
    end
end

@inline function _make_roots_function(
    ::Type{M},
    param_set::APS,
    ::pθ_li,
    p,
    θ_li,
    q_tot,
) where {M}
    return _T -> begin
        T_val = ReLU(_T)
        (_q_liq, _q_ice) = _condensate_partition_from_p(param_set, T_val, p, q_tot)
        liquid_ice_pottemp_given_pressure(param_set, T_val, p, q_tot, _q_liq, _q_ice) - θ_li
    end
end

# ρθ_li formulation
@inline function _make_roots_function(
    ::Type{RS.NewtonsMethod},
    param_set::APS,
    ::ρθ_li,
    ρ,
    θ_li,
    q_tot,
)
    return _T -> begin
        T_val = ReLU(_T)
        (_q_liq, _q_ice) = condensate_partition(param_set, T_val, ρ, q_tot)
        f = liquid_ice_pottemp(param_set, T_val, ρ, q_tot, _q_liq, _q_ice) - θ_li
        (f, ∂θ_li_∂T_sat_ρ(param_set, T_val, ρ, q_tot))
    end
end

@inline function _make_roots_function(
    ::Type{M},
    param_set::APS,
    ::ρθ_li,
    ρ,
    θ_li,
    q_tot,
) where {M}
    return _T -> begin
        T_val = ReLU(_T)
        (_q_liq, _q_ice) = condensate_partition(param_set, T_val, ρ, q_tot)
        liquid_ice_pottemp(param_set, T_val, ρ, q_tot, _q_liq, _q_ice) - θ_li
    end
end

# pρ formulation
@inline function _make_roots_function(
    ::Type{RS.NewtonsMethod},
    param_set::APS,
    ::pρ,
    p,
    ρ,
    q_tot,
)
    return _T -> begin
        T_val = ReLU(_T)
        (_q_liq, _q_ice) = condensate_partition(param_set, T_val, ρ, q_tot)
        f = air_pressure(param_set, T_val, ρ, q_tot, _q_liq, _q_ice) - p
        (f, ∂p_∂T_sat_ρ(param_set, T_val, ρ, q_tot))
    end
end

@inline function _make_roots_function(
    ::Type{M},
    param_set::APS,
    ::pρ,
    p,
    ρ,
    q_tot,
) where {M}
    return _T -> begin
        T_val = ReLU(_T)
        (_q_liq, _q_ice) = condensate_partition(param_set, T_val, ρ, q_tot)
        air_pressure(param_set, T_val, ρ, q_tot, _q_liq, _q_ice) - p
    end
end

"""
    _phase_partition_from_T_p(param_set, T, p, q_tot)

Helper to compute equilibrium phase partition given temperature, pressure, and total humidity.
Returns `(ρ, q_liq, q_ice)` tuple.
"""
@inline function _phase_partition_from_T_p(
    param_set::APS,
    T,
    p,
    q_tot,
    clamped::Val = Val(true),
)
    (q_liq, q_ice) = _condensate_partition_from_p(param_set, T, p, q_tot, clamped)
    ρ = air_density(param_set, T, p, q_tot, q_liq, q_ice)
    return (ρ, q_liq, q_ice)
end

"""
    _condensate_partition_from_p(param_set, T, p, q_tot, clamped = Val(true))

Internal function. Equilibrium `(q_liq, q_ice)` at a given temperature and *pressure*.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `p`: air pressure [Pa]
 - `q_tot`: total specific humidity [kg/kg]
 - `clamped`: `Val(false)` continues the saturated branch to negative condensate,
   see `_clamp_excess`

# Returns
 - `(q_liq, q_ice)`: liquid and ice specific humidities [kg/kg]

This is the pressure analogue of [`condensate_partition`](@ref). Obtaining the saturation
specific humidity from the pressure directly, via
[`q_vap_saturation_from_pressure`](@ref), keeps the result self-consistent: computing a
density first would require a phase partition that is not yet known, and assuming no
condensate there biases the saturation humidity high.
"""
@inline function _condensate_partition_from_p(
    param_set::APS,
    T,
    p,
    q_tot,
    clamped::Val = Val(true),
)
    λ = liquid_fraction_ramp(param_set, T)
    p_v_sat = saturation_vapor_pressure_mixture(param_set, T, λ)
    q_vap_sat = q_vap_saturation_from_pressure_calc(param_set, q_tot, p, p_v_sat)
    q_c = _clamp_excess(clamped, q_tot - q_vap_sat)
    return (λ * q_c, (1 - λ) * q_c)
end

"""
    _internal_energy_sat_from_p(param_set, T, p, q_tot, clamped = Val(true))

Internal function. Equilibrium internal energy at a given temperature and pressure.

Pressure analogue of [`internal_energy_sat`](@ref); see
`_condensate_partition_from_p` for why the partition is taken from the pressure
rather than from a condensate-free density.
"""
@inline function _internal_energy_sat_from_p(
    param_set::APS,
    T,
    p,
    q_tot,
    clamped::Val = Val(true),
)
    (q_liq, q_ice) = _condensate_partition_from_p(param_set, T, p, q_tot, clamped)
    return internal_energy(param_set, T, q_tot, q_liq, q_ice)
end

"""
    _enthalpy_sat_from_p(param_set, T, p, q_tot, clamped = Val(true))

Internal function. Equilibrium specific enthalpy at a given temperature and pressure.

Pressure analogue of [`enthalpy_sat`](@ref).
"""
@inline function _enthalpy_sat_from_p(
    param_set::APS,
    T,
    p,
    q_tot,
    clamped::Val = Val(true),
)
    (q_liq, q_ice) = _condensate_partition_from_p(param_set, T, p, q_tot, clamped)
    return enthalpy(param_set, T, q_tot, q_liq, q_ice)
end

"""
    _saturation_derivative_vars_p(param_set, T, p, q_tot, clamped = Val(true))

Internal function. Phase partition and its temperature derivatives at fixed *pressure*.

# Returns
A `NamedTuple` with the same fields as `_saturation_derivative_vars`, plus
`q_vap_sat` and `p_v_sat`.

Differentiating `q_v^* = ε (1 - q_tot) p_v^* / (p - p_v^*)` at fixed `p` gives

    ∂q_v^*/∂T|_p = q_v^* · (∂ln p_v^*/∂T) · p / (p - p_v^*),

where the `p / (p - p_v^*)` factor comes from the saturation vapor pressure appearing in the
denominator as well. It approaches 1 when `p_v^* ≪ p` but grows in warm, moist air.
"""
@inline function _saturation_derivative_vars_p(
    param_set::APS,
    T,
    p,
    q_tot,
    clamped::Val = Val(true),
)
    FT = eltype(param_set)
    R_v = TP.R_v(param_set)

    λ = liquid_fraction_ramp(param_set, T)
    ∂λ_∂T = ∂λ_∂T_ramp(param_set, T)
    p_v_sat = saturation_vapor_pressure_mixture(param_set, T, λ)
    q_vap_sat = q_vap_saturation_from_pressure_calc(param_set, q_tot, p, p_v_sat)

    q_c = _clamp_excess(clamped, q_tot - q_vap_sat)
    q_liq = λ * q_c
    q_ice = (1 - λ) * q_c

    # Logarithmic derivative of the mixed-phase saturation vapor pressure
    L = latent_heat_mixed(param_set, T, λ)
    ∂lnp_∂λ = log_saturation_vapor_pressure_ratio(param_set, T)
    ∂lnp_v_sat_∂T = L / (R_v * T^2) + ∂lnp_∂λ * ∂λ_∂T

    # Guard the denominator the same way q_vap_saturation_from_pressure_calc does, so that
    # p approaching p_v_sat cannot produce a division by zero inside a kernel. In that
    # regime that function returns the constant 1, so the derivative of what it actually
    # computes is zero — using 1 here would differentiate a value it never returns.
    Δp = p - p_v_sat
    amplification = ifelse(Δp ≥ ϵ_numerics(FT), p / Δp, zero(Δp))
    ∂qvs_∂T = q_vap_sat * ∂lnp_v_sat_∂T * amplification

    ∂q_liq_∂T = ∂λ_∂T * q_c + λ * (-∂qvs_∂T)
    ∂q_ice_∂T = -∂λ_∂T * q_c + (1 - λ) * (-∂qvs_∂T)

    return (;
        λ,
        q_c,
        q_liq,
        q_ice,
        ∂λ_∂T,
        ∂qvs_∂T,
        ∂q_liq_∂T,
        ∂q_ice_∂T,
        q_vap_sat,
        p_v_sat,
    )
end

"""
    _find_zero_and_convergence(
        roots_func,
        numerical_method,
        solution_type,
        tol,
        maxiter,
    )

Helper function to find the root of `roots_func` using `numerical_method`.
Returns `(root, converged)` tuple.

Used by `saturation_adjustment` functions to handle common solver logic.
"""
@inline function _find_zero_and_convergence(
    roots_func,
    numerical_method,
    solution_type,
    tol,
    maxiter,
)
    sol =
        RS.find_zero(roots_func, numerical_method, solution_type, tol, maxiter)
    DataCollection.log_meta(sol)
    return (sol.root, sol.converged)
end

"""
    _saturation_adjustment_generic(
        ::Type{Method},
        param_set,
        indep_vars,
        var₁,
        var₂,
        q_tot,
        maxiter,
        relative_temperature_tol,
        T_guess,
    )

Generic kernel for saturation adjustment. Handles unsaturated check, solver initialization,
and root-finding for all `saturation_adjustment` variants.

# Arguments
- `Method`: Root-solving method type (e.g., `RS.NewtonsMethod`, `RS.SecantMethod`).
- `param_set`: Thermodynamics parameter set.
- `indep_vars`: [`IndepVars`](@ref) singleton naming the formulation.
- `var₁`, `var₂`: the formulation's two independent variables.
- `q_tot`: Total specific humidity.
- `maxiter`: Maximum iterations for the solver.
- `relative_temperature_tol`: Relative tolerance for temperature solution.
- `T_guess`: Optional initial temperature guess.

# Returns
- `(T, converged)`: Temperature and convergence flag.

What varies between formulations is supplied by `_temperature_unsaturated`,
`_temperature_all_ice`, `_q_vap_sat_at`, and `_make_roots_function`, all of which dispatch
on `indep_vars`.
"""
@inline function _saturation_adjustment_generic(
    ::Type{Method},  # RootSolvers MethodType
    param_set::APS,
    indep_vars::IndepVars,
    var₁,
    var₂,
    q_tot,
    maxiter,
    relative_temperature_tol,
    T_guess,
) where {Method}
    _T_min = TP.T_min(param_set)
    tol =
        relative_temperature_tol isa Real ?
        RS.RelativeSolutionTolerance(relative_temperature_tol) :
        relative_temperature_tol

    # Unsaturated check: the no-condensate temperature solves the problem exactly whenever
    # the state turns out to be subsaturated there, so it is returned unmodified. Clamping
    # it would replace an exact answer with the clamp.
    T_unsat = _temperature_unsaturated(param_set, indep_vars, var₁, var₂, q_tot)

    q_v_sat = _q_vap_sat_at(param_set, indep_vars, var₁, var₂, T_unsat, q_tot)
    if q_tot <= q_v_sat
        return (T_unsat, true)
    end

    # Saturated case: solve for T
    T_ice = _temperature_all_ice(param_set, indep_vars, var₁, var₂, q_tot)
    roots_func =
        _make_roots_function(Method, param_set, indep_vars, var₁, var₂, q_tot)

    # Initialize solver (logic merged from config_sa_method.jl). The lower bracket is
    # clamped to `T_min`: that is a bound on the search, not on the answer.
    solver =
        _make_sa_solver(Method, param_set, max(_T_min, T_unsat), T_ice, T_guess)

    # `solution_type()` rather than a hard-coded `RS.CompactSolution()`, so that redefining
    # it to `RS.VerboseSolution()` actually reaches the solver and lets `DataCollection`
    # gather iteration statistics.
    (T, converged) = _find_zero_and_convergence(
        roots_func,
        solver,
        solution_type(),
        tol,
        maxiter,
    )

    return (T, converged)
end

# -------------------------------
# Derivatives for Newton's method
# -------------------------------

"""
    _saturation_derivative_vars(param_set, T, ρ, q_tot, q_vap_sat, ::Val{:ρ}, clamped)

Helper to compute common intermediate variables for saturation derivatives at fixed
density. The fixed-pressure counterpart is `_saturation_derivative_vars_p`, which takes
`p` rather than `ρ` so that the phase partition stays consistent with the pressure.

Returns a named tuple with phase partition and derivative information:
- `λ`, `q_c`, `q_liq`, `q_ice`: Phase partition variables
- `∂λ_∂T`, `∂qvs_∂T`: Temperature derivatives of liquid fraction and saturation humidity
- `∂q_liq_∂T`, `∂q_ice_∂T`: Phase partition derivatives

`clamped` must match the setting used for the residual being differentiated: with
`Val(false)` the saturation excess is allowed to go negative, so that the derivative
describes the same analytic continuation the fixed-iteration solvers iterate on.
"""
@inline function _saturation_derivative_vars(
    param_set::APS,
    T,
    ρ,
    q_tot,
    q_vap_sat,
    ::Val{:ρ},
    clamped::Val = Val(true),
)
    λ = liquid_fraction_ramp(param_set, T)
    q_c = saturation_excess(param_set, T, ρ, q_tot, clamped)
    q_liq = λ * q_c
    q_ice = (1 - λ) * q_c

    ∂λ_∂T = ∂λ_∂T_ramp(param_set, T)

    # ∂q_vap_sat/∂T at fixed ρ (includes -1/T term)
    ∂qvs_∂T = ∂q_vap_sat_∂T(param_set, T, ρ)

    # Phase partition derivatives
    ∂q_liq_∂T = ∂λ_∂T * q_c + λ * (-∂qvs_∂T)
    ∂q_ice_∂T = -∂λ_∂T * q_c + (1 - λ) * (-∂qvs_∂T)

    return (; λ, q_c, q_liq, q_ice, ∂λ_∂T, ∂qvs_∂T, ∂q_liq_∂T, ∂q_ice_∂T)
end

"""
    _select_sat_branch(unsat_branch, q_tot, q_vap_sat, x_unsat, x_sat)

Internal function. Choose between the unsaturated and saturated forms of a temperature
derivative.

With `Val(true)` (the default for the public derivative functions) this returns the exact
derivative, which is discontinuous at the saturation boundary: below it the condensate
terms are absent.

With `Val(false)` it returns the saturated form everywhere, i.e. the analytic continuation
of the saturated branch. That continuation is larger than the unsaturated derivative, so
Newton steps taken with it are damped rather than amplified when an iterate overshoots into
the unsaturated region. The fixed-iteration solvers use it to avoid the limit cycle that the
discontinuity would otherwise produce (see `saturation_adjustment_fixed_iters`).

Dispatching on `Val` keeps the choice a compile-time constant, so no branch reaches the GPU.
"""
@inline _select_sat_branch(::Val{true}, q_tot, q_vap_sat, x_unsat, x_sat) =
    ifelse(q_tot <= q_vap_sat, x_unsat, x_sat)

@inline _select_sat_branch(::Val{false}, q_tot, q_vap_sat, x_unsat, x_sat) = x_sat

"""
    ∂e_int_∂T_sat_ρ(param_set, T, ρ, q_tot, unsat_branch = Val(true))

Derivative of `internal_energy_sat` with respect to temperature at fixed density.

Uses `∂q_vap_sat/∂T|_ρ` from Clausius-Clapeyron, which includes the `-1/T` term
from the density dependence of saturation vapor pressure.

Passing `Val(false)` for `unsat_branch` returns the saturated form even below saturation;
see `_select_sat_branch`.
"""
@inline function ∂e_int_∂T_sat_ρ(
    param_set::APS,
    T,
    ρ,
    q_tot,
    unsat_branch::Val = Val(true),
)
    q_vap_sat = q_vap_saturation(param_set, T, ρ)
    vars = _saturation_derivative_vars(
        param_set,
        T,
        ρ,
        q_tot,
        q_vap_sat,
        Val(:ρ),
        unsat_branch,
    )

    (c_unsat, de_dT_sat) =
        _∂energy_∂T_sat(param_set, T, q_tot, vars, cv_m, internal_energy_vapor,
            internal_energy_liquid, internal_energy_ice)

    return _select_sat_branch(unsat_branch, q_tot, q_vap_sat, c_unsat, de_dT_sat)
end

"""
    _∂energy_∂T_sat(param_set, T, q_tot, vars, c_mixture, u_vap, u_liq, u_ice)

Internal function. Temperature derivative of a mixture energy along the saturated branch.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `q_tot`: total specific humidity [kg/kg]
 - `vars`: phase partition and its derivatives, from `_saturation_derivative_vars` or
   `_saturation_derivative_vars_p`
 - `c_mixture`: mixture heat capacity, `cv_m` or `cp_m`
 - `u_vap`, `u_liq`, `u_ice`: the corresponding component energies, e.g.
   `internal_energy_vapor` or `enthalpy_vapor`

# Returns
 - `(c_unsat, d_sat)`: the unsaturated derivative (just the heat capacity) and the saturated
   one [J/(kg·K)]

Internal energy and enthalpy differ only in which heat capacity and component energies are
used, so both are expressed here. Along the saturated branch the derivative is the mixture
heat capacity plus the energy carried by the phases whose amounts are changing,
`Σᵢ uᵢ ∂qᵢ/∂T`; below saturation nothing changes phase and only the heat capacity remains.
"""
@inline function _∂energy_∂T_sat(
    param_set::APS,
    T,
    q_tot,
    vars,
    c_mixture::F,
    u_vap::Fv,
    u_liq::Fl,
    u_ice::Fi,
) where {F, Fv, Fl, Fi}
    c_unsat = c_mixture(param_set, q_tot, zero(q_tot), zero(q_tot))
    c_sat = c_mixture(param_set, q_tot, vars.q_liq, vars.q_ice)
    d_sat =
        c_sat +
        u_vap(param_set, T) * vars.∂qvs_∂T +
        u_liq(param_set, T) * vars.∂q_liq_∂T +
        u_ice(param_set, T) * vars.∂q_ice_∂T
    return (c_unsat, d_sat)
end

"""
    ∂e_int_∂T_sat_p(param_set, T, p, q_tot)

Derivative of `internal_energy_sat` with respect to temperature at fixed pressure.

Uses `∂q_vap_sat/∂T|_p = q_vap_sat * L / (R_v T²)` (Clausius-Clapeyron at constant
pressure), which differs from the constant-density form used in [`∂e_int_∂T_sat_ρ`](@ref)
by an additional `q_vap_sat / T` term.
"""
@inline function ∂e_int_∂T_sat_p(
    param_set::APS,
    T,
    p,
    q_tot,
    unsat_branch::Val = Val(true),
)
    vars = _saturation_derivative_vars_p(param_set, T, p, q_tot, unsat_branch)

    (c_unsat, de_dT_sat) =
        _∂energy_∂T_sat(param_set, T, q_tot, vars, cv_m, internal_energy_vapor,
            internal_energy_liquid, internal_energy_ice)

    return _select_sat_branch(unsat_branch, q_tot, vars.q_vap_sat, c_unsat, de_dT_sat)
end

"""
    ∂h_∂T_sat_p(param_set, T, p, q_tot)

Derivative of `enthalpy_sat` with respect to temperature at fixed pressure.

Structured identically to [`∂e_int_∂T_sat_p`](@ref) but with component enthalpies
(`cp_m` instead of `cv_m`, `enthalpy_vapor` instead of `internal_energy_vapor`, etc.).
"""
@inline function ∂h_∂T_sat_p(
    param_set::APS,
    T,
    p,
    q_tot,
    unsat_branch::Val = Val(true),
)
    vars = _saturation_derivative_vars_p(param_set, T, p, q_tot, unsat_branch)

    (c_unsat, dh_dT_sat) = _∂energy_∂T_sat(param_set, T, q_tot, vars, cp_m,
        enthalpy_vapor, enthalpy_liquid, enthalpy_ice)

    return _select_sat_branch(unsat_branch, q_tot, vars.q_vap_sat, c_unsat, dh_dT_sat)
end

"""
    ∂θ_li_∂T_sat_p(param_set, T, p, q_tot)

Derivative of `liquid_ice_pottemp_given_pressure` (at saturation equilibrium)
with respect to temperature at fixed pressure.

Uses the product rule on `θ_li = θ * (1 - L_c / (cp_m T))`, differentiating
each factor through the T-dependent phase partition.
"""
@inline function ∂θ_li_∂T_sat_p(
    param_set::APS,
    T,
    p,
    q_tot,
    unsat_branch::Val = Val(true),
)
    vars = _saturation_derivative_vars_p(param_set, T, p, q_tot, unsat_branch)
    st = _θ_li_derivative_state(param_set, T, p, q_tot, vars)

    # At fixed p only α varies with T: ∂θ/∂T = θ (1/T + ln(p₀/p) ∂α/∂T)
    ∂θ_∂T = st.θ * (1 / T - st.ln_p_over_p0 * st.∂α_∂T)
    dθ_li_dT_sat = ∂θ_∂T * st.F + st.θ * st.∂F_∂T

    # Unsaturated: α is constant, L_c = 0, F = 1, so ∂θ_li/∂T = θ/T
    dθ_li_dT_unsat = st.θ / T

    return _select_sat_branch(
        unsat_branch,
        q_tot,
        vars.q_vap_sat,
        dθ_li_dT_unsat,
        dθ_li_dT_sat,
    )
end

"""
    _θ_li_derivative_state(param_set, T, p, q_tot, vars)

Internal function. State and partial derivatives shared by the two liquid-ice potential
temperature derivatives.

# Arguments
 - `param_set`: thermodynamics parameter set, see [`Thermodynamics`](@ref)
 - `T`: temperature [K]
 - `p`: air pressure [Pa]; at fixed density the caller passes `ρ R_m T`
 - `q_tot`: total specific humidity [kg/kg]
 - `vars`: phase partition and its derivatives

# Returns
A `NamedTuple` with `θ`, `α`, `F`, `ln_p_over_p0`, `∂α_∂T`, `∂R_m_∂T`, `∂cp_m_∂T`, `∂F_∂T`
and `cp_m`.

Everything here is common to the fixed-pressure and fixed-density cases; what differs is
only `∂θ/∂T`, since at fixed density `p` itself varies with temperature. `ln_p_over_p0` is
returned because both the Exner function and `∂θ/∂T` need it, and computing it once avoids
a second logarithm.
"""
@inline function _θ_li_derivative_state(param_set::APS, T, p, q_tot, vars)
    R_v = TP.R_v(param_set)
    cp_v = TP.cp_v(param_set)
    cp_l = TP.cp_l(param_set)
    cp_i = TP.cp_i(param_set)
    LH_v0 = TP.LH_v0(param_set)
    LH_s0 = TP.LH_s0(param_set)
    p0 = TP.p_ref_theta(param_set)

    R_m = gas_constant_air(param_set, q_tot, vars.q_liq, vars.q_ice)
    _cp_m = cp_m(param_set, q_tot, vars.q_liq, vars.q_ice)
    α = R_m / _cp_m

    # Π = (p/p₀)^α; keeping the logarithm lets ∂θ/∂T reuse it below.
    ln_p_over_p0 = log(p / p0)
    Π = exp(α * ln_p_over_p0)
    θ = T / Π

    L_c = humidity_weighted_latent_heat(param_set, vars.q_liq, vars.q_ice)
    F = 1 - L_c / (_cp_m * T)

    ∂R_m_∂T = R_v * vars.∂qvs_∂T
    ∂cp_m_∂T = (cp_l - cp_v) * vars.∂q_liq_∂T + (cp_i - cp_v) * vars.∂q_ice_∂T

    # ∂α/∂T = (∂R_m cp_m - R_m ∂cp_m) / cp_m²
    ∂α_∂T = (∂R_m_∂T * _cp_m - R_m * ∂cp_m_∂T) / _cp_m^2

    ∂L_c_∂T = LH_v0 * vars.∂q_liq_∂T + LH_s0 * vars.∂q_ice_∂T

    # ∂F/∂T = -1/(cp_m T) (∂L_c/∂T - L_c (1/T + ∂cp_m/∂T / cp_m))
    ∂F_∂T = -1 / (_cp_m * T) * (∂L_c_∂T - L_c * (1 / T + ∂cp_m_∂T / _cp_m))

    return (; θ, α, F, ln_p_over_p0, ∂α_∂T, ∂R_m_∂T, ∂cp_m_∂T, ∂F_∂T, cp_m = _cp_m)
end

"""
    ∂θ_li_∂T_sat_ρ(param_set, T, ρ, q_tot)

Derivative of `liquid_ice_pottemp` (at saturation equilibrium)
with respect to temperature at fixed density.

Structured like [`∂θ_li_∂T_sat_p`](@ref), but accounts for `p = ρ R_m T` varying
with `T`. This introduces two corrections to `∂θ/∂T`:
  - an extra `-α/T` from `∂ln(p)/∂T|_ρ = 1/T + ∂R_m/∂T / R_m`, and
  - an extra `-∂R_m/∂T / cp_m` from the same.

The `∂q_vap_sat/∂T` also differs: at fixed ρ it carries an extra `-q_vap_sat / T`
relative to the fixed-p Clausius–Clapeyron form.
"""
@inline function ∂θ_li_∂T_sat_ρ(
    param_set::APS,
    T,
    ρ,
    q_tot,
    unsat_branch::Val = Val(true),
)
    q_vap_sat = q_vap_saturation(param_set, T, ρ)
    vars = _saturation_derivative_vars(
        param_set,
        T,
        ρ,
        q_tot,
        q_vap_sat,
        Val(:ρ),
        unsat_branch,
    )

    R_m = gas_constant_air(param_set, q_tot, vars.q_liq, vars.q_ice)
    p = ρ * R_m * T  # ideal gas at fixed ρ
    st = _θ_li_derivative_state(param_set, T, p, q_tot, vars)

    # At fixed ρ, p varies with T as well, contributing -α/T and -∂R_m/∂T / cp_m
    # relative to the fixed-pressure case.
    ∂θ_∂T =
        st.θ * (
            (1 - st.α) / T - st.ln_p_over_p0 * st.∂α_∂T -
            st.∂R_m_∂T / st.cp_m
        )
    dθ_li_dT_sat = ∂θ_∂T * st.F + st.θ * st.∂F_∂T

    # Unsaturated: q_liq = q_ice = 0, L_c = 0, F = 1, R_m constant
    # θ = (p0/(ρ R_m))^α · T^(1-α)  =>  dθ/dT = (1-α)·θ/T
    dθ_li_dT_unsat = st.θ * (1 - st.α) / T

    return _select_sat_branch(
        unsat_branch,
        q_tot,
        q_vap_sat,
        dθ_li_dT_unsat,
        dθ_li_dT_sat,
    )
end

"""
    ∂p_∂T_sat_ρ(param_set, T, ρ, q_tot)

Derivative of `air_pressure` (at saturation equilibrium) with respect to
temperature at fixed density.

From `p = ρ R_m T` and `∂R_m/∂T = R_v · ∂q_vap_sat/∂T`:

    ∂p/∂T|_ρ = ρ · (R_m + T · ∂R_m/∂T)

Unsaturated: `R_m` is constant, so `∂p/∂T = ρ R_m`.
"""
@inline function ∂p_∂T_sat_ρ(
    param_set::APS,
    T,
    ρ,
    q_tot,
    unsat_branch::Val = Val(true),
)
    R_v = TP.R_v(param_set)
    q_vap_sat = q_vap_saturation(param_set, T, ρ)

    # Unsaturated: R_m constant
    R_m_unsat = gas_constant_air(param_set, q_tot, zero(q_tot), zero(q_tot))
    dp_dT_unsat = ρ * R_m_unsat

    vars = _saturation_derivative_vars(
        param_set,
        T,
        ρ,
        q_tot,
        q_vap_sat,
        Val(:ρ),
        unsat_branch,
    )
    R_m = gas_constant_air(param_set, q_tot, vars.q_liq, vars.q_ice)

    # ∂R_m/∂T
    ∂R_m_∂T = R_v * vars.∂qvs_∂T

    # p = ρ R_m T  =>  ∂p/∂T = ρ(R_m + T · ∂R_m/∂T)
    dp_dT_sat = ρ * (R_m + T * ∂R_m_∂T)

    return _select_sat_branch(unsat_branch, q_tot, q_vap_sat, dp_dT_unsat, dp_dT_sat)
end
