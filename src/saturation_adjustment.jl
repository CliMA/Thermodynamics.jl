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

"""
    saturation_adjustment(
        ::Type{M},  # RS.RootSolvingMethod type
        param_set,
        ::ρe,
        ρ::Real,
        e_int::Real,
        q_tot::Real,
        maxiter::Int,
        tol,
        [T_guess::Union{Nothing, Real} = nothing],
        [forced_fixed_iters::Bool = false]
    )

Compute the saturation equilibrium temperature `T` and phase partition `(q_liq, q_ice)`
given density `ρ`, internal energy `e_int`, and total specific humidity `q_tot`.

# Arguments
- `M`: Root-solving method type from RS.jl. Use `RS.NewtonsMethod`,
  `RS.SecantMethod`, `RS.BrentsMethod`, or `RS.NewtonsMethodAD`.
- `param_set`: Thermodynamics parameter set, see [`Thermodynamics`](@ref).
- `ρ`: Density of moist air [kg/m³].
- `e_int`: Specific internal energy [J/kg].
- `q_tot`: Total specific humidity [kg/kg].
- `maxiter`: Maximum iterations for the solver [dimensionless integer].
- `tol`: Relative tolerance for the temperature solution (or a `RS.RelativeSolutionTolerance`).
- `T_guess`: Optional initial guess for the temperature [K]. Defaults to `nothing`.
- `forced_fixed_iters`: Optional boolean to force a fixed number of iterations (`maxiter`)
  without checking for convergence. Useful for GPU optimization to avoid branch divergence.
  When `true`, `T_guess` and `tol` are ignored. Defaults to `false`.

# Returns
- `NamedTuple` `(; T, q_liq, q_ice, converged)`:
    - `T`: Temperature [K]
    - `q_liq`: Liquid specific humidity [kg/kg]
    - `q_ice`: Ice specific humidity [kg/kg]
    - `converged`: Boolean flag indicating if the solver converged

# Notes
- This function solves for `T` such that `e_int = internal_energy_sat(param_set, T, ρ, q_tot)` using
  root-finding, then computes `(q_liq, q_ice)` from [`condensate_partition`](@ref).
- `NewtonsMethod` is recommended for all formulations (analytic derivatives available).
- **GPU broadcasting**: Pass `forced_fixed_iters` as a positional Bool.
"""
function saturation_adjustment(
    ::Type{M},  # RS.AbstractMethod type
    param_set::APS,
    ::ρe,
    ρ,
    e_int,
    q_tot,
    maxiter::Int,
    tol,
    T_guess = nothing,
    forced_fixed_iters::Bool = false,
) where {M}
    if forced_fixed_iters
        return saturation_adjustment_fixed_iters(param_set, ρe(), ρ, e_int, q_tot, maxiter)
    end
    temp_unsat_func = (param_set, e_int, q_tot) -> air_temperature(param_set, e_int, q_tot)

    q_sat_unsat_func = (param_set, T, q_tot) -> q_vap_saturation(param_set, T, ρ)

    # For T_ice upper bound: all water is ice
    temp_ice_func =
        (param_set, e_int, q_tot) ->
            air_temperature(param_set, ρe(), e_int, q_tot, zero(q_tot), q_tot)

    # Root function (supports analytic derivative if M=NewtonsMethod)
    roots_func = _make_roots_function(M, param_set, ρe(), ρ, e_int, q_tot)

    (T, converged) = _saturation_adjustment_generic(
        M,
        param_set,
        e_int,
        q_tot,
        maxiter,
        tol,
        T_guess,
        temp_unsat_func,
        q_sat_unsat_func,
        temp_ice_func,
        roots_func,
    )

    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)
    return (; T, q_liq, q_ice, converged)
end

"""
    saturation_adjustment(
        ::Type{M},  # RS.RootSolvingMethod type
        param_set,
        ::pe,
        p::Real,
        e_int::Real,
        q_tot::Real,
        maxiter::Int,
        tol,
        [T_guess::Union{Nothing, Real} = nothing],
        [forced_fixed_iters::Bool = false]
    )

Compute the saturation equilibrium temperature `T` and phase partition `(q_liq, q_ice)`
given pressure `p`, specific internal energy `e_int`, and total specific humidity `q_tot`.

# Arguments
- `M`: Root-solving method type from `RootSolvers.jl`. Supported types:
  `RS.SecantMethod`, `RS.BrentsMethod`, `RS.NewtonsMethod`, `RS.NewtonsMethodAD`.
- `param_set`: Thermodynamics parameter set, see [`Thermodynamics`](@ref).
- `p`: Pressure of moist air [Pa].
- `e_int`: Specific internal energy [J/kg].
- `q_tot`: Total specific humidity [kg/kg].
- `maxiter`: Maximum iterations for the solver [dimensionless integer].
- `tol`: Relative tolerance for the temperature solution (or a `RS.RelativeSolutionTolerance`).
- `T_guess`: Optional initial guess for the temperature [K]. Defaults to `nothing`.
- `forced_fixed_iters`: Optional boolean to force a fixed number of iterations (`maxiter`)
  without checking for convergence. Useful for GPU optimization to avoid branch divergence.
  When `true`, `T_guess` and `tol` are ignored. Defaults to `false`.

# Returns
- `NamedTuple` `(; T, q_liq, q_ice, converged)`:
    - `T`: Temperature [K]
    - `q_liq`: Liquid specific humidity [kg/kg]
    - `q_ice`: Ice specific humidity [kg/kg]
    - `converged`: Boolean flag indicating if the solver converged

# Notes
- **GPU broadcasting**: Pass `forced_fixed_iters` as a positional Bool.
"""
function saturation_adjustment(
    ::Type{M},  # RS.AbstractMethod type
    param_set::APS,
    ::pe,
    p,
    e_int,
    q_tot,
    maxiter::Int,
    tol,
    T_guess = nothing,
    forced_fixed_iters::Bool = false,
) where {M}
    if forced_fixed_iters
        return saturation_adjustment_fixed_iters(param_set, pe(), p, e_int, q_tot, maxiter)
    end
    q_sat_func =
        (param_set, T, q_tot) ->
            q_vap_saturation_from_pressure(param_set, q_tot, p, T)

    temp_ice_func =
        (param_set, e_int, q_tot) ->
            air_temperature(param_set, e_int, q_tot, zero(q_tot), q_tot)

    temp_unsat_func =
        (param_set, e_int, q_tot) -> air_temperature(param_set, e_int, q_tot)

    roots_func = _make_roots_function(M, param_set, pe(), p, e_int, q_tot)

    (T, converged) = _saturation_adjustment_generic(
        M,
        param_set,
        e_int,
        q_tot,
        maxiter,
        tol,
        T_guess,
        temp_unsat_func,
        q_sat_func,
        temp_ice_func,
        roots_func,
    )

    (ρ, q_liq, q_ice) = _phase_partition_from_T_p(param_set, T, p, q_tot)
    return (; T, q_liq, q_ice, converged)
end

"""
    saturation_adjustment(
        ::Type{M},  # RS.RootSolvingMethod type
        param_set,
        ::ph,
        p::Real,
        h::Real,
        q_tot::Real,
        maxiter::Int,
        tol,
        [T_guess::Union{Nothing, Real} = nothing],
        [forced_fixed_iters::Bool = false]
    )

Compute the saturation equilibrium temperature `T` and phase partition `(q_liq, q_ice)`
given pressure `p`, specific enthalpy `h`, and total specific humidity `q_tot`.

Returns a `NamedTuple` `(; T, q_liq, q_ice, converged)`.

# Arguments
- `M`: Root-solving method type from `RootSolvers.jl`. Supported types:
  `RS.SecantMethod`, `RS.BrentsMethod`, `RS.NewtonsMethod`, `RS.NewtonsMethodAD`.
- `param_set`: Thermodynamics parameter set, see [`Thermodynamics`](@ref).
- `p`: Pressure of moist air [Pa].
- `h`: Specific enthalpy [J/kg].
- `q_tot`: Total specific humidity [kg/kg].
- `maxiter`: Maximum iterations for the solver [dimensionless integer].
- `tol`: Relative tolerance for the temperature solution (or a `RS.RelativeSolutionTolerance`).
- `T_guess`: Optional initial guess for the temperature [K]. Defaults to `nothing`.
- `forced_fixed_iters`: Optional boolean to force a fixed number of iterations (`maxiter`)
  without checking for convergence. Useful for GPU optimization to avoid branch divergence.
  When `true`, `T_guess` and `tol` are ignored. Defaults to `false`.

# Returns
- `NamedTuple` `(; T, q_liq, q_ice, converged)`:
    - `T`: Temperature [K]
    - `q_liq`: Liquid specific humidity [kg/kg]
    - `q_ice`: Ice specific humidity [kg/kg]
    - `converged`: Boolean flag indicating if the solver converged

# Notes
- **GPU broadcasting**: Pass `forced_fixed_iters` as a positional Bool.
"""
function saturation_adjustment(
    ::Type{M},  # RS.AbstractMethod type
    param_set::APS,
    ::ph,
    p,
    h,
    q_tot,
    maxiter::Int,
    tol,
    T_guess = nothing,
    forced_fixed_iters::Bool = false,
) where {M}
    if forced_fixed_iters
        return saturation_adjustment_fixed_iters(param_set, ph(), p, h, q_tot, maxiter)
    end
    q_sat_func =
        (param_set, T, q_tot) ->
            q_vap_saturation_from_pressure(param_set, q_tot, p, T)

    temp_unsat_func =
        (param_set, h, q_tot) ->
            air_temperature(param_set, ph(), h, q_tot, 0, 0)

    temp_ice_func =
        (param_set, h, q_tot) ->
            air_temperature(param_set, ph(), h, q_tot, zero(q_tot), q_tot)

    roots_func = _make_roots_function(M, param_set, ph(), p, h, q_tot)

    (T, converged) = _saturation_adjustment_generic(
        M,
        param_set,
        h,
        q_tot,
        maxiter,
        tol,
        T_guess,
        temp_unsat_func,
        q_sat_func,
        temp_ice_func,
        roots_func,
    )

    (ρ, q_liq, q_ice) = _phase_partition_from_T_p(param_set, T, p, q_tot)
    return (; T, q_liq, q_ice, converged)
end

"""
    saturation_adjustment(
        ::Type{M},  # RS.RootSolvingMethod type
        param_set,
        ::pθ_li,
        p::Real,
        θ_li::Real,
        q_tot::Real,
        maxiter::Int,
        tol,
        [T_guess::Union{Nothing, Real} = nothing],
        [forced_fixed_iters::Bool = false]
    )

Compute the saturation equilibrium temperature `T` and phase partition `(q_liq, q_ice)`
given pressure `p`, liquid-ice potential temperature `θ_li`, and total specific humidity `q_tot`.

# Arguments
- `M`: Root-solving method type from `RootSolvers.jl`. Supported types:
  `RS.SecantMethod`, `RS.BrentsMethod`, `RS.NewtonsMethod`, `RS.NewtonsMethodAD`.
- `param_set`: Thermodynamics parameter set, see [`Thermodynamics`](@ref).
- `p`: Pressure of moist air [Pa].
- `θ_li`: Liquid-ice potential temperature [K].
- `q_tot`: Total specific humidity [kg/kg].
- `maxiter`: Maximum iterations for the solver [dimensionless integer].
- `tol`: Relative tolerance for the temperature solution (or a `RS.RelativeSolutionTolerance`).
- `T_guess`: Optional initial guess for the temperature [K]. Defaults to `nothing`.
- `forced_fixed_iters`: Optional boolean to force a fixed number of iterations (`maxiter`)
  without checking for convergence. Useful for GPU optimization to avoid branch divergence.
  When `true`, `T_guess` and `tol` are ignored. Defaults to `false`.

# Returns
- `NamedTuple` `(; T, q_liq, q_ice, converged)`:
    - `T`: Temperature [K]
    - `q_liq`: Liquid specific humidity [kg/kg]
    - `q_ice`: Ice specific humidity [kg/kg]
    - `converged`: Boolean flag indicating if the solver converged

# Notes
- **GPU broadcasting**: Pass `forced_fixed_iters` as a positional Bool.
"""
function saturation_adjustment(
    ::Type{M},  # RS.AbstractMethod type
    param_set::APS,
    ::pθ_li,
    p,
    θ_li,
    q_tot,
    maxiter::Int,
    tol,
    T_guess = nothing,
    forced_fixed_iters::Bool = false,
) where {M}
    if forced_fixed_iters
        return saturation_adjustment_fixed_iters(
            param_set,
            pθ_li(),
            p,
            θ_li,
            q_tot,
            maxiter,
        )
    end
    temp_unsat_func =
        (param_set, θ_li, q_tot) ->
            air_temperature(param_set, pθ_li(), p, θ_li, q_tot)

    q_sat_func =
        (param_set, T, q_tot) ->
            q_vap_saturation_from_pressure(param_set, q_tot, p, T)

    temp_ice_func =
        (param_set, θ_li, q_tot) ->
            air_temperature(
                param_set,
                pθ_li(),
                p,
                θ_li,
                q_tot,
                zero(q_tot),
                q_tot,
            )

    roots_func = _make_roots_function(M, param_set, pθ_li(), p, θ_li, q_tot)

    (T, converged) = _saturation_adjustment_generic(
        M,
        param_set,
        θ_li,
        q_tot,
        maxiter,
        tol,
        T_guess,
        temp_unsat_func,
        q_sat_func,
        temp_ice_func,
        roots_func,
    )

    # Compute equilibrium phase partition
    (ρ, q_liq, q_ice) = _phase_partition_from_T_p(param_set, T, p, q_tot)
    return (; T, q_liq, q_ice, converged)
end

"""
    saturation_adjustment(
        ::Type{M},  # RS.RootSolvingMethod type
        param_set,
        ::ρθ_li,
        ρ::Real,
        θ_li::Real,
        q_tot::Real,
        maxiter::Int,
        tol,
        [T_guess::Union{Nothing, Real} = nothing],
        [forced_fixed_iters::Bool = false]
    )

Compute the saturation equilibrium temperature `T` and phase partition `(q_liq, q_ice)`
given density `ρ`, liquid-ice potential temperature `θ_li`, and total specific humidity `q_tot`.

# Arguments
- `M`: Root-solving method type from `RootSolvers.jl`. Supported types:
  `RS.SecantMethod`, `RS.BrentsMethod`, `RS.NewtonsMethod`, `RS.NewtonsMethodAD`.
- `param_set`: Thermodynamics parameter set, see [`Thermodynamics`](@ref).
- `ρ`: Density of moist air [kg/m³].
- `θ_li`: Liquid-ice potential temperature [K].
- `q_tot`: Total specific humidity [kg/kg].
- `maxiter`: Maximum iterations for the solver [dimensionless integer].
- `tol`: Relative tolerance for the temperature solution (or a `RS.RelativeSolutionTolerance`).
- `T_guess`: Optional initial guess for the temperature [K]. Defaults to `nothing`.
- `forced_fixed_iters`: Optional boolean to force a fixed number of iterations (`maxiter`)
  without checking for convergence. Useful for GPU optimization to avoid branch divergence.
  When `true`, `T_guess` and `tol` are ignored. Defaults to `false`.

# Returns
- `NamedTuple` `(; T, q_liq, q_ice, converged)`:
    - `T`: Temperature [K]
    - `q_liq`: Liquid specific humidity [kg/kg]
    - `q_ice`: Ice specific humidity [kg/kg]
    - `converged`: Boolean flag indicating if the solver converged

# Notes
- **GPU broadcasting**: Pass `forced_fixed_iters` as a positional Bool.
"""
function saturation_adjustment(
    ::Type{M},  # RS.AbstractMethod type
    param_set::APS,
    ::ρθ_li,
    ρ,
    θ_li,
    q_tot,
    maxiter::Int,
    tol,
    T_guess = nothing,
    forced_fixed_iters::Bool = false,
) where {M}
    if forced_fixed_iters
        return saturation_adjustment_fixed_iters(
            param_set,
            ρθ_li(),
            ρ,
            θ_li,
            q_tot,
            maxiter,
        )
    end
    temp_unsat_func =
        (param_set, θ_li, q_tot) ->
            air_temperature(param_set, ρθ_li(), ρ, θ_li, q_tot)

    q_sat_func = (param_set, T, q_tot) ->
        q_vap_saturation(param_set, T, ρ)

    temp_ice_func =
        (param_set, θ_li, q_tot) ->
            air_temperature(
                param_set,
                ρθ_li(),
                ρ,
                θ_li,
                q_tot,
                zero(q_tot),
                q_tot,
            )

    roots_func = _make_roots_function(M, param_set, ρθ_li(), ρ, θ_li, q_tot)

    (T, converged) = _saturation_adjustment_generic(
        M,
        param_set,
        θ_li,
        q_tot,
        maxiter,
        tol,
        T_guess,
        temp_unsat_func,
        q_sat_func,
        temp_ice_func,
        roots_func,
    )

    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)
    return (; T, q_liq, q_ice, converged)
end

"""
    saturation_adjustment(
        ::Type{M},  # RS.RootSolvingMethod type
        param_set,
        ::pρ,
        p::Real,
        ρ::Real,
        q_tot::Real,
        maxiter::Int,
        tol,
        [T_guess::Union{Nothing, Real} = nothing],
        [forced_fixed_iters::Bool = false]
    )

Compute the saturation equilibrium temperature `T` and phase partition `(q_liq, q_ice)`
given pressure `p`, density `ρ`, and total specific humidity `q_tot`.

# Arguments
- `M`: Root-solving method type from `RootSolvers.jl`. Supported types:
  `RS.SecantMethod`, `RS.BrentsMethod`, `RS.NewtonsMethod`, `RS.NewtonsMethodAD`.
- `param_set`: Thermodynamics parameter set, see [`Thermodynamics`](@ref).
- `p`: Pressure of moist air [Pa].
- `ρ`: Density of moist air [kg/m³].
- `q_tot`: Total specific humidity [kg/kg].
- `maxiter`: Maximum iterations for the solver [dimensionless integer].
- `tol`: Relative tolerance for the temperature solution (or a `RS.RelativeSolutionTolerance`).
- `T_guess`: Optional initial guess for the temperature [K]. Defaults to `nothing`.
- `forced_fixed_iters`: Optional boolean to force a fixed number of iterations (`maxiter`)
  without checking for convergence. Useful for GPU optimization to avoid branch divergence.
  When `true`, `T_guess` and `tol` are ignored. Defaults to `false`.

# Returns
- `NamedTuple` `(; T, q_liq, q_ice, converged)`:
    - `T`: Temperature [K]
    - `q_liq`: Liquid specific humidity [kg/kg]
    - `q_ice`: Ice specific humidity [kg/kg]
    - `converged`: Boolean flag indicating if the solver converged

# Notes
- **GPU broadcasting**: Pass `forced_fixed_iters` as a positional Bool.
"""
function saturation_adjustment(
    ::Type{M},  # RS.AbstractMethod type
    param_set::APS,
    ::pρ,
    p,
    ρ,
    q_tot,
    maxiter::Int,
    tol,
    T_guess = nothing,
    forced_fixed_iters::Bool = false,
) where {M}
    if forced_fixed_iters
        return saturation_adjustment_fixed_iters(param_set, pρ(), p, ρ, q_tot, maxiter)
    end
    temp_unsat_func =
        (param_set, p, q_tot) ->
            air_temperature(param_set, pρ(), p, ρ, q_tot)

    q_sat_func = (param_set, T, q_tot) ->
        q_vap_saturation(param_set, T, ρ)

    temp_ice_func =
        (param_set, p, q_tot) ->
            air_temperature(param_set, pρ(), p, ρ, q_tot, zero(q_tot), q_tot)

    roots_func = _make_roots_function(M, param_set, pρ(), p, ρ, q_tot)

    (T, converged) = _saturation_adjustment_generic(
        M,
        param_set,
        p, # thermo_var (target p)
        q_tot,
        maxiter,
        tol,
        T_guess,
        temp_unsat_func,
        q_sat_func,
        temp_ice_func,
        roots_func,
    )

    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)
    return (; T, q_liq, q_ice, converged)
end

# ---------------------------------------------
# Convenience methods with reasonable defaults
# ---------------------------------------------

"""
    saturation_adjustment(
        param_set,
        ::ρe,
        ρ,
        e_int,
        q_tot;
        maxiter::Int = 2,
    )

Convenience method for `ρe` formulation with reasonable GPU-optimized defaults.

Runs a fixed number of safeguarded Newton iterations (see
`saturation_adjustment_fixed_iters`) rather than a convergence-tested solve, so
every GPU lane executes the same instructions.

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

For more control over solver parameters, use the full signature with explicit method type.

# Returns
- `NamedTuple` `(; T, q_liq, q_ice)`. Use the full signature if you need the `converged`
  flag.
"""
function saturation_adjustment(
    param_set::APS,
    ::ρe,
    ρ,
    e_int,
    q_tot;
    maxiter::Int = 2,
)
    sa_result = saturation_adjustment_fixed_iters(param_set, ρe(), ρ, e_int, q_tot, maxiter)
    return (; sa_result.T, sa_result.q_liq, sa_result.q_ice)
end

"""
    saturation_adjustment(
        param_set,
        ::pe,
        p,
        e_int,
        q_tot;
        maxiter::Int = 2,
    )

Convenience method for `pe` formulation with reasonable GPU-optimized defaults.

Runs a fixed number of safeguarded Newton iterations (see
`saturation_adjustment_fixed_iters`) rather than a convergence-tested solve, so
every GPU lane executes the same instructions.

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

For more control over solver parameters, use the full signature with explicit method type.

# Returns
- `NamedTuple` `(; T, q_liq, q_ice)`. Use the full signature if you need the `converged`
  flag.
"""
function saturation_adjustment(
    param_set::APS,
    ::pe,
    p,
    e_int,
    q_tot;
    maxiter::Int = 2,
)
    sa_result = saturation_adjustment_fixed_iters(param_set, pe(), p, e_int, q_tot, maxiter)
    return (; sa_result.T, sa_result.q_liq, sa_result.q_ice)
end

"""
    saturation_adjustment(
        param_set,
        ::ph,
        p,
        h,
        q_tot;
        maxiter::Int = 2,
    )

Convenience method for `ph` formulation with reasonable GPU-optimized defaults.

Runs a fixed number of safeguarded Newton iterations (see
`saturation_adjustment_fixed_iters`) rather than a convergence-tested solve, so
every GPU lane executes the same instructions.

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

For more control over solver parameters, use the full signature with explicit method type.

# Returns
- `NamedTuple` `(; T, q_liq, q_ice)`. Use the full signature if you need the `converged`
  flag.
"""
function saturation_adjustment(
    param_set::APS,
    ::ph,
    p,
    h,
    q_tot;
    maxiter::Int = 2,
)
    sa_result = saturation_adjustment_fixed_iters(param_set, ph(), p, h, q_tot, maxiter)
    return (; sa_result.T, sa_result.q_liq, sa_result.q_ice)
end

"""
    saturation_adjustment(
        param_set,
        ::pθ_li,
        p,
        θ_li,
        q_tot;
        maxiter::Int = 2,
    )

Convenience method for `pθ_li` formulation with reasonable GPU-optimized defaults.

Runs a fixed number of safeguarded Newton iterations (see
`saturation_adjustment_fixed_iters`) rather than a convergence-tested solve, so
every GPU lane executes the same instructions.

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

For more control over solver parameters, use the full signature with explicit method type.

# Returns
- `NamedTuple` `(; T, q_liq, q_ice)`. Use the full signature if you need the `converged`
  flag.
"""
function saturation_adjustment(
    param_set::APS,
    ::pθ_li,
    p,
    θ_li,
    q_tot;
    maxiter::Int = 2,
)
    sa_result =
        saturation_adjustment_fixed_iters(param_set, pθ_li(), p, θ_li, q_tot, maxiter)
    return (; sa_result.T, sa_result.q_liq, sa_result.q_ice)
end

"""
    saturation_adjustment(
        param_set,
        ::ρθ_li,
        ρ,
        θ_li,
        q_tot;
        maxiter::Int = 2,
    )

Convenience method for `ρθ_li` formulation with reasonable GPU-optimized defaults.

Runs a fixed number of safeguarded Newton iterations (see
`saturation_adjustment_fixed_iters`) rather than a convergence-tested solve, so
every GPU lane executes the same instructions.

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

For more control over solver parameters, use the full signature with explicit method type.

# Returns
- `NamedTuple` `(; T, q_liq, q_ice)`. Use the full signature if you need the `converged`
  flag.
"""
function saturation_adjustment(
    param_set::APS,
    ::ρθ_li,
    ρ,
    θ_li,
    q_tot;
    maxiter::Int = 2,
)
    sa_result =
        saturation_adjustment_fixed_iters(param_set, ρθ_li(), ρ, θ_li, q_tot, maxiter)
    return (; sa_result.T, sa_result.q_liq, sa_result.q_ice)
end

"""
    saturation_adjustment(
        param_set,
        ::pρ,
        p,
        ρ,
        q_tot;
        maxiter::Int = 2,
    )

Convenience method for `pρ` formulation with reasonable GPU-optimized defaults.

Runs a fixed number of safeguarded Newton iterations (see
`saturation_adjustment_fixed_iters`) rather than a convergence-tested solve, so
every GPU lane executes the same instructions.

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

For more control over solver parameters, use the full signature with explicit method type.

# Returns
- `NamedTuple` `(; T, q_liq, q_ice)`. Use the full signature if you need the `converged`
  flag.
"""
function saturation_adjustment(
    param_set::APS,
    ::pρ,
    p,
    ρ,
    q_tot;
    maxiter::Int = 2,
)
    sa_result = saturation_adjustment_fixed_iters(param_set, pρ(), p, ρ, q_tot, maxiter)
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
 - `T`: current temperature iterate [K]
 - `ΔT_raw`: unsafeguarded Newton increment [K]

# Returns
 - `(T_new, ΔT)`: updated temperature and the increment actually applied [K]

The increment is limited to `ΔT_max` and the result is kept at or above `T_init_min`. Both
guards are branchless. Limiting matters because the residual's derivative changes sharply
across the saturation boundary, so an unlimited step can traverse the whole physical
temperature range in one iteration and land somewhere with no useful information; keeping
the iterate positive additionally prevents `log` of a negative temperature downstream.
"""
@inline function _newton_update(param_set::APS, T, ΔT_raw)
    FT = eltype(param_set)
    T_init_min = TP.T_init_min(param_set)
    # Wide enough never to bind near a solution, small enough to keep a diverging
    # iterate within the range where the saturation functions remain informative.
    ΔT_max = FT(50)
    ΔT = clamp(ΔT_raw, -ΔT_max, ΔT_max)
    T_new = max(T_init_min, T + ΔT)
    return (T_new, T_new - T)
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
 - `ΔT`: increment applied by the final iteration [K]

# Returns
 - `converged`: `true` when the last increment left the temperature essentially unchanged

This is a test on the iteration itself, not on a residual: the fixed-iteration solvers take
a set number of steps and never evaluate a stopping criterion. A small final increment means
the iteration has settled; a large one means `maxiter` ran out before it did.
"""
@inline function _fixed_iters_converged(T, ΔT)
    rtol = sqrt(eps(typeof(T)))
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
    T_init_min = TP.T_init_min(param_set)
    T_0 = max(T_init_min, T_unsat)
    saturated = _is_saturated(param_set, T_0, ρ, q_tot)

    T = T_0
    ΔT = zero(T)
    for _ in 1:maxiter
        e_val = internal_energy_sat(param_set, T, ρ, q_tot, Val(false))
        de_int_dT = ∂e_int_∂T_sat_ρ(param_set, T, ρ, q_tot, Val(false))
        (T, ΔT) = _newton_update(param_set, T, (e_int - e_val) / de_int_dT)
    end
    (T, ΔT) = _select_solution(saturated, T, ΔT, T_0)

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
    T_init_min = TP.T_init_min(param_set)
    T_0 = max(T_init_min, T_unsat)
    saturated = _is_saturated_from_p(param_set, T_0, p, q_tot)

    T = T_0
    ΔT = zero(T)
    for _ in 1:maxiter
        (q_liq, q_ice) = _condensate_partition_from_p(param_set, T, p, q_tot, Val(false))
        e_val = internal_energy(param_set, T, q_tot, q_liq, q_ice)
        de_int_dT = ∂e_int_∂T_sat_p(param_set, T, p, q_tot, Val(false))
        (T, ΔT) = _newton_update(param_set, T, (e_int - e_val) / de_int_dT)
    end
    (T, ΔT) = _select_solution(saturated, T, ΔT, T_0)

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
    T_init_min = TP.T_init_min(param_set)
    T_0 = max(T_init_min, T_unsat)
    saturated = _is_saturated_from_p(param_set, T_0, p, q_tot)

    T = T_0
    ΔT = zero(T)
    for _ in 1:maxiter
        (q_liq, q_ice) = _condensate_partition_from_p(param_set, T, p, q_tot, Val(false))
        h_val = enthalpy(param_set, T, q_tot, q_liq, q_ice)
        dh_dT = ∂h_∂T_sat_p(param_set, T, p, q_tot, Val(false))
        (T, ΔT) = _newton_update(param_set, T, (h - h_val) / dh_dT)
    end
    (T, ΔT) = _select_solution(saturated, T, ΔT, T_0)

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
    T_init_min = TP.T_init_min(param_set)
    T_0 = max(T_init_min, T_unsat)
    saturated = _is_saturated_from_p(param_set, T_0, p, q_tot)

    T = T_0
    ΔT = zero(T)
    for _ in 1:maxiter
        (q_liq, q_ice) = _condensate_partition_from_p(param_set, T, p, q_tot, Val(false))
        θ_li_val = liquid_ice_pottemp_given_pressure(param_set, T, p, q_tot, q_liq, q_ice)
        dθ_li_dT = ∂θ_li_∂T_sat_p(param_set, T, p, q_tot, Val(false))
        (T, ΔT) = _newton_update(param_set, T, (θ_li - θ_li_val) / dθ_li_dT)
    end
    (T, ΔT) = _select_solution(saturated, T, ΔT, T_0)

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
    T_init_min = TP.T_init_min(param_set)
    T_0 = max(T_init_min, T_unsat)
    saturated = _is_saturated(param_set, T_0, ρ, q_tot)

    T = T_0
    ΔT = zero(T)
    for _ in 1:maxiter
        (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot, Val(false))
        θ_li_val = liquid_ice_pottemp(param_set, T, ρ, q_tot, q_liq, q_ice)
        dθ_li_dT = ∂θ_li_∂T_sat_ρ(param_set, T, ρ, q_tot, Val(false))
        (T, ΔT) = _newton_update(param_set, T, (θ_li - θ_li_val) / dθ_li_dT)
    end
    (T, ΔT) = _select_solution(saturated, T, ΔT, T_0)

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
    T_init_min = TP.T_init_min(param_set)
    T_0 = max(T_init_min, T_unsat)
    saturated = _is_saturated(param_set, T_0, ρ, q_tot)

    T = T_0
    ΔT = zero(T)
    for _ in 1:maxiter
        (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot, Val(false))
        p_val = air_pressure(param_set, T, ρ, q_tot, q_liq, q_ice)
        dp_dT = ∂p_∂T_sat_ρ(param_set, T, ρ, q_tot, Val(false))
        (T, ΔT) = _newton_update(param_set, T, (p - p_val) / dp_dT)
    end
    (T, ΔT) = _select_solution(saturated, T, ΔT, T_0)

    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)
    return (; T, q_liq, q_ice, converged = _fixed_iters_converged(T, ΔT))
end

# ---------------------------------------------
# Internal helpers
# ---------------------------------------------

"""
    bound_upper_temperature(param_set, T_lo, T_hi)

Internal function. Bounds the upper temperature guess `T_hi` for bracket methods.

Returns `T_hi` bounded by `T_max`, ensuring it is at least `T_lo + 0.1` for
valid numerical initialization.
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
  otherwise `max(T_init_min, T_unsat)`.
- For bracket methods (`SecantMethod`, `BrentsMethod`): Constructs bracket `[T_lo, T_hi]`
  where `T_hi` is bounded by `T_ice` and `T_max`.
"""
@inline function _make_sa_solver(
    ::Type{RS.NewtonsMethod},
    param_set::APS,
    T_unsat,
    T_ice,
    T_guess,
)
    T_init_min = TP.T_init_min(param_set)
    T_init = T_guess isa Nothing ? max(T_init_min, T_unsat) : T_guess
    return RS.NewtonsMethod(T_init)
end

@inline function _make_sa_solver(
    ::Type{RS.NewtonsMethodAD},
    param_set::APS,
    T_unsat,
    T_ice,
    T_guess,
)
    T_init_min = TP.T_init_min(param_set)
    T_init = T_guess isa Nothing ? max(T_init_min, T_unsat) : T_guess
    return RS.NewtonsMethodAD(T_init)
end

@inline function _make_sa_solver(
    ::Type{RS.SecantMethod},
    param_set::APS,
    T_unsat,
    T_ice,
    T_guess,
)
    T_init_min = TP.T_init_min(param_set)
    T_lo = T_guess isa Nothing ? max(T_init_min, T_unsat) : max(T_init_min, T_guess)
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
    T_init_min = TP.T_init_min(param_set)
    # BrentsMethod requires strict bracketing - ignore T_guess
    T_lo = max(T_init_min, T_unsat)
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
    return T -> begin
        (_q_liq, _q_ice) = condensate_partition(param_set, T, ρ, q_tot)
        liquid_ice_pottemp(param_set, T, ρ, q_tot, _q_liq, _q_ice) - θ_li
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
    return T -> begin
        (_q_liq, _q_ice) = condensate_partition(param_set, T, ρ, q_tot)
        air_pressure(param_set, T, ρ, q_tot, _q_liq, _q_ice) - p
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
    # p approaching p_v_sat cannot produce a division by zero inside a kernel.
    Δp = p - p_v_sat
    amplification = ifelse(Δp ≥ ϵ_numerics(FT), p / Δp, one(Δp))
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
function _find_zero_and_convergence(
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
        thermo_var,
        q_tot,
        maxiter,
        relative_temperature_tol,
        T_guess,
        temp_unsat_func,
        q_sat_unsat_func,
        temp_ice_func,
        roots_func,
    )

Generic kernel for saturation adjustment. Handles unsaturated check, solver initialization,
and root-finding for all `saturation_adjustment` variants.

# Arguments
- `Method`: Root-solving method type (e.g., `RS.NewtonsMethod`, `RS.SecantMethod`).
- `param_set`: Thermodynamics parameter set.
- `thermo_var`: Thermodynamic variable to match (e.g., `e_int`, `h`, `p`, `θ_li`).
- `q_tot`: Total specific humidity.
- `maxiter`: Maximum iterations for the solver.
- `relative_temperature_tol`: Relative tolerance for temperature solution.
- `T_guess`: Optional initial temperature guess.
- `temp_unsat_func`: Closure `(param_set, thermo_var, q_tot) -> T_unsat`.
- `q_sat_func`: Closure `(param_set, T, q_tot) -> q_vap_sat`.
- `temp_ice_func`: Closure `(param_set, thermo_var, q_tot) -> T_ice` for upper bound.
- `roots_func`: Closure `(T) -> residual` or `(T) -> (residual, derivative)`.

# Returns
- `(T, converged)`: Temperature and convergence flag.
"""
@inline function _saturation_adjustment_generic(
    ::Type{Method},  # RootSolvers MethodType
    param_set::APS,
    thermo_var, # This can be e_int or h
    q_tot,
    maxiter,
    relative_temperature_tol,
    T_guess,
    # The following are functions passed in to customize the behavior
    temp_unsat_func,     # (param_set, thermo_var, q_tot) -> T
    q_sat_func,          # (param_set, T, q_tot) -> q_sat (for unsaturated check)
    temp_ice_func,       # (param_set, thermo_var, q_tot) -> T_ice (for upper bound)
    roots_func,          # (T) -> val or (val, deriv)
) where {Method}
    _T_min = TP.T_min(param_set)
    tol =
        relative_temperature_tol isa Real ?
        RS.RelativeSolutionTolerance(relative_temperature_tol) :
        relative_temperature_tol

    # Encapsulated "Unsaturated Check" logic
    T_unsat = max(_T_min, temp_unsat_func(param_set, thermo_var, q_tot))

    q_v_sat = q_sat_func(param_set, T_unsat, q_tot)
    if q_tot <= q_v_sat
        return (T_unsat, true)
    end

    # Saturated case: solve for T
    T_ice = temp_ice_func(param_set, thermo_var, q_tot)

    # Initialize solver (logic merged from config_sa_method.jl)
    solver = _make_sa_solver(Method, param_set, T_unsat, T_ice, T_guess)

    (T, converged) = _find_zero_and_convergence(
        roots_func,
        solver,
        RS.CompactSolution(),
        tol,
        maxiter,
    )

    return (T, converged)
end

# -------------------------------
# Derivatives for Newton's method
# -------------------------------

"""
    _saturation_derivative_vars(param_set, T, ρ, q_tot, q_vap_sat, ::Val{:ρ})

Helper to compute common intermediate variables for saturation derivatives at fixed
density. The fixed-pressure counterpart is `_saturation_derivative_vars_p`, which takes
`p` rather than `ρ` so that the phase partition stays consistent with the pressure.

Returns a named tuple with phase partition and derivative information:
- `λ`, `q_c`, `q_liq`, `q_ice`: Phase partition variables
- `∂λ_∂T`, `∂qvs_∂T`: Temperature derivatives of liquid fraction and saturation humidity
- `∂q_liq_∂T`, `∂q_ice_∂T`: Phase partition derivatives
"""
@inline function _saturation_derivative_vars(
    param_set::APS,
    T,
    ρ,
    q_tot,
    q_vap_sat,
    ::Val{:ρ},
)
    λ = liquid_fraction_ramp(param_set, T)
    q_c = saturation_excess(param_set, T, ρ, q_tot)
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
    cvm_unsat = cv_m(param_set, q_tot, zero(q_tot), zero(q_tot))

    vars = _saturation_derivative_vars(param_set, T, ρ, q_tot, q_vap_sat, Val(:ρ))

    # Component internal energies
    e_vap = internal_energy_vapor(param_set, T)
    e_liq = internal_energy_liquid(param_set, T)
    e_ice = internal_energy_ice(param_set, T)

    # Full derivative: cv_m + Σ(e_i * ∂q_i/∂T)
    cvm_sat = cv_m(param_set, q_tot, vars.q_liq, vars.q_ice)
    de_dT_sat =
        cvm_sat + e_vap * vars.∂qvs_∂T + e_liq * vars.∂q_liq_∂T + e_ice * vars.∂q_ice_∂T

    return _select_sat_branch(unsat_branch, q_tot, q_vap_sat, cvm_unsat, de_dT_sat)
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
    cvm_unsat = cv_m(param_set, q_tot, zero(q_tot), zero(q_tot))

    vars = _saturation_derivative_vars_p(param_set, T, p, q_tot)

    # Component internal energies
    e_vap = internal_energy_vapor(param_set, T)
    e_liq = internal_energy_liquid(param_set, T)
    e_ice = internal_energy_ice(param_set, T)

    # Full derivative: cv_m + Σ(e_i * ∂q_i/∂T)
    cvm_sat = cv_m(param_set, q_tot, vars.q_liq, vars.q_ice)
    de_dT_sat =
        cvm_sat + e_vap * vars.∂qvs_∂T + e_liq * vars.∂q_liq_∂T + e_ice * vars.∂q_ice_∂T

    return _select_sat_branch(unsat_branch, q_tot, vars.q_vap_sat, cvm_unsat, de_dT_sat)
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
    cpm_unsat = cp_m(param_set, q_tot, zero(q_tot), zero(q_tot))

    vars = _saturation_derivative_vars_p(param_set, T, p, q_tot)

    # Component enthalpies
    h_vap = enthalpy_vapor(param_set, T)
    h_liq = enthalpy_liquid(param_set, T)
    h_ice = enthalpy_ice(param_set, T)

    # Full derivative: cp_m + Σ(h_i * ∂q_i/∂T)
    cpm_sat = cp_m(param_set, q_tot, vars.q_liq, vars.q_ice)
    dh_dT_sat =
        cpm_sat + h_vap * vars.∂qvs_∂T + h_liq * vars.∂q_liq_∂T + h_ice * vars.∂q_ice_∂T

    return _select_sat_branch(unsat_branch, q_tot, vars.q_vap_sat, cpm_unsat, dh_dT_sat)
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
    vars = _saturation_derivative_vars_p(param_set, T, p, q_tot)
    q_vap_sat = vars.q_vap_sat

    # Current thermodynamic state
    R_m = gas_constant_air(param_set, q_tot, vars.q_liq, vars.q_ice)
    _cp_m = cp_m(param_set, q_tot, vars.q_liq, vars.q_ice)
    p0 = TP.p_ref_theta(param_set)
    α = R_m / _cp_m
    Π = fast_power(p / p0, α)
    θ = T / Π
    L_c = humidity_weighted_latent_heat(param_set, vars.q_liq, vars.q_ice)
    F = 1 - L_c / (_cp_m * T)

    # Unsaturated: α is constant, L_c = 0, F = 1, so ∂θ_li/∂T = θ/T
    dθ_li_dT_unsat = θ / T

    # Parameters
    R_v = TP.R_v(param_set)
    cp_v = TP.cp_v(param_set)
    cp_l = TP.cp_l(param_set)
    cp_i = TP.cp_i(param_set)
    LH_v0 = TP.LH_v0(param_set)
    LH_s0 = TP.LH_s0(param_set)

    # ∂R_m/∂T and ∂cp_m/∂T
    ∂R_m_∂T = R_v * vars.∂qvs_∂T
    ∂cp_m_∂T = (cp_l - cp_v) * vars.∂q_liq_∂T + (cp_i - cp_v) * vars.∂q_ice_∂T

    # ∂α/∂T = (∂R_m * cp_m - R_m * ∂cp_m) / cp_m²
    ∂α_∂T = (∂R_m_∂T * _cp_m - R_m * ∂cp_m_∂T) / _cp_m^2

    # ∂θ/∂T = θ * (1/T + ln(p₀/p) * ∂α/∂T)
    ∂θ_∂T = θ * (1 / T + log(p0 / p) * ∂α_∂T)

    # ∂L_c/∂T
    ∂L_c_∂T = LH_v0 * vars.∂q_liq_∂T + LH_s0 * vars.∂q_ice_∂T

    # ∂F/∂T = -1/(cp_m*T) * (∂L_c/∂T - L_c * (1/T + ∂cp_m/∂T / cp_m))
    ∂F_∂T = -1 / (_cp_m * T) * (∂L_c_∂T - L_c * (1 / T + ∂cp_m_∂T / _cp_m))

    # Product rule: ∂(θ * F)/∂T
    dθ_li_dT_sat = ∂θ_∂T * F + θ * ∂F_∂T

    return _select_sat_branch(
        unsat_branch,
        q_tot,
        q_vap_sat,
        dθ_li_dT_unsat,
        dθ_li_dT_sat,
    )
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

    vars = _saturation_derivative_vars(param_set, T, ρ, q_tot, q_vap_sat, Val(:ρ))

    # Current thermodynamic state
    R_m = gas_constant_air(param_set, q_tot, vars.q_liq, vars.q_ice)
    _cp_m = cp_m(param_set, q_tot, vars.q_liq, vars.q_ice)
    p0 = TP.p_ref_theta(param_set)
    p = ρ * R_m * T  # ideal gas at fixed ρ
    α = R_m / _cp_m
    Π = fast_power(p / p0, α)
    θ = T / Π
    L_c = humidity_weighted_latent_heat(param_set, vars.q_liq, vars.q_ice)
    F = 1 - L_c / (_cp_m * T)

    # Unsaturated: q_liq = q_ice = 0, L_c = 0, F = 1, R_m constant
    # θ = (p0/(ρ R_m))^α · T^(1-α)  =>  dθ/dT = (1-α)·θ/T
    dθ_li_dT_unsat = θ * (1 - α) / T

    # Parameters
    R_v = TP.R_v(param_set)
    cp_v = TP.cp_v(param_set)
    cp_l = TP.cp_l(param_set)
    cp_i = TP.cp_i(param_set)
    LH_v0 = TP.LH_v0(param_set)
    LH_s0 = TP.LH_s0(param_set)

    # ∂R_m/∂T and ∂cp_m/∂T
    ∂R_m_∂T = R_v * vars.∂qvs_∂T
    ∂cp_m_∂T = (cp_l - cp_v) * vars.∂q_liq_∂T + (cp_i - cp_v) * vars.∂q_ice_∂T

    # ∂α/∂T = (∂R_m · cp_m - R_m · ∂cp_m) / cp_m²
    ∂α_∂T = (∂R_m_∂T * _cp_m - R_m * ∂cp_m_∂T) / _cp_m^2

    # ∂θ/∂T at fixed ρ: extra -α/T and -∂R_m/∂T/cp_m vs. fixed p
    ∂θ_∂T = θ * ((1 - α) / T + log(p0 / p) * ∂α_∂T - ∂R_m_∂T / _cp_m)

    # ∂L_c/∂T
    ∂L_c_∂T = LH_v0 * vars.∂q_liq_∂T + LH_s0 * vars.∂q_ice_∂T

    # ∂F/∂T = -1/(cp_m·T) · (∂L_c/∂T - L_c · (1/T + ∂cp_m/∂T / cp_m))
    ∂F_∂T = -1 / (_cp_m * T) * (∂L_c_∂T - L_c * (1 / T + ∂cp_m_∂T / _cp_m))

    # Product rule: ∂(θ · F)/∂T
    dθ_li_dT_sat = ∂θ_∂T * F + θ * ∂F_∂T

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

    vars = _saturation_derivative_vars(param_set, T, ρ, q_tot, q_vap_sat, Val(:ρ))
    R_m = gas_constant_air(param_set, q_tot, vars.q_liq, vars.q_ice)

    # ∂R_m/∂T
    ∂R_m_∂T = R_v * vars.∂qvs_∂T

    # p = ρ R_m T  =>  ∂p/∂T = ρ(R_m + T · ∂R_m/∂T)
    dp_dT_sat = ρ * (R_m + T * ∂R_m_∂T)

    return _select_sat_branch(unsat_branch, q_tot, q_vap_sat, dp_dT_unsat, dp_dT_sat)
end
