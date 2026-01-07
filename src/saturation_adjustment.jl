# Saturation adjustment functions for various combinations of input variables

import RootSolvers as RS

export saturation_adjustment

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
        [T_guess::Union{Nothing, Real}]
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
- `T_guess`: Optional initial guess for the temperature [K].

# Returns
- `NamedTuple` `(; T, q_liq, q_ice, converged)`:
    - `T`: Temperature [K]
    - `q_liq`: Liquid specific humidity [kg/kg]
    - `q_ice`: Ice specific humidity [kg/kg]
    - `converged`: Boolean flag indicating if the solver converged

# Notes
- This function solves for `T` such that `e_int = internal_energy_sat(param_set, T, ρ, q_tot)` using
  root-finding, then computes `(q_liq, q_ice)` from [`condensate_partition`](@ref).
- For `ρe` formulation, `NewtonsMethod` is recommended (fast + analytic derivative available).
- For other formulations, `SecantMethod` or `BrentsMethod` are recommended.
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
) where {M}
    T_init_min = TP.T_init_min(param_set)
    sol_tol = tol isa Real ? RS.RelativeSolutionTolerance(tol) : tol

    # Temperature for unsaturated case (always computed)
    T_unsat = max(T_init_min, air_temperature(param_set, e_int, q_tot))

    # Check if unsaturated
    if q_tot <= q_vap_saturation(param_set, T_unsat, ρ)
        return (;
            T = T_unsat,
            q_liq = zero(T_unsat),
            q_ice = zero(T_unsat),
            converged = true,
        )
    end

    # Root function: e_int - internal_energy_sat(T, ρ, q_tot) = 0
    roots = _make_roots_function(M, param_set, ρ, e_int, q_tot)

    # Construct numerical method based on type
    numerical_method = sa_numerical_method(
        M,
        param_set,
        ρe(),
        ρ,
        e_int,
        q_tot,
        T_guess, # Use T_guess here, it will handle initialization internally
    )

    # Solve for temperature and check convergence
    (T, converged) = _find_zero_and_convergence(
        roots,
        numerical_method,
        solution_type(),
        sol_tol,
        maxiter,
    )

    # Compute equilibrium phase partition
    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)
    return (; T, q_liq, q_ice, converged)
end

"""
    saturation_adjustment(
        ::Type{M},  # RS.RootSolvingMethod type
        param_set,
        ::pe,
        p,
        e_int,
        q_tot,
        maxiter::Int,
        tol,
        [T_guess::Union{Nothing, Real}]
    )

Compute the saturation equilibrium temperature `T` and phase partition `(q_liq, q_ice)`
given pressure `p`, specific internal energy `e_int`, and total specific humidity `q_tot`.

Returns a `NamedTuple` `(; T, q_liq, q_ice, converged)`.
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
) where {M}
    e_int_sat_given_p =
        T ->
            internal_energy_sat(param_set, T, air_density(param_set, T, p, q_tot), q_tot)

    q_sat_unsat_p =
        (param_set, T, q_tot) ->
            q_vap_saturation(param_set, T, air_density(param_set, T, p, q_tot))

    # Construct numerical method for the root solver
    numerical_method = sa_numerical_method(
        M,
        param_set,
        pe(),
        p,
        e_int,
        q_tot,
        T_guess,
    )

    temp_from_energy_unsat =
        (param_set, e_int, q_tot) -> air_temperature(param_set, e_int, q_tot)

    (T, converged) = _saturation_adjustment_generic(
        numerical_method,
        param_set,
        e_int,
        q_tot,
        maxiter,
        tol,
        temp_from_energy_unsat,
        q_sat_unsat_p,
        e_int_sat_given_p,
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
        [T_guess::Union{Nothing, Real}]
    )

Compute the saturation equilibrium temperature `T` and phase partition `(q_liq, q_ice)`
given pressure `p`, specific enthalpy `h`, and total specific humidity `q_tot`.

Returns a `NamedTuple` `(; T, q_liq, q_ice, converged)`.

# Arguments
- `M`: The numerical method type for root-finding. Supported types:
  `RS.SecantMethod`, `RS.BrentsMethod`, `RS.NewtonsMethod`, `RS.NewtonsMethodAD`.
- `param_set`: Thermodynamics parameter set.
- `p`: Pressure of moist air.
- `h`: Specific enthalpy.
- `q_tot`: Total specific humidity.
- `maxiter`: Maximum iterations for the solver.
- `tol`: Relative tolerance for the temperature solution (or a `RS.RelativeSolutionTolerance`).
- `T_guess`: Optional initial guess for the temperature.
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
) where {M}
    h_sat_given_p =
        T ->
            enthalpy_sat(param_set, T, air_density(param_set, T, p, q_tot), q_tot)

    q_sat_unsat_p =
        (param_set, T, q_tot) ->
            q_vap_saturation(param_set, T, air_density(param_set, T, p, q_tot))

    # Unsaturated temperature estimate from (h, q_tot) with zero condensate.
    temp_from_hq_unsat =
        (param_set, h, q_tot) ->
            air_temperature(param_set, ph(), h, q_tot, 0, 0)

    # Construct numerical method for the root solver
    numerical_method = sa_numerical_method(
        M,
        param_set,
        ph(),
        p,
        h,
        q_tot,
        T_guess,
    )

    (T, converged) = _saturation_adjustment_generic(
        numerical_method,
        param_set,
        h,
        q_tot,
        maxiter,
        tol,
        temp_from_hq_unsat,
        q_sat_unsat_p,
        h_sat_given_p,
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
        θ_liq_ice::Real,
        q_tot::Real,
        maxiter::Int,
        tol,
        [T_guess::Union{Nothing, Real}]
    )

Compute the saturation equilibrium temperature `T` and phase partition `(q_liq, q_ice)`
given pressure `p`, liquid-ice potential temperature `θ_liq_ice`, and total specific humidity `q_tot`.

Returns a `NamedTuple` `(; T, q_liq, q_ice, converged)`.
"""
function saturation_adjustment(
    ::Type{M},  # RS.AbstractMethod type
    param_set::APS,
    ::pθ_li,
    p,
    θ_liq_ice,
    q_tot,
    maxiter::Int,
    tol,
    T_guess = nothing,
) where {M}
    θ_liq_ice_sat_given_p =
        T -> begin
            _ρ = air_density(param_set, T, p, q_tot)
            (_q_liq, _q_ice) = condensate_partition(param_set, T, _ρ, q_tot)
            liquid_ice_pottemp_given_pressure(param_set, T, p, q_tot, _q_liq, _q_ice)
        end

    temp_from_θ_liq_ice_func =
        (param_set, θ_liq_ice, q_tot) ->
            air_temperature(param_set, pθ_li(), p, θ_liq_ice, q_tot)

    q_sat_unsat_p =
        (param_set, T, q_tot) ->
            q_vap_saturation(param_set, T, air_density(param_set, T, p, q_tot))

    # Construct numerical method for the root solver
    numerical_method = sa_numerical_method(
        M,
        param_set,
        pθ_li(),
        p,
        θ_liq_ice,
        q_tot,
        T_guess,
    )

    (T, converged) = _saturation_adjustment_generic(
        numerical_method,
        param_set,
        θ_liq_ice,
        q_tot,
        maxiter,
        tol,
        temp_from_θ_liq_ice_func,
        q_sat_unsat_p,
        θ_liq_ice_sat_given_p,
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
        θ_liq_ice::Real,
        q_tot::Real,
        maxiter::Int,
        tol,
        [T_guess::Union{Nothing, Real}]
    )

Compute the saturation equilibrium temperature `T` and phase partition `(q_liq, q_ice)`
given density `ρ`, liquid-ice potential temperature `θ_liq_ice`, and total specific humidity `q_tot`.

Returns a `NamedTuple` `(; T, q_liq, q_ice, converged)`.
"""
function saturation_adjustment(
    ::Type{M},  # RS.AbstractMethod type
    param_set::APS,
    ::ρθ_li,
    ρ,
    θ_liq_ice,
    q_tot,
    maxiter::Int,
    tol,
    T_guess = nothing,
) where {M}
    θ_liq_ice_sat_given_ρ =
        T -> begin
            (_q_liq, _q_ice) = condensate_partition(param_set, T, ρ, q_tot)
            liquid_ice_pottemp(param_set, T, ρ, q_tot, _q_liq, _q_ice)
        end

    temp_from_θ_liq_ice_func =
        (param_set, θ_liq_ice, q_tot) ->
            air_temperature(param_set, ρθ_li(), ρ, θ_liq_ice, q_tot)

    q_sat_unsat_ρ = (param_set, T, q_tot) ->
        q_vap_saturation(param_set, T, ρ)

    # Construct numerical method for the root solver
    numerical_method = sa_numerical_method(
        M,
        param_set,
        ρθ_li(),
        ρ,
        θ_liq_ice,
        q_tot,
        T_guess,
    )

    (T, converged) = _saturation_adjustment_generic(
        numerical_method,
        param_set,
        θ_liq_ice,
        q_tot,
        maxiter,
        tol,
        temp_from_θ_liq_ice_func,
        q_sat_unsat_ρ,
        θ_liq_ice_sat_given_ρ,
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
        [T_guess::Union{Nothing, Real}]
    )

Compute the saturation equilibrium temperature `T` and phase partition `(q_liq, q_ice)`
given pressure `p`, density `ρ`, and total specific humidity `q_tot`.

Returns a `NamedTuple` `(; T, q_liq, q_ice, converged)`.
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
) where {M}
    pressure_sat_given_ρ =
        T -> begin
            (_q_liq, _q_ice) = condensate_partition(param_set, T, ρ, q_tot)
            air_pressure(param_set, T, ρ, q_tot, _q_liq, _q_ice)
        end

    temp_from_pρq_func =
        (param_set, p, q_tot) ->
            air_temperature(param_set, pρ(), p, ρ, q_tot)

    q_sat_unsat_ρ = (param_set, T, q_tot) ->
        q_vap_saturation(param_set, T, ρ)

    # Construct numerical method for the root solver
    numerical_method = sa_numerical_method(
        M,
        param_set,
        pρ(),
        p,
        ρ,
        q_tot,
        T_guess,
    )

    # Using generic helper with p as target thermo_var
    (T, converged) = _saturation_adjustment_generic(
        numerical_method,
        param_set,
        p, # thermo_var (target p)
        q_tot,
        maxiter,
        tol,
        temp_from_pρq_func,
        q_sat_unsat_ρ,
        pressure_sat_given_ρ,
    )

    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)
    return (; T, q_liq, q_ice, converged)
end

# ---------------------------------------------
# Helper functions 
# ---------------------------------------------

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
temperature-dependent liquid fraction (see [`liquid_fraction`](@ref)) and saturation 
excess (see [`saturation_excess`](@ref)).
"""
@inline function internal_energy_sat(param_set::APS, T, ρ, q_tot)
    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)
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
temperature-dependent liquid fraction (see [`liquid_fraction`](@ref)) and saturation 
excess (see [`saturation_excess`](@ref)).
"""
@inline function enthalpy_sat(param_set::APS, T, ρ, q_tot)
    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)
    return enthalpy(param_set, T, q_tot, q_liq, q_ice)
end



"""
    _make_roots_function(::Type{M}, param_set, ρ, e_int, q_tot)

Helper function to create the root function for Newton's method (with derivative) or
other methods (without derivative), dispatching on the method type.
"""
@inline function _make_roots_function(
    ::Type{RS.NewtonsMethod},
    param_set::APS,
    ρ,
    e_int,
    q_tot,
)
    return _T -> begin
        T_val = ReLU(_T)
        f = e_int - internal_energy_sat(param_set, T_val, ρ, q_tot)
        (f, -∂e_int_∂T_sat(param_set, T_val, ρ, q_tot))
    end
end

@inline function _make_roots_function(
    ::Type{M},
    param_set::APS,
    ρ,
    e_int,
    q_tot,
) where {M}
    return _T -> begin
        T_val = ReLU(_T)
        e_int - internal_energy_sat(param_set, T_val, ρ, q_tot)
    end
end

"""
    _phase_partition_from_T_p(param_set, T, p, q_tot)

Helper to compute equilibrium phase partition given temperature, pressure, and total humidity.
Returns `(ρ, q_liq, q_ice)` tuple.
"""
@inline function _phase_partition_from_T_p(param_set::APS, T, p, q_tot)
    ρ = air_density(param_set, T, p, q_tot)
    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)
    return (ρ, q_liq, q_ice)
end

"""
    _find_zero_with_convergence_check(
        roots_func,
        numerical_method,
        solution_type,
        tol,
        maxiter,
        warning_func,
        warning_args...,
    )

Helper function to find the root of `roots_func` using `numerical_method`,
checking for convergence and issuing warnings or errors if necessary.

Used by `saturation_adjustment_*` functions to handle common solver logic.
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
        numerical_method,
        param_set,
        thermo_var,
        q_tot,
        maxiter,
        relative_temperature_tol,
        temp_from_var_unsat_func,
        q_sat_unsat_func,
        sat_val_func,
        warning_func,
        warning_args...,
    )

Generic kernel for saturation adjustment.

Encapsulates unsaturated check logic and root-finding for saturation adjustments.
Arguments `temp_from_var_unsat_func`, `q_sat_unsat_func`, and `sat_val_func`
are closures that capture any specific independent variables (like p or ρ).
"""
@inline function _saturation_adjustment_generic(
    numerical_method,  # Pre-constructed RootSolvers method instance
    param_set::APS,
    thermo_var, # This can be e_int or h
    q_tot,
    maxiter,
    relative_temperature_tol,
    # The following are functions passed in to customize the behavior
    temp_from_var_unsat_func,   # (param_set, thermo_var, q_tot) -> T
    q_sat_unsat_func,           # (param_set, T, q_tot) -> q_sat (for unsaturated check)
    sat_val_func,               # (T) -> val (to match thermo_var)
)
    _T_min = TP.T_min(param_set)
    tol =
        relative_temperature_tol isa Real ?
        RS.RelativeSolutionTolerance(relative_temperature_tol) :
        relative_temperature_tol

    # Encapsulated "Unsaturated Check" logic
    T_unsat = max(_T_min, temp_from_var_unsat_func(param_set, thermo_var, q_tot))

    q_v_sat = q_sat_unsat_func(param_set, T_unsat, q_tot)
    if q_tot <= q_v_sat
        return (T_unsat, true)
    end

    # Saturated case: find the root
    function roots(T)
        T_safe = max(T, _T_min)
        return sat_val_func(T_safe) - thermo_var
    end

    return _find_zero_and_convergence(
        roots,
        numerical_method,
        RS.CompactSolution(),
        tol,
        maxiter,
    )
end

# -------------------------------
# Derivatives for Newton's method
# -------------------------------

"""
    ∂e_int_∂T(param_set, T, e_int, ρ, q_tot, ...)

Derivative of internal energy with respect to temperature at saturation.

# Arguments
- `param_set`: Parameter set.
- `T`: Temperature.
- `e_int`: Internal energy (unused, for interface consistency).
- `ρ`: Density.
- `q_tot`: Total specific humidity.
"""
@inline function ∂e_int_∂T(
    param_set::APS,
    T,
    e_int,
    ρ,
    q_tot,
    λ = liquid_fraction(param_set, T),
    p_vap_sat = saturation_vapor_pressure(param_set, T),
)
    FT = eltype(param_set)
    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)
    q_cond = q_liq + q_ice
    q_vap_sat = q_vap_from_p_vap(param_set, T, ρ, p_vap_sat)
    L = latent_heat_mixed(param_set, T, λ)

    # ∂λ/∂T for λ = ((T - Tⁱ) / (Tᶠ - Tⁱ))^n
    Tᶠ = TP.T_freeze(param_set)
    Tⁱ = TP.T_icenuc(param_set)
    n = TP.pow_icenuc(param_set)
    ∂λ_∂T = ifelse(
        Tⁱ < T < Tᶠ,
        n / (Tᶠ - Tⁱ) * ((T - Tⁱ) / (Tᶠ - Tⁱ))^(n - 1),
        FT(0),
    )

    # ∂q_vap_sat/∂T from Clausius-Clapeyron
    _∂q_vap_sat_∂T = ∂q_vap_sat_∂T(param_set, λ, T, q_vap_sat, L)

    # Derivatives of phase fractions (when saturated, ∂q_c/∂T = -∂q_vap_sat/∂T)
    # q_c = q_tot - q_vap_sat
    ∂q_c_∂T = -_∂q_vap_sat_∂T
    ∂q_liq_∂T = ∂λ_∂T * q_cond + λ * ∂q_c_∂T
    ∂q_ice_∂T = -∂λ_∂T * q_cond + (1 - λ) * ∂q_c_∂T
    ∂q_vap_∂T = _∂q_vap_sat_∂T  # = -∂q_liq_∂T - ∂q_ice_∂T

    # Component internal energies
    e_vap = internal_energy_vapor(param_set, T)
    e_liq = internal_energy_liquid(param_set, T)
    e_ice = internal_energy_ice(param_set, T)

    # Full derivative: ∂e_int/∂T = cv_m + Σ(e_i * ∂q_i/∂T)
    cvm = cv_m(param_set, q_tot, q_liq, q_ice)
    return cvm + e_vap * ∂q_vap_∂T + e_liq * ∂q_liq_∂T + e_ice * ∂q_ice_∂T
end

"""
    ∂e_int_∂T_sat(param_set, T, ρ, q_tot)

Helper to compute `∂e_int/∂T` at saturation for Newton's method,
calculating necessary intermediate variables.
"""
@inline function ∂e_int_∂T_sat(param_set::APS, T, ρ, q_tot)
    p_vap_sat = saturation_vapor_pressure(param_set, T)
    λ = liquid_fraction(param_set, T)
    e_int_sat = internal_energy_sat(param_set, T, ρ, q_tot)
    return ∂e_int_∂T(param_set, T, e_int_sat, ρ, q_tot, λ, p_vap_sat)
end

"""
    ∂q_vap_sat_∂T(param_set, ts)
    ∂q_vap_sat_∂T(param_set, λ, T, q_vap_sat, [L])

Derivative of saturation vapor specific humidity with respect to temperature.

Computed via the Clausius-Clapeyron relation: `∂q_sat/∂T = q_sat * (L / (Rv * T^2) - 1 / T)`.
"""
@inline function ∂q_vap_sat_∂T(
    param_set::APS,
    λ,
    T,
    q_vap_sat,
    L = latent_heat_mixed(param_set, T, λ),
)
    R_v = TP.R_v(param_set)
    return q_vap_sat * (L / (R_v * T^2) - 1 / T)
end
