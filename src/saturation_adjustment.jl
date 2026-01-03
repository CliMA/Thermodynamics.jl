# Saturation adjustment functions for various combinations of input variables

export saturation_adjustment_ρeq
export saturation_adjustment_phq

"""
    saturation_adjustment_ρeq(
        sat_adjust_method::Type,
        param_set::APS,
        ρ::Real,
        e_int::Real,
        q_tot::Real,
        maxiter::Int,
        tol::Real,
        [T_guess::Union{Nothing, Real}]
    )

Compute the saturation equilibrium temperature `T` and phase partition `(q_liq, q_ice)`
given density `ρ`, internal energy `e_int`, and total specific humidity `q_tot`.

Returns a tuple `(T, q_liq, q_ice)`.

# Arguments
- `sat_adjust_method`: The numerical method for root-finding. Supported:
  `SecantMethod`, `BrentsMethod`, `NewtonsMethod`, `NewtonsMethodAD`.
- `param_set`: An `AbstractParameterSet` containing thermodynamic parameters.
- `ρ`: Density of moist air.
- `e_int`: Specific internal energy.
- `q_tot`: Total specific humidity.
- `maxiter`: Maximum iterations for the solver.
- `tol`: Relative tolerance for the temperature solution.
- `T_guess`: Optional initial guess for the temperature.
"""
@inline function saturation_adjustment_ρeq(
    ::Type{sat_adjust_method},
    param_set::APS,
    ρ,
    e_int,
    q_tot,
    maxiter::Int,
    tol::Real,
    T_guess = nothing,
) where {sat_adjust_method}
    T_init_min = TP.T_init_min(param_set)
    sol_tol = RS.RelativeSolutionTolerance(tol)

    # Temperature for unsaturated case (always computed)
    T_unsat = max(T_init_min, air_temperature(param_set, e_int, q_tot))

    # Check if unsaturated
    if q_tot <= q_vap_saturation(param_set, T_unsat, ρ)
        return (T_unsat, zero(T_unsat), zero(T_unsat))
    end

    # Root function: e_int - internal_energy_sat(T, ρ, q_tot) = 0
    # For NewtonsMethod, also return the derivative
    @inline function roots(_T)
        T = ReLU(_T)
        f = e_int - internal_energy_sat(param_set, T, ρ, q_tot)
        if sat_adjust_method <: RS.NewtonsMethod
            return (f, -∂e_int_∂T_sat(param_set, T, ρ, q_tot))
        else
            return f
        end
    end

    # Construct numerical method based on type
    numerical_method = sa_numerical_method(
        sat_adjust_method,
        param_set,
        ρ,
        e_int,
        q_tot,
        T_guess, # Use T_guess here, it will handle initialization internally
    )

    # Solve for temperature and check convergence
    T = _find_zero_with_convergence_check(
        roots,
        numerical_method,
        solution_type(),
        sol_tol,
        maxiter,
        print_warning_ρeq,
        sat_adjust_method,
        ρ,
        e_int,
        q_tot,
    )

    # Compute equilibrium phase partition
    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)

    return (T, q_liq, q_ice)
end

# Default to SecantMethod when no root solver method is specified
@inline saturation_adjustment_ρeq(
    param_set::APS,
    ρ,
    e_int,
    q_tot,
    maxiter::Int,
    tol::Real,
    T_guess = nothing,
) = saturation_adjustment_ρeq(
    SecantMethod,
    param_set,
    ρ,
    e_int,
    q_tot,
    maxiter,
    tol,
    T_guess,
)

"""
    saturation_adjustment_phq(
        sat_adjust_method::Type,
        param_set::APS,
        p::Real,
        h::Real,
        q_tot::Real,
        maxiter::Int,
        tol::Real,
        [T_guess::Union{Nothing, Real}]
    )

Compute the saturation equilibrium temperature `T` and phase partition `(q_liq, q_ice)`
given pressure `p`, specific enthalpy `h`, and total specific humidity `q_tot`.

Returns a tuple `(T, q_liq, q_ice)`.

# Arguments
- `sat_adjust_method`: The numerical method for root-finding. Supported types:
  `SecantMethod`, `BrentsMethod`, `NewtonsMethod`, `NewtonsMethodAD`.
- `param_set`: An `AbstractParameterSet` containing thermodynamic parameters.
- `p`: Pressure of moist air.
- `h`: Specific enthalpy.
- `q_tot`: Total specific humidity.
- `maxiter`: Maximum iterations for the solver.
- `tol`: Relative tolerance for the temperature solution.
- `T_guess`: Optional initial guess for the temperature.
"""
@inline function saturation_adjustment_phq(
    ::Type{sat_adjust_method},
    param_set::APS,
    p,
    h,
    q_tot,
    maxiter::Int,
    tol::Real,
    T_guess = nothing,
) where {sat_adjust_method}
    
    @inline h_sat_given_p(T, param_set, p, q_tot) =
        enthalpy_sat(
            param_set,
            T,
            air_density(param_set, T, p, q_tot),
            q_tot,
        )

    T = _saturation_adjustment_p_thermo_q(
        sat_adjust_method,
        param_set,
        p,
        h,
        q_tot,
        maxiter,
        tol,
        T_guess,
        air_temperature_given_hq, 
        h_sat_given_p,
        sa_numerical_method_phq,
        print_warning_hpq,             
        sat_adjust_method,
        h,
        p,
        q_tot,
        T_guess, 
    )
    
    # Compute equilibrium phase partition
    ρ = air_density(param_set, T, p, q_tot)
    (q_liq, q_ice) = condensate_partition(param_set, T, ρ, q_tot)
    
    return (T, q_liq, q_ice)
end

# ---------------------------------------------
# Helper functions 
# ---------------------------------------------

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
function _find_zero_with_convergence_check(
    roots_func,
    numerical_method,
    solution_type,
    tol,
    maxiter,
    warning_func,
    warning_args...,
)
    sol =
        RS.find_zero(roots_func, numerical_method, solution_type, tol, maxiter)

    DataCollection.log_meta(sol)
    if !sol.converged
        if print_warning()
            # Pass `sol.root` and `tol.tol` to match original warning signatures
            warning_func(warning_args..., sol.root, maxiter, tol.tol)
        end
        if error_on_non_convergence()
            error("Failed to converge with printed set of inputs.")
        end
    end
    return sol.root
end

"""
    _saturation_adjustment_p_thermo_q(
        sat_method,
        param_set,
        p,
        thermo_var,
        q_tot,
        maxiter,
        relative_temperature_tol,
        T_guess,
        temp_from_var_func,
        sat_val_func,
        numerical_method_constructor,
        warning_func,
        warning_args...,
    )

Private kernel for saturation adjustment given pressure `p`, total humidity `q_tot`,
and a thermodynamic variable `thermo_var` (e.g., specific enthalpy or internal energy).

Encapsulates unsaturated check logic and root-finding for pressure-based adjustments.
"""
@inline function _saturation_adjustment_p_thermo_q(
    sat_method,
    param_set::APS,
    p,
    thermo_var, # This can be e_int or h
    q_tot,
    maxiter,
    relative_temperature_tol,
    T_guess,
    # The following are functions passed in to customize the behavior
    temp_from_var_func, # Function to compute T from the thermo_var
    sat_val_func,       # Function to compute the saturated value of the thermo_var
    numerical_method_constructor,
    warning_func,
    warning_args...,
)
    _T_min = TP.T_min(param_set)
    tol = RS.RelativeSolutionTolerance(relative_temperature_tol)

    # Encapsulated "Unsaturated Check" logic
    T_1 = max(
        _T_min,
        temp_from_var_func(param_set, thermo_var, q_tot),
    )
    @inline ρ_T(T) = air_density(param_set, T, p, q_tot)
    ρ_1 = ρ_T(T_1)
    q_v_sat = q_vap_saturation(param_set, T_1, ρ_1)
    if q_tot <= q_v_sat
        return T_1
    end

    # Saturated case: find the root
    @inline roots(T) =
        sat_val_func(T, param_set, p, q_tot) - thermo_var

    numerical_method = numerical_method_constructor(
        sat_method,
        param_set,
        p,
        thermo_var,
        q_tot,
        T_guess,
    )

    return _find_zero_with_convergence_check(
        roots,
        numerical_method,
        RS.CompactSolution(),
        tol,
        maxiter,
        warning_func,
        warning_args...,
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
    ∂λ_∂T = ifelse(Tⁱ < T < Tᶠ, n / (Tᶠ - Tⁱ) * ((T - Tⁱ) / (Tᶠ - Tⁱ))^(n - 1), FT(0))

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
@inline function ∂e_int_∂T_sat(
    param_set::APS,
    T,
    ρ,
    q_tot,
)
    p_vap_sat = saturation_vapor_pressure(param_set, T)
    λ = liquid_fraction(param_set, T)
    e_int_sat = internal_energy_sat(param_set, T, ρ, q_tot)
    return ∂e_int_∂T(
        param_set,
        T,
        e_int_sat,
        ρ,
        q_tot,
        λ,
        p_vap_sat,
    )
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

@inline function ∂q_vap_sat_∂T(param_set::APS, ts::ThermodynamicState)
    λ = liquid_fraction(param_set, ts)
    T = air_temperature(param_set, ts)
    q_vap_sat = vapor_specific_humidity(param_set, ts)
    return ∂q_vap_sat_∂T(param_set, λ, T, q_vap_sat)
end

