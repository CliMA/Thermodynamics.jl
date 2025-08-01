# Saturation adjustment functions for various combinations of input variables

"""
    saturation_adjustment(
        sat_adjust_method,
        param_set,
        e_int,
        ρ,
        q_tot,
        phase_type,
        maxiter,
        relative_temperature_tol,
        T_guess,
    )

Computes the saturation equilibrium temperature given internal energy `e_int`,
density `ρ`, and total specific humidity `q_tot`.

This function finds the temperature `T` that satisfies the root equation:
`e_int - internal_energy_sat(T, ρ, q_tot) = 0`.
It is the most common entry point for saturation adjustment.

# Arguments
- `sat_adjust_method`: The numerical method for root-finding (e.g., `NewtonsMethod`, `SecantMethod`).
- `param_set`: An `AbstractParameterSet` containing thermodynamic parameters.
- `e_int`: Specific internal energy.
- `ρ`: Density of moist air.
- `q_tot`: Total specific humidity (vapor + condensate).
- `phase_type`: A thermodynamic phase type (`PhaseEquil`, etc.).
- `maxiter`: Maximum iterations for the solver.
- `relative_temperature_tol`: Relative tolerance for the temperature solution.
- `T_guess`: An initial guess for the temperature.
"""
@inline function saturation_adjustment(
    ::Type{sat_adjust_method},
    param_set::APS,
    e_int::FT,
    ρ::FT,
    q_tot::FT,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::Real,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, sat_adjust_method, phase_type <: PhaseEquil}
    _T_min = TP.T_min(param_set)
    tol = RS.RelativeSolutionTolerance(relative_temperature_tol)

    # Unsaturated Check
    T_1 = max(_T_min, air_temperature(param_set, e_int, PhasePartition(q_tot)))
    if q_tot <= q_vap_saturation(param_set, T_1, ρ, phase_type)
        return T_1
    end

    # Saturated case (with special handling for Newton's method's derivative)
    @inline function roots(_T)
        T = ReLU(_T)
        e_int_sat_val = internal_energy_sat(param_set, T, ρ, q_tot, phase_type)
        f = e_int_sat_val - e_int

        # NewtonsMethod needs to return the derivative but other methods don't, 
        # which makes this function difficult to merge with others
        if sat_adjust_method <: RS.NewtonsMethod
            return (f, ∂e_int_∂T_sat(param_set, T, ρ, q_tot, phase_type))
        else
            return f
        end
    end

    numerical_method = sa_numerical_method(
        sat_adjust_method,
        param_set,
        ρ,
        e_int,
        q_tot,
        phase_type,
        T_guess,
    )

    return _find_zero_with_convergence_check(
        roots,
        numerical_method,
        solution_type(),
        tol,
        maxiter,
        print_warning_ρeq,
        sat_adjust_method,
        ρ,
        e_int,
        q_tot,
    )
end

"""
    saturation_adjustment_given_peq(
        sat_adjust_method,
        param_set,
        p,
        e_int,
        q_tot,
        ...
    )

Computes the saturation equilibrium temperature given pressure `p`, internal energy `e_int`,
and total specific humidity `q_tot`.

This function finds the temperature `T` that satisfies the root equation:
`e_int - internal_energy_sat(T, ρ(T, p), q_tot) = 0`.

See also [`saturation_adjustment`](@ref).
"""
@inline function saturation_adjustment_given_peq(
    ::Type{sat_adjust_method},
    param_set::APS,
    p::FT,
    e_int::FT,
    q_tot::FT,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::Real,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, sat_adjust_method, phase_type <: PhaseEquil}
    # Define how to compute saturated internal energy given p
    @inline e_int_sat_given_p(T, param_set, p, q_tot, phase_type) =
        internal_energy_sat(
            param_set,
            T,
            air_density(param_set, T, p, PhasePartition(q_tot)),
            q_tot,
            phase_type,
        )

    return _saturation_adjustment_p_thermo_q(
        sat_adjust_method,
        param_set,
        p,
        e_int,
        q_tot,
        phase_type,
        maxiter,
        relative_temperature_tol,
        T_guess,
        air_temperature,             # temp_from_var_func for e_int
        e_int_sat_given_p,           # sat_val_func
        sa_numerical_method_peq,     # constructor for numerical method
        print_warning_peq,           # warning function
        sat_adjust_method,           # arguments for the warning
        p,
        e_int,
        q_tot,
    )
end

"""
    saturation_adjustment_given_phq(
        sat_adjust_method,
        param_set,
        p,
        h,
        q_tot,
        ...
    )

Computes the saturation equilibrium temperature given pressure `p`, specific enthalpy `h`,
and total specific humidity `q_tot`.

This function finds the temperature `T` that satisfies the root equation:
`h - specific_enthalpy_sat(T, ρ(T, p), q_tot) = 0`.

See also [`saturation_adjustment`](@ref).
"""
@inline function saturation_adjustment_given_phq(
    ::Type{sat_adjust_method},
    param_set::APS,
    p::FT,
    h::FT,
    q_tot::FT,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::Real,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, sat_adjust_method, phase_type <: PhaseEquil}
    # Define how to compute saturated enthalpy given p
    @inline h_sat_given_p(T, param_set, p, q_tot, phase_type) =
        specific_enthalpy_sat(
            param_set,
            T,
            air_density(param_set, T, p, PhasePartition(q_tot)),
            q_tot,
            phase_type,
        )

    return _saturation_adjustment_p_thermo_q(
        sat_adjust_method,
        param_set,
        p,
        h,
        q_tot,
        phase_type,
        maxiter,
        relative_temperature_tol,
        T_guess,
        air_temperature_from_enthalpy, # temp_from_var_func for h
        h_sat_given_p,                 # sat_val_func
        sa_numerical_method_phq,       # constructor for numerical method
        print_warning_hpq,             # warning function
        sat_adjust_method,
        h,
        p,
        q_tot,
        T_guess, # arguments for the warning
    )
end

"""
    saturation_adjustment_ρpq(
        sat_adjust_method,
        param_set,
        ρ,
        p,
        q_tot,
        ...
    )

Computes the saturation equilibrium temperature given density `ρ`, pressure `p`,
and total specific humidity `q_tot`.

This function finds the temperature `T` that satisfies the consistency equation
`T - air_temperature_given_ρp(p, ρ, q_sat(T)) = 0`.

See also [`saturation_adjustment`](@ref).
"""
@inline function saturation_adjustment_ρpq(
    ::Type{sat_adjust_method},
    param_set::APS,
    ρ::FT,
    p::FT,
    q_tot::FT,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::Real,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, sat_adjust_method, phase_type <: PhaseEquil}
    tol = RS.RelativeSolutionTolerance(relative_temperature_tol)
    # Use `oftype` to preserve diagonalized type signatures:
    @inline roots(T) =
        T - air_temperature_given_ρp(
            param_set,
            oftype(T, p),
            oftype(T, ρ),
            PhasePartition_equil(
                param_set,
                T,
                oftype(T, ρ),
                oftype(T, q_tot),
                phase_type,
            ),
        )

    numerical_method = sa_numerical_method_ρpq(
        sat_adjust_method,
        param_set,
        ρ,
        p,
        q_tot,
        phase_type,
        T_guess,
    )

    return _find_zero_with_convergence_check(
        roots,
        numerical_method,
        RS.CompactSolution(),
        tol,
        maxiter,
        print_warning_ρpq,
        sat_adjust_method,
        ρ,
        p,
        q_tot,
    )
end

"""
    saturation_adjustment_given_ρθq(
        param_set,
        ρ,
        θ_liq_ice,
        q_tot,
        ...
    )

Computes the saturation equilibrium temperature given density `ρ`, liquid-ice potential
temperature `θ_liq_ice`, and total specific humidity `q_tot`.

This function finds the temperature `T` that satisfies the root equation:
`θ_liq_ice - liquid_ice_pottemp_sat(T, ρ, q_tot) = 0`. 
It uses the `SecantMethod`.

See also [`saturation_adjustment`](@ref).
"""
@inline function saturation_adjustment_given_ρθq(
    param_set::APS,
    ρ::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    ::Type{phase_type},
    maxiter::Int,
    tol::RS.AbstractTolerance,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, phase_type <: PhaseEquil}
    # Define the specific functions for this saturation adjustment method
    @inline T_1_func(q_pt) =
        air_temperature_given_ρθq(param_set, ρ, θ_liq_ice, q_pt)

    @inline q_vap_sat_func(T) = q_vap_saturation(param_set, T, ρ, phase_type)

    @inline residual_func(T) =
        liquid_ice_pottemp_sat(param_set, ReLU(T), ρ, phase_type, q_tot) -
        θ_liq_ice

    @inline numerical_method_constructor(T_1) = begin
        T_2 = air_temperature_given_ρθq(
            param_set,
            ρ,
            θ_liq_ice,
            PhasePartition(q_tot, FT(0), q_tot),
        ) # Assume all ice
        T_2 = bound_upper_temperature(T_1, T_2)
        return RS.SecantMethod(TP.T_init_min(param_set), T_2)
    end

    return _saturation_adjustment_θ_q_kernel(
        param_set,
        θ_liq_ice,
        q_tot,
        phase_type,
        maxiter,
        tol,
        T_1_func,           # T_1_func
        q_vap_sat_func,     # q_vap_sat_func
        residual_func,       # residual_func
        numerical_method_constructor, # numerical_method_constructor
        print_warning_ρθq,  # warning function
        RS.SecantMethod,
        ρ,
        θ_liq_ice,
        q_tot, # arguments for the warning
    )
end

"""
    saturation_adjustment_given_pθq(
        sat_adjust_method,
        param_set,
        p,
        θ_liq_ice,
        q_tot,
        ...
    )

Computes the saturation equilibrium temperature given pressure `p`, liquid-ice potential
temperature `θ_liq_ice`, and total specific humidity `q_tot`.

This function finds the temperature `T` that satisfies the root equation:
`θ_liq_ice - liquid_ice_pottemp_given_pressure(T, p, q_sat(T)) = 0`.

See also [`saturation_adjustment`](@ref).
"""
@inline function saturation_adjustment_given_pθq(
    ::Type{sat_adjust_method},
    param_set::APS,
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::Real,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, sat_adjust_method, phase_type <: PhaseEquil}
    tol = RS.RelativeSolutionTolerance(relative_temperature_tol)

    # Define the specific functions for this saturation adjustment method
    @inline T_1_func(q_pt) =
        air_temperature_given_pθq(param_set, p, θ_liq_ice, q_pt)

    @inline q_vap_sat_func(T) =
        q_vap_saturation_from_pressure(param_set, q_tot, p, T, phase_type)

    @inline residual_func(T) = begin
        q = PhasePartition(oftype(T, 0))
        λ = liquid_fraction(param_set, T, phase_type, q)
        q_pt = PhasePartition_equil_given_p(
            param_set,
            T,
            oftype(T, p),
            oftype(T, q_tot),
            phase_type,
            λ,
        )
        return oftype(T, θ_liq_ice) - liquid_ice_pottemp_given_pressure(
            param_set,
            T,
            oftype(T, p),
            q_pt,
        )
    end

    @inline numerical_method_constructor(T_1) = sa_numerical_method_pθq(
        sat_adjust_method,
        param_set,
        p,
        θ_liq_ice,
        q_tot,
        phase_type,
        T_guess,
    )

    return _saturation_adjustment_θ_q_kernel(
        param_set,
        θ_liq_ice,
        q_tot,
        phase_type,
        maxiter,
        tol,
        T_1_func,           # T_1_func
        q_vap_sat_func,     # q_vap_sat_func
        residual_func,       # residual_func
        numerical_method_constructor, # numerical_method_constructor
        print_warning_pθq,  # warning function
        sat_adjust_method,
        p,
        θ_liq_ice,
        q_tot, # arguments for the warning
    )
end

### Helper functions ###

"""
    ΔT_min(::Type{FT})

Minimum temperature interval for the Secant method bracketing.
"""
@inline ΔT_min(::Type{FT}) where {FT} = FT(3)

"""
    ΔT_max(::Type{FT})

Maximum temperature interval for the Secant method bracketing.
"""
@inline ΔT_max(::Type{FT}) where {FT} = FT(10)

"""
    bound_upper_temperature(T_1, T_2)

Bounds the upper temperature guess `T_2` for the Secant method
to prevent divergence.
"""
@inline function bound_upper_temperature(T_1::FT, T_2::FT) where {FT <: Real}
    T_2 = max(T_1 + ΔT_min(FT), T_2)
    return min(T_1 + ΔT_max(FT), T_2)
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

A generic helper function to find the root of `roots_func` using `numerical_method`,
and handle convergence checks and warnings.

This function abstracts the common logic of calling the root solver, logging metadata,
and issuing warnings or errors if the solver fails to converge. It is used by all
`saturation_adjustment_*` functions to reduce code duplication.
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
        sat_adjust_method,
        param_set,
        p,
        thermo_var,
        q_tot,
        phase_type,
        maxiter,
        relative_temperature_tol,
        T_guess,
        temp_from_var_func,
        sat_val_func,
        numerical_method_constructor,
        warning_func,
        warning_args...,
    )

A private kernel function that consolidates the logic for saturation adjustment when
pressure, total humidity, and a thermodynamic variable (like internal energy or enthalpy) are given.

This function encapsulates the common unsaturated check logic and root finding for
pressure-based saturation adjustment methods.
"""
@inline function _saturation_adjustment_p_thermo_q(
    sat_method,
    param_set::APS,
    p::FT,
    thermo_var::FT, # This can be e_int or h
    q_tot::FT,
    ::Type{phase_type},
    maxiter,
    relative_temperature_tol,
    T_guess,
    # The following are functions passed in to customize the behavior
    temp_from_var_func, # Function to compute T from the thermo_var
    sat_val_func,       # Function to compute the saturated value of the thermo_var
    numerical_method_constructor,
    warning_func,
    warning_args...,
) where {FT <: Real, phase_type <: PhaseEquil}
    _T_min = TP.T_min(param_set)
    tol = RS.RelativeSolutionTolerance(relative_temperature_tol)

    # Encapsulated "Unsaturated Check" logic
    T_1 = max(
        _T_min,
        temp_from_var_func(param_set, thermo_var, PhasePartition(q_tot)),
    )
    @inline ρ_T(T) = air_density(param_set, T, p, PhasePartition(q_tot))
    ρ_1 = ρ_T(T_1)
    q_v_sat = q_vap_saturation(param_set, T_1, ρ_1, phase_type)
    if q_tot <= q_v_sat
        return T_1
    end

    # Saturated case: find the root
    @inline roots(T) =
        sat_val_func(T, param_set, p, q_tot, phase_type) - thermo_var

    numerical_method = numerical_method_constructor(
        sat_method,
        param_set,
        p,
        thermo_var,
        q_tot,
        phase_type,
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

"""
    _saturation_adjustment_θ_q_kernel(
        param_set,
        θ_liq_ice,
        q_tot,
        phase_type,
        maxiter,
        tol,
        # Functional arguments to customize behavior
        T_1_func,
        q_vap_sat_func,
        residual_func,
        numerical_method_constructor,
        warning_func,
        warning_args...,
    )

A generic private kernel for saturation adjustment given liquid-ice potential temperature `θ_liq_ice`.
It encapsulates the common workflow for both pressure-based (`pθq`) and density-based (`ρθq`) methods.

# Customization via Functional Arguments
- `T_1_func`: A function `(q_pt) -> T` that computes the initial "all-vapor" temperature.
- `q_vap_sat_func`: A function `(T) -> q_v_sat` that computes the saturation specific humidity.
- `residual_func`: A function `(T) -> F` that computes the residual for the root-finder.
- `numerical_method_constructor`: A function `(T_1) -> method` that returns a configured numerical
  method instance (e.g., `SecantMethod` or a generic method from a helper).
- `warning_func`, `warning_args...`: The warning function and its arguments for convergence failure.
"""
@inline function _saturation_adjustment_θ_q_kernel(
    param_set::APS,
    θ_liq_ice::FT,
    q_tot::FT,
    ::Type{phase_type},
    maxiter::Int,
    tol::RS.AbstractTolerance,
    # --- Functional arguments ---
    T_1_func,
    q_vap_sat_func,
    residual_func,
    numerical_method_constructor,
    warning_func,
    warning_args...,
) where {FT <: Real, phase_type <: PhaseEquil}
    _T_min = TP.T_min(param_set)

    # 1. Unsaturated Check (logic is now generic)
    T_1 = max(_T_min, T_1_func(PhasePartition(q_tot))) # Assume all vapor
    if q_tot <= q_vap_sat_func(T_1)
        return T_1
    end

    # 2. Saturated case: setup and solve
    @inline roots(T) = residual_func(T)
    numerical_method = numerical_method_constructor(T_1)

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

"""
    ∂e_int_∂T(param_set, T, e_int, ρ, q_tot, phase_type, ...)

The derivative of internal energy with respect to temperature at saturation.
Note that the `e_int` argument is a placeholder for interface consistency 
and is not used.
"""
@inline function ∂e_int_∂T(
    param_set::APS,
    T::FT,
    e_int::FT,
    ρ::FT,
    q_tot::FT,
    ::Type{phase_type},
    λ = liquid_fraction(param_set, T, phase_type),
    p_vap_sat = saturation_vapor_pressure(
        param_set,
        phase_type,
        T,
        PhasePartition(FT(0)),
        λ,
    ),
    q = PhasePartition_equil(param_set, T, ρ, q_tot, p_vap_sat, λ),
    cvm = cv_m(param_set, q),
) where {FT <: Real, phase_type <: PhaseEquil}
    T_0 = TP.T_0(param_set)
    cv_v = TP.cv_v(param_set)
    cv_l = TP.cv_l(param_set)
    cv_i = TP.cv_i(param_set)
    e_int_v0 = TP.e_int_v0(param_set)
    e_int_i0 = TP.e_int_i0(param_set)

    q_c = condensate_specific_humidity(q)
    q_vap_sat = q_vap_from_p_vap(param_set, T, ρ, p_vap_sat)
    L = latent_heat_mixed(param_set, T, λ)

    Tᶠ = TP.T_freeze(param_set)
    Tⁱ = TP.T_icenuc(param_set)
    n = TP.pow_icenuc(param_set)
    ∂λ_∂T = ifelse(Tⁱ < T < Tᶠ, n * (1 / (Tᶠ - Tⁱ))^n * T^(n - 1), FT(0))

    _∂q_vap_sat_∂T = ∂q_vap_sat_∂T(param_set, λ, T, q_vap_sat, L)
    dcvm_dq_vap = cv_v - λ * cv_l - (1 - λ) * cv_i
    return cvm +
           (e_int_v0 + (1 - λ) * e_int_i0 + (T - T_0) * dcvm_dq_vap) *
           _∂q_vap_sat_∂T +
           q_c * e_int_i0 * ∂λ_∂T
end

"""
    ∂e_int_∂T_sat(param_set, T, ρ, q_tot, phase_type)

Helper function to compute the derivative of internal energy with respect to
temperature at saturation, for use in Newton's method. It computes all necessary
intermediate variables.
"""
@inline function ∂e_int_∂T_sat(
    param_set::APS,
    T::FT,
    ρ::FT,
    q_tot::FT,
    ::Type{phase_type},
) where {FT <: Real, phase_type <: PhaseEquil}
    λ = liquid_fraction(param_set, T, phase_type)
    p_vap_sat = saturation_vapor_pressure(
        param_set,
        phase_type,
        T,
        PhasePartition(FT(0)),
        λ,
    )
    q = PhasePartition_equil(param_set, T, ρ, q_tot, p_vap_sat, λ)
    cvm = cv_m(param_set, q)
    e_int_sat = internal_energy_sat(param_set, T, ρ, q_tot, phase_type)
    return ∂e_int_∂T(
        param_set,
        T,
        e_int_sat,
        ρ,
        q_tot,
        phase_type,
        λ,
        p_vap_sat,
        q,
        cvm,
    )
end

"""
    ∂q_vap_sat_∂T(param_set, ts)
    ∂q_vap_sat_∂T(param_set, λ, T, q_vap_sat)

The derivative of the saturation vapor specific humidity with respect to 
temperature, according to the Clausius Clapeyron relation and the definition 
of specific humidity `q = p_v/(ρ R_v T)` in terms of vapor pressure `p_v`.

It is computed either given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ts` ThermodynamicState

 or given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `λ` liquid fraction
 - `T` temperature
 - `q_vap_sat` saturation vapor specific humidity
"""
@inline function ∂q_vap_sat_∂T(
    param_set::APS,
    λ::FT,
    T::FT,
    q_vap_sat::FT,
    L = latent_heat_mixed(param_set, T, λ),
) where {FT <: Real}
    R_v = TP.R_v(param_set)
    return q_vap_sat * (L / (R_v * T^2) - 1 / T)
end

@inline function ∂q_vap_sat_∂T(param_set::APS, ts::ThermodynamicState)
    λ = liquid_fraction(param_set, ts)
    T = air_temperature(param_set, ts)
    q_vap_sat = vapor_specific_humidity(param_set, ts)
    return ∂q_vap_sat_∂T(param_set, λ, T, q_vap_sat)
end

"""
    virt_temp_from_RH(param_set, T, ρ, RH, phase_type)

Computes the virtual temperature from temperature `T`, density `ρ`, and
relative humidity `RH`.
"""
@inline function virt_temp_from_RH(
    param_set::APS,
    T::FT,
    ρ::FT,
    RH::FT,
    ::Type{phase_type},
) where {FT <: Real, phase_type <: ThermodynamicState}
    q_tot = RH * q_vap_saturation(param_set, T, ρ, phase_type)
    q_pt = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    return virtual_temperature(param_set, T, q_pt)
end
