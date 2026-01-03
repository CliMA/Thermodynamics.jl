#-------------------------------
# Older methods
#-------------------------------

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
    e_int,
    ρ,
    q_tot,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::Real,
    T_guess = nothing,
) where {sat_adjust_method, phase_type <: PhaseEquil}
    (T, _, _) = saturation_adjustment_ρeq(
        sat_adjust_method,
        param_set,
        ρ,
        e_int,
        q_tot,
        maxiter,
        relative_temperature_tol,
        T_guess,
    )
    return T
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
    p,
    e_int,
    q_tot,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::Real,
    T_guess = nothing,
) where {sat_adjust_method, phase_type <: PhaseEquil}
    # Define how to compute saturated internal energy given p
    @inline e_int_sat_given_p(T, param_set, p, q_tot) =
        internal_energy_sat(
            param_set,
            T,
            air_density(param_set, T, p, q_tot),
            q_tot,
        )

    return _saturation_adjustment_p_thermo_q(
        sat_adjust_method,
        param_set,
        p,
        e_int,
        q_tot,
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
`h - enthalpy_sat(T, ρ(T, p), q_tot) = 0`.

See also [`saturation_adjustment`](@ref).
"""
@inline function saturation_adjustment_given_phq(
    ::Type{sat_adjust_method},
    param_set::APS,
    p,
    h,
    q_tot,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::Real,
    T_guess = nothing,
) where {sat_adjust_method, phase_type <: PhaseEquil}
    # Define how to compute saturated enthalpy given p
    @inline h_sat_given_p(T, param_set, p, q_tot) =
        enthalpy_sat(
            param_set,
            T,
            air_density(param_set, T, p, q_tot),
            q_tot,
        )

    return _saturation_adjustment_p_thermo_q(
        sat_adjust_method,
        param_set,
        p,
        h,
        q_tot,
        maxiter,
        relative_temperature_tol,
        T_guess,
        air_temperature_given_hq, # temp_from_var_func for h
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
`T - air_temperature_given_pρq(p, ρ, q_sat(T)) = 0`.

See also [`saturation_adjustment`](@ref).
"""
@inline function saturation_adjustment_ρpq(
    ::Type{sat_adjust_method},
    param_set::APS,
    ρ,
    p,
    q_tot,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::Real,
    T_guess = nothing,
) where {sat_adjust_method, phase_type <: PhaseEquil}
    tol = RS.RelativeSolutionTolerance(relative_temperature_tol)
    # Use `oftype` to preserve diagonalized type signatures:
    @inline roots(T) =
        T - air_temperature_given_pρq(
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
    ρ,
    θ_liq_ice,
    q_tot,
    ::Type{phase_type},
    maxiter::Int,
    tol::RS.AbstractTolerance,
    T_guess = nothing,
) where {phase_type <: PhaseEquil}
    FT = eltype(param_set)
    # Define the specific functions for this saturation adjustment method
    @inline T_1_func(args...) =
        air_temperature_given_ρθq(param_set, ρ, θ_liq_ice, args...)

    @inline q_vap_sat_func(T) = q_vap_saturation(param_set, T, ρ)

    @inline residual_func(T) =
        liquid_ice_pottemp_sat(param_set, ReLU(T), ρ, PhaseEquil, q_tot) -
        θ_liq_ice

    @inline numerical_method_constructor(T_1) = begin
        T_2 = air_temperature_given_ρθq(
            param_set,
            ρ,
            θ_liq_ice,
            q_tot,
            FT(0),
            q_tot,
        ) # Assume all ice
        T_2 = bound_upper_temperature(param_set, T_1, T_2)
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
    p,
    θ_liq_ice,
    q_tot,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::Real,
    T_guess = nothing,
) where {sat_adjust_method, phase_type <: PhaseEquil}
    tol = RS.RelativeSolutionTolerance(relative_temperature_tol)

    # Define the specific functions for this saturation adjustment method
    @inline T_1_func(args...) =
        air_temperature_given_pθq(param_set, p, θ_liq_ice, args...)

    @inline q_vap_sat_func(T) =
        q_vap_saturation_from_pressure(param_set, q_tot, p, T)

    @inline residual_func(T) = begin
        λ = liquid_fraction(param_set, T)
        q_pt = PhasePartition_equil_given_p(
            param_set,
            T,
            oftype(T, p),
            oftype(T, q_tot),
            PhaseEquil,
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
    θ_liq_ice,
    q_tot,
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
) where {phase_type <: PhaseEquil}
    _T_min = TP.T_min(param_set)

    # 1. Unsaturated Check (logic is now generic)
    T_1 = max(_T_min, T_1_func(q_tot)) # Assume all vapor
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
    virt_temp_from_RH(param_set, T, ρ, RH, phase_type)

Computes the virtual temperature from temperature `T`, density `ρ`, and
relative humidity `RH`.
"""
@inline function virt_temp_from_RH(
    param_set::APS,
    T,
    ρ,
    RH,
    ::Type{phase_type},
) where {phase_type <: ThermodynamicState}
    q_tot = RH * q_vap_saturation(param_set, T, ρ, phase_type)
    q_pt = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    return virtual_temperature(param_set, T, q_pt)
end
