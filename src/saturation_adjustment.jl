# Saturation adjustment functions for various combinations of input variables

# TODO: Remove catches around freezing temperature (given we have continuous phase partition)

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

Compute the temperature that is consistent with

 - `sat_adjust_method` the numerical method to use.
    See the [`Thermodynamics`](@ref) for options.
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `e_int` internal energy
 - `ρ` (moist-)air density
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
 - `maxiter` maximum iterations for non-linear equation solve
 - `relative_temperature_tol` relative temperature tolerance
 - `T_guess` initial temperature guess

by finding the root of

`e_int - internal_energy_sat(param_set, T, ρ, q_tot, phase_type) = 0`

using the given numerical method `sat_adjust_method`.

See also [`saturation_adjustment`](@ref).
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

    T_1 = max(_T_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    if T_1 ≥ _T_min
        q_v_sat = q_vap_saturation(param_set, T_1, ρ, phase_type)
        unsaturated = q_tot <= q_v_sat
        if unsaturated
            return T_1
        end
    end
    _T_freeze = TP.T_freeze(param_set)
    @inline e_int_sat(T) =
        internal_energy_sat(param_set, ReLU(T), ρ, q_tot, phase_type)
    temperature_tol = _T_freeze * relative_temperature_tol
    e_int_upper = e_int_sat(_T_freeze + temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    if e_int < e_int_upper
        e_int_lower = e_int_sat(_T_freeze - temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
        if e_int_lower < e_int # < e_int_upper
            return _T_freeze
        end
    end
    @inline function roots(_T) # ff′
        T = ReLU(_T)
        return ifelse(
            sat_adjust_method <: RS.NewtonsMethod,
            begin
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
                f′ = ∂e_int_∂T(
                    param_set,
                    T,
                    e_int,
                    ρ,
                    q_tot,
                    phase_type,
                    λ,
                    p_vap_sat,
                    q,
                    cvm,
                )
                _e_int_sat =
                    internal_energy_sat(param_set, T, ρ, q_tot, phase_type, q)
                f = _e_int_sat - e_int
                (f, f′)
            end,
            e_int_sat(T) - e_int,
        )
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
        e_int,
        p,
        q_tot,
        phase_type,
        maxiter,
        relative_temperature_tol,
        T_guess,
    )

Compute the temperature that is consistent with

 - `sat_adjust_method` the numerical method to use.
    See the [`Thermodynamics`](@ref) for options.
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `e_int` internal energy
 - `p` air pressure
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
 - `maxiter` maximum iterations for non-linear equation solve
 - `relative_temperature_tol` relative temperature tolerance
 - `T_guess` initial temperature guess

by finding the root of

`e_int - internal_energy_sat(param_set, T, ρ(T), q_tot, phase_type) = 0`

where `ρ(T) = air_density(param_set, T, p, PhasePartition(q_tot))`

using the given numerical method `sat_adjust_method`.

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
    _T_min = TP.T_min(param_set)
    tol = RS.RelativeSolutionTolerance(relative_temperature_tol)

    T_1 = max(_T_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    @inline ρ_T(T) = air_density(param_set, T, p, PhasePartition(q_tot))
    ρ_1 = ρ_T(T_1)
    q_v_sat = q_vap_saturation(param_set, T_1, ρ_1, phase_type)
    unsaturated = q_tot <= q_v_sat
    if unsaturated && T_1 ≥ _T_min
        return T_1
    end
    _T_freeze = TP.T_freeze(param_set)
    @inline e_int_sat(T) =
        internal_energy_sat(param_set, ReLU(T), ρ_T(T), q_tot, phase_type)

    temperature_tol = _T_freeze * relative_temperature_tol
    e_int_upper = e_int_sat(_T_freeze + temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    e_int_lower = e_int_sat(_T_freeze - temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    if e_int_lower < e_int < e_int_upper
        return _T_freeze
    end
    @inline roots(T) = e_int_sat(T) - e_int

    numerical_method = sa_numerical_method_peq(
        sat_adjust_method,
        param_set,
        p,
        e_int,
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
        print_warning_peq,
        sat_adjust_method,
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
        phase_type,
        maxiter,
        relative_temperature_tol
        T_guess,
    )

Compute the temperature that is consistent with

 - `sat_adjust_method` the numerical method to use.
    See the [`Thermodynamics`](@ref) for options.
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `h` specific enthalpy
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
 - `maxiter` maximum iterations for non-linear equation solve
 - `relative_temperature_tol` relative temperature tolerance
 - `T_guess` initial temperature guess

by finding the root of

`h - specific_enthalpy_sat(param_set, T, ρ(T), q_tot, phase_type) = 0`

where `ρ(T) = air_density(param_set, T, p, PhasePartition(q_tot))`

using the given numerical method `sat_adjust_method`.

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
    _T_min = TP.T_min(param_set)
    tol = RS.RelativeSolutionTolerance(relative_temperature_tol)

    T_1 = max(
        _T_min,
        air_temperature_from_enthalpy(param_set, h, PhasePartition(q_tot)),
    ) # Assume all vapor
    @inline ρ_T(T) = air_density(param_set, T, p, PhasePartition(q_tot))
    ρ_1 = ρ_T(T_1)
    q_v_sat = q_vap_saturation(param_set, T_1, ρ_1, phase_type)
    unsaturated = q_tot <= q_v_sat
    if unsaturated && T_1 ≥ _T_min
        return T_1
    end
    _T_freeze = TP.T_freeze(param_set)
    @inline h_sat(T) =
        specific_enthalpy_sat(param_set, ReLU(T), ρ_T(T), q_tot, phase_type)

    temperature_tol = _T_freeze * relative_temperature_tol
    h_upper = h_sat(_T_freeze + temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    h_lower = h_sat(_T_freeze - temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    if h_lower < h < h_upper
        return _T_freeze
    end
    @inline roots(T) = h_sat(T) - h

    numerical_method = sa_numerical_method_phq(
        sat_adjust_method,
        param_set,
        p,
        h,
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
        print_warning_hpq,
        sat_adjust_method,
        h,
        p,
        q_tot,
        T_guess,
    )
end

"""
    saturation_adjustment_ρpq(
        sat_adjust_method,
        param_set,
        ρ,
        p,
        q_tot,
        phase_type,
        maxiter,
        relative_temperature_tol,
        T_guess,
    )
Compute the temperature that is consistent with
 - `sat_adjust_method` the numerical method to use.
    See the [`Thermodynamics`](@ref) for options.
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `p` pressure
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
 - `maxiter` maximum iterations for non-linear equation solve
 - `relative_temperature_tol` relative temperature tolerance
 - `T_guess` initial temperature guess
by finding the root of

```
T - air_temperature_given_ρp(
        param_set,
        p,
        ρ,
        PhasePartition_equil(param_set, T, ρ, q_tot, phase_type),
    )
```
using Newtons method using ForwardDiff.
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
    ΔT_min(::Type{FT})

Minimum interval for saturation adjustment using Secant method
"""
@inline ΔT_min(::Type{FT}) where {FT} = FT(3)

"""
    ΔT_max(::Type{FT})

Maximum interval for saturation adjustment using Secant method
"""
@inline ΔT_max(::Type{FT}) where {FT} = FT(10)

"""
    bound_upper_temperature(T_1, T_2)

Bounds the upper temperature, `T_2`, for
saturation adjustment using Secant method
"""
@inline function bound_upper_temperature(T_1::FT, T_2::FT) where {FT <: Real}
    T_2 = max(T_1 + ΔT_min(FT), T_2)
    return min(T_1 + ΔT_max(FT), T_2)
end

"""
    saturation_adjustment_given_ρθq(
        param_set,
        ρ,
        θ_liq_ice,
        q_tot,
        phase_type,
        maxiter,
        tol,
        T_guess,
    )

Compute the temperature `T` that is consistent with

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
 - `maxiter` maximum iterations for non-linear equation solve
 - `tol` absolute tolerance for saturation adjustment iterations. Can be one of:
    - `RelativeSolutionTolerance()` to stop when `|x_2 - x_1| < tol`
    - `ResidualTolerance()` to stop when `|f(x)| < tol`
    - `RelativeRelativeSolutionTolerance()` to stop when `|x_2 - x_1|/x_1 < tol`
 - `T_guess` initial temperature guess

by finding the root of

`θ_{liq_ice} - liquid_ice_pottemp_sat(param_set, T, ρ, phase_type, q_tot) = 0`

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
    _T_min = TP.T_min(param_set)
    T_init_min = TP.T_init_min(param_set)
    @inline air_temp(q) = air_temperature_given_ρθq(param_set, ρ, θ_liq_ice, q)
    T_1 = max(_T_min, air_temp(PhasePartition(q_tot))) # Assume all vapor
    q_v_sat = q_vap_saturation(param_set, T_1, ρ, phase_type)
    unsaturated = q_tot <= q_v_sat
    if unsaturated && T_1 ≥ _T_min
        return T_1
    end
    T_2 = air_temp(PhasePartition(q_tot, FT(0), q_tot)) # Assume all ice
    T_2 = bound_upper_temperature(T_1, T_2)
    @inline roots(T) =
        liquid_ice_pottemp_sat(param_set, ReLU(T), ρ, phase_type, q_tot) -
        θ_liq_ice

    numerical_method = RS.SecantMethod(T_init_min, T_2)

    return _find_zero_with_convergence_check(
        roots,
        numerical_method,
        RS.CompactSolution(),
        tol,
        maxiter,
        print_warning_ρθq,
        RS.SecantMethod,
        ρ,
        θ_liq_ice,
        q_tot,
    )
end

"""
    saturation_adjustment_given_pθq(
        sat_adjust_method,
        param_set,
        p,
        θ_liq_ice,
        q_tot,
        phase_type,
        maxiter,
        relative_temperature_tol,
        T_guess
    )

Compute the temperature `T` that is consistent with

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
 - `relative_temperature_tol` relative temperature tolerance
 - `maxiter` maximum iterations for non-linear equation solve
 - `sat_adjust_method` the numerical method to use.
 - `T_guess` initial temperature guess

by finding the root of

`θ_{liq_ice} - liquid_ice_pottemp_given_pressure(param_set, T, p, phase_type, q_tot) = 0`

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
    T_min = TP.T_min(param_set)
    T_freeze = TP.T_freeze(param_set)
    @inline air_temp(q) = air_temperature_given_pθq(param_set, p, θ_liq_ice, q)
    @inline function θ_liq_ice_closure(T)
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
        return liquid_ice_pottemp_given_pressure(
            param_set,
            T,
            oftype(T, p),
            q_pt,
        )
    end
    @inline q_vap_sat(T) =
        q_vap_saturation_from_pressure(param_set, q_tot, p, T, phase_type)
    T_1 = max(T_min, air_temp(PhasePartition(q_tot))) # Assume all vapor
    q_v_sat_1 = q_vap_sat(T_1)
    unsaturated = q_tot <= q_v_sat_1
    if unsaturated && T_1 ≥ T_min
        return T_1
    end
    temperature_tol = T_freeze * relative_temperature_tol
    θ_liq_ice_upper = θ_liq_ice_closure(T_freeze + temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    θ_liq_ice_lower = θ_liq_ice_closure(T_freeze - temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    if θ_liq_ice_lower < θ_liq_ice < θ_liq_ice_upper
        return T_freeze
    end
    roots(T) = oftype(T, θ_liq_ice) - θ_liq_ice_closure(T)

    numerical_method = sa_numerical_method_pθq(
        sat_adjust_method,
        param_set,
        p,
        θ_liq_ice,
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
        print_warning_pθq,
        sat_adjust_method,
        p,
        θ_liq_ice,
        q_tot,
    )
end

# Derivative of the internal energy with respect to temperature, needed 
# for Newton's method in saturation adjustment
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

    q_c = condensate_shum(q)
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

# Helper function for saturation adjustment when virtual temperature 
# and relative humidity are given
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
