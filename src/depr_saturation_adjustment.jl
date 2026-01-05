"""
    saturation_adjustment_ρeq(args...)

Deprecated. Use [`saturation_adjustment`](@ref) with `ρeq()` dispatch instead.
"""
@inline function saturation_adjustment_ρeq(
    sat_adjust_method::Type,
    param_set::APS,
    ρ,
    e_int,
    q_tot,
    maxiter::Int,
    tol::Real,
    T_guess = nothing,
)
    return saturation_adjustment(
        sat_adjust_method,
        param_set,
        ρeq(),
        ρ,
        e_int,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
end

@inline function saturation_adjustment_ρeq(
    param_set::APS,
    ρ,
    e_int,
    q_tot,
    maxiter::Int,
    tol::Real,
    T_guess = nothing,
)
    return saturation_adjustment(
        SecantMethod,
        param_set,
        ρeq(),
        ρ,
        e_int,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
end

"""
    saturation_adjustment_phq(args...)

Deprecated. Use [`saturation_adjustment`](@ref) with `phq()` dispatch instead.
"""
@inline function saturation_adjustment_phq(
    sat_adjust_method::Type,
    param_set::APS,
    p,
    h,
    q_tot,
    maxiter::Int,
    tol::Real,
    T_guess = nothing,
)
    return saturation_adjustment(
        sat_adjust_method,
        param_set,
        phq(),
        p,
        h,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
end

"""
    saturation_adjustment(args...)

Legacy wrapper for `saturation_adjustment` with `phase_type` argument.
Forwards to `saturation_adjustment(..., ::ρeq, ...)` with arguments `(ρ, e_int)` swapped.
"""
@inline function saturation_adjustment(
    sat_adjust_method::Type,
    param_set::APS,
    e_int::Number,
    ρ::Number,
    q_tot::Number,
    phase_type::Type,
    maxiter::Int,
    tol::Real,
    T_guess = nothing,
)
    (T, _...) = saturation_adjustment(
        sat_adjust_method,
        param_set,
        ρeq(),
        ρ,
        e_int,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
    return T
end

"""
    saturation_adjustment_given_phq(args...)

Deprecated. Use [`saturation_adjustment`](@ref) with `phq()` dispatch instead.
"""
@inline function saturation_adjustment_given_phq(
    sat_adjust_method::Type,
    param_set::APS,
    p,
    h,
    q_tot,
    phase_type::Type,
    maxiter::Int,
    tol::Real,
    T_guess = nothing,
)
    (T, _...) = saturation_adjustment(
        sat_adjust_method,
        param_set,
        phq(),
        p,
        h,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
    return T
end

"""
    saturation_adjustment_given_peq(args...)

Deprecated. Use [`saturation_adjustment`](@ref) with `peq()` dispatch instead.
"""
@inline function saturation_adjustment_given_peq(
    sat_adjust_method::Type,
    param_set::APS,
    p,
    e_int,
    q_tot,
    phase_type::Type,
    maxiter::Int,
    tol,
    T_guess = nothing,
)
    (T, _...) = saturation_adjustment(
        sat_adjust_method,
        param_set,
        peq(),
        p,
        e_int,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
    return T
end

"""
    saturation_adjustment_given_pθq(args...)

Deprecated. Use [`saturation_adjustment`](@ref) with `pθ_liq_ice_q()` dispatch instead.
"""
@inline function saturation_adjustment_given_pθq(
    sat_adjust_method::Type,
    param_set::APS,
    p::Number,
    θ_liq_ice::Number,
    q_tot::Number,
    phase_type::Type,
    maxiter::Int,
    tol::Real,
    T_guess = nothing,
)
    (T, _...) = saturation_adjustment(
        sat_adjust_method,
        param_set,
        pθ_liq_ice_q(),
        p,
        θ_liq_ice,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
    return T
end

"""
    saturation_adjustment_given_ρθq(args...)

Deprecated. Use [`saturation_adjustment`](@ref) with `ρθ_liq_ice_q()` dispatch instead.
"""
@inline function saturation_adjustment_given_ρθq(
    param_set::APS,
    ρ::Number,
    θ_liq_ice::Number,
    q_tot::Number,
    phase_type::Type,
    maxiter::Int,
    tol,
    T_guess = nothing,
    sat_adjust_method::Type = RS.SecantMethod,
)
    (T, _...) = saturation_adjustment(
        sat_adjust_method,
        param_set,
        ρθ_liq_ice_q(),
        ρ,
        θ_liq_ice,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
    return T
end

"""
    saturation_adjustment_ρpq(args...)

Deprecated. Use [`saturation_adjustment`](@ref) with `pρq()` dispatch instead.
"""
@inline function saturation_adjustment_ρpq(
    sat_adjust_method::Type,
    param_set::APS,
    ρ::Number,
    p::Number,
    q_tot::Number,
    phase_type::Type,
    maxiter::Int,
    tol::Real,
    T_guess = nothing,
)
    (T, _...) = saturation_adjustment(
        sat_adjust_method,
        param_set,
        pρq(),
        p,
        ρ,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
    return T
end
