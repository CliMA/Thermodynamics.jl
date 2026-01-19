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
    res = saturation_adjustment(
        sat_adjust_method,
        param_set,
        ρe(),
        ρ,
        e_int,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
    return (res.T, res.q_liq, res.q_ice)
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
    res = saturation_adjustment(
        SecantMethod,
        param_set,
        ρe(),
        ρ,
        e_int,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
    return (res.T, res.q_liq, res.q_ice)
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
    res = saturation_adjustment(
        sat_adjust_method,
        param_set,
        ph(),
        p,
        h,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
    return (res.T, res.q_liq, res.q_ice)
end

"""
    saturation_adjustment(args...)

Legacy wrapper for `saturation_adjustment` with `phase_type` argument.
Forwards to `saturation_adjustment(..., ::ρeq, ...)` with arguments `(ρ, e_int)` swapped.
"""
@inline function saturation_adjustment(
    sat_adjust_method::Type,
    param_set::APS,
    e_int::Real,
    ρ::Real,
    q_tot::Real,
    phase_type::Type,
    maxiter::Int,
    tol::Real,
    T_guess = nothing,
)
    res = saturation_adjustment(
        sat_adjust_method,
        param_set,
        ρe(),
        ρ,
        e_int,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
    T = res.T
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
    res = saturation_adjustment(
        sat_adjust_method,
        param_set,
        ph(),
        p,
        h,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
    T = res.T
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
    res = saturation_adjustment(
        sat_adjust_method,
        param_set,
        pe(),
        p,
        e_int,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
    T = res.T
    return T
end

"""
    saturation_adjustment_given_pθq(args...)

Deprecated. Use [`saturation_adjustment`](@ref) with `pθ_li_q()` dispatch instead.
"""
@inline function saturation_adjustment_given_pθq(
    sat_adjust_method::Type,
    param_set::APS,
    p::Real,
    θ_li::Real,
    q_tot::Real,
    phase_type::Type,
    maxiter::Int,
    tol::Real,
    T_guess = nothing,
)
    res = saturation_adjustment(
        sat_adjust_method,
        param_set,
        pθ_li(),
        p,
        θ_li,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
    T = res.T
    return T
end

"""
    saturation_adjustment_given_ρθq(args...)

Deprecated. Use [`saturation_adjustment`](@ref) with `ρθ_li_q()` dispatch instead.
"""
@inline function saturation_adjustment_given_ρθq(
    param_set::APS,
    ρ::Real,
    θ_li::Real,
    q_tot::Real,
    phase_type::Type,
    maxiter::Int,
    tol,
    T_guess = nothing,
    sat_adjust_method::Type = RS.SecantMethod,
)
    res = saturation_adjustment(
        sat_adjust_method,
        param_set,
        ρθ_li(),
        ρ,
        θ_li,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
    T = res.T
    return T
end

"""
    saturation_adjustment_ρpq(args...)

Deprecated. Use [`saturation_adjustment`](@ref) with `pρq()` dispatch instead.
"""
@inline function saturation_adjustment_ρpq(
    sat_adjust_method::Type,
    param_set::APS,
    ρ::Real,
    p::Real,
    q_tot::Real,
    phase_type::Type,
    maxiter::Int,
    tol::Real,
    T_guess = nothing,
)
    res = saturation_adjustment(
        sat_adjust_method,
        param_set,
        pρ(),
        p,
        ρ,
        q_tot,
        maxiter,
        tol,
        T_guess,
    )
    T = res.T
    return T
end
