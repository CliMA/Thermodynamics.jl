# Configures saturation adjustment numerical methods

#####
##### Common solver logic
#####

@inline function _make_sa_solver(
    method_selector,
    param_set::APS,
    T_unsat,
    T_ice,
    T_guess,
)
    T_init_min = TP.T_init_min(param_set)

    if method_selector isa RootSolvers.NewtonsSelector
        T_init = T_guess isa Nothing ? max(T_init_min, T_unsat) : T_guess
        return RootSolvers.NewtonsMethod(T_init)
    elseif method_selector isa RootSolvers.NewtonsADSelector
        T_init = T_guess isa Nothing ? max(T_init_min, T_unsat) : T_guess
        return RootSolvers.NewtonsMethodAD(T_init)
    elseif method_selector isa RootSolvers.SecantSelector
        T_lo = T_guess isa Nothing ? max(T_init_min, T_unsat) : max(T_init_min, T_guess)
        T_hi = bound_upper_temperature(param_set, T_lo, T_ice)
        return RootSolvers.SecantMethod(T_lo, T_hi)
    elseif method_selector isa RootSolvers.BrentsSelector
        # BrentsMethod requires strict bracketing - ignore T_guess
        T_lo = max(T_init_min, T_unsat)
        T_hi = bound_upper_temperature(param_set, T_lo, T_ice)
        return RootSolvers.BrentsMethod(T_lo, T_hi)
    else
        error("Unknown method selector type: $(typeof(method_selector))")
    end
end

#####
##### Thermodynamic variable inputs: ρ, e_int, q_tot
#####

@inline function sa_numerical_method(
    method_selector,  # RootSolvers.AbstractMethodSelector
    param_set::APS,
    ::ρe,
    ρ,
    e_int,
    q_tot,
    T_guess,
)
    FT = eltype(param_set)
    T_unsat = air_temperature(param_set, ρe(), e_int, q_tot)
    T_ice = air_temperature(param_set, ρe(), e_int, q_tot, FT(0), q_tot)
    return _make_sa_solver(method_selector, param_set, T_unsat, T_ice, T_guess)
end

#####
##### Thermodynamic variable inputs: ρ, p, q_tot
#####

# Takes p, ρ order to match IndepVars and air_temperature signature
@inline function sa_numerical_method(
    method_selector,  # RootSolvers.AbstractMethodSelector
    param_set::APS,
    ::pρ,
    p,
    ρ,
    q_tot,
    T_guess,
)
    FT = eltype(param_set)
    T_unsat = air_temperature(param_set, pρ(), p, ρ, q_tot)
    T_ice = air_temperature(param_set, pρ(), p, ρ, q_tot, FT(0), q_tot)
    return _make_sa_solver(method_selector, param_set, T_unsat, T_ice, T_guess)
end

#####
##### Thermodynamic variable inputs: p, e_int, q_tot
#####

@inline function sa_numerical_method(
    method_selector,  # RootSolvers.AbstractMethodSelector
    param_set::APS,
    ::pe,
    p,
    e_int,
    q_tot,
    T_guess,
)
    FT = eltype(param_set)
    T_unsat = air_temperature(param_set, e_int, q_tot)
    T_ice = air_temperature(param_set, e_int, q_tot, FT(0), q_tot)
    return _make_sa_solver(method_selector, param_set, T_unsat, T_ice, T_guess)
end


#####
##### Thermodynamic variable inputs: p, h, q_tot
#####

@inline function sa_numerical_method(
    method_selector,  # RootSolvers.AbstractMethodSelector
    param_set::APS,
    ::ph,
    p,
    h,
    q_tot,
    T_guess,
)
    FT = eltype(param_set)
    T_unsat = air_temperature(param_set, ph(), h, q_tot)
    T_ice = air_temperature(param_set, ph(), h, q_tot, FT(0), q_tot)
    return _make_sa_solver(method_selector, param_set, T_unsat, T_ice, T_guess)
end

#####
##### Thermodynamic variable inputs: p, θ_liq_ice, q_tot
#####

@inline function sa_numerical_method(
    method_selector,  # RootSolvers.AbstractMethodSelector
    param_set::APS,
    ::pθ_li,
    p,
    θ_liq_ice,
    q_tot,
    T_guess,
)
    FT = eltype(param_set)
    T_unsat = air_temperature(param_set, pθ_li(), p, θ_liq_ice, q_tot)
    T_ice = air_temperature(
        param_set,
        pθ_li(),
        p,
        θ_liq_ice,
        q_tot,
        FT(0),
        q_tot,
    )
    return _make_sa_solver(method_selector, param_set, T_unsat, T_ice, T_guess)
end

#####
##### Thermodynamic variable inputs: ρ, θ_liq_ice, q_tot
#####

@inline function sa_numerical_method(
    method_selector,  # RootSolvers.AbstractMethodSelector
    param_set::APS,
    ::ρθ_li,
    ρ,
    θ_liq_ice,
    q_tot,
    T_guess,
)
    FT = eltype(param_set)
    T_unsat = air_temperature(param_set, ρθ_li(), ρ, θ_liq_ice, q_tot)
    T_ice = air_temperature(
        param_set,
        ρθ_li(),
        ρ,
        θ_liq_ice,
        q_tot,
        FT(0),
        q_tot,
    )
    return _make_sa_solver(method_selector, param_set, T_unsat, T_ice, T_guess)
end

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
    # Ensure T_hi > T_lo for numerical initialization (bracket width / finite difference)
    return max(T_lo + FT(0.1), T_hi_phys)
end
