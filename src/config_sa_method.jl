# Configures saturation adjustment numerical methods

#####
##### Common solver logic
#####

@inline function _make_sa_solver(
    ::Type{NM},
    param_set::APS,
    T_unsat,
    T_ice,
    T_guess,
) where {NM <: Union{RS.NewtonsMethod, RS.NewtonsMethodAD}}
    T_init_min = TP.T_init_min(param_set)
    T_init = if T_guess isa Nothing
        max(T_init_min, T_unsat)
    else
        T_guess
    end
    return NM(T_init)
end

@inline function _make_sa_solver(
    ::Type{NM},
    param_set::APS,
    T_unsat,
    T_ice,
    T_guess,
) where {NM <: RS.SecantMethod}
    T_init_min = TP.T_init_min(param_set)
    if T_guess isa Nothing
        T_1 = max(T_init_min, T_unsat)
    else
        T_1 = max(T_init_min, T_guess)
    end
    T_2 = bound_upper_temperature(param_set, T_1, T_ice)
    return RS.SecantMethod(T_1, T_2)
end

@inline function _make_sa_solver(
    ::Type{NM},
    param_set::APS,
    T_unsat,
    T_ice,
    T_guess,
) where {NM <: RS.BrentsMethod}
    T_init_min = TP.T_init_min(param_set)
    # BrentsMethod requires strict bracketing - ignore T_guess
    T_1 = max(T_init_min, T_unsat)
    T_2 = bound_upper_temperature(param_set, T_1, T_ice)
    return RS.BrentsMethod(T_1, T_2)
end

#####
##### Thermodynamic variable inputs: ρ, e_int, q_tot
#####

@inline function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ::ρeq,
    ρ,
    e_int,
    q_tot,
    T_guess,
) where {NM}
    FT = eltype(param_set)
    T_unsat = air_temperature(param_set, ρeq(), e_int, q_tot)
    T_ice = air_temperature(param_set, ρeq(), e_int, q_tot, FT(0), q_tot)
    return _make_sa_solver(NM, param_set, T_unsat, T_ice, T_guess)
end

#####
##### Thermodynamic variable inputs: ρ, p, q_tot
#####

# Takes p, ρ order to match IndepVars and air_temperature signature
@inline function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ::pρq,
    p,
    ρ,
    q_tot,
    T_guess,
) where {NM}
    FT = eltype(param_set)
    T_unsat = air_temperature(param_set, pρq(), p, ρ, q_tot)
    T_ice = air_temperature(param_set, pρq(), p, ρ, q_tot, FT(0), q_tot)
    return _make_sa_solver(NM, param_set, T_unsat, T_ice, T_guess)
end

#####
##### Thermodynamic variable inputs: p, e_int, q_tot
#####

@inline function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ::peq,
    p,
    e_int,
    q_tot,
    T_guess,
) where {NM}
    FT = eltype(param_set)
    T_unsat = air_temperature(param_set, e_int, q_tot)
    T_ice = air_temperature(param_set, e_int, q_tot, FT(0), q_tot)
    return _make_sa_solver(NM, param_set, T_unsat, T_ice, T_guess)
end


#####
##### Thermodynamic variable inputs: p, h, q_tot
#####

@inline function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ::phq,
    p,
    h,
    q_tot,
    T_guess,
) where {NM}
    FT = eltype(param_set)
    T_unsat = air_temperature(param_set, phq(), h, q_tot)
    T_ice = air_temperature(param_set, phq(), h, q_tot, FT(0), q_tot)
    return _make_sa_solver(NM, param_set, T_unsat, T_ice, T_guess)
end

#####
##### Thermodynamic variable inputs: p, θ_liq_ice, q_tot
#####

@inline function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ::pθ_liq_ice_q,
    p,
    θ_liq_ice,
    q_tot,
    T_guess,
) where {NM}
    FT = eltype(param_set)
    T_unsat = air_temperature(param_set, pθ_liq_ice_q(), p, θ_liq_ice, q_tot)
    T_ice = air_temperature(
        param_set,
        pθ_liq_ice_q(),
        p,
        θ_liq_ice,
        q_tot,
        FT(0),
        q_tot,
    )
    return _make_sa_solver(NM, param_set, T_unsat, T_ice, T_guess)
end

#####
##### Thermodynamic variable inputs: ρ, θ_liq_ice, q_tot
#####

@inline function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ::ρθ_liq_ice_q,
    ρ,
    θ_liq_ice,
    q_tot,
    T_guess,
) where {NM}
    FT = eltype(param_set)
    T_unsat = air_temperature(param_set, ρθ_liq_ice_q(), ρ, θ_liq_ice, q_tot)
    T_ice = air_temperature(
        param_set,
        ρθ_liq_ice_q(),
        ρ,
        θ_liq_ice,
        q_tot,
        FT(0),
        q_tot,
    )
    return _make_sa_solver(NM, param_set, T_unsat, T_ice, T_guess)
end

"""
    bound_upper_temperature(param_set, T_1, T_2)

Internal function. Bounds the upper temperature guess `T_2` for bracket methods
to ensure numerical stability and physical validity.

Returns T_2 bounded to [T_1 + 3K, min(T_max, T_1 + 10K)] to ensure:
- Minimum bracket width of 3K for numerical stability
- Maximum bracket width of 10K to maintain accuracy
- Never exceeds physical T_max
"""
@inline function bound_upper_temperature(param_set, T_1, T_2)
    T_max = TP.T_max(param_set)
    T_2_bounded = max(T_1 + 3, T_2)
    return min(min(T_max, T_1 + 10), T_2_bounded)
end
