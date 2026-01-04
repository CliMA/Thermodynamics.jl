# Configures saturation adjustment numerical methods

#####
##### Thermodynamic variable inputs: ρ, e_int, q_tot
#####
@inline function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ρ,
    e_int,
    q_tot,
    T_guess,
) where {NM <: Union{RS.NewtonsMethod, RS.NewtonsMethodAD}}
    T_init_min = TP.T_init_min(param_set)
    T_init = if T_guess isa Nothing
        max(T_init_min, air_temperature(param_set, e_int, q_tot)) # Assume all vapor
    else
        T_guess
    end
    return NM(T_init)
end

@inline function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    e_int,
    q_tot,
    T_guess,
) where {FT, NM <: RS.SecantMethod}
    T_init_min = TP.T_init_min(param_set)
    if T_guess isa Nothing
        # Physical bracketing: all-vapor and all-ice bounds
        T_1 = max(T_init_min, air_temperature(param_set, e_int, q_tot))
        T_2 = air_temperature(param_set, e_int, q_tot, FT(0), q_tot)
    else
        # Use T_guess as lower bracket, all-ice as upper
        T_1 = max(T_init_min, T_guess)
        T_2 = air_temperature(param_set, e_int, q_tot, FT(0), q_tot)
    end
    T_2 = bound_upper_temperature(param_set, T_1, T_2)
    return RS.SecantMethod(T_1, T_2)
end

@inline function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    e_int,
    q_tot,
    T_guess,
) where {FT, NM <: RS.BrentsMethod}
    T_init_min = TP.T_init_min(param_set)
    # BrentsMethod requires strict bracketing - ignore T_guess
    # Physical bracketing: all-vapor and all-ice bounds
    T_1 = max(T_init_min, air_temperature(param_set, e_int, q_tot))
    T_2 = air_temperature(param_set, e_int, q_tot, FT(0), q_tot)
    T_2 = bound_upper_temperature(param_set, T_1, T_2)
    return RS.BrentsMethod(T_1, T_2)
end

#####
##### Thermodynamic variable inputs: ρ, p, q_tot
#####

@inline function sa_numerical_method_ρpq(
    ::Type{NM},
    param_set::APS,
    ρ,
    p,
    q_tot,
    T_guess,
) where {NM <: RS.NewtonsMethodAD}
    # Use scalar q_tot, assume all vapor for init
    T_init = if T_guess isa Nothing
        air_temperature_given_pρq(param_set, p, ρ, q_tot)
    else
        T_guess
    end
    return RS.NewtonsMethodAD(T_init)
end

@inline function sa_numerical_method_ρpq(
    ::Type{NM},
    param_set::APS,
    ρ,
    p,
    q_tot,
    T_guess,
) where {NM <: RS.BrentsMethod}
    T_init_min = TP.T_init_min(param_set)
    # BrentsMethod requires strict bracketing - ignore T_guess
    # Physical bracketing: all-vapor assumption for lower bound
    T_1 = max(T_init_min, air_temperature_given_pρq(param_set, p, ρ, q_tot))
    # Upper bound: conservative +10K for stability
    T_2 = T_1 + 10
    return RS.BrentsMethod(T_1, T_2)
end

#####
##### Thermodynamic variable inputs: p, e_int, q_tot
#####

@inline function sa_numerical_method_peq(
    ::Type{NM},
    param_set::APS,
    p,
    e_int,
    q_tot,
    T_guess,
) where {NM <: RS.NewtonsMethodAD}
    T_init_min = TP.T_init_min(param_set)
    T_init = if T_guess isa Nothing
        max(T_init_min, air_temperature(param_set, e_int, q_tot)) # Assume all vapor
    else
        T_guess
    end
    return RS.NewtonsMethodAD(T_init)
end

@inline function sa_numerical_method_peq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    e_int,
    q_tot,
    T_guess,
) where {FT, NM <: RS.SecantMethod}
    T_init_min = TP.T_init_min(param_set)
    if T_guess isa Nothing
        T_1 = max(T_init_min, air_temperature(param_set, e_int, q_tot))
        T_2 = air_temperature(param_set, e_int, q_tot, FT(0), q_tot)
    else
        T_1 = max(T_init_min, T_guess)
        T_2 = air_temperature(param_set, e_int, q_tot, FT(0), q_tot)
    end
    T_2 = bound_upper_temperature(param_set, T_1, T_2)
    return RS.SecantMethod(T_1, T_2)
end

#####
##### Thermodynamic variable inputs: p, h, q_tot
#####

@inline function sa_numerical_method_phq(
    ::Type{NM},
    param_set::APS,
    p,
    h,
    q_tot,
    T_guess,
) where {NM <: RS.NewtonsMethodAD}
    T_init_min = TP.T_init_min(param_set)
    T_init = if T_guess isa Nothing # Assume all vapor
        max(
            T_init_min,
            air_temperature_given_hq(param_set, h, q_tot),
        )
    else
        T_guess
    end
    return RS.NewtonsMethodAD(T_init)
end

@inline function sa_numerical_method_phq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    h,
    q_tot,
    T_guess,
) where {FT, NM <: RS.SecantMethod}
    T_init_min = TP.T_init_min(param_set)
    if T_guess isa Nothing
        T_1 = max(T_init_min, air_temperature_given_hq(param_set, h, q_tot))
        T_2 = air_temperature_given_hq(param_set, h, q_tot, FT(0), q_tot) # Assume all ice
    else
        T_1 = max(T_init_min, T_guess)
        T_2 = air_temperature_given_hq(param_set, h, q_tot, FT(0), q_tot)
    end
    T_2 = bound_upper_temperature(param_set, T_1, T_2)
    return RS.SecantMethod(T_1, T_2)
end

@inline function sa_numerical_method_phq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    h,
    q_tot,
    T_guess,
) where {FT, NM <: RS.BrentsMethod}
    T_init_min = TP.T_init_min(param_set)
    # BrentsMethod requires strict bracketing - ignore T_guess
    T_2 = air_temperature_given_hq(param_set, h, q_tot, FT(0), q_tot) # Assume all ice
    T_1 = max(
        T_init_min,
        air_temperature_given_hq(param_set, h, q_tot),
    ) # Assume all vapor
    T_2 = bound_upper_temperature(param_set, T_1, T_2)
    return RS.BrentsMethod(T_1, T_2)
end

#####
##### Thermodynamic variable inputs: p, θ_liq_ice, q_tot
#####

@inline function sa_numerical_method_pθq(
    ::Type{NM},
    param_set::APS,
    p,
    θ_liq_ice,
    q_tot,
    T_guess,
) where {NM <: RS.BrentsMethod}
    T_init_min = TP.T_init_min(param_set)
    # BrentsMethod requires strict bracketing - ignore T_guess
    @inline air_temp(args...) = air_temperature_given_pθq(param_set, p, θ_liq_ice, args...)
    T_1 = max(T_init_min, air_temp(q_tot)) # Assume all vapor
    T_2 = T_1 + 10
    T_1 = T_1 - 10
    return RS.BrentsMethod(T_1, T_2)
end

@inline function sa_numerical_method_pθq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    θ_liq_ice,
    q_tot,
    T_guess,
) where {FT, NM <: RS.SecantMethod}
    T_init_min = TP.T_init_min(param_set)
    @inline air_temp(args...) = air_temperature_given_pθq(param_set, p, θ_liq_ice, args...)
    if T_guess isa Nothing
        T_1 = max(T_init_min, air_temp(q_tot))
        T_2 = air_temp(q_tot, FT(0), q_tot)
    else
        T_1 = max(T_init_min, T_guess)
        T_2 = air_temp(q_tot, FT(0), q_tot)
    end
    T_2 = bound_upper_temperature(param_set, T_1, T_2)
    return RS.SecantMethod(T_1, T_2)
end

@inline function sa_numerical_method_pθq(
    ::Type{NM},
    param_set::APS,
    p,
    θ_liq_ice,
    q_tot,
    T_guess,
) where {NM <: RS.NewtonsMethodAD}
    T_init_min = TP.T_init_min(param_set)
    @inline air_temp(args...) = air_temperature_given_pθq(param_set, p, θ_liq_ice, args...)
    T_init = if T_guess isa Nothing
        max(T_init_min, air_temp(q_tot)) # Assume all vapor
    else
        T_guess
    end
    return RS.NewtonsMethodAD(T_init)
end

#####
##### Thermodynamic variable inputs: ρ, θ_liq_ice, q_tot
#####

@inline function sa_numerical_method_ρθq(
    ::Type{NM},
    param_set::APS,
    ρ,
    θ_liq_ice,
    q_tot,
    T_guess,
) where {NM <: RS.BrentsMethod}
    T_init_min = TP.T_init_min(param_set)
    # BrentsMethod requires strict bracketing - ignore T_guess
    @inline air_temp(args...) = air_temperature_given_ρθq(param_set, ρ, θ_liq_ice, args...)
    T_1 = max(T_init_min, air_temp(q_tot)) # Assume all vapor
    T_2 = T_1 + 10
    T_1 = T_1 - 10
    return RS.BrentsMethod(T_1, T_2)
end

@inline function sa_numerical_method_ρθq(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    θ_liq_ice,
    q_tot,
    T_guess,
) where {FT, NM <: RS.SecantMethod}
    T_init_min = TP.T_init_min(param_set)
    @inline air_temp(args...) = air_temperature_given_ρθq(param_set, ρ, θ_liq_ice, args...)
    if T_guess isa Nothing
        T_1 = max(T_init_min, air_temp(q_tot))
        T_2 = air_temp(q_tot, FT(0), q_tot)
    else
        T_1 = max(T_init_min, T_guess)
        T_2 = air_temp(q_tot, FT(0), q_tot)
    end
    T_2 = bound_upper_temperature(param_set, T_1, T_2)
    return RS.SecantMethod(T_1, T_2)
end

@inline function sa_numerical_method_ρθq(
    ::Type{NM},
    param_set::APS,
    ρ,
    θ_liq_ice,
    q_tot,
    T_guess,
) where {NM <: RS.NewtonsMethodAD}
    T_init_min = TP.T_init_min(param_set)
    @inline air_temp(args...) = air_temperature_given_ρθq(param_set, ρ, θ_liq_ice, args...)
    T_init = if T_guess isa Nothing
        max(T_init_min, air_temp(q_tot)) # Assume all vapor
    else
        T_guess
    end
    return RS.NewtonsMethodAD(T_init)
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