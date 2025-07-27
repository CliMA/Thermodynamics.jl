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
    ::Type{phase_type},
    T_guess,
) where {NM <: RS.NewtonsMethod, phase_type <: PhaseEquil}
    T_init_min = TP.T_init_min(param_set)
    T_init = if T_guess isa Nothing
        max(T_init_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    else
        T_guess
    end
    return RS.NewtonsMethod(T_init)
end

@inline function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ρ,
    e_int,
    q_tot,
    ::Type{phase_type},
    T_guess,
) where {NM <: RS.NewtonsMethodAD, phase_type <: PhaseEquil}
    T_init_min = TP.T_init_min(param_set)
    T_init = if T_guess isa Nothing
        max(T_init_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    else
        T_guess
    end
    return RS.NewtonsMethodAD(T_init)
end

@inline function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    e_int,
    q_tot,
    ::Type{phase_type},
    T_guess,
) where {FT, NM <: RS.SecantMethod, phase_type <: PhaseEquil}
    T_init_min = TP.T_init_min(param_set)
    q_pt = PhasePartition(q_tot, FT(0), q_tot) # Assume all ice
    T_2 = air_temperature(param_set, e_int, q_pt)
    T_1 = max(
        T_init_min,
        air_temperature(param_set, e_int, PhasePartition(q_tot)),
    ) # Assume all vapor
    T_2 = bound_upper_temperature(T_1, T_2)
    return RS.SecantMethod(T_1, T_2)
end

@inline function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    e_int,
    q_tot,
    ::Type{phase_type},
    T_guess,
) where {FT, NM <: RS.RegulaFalsiMethod, phase_type <: PhaseEquil}
    T_init_min = TP.T_init_min(param_set)
    q_pt = PhasePartition(q_tot, FT(0), q_tot) # Assume all ice
    T_2 = air_temperature(param_set, e_int, q_pt)
    T_1 = max(
        T_init_min,
        air_temperature(param_set, e_int, PhasePartition(q_tot)),
    ) # Assume all vapor
    T_2 = bound_upper_temperature(T_1, T_2)
    return RS.RegulaFalsiMethod(T_1, T_2)
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
    ::Type{phase_type},
    T_guess,
) where {NM <: RS.NewtonsMethodAD, phase_type <: PhaseEquil}
    q_pt = PhasePartition(q_tot)
    T_init = if T_guess isa Nothing
        air_temperature_given_ρp(param_set, p, ρ, q_pt)
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
    ::Type{phase_type},
    T_guess,
) where {NM <: RS.RegulaFalsiMethod, phase_type <: PhaseEquil}
    q_pt = PhasePartition(q_tot)
    T_1 = air_temperature_given_ρp(param_set, p, ρ, q_pt) - 5
    T_2 = air_temperature_given_ρp(param_set, p, ρ, q_pt) + 5
    return RS.RegulaFalsiMethod(T_1, T_2)
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
    ::Type{phase_type},
    T_guess,
) where {NM <: RS.NewtonsMethodAD, phase_type <: PhaseEquil}
    T_init_min = TP.T_init_min(param_set)
    T_init = if T_guess isa Nothing
        max(T_init_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
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
    ::Type{phase_type},
    T_guess,
) where {FT, NM <: RS.SecantMethod, phase_type <: PhaseEquil}
    T_init_min = TP.T_init_min(param_set)
    q_pt = PhasePartition(q_tot, FT(0), q_tot) # Assume all ice
    T_2 = air_temperature(param_set, e_int, q_pt)
    T_1 = max(
        T_init_min,
        air_temperature(param_set, e_int, PhasePartition(q_tot)),
    ) # Assume all vapor
    T_2 = bound_upper_temperature(T_1, T_2)
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
    ::Type{phase_type},
    T_guess,
) where {NM <: RS.NewtonsMethodAD, phase_type <: PhaseEquil}
    T_init_min = TP.T_init_min(param_set)
    T_init = if T_guess isa Nothing # Assume all vapor
        max(
            T_init_min,
            air_temperature_from_enthalpy(param_set, h, PhasePartition(q_tot)),
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
    ::Type{phase_type},
    T_guess,
) where {FT, NM <: RS.SecantMethod, phase_type <: PhaseEquil}
    T_init_min = TP.T_init_min(param_set)
    q_pt = PhasePartition(q_tot, FT(0), q_tot) # Assume all ice
    T_2 = air_temperature_from_enthalpy(param_set, h, q_pt)
    T_1 = max(
        T_init_min,
        air_temperature_from_enthalpy(param_set, h, PhasePartition(q_tot)),
    ) # Assume all vapor
    T_2 = bound_upper_temperature(T_1, T_2)
    return RS.SecantMethod(T_1, T_2)
end

@inline function sa_numerical_method_phq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    h,
    q_tot,
    ::Type{phase_type},
    T_guess,
) where {FT, NM <: RS.RegulaFalsiMethod, phase_type <: PhaseEquil}
    T_init_min = TP.T_init_min(param_set)
    q_pt = PhasePartition(q_tot, FT(0), q_tot) # Assume all ice
    T_2 = air_temperature_from_enthalpy(param_set, h, q_pt)
    T_1 = max(
        T_init_min,
        air_temperature_from_enthalpy(param_set, h, PhasePartition(q_tot)),
    ) # Assume all vapor
    T_2 = bound_upper_temperature(T_1, T_2)
    return RS.RegulaFalsiMethod(T_1, T_2)
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
    ::Type{phase_type},
    T_guess,
) where {NM <: RS.RegulaFalsiMethod, phase_type <: PhaseEquil}
    _T_init_min = TP.T_init_min(param_set)
    _T_max = TP.T_max(param_set)
    @inline air_temp(q) = air_temperature_given_pθq(param_set, p, θ_liq_ice, q)
    T_1 = max(_T_init_min, air_temp(PhasePartition(q_tot))) # Assume all vapor
    T_2 = T_1 + 10
    T_1 = T_1 - 10
    return RS.RegulaFalsiMethod(T_1, T_2)
end

@inline function sa_numerical_method_pθq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    θ_liq_ice,
    q_tot,
    ::Type{phase_type},
    T_guess,
) where {FT, NM <: RS.SecantMethod, phase_type <: PhaseEquil}
    _T_init_min = TP.T_init_min(param_set)
    @inline air_temp(q) = air_temperature_given_pθq(param_set, p, θ_liq_ice, q)
    T_1 = max(_T_init_min, air_temp(PhasePartition(q_tot))) # Assume all vapor
    T_2 = air_temp(PhasePartition(q_tot, FT(0), q_tot)) # Assume all ice
    T_2 = bound_upper_temperature(T_1, T_2)
    return RS.SecantMethod(T_1, T_2)
end

@inline function sa_numerical_method_pθq(
    ::Type{NM},
    param_set::APS,
    p,
    θ_liq_ice,
    q_tot,
    ::Type{phase_type},
    T_guess,
) where {NM <: RS.NewtonsMethodAD, phase_type <: PhaseEquil}
    T_init_min = TP.T_init_min(param_set)
    @inline air_temp(q) = air_temperature_given_pθq(param_set, p, θ_liq_ice, q)
    T_init = if T_guess isa Nothing
        max(T_init_min, air_temp(PhasePartition(q_tot))) # Assume all vapor
    else
        T_guess
    end
    return RS.NewtonsMethodAD(T_init)
end
