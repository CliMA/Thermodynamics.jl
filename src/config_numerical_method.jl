# These functions (variants of sa_numerical_method)
# return an instance of a numerical method to solve
# saturation adjustment, for different combinations
# of thermodynamic variable inputs.

# KA.@print only accepts literal strings, so we must
# branch to print which method is being used.
@inline function print_numerical_method(
    ::Type{sat_adjust_method},
) where {sat_adjust_method}
    if sat_adjust_method <: RS.NewtonsMethod
        KA.@print("    Method=NewtonsMethod")
    elseif sat_adjust_method <: RS.NewtonsMethodAD
        KA.@print("    Method=NewtonsMethodAD")
    elseif sat_adjust_method <: RS.SecantMethod
        KA.@print("    Method=SecantMethod")
    elseif sat_adjust_method <: RS.RegulaFalsiMethod
        KA.@print("    Method=RegulaFalsiMethod")
    else
        error("Unsupported numerical method")
    end
end

@inline function print_T_guess(
    ::Type{sat_adjust_method},
    T_guess::Real,
) where {sat_adjust_method}
    if sat_adjust_method <: RS.NewtonsMethod
        KA.@print(", T_guess=", T_guess)
    elseif sat_adjust_method <: RS.NewtonsMethodAD
        KA.@print(", T_guess=", T_guess)
    end
end

@inline function print_T_guess(
    ::Type{sat_adjust_method},
    T_guess::Nothing,
) where {sat_adjust_method}
    if sat_adjust_method <: RS.NewtonsMethod
        KA.@print(", T_guess=nothing")
    elseif sat_adjust_method <: RS.NewtonsMethodAD
        KA.@print(", T_guess=nothing")
    end
end

#####
##### Thermodynamic variable inputs: ρ, e_int, q_tot
#####
@inline function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    e_int::FT,
    q_tot::FT,
    ::Type{phase_type},
    T_guess::Union{FT, Nothing},
) where {FT, NM <: RS.NewtonsMethod, phase_type <: PhaseEquil}
    T_min = TP.T_min(param_set)
    T_init = if T_guess isa Nothing
        max(T_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    else
        T_guess
    end
    return RS.NewtonsMethod(T_init)
end

@inline function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    e_int::FT,
    q_tot::FT,
    ::Type{phase_type},
    T_guess::FT,
) where {FT, NM <: RS.NewtonsMethodAD, phase_type <: PhaseEquil}
    T_min = TP.T_min(param_set)
    T_init = if T_guess isa Nothing
        max(T_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    else
        T_guess
    end
    return RS.NewtonsMethodAD(T_init)
end

@inline function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    e_int::FT,
    q_tot::FT,
    ::Type{phase_type},
    T_guess::Union{FT, Nothing},
) where {FT, NM <: RS.SecantMethod, phase_type <: PhaseEquil}
    T_min = TP.T_min(param_set)
    q_pt = PhasePartition(q_tot, FT(0), q_tot) # Assume all ice
    T_2 = air_temperature(param_set, e_int, q_pt)
    T_1 = max(T_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    T_2 = bound_upper_temperature(T_1, T_2)
    return RS.SecantMethod(T_1, T_2)
end

@inline function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    e_int::FT,
    q_tot::FT,
    ::Type{phase_type},
    T_guess::Union{FT, Nothing},
) where {FT, NM <: RS.RegulaFalsiMethod, phase_type <: PhaseEquil}
    T_min = TP.T_min(param_set)
    q_pt = PhasePartition(q_tot, FT(0), q_tot) # Assume all ice
    T_2 = air_temperature(param_set, e_int, q_pt)
    T_1 = max(T_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    T_2 = bound_upper_temperature(T_1, T_2)
    return RS.RegulaFalsiMethod(T_1, T_2)
end

#####
##### Thermodynamic variable inputs: ρ, p, q_tot
#####

@inline function sa_numerical_method_ρpq(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    p::FT,
    q_tot::FT,
    ::Type{phase_type},
    T_guess::Union{FT, Nothing},
) where {FT, NM <: RS.NewtonsMethodAD, phase_type <: PhaseEquil}
    q_pt = PhasePartition(q_tot)
    T_init = if T_guess isa Nothing
        air_temperature_from_ideal_gas_law(param_set, p, ρ, q_pt)
    else
        T_guess
    end
    return RS.NewtonsMethodAD(T_init)
end

@inline function sa_numerical_method_ρpq(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    p::FT,
    q_tot::FT,
    ::Type{phase_type},
    T_guess::Union{FT, Nothing},
) where {FT, NM <: RS.RegulaFalsiMethod, phase_type <: PhaseEquil}
    q_pt = PhasePartition(q_tot)
    T_1 = air_temperature_from_ideal_gas_law(param_set, p, ρ, q_pt) - 5
    T_2 = air_temperature_from_ideal_gas_law(param_set, p, ρ, q_pt) + 5
    return RS.RegulaFalsiMethod(T_1, T_2)
end

#####
##### Thermodynamic variable inputs: p, e_int, q_tot
#####

@inline function sa_numerical_method_peq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    e_int::FT,
    q_tot::FT,
    ::Type{phase_type},
    T_guess::Union{FT, Nothing},
) where {FT, NM <: RS.NewtonsMethodAD, phase_type <: PhaseEquil}
    T_min = TP.T_min(param_set)
    T_init = if T_guess isa Nothing
        max(T_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    else
        T_guess
    end
    return RS.NewtonsMethodAD(T_init)
end

@inline function sa_numerical_method_peq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    e_int::FT,
    q_tot::FT,
    ::Type{phase_type},
    T_guess::Union{FT, Nothing},
) where {FT, NM <: RS.SecantMethod, phase_type <: PhaseEquil}
    T_min = TP.T_min(param_set)
    q_pt = PhasePartition(q_tot, FT(0), q_tot) # Assume all ice
    T_2 = air_temperature(param_set, e_int, q_pt)
    T_1 = max(T_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    T_2 = bound_upper_temperature(T_1, T_2)
    return RS.SecantMethod(T_1, T_2)
end

#####
##### Thermodynamic variable inputs: p, h, q_tot
#####

@inline function sa_numerical_method_phq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    h::FT,
    q_tot::FT,
    ::Type{phase_type},
    T_guess::Union{FT, Nothing},
) where {FT, NM <: RS.NewtonsMethodAD, phase_type <: PhaseEquil}
    T_min = TP.T_min(param_set)
    T_init = if T_guess isa Nothing # Assume all vapor
        max(
            T_min,
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
    h::FT,
    q_tot::FT,
    ::Type{phase_type},
    T_guess::Union{FT, Nothing},
) where {FT, NM <: RS.SecantMethod, phase_type <: PhaseEquil}
    T_min = TP.T_min(param_set)
    q_pt = PhasePartition(q_tot, FT(0), q_tot) # Assume all ice
    T_2 = air_temperature_from_enthalpy(param_set, h, q_pt)
    T_1 = max(
        T_min,
        air_temperature_from_enthalpy(param_set, h, PhasePartition(q_tot)),
    ) # Assume all vapor
    T_2 = bound_upper_temperature(T_1, T_2)
    return RS.SecantMethod(T_1, T_2)
end

@inline function sa_numerical_method_phq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    h::FT,
    q_tot::FT,
    ::Type{phase_type},
    T_guess::Union{FT, Nothing},
) where {FT, NM <: RS.RegulaFalsiMethod, phase_type <: PhaseEquil}
    T_min = TP.T_min(param_set)
    q_pt = PhasePartition(q_tot, FT(0), q_tot) # Assume all ice
    T_2 = air_temperature_from_enthalpy(param_set, h, q_pt)
    T_1 = max(
        T_min,
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
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    ::Type{phase_type},
    T_guess::Union{FT, Nothing},
) where {FT, NM <: RS.RegulaFalsiMethod, phase_type <: PhaseEquil}
    _T_min = TP.T_min(param_set)
    _T_max = TP.T_max(param_set)
    @inline air_temp(q) = air_temperature_given_pθq(param_set, p, θ_liq_ice, q)
    T_1 = max(_T_min, air_temp(PhasePartition(q_tot))) # Assume all vapor
    T_2 = T_1 + 10
    T_1 = T_1 - 10
    return RS.RegulaFalsiMethod(T_1, T_2)
end

@inline function sa_numerical_method_pθq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    ::Type{phase_type},
    T_guess::Union{FT, Nothing},
) where {FT, NM <: RS.SecantMethod, phase_type <: PhaseEquil}
    _T_min = TP.T_min(param_set)
    @inline air_temp(q) = air_temperature_given_pθq(param_set, p, θ_liq_ice, q)
    T_1 = max(_T_min, air_temp(PhasePartition(q_tot))) # Assume all vapor
    T_2 = air_temp(PhasePartition(q_tot, FT(0), q_tot)) # Assume all ice
    T_2 = bound_upper_temperature(T_1, T_2)
    return RS.SecantMethod(T_1, T_2)
end

@inline function sa_numerical_method_pθq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    ::Type{phase_type},
    T_guess::Union{FT, Nothing},
) where {FT, NM <: RS.NewtonsMethodAD, phase_type <: PhaseEquil}
    T_min = TP.T_min(param_set)
    @inline air_temp(q) = air_temperature_given_pθq(param_set, p, θ_liq_ice, q)
    T_init = if T_guess isa Nothing
        max(T_min, air_temp(PhasePartition(q_tot))) # Assume all vapor
    else
        T_guess
    end
    return RS.NewtonsMethodAD(T_init)
end
