# These functions (variants of sa_numerical_method)
# return an instance of a numerical method to solve
# saturation adjustment, for different combinations
# of thermodynamic variable inputs.

# KA.@print only accepts literal strings, so we must
# branch to print which method is being used.
function print_numerical_method(
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

#####
##### Thermodynamic variable inputs: ρ, e_int, q_tot
#####
function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    e_int::FT,
    q_tot::FT,
    phase_type::Type{<:PhaseEquil},
) where {FT, NM <: RS.NewtonsMethod}
    T_min::FT = CPP.T_min(param_set)
    T_init =
        max(T_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    return RS.NewtonsMethod(
        T_init,
        T_ -> ∂e_int_∂T(param_set, heavisided(T_), e_int, ρ, q_tot, phase_type),
    )
end

function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    e_int::FT,
    q_tot::FT,
    phase_type::Type{<:PhaseEquil},
) where {FT, NM <: RS.NewtonsMethodAD}
    T_min::FT = CPP.T_min(param_set)
    T_init =
        max(T_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    return RS.NewtonsMethodAD(T_init)
end

function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    e_int::FT,
    q_tot::FT,
    phase_type::Type{<:PhaseEquil},
) where {FT, NM <: RS.SecantMethod}
    T_min::FT = CPP.T_min(param_set)
    q_pt = PhasePartition(q_tot, FT(0), q_tot) # Assume all ice
    T_2 = air_temperature(param_set, e_int, q_pt)
    T_1 = max(T_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    T_2 = bound_upper_temperature(T_1, T_2)
    return RS.SecantMethod(T_1, T_2)
end

function sa_numerical_method(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    e_int::FT,
    q_tot::FT,
    phase_type::Type{<:PhaseEquil},
) where {FT, NM <: RS.RegulaFalsiMethod}
    T_min::FT = CPP.T_min(param_set)
    q_pt = PhasePartition(q_tot, FT(0), q_tot) # Assume all ice
    T_2 = air_temperature(param_set, e_int, q_pt)
    T_1 = max(T_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    T_2 = bound_upper_temperature(T_1, T_2)
    return RS.RegulaFalsiMethod(T_1, T_2)
end

#####
##### Thermodynamic variable inputs: ρ, p, q_tot
#####

function sa_numerical_method_ρpq(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    p::FT,
    q_tot::FT,
    phase_type::Type{<:PhaseEquil},
) where {FT, NM <: RS.NewtonsMethodAD}
    q_pt = PhasePartition(q_tot)
    T_init = air_temperature_from_ideal_gas_law(param_set, p, ρ, q_pt)
    return RS.NewtonsMethodAD(T_init)
end

function sa_numerical_method_ρpq(
    ::Type{NM},
    param_set::APS,
    ρ::FT,
    p::FT,
    q_tot::FT,
    phase_type::Type{<:PhaseEquil},
) where {FT, NM <: RS.RegulaFalsiMethod}
    q_pt = PhasePartition(q_tot)
    T_1 = air_temperature_from_ideal_gas_law(param_set, p, ρ, q_pt) - 5
    T_2 = air_temperature_from_ideal_gas_law(param_set, p, ρ, q_pt) + 5
    return RS.RegulaFalsiMethod(T_1, T_2)
end

#####
##### Thermodynamic variable inputs: p, e_int, q_tot
#####

function sa_numerical_method_peq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    e_int::FT,
    q_tot::FT,
    phase_type::Type{<:PhaseEquil},
) where {FT, NM <: RS.NewtonsMethodAD}
    T_min::FT = CPP.T_min(param_set)
    T_init =
        max(T_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    return RS.NewtonsMethodAD(T_init)
end

function sa_numerical_method_peq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    e_int::FT,
    q_tot::FT,
    phase_type::Type{<:PhaseEquil},
) where {FT, NM <: RS.SecantMethod}
    T_min::FT = CPP.T_min(param_set)
    q_pt = PhasePartition(q_tot, FT(0), q_tot) # Assume all ice
    T_2 = air_temperature(param_set, e_int, q_pt)
    T_1 = max(T_min, air_temperature(param_set, e_int, PhasePartition(q_tot))) # Assume all vapor
    T_2 = bound_upper_temperature(T_1, T_2)
    return RS.SecantMethod(T_1, T_2)
end

#####
##### Thermodynamic variable inputs: p, θ_liq_ice, q_tot
#####

function sa_numerical_method_pθq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    phase_type::Type{<:PhaseEquil},
) where {FT, NM <: RS.RegulaFalsiMethod}
    _T_min::FT = CPP.T_min(param_set)
    _T_max::FT = CPP.T_max(param_set)
    air_temp(q) = air_temperature_given_pθq(param_set, p, θ_liq_ice, q)
    T_1 = max(_T_min, air_temp(PhasePartition(q_tot))) # Assume all vapor
    T_2 = T_1 + 10
    T_1 = T_1 - 10
    return RS.RegulaFalsiMethod(T_1, T_2)
end

function sa_numerical_method_pθq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    phase_type::Type{<:PhaseEquil},
) where {FT, NM <: RS.SecantMethod}
    _T_min::FT = CPP.T_min(param_set)
    air_temp(q) = air_temperature_given_pθq(param_set, p, θ_liq_ice, q)
    T_1 = max(_T_min, air_temp(PhasePartition(q_tot))) # Assume all vapor
    T_2 = air_temp(PhasePartition(q_tot, FT(0), q_tot)) # Assume all ice
    T_2 = bound_upper_temperature(T_1, T_2)
    return RS.SecantMethod(T_1, T_2)
end

function sa_numerical_method_pθq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    phase_type::Type{<:PhaseEquil},
) where {FT, NM <: RS.NewtonsMethodAD}
    T_min::FT = CPP.T_min(param_set)
    air_temp(q) = air_temperature_given_pθq(param_set, p, θ_liq_ice, q)
    T_init = max(T_min, air_temp(PhasePartition(q_tot))) # Assume all vapor
    return RS.NewtonsMethodAD(T_init)
end
