using Test
using JET
using Thermodynamics
import ClimaParams as CP
import Thermodynamics.Parameters as TP
import RootSolvers as RS
import Thermodynamics as TD

include(joinpath(@__DIR__, "TestedProfiles.jl"))
using .TestedProfiles

# Force compiled pruning of warning paths for JET
TD.print_warning() = false

#####
##### Input Generation (in part replicated from perf/common_micro_bm.jl)
#####

function functional_inputs(param_set, ::Type{FT}; n = 1500) where {FT}
    ps = TestedProfiles.PhaseEquilProfiles(param_set, Array{FT})
    n = min(n, length(ps.z))
    sl = 1:n

    T = ps.T[sl]
    ρ = ps.ρ[sl]
    p = ps.p[sl]
    q_tot = ps.q_tot[sl]
    q_liq = ps.q_liq[sl]
    q_ice = ps.q_ice[sl]

    e_int_ρ = TD.internal_energy_sat.(Ref(param_set), T, ρ, q_tot)
    θ_ρ = TD.liquid_ice_pottemp.(Ref(param_set), T, ρ, q_tot, q_liq, q_ice)

    ρ_p = TD.air_density.(Ref(param_set), T, p, q_tot)
    q_liq_p, q_ice_p = let qi = TD.condensate_partition.(Ref(param_set), T, ρ_p, q_tot)
        (first.(qi), last.(qi))
    end
    e_int_p = TD.internal_energy_sat.(Ref(param_set), T, ρ_p, q_tot)
    h_p = TD.enthalpy_sat.(Ref(param_set), T, ρ_p, q_tot)
    θ_p =
        TD.liquid_ice_pottemp_given_pressure.(Ref(param_set), T, p, q_tot, q_liq_p, q_ice_p)

    return (;
        T,
        ρ,
        p,
        q_tot,
        q_liq,
        q_ice,
        e_int_ρ,
        θ_ρ,
        e_int_p,
        h_p,
        θ_p,
        q_liq_p,
        q_ice_p,
    )
end

function find_dry(inputs)
    i = findfirst(==(zero(eltype(inputs.q_tot))), inputs.q_tot)
    isnothing(i) && error("Dry index not found")
    return i
end

function find_moist_index(inputs)
    i = findfirst(>(0.01), inputs.q_tot)
    isnothing(i) && error("Moist index not found")
    return i
end

function find_saturated_index(inputs; use_p_based::Bool)
    q_liq = use_p_based ? inputs.q_liq_p : inputs.q_liq
    q_ice = use_p_based ? inputs.q_ice_p : inputs.q_ice
    i = findfirst(i -> (q_liq[i] + q_ice[i]) > 0, eachindex(q_liq))
    isnothing(i) && error("Saturated index not found")
    return i
end

use_p_based(::Type{TD.pe}) = true
use_p_based(::Type{TD.ph}) = true
use_p_based(::Type{TD.pθ_li}) = true
use_p_based(::Type) = false

function conditions_index(inputs, sym, ftype::Type)
    if sym == :dry
        return find_dry(inputs)
    elseif sym == :moist
        return find_moist_index(inputs)
    elseif sym == :saturated
        return find_saturated_index(inputs; use_p_based = use_p_based(ftype))
    else
        error("Bad sym given")
    end
end

get_kwargs(x, ::Type{TD.ρe}) = (; ρ = x.ρ, e_int = x.e_int_ρ, q_tot = x.q_tot)
get_kwargs(x, ::Type{TD.pe}) = (; p = x.p, e_int = x.e_int_p, q_tot = x.q_tot)
get_kwargs(x, ::Type{TD.ph}) = (; p = x.p, h = x.h_p, q_tot = x.q_tot)
get_kwargs(x, ::Type{TD.pρ}) = (; p = x.p, ρ = x.ρ, q_tot = x.q_tot)
get_kwargs(x, ::Type{TD.ρθ_li}) = (; ρ = x.ρ, θ_liq_ice = x.θ_ρ, q_tot = x.q_tot)
get_kwargs(x, ::Type{TD.pθ_li}) = (; p = x.p, θ_liq_ice = x.θ_p, q_tot = x.q_tot)

conditions(::Type) = (:dry, :saturated)

function sample_args(inputs, param_set, sym, ftype::Type)
    i = conditions_index(inputs, sym, ftype)
    kwargs = get_kwargs(inputs, ftype)
    return getindex.(values(kwargs), i)
end

solver_for(::Type{TD.ρe}) = RS.NewtonsMethod
solver_for(::Type) = RS.SecantMethod

#####
##### JET Test Logic
#####

@noinline function check_saturation_adjustment(
    solver,
    param_set,
    method,
    args,
    maxiter,
    tol,
)
    TD.saturation_adjustment(solver, param_set, method, args..., maxiter, tol)
end

function jet_thermo_states(::Type{FT}) where {FT}
    param_set = TP.ThermodynamicsParameters(FT)
    inputs = functional_inputs(param_set, FT)

    @testset "JET optimization tests" begin
        # Test for allocation-free / type-stable saturation adjustment
        for F in (TD.ρe, TD.pe, TD.ph, TD.pρ, TD.ρθ_li, TD.pθ_li)
            for cond in conditions(F)
                args = sample_args(inputs, param_set, cond, F)
                solver = solver_for(F)
                JET.@test_opt check_saturation_adjustment(
                    solver,
                    param_set,
                    F(),
                    args,
                    40,
                    FT(1e-10),
                )
            end
        end
    end
end

# Run the tests
jet_thermo_states(Float32)
