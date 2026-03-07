import StatsBase
import PrettyTables
import OrderedCollections
using Test
import BenchmarkTools

import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import ClimaParams as CP
import RootSolvers as RS

include(joinpath(@__DIR__, "..", "test", "TestedProfiles.jl"))
using .TestedProfiles

#####
##### Finding indexes in profiles satisfying certain conditions
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

#=
    conditions_index

Find an index in the given `profiles` that satisfies
some condition. For example, find the index that we're
sure that saturation adjustment is performed for a given
(profiles.p[i], profiles.θ_li[i], profiles.q_tot[i]).
=#
use_p_based(::Type{TD.peq}) = true
use_p_based(::Type{TD.phq}) = true
use_p_based(::Type{TD.pθ_li_q}) = true
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

function sample_args(inputs, param_set, sym, ftype::Type)
    i = conditions_index(inputs, sym, ftype)
    kwargs = get_kwargs(inputs, ftype)
    return getindex.(values(kwargs), i)
end

#####
##### BenchmarkTools's trial utils
#####

get_summary(trial) = (;
    # Using some BenchmarkTools internals :/
    mem = BenchmarkTools.prettymemory(trial.memory),
    nalloc = trial.allocs,
    t_min = BenchmarkTools.prettytime(minimum(trial.times)),
    t_max = BenchmarkTools.prettytime(maximum(trial.times)),
    t_mean = BenchmarkTools.prettytime(StatsBase.mean(trial.times)),
    t_med = BenchmarkTools.prettytime(StatsBase.median(trial.times)),
    n_samples = length(trial),
)

function tabulate_summary(summary)
    summary_keys = collect(keys(summary))
    mem = map(k -> summary[k].mem, summary_keys)
    nalloc = map(k -> summary[k].nalloc, summary_keys)
    t_mean = map(k -> summary[k].t_mean, summary_keys)
    t_min = map(k -> summary[k].t_min, summary_keys)
    t_max = map(k -> summary[k].t_max, summary_keys)
    t_med = map(k -> summary[k].t_med, summary_keys)
    n_samples = map(k -> summary[k].n_samples, summary_keys)

    table_data = hcat(
        string.(collect(keys(summary))),
        mem,
        nalloc,
        t_min,
        t_max,
        t_mean,
        t_med,
        n_samples,
    )

    # Create a single header row merging the two-row header logic
    header_row = [
    "IndepVars (+conditions)" "Memory estimate" "allocs estimate" "Time min" "Time max" "Time mean" "Time median" "N-samples"
]

    final_table = vcat(header_row, table_data)

    println(
        "Summary Table (Columns: IndepVars, Memory, Allocations, Time(min, max, mean, median), N-samples):",
    )
    display(final_table)
end

#####
##### Constructor-specific configurations
#####

get_kwargs(x, ::Type{TD.ρeq}) = (; ρ = x.ρ, e_int = x.e_int_ρ, q_tot = x.q_tot)
get_kwargs(x, ::Type{TD.peq}) = (; p = x.p, e_int = x.e_int_p, q_tot = x.q_tot)
get_kwargs(x, ::Type{TD.phq}) = (; p = x.p, h = x.h_p, q_tot = x.q_tot)
get_kwargs(x, ::Type{TD.pρq}) = (; p = x.p, ρ = x.ρ, q_tot = x.q_tot)
get_kwargs(x, ::Type{TD.ρθ_li_q}) = (; ρ = x.ρ, θ_li = x.θ_ρ, q_tot = x.q_tot)
get_kwargs(x, ::Type{TD.pθ_li_q}) = (; p = x.p, θ_li = x.θ_p, q_tot = x.q_tot)

conditions(::Type) = (:dry, :saturated)

solver_for(::Type{TD.ρeq}) = RS.NewtonsMethod
solver_for(::Type) = RS.SecantMethod
