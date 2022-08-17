import Thermodynamics
import StatsBase
import PrettyTables
import OrderedCollections
const TD = Thermodynamics
const TP = TD.Parameters
using JET
using Test

import UnPack
import BenchmarkTools

import CLIMAParameters
const CP = CLIMAParameters
function get_parameter_set(::Type{FT}) where {FT}
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    aliases = string.(fieldnames(TP.ThermodynamicsParameters))
    param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    param_set = TP.ThermodynamicsParameters{FT}(; param_pairs...)
    # logfilepath = joinpath(@__DIR__, "logfilepath_$FT.toml")
    # CP.log_parameter_information(toml_dict, logfilepath)
    return param_set
end

#####
##### Finding indexes in profiles satisfying certain conditions
#####

function find_freezing_index(profiles, param_set)
    i = findfirst(T -> T === TP.T_freeze(param_set), profiles.T)
    isnothing(i) && error("Freezing index not found")
    return i
end

function find_dry(profiles)
    i = findfirst(
        q_tot -> q_tot === Float64(0) || q_tot === Float32(0),
        profiles.q_tot,
    )
    isnothing(i) && error("Dry index not found")
    return i
end

function find_moist_index(profiles)
    i = findfirst(q_tot -> q_tot > 0.01, profiles.q_tot)
    isnothing(i) && error("Moist index not found")
    return i
end

function find_sat_adjust_index(profiles, param_set, f::TDC) where {TDC}
    kwargs = get_kwargs(profiles, f)
    ts = f.(param_set, values(kwargs)...) # TDC is the thermo constructor
    T_sa = TD.air_temperature.(param_set, ts)
    Z = zip(T_sa, profiles.T)
    i = findfirst(x -> abs(x[1] - x[2]) > 0.0001, collect(Z))
    isnothing(i) && error("Saturation adjustment index not found")
    # @info "T_error for i_sat_adjust = $(abs(T_sa[i] - T[i]))"
    return i
end

#=
    conditions_index

Find an index in the given `profiles` that satisfies
some condition. For example, find the index that we're
sure that saturation adjustment is performed for a given
(profiles.p[i], profiles.θ_liq_ice[i], profiles.q_tot[i]).
=#
function conditions_index(profiles, param_set, sym, constructor)
    i = if sym == :dry
        find_dry(profiles)
    elseif sym == :freezing
        find_freezing_index(profiles, param_set)
    elseif sym == :sat_adjust
        find_sat_adjust_index(profiles, param_set, constructor)
    elseif sym == :moist
        find_moist_index(profiles)
    else
        error("Bad sym given")
    end
    return i
end

function sample_args(profiles, param_set, sym, constructor)
    i = conditions_index(profiles, param_set, sym, constructor)
    kwargs = get_kwargs(profiles, constructor)
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

    header = (
        [
            "Constructor",
            "Memory",
            "allocs",
            "Time",
            "Time",
            "Time",
            "Time",
            "N-samples",
        ],
        [
            "(+conditions)",
            "estimate",
            "estimate",
            "min",
            "max",
            "mean",
            "median",
            "",
        ],
    )

    PrettyTables.pretty_table(
        table_data;
        header,
        crop = :none,
        alignment = vcat(:l, repeat([:r], length(header[1]) - 1)),
    )
end

#####
##### Constructor-specific configurations
#####

# unpack variables from profiles into NamedTuple:
up(profiles, syms) = (; zip(syms, getproperty.(Ref(profiles), syms))...)
get_kwargs(p, ::typeof(TD.PhaseEquil_ρeq)) = up(p, :(ρ, e_int, q_tot).args)
get_kwargs(p, ::typeof(TD.PhaseEquil_ρTq)) = up(p, :(p, T, q_tot).args)
get_kwargs(p, ::typeof(TD.PhaseEquil_pTq)) = up(p, :(p, T, q_tot).args)
get_kwargs(p, ::typeof(TD.PhaseEquil_peq)) = up(p, :(p, e_int, q_tot).args)
get_kwargs(p, ::typeof(TD.PhaseEquil_phq)) = up(p, :(p, h, q_tot).args)
get_kwargs(p, ::typeof(TD.PhaseEquil_ρθq)) = up(p, :(ρ, θ_liq_ice, q_tot).args)
get_kwargs(p, ::typeof(TD.PhaseEquil_pθq)) = up(p, :(p, θ_liq_ice, q_tot).args)
get_kwargs(p, ::typeof(TD.PhaseEquil_ρpq)) = up(p, :(ρ, p, q_tot).args)

# Conditions to perform microbenchmarks and JET tests:
# note: No freezing points exist in
#       TD.TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)!
#       so we're not testing performance of these branches.

conditions(::typeof(TD.PhaseEquil_ρeq)) = (:dry, :sat_adjust)
conditions(::typeof(TD.PhaseEquil_ρTq)) = (:dry, :moist)
conditions(::typeof(TD.PhaseEquil_pTq)) = (:dry, :moist) # no sat adjust exists!
conditions(::typeof(TD.PhaseEquil_peq)) = (:dry, :sat_adjust)
conditions(::typeof(TD.PhaseEquil_phq)) = (:dry, :sat_adjust)
conditions(::typeof(TD.PhaseEquil_ρθq)) = (:dry, :sat_adjust)
conditions(::typeof(TD.PhaseEquil_pθq)) = (:dry, :sat_adjust)
conditions(::typeof(TD.PhaseEquil_ρpq)) = (:dry, :sat_adjust)
