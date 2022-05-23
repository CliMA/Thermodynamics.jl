import Thermodynamics
import StatsBase
import PrettyTables
import OrderedCollections
const TD = Thermodynamics
const ICP = TD.InternalClimaParams

import UnPack
import BenchmarkTools

import CLIMAParameters
const CP = CLIMAParameters

struct EarthParameterSet <: CP.AbstractEarthParameterSet end
const param_set = EarthParameterSet()

function find_freezing_index(profiles)
    i = findfirst(T -> T === ICP.T_freeze(param_set), profiles.T)
    isnothing(i) && error("Freezing index not found")
    return i
end

function find_dry(profiles)
    i = findfirst(q_tot -> q_tot === 0.0, profiles.q_tot)
    isnothing(i) && error("Dry index not found")
    return i
end

function find_moist_index(profiles)
    i = findfirst(q_tot -> q_tot > 0.01, profiles.q_tot)
    isnothing(i) && error("Moist index not found")
    return i
end

function find_sat_adjust_index(profiles, f::TDC) where {TDC}
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
function conditions_index(profiles, sym, constructor)
    i = if sym == :dry
        find_dry(profiles)
    elseif sym == :freezing
        find_freezing_index(profiles)
    elseif sym == :sat_adjust
        find_sat_adjust_index(profiles, constructor)
    elseif sym == :moist
        find_moist_index(profiles)
    else
        error("Bad sym given")
    end
    return i
end

function sample_args(profiles, sym, constructor)
    i = conditions_index(profiles, sym, constructor)
    kwargs = get_kwargs(profiles, constructor)
    return getindex.(values(kwargs), i)
end

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

    println()
    PrettyTables.pretty_table(
        table_data;
        header,
        crop = :none,
        alignment = vcat(:l, repeat([:r], length(header[1]) - 1)),
    )
end
