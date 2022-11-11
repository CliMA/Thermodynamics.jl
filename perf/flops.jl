include("common.jl")
import GFlops
import OrderedCollections
import PrettyTables

function total_flops(flops)
    n_flops = Dict()
    strip_nums(x) =
        strip("$x", ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'])
    for pn in strip_nums.(propertynames(flops))
        n_flops[pn] = 0
    end
    for pn in propertynames(flops)
        n_flops[strip_nums(pn)] += getproperty(flops, pn)
    end
    return n_flops
end

function summarize_flops(flops)
    fnames = collect(keys(flops))
    flop_types = keys(flops[first(fnames)])
    header = (["Function", flop_types...],)
    flop_vals = map(collect(flop_types)) do ft
        map(x -> flops[x][ft], fnames)
    end
    table_data = hcat(fnames, flop_vals...)

    PrettyTables.pretty_table(
        table_data;
        header,
        crop = :none,
        alignment = vcat(:l, repeat([:c], length(header[1]) - 1)),
    )
end

#! format: off
@testset "Thermodynamics - Flops" begin

    flops = OrderedCollections.OrderedDict()
    ts = TD.PhaseEquil_ρeq(param_set, ρ[1], e_int[1], q_tot[1])

    flops["PhasePartition"] = total_flops(GFlops.@count_ops TD.PhasePartition($param_set, $ts))
    flops["liquid_fraction"] = total_flops(GFlops.@count_ops TD.liquid_fraction($param_set, $ts))
    flops["q_vap_saturation"] = total_flops(GFlops.@count_ops TD.q_vap_saturation($param_set, $ts))

    summarize_flops(flops)
end
#! format: on
