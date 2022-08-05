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
    header = (["Function", keys(flops[first(fnames)])...],)
    flops_neg = map(x -> flops[x]["neg"], fnames)
    flops_add = map(x -> flops[x]["add"], fnames)
    flops_sub = map(x -> flops[x]["sub"], fnames)
    flops_muladd = map(x -> flops[x]["muladd"], fnames)
    flops_fma = map(x -> flops[x]["fma"], fnames)
    flops_mul = map(x -> flops[x]["mul"], fnames)
    flops_div = map(x -> flops[x]["div"], fnames)
    flops_rem = map(x -> flops[x]["rem"], fnames)
    flops_sqrt = map(x -> flops[x]["sqrt"], fnames)
    flops_abs = map(x -> flops[x]["abs"], fnames)
    table_data = hcat(
        fnames,
        flops_neg,
        flops_add,
        flops_sub,
        flops_muladd,
        flops_fma,
        flops_mul,
        flops_div,
        flops_rem,
        flops_sqrt,
        flops_abs,
    )

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
