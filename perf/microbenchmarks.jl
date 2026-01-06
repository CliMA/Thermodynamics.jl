include("common_micro_bm.jl")

function benchmark_thermo_states(::Type{FT}) where {FT}
    summary = OrderedCollections.OrderedDict()
    ArrayType = Array{FT}
    param_set = TP.ThermodynamicsParameters(FT)
    inputs = functional_inputs(param_set, FT)

    for F in (TD.ρeq, TD.peq, TD.phq, TD.pρq, TD.ρθ_liq_ice_q, TD.pθ_liq_ice_q)
        for cond in conditions(F)
            args = sample_args(inputs, param_set, cond, F)
            solver = solver_for(F)
            trial = BenchmarkTools.@benchmark TD.saturation_adjustment(
                $solver,
                $param_set,
                $(F)(),
                $args...,
                40,
                FT(1e-10),
            )
            summary[Symbol(nameof(F), :_, cond)] = get_summary(trial)
        end
    end
    @info "Microbenchmarks for $FT"
    tabulate_summary(summary)
end

benchmark_thermo_states(Float64)
benchmark_thermo_states(Float32)
