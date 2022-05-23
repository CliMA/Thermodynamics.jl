include("common_micro_bm.jl")

function benchmark_thermo_states()
    summary = OrderedCollections.OrderedDict()
    ArrayType = Array{Float64}
    profiles = TD.TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)

    for C in (
        TD.PhaseEquil_ρeq,
        TD.PhaseEquil_ρTq,
        TD.PhaseEquil_pTq,
        TD.PhaseEquil_peq,
        TD.PhaseEquil_phq,
        TD.PhaseEquil_ρθq,
        TD.PhaseEquil_pθq,
        TD.PhaseEquil_ρpq,
    )
        for cond in conditions(C)
            args = sample_args(profiles, cond, C)
            trial = BenchmarkTools.@benchmark $C(param_set, $args...)
            summary[Symbol(nameof(C), :_, cond)] = get_summary(trial)
        end
    end

    tabulate_summary(summary)
end

benchmark_thermo_states()
