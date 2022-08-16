include("common_micro_bm.jl")

function benchmark_thermo_states(::Type{FT}) where {FT}
    summary = OrderedCollections.OrderedDict()
    ArrayType = Array{FT}
    param_set = get_parameter_set(FT)
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
            args = sample_args(profiles, param_set, cond, C)
            trial = BenchmarkTools.@benchmark $C($param_set, $args...)
            summary[Symbol(nameof(C), :_, cond)] = get_summary(trial)
        end
    end
    @info "Microbenchmarks for $FT"
    tabulate_summary(summary)
end

benchmark_thermo_states(Float64)
benchmark_thermo_states(Float32)
