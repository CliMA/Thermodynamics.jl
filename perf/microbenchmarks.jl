include("common_micro_bm.jl")

# unpack variables from profiles into NamedTuple:
up(profiles, syms) = (; zip(syms, getproperty.(Ref(profiles), syms))...)
get_kwargs(p, ::typeof(TD.PhaseEquil_ρeq)) = up(p, :(ρ, e_int, q_tot).args)
get_kwargs(p, ::typeof(TD.PhaseEquil_pθq)) = up(p, :(p, θ_liq_ice, q_tot).args)
get_kwargs(p, ::typeof(TD.PhaseEquil_peq)) = up(p, :(p, e_int, q_tot).args)
get_kwargs(p, ::typeof(TD.PhaseEquil_pTq)) = up(p, :(p, T, q_tot).args)

conditions(::typeof(TD.PhaseEquil_ρeq)) = (:dry, :sat_adjust)
conditions(::typeof(TD.PhaseEquil_peq)) = (:dry, :sat_adjust)
conditions(::typeof(TD.PhaseEquil_pθq)) = (:dry, :sat_adjust)
conditions(::typeof(TD.PhaseEquil_pTq)) = (:dry, :moist) # no sat adjust exists!

# No freezing points exist in
# TD.TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)!
# so we're not testing performance of these branches.
function benchmark_thermo_states()
    summary = OrderedCollections.OrderedDict()
    ArrayType = Array{Float64}
    profiles = TD.TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)

    for C in (
        TD.PhaseEquil_ρeq,
        TD.PhaseEquil_peq,
        TD.PhaseEquil_pθq,
        TD.PhaseEquil_pTq,
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
