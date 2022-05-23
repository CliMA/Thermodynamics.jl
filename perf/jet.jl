include("common_micro_bm.jl")

# JET does not like KernelAbstractions.@print calls,
# so we disable them here.
TD.print_warning() = false

function jet_thermo_states()
    ArrayType = Array{Float64}
    profiles = TD.TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)

    @testset "JET tests" begin
        for C in (
            TD.PhaseEquil_ρeq,
            TD.PhaseEquil_peq,
            TD.PhaseEquil_pθq,
            TD.PhaseEquil_pTq,
        )
            for cond in conditions(C)
                args = sample_args(profiles, cond, C)
                JET.@test_opt C(param_set, args...)
            end
        end
    end
end

jet_thermo_states()
