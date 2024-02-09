include("common_micro_bm.jl")

function jet_thermo_states(::Type{FT}) where {FT}
    param_set = TP.ThermodynamicsParameters(FT)
    ArrayType = Array{FT}
    profiles = TD.TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)

    @testset "JET tests" begin
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
                JET.@test_opt C(TD.NullLogger(), param_set, args...)
            end
        end
    end
end

jet_thermo_states(Float32)
