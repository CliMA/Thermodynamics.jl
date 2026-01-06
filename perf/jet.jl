include("common_micro_bm.jl")

# JET does not like KernelAbstractions.@print calls,
# so we disable them here.
TD.print_warning() = false

function jet_thermo_states(::Type{FT}) where {FT}
    param_set = TP.ThermodynamicsParameters(FT)
    inputs = functional_inputs(param_set, FT)

    @testset "JET tests" begin
        for F in (TD.ρeq, TD.peq, TD.phq, TD.pρq, TD.ρθ_liq_ice_q, TD.pθ_liq_ice_q)
            for cond in conditions(F)
                args = sample_args(inputs, param_set, cond, F)
                solver = solver_for(F)
                JET.@test_opt TD.saturation_adjustment(
                    solver,
                    param_set,
                    F(),
                    args...,
                    40,
                    FT(1e-10),
                )
            end
        end
    end
end

jet_thermo_states(Float32)
