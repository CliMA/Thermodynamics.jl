using Pkg.Artifacts

using ArtifactWrappers

# Get dycoms dataset folder:
dycoms_dataset = ArtifactWrapper(
    @__DIR__,
    "dycoms",
    ArtifactFile[ArtifactFile(
        url = "https://caltech.box.com/shared/static/bxau6i46y6ikxn2sy9krgz0sw5vuptfo.nc",
        filename = "test_data_PhaseEquil.nc",
    ),],
)
dycoms_dataset_path = get_data_folder(dycoms_dataset)


@testset "Data tests" begin
    FT = Float64
    param_set = TP.ThermodynamicsParameters(FT)
    data = joinpath(dycoms_dataset_path, "test_data_PhaseEquil.nc")
    ds_PhaseEquil = Dataset(data, "r")
    e_int = Array{FT}(ds_PhaseEquil["e_int"][:])
    ρ = Array{FT}(ds_PhaseEquil["ρ"][:])
    q_tot = Array{FT}(ds_PhaseEquil["q_tot"][:])

    ts = PhaseEquil_ρeq.(Ref(param_set), ρ, e_int, q_tot, 4)
    # ts = PhaseEquil_ρeq.(Ref(param_set), ρ, e_int, q_tot, 3) # Fails
end

@testset "pθq data-driven tests" begin
    param_set = TP.ThermodynamicsParameters(Float64)
    #! format: off
    pθq_broken = [
        # (; p = , θ_liq_ice = , q_tot = ),
        # (; p = , θ_liq_ice = , q_tot = ),
        # (; p = , θ_liq_ice = , q_tot = ),
        # (; p = , θ_liq_ice = , q_tot = ),
        # (; p = , θ_liq_ice = , q_tot = ),
    ]
    #! format: on
    # config = (50, 1e-3, RootSolvers.RegulaFalsiMethod)
    # config = (50, 1e-3, RootSolvers.NewtonMethodAD)
    config = ()
    if !isempty(pθq_broken)
        p_arr = getproperty.(pθq_broken, :p)
        θ_liq_ice_arr = getproperty.(pθq_broken, :θ_liq_ice)
        q_tot_arr = getproperty.(pθq_broken, :q_tot)
        for (p, θ_liq_ice, q_tot) in zip(p_arr, θ_liq_ice_arr, q_tot_arr)
            @test_broken begin
                ts = PhaseEquil_pθq(param_set, p, θ_liq_ice, q_tot, config...)
                air_pressure(param_set, ts) == p
            end
        end
    end

    #! format: off
    pθq = [
        (; p = 82307.30319719888, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 88357.42589002676, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 81090.35731696963, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 75461.41343839701, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 71158.23329080557, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 68032.73180723468, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 59906.16044902247, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 51740.77281055945, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 50056.7894012541, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 48806.70567178993, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 27564.73889538213, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
    ]
    #! format: on
    if !isempty(pθq)
        p_arr = getproperty.(pθq, :p)
        θ_liq_ice_arr = getproperty.(pθq, :θ_liq_ice)
        q_tot_arr = getproperty.(pθq, :q_tot)
        for (p, θ_liq_ice, q_tot) in zip(p_arr, θ_liq_ice_arr, q_tot_arr)
            @test begin
                ts = PhaseEquil_pθq(param_set, p, θ_liq_ice, q_tot, config...)
                air_pressure(param_set, ts) == p
            end
        end
    end
end
