using Pkg.Artifacts

using ArtifactWrappers

# Get dycoms dataset folder:
dycoms_dataset = ArtifactWrapper(
    @__DIR__,
    isempty(get(ENV, "CI", "")),
    "dycoms",
    ArtifactFile[ArtifactFile(
        url = "https://caltech.box.com/shared/static/bxau6i46y6ikxn2sy9krgz0sw5vuptfo.nc",
        filename = "test_data_PhaseEquil.nc",
    ),],
)
dycoms_dataset_path = get_data_folder(dycoms_dataset)


@testset "Data tests" begin
    FT = Float64

    data = joinpath(dycoms_dataset_path, "test_data_PhaseEquil.nc")
    ds_PhaseEquil = Dataset(data, "r")
    e_int = Array{FT}(ds_PhaseEquil["e_int"][:])
    ρ = Array{FT}(ds_PhaseEquil["ρ"][:])
    q_tot = Array{FT}(ds_PhaseEquil["q_tot"][:])

    ts = PhaseEquil_ρeq.(Ref(param_set), ρ, e_int, q_tot, 4)
    # ts = PhaseEquil_ρeq.(Ref(param_set), ρ, e_int, q_tot, 3) # Fails
end
