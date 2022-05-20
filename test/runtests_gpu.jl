if !haskey(ENV, "BUILDKITE")
    import Pkg
    Pkg.develop(Pkg.PackageSpec(; path = dirname(@__DIR__)))
end
using Test

using KernelAbstractions
const KA = KernelAbstractions

import CUDAKernels
const CK = CUDAKernels

import Thermodynamics
const TD = Thermodynamics

import UnPack

using Random
using LinearAlgebra

import RootSolvers
const RS = RootSolvers

import CLIMAParameters
const CP = CLIMAParameters

const ICP = TD.InternalClimaParams

function get_parameter_set(::Type{FT}) where {FT}
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    aliases = string.(fieldnames(ICP.ThermodynamicsParameters))
    param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    param_set = ICP.ThermodynamicsParameters{FT}(; param_pairs...)
    logfilepath = joinpath(@__DIR__, "logfilepath_$FT.toml")
    CP.log_parameter_information(toml_dict, logfilepath)
    return param_set
end

const param_set_Float64 = get_parameter_set(Float64)
const param_set_Float32 = get_parameter_set(Float32)
parameter_set(::Type{Float64}) = param_set_Float64
parameter_set(::Type{Float32}) = param_set_Float32


if get(ARGS, 1, "Array") == "CuArray"
    import CUDA
    ArrayType = CUDA.CuArray
    CUDA.allowscalar(false)
    device(::Type{T}) where {T <: CUDA.CuArray} = CK.CUDADevice()
else
    ArrayType = Array
    device(::Type{T}) where {T <: Array} = CK.CPU()
end

@show ArrayType

@kernel function test_thermo_kernel!(
    param_set,
    dst::AbstractArray{FT},
    e_int,
    ρ,
    p,
    q_tot,
) where {FT}
    i = @index(Group, Linear)
    @inbounds begin

        ts = TD.PhaseEquil_ρeq(param_set, FT(ρ[i]), FT(e_int[i]), FT(q_tot[i]))
        dst[1, i] = TD.air_temperature(param_set, ts)

        ts_ρpq = TD.PhaseEquil_ρpq(
            param_set,
            FT(ρ[i]),
            FT(p[i]),
            FT(q_tot[i]),
            true,
            100,
            RS.RegulaFalsiMethod,
        )
        dst[2, i] = TD.air_temperature(param_set, ts_ρpq)
    end
end

# Since we use `rand` to generate the ProfileSet,
# just initialize on the CPU, and provide convert
# function to move arrays to the GPU.
convert_profile_set(ps::TD.TestedProfiles.ProfileSet, ArrayType, slice) =
    TD.TestedProfiles.ProfileSet(
        ArrayType(ps.z[slice]),
        ArrayType(ps.T[slice]),
        ArrayType(ps.p[slice]),
        ArrayType(ps.RS[slice]),
        ArrayType(ps.e_int[slice]),
        ArrayType(ps.h[slice]),
        ArrayType(ps.ρ[slice]),
        ArrayType(ps.θ_liq_ice[slice]),
        ArrayType(ps.q_tot[slice]),
        ArrayType(ps.q_liq[slice]),
        ArrayType(ps.q_ice[slice]),
        TD.PhasePartition.(ps.q_tot[slice], ps.q_liq[slice], ps.q_ice[slice]),
        ArrayType(ps.RH[slice]),
        ArrayType(ps.e_pot[slice]),
        ArrayType(ps.u[slice]),
        ArrayType(ps.v[slice]),
        ArrayType(ps.w[slice]),
        ArrayType(ps.e_kin[slice]),
        ps.phase_type,
    )

@testset "Thermodynamics - kernels" begin
    FT = Float32
    param_set = parameter_set(FT)
    dev = device(ArrayType)
    profiles = TD.TestedProfiles.PhaseEquilProfiles(param_set, Array)
    slice = Colon()
    profiles = convert_profile_set(profiles, ArrayType, slice)

    n_profiles = length(profiles.z)
    n_vars = length(propertynames(profiles))
    d_dst = ArrayType(Array{FT}(undef, 2, n_profiles))
    fill!(d_dst, 0)

    UnPack.@unpack e_int, ρ, p, q_tot = profiles

    work_groups = (1,)
    ndrange = (n_profiles,)
    kernel! = test_thermo_kernel!(dev, work_groups)
    event = kernel!(param_set, d_dst, e_int, ρ, p, q_tot, ndrange = ndrange)
    wait(dev, event)

    ts_correct =
        TD.PhaseEquil_ρeq.(param_set, Array(ρ), Array(e_int), Array(q_tot))
    @test all(Array(d_dst)[1, :] .≈ TD.air_temperature.(param_set, ts_correct))

    ts_correct =
        TD.PhaseEquil_ρpq.(
            param_set,
            Array(ρ),
            Array(p),
            Array(q_tot),
            true,
            100,
            RS.RegulaFalsiMethod,
        )
    @test all(Array(d_dst)[2, :] .≈ TD.air_temperature.(param_set, ts_correct))

end

rm(joinpath(@__DIR__, "logfilepath_Float32.toml"); force = true)
rm(joinpath(@__DIR__, "logfilepath_Float64.toml"); force = true)
