using Test
using KernelAbstractions
using CUDAKernels
using Thermodynamics
using Thermodynamics.TemperatureProfiles
using Thermodynamics.TestedProfiles
using Thermodynamics.TestedProfiles: ProfileSet
using UnPack
using Random
using RootSolvers
const TD = Thermodynamics

using LinearAlgebra
using CLIMAParameters
using CLIMAParameters.Planet

struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

if get(ARGS, 1, "Array") == "CuArray"
    using CUDA
    ArrayType = CUDA.CuArray
    CUDA.allowscalar(false)
    device(::Type{T}) where {T <: CuArray} = CUDADevice()
else
    ArrayType = Array
    device(::Type{T}) where {T <: Array} = CPU()
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

        ts = PhaseEquil_ρeq(param_set, FT(ρ[i]), FT(e_int[i]), FT(q_tot[i]))
        dst[1, i] = air_temperature(ts)

        ts_ρpq = PhaseEquil_ρpq(
            param_set,
            FT(ρ[i]),
            FT(p[i]),
            FT(q_tot[i]),
            true,
            100,
            RegulaFalsiMethod,
        )
        dst[2, i] = air_temperature(ts_ρpq)
    end
end

# Since we use `rand` to generate the ProfileSet,
# just initialize on the CPU, and provide convert
# function to move arrays to the GPU.
convert_profile_set(ps::ProfileSet, ArrayType, slice) = ProfileSet(
    ArrayType(ps.z[slice]),
    ArrayType(ps.T[slice]),
    ArrayType(ps.p[slice]),
    ArrayType(ps.RS[slice]),
    ArrayType(ps.e_int[slice]),
    ArrayType(ps.ρ[slice]),
    ArrayType(ps.θ_liq_ice[slice]),
    ArrayType(ps.q_tot[slice]),
    ArrayType(ps.q_liq[slice]),
    ArrayType(ps.q_ice[slice]),
    PhasePartition.(ps.q_tot[slice], ps.q_liq[slice], ps.q_ice[slice]),
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
    dev = device(ArrayType)
    profiles = TestedProfiles.PhaseEquilProfiles(param_set, Array)
    slice = Colon()
    profiles = convert_profile_set(profiles, ArrayType, slice)

    n_profiles = length(profiles.z)
    n_vars = length(propertynames(profiles))
    d_dst = ArrayType(Array{FT}(undef, 2, n_profiles))
    fill!(d_dst, 0)

    @unpack e_int, ρ, p, q_tot = profiles

    work_groups = (1,)
    ndrange = (n_profiles,)
    kernel! = test_thermo_kernel!(dev, work_groups)
    event = kernel!(param_set, d_dst, e_int, ρ, p, q_tot, ndrange = ndrange)
    wait(dev, event)

    ts_correct =
        PhaseEquil_ρeq.(param_set, Array(ρ), Array(e_int), Array(q_tot))
    @test all(Array(d_dst)[1, :] .≈ air_temperature.(ts_correct))

    ts_correct =
        PhaseEquil_ρpq.(
            param_set,
            Array(ρ),
            Array(p),
            Array(q_tot),
            true,
            100,
            RegulaFalsiMethod,
        )
    @test all(Array(d_dst)[2, :] .≈ air_temperature.(ts_correct))

end
