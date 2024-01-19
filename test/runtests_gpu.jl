#=
using Revise
push!(ARGS, "CuArray")
include("test/runtests_gpu.jl")
=#
using Test

using KernelAbstractions
import KernelAbstractions as KA
using Random
using LinearAlgebra
import RootSolvers as RS

import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import CLIMAParameters as CP

if get(ARGS, 1, "Array") == "CuArray"
    import CUDA
    ArrayType = CUDA.CuArray
    CUDA.allowscalar(false)
else
    ArrayType = Array
end

const param_set_Float64 = TP.ThermodynamicsParameters(Float64)
const param_set_Float32 = TP.ThermodynamicsParameters(Float32)
parameter_set(::Type{Float64}) = param_set_Float64
parameter_set(::Type{Float32}) = param_set_Float32

@show ArrayType

@kernel function test_thermo_kernel!(
    param_set,
    dst::AbstractArray{FT},
    e_int,
    ρ,
    p,
    q_tot,
) where {FT}
    i = @index(Global)
    @inbounds begin

        param_set = parameter_set(FT)
        ts = TD.PhaseEquil_ρeq(param_set, FT(ρ[i]), FT(e_int[i]), FT(q_tot[i]))
        dst[1, i] = TD.air_temperature(param_set, ts)

        ts_ρpq = TD.PhaseEquil_ρpq(
            param_set,
            FT(ρ[i]),
            FT(p[i]),
            FT(q_tot[i]),
            true,
            100,
            FT(sqrt(eps(FT))),
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
    profiles = TD.TestedProfiles.PhaseEquilProfiles(param_set, Array)
    slice = Colon()
    profiles = convert_profile_set(profiles, ArrayType, slice)

    n_profiles = length(profiles.z)
    n_vars = length(propertynames(profiles))
    d_dst = ArrayType(Array{FT}(undef, 2, n_profiles))
    fill!(d_dst, 0)

    (; e_int, ρ, p, q_tot) = profiles

    ndrange = (n_profiles,)
    backend = KA.get_backend(d_dst)
    kernel! = test_thermo_kernel!(backend)
    kernel!(param_set, d_dst, e_int, ρ, p, q_tot; ndrange = ndrange)
    KA.synchronize(backend)

    ts_cpu =
        TD.PhaseEquil_ρeq.(
            param_set,
            Array{FT}(ρ),
            Array{FT}(e_int),
            Array{FT}(q_tot),
        )
    @test all(Array(d_dst)[1, :] .≈ TD.air_temperature.(param_set, ts_cpu))

    ts_correct =
        TD.PhaseEquil_ρpq.(
            param_set,
            Array(ρ),
            Array(p),
            Array(q_tot),
            true,
            100,
            sqrt(eps()),
            RS.RegulaFalsiMethod,
        )
    @test all(Array(d_dst)[2, :] .≈ TD.air_temperature.(param_set, ts_correct))

end

rm(joinpath(@__DIR__, "logfilepath_Float32.toml"); force = true)
rm(joinpath(@__DIR__, "logfilepath_Float64.toml"); force = true)
