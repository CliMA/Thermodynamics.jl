#=
using Revise; "CuArray" in ARGS || push!(ARGS, "CuArray"); include("perf/kernel_bm.jl")
using Revise; include("perf/kernel_bm.jl")
=#
using Test

using KernelAbstractions
using BenchmarkTools
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

        dst[1, i] = TD.air_pressure(param_set, ts)
        dst[2, i] = TD.air_temperature(param_set, ts)
        dst[3, i] = TD.air_density(param_set, ts)
        dst[4, i] = TD.specific_volume(param_set, ts)
        dst[5, i] = TD.soundspeed_air(param_set, ts)
        dst[6, i] = TD.total_specific_humidity(param_set, ts)
        dst[7, i] = TD.liquid_specific_humidity(param_set, ts)
        dst[8, i] = TD.ice_specific_humidity(param_set, ts)
        dst[9, i] = TD.vapor_specific_humidity(param_set, ts)
        dst[10, i] = TD.total_energy(param_set, ts, FT(0), FT(0))
        dst[11, i] = TD.internal_energy(param_set, ts)
        dst[12, i] = TD.internal_energy_sat(param_set, ts)
        dst[13, i] = TD.internal_energy_dry(param_set, ts)
        dst[14, i] = TD.internal_energy_vapor(param_set, ts)
        dst[15, i] = TD.internal_energy_liquid(param_set, ts)
        dst[16, i] = TD.internal_energy_ice(param_set, ts)
        dst[17, i] = TD.cp_m(param_set, ts)
        dst[18, i] = TD.cv_m(param_set, ts)
        dst[19, i] = TD.gas_constant_air(param_set, ts)
        # dst[20, i] = TD.gas_constants(param_set, ts)
        dst[21, i] = TD.latent_heat_vapor(param_set, ts)
        dst[22, i] = TD.latent_heat_sublim(param_set, ts)
        dst[23, i] = TD.latent_heat_fusion(param_set, ts)
        dst[24, i] = TD.latent_heat_liq_ice(param_set, ts)
        dst[25, i] = TD.saturation_vapor_pressure(param_set, ts, TD.Liquid())
        # dst[26, i] = TD.q_vap_saturation_generic(param_set, ts)
        dst[27, i] = TD.q_vap_saturation(param_set, ts)
        dst[28, i] = TD.q_vap_saturation_liquid(param_set, ts)
        dst[29, i] = TD.q_vap_saturation_ice(param_set, ts)
        dst[30, i] = TD.saturation_excess(param_set, ts)
        dst[31, i] = TD.supersaturation(param_set, ts, TD.Liquid())
        dst[32, i] = TD.liquid_fraction(param_set, ts)
        dst[33, i] = TD.PhasePartition_equil(param_set, ts).tot
        dst[34, i] = TD.dry_pottemp(param_set, ts)
        dst[35, i] = TD.virtual_pottemp(param_set, ts)
        dst[36, i] = TD.virtual_dry_static_energy(param_set, ts, FT(0))
        dst[37, i] = TD.exner(param_set, ts)
        # dst[38, i] = TD.shum_to_mixing_ratio(param_set, ts)
        dst[39, i] = TD.mixing_ratios(param_set, ts).tot
        dst[40, i] = TD.vol_vapor_mixing_ratio(param_set, ts)
        dst[41, i] = TD.liquid_ice_pottemp(param_set, ts)
        dst[42, i] = TD.liquid_ice_pottemp_sat(param_set, ts)
        dst[43, i] = TD.relative_humidity(param_set, ts)
        dst[44, i] = TD.virtual_temperature(param_set, ts)
        dst[45, i] = TD.condensate(param_set, ts)
        dst[46, i] = TD.has_condensate(param_set, ts)
        dst[47, i] = TD.specific_enthalpy(param_set, ts)
        dst[48, i] = TD.total_specific_enthalpy(param_set, ts, FT(0))
        dst[49, i] = TD.moist_static_energy(param_set, ts, FT(0))
        dst[50, i] = TD.specific_entropy(param_set, ts)
        dst[51, i] = TD.saturated(param_set, ts)

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

function test_thermo!(param_set, d_dst, profiles)
    (; e_int, ρ, p, q_tot) = profiles
    n_profiles = length(profiles.z)
    ndrange = (n_profiles,)
    backend = KA.get_backend(d_dst)
    kernel! = test_thermo_kernel!(backend)
    kernel!(param_set, d_dst, e_int, ρ, p, q_tot; ndrange = ndrange)
    KA.synchronize(backend)
    return nothing
end

@testset "Thermodynamics - kernels" begin
    FT = Float32
    param_set = parameter_set(FT)
    profiles = TD.TestedProfiles.PhaseEquilProfiles(param_set, Array)
    slice = Colon()
    profiles = convert_profile_set(profiles, ArrayType, slice)

    n_profiles = length(profiles.z)
    n_vars = length(propertynames(profiles))
    d_dst = ArrayType(Array{FT}(undef, 51, n_profiles))
    fill!(d_dst, 0)

    test_thermo!(param_set, d_dst, profiles) # compile first
    trial =
        BenchmarkTools.@benchmark test_thermo!($param_set, $d_dst, $profiles)
    show(stdout, MIME("text/plain"), trial)

    (; e_int, ρ, p, q_tot) = profiles
    # Test
    ts_cpu =
        TD.PhaseEquil_ρeq.(
            param_set,
            Array{FT}(ρ),
            Array{FT}(e_int),
            Array{FT}(q_tot),
        )
    @test all(Array(d_dst)[1, :] .≈ TD.air_pressure.(param_set, ts_cpu))
    @test all(Array(d_dst)[2, :] .≈ TD.air_temperature.(param_set, ts_cpu))

end

rm(joinpath(@__DIR__, "logfilepath_Float32.toml"); force = true)
rm(joinpath(@__DIR__, "logfilepath_Float64.toml"); force = true)
