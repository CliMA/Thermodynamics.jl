#=
using Revise; "CuArray" in ARGS || push!(ARGS, "CuArray"); include("perf/kernel_bm.jl")
using Revise; include("perf/kernel_bm.jl")
=#
using Test

using KernelAbstractions
import BenchmarkTools as BMT
import KernelAbstractions as KA
using Random
using LinearAlgebra
import RootSolvers as RS

import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import ClimaParams as CP

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

@kernel function test_thermo_ka_kernel!(
    param_set,
    dst::AbstractArray{FT},
    e_int,
    ρ,
    q_tot,
) where {FT}
    i = @index(Global)
    @inbounds begin

        param_set = parameter_set(FT)
        ts = TD.PhaseEquil_ρeq(param_set, FT(ρ[i]), FT(e_int[i]), FT(q_tot[i]))

        s = FT(0)
        s += TD.air_pressure(param_set, ts)
        s += TD.air_temperature(param_set, ts)
        s += TD.air_density(param_set, ts)
        s += TD.specific_volume(param_set, ts)
        s += TD.soundspeed_air(param_set, ts)
        s += TD.total_specific_humidity(param_set, ts)
        s += TD.liquid_specific_humidity(param_set, ts)
        s += TD.ice_specific_humidity(param_set, ts)
        s += TD.vapor_specific_humidity(param_set, ts)
        s += TD.total_energy(param_set, ts, FT(0), FT(0))
        s += TD.internal_energy(param_set, ts)
        s += TD.internal_energy_sat(param_set, ts)
        s += TD.internal_energy_dry(param_set, ts)
        s += TD.internal_energy_vapor(param_set, ts)
        s += TD.internal_energy_liquid(param_set, ts)
        s += TD.internal_energy_ice(param_set, ts)
        s += TD.cp_m(param_set, ts)
        s += TD.cv_m(param_set, ts)
        s += TD.gas_constant_air(param_set, ts)
        # s += TD.gas_constants(param_set, ts)
        s += TD.latent_heat_vapor(param_set, ts)
        s += TD.latent_heat_sublim(param_set, ts)
        s += TD.latent_heat_fusion(param_set, ts)
        s += TD.latent_heat_liq_ice(param_set, ts)
        s += TD.saturation_vapor_pressure(param_set, ts, TD.Liquid())
        # s += = TD.q_vap_saturation_generic(param_set, ts)
        s += TD.q_vap_saturation(param_set, ts)
        s += TD.q_vap_saturation_liquid(param_set, ts)
        s += TD.q_vap_saturation_ice(param_set, ts)
        s += TD.saturation_excess(param_set, ts)
        s += TD.supersaturation(param_set, ts, TD.Liquid())
        s += TD.liquid_fraction(param_set, ts)
        s += TD.PhasePartition_equil(param_set, ts).tot
        s += TD.dry_pottemp(param_set, ts)
        s += TD.virtual_pottemp(param_set, ts)
        s += TD.virtual_dry_static_energy(param_set, ts, FT(0))
        s += TD.exner(param_set, ts)
        # s += TD.shum_to_mixing_ratio(param_set, ts)
        s += TD.mixing_ratios(param_set, ts).tot
        s += TD.vol_vapor_mixing_ratio(param_set, ts)
        s += TD.liquid_ice_pottemp(param_set, ts)
        s += TD.liquid_ice_pottemp_sat(param_set, ts)
        s += TD.relative_humidity(param_set, ts)
        s += TD.virtual_temperature(param_set, ts)
        s += TD.condensate_shum(param_set, ts)
        s += TD.has_condensate(param_set, ts)
        s += TD.specific_enthalpy(param_set, ts)
        s += TD.total_specific_enthalpy(param_set, ts, FT(0))
        s += TD.moist_static_energy(param_set, ts, FT(0))
        s += TD.specific_entropy(param_set, ts)
        s += TD.saturated(param_set, ts)
        dst[i] = s

    end
end

function test_thermo_bc_kernel!(x, param_set, profiles)
    (; e_int, ρ, q_tot) = profiles
    @. x = test_thermo_bc_kernel(param_set, e_int, ρ, q_tot)
    return nothing
end

function test_thermo_bc_kernel(
    param_set,
    e_int::FT,
    ρ::FT,
    q_tot::FT,
) where {FT}
    param_set = parameter_set(FT)
    ts = TD.PhaseEquil_ρeq(param_set, FT(ρ), FT(e_int), FT(q_tot))

    s = FT(0)
    s += TD.air_pressure(param_set, ts)
    s += TD.air_temperature(param_set, ts)
    s += TD.air_density(param_set, ts)
    s += TD.specific_volume(param_set, ts)
    s += TD.soundspeed_air(param_set, ts)
    s += TD.total_specific_humidity(param_set, ts)
    s += TD.liquid_specific_humidity(param_set, ts)
    s += TD.ice_specific_humidity(param_set, ts)
    s += TD.vapor_specific_humidity(param_set, ts)
    s += TD.total_energy(param_set, ts, FT(0), FT(0))
    s += TD.internal_energy(param_set, ts)
    s += TD.internal_energy_sat(param_set, ts)
    s += TD.internal_energy_dry(param_set, ts)
    s += TD.internal_energy_vapor(param_set, ts)
    s += TD.internal_energy_liquid(param_set, ts)
    s += TD.internal_energy_ice(param_set, ts)
    s += TD.cp_m(param_set, ts)
    s += TD.cv_m(param_set, ts)
    s += TD.gas_constant_air(param_set, ts)
    # s += TD.gas_constants(param_set, ts)
    s += TD.latent_heat_vapor(param_set, ts)
    s += TD.latent_heat_sublim(param_set, ts)
    s += TD.latent_heat_fusion(param_set, ts)
    s += TD.latent_heat_liq_ice(param_set, ts)
    s += TD.saturation_vapor_pressure(param_set, ts, TD.Liquid())
    # s += = TD.q_vap_saturation_generic(param_set, ts)
    s += TD.q_vap_saturation(param_set, ts)
    s += TD.q_vap_saturation_liquid(param_set, ts)
    s += TD.q_vap_saturation_ice(param_set, ts)
    s += TD.saturation_excess(param_set, ts)
    s += TD.supersaturation(param_set, ts, TD.Liquid())
    s += TD.liquid_fraction(param_set, ts)
    s += TD.PhasePartition_equil(param_set, ts).tot
    s += TD.dry_pottemp(param_set, ts)
    s += TD.virtual_pottemp(param_set, ts)
    s += TD.virtual_dry_static_energy(param_set, ts, FT(0))
    s += TD.exner(param_set, ts)
    # s += TD.shum_to_mixing_ratio(param_set, ts)
    s += TD.mixing_ratios(param_set, ts).tot
    s += TD.vol_vapor_mixing_ratio(param_set, ts)
    s += TD.liquid_ice_pottemp(param_set, ts)
    s += TD.liquid_ice_pottemp_sat(param_set, ts)
    s += TD.relative_humidity(param_set, ts)
    s += TD.virtual_temperature(param_set, ts)
    s += TD.condensate_shum(param_set, ts)
    s += TD.has_condensate(param_set, ts)
    s += TD.specific_enthalpy(param_set, ts)
    s += TD.total_specific_enthalpy(param_set, ts, FT(0))
    s += TD.moist_static_energy(param_set, ts, FT(0))
    s += TD.specific_entropy(param_set, ts)
    s += TD.saturated(param_set, ts)
    return s
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
    (; e_int, ρ, q_tot) = profiles
    n_profiles = length(profiles.z)
    ndrange = (n_profiles,)
    backend = KA.get_backend(d_dst)
    kernel! = test_thermo_ka_kernel!(backend)
    kernel!(param_set, d_dst, e_int, ρ, q_tot; ndrange = ndrange)
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
    d_dst = ArrayType(zeros(FT, n_profiles))
    fill!(d_dst, 0)

    test_thermo!(param_set, d_dst, profiles) # compile first
    trial = if d_dst isa Array
        BMT.@benchmark test_thermo!($param_set, $d_dst, $profiles)
    else
        BMT.@benchmark CUDA.@sync test_thermo!($param_set, $d_dst, $profiles)
    end
    show(stdout, MIME("text/plain"), trial)
    println("")

    @test !any(z -> iszero(z), Array(d_dst))

    x = ArrayType(zeros(FT, n_profiles))
    test_thermo_bc_kernel!(x, param_set, profiles) # compile first
    trial = if x isa Array
        BMT.@benchmark test_thermo_bc_kernel!($x, $param_set, $profiles)
    else
        BMT.@benchmark CUDA.@sync test_thermo_bc_kernel!($x, $param_set, $profiles)
    end
    show(stdout, MIME("text/plain"), trial)
    println("")

    @test !any(z -> iszero(z), Array(x))
end

rm(joinpath(@__DIR__, "logfilepath_Float32.toml"); force = true)
rm(joinpath(@__DIR__, "logfilepath_Float64.toml"); force = true)
