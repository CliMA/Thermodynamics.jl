using Test

using KernelAbstractions
import KernelAbstractions as KA
using Random
using LinearAlgebra
import RootSolvers as RS

import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import ClimaParams as CP

include("TestedProfiles.jl")

# Determine array type based on command line argument
# This allows testing on both CPU (Array) and GPU (CuArray)
if get(ARGS, 1, "Array") == "CuArray"
    import CUDA
    ArrayType = CUDA.CuArray
    CUDA.allowscalar(false)  # Prevent scalar operations on GPU for performance
else
    ArrayType = Array
end

# Pre-allocate parameter sets for different floating-point types
const param_set_Float64 = TP.ThermodynamicsParameters(Float64)
const param_set_Float32 = TP.ThermodynamicsParameters(Float32)

"""
    parameter_set(::Type{FT})

Get the appropriate parameter set for a given floating-point type.

# Arguments
- `FT`: Floating-point type (Float32 or Float64)

# Returns
- Thermodynamic parameter set for the specified type
"""
parameter_set(::Type{Float64}) = param_set_Float64
parameter_set(::Type{Float32}) = param_set_Float32

@show ArrayType

"""
    test_thermo_kernel!(param_set, dst, e_int, ρ, p, q_tot)

GPU kernel for testing thermodynamic state constructors.
This kernel creates thermodynamic states from different input combinations
and computes air temperature to verify consistency.

# Arguments
- `param_set`: Thermodynamic parameter set
- `dst`: Destination array for storing computed temperatures
- `e_int`: Internal energy array
- `ρ`: Density array
- `p`: Pressure array
- `q_tot`: Total specific humidity array

# Details
The kernel tests two different thermodynamic state constructors:
1. PhaseEquil_ρeq: Uses density, internal energy, and humidity
2. PhaseEquil_ρpq: Uses density, pressure, and humidity with saturation adjustment
"""
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
        # Get parameter set for current floating-point type
        param_set = parameter_set(FT)

        # Test 1: PhaseEquil_ρeq constructor
        # Uses density, internal energy, and total humidity
        ts = TD.PhaseEquil_ρeq(param_set, FT(ρ[i]), FT(e_int[i]), FT(q_tot[i]))
        dst[1, i] = TD.air_temperature(param_set, ts)

        # Test 2: PhaseEquil_ρpq constructor with saturation adjustment
        # Uses density, pressure, and total humidity with iterative convergence
        ts_ρpq = TD.PhaseEquil_ρpq(
            param_set,
            FT(ρ[i]),
            FT(p[i]),
            FT(q_tot[i]),
            true,                    # Enable saturation adjustment
            100,                     # Maximum iterations
            FT(sqrt(eps(FT))),       # Tolerance
            RS.BrentsMethod,    # Root finding method
        )
        dst[2, i] = TD.air_temperature(param_set, ts_ρpq)
    end
end

"""
    convert_profile_set(ps, ArrayType, slice)

Convert a ProfileSet to use a different array type (e.g., for GPU computation).
This function moves arrays to the target device while preserving the structure.

# Arguments
- `ps`: ProfileSet to convert
- `ArrayType`: Target array type (Array or CuArray)
- `slice`: Index slice to apply to arrays

# Returns
- New ProfileSet with arrays converted to the target type
"""
convert_profile_set(ps::TestedProfiles.ProfileSet, ArrayType, slice) =
    TestedProfiles.ProfileSet(
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
        ArrayType(ps.RH[slice]),
        ArrayType(ps.e_pot[slice]),
        ArrayType(ps.u[slice]),
        ArrayType(ps.v[slice]),
        ArrayType(ps.w[slice]),
        ArrayType(ps.e_kin[slice]),
    )

@testset "Thermodynamics - kernels" begin
    FT = Float32  # Use Float32 for GPU compatibility
    param_set = parameter_set(FT)

    # Load test profiles and convert to target array type
    profiles = TestedProfiles.PhaseEquilProfiles(param_set, Array)
    slice = Colon()  # Use all profiles
    profiles = convert_profile_set(profiles, ArrayType, slice)

    # Set up kernel execution
    n_profiles = length(profiles.z)
    n_vars = length(propertynames(profiles))
    d_dst = ArrayType(Array{FT}(undef, 2, n_profiles))
    fill!(d_dst, 0)

    # Extract required variables from profiles
    (; e_int, ρ, p, q_tot) = profiles

    # Execute kernel
    ndrange = (n_profiles,)
    backend = KA.get_backend(d_dst)
    kernel! = test_thermo_kernel!(backend)
    kernel!(param_set, d_dst, e_int, ρ, p, q_tot; ndrange = ndrange)
    KA.synchronize(backend)

    # Test 1: Verify PhaseEquil_ρeq results
    # Compare GPU results with CPU reference implementation
    ts_cpu =
        TD.PhaseEquil_ρeq.(
            param_set,
            Array{FT}(ρ),
            Array{FT}(e_int),
            Array{FT}(q_tot),
        )
    @test all(Array(d_dst)[1, :] .≈ TD.air_temperature.(param_set, ts_cpu))

    # Test 2: Verify PhaseEquil_ρpq results
    # Compare GPU results with CPU reference implementation including saturation adjustment
    ts_correct =
        TD.PhaseEquil_ρpq.(
            param_set,
            Array(ρ),
            Array(p),
            Array(q_tot),
            true,                    # Enable saturation adjustment
            100,                     # Maximum iterations
            sqrt(eps()),            # Tolerance
            RS.BrentsMethod,   # Root finding method
        )
    @test all(Array(d_dst)[2, :] .≈ TD.air_temperature.(param_set, ts_correct))
end

# Clean up any temporary files created during testing
rm(joinpath(@__DIR__, "logfilepath_Float32.toml"); force = true)
rm(joinpath(@__DIR__, "logfilepath_Float64.toml"); force = true)
