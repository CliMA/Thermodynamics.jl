using Test
import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import Thermodynamics.TemperatureProfiles as TDTP
using ForwardDiff

import ClimaParams as CP

"""
    test_temperature_profile_consistency(profile, param_set, z, _p_ref, _grav, _R_d)

Test that a temperature profile is physically consistent by verifying:
1. Surface pressure matches the reference pressure
2. Virtual temperature can be reconstructed from pressure gradient using hydrostatic balance

# Arguments
- `profile`: Temperature profile function to test
- `param_set`: Thermodynamic parameter set
- `z`: Array of height levels
- `_p_ref`: Reference pressure (typically MSLP)
- `_grav`: Gravitational acceleration
- `_R_d`: Gas constant for dry air
"""
function test_temperature_profile_consistency(
    profile,
    param_set,
    z,
    _p_ref,
    _grav,
    _R_d,
)
    # Evaluate profile at all height levels
    args = profile.(Ref(param_set), z)
    T_virt = first.(args)
    p = last.(args)

    # Test 1: Surface pressure should equal reference pressure
    mask_z_0 = z .≈ 0
    @test all(p[mask_z_0] .≈ _p_ref)

    # Test 2: Virtual temperature can be reconstructed from pressure gradient
    # This verifies hydrostatic balance consistency
    function log_p_over_p_ref(_z)
        _, _p = profile(param_set, _z)
        return log(_p / _p_ref)
    end

    ∇log_p_over_p_ref = _z -> ForwardDiff.derivative(log_p_over_p_ref, _z)

    # Reconstruct virtual temperature using hydrostatic balance:
    # T_virt = -g / (R_d * d(log(p/p_ref))/dz)
    T_virt_rec = -_grav ./ (_R_d .* ∇log_p_over_p_ref.(z))
    @test all(T_virt_rec .≈ T_virt)
end

@testset "TemperatureProfiles - DecayingTemperatureProfile" begin
    for FT in [Float32, Float64]
        param_set = TD.Parameters.ThermodynamicsParameters(FT)
        _grav = FT(TP.grav(param_set))
        _R_d = FT(TP.R_d(param_set))
        _p_ref = FT(TP.MSLP(param_set))

        # Define height range for testing (0 to 25 km)
        z = collect(range(FT(0), stop = FT(25e3), length = 100))

        # Number of test cases for each parameter
        n = 7

        # Define realistic parameter ranges for atmospheric temperature profiles
        _T_virt_surf_range = collect(range(FT(260), stop = FT(320), length = n))  # Surface temperature (K)
        _T_min_ref_range = collect(range(FT(170), stop = FT(230), length = n))    # Minimum temperature (K)
        _H_t_range = collect(range(FT(6e3), stop = FT(13e3), length = n))         # Scale height (m)

        # Test different temperature profile types
        profiles = []
        for (_T_virt_surf, _T_min_ref, _H_t) in
            zip(_T_virt_surf_range, _T_min_ref_range, _H_t_range)

            # Create different types of temperature profiles
            profiles = [
                # Decaying temperature profile (most realistic for troposphere)
                TDTP.DecayingTemperatureProfile{FT}(
                    param_set,
                    _T_virt_surf,
                    _T_min_ref,
                    _H_t,
                ),
                # Dry adiabatic profile (isentropic)
                TDTP.DryAdiabaticProfile{FT}(
                    param_set,
                    _T_virt_surf,
                    _T_min_ref,
                ),
                # Isothermal profile (constant temperature)
                TDTP.IsothermalProfile(param_set, _T_virt_surf),
                # Generic isothermal profile
                TDTP.IsothermalProfile(param_set, FT),
            ]

            # Test each profile for physical consistency
            for profile in profiles
                test_temperature_profile_consistency(
                    profile,
                    param_set,
                    z,
                    _p_ref,
                    _grav,
                    _R_d,
                )
            end
        end
    end
end

# Clean up any temporary files created during testing
rm(joinpath(@__DIR__, "logfilepath_Float32.toml"); force = true)
rm(joinpath(@__DIR__, "logfilepath_Float64.toml"); force = true)
