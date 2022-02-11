using Test
using Thermodynamics.TemperatureProfiles
using Thermodynamics

import CLIMAParameters
#using CLIMAParameters.Planet

using ForwardDiff

#struct EarthParameterSet <: AbstractEarthParameterSet end
#const param_set = EarthParameterSet()

param_set = Thermodynamics.ThermodynamicsParameters(full_parameter_set)

CLIMAParameters.write_log_file(full_parameter_set,"test_log_file.toml")

@testset "TemperatureProfiles - DecayingTemperatureProfile" begin
    for FT in [Float32, Float64]
        _grav = FT(param_set.grav)
        _R_d = FT(param_set.R_d)
        _MSLP = FT(param_set.MSLP)

        z = collect(range(FT(0), stop = FT(25e3), length = 100))
        n = 7
        _T_virt_surf_range = collect(range(FT(260), stop = FT(320), length = n))
        _T_min_ref_range = collect(range(FT(170), stop = FT(230), length = n))
        _H_t_range = collect(range(FT(6e3), stop = FT(13e3), length = n))

        profiles = []
        for (_T_virt_surf, _T_min_ref, _H_t) in
            zip(_T_virt_surf_range, _T_min_ref_range, _H_t_range)
            profiles = [
                DecayingTemperatureProfile{FT}(
                    param_set,
                    _T_virt_surf,
                    _T_min_ref,
                    _H_t,
                ),
                DryAdiabaticProfile{FT}(param_set, _T_virt_surf, _T_min_ref),
                IsothermalProfile(param_set, _T_virt_surf),
                IsothermalProfile(param_set, FT),
            ]

            for profile in profiles
                args = profile.(Ref(param_set), z)
                T_virt = first.(args)
                p = last.(args)

                # Test that surface pressure is equal to the
                # specified boundary condition (MSLP)
                mask_z_0 = z .≈ 0
                @test all(p[mask_z_0] .≈ _MSLP)

                function log_p_over_MSLP(_z)
                    _, _p = profile(param_set, _z)
                    return log(_p / _MSLP)
                end
                ∇log_p_over_MSLP =
                    _z -> ForwardDiff.derivative(log_p_over_MSLP, _z)
                # Uses density computed from pressure derivative
                # and ideal gas law to reconstruct virtual temperature
                T_virt_rec = -_grav ./ (_R_d .* ∇log_p_over_MSLP.(z))
                @test all(T_virt_rec .≈ T_virt)
            end
        end

    end
end
