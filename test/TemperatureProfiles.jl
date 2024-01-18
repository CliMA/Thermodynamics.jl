using Test
import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import Thermodynamics.TemperatureProfiles as TDTP
using ForwardDiff

import CLIMAParameters as CP

@testset "TemperatureProfiles - DecayingTemperatureProfile" begin
    for FT in [Float32, Float64]
        param_set = TP.ThermodynamicsParameters(FT)
        _grav = FT(TP.grav(param_set))
        _R_d = FT(TP.R_d(param_set))
        _p_ref = FT(TP.MSLP(param_set))

        z = collect(range(FT(0), stop = FT(25e3), length = 100))
        n = 7
        _T_virt_surf_range = collect(range(FT(260), stop = FT(320), length = n))
        _T_min_ref_range = collect(range(FT(170), stop = FT(230), length = n))
        _H_t_range = collect(range(FT(6e3), stop = FT(13e3), length = n))

        profiles = []
        for (_T_virt_surf, _T_min_ref, _H_t) in
            zip(_T_virt_surf_range, _T_min_ref_range, _H_t_range)
            profiles = [
                TDTP.DecayingTemperatureProfile{FT}(
                    param_set,
                    _T_virt_surf,
                    _T_min_ref,
                    _H_t,
                ),
                TDTP.DryAdiabaticProfile{FT}(
                    param_set,
                    _T_virt_surf,
                    _T_min_ref,
                ),
                TDTP.IsothermalProfile(param_set, _T_virt_surf),
                TDTP.IsothermalProfile(param_set, FT),
            ]

            for profile in profiles
                args = profile.(Ref(param_set), z)
                T_virt = first.(args)
                p = last.(args)

                # Test that surface pressure is equal to the
                # specified boundary condition p_ref (by default set to MSLP)
                mask_z_0 = z .≈ 0
                @test all(p[mask_z_0] .≈ _p_ref)

                function log_p_over_p_ref(_z)
                    _, _p = profile(param_set, _z)
                    return log(_p / _p_ref)
                end
                ∇log_p_over_p_ref =
                    _z -> ForwardDiff.derivative(log_p_over_p_ref, _z)
                # Uses density computed from pressure derivative
                # and ideal gas law to reconstruct virtual temperature
                T_virt_rec = -_grav ./ (_R_d .* ∇log_p_over_p_ref.(z))
                @test all(T_virt_rec .≈ T_virt)
            end
        end

    end
end

rm(joinpath(@__DIR__, "logfilepath_Float32.toml"); force = true)
rm(joinpath(@__DIR__, "logfilepath_Float64.toml"); force = true)
