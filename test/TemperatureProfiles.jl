using Test
import Thermodynamics
const TD = Thermodynamics
const TP = TD.Parameters
const TDTP = TD.TemperatureProfiles
using ForwardDiff

import CLIMAParameters
const CP = CLIMAParameters

function get_parameter_set(::Type{FT}) where {FT}
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    aliases = string.(fieldnames(TP.ThermodynamicsParameters))
    param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    param_set = TP.ThermodynamicsParameters{FT}(; param_pairs...)
    logfilepath = joinpath(@__DIR__, "logfilepath_$FT.toml")
    CP.log_parameter_information(toml_dict, logfilepath)
    return param_set
end

const param_set_Float64 = get_parameter_set(Float64)
const param_set_Float32 = get_parameter_set(Float32)
parameter_set(::Type{Float64}) = param_set_Float64
parameter_set(::Type{Float32}) = param_set_Float32

@testset "TemperatureProfiles - DecayingTemperatureProfile" begin
    for FT in [Float32, Float64]
        param_set = parameter_set(FT)
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
