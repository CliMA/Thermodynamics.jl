module CreateParametersExt

import Thermodynamics.Parameters.ThermodynamicsParameters
import ClimaParams as CP

ThermodynamicsParameters(::Type{FT}) where {FT <: Real} =
    ThermodynamicsParameters(CP.create_toml_dict(FT))

function ThermodynamicsParameters(toml_dict::CP.ParamDict{FT}) where {FT}
    name_map = (;
        :temperature_min_reference => :T_min_ref,
        :entropy_water_vapor => :entropy_water_vapor,
        :entropy_dry_air => :entropy_dry_air,
        :potential_temperature_reference_pressure => :p_ref_theta,
        :entropy_reference_temperature => :entropy_reference_temperature,
        :temperature_saturation_adjustment_max => :T_max,
        :pow_icenuc => :pow_icenuc,
        :temperature_triple_point => :T_triple,
        :gas_constant_dry_air => :R_d,
        :gas_constant_vapor => :R_v,
        :isobaric_specific_heat_dry_air => :cp_d,
        :pressure_triple_point => :press_triple,
        :thermodynamics_temperature_reference => :T_0,
        :temperature_water_freeze => :T_freeze,
        :isobaric_specific_heat_ice => :cp_i,
        :latent_heat_sublimation_at_reference => :LH_s0,
        :isobaric_specific_heat_vapor => :cp_v,
        :mean_sea_level_pressure => :MSLP,
        :isobaric_specific_heat_liquid => :cp_l,
        :latent_heat_vaporization_at_reference => :LH_v0,
        :temperature_saturation_adjustment_min => :T_min,
        :temperature_saturation_adjustment_init_min => :T_init_min,
        :temperature_surface_reference => :T_surf_ref,
        :gravitational_acceleration => :grav,
        :temperature_homogenous_nucleation => :T_icenuc,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "Thermodynamics")
    return ThermodynamicsParameters{FT}(; parameters...)
end

end
