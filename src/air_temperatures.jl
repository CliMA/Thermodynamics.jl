export air_temperature
export potential_temperature
export potential_temperature_given_pressure
export virtual_temperature
export virtual_pottemp
export liquid_ice_pottemp
export liquid_ice_pottemp_given_pressure

"""
    air_temperature(param_set, e_int, q_tot=0, q_liq=0, q_ice=0)

The air temperature.

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `e_int`: specific internal energy of moist air [J/kg]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `T`: air temperature [K]

If the specific humidities are not given, the result is for dry air.
This method inverts [`internal_energy`](@ref) by solving for `T` given `e_int`.
"""
@inline function air_temperature(
    param_set::APS,
    e_int::Real,     # Number type needed to disambiguate from deprecated methods that are still there
    q_tot::Real = 0,
    q_liq::Real = 0,
    q_ice::Real = 0,
)
    return air_temperature(param_set, ρe(), e_int, q_tot, q_liq, q_ice)
end

@inline function air_temperature(
    param_set::APS,
    ::Union{ρe, pe},
    e_int::Real,
    q_tot::Real = 0,
    q_liq::Real = 0,
    q_ice::Real = 0,
)
    T_0 = TP.T_0(param_set)
    R_d = TP.R_d(param_set)
    e_int_v0 = TP.e_int_v0(param_set)
    e_int_i0 = TP.e_int_i0(param_set)
    cvm = cv_m(param_set, q_tot, q_liq, q_ice)
    q_vap = vapor_specific_humidity(q_tot, q_liq, q_ice)
    return T_0 +
           (
        e_int - q_vap * e_int_v0 + q_ice * e_int_i0 + (1 - q_tot) * R_d * T_0
    ) / cvm
end

"""
    air_temperature(param_set, ::ph, h, q_tot=0, q_liq=0, q_ice=0)

The air temperature.

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `h`: specific enthalpy of moist air [J/kg]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `T`: air temperature [K]

If the specific humidities are not given, the result is for dry air.
This method inverts [`enthalpy`](@ref) by solving for `T` given `h`.
"""
@inline function air_temperature(
    param_set::APS,
    ::ph,
    h::Real,
    q_tot::Real = 0,
    q_liq::Real = 0,
    q_ice::Real = 0,
)
    cpm = cp_m(param_set, q_tot, q_liq, q_ice)
    T_0 = TP.T_0(param_set)
    LH_v0 = TP.LH_v0(param_set)
    LH_f0 = TP.LH_f0(param_set)
    q_vap = vapor_specific_humidity(q_tot, q_liq, q_ice)
    return T_0 + (h - q_vap * LH_v0 + q_ice * LH_f0) / cpm
end

"""
    air_temperature(param_set, ::pρ, p, ρ, q_tot=0, q_liq=0, q_ice=0)

The air temperature.

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `p`: air pressure [Pa]
 - `ρ`: air density [kg/m³]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `T`: air temperature [K]

If the specific humidities are not given, the result is for dry air.
This directly applies the ideal gas law: `T = p / (R_m ρ)`.
"""
@inline function air_temperature(
    param_set::APS,
    ::pρ,
    p::Real,
    ρ::Real,
    q_tot::Real = 0,
    q_liq::Real = 0,
    q_ice::Real = 0,
)
    R_m = gas_constant_air(param_set, q_tot, q_liq, q_ice)
    return p / (R_m * ρ)
end

"""
    air_temperature(
        param_set,
        ::pθ_li,
        p,
        θ_li,
        q_tot=0,
        q_liq=0,
        q_ice=0
    )

The air temperature obtained by inverting the liquid-ice potential temperature.

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `p`: pressure [Pa]
 - `θ_li`: liquid-ice potential temperature [K]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `T`: air temperature [K]

If the specific humidities are not given, `θ_li` is assumed to be the dry-air
potential temperature.
This inverts [`liquid_ice_pottemp_given_pressure`](@ref) by solving for `T`.
"""
@inline function air_temperature(
    param_set::APS,
    ::pθ_li,
    p::Real,
    θ_li::Real,
    q_tot::Real = 0,
    q_liq::Real = 0,
    q_ice::Real = 0,
)
    cpm = cp_m(param_set, q_tot, q_liq, q_ice)
    return θ_li * exner_given_pressure(param_set, p, q_tot, q_liq, q_ice) +
           humidity_weighted_latent_heat(param_set, q_liq, q_ice) / cpm
end


"""
    air_temperature(param_set, ::ρθ_li, ρ, θ_li, q_tot=0, q_liq=0, q_ice=0)

The air temperature obtained by inverting the liquid-ice potential temperature.

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `ρ`: (moist-)air density [kg/m³]
 - `θ_li`: liquid-ice potential temperature [K]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `T`: air temperature [K]

If the specific humidities are not given, the results are for dry air.

# Method
This function uses a second-order Taylor expansion to avoid an iterative inversion.
The unsaturated temperature `T_unsat` is first computed assuming the ideal gas law with
potential temperature, then latent heat corrections are applied.
"""
@inline function air_temperature(
    param_set::APS,
    ::ρθ_li,
    ρ::Real,
    θ_li::Real,
    q_tot::Real = 0,
    q_liq::Real = 0,
    q_ice::Real = 0,
)
    # Second-order Taylor expansion around unsaturated temperature
    p_ref = TP.p_ref_theta(param_set)
    cvm = cv_m(param_set, q_tot, q_liq, q_ice)
    cpm = cp_m(param_set, q_tot, q_liq, q_ice)
    R_m = gas_constant_air(param_set, q_tot, q_liq, q_ice)
    κ = 1 - cvm / cpm

    # Unsaturated temperature corresponding to (ρ, θ_li) in the dry/moist EOS sense
    T_unsat = (ρ * R_m * θ_li / p_ref)^(R_m / cvm) * θ_li

    # Latent-heat correction (humidity-weighted at reference temperature)
    L_q = humidity_weighted_latent_heat(param_set, q_liq, q_ice)
    ΔT₁ = L_q / cvm
    ΔT₂ = -κ / (2 * T_unsat) * (L_q / cvm)^2
    return T_unsat + ΔT₁ + ΔT₂
end


"""
    potential_temperature(param_set, T, ρ, q_tot=0, q_liq=0, q_ice=0)

The potential temperature.

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T`: temperature [K]
 - `ρ`: (moist-)air density [kg/m³]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `θ`: potential temperature [K]

If the specific humidities are not given, the result is for dry air.

Note: if any of `q_tot`, `q_liq`, or `q_ice` are nonzero, the Exner function uses the
moist exponent `R_m/cp_m` (reducing to the dry exponent `R_d/cp_d` in the dry limit).
"""
@inline function potential_temperature(
    param_set::APS,
    T,
    ρ,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    return T / exner(param_set, T, ρ, q_tot, q_liq, q_ice)
end

"""
    potential_temperature_given_pressure(param_set, T, p, q_tot=0, q_liq=0, q_ice=0)

The potential temperature.

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `p` pressure
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is for dry air.

Note: if any of `q_tot`, `q_liq`, or `q_ice` are nonzero, the Exner function uses the
moist exponent `R_m/cp_m` (reducing to the dry exponent `R_d/cp_d` in the dry limit).
"""
@inline function potential_temperature_given_pressure(
    param_set::APS,
    T,
    p,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    return T / exner_given_pressure(param_set, p, q_tot, q_liq, q_ice)
end

"""
    virtual_temperature(param_set, T, q_tot=0, q_liq=0, q_ice=0)

The virtual temperature.

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T`: temperature [K]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `T_v`: virtual temperature [K]

If the specific humidities are not given, the result is for dry air.
The virtual temperature is defined such that dry air at `T_v` has the same density
as moist air at `T`, i.e., `T_v = T (R_m / R_d)`.
"""
@inline function virtual_temperature(
    param_set::APS,
    T,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    R_d = TP.R_d(param_set)
    R_m = gas_constant_air(param_set, q_tot, q_liq, q_ice)
    return T * R_m / R_d
end

"""
    virtual_pottemp(param_set, T, ρ, q_tot=0, q_liq=0, q_ice=0)

The virtual potential temperature.

 - `param_set` thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `q_tot` total specific humidity
 - `q_liq` liquid specific humidity
 - `q_ice` ice specific humidity

If the specific humidities are not given, the result is for dry air.
"""
@inline function virtual_pottemp(
    param_set::APS,
    T,
    ρ,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    R_d = TP.R_d(param_set)
    R_m = gas_constant_air(param_set, q_tot, q_liq, q_ice)
    return potential_temperature(param_set, T, ρ, q_tot, q_liq, q_ice) * R_m / R_d
end

"""
    liquid_ice_pottemp(param_set, T, ρ, q_tot=0, q_liq=0, q_ice=0)

The liquid-ice potential temperature.

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T`: temperature [K]
 - `ρ`: (moist-)air density [kg/m³]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `θ_li`: liquid-ice potential temperature [K]

If the specific humidities are not given, the result is for dry air.

# Reference
Betts (1973), "Non-precipitating cumulus convection and its parameterization," 
*Quarterly Journal of the Royal Meteorological Society*, **99**, 178-196,
doi:[10.1002/qj.49709941915](https://doi.org/10.1002/qj.49709941915).
"""
@inline function liquid_ice_pottemp(
    param_set::APS,
    T,
    ρ,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    return liquid_ice_pottemp_given_pressure(
        param_set,
        T,
        air_pressure(param_set, T, ρ, q_tot, q_liq, q_ice),
        q_tot,
        q_liq,
        q_ice,
    )
end

"""
    liquid_ice_pottemp_given_pressure(param_set, T, p, q_tot=0, q_liq=0, q_ice=0)

The liquid-ice potential temperature.

# Arguments
 - `param_set`: thermodynamics parameter set, see the [`Thermodynamics`](@ref) for more details
 - `T`: temperature [K]
 - `p`: pressure [Pa]
 - `q_tot`: total specific humidity [kg/kg]
 - `q_liq`: liquid specific humidity [kg/kg]
 - `q_ice`: ice specific humidity [kg/kg]

# Returns
 - `θ_li`: liquid-ice potential temperature [K]

If the specific humidities are not given, the result is for dry air.
The latent heats of phase transitions are approximated as constants.

# Reference
Betts (1973), "Non-precipitating cumulus convection and its parameterization," 
*Quarterly Journal of the Royal Meteorological Society*, **99**, 178-196,
doi:[10.1002/qj.49709941915](https://doi.org/10.1002/qj.49709941915).
"""
@inline function liquid_ice_pottemp_given_pressure(
    param_set::APS,
    T,
    p,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    # liquid-ice potential temperature, approximating latent heats
    # of phase transitions as constants
    cpm = cp_m(param_set, q_tot, q_liq, q_ice)
    return potential_temperature_given_pressure(param_set, T, p, q_tot, q_liq, q_ice) * (
        1 - humidity_weighted_latent_heat(param_set, q_liq, q_ice) / (cpm * T)
    )
end
