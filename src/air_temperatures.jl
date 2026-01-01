# Temperature functions
export virtual_temperature

"""
    virtual_temperature(param_set, T, q_tot=0, q_liq=0, q_ice=0)

The virtual temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `q_tot` total specific humidity of water
 - `q_liq` specific humidity of liquid
 - `q_ice` specific humidity of ice

If the specific humidities are not given, the result is for dry air.
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
