export air_temperature_given_hq
export air_temperature_given_pρq
export air_temperature_given_pθq
export air_temperature_given_ρθq

"""
    air_temperature_given_hq(param_set, h, q_tot=0, q_liq=0, q_ice=0)

!!! warn "Deprecated"
    This function is deprecated and will be removed in a future release.
    Please use [`air_temperature`](@ref) with the [`phq`](@ref) type instead.
"""
@inline function air_temperature_given_hq(
    param_set::APS,
    h,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    return air_temperature(param_set, ph(), h, q_tot, q_liq, q_ice)
end

"""
    air_temperature_given_pρq(param_set, p, ρ, q_tot=0, q_liq=0, q_ice=0)

!!! warn "Deprecated"
    This function is deprecated and will be removed in a future release.
    Please use [`air_temperature`](@ref) with the [`pρq`](@ref) type instead.
"""
@inline function air_temperature_given_pρq(
    param_set::APS,
    p,
    ρ,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    return air_temperature(param_set, pρ(), p, ρ, q_tot, q_liq, q_ice)
end

"""
    air_temperature_given_pθq(param_set, p, θ_liq_ice, q_tot=0, q_liq=0, q_ice=0)

!!! warn "Deprecated"
    This function is deprecated and will be removed in a future release.
    Please use [`air_temperature`](@ref) with the [`pθ_liq_ice_q`](@ref) type instead.
"""
@inline function air_temperature_given_pθq(
    param_set::APS,
    p,
    θ_liq_ice,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    return air_temperature(
        param_set,
        pθ_li(),
        p,
        θ_liq_ice,
        q_tot,
        q_liq,
        q_ice,
    )
end

"""
    air_temperature_given_ρθq(param_set, ρ, θ_liq_ice, q_tot=0, q_liq=0, q_ice=0)

!!! warn "Deprecated"
    This function is deprecated and will be removed in a future release.
    Please use [`air_temperature`](@ref) with the [`ρθ_liq_ice_q`](@ref) type instead.
"""
@inline function air_temperature_given_ρθq(
    param_set::APS,
    ρ,
    θ_liq_ice,
    q_tot = 0,
    q_liq = 0,
    q_ice = 0,
)
    return air_temperature(
        param_set,
        ρθ_li(),
        ρ,
        θ_liq_ice,
        q_tot,
        q_liq,
        q_ice,
    )
end
