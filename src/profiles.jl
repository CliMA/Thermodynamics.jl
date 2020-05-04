#### Reference state profiles

"""
    fixed_lapse_rate_ref_state(
        param_set::APS,
        z::FT,
        T_surface::FT,
        T_min::FT,
        ) where {FT <: AbstractFloat}

Fixed lapse rate hydrostatic reference state
"""
function fixed_lapse_rate_ref_state(
    param_set::APS,
    z::FT,
    T_surface::FT,
    T_min::FT,
) where {FT <: AbstractFloat}
    _grav::FT = grav(param_set)
    _cp_d::FT = cp_d(param_set)
    _R_d::FT = R_d(param_set)
    _MSLP::FT = MSLP(param_set)
    Γ = _grav / _cp_d
    z_tropopause = (T_surface - T_min) / Γ
    H_min = _R_d * T_min / _grav
    T = max(T_surface - Γ * z, T_min)
    p = _MSLP * (T / T_surface)^(_grav / (_R_d * Γ))
    T == T_min && (p = p * exp(-(z - z_tropopause) / H_min))
    ρ = p / (_R_d * T)
    return T, p, ρ
end

"""
    tested_profiles(param_set::APS, n::Int, ::Type{FT}) where {FT}

A range of input arguments to thermodynamic state constructors
 - `param_set` an `AbstractParameterSet`, see the [`MoistThermodynamics`](@ref) for more details
 - `e_int` internal energy
 - `ρ` (moist-)air density
 - `q_tot` total specific humidity
 - `q_pt` phase partition
 - `T` air temperature
 - `θ_liq_ice` liquid-ice potential temperature
that are tested for convergence in saturation adjustment.

Note that the output vectors are of size ``n*n_RH``, and they
should span the input arguments to all of the constructors.
"""
function tested_profiles(param_set::APS, n::Int, ::Type{FT}) where {FT}
    n_RS1 = 10
    n_RS2 = 20
    n_RS = n_RS1 + n_RS2
    z_range = range(FT(0), stop = FT(2.5e4), length = n)
    relative_sat1 = range(FT(0), stop = FT(1), length = n_RS1)
    relative_sat2 = range(FT(1), stop = FT(1.02), length = n_RS2)
    relative_sat = [relative_sat1..., relative_sat2...]
    T_min = FT(150)
    T_surface = FT(350)

    args =
        fixed_lapse_rate_ref_state.(
            Ref(param_set),
            z_range,
            Ref(T_surface),
            Ref(T_min),
        )
    T, p, ρ = getindex.(args, 1), getindex.(args, 2), getindex.(args, 3)

    p = collect(Iterators.flatten([p for RS in 1:n_RS]))
    ρ = collect(Iterators.flatten([ρ for RS in 1:n_RS]))
    T = collect(Iterators.flatten([T for RS in 1:n_RS]))
    relative_sat = collect(Iterators.flatten([relative_sat for RS in 1:n]))

    # Additional variables
    q_sat = q_vap_saturation.(Ref(param_set), T, ρ)
    q_tot = min.(relative_sat .* q_sat, FT(1))
    q_pt = PhasePartition_equil.(Ref(param_set), T, ρ, q_tot)
    e_int = internal_energy.(Ref(param_set), T, q_pt)
    θ_liq_ice = liquid_ice_pottemp.(Ref(param_set), T, ρ, q_pt)
    return e_int, ρ, q_tot, q_pt, T, p, θ_liq_ice
end
