export air_temperature
export air_temperature_given_ρp
export dry_pottemp
export virtual_pottemp
export liquid_ice_pottemp
export liquid_ice_pottemp_sat
export virtual_temperature

"""
    air_temperature(param_set, e_int[, q::PhasePartition])

The air temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `e_int` specific internal energy of moist air

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air.
"""
@inline function air_temperature(
    param_set::APS,
    e_int,
    q::PhasePartition = q_pt_0(param_set),
    cvm::Number = cv_m(param_set, q),
)
    T_0 = TP.T_0(param_set)
    R_d = TP.R_d(param_set)
    e_int_v0 = TP.e_int_v0(param_set)
    e_int_i0 = TP.e_int_i0(param_set)
    return T_0 +
           (
        e_int - (q.tot - q.liq - q.ice) * e_int_v0 +
        q.ice * e_int_i0 +
        (1 - q.tot) * R_d * T_0
    ) / cvm
end

"""
    air_temperature_from_enthalpy(param_set, h[, q::PhasePartition]s)

The air temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `h` specific enthalpy of moist air

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air.
"""
@inline function air_temperature_from_enthalpy(
    param_set::APS,
    h,
    q::PhasePartition = q_pt_0(param_set),
)
    cp_m_ = cp_m(param_set, q)
    T_0 = TP.T_0(param_set)
    LH_v0 = TP.LH_v0(param_set)
    LH_f0 = TP.LH_f0(param_set)
    return T_0 + (h - (q.tot - q.liq - q.ice) * LH_v0 + q.ice * LH_f0) / cp_m_
end

"""
    air_temperature_given_ρp(param_set, p, ρ[, q::PhasePartition])

The air temperature, where

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `ρ` air density

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air.
"""
@inline function air_temperature_given_ρp(
    param_set::APS,
    p,
    ρ,
    q::PhasePartition = q_pt_0(param_set),
)
    R_m = gas_constant_air(param_set, q)
    return p / (R_m * ρ)
end

"""
    dry_pottemp(param_set, T, ρ[, q::PhasePartition])

The dry potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air.
"""
@inline function dry_pottemp(
    param_set::APS,
    T,
    ρ,
    q::PhasePartition = q_pt_0(param_set),
    cpm = cp_m(param_set, q),
)
    return T / exner(param_set, T, ρ, q, cpm)
end

"""
    dry_pottemp_given_pressure(param_set, T, p[, q::PhasePartition])

The dry potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `p` pressure

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air, i.e., using the adiabatic exponent 
 for dry air.
"""
@inline function dry_pottemp_given_pressure(
    param_set::APS,
    T,
    p,
    q::PhasePartition = q_pt_0(param_set),
    cpm = cp_m(param_set, q),
)
    return T / exner_given_pressure(param_set, p, q, cpm)
end

"""
    latent_heat_liq_ice(param_set[, q::PhasePartition])

Specific-humidity weighted sum of latent heats of liquid and ice evaluated at reference temperature 
`T_0`, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q` [`PhasePartition`](@ref). 

When `q` is not provided, `latent_heat_liq_ice` is zero.

 This is used in the evaluation of the liquid-ice potential temperature.
"""
@inline function latent_heat_liq_ice(
    param_set::APS,
    q::PhasePartition = q_pt_0(param_set),
)
    LH_v0 = TP.LH_v0(param_set)
    LH_s0 = TP.LH_s0(param_set)
    return LH_v0 * q.liq + LH_s0 * q.ice
end

"""
    liquid_ice_pottemp_given_pressure(param_set, T, p[, q::PhasePartition])

The liquid-ice potential temperature, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `p` pressure

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the result is the dry-air potential temperature.
"""
@inline function liquid_ice_pottemp_given_pressure(
    param_set::APS,
    T,
    p,
    q::PhasePartition = q_pt_0(param_set),
    cpm = cp_m(param_set, q),
)
    # liquid-ice potential temperature, approximating latent heats
    # of phase transitions as constants
    return dry_pottemp_given_pressure(param_set, T, p, q, cpm) *
           (1 - latent_heat_liq_ice(param_set, q) / (cpm * T))
end

"""
    liquid_ice_pottemp(param_set, T, ρ[, q::PhasePartition])

The liquid-ice potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the result is the dry-air potential temperature.
"""
@inline function liquid_ice_pottemp(
    param_set::APS,
    T,
    ρ,
    q::PhasePartition = q_pt_0(param_set),
    cpm = cp_m(param_set, q),
)
    return liquid_ice_pottemp_given_pressure(
        param_set,
        T,
        air_pressure(param_set, T, ρ, q),
        q,
        cpm,
    )
end

"""
    air_temperature_given_pθq(
        param_set,
        p,
        θ_liq_ice,
        [q::PhasePartition]
    )

The air temperature obtained by inverting the liquid-ice potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `θ_liq_ice` liquid-ice potential temperature
 
and, optionally,
 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the `θ_liq_ice` is assumed to be the dry-air potential temperature.
"""
@inline function air_temperature_given_pθq(
    param_set::APS,
    p,
    θ_liq_ice,
    q::PhasePartition = q_pt_0(param_set),
    cpm = cp_m(param_set, q),
)
    return θ_liq_ice * exner_given_pressure(param_set, p, q, cpm) +
           latent_heat_liq_ice(param_set, q) / cpm
end

"""
    air_temperature_given_ρθq(param_set, ρ, θ_liq_ice[, q::PhasePartition])

The air temperature obtained by inverting the liquid-ice potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `θ_liq_ice` liquid-ice potential temperature

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the results are for dry air.
"""
@inline function air_temperature_given_ρθq(
    param_set::APS,
    ρ,
    θ_liq_ice,
    q::PhasePartition = q_pt_0(param_set),
)

    p0 = TP.p_ref_theta(param_set)
    cvm = cv_m(param_set, q)
    cpm = cp_m(param_set, q)
    R_m = gas_constant_air(param_set, q)
    κ = 1 - cvm / cpm
    T_u = (ρ * R_m * θ_liq_ice / p0)^(R_m / cvm) * θ_liq_ice
    T_1 = latent_heat_liq_ice(param_set, q) / cvm
    T_2 = -κ / (2 * T_u) * (latent_heat_liq_ice(param_set, q) / cvm)^2
    return T_u + T_1 + T_2
end

"""
    liquid_ice_pottemp_sat(param_set, T, ρ, phase_type[, q::PhasePartition, cpm])

The saturated liquid ice potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `phase_type` a thermodynamic state type

and, optionally,

 - `q` [`PhasePartition`](@ref). 

When `q` is not provided, the air assumed to be dry.
"""
@inline function liquid_ice_pottemp_sat(
    param_set::APS,
    T,
    ρ,
    ::Type{phase_type},
    q::PhasePartition = q_pt_0(param_set),
    cpm = cp_m(param_set, q),
) where {phase_type <: ThermodynamicState}
    q_v_sat = q_vap_saturation(param_set, T, ρ, phase_type, q)
    return liquid_ice_pottemp(param_set, T, ρ, PhasePartition(q_v_sat), cpm)
end

"""
    liquid_ice_pottemp_sat(param_set, T, ρ, phase_type, q_tot)

The saturated liquid ice potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density
 - `phase_type` a thermodynamic state type
 - `q_tot` total specific humidity
"""
@inline function liquid_ice_pottemp_sat(
    param_set::APS,
    T,
    ρ,
    ::Type{phase_type},
    q_tot,
) where {phase_type <: ThermodynamicState}
    q = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    cpm = cp_m(param_set, q)
    return liquid_ice_pottemp(param_set, T, ρ, q, cpm)
end

"""
    virtual_pottemp(param_set, T, ρ[, q::PhasePartition])

The virtual potential temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `ρ` (moist-)air density

and, optionally,

 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the result is the dry-air potential temperature.

The virtual potential temperature is defined as `θ_v = (R_m/R_d) * θ`, where `θ` is the
potential temperature. It is the potential temperature a dry air parcel would need to have to
have the same density as the moist air parcel at the same pressure.
"""
@inline function virtual_pottemp(
    param_set::APS,
    T,
    ρ,
    q::PhasePartition = q_pt_0(param_set),
    cpm = cp_m(param_set, q),
)
    R_d = TP.R_d(param_set)
    return gas_constant_air(param_set, q) / R_d *
           dry_pottemp(param_set, T, ρ, q, cpm)
end

"""
    virtual_temperature(param_set, T[, q::PhasePartition])

The virtual temperature, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature

and, optionally,
 
 - `q` [`PhasePartition`](@ref). 
 
When `q` is not provided, the result is the regular temperature. 

The virtual temperature is defined as `T_v = (R_m/R_d) * T`. It is the temperature
a dry air parcel would need to have to have the same density as the moist air parcel
at the same pressure.
"""
@inline function virtual_temperature(
    param_set::APS,
    T,
    q::PhasePartition = q_pt_0(param_set),
)
    R_d = TP.R_d(param_set)
    return gas_constant_air(param_set, q) / R_d * T
end

"""
    temperature_and_humidity_given_TᵥρRH(param_set, T_virt, ρ, RH)

The air temperature and total specific humidity `q_tot`, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T_virt` virtual temperature
 - `ρ` air density
 - `RH` relative humidity
 - `phase_type` a thermodynamic state type
"""
@inline function temperature_and_humidity_given_TᵥρRH(
    param_set::APS{FT},
    T_virt,
    ρ,
    RH,
    ::Type{phase_type},
    maxiter::Int = 100,
    tol::RS.AbstractTolerance = RS.ResidualTolerance(sqrt(eps(FT))),
) where {FT, phase_type <: ThermodynamicState}

    T_init_min = TP.T_init_min(param_set)
    _T_max = T_virt
    @inline roots(T) =
        T_virt - virt_temp_from_RH(param_set, ReLU(T), ρ, RH, phase_type)
    sol = RS.find_zero(
        roots,
        RS.SecantMethod(T_init_min, _T_max),
        RS.CompactSolution(),
        tol,
        maxiter,
    )
    if !sol.converged
        if print_warning()
            print_warning_TᵥρRH(
                RS.SecantMethod,
                T_virt,
                RH,
                ρ,
                sol.root,
                maxiter,
                tol.tol,
            )
        end
        if error_on_non_convergence()
            error("Failed to converge with printed set of inputs.")
        end
    end
    T = sol.root

    # Re-compute specific humidity and phase partitioning
    # given the temperature
    q_tot = RH * q_vap_saturation(param_set, T, ρ, phase_type)
    q_pt = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    return (T, q_pt)

end

"""
    air_temperature_given_ρθq_nonlinear(param_set, ρ, θ_liq_ice, maxiter, tol, q::PhasePartition)

Computes temperature `T`, given

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `ρ` (moist-)air density
 - `θ_liq_ice` liquid-ice potential temperature

and, optionally,
 - `maxiter` maximum iterations for non-linear equation solve
- `tol` absolute tolerance for non-linear equation iterations. Can be one of:
    - `RelativeSolutionTolerance()` to stop when `|x_2 - x_1|/x_1 < tol`
    - `ResidualTolerance()` to stop when `|f(x)| < tol`
 - `q` [`PhasePartition`](@ref).When `q` is not provided, the results are for dry air,

The temperature `T` is found by finding the root of
`T - air_temperature_given_pθq(param_set,
                               air_pressure(param_set, T, ρ, q),
                               θ_liq_ice,
                               q) = 0`
"""
@inline function air_temperature_given_ρθq_nonlinear(
    param_set::APS,
    ρ,
    θ_liq_ice,
    maxiter::Int,
    tol::RS.AbstractTolerance,
    q::PhasePartition = q_pt_0(param_set),
)
    T_init_min = TP.T_init_min(param_set)
    _T_max = TP.T_max(param_set)
    @inline roots(T) =
        T - air_temperature_given_pθq(
            param_set,
            air_pressure(param_set, ReLU(T), ρ, q),
            θ_liq_ice,
            q,
        )
    sol = RS.find_zero(
        roots,
        RS.SecantMethod(T_init_min, _T_max),
        RS.CompactSolution(),
        tol,
        maxiter,
    )
    if !sol.converged
        if print_warning()
            print_warning_ρθq_nonlinear(
                RS.SecantMethod,
                θ_liq_ice,
                ρ,
                q.tot,
                q.liq,
                q.ice,
                sol.root,
                maxiter,
                tol.tol,
            )
        end
        if error_on_non_convergence()
            error("Failed to converge with printed set of inputs.")
        end
    end
    return sol.root
end
