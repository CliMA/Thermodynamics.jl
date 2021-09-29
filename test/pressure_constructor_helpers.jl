#=
These helper functions ensure that the density-formulation
and pressure formulation are consistent with one another.
The reason we have a separate pressure formulation in `src/`
is because it is more performant--but both strategies should
produce the same results within machine precision.
=#


#=
    air_density_equil(
        param_set,
        T,
        p,
        q_tot,
    )
Compute the (equilibrium-)air density `ρ` that is consistent with
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `p` pressure
 - `q_tot` total specific humidity
by finding the root of
```
ρ - air_density(
    param_set,
    T,
    p,
    PhasePartition_equil(param_set, T, ρ, q_tot, PhaseEquil)
)
```
This air density is in sync with that computed from `air_density`,
which relies on the full `PhasePartition`.
See also [`saturation_adjustment`](@ref).
!!! warn
    This is an expensive function, and has internalized
    iterative solver parameters to simplify the interface.
=#
function air_density_equil(
    param_set::APS,
    T::FT,
    p::FT,
    q_tot::FT,
) where {FT <: Real}
    # Assume all liquid
    ρ_init = TD.air_density(param_set, T, p, TD.PhasePartition(q_tot))
    phase_type = PhaseEquil
    maxiter = 50
    tol = SolutionTolerance(eps(FT))
    sol = find_zero(
        ρ -> begin
            _q_tot = oftype(ρ, q_tot)
            _p = oftype(ρ, p)

            q_pt = TD.PhasePartition_equil(
                param_set,
                oftype(ρ, T),
                ρ,
                oftype(ρ, q_tot),
                phase_type,
            )
            # TODO: is one of these equations preferred over the other?
            # (e.g., physically / numerically)
            # T - air_temperature_from_ideal_gas_law(param_set, _p, ρ, q_pt)
            ρ - TD.air_density(param_set, oftype(ρ, T), _p, q_pt)
        end,
        NewtonsMethodAD(ρ_init),
        CompactSolution(),
        tol,
        maxiter,
    )
    if !sol.converged
        if TD.print_warning()
            @print("-----------------------------------------\n")
            @print("maxiter reached in air_density_equil:\n")
            @print("    Method=NewtonsMethodAD")
            @print(", T=", T)
            @print(", p=", p)
            @print(", q_tot=", q_tot)
            @print(", ρ=", sol.root)
            @print(", maxiter=", maxiter)
            @print(", tol=", tol.tol, "\n")
        end
        if TD.error_on_non_convergence()
            error("Failed to converge with printed set of inputs.")
        end
    end
    return sol.root
end

function saturation_adjustment_given_pθq_unified(
    param_set::APS,
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    phase_type::Type{<:PhaseEquil},
    maxiter::Int,
    tol::RS.AbstractTolerance,
) where {FT <: Real}
    _T_min::FT = T_min(param_set)
    air_temp(q) = TD.air_temperature_given_pθq(param_set, p, θ_liq_ice, q)
    T_1 = air_temp(TD.PhasePartition(q_tot)) # Assume all vapor
    ρ_T(T) = air_density_equil(param_set, TD.heavisided(T), p, q_tot)
    ρ = ρ_T(T_1)
    q_v_sat = TD.q_vap_saturation(param_set, T_1, ρ, phase_type)
    unsaturated = q_tot <= q_v_sat
    if unsaturated && T_1 > _T_min
        return T_1
    end
    T_2 = air_temp(TD.PhasePartition(q_tot, FT(0), q_tot)) # Assume all ice
    T_2 = TD.bound_upper_temperature(T_1, T_2)
    sol = find_zero(
        T ->
            liquid_ice_pottemp_sat(
                param_set,
                TD.heavisided(T),
                ρ_T(T),
                phase_type,
                q_tot,
            ) - θ_liq_ice,
        SecantMethod(T_1, T_2),
        CompactSolution(),
        tol,
        maxiter,
    )
    if !sol.converged
        if TD.print_warning()
            @print("-----------------------------------------\n")
            @print("maxiter reached in saturation_adjustment_given_pθq_unified:\n")
            @print("    Method=SecantMethod")
            @print(", p=", p)
            @print(", θ_liq_ice=", θ_liq_ice)
            @print(", q_tot=", q_tot)
            @print(", T=", sol.root)
            @print(", maxiter=", maxiter)
            @print(", tol=", tol.tol, "\n")
        end
        if TD.error_on_non_convergence()
            error("Failed to converge with printed set of inputs.")
        end
    end
    return sol.root
end

#=
This constructor re-computes pressure on the fly,
and uses all of the density-formulated equations
=#
function PhaseEquil_pθq_unified(
    param_set::APS,
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    maxiter::IT = nothing,
    temperature_tol::FTT = nothing,
) where {FT <: Real, IT <: TD.ITERTYPE, FTT <: TD.TOLTYPE(FT)}
    # NOTE THAT maxiter and temperature_tol should be in sync
    # with the PhaseEquil_pθq in `src/`
    maxiter === nothing && (maxiter = 50)
    temperature_tol === nothing && (temperature_tol = FT(1e-3))
    phase_type = PhaseEquil
    q_tot_safe = clamp(q_tot, FT(0), FT(1))
    tol = ResidualTolerance(temperature_tol)
    T = saturation_adjustment_given_pθq_unified(
        param_set,
        p,
        θ_liq_ice,
        q_tot_safe,
        phase_type,
        maxiter,
        tol,
    )
    ρ = air_density_equil(param_set, T, p, q_tot_safe)
    q_pt = TD.PhasePartition_equil(param_set, T, ρ, q_tot_safe, phase_type)
    e_int = internal_energy(param_set, T, q_pt)
    return PhaseEquil{FT, typeof(param_set)}(
        param_set,
        ρ,
        p,
        e_int,
        q_tot_safe,
        T,
    )
end
