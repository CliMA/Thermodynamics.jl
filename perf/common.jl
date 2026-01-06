using Test

import BenchmarkTools
import RootSolvers as RS

import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import ClimaParams as CP

const FT = Float64
const param_set = TP.ThermodynamicsParameters(FT)

include(joinpath(@__DIR__, "..", "test", "TestedProfiles.jl"))
using .TestedProfiles

"""
    functional_inputs(param_set, ::Type{FT}; n=256)

Construct equilibrium-consistent scalar inputs for functional benchmarks, derived from
`test/TestedProfiles.jl` (no deprecated PhasePartition/state types).
"""
function functional_inputs(param_set, ::Type{FT}; n = 256) where {FT}
    ps = TestedProfiles.PhaseEquilProfiles(param_set, Array{FT})
    n = min(n, length(ps.z))
    sl = 1:n

    T = ps.T[sl]
    ρ = ps.ρ[sl]
    p = ps.p[sl]
    q_tot = ps.q_tot[sl]
    q_liq = ps.q_liq[sl]
    q_ice = ps.q_ice[sl]

    # ρ-based thermodynamic variables (consistent with (T,ρ,q_tot))
    e_int_ρ = TD.internal_energy_sat.(Ref(param_set), T, ρ, q_tot)
    p_ρ = TD.air_pressure.(Ref(param_set), T, ρ, q_tot, q_liq, q_ice)
    θ_ρ = TD.liquid_ice_pottemp.(Ref(param_set), T, ρ, q_tot, q_liq, q_ice)

    # p-based thermodynamic variables (consistent with (T,p,q_tot), using ρ(T,p,q_tot))
    ρ_p = TD.air_density.(Ref(param_set), T, p, q_tot)
    q_liq_p, q_ice_p = let qi = TD.condensate_partition.(Ref(param_set), T, ρ_p, q_tot)
        (first.(qi), last.(qi))
    end
    e_int_p = TD.internal_energy_sat.(Ref(param_set), T, ρ_p, q_tot)
    h_p = TD.enthalpy_sat.(Ref(param_set), T, ρ_p, q_tot)
    θ_p =
        TD.liquid_ice_pottemp_given_pressure.(Ref(param_set), T, p, q_tot, q_liq_p, q_ice_p)

    return (;
        T,
        ρ,
        p,
        q_tot,
        q_liq,
        q_ice,
        e_int_ρ,
        p_ρ,
        θ_ρ,
        e_int_p,
        h_p,
        θ_p,
        q_liq_p,
        q_ice_p,
    )
end

const inputs = functional_inputs(param_set, FT)

function sa_ρeq_loop(param_set, ρ, e_int, q_tot; maxiter = 40, tol = FT(1e-10))
    s = zero(FT)
    @inbounds for i in eachindex(q_tot)
        (T, ql, qi) = TD.saturation_adjustment(
            RS.NewtonsMethod,
            param_set,
            TD.ρeq(),
            ρ[i],
            e_int[i],
            q_tot[i],
            maxiter,
            tol,
        )
        s += T + ql + qi
    end
    return s
end

function sa_pθ_loop(param_set, p, θ_liq_ice, q_tot; maxiter = 40, tol = FT(1e-10))
    s = zero(FT)
    @inbounds for i in eachindex(q_tot)
        (T, ql, qi) = TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            TD.pθ_liq_ice_q(),
            p[i],
            θ_liq_ice[i],
            q_tot[i],
            maxiter,
            tol,
        )
        s += T + ql + qi
    end
    return s
end

function eos_loop(param_set, T, ρ, q_tot, q_liq, q_ice)
    s = zero(FT)
    @inbounds for i in eachindex(q_tot)
        p = TD.air_pressure(param_set, T[i], ρ[i], q_tot[i], q_liq[i], q_ice[i])
        ρ_chk = TD.air_density(param_set, T[i], p, q_tot[i], q_liq[i], q_ice[i])
        Π = TD.exner_given_pressure(param_set, p, q_tot[i], q_liq[i], q_ice[i])
        s += p + ρ_chk + Π
    end
    return s
end
