import ForwardDiff
import ClimaParams as CP

import Thermodynamics as TD
include("../../test/TestedProfiles.jl")
import .TestedProfiles
import Thermodynamics.Parameters as TP

ArrayType = Array{Float64}
FT = eltype(ArrayType)

param_set = TP.ThermodynamicsParameters(FT)
profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
(; T, p, e_int, ρ, θ_liq_ice) = profiles
(; q_tot, q_liq, q_ice, RH, e_kin, e_pot) = profiles

k = findfirst(q -> q > 0.01, q_tot) # test for one value with q_tot above some threshhold
ts_sol = TD.PhaseEquil_ρTq(param_set, ρ[k], T[k], q_tot[k])

function q_vap_sat(_T::FT) where {FT}
    _ρ = TD.air_density(param_set, ts_sol)
    _q_tot = TD.total_specific_humidity(param_set, ts_sol)
    _phase_type = TD.PhaseEquil{FT}
    _q_pt = TD.PhasePartition_equil(
        param_set,
        _T,
        oftype(_T, _ρ),
        oftype(_T, _q_tot),
        _phase_type,
    )
    return TD.q_vap_saturation(
        param_set,
        _T,
        oftype(_T, _ρ),
        _phase_type,
        _q_pt,
    )
end

function ∂q_vap_sat_∂T_vs_T(_T::FT) where {FT}
    _ρ = TD.air_density(param_set, ts_sol)
    _q_tot = TD.total_specific_humidity(param_set, ts_sol)
    _λ = TD.liquid_fraction(param_set, ts_sol)
    _phase_type = TD.PhaseEquil{FT}
    _q_pt = TD.PhasePartition_equil(
        param_set,
        _T,
        oftype(_T, _ρ),
        oftype(_T, _q_tot),
        _phase_type,
    )
    _q_vap_sat = TD.q_vap_saturation(param_set, _T, _ρ, _phase_type, _q_pt)
    return TD.∂q_vap_sat_∂T(
        param_set,
        oftype(_T, _λ),
        _T,
        oftype(_T, _q_vap_sat),
    )
end

∂q_vap_sat_∂T_fd = _T -> ForwardDiff.derivative(q_vap_sat, _T)

_ρ = TD.air_density(param_set, ts_sol)
_q_tot = TD.total_specific_humidity(param_set, ts_sol)

T_sorted = sort(T)
import Plots
p1 = Plots.plot()
Plots.plot!(
    T_sorted,
    ∂q_vap_sat_∂T_fd.(T_sorted);
    label = "∂qvsat_∂T ForwardDiff",
    yaxis = :log,
)
Plots.plot!(
    T_sorted,
    ∂q_vap_sat_∂T_vs_T.(T_sorted);
    label = "∂qvsat_∂T Analytic",
    yaxis = :log,
)
Plots.plot!(; xlabel = "T [K]", legend = :topleft)

p2 = Plots.plot()
Plots.plot!(
    T_sorted,
    ∂q_vap_sat_∂T_fd.(T_sorted) .- ∂q_vap_sat_∂T_vs_T.(T_sorted);
    label = "error",
)
Plots.plot!(; xlabel = "T [K]", legend = :topleft)

p3 = Plots.plot()
Plots.plot!(T_sorted, q_vap_sat.(T_sorted); label = "q_vap_sat")
Plots.plot!(; xlabel = "T [K]")
ρq_vals = "ρ=$_ρ, q_tot=$_q_tot"
Plots.plot(
    p1,
    p2,
    p3;
    layout = Plots.grid(3, 1),
    plot_title = "$ρq_vals",
    titlefontsizes = 5,
)

Plots.savefig("Clausius_Clapeyron.svg")
