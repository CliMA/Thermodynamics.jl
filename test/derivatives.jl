using Test

import ForwardDiff

using Thermodynamics
import Thermodynamics as TD
import Thermodynamics.Parameters as TP
using Thermodynamics.TemperatureProfiles
using Thermodynamics.TestedProfiles
import ClimaParams as CP

FT = Float32
param_set = TP.ThermodynamicsParameters(FT)

_cp_d = FT(TP.cp_d(param_set))
_cp_v = FT(TP.cp_v(param_set))
_cp_l = FT(TP.cp_l(param_set))
_cp_i = FT(TP.cp_i(param_set))
_cv_d = FT(TP.cv_d(param_set))
_cv_v = FT(TP.cv_v(param_set))
_cv_l = FT(TP.cv_l(param_set))
_cv_i = FT(TP.cv_i(param_set))
_R_d = FT(TP.R_d(param_set))
_R_v = FT(TP.R_v(param_set))

_T_0 = FT(TP.T_0(param_set))
_e_int_v0 = FT(TP.e_int_v0(param_set))
_e_int_i0 = FT(TP.e_int_i0(param_set))

Δcv_v = _cv_v - _cv_d
ΔR_v = _R_v - _R_d

ρ = FT(1.0)
T = FT(280.0)
q_tot = FT(-1e-12)
#q_pt = PhasePartition(ForwardDiff.Dual(q_tot), FT(0), FT(0))

e_int = TD.internal_energy(param_set, T, PhasePartition(q_tot))
ts = TD.PhaseEquil_ρeq(param_set, ρ, e_int, q_tot)
ts_nonequil = TD.PhaseNonEquil(param_set, e_int, ρ, PhasePartition(q_tot))

ts_dual = TD.PhaseEquil_ρeq(param_set, ρ, e_int, ForwardDiff.Dual(q_tot, FT(1.0)))
@info "PhaseEquil"
#@show ts_dual
@show ForwardDiff.value(ts_dual.p)
@show ForwardDiff.partials(ts_dual.p)
# @show ForwardDiff.value(ts_dual.T)
# @show ForwardDiff.partials(ts_dual.T)

@info "Expected PhaseEquil"

∂e_int_∂q_tot = _T_0 * (Δcv_v - _R_d) - _e_int_v0
kappa_m = TD.gas_constant_air(param_set, ts) / TD.cv_m(param_set, ts)
∂kappa_m∂q_tot = (
    ΔR_v * TD.cv_m(param_set, ts) -
    Δcv_v * TD.gas_constant_air(param_set, ts)
) / abs2(TD.cv_m(param_set, ts))

∂p_∂q_tot_expected = kappa_m * ∂e_int_∂q_tot +
    ∂kappa_m∂q_tot * (
        _cp_d * _T_0 + e_int +
        ∂e_int_∂q_tot * q_tot
    )
@show ∂p_∂q_tot_expected

term1 = -(_T_0 * _R_d) - _e_int_v0
∂T_∂q_tot_expected = 1 / TD.cv_m(param_set, ts) * term1 -
    (Δcv_v / abs2(TD.cv_m(param_set, ts))) * (
        _R_d * _T_0 + e_int + term1 * q_tot
    )
# @show ∂T_∂q_tot_expected
p_nonequil = TD.air_pressure(param_set, TD.PhaseNonEquil(param_set, e_int, ρ, PhasePartition(ForwardDiff.Dual(q_tot, FT(1.0)))))
t_nonequil = TD.air_temperature(param_set, TD.PhaseNonEquil(param_set, e_int, ρ, PhasePartition(ForwardDiff.Dual(q_tot, FT(1.0)))))
@info "PhaseNonEquil"
@show ForwardDiff.value(p_nonequil)
@show ForwardDiff.partials(p_nonequil)
@show ForwardDiff.value(t_nonequil)
@show ForwardDiff.partials(t_nonequil)

@info "Expected PhaseNonEquil"

kappa_m = TD.gas_constant_air(param_set, ts_nonequil) / TD.cv_m(param_set, ts_nonequil)
∂kappa_m∂q_tot = (
    ΔR_v * TD.cv_m(param_set, ts_nonequil) -
    Δcv_v * TD.gas_constant_air(param_set, ts_nonequil)
) / abs2(TD.cv_m(param_set, ts_nonequil))

∂p_∂q_tot_expected = kappa_m * ∂e_int_∂q_tot +
    ∂kappa_m∂q_tot * (
        _cp_d * _T_0 + e_int +
        ∂e_int_∂q_tot * q_tot
    )
@show ∂p_∂q_tot_expected

∂T_∂q_tot_expected = 1 / TD.cv_m(param_set, ts_nonequil) * term1 -
    (Δcv_v / abs2(TD.cv_m(param_set, ts_nonequil))) * (
        _R_d * _T_0 + e_int + term1 * q_tot
    )
@show ∂T_∂q_tot_expected

nothing