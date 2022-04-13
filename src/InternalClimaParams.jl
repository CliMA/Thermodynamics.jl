module InternalClimaParams

import CLIMAParameters
const CP = CLIMAParameters
const CPP = CP.Planet
const APS = CP.AbstractParameterSet

T_min(param_set::APS) = CPP.T_min(param_set)
T_max(param_set::APS) = CPP.T_max(param_set)
MSLP(param_set::APS) = CPP.MSLP(param_set)
R_d(param_set::APS) = CPP.R_d(param_set)
cp_d(param_set::APS) = CPP.cp_d(param_set)
kappa_d(param_set::APS) = CPP.kappa_d(param_set)
molmass_ratio(param_set::APS) = CPP.molmass_ratio(param_set)
cp_v(param_set::APS) = CPP.cp_v(param_set)
cp_l(param_set::APS) = CPP.cp_l(param_set)
cp_i(param_set::APS) = CPP.cp_i(param_set)
cv_d(param_set::APS) = CPP.cv_d(param_set)
cv_v(param_set::APS) = CPP.cv_v(param_set)
cv_l(param_set::APS) = CPP.cv_l(param_set)
cv_i(param_set::APS) = CPP.cv_i(param_set)
T_0(param_set::APS) = CPP.T_0(param_set)
e_int_v0(param_set::APS) = CPP.e_int_v0(param_set)
e_int_i0(param_set::APS) = CPP.e_int_i0(param_set)
LH_v0(param_set::APS) = CPP.LH_v0(param_set)
LH_s0(param_set::APS) = CPP.LH_s0(param_set)
LH_f0(param_set::APS) = CPP.LH_f0(param_set)
press_triple(param_set::APS) = CPP.press_triple(param_set)
R_v(param_set::APS) = CPP.R_v(param_set)
T_triple(param_set::APS) = CPP.T_triple(param_set)
T_freeze(param_set::APS) = CPP.T_freeze(param_set)
T_icenuc(param_set::APS) = CPP.T_icenuc(param_set)
pow_icenuc(param_set::APS) = CPP.pow_icenuc(param_set)
entropy_reference_temperature(param_set::APS) =
    CPP.entropy_reference_temperature(param_set)
entropy_dry_air(param_set::APS) = CPP.entropy_dry_air(param_set)
entropy_water_vapor(param_set::APS) = CPP.entropy_water_vapor(param_set)
T_surf_ref(param_set::APS) = CPP.T_surf_ref(param_set)
T_min_ref(param_set::APS) = CPP.T_min_ref(param_set)
grav(param_set::APS) = CPP.grav(param_set)

end
