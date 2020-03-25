module Parameters

export gas_constant,
       astro_unit,
       molmass_dryair,
       R_d,
       kappa_d,
       cp_d,
       cv_d,
       molmass_water,
       molmass_ratio,
       R_v,
       cp_v,
       cp_l,
       cp_i,
       cv_v,
       cv_l,
       cv_i,
       T_freeze,
       T_min,
       T_max,
       T_icenuc,
       T_triple,
       T_0,
       LH_v0,
       LH_s0,
       LH_f0,
       e_int_v0,
       e_int_i0,
       press_triple,
       grav,
       MSLP

# Universal constants
function gas_constant end
function astro_unit end

# Properties of dry air
function molmass_dryair end
function R_d end
function kappa_d end
function cp_d end
function cv_d end

# Properties of water
function molmass_water end
function molmass_ratio end
function R_v end
function cp_v end
function cp_l end
function cp_i end
function cv_v end
function cv_l end
function cv_i end
function T_freeze end
function T_min end
function T_max end
function T_icenuc end
function T_triple end
function T_0 end
function LH_v0 end
function LH_s0 end
function LH_f0 end
function e_int_v0 end
function e_int_i0 end
function press_triple end

# Planetary parameters
function grav end
function MSLP end

end