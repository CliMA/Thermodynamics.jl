
# Universal constants
MoistThermodynamics.Parameters.gas_constant() = 8.3144598
MoistThermodynamics.Parameters.astro_unit() = 1.4959787e11

# Properties of dry air
MoistThermodynamics.Parameters.molmass_dryair(ps::ParameterSet) = 28.97e-3
MoistThermodynamics.Parameters.R_d(ps::ParameterSet)            = gas_constant() / molmass_dryair(ps)
MoistThermodynamics.Parameters.kappa_d(ps::ParameterSet)        = 2 / 7
MoistThermodynamics.Parameters.cp_d(ps::ParameterSet)           = R_d(ps) / kappa_d(ps)
MoistThermodynamics.Parameters.cv_d(ps::ParameterSet)           = cp_d(ps) - R_d(ps)

# Properties of water
MoistThermodynamics.Parameters.molmass_water(ps::ParameterSet)  = 18.01528e-3
MoistThermodynamics.Parameters.molmass_ratio(ps::ParameterSet)  = molmass_dryair(ps) / molmass_water(ps)
MoistThermodynamics.Parameters.R_v(ps::ParameterSet)            = gas_constant() / molmass_water(ps)
MoistThermodynamics.Parameters.cp_v(ps::ParameterSet)           = 1859
MoistThermodynamics.Parameters.cp_l(ps::ParameterSet)           = 4181
MoistThermodynamics.Parameters.cp_i(ps::ParameterSet)           = 2100
MoistThermodynamics.Parameters.cv_v(ps::ParameterSet)           = cp_v(ps) - R_v(ps)
MoistThermodynamics.Parameters.cv_l(ps::ParameterSet)           = cp_l(ps)
MoistThermodynamics.Parameters.cv_i(ps::ParameterSet)           = cp_i(ps)
MoistThermodynamics.Parameters.T_freeze(ps::ParameterSet)       = 273.15
MoistThermodynamics.Parameters.T_min(ps::ParameterSet)          = 150.0
MoistThermodynamics.Parameters.T_max(ps::ParameterSet)          = 1000.0
MoistThermodynamics.Parameters.T_icenuc(ps::ParameterSet)       = 233.00
MoistThermodynamics.Parameters.T_triple(ps::ParameterSet)       = 273.16
MoistThermodynamics.Parameters.T_0(ps::ParameterSet)            = T_triple(ps)
MoistThermodynamics.Parameters.LH_v0(ps::ParameterSet)          = 2.5008e6
MoistThermodynamics.Parameters.LH_s0(ps::ParameterSet)          = 2.8344e6
MoistThermodynamics.Parameters.LH_f0(ps::ParameterSet)          = LH_s0(ps) - LH_v0(ps)
MoistThermodynamics.Parameters.e_int_v0(ps::ParameterSet)       = LH_v0(ps) - R_v(ps) * T_0(ps)
MoistThermodynamics.Parameters.e_int_i0(ps::ParameterSet)       = LH_f0(ps)
MoistThermodynamics.Parameters.press_triple(ps::ParameterSet)   = 611.657

# Planetary parameters
MoistThermodynamics.Parameters.grav(ps::ParameterSet)           = 9.81
MoistThermodynamics.Parameters.MSLP(ps::ParameterSet)           = 1.01325e5
