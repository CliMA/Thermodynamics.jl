"""
# Common Parameters, Includes, and Functions for Thermodynamics Tests

This file contains common thermodynamic parameter definitions and tolerances
that are shared across multiple test files to avoid duplication.
"""

using Test
using Thermodynamics
import Thermodynamics as TD
import Thermodynamics.Parameters as TP
using Thermodynamics.TestedProfiles
import RootSolvers as RS
using LinearAlgebra
import ForwardDiff

# Test both Float32 and Float64 for type stability
array_types = [Array{Float32}, Array{Float64}]

# Reference parameter set for Float64
param_set_Float64 = TP.ThermodynamicsParameters(Float64)
param_set_Float32 = TP.ThermodynamicsParameters(Float32)

# Saturation adjustment tolerance (relative change of temperature between consecutive iterations)
rtol_temperature = 1e-4

# Tolerances for tested quantities:
atol_temperature = 0.3   # Expected absolute temperature accuracy
atol_energy_temperature = TP.cv_d(param_set_Float64) * atol_temperature  # Expected absolute energy accuracy due to temperature accuracy
rtol_humidity = 1e-2     # Relative accuracy of specific humidity (for energy tolerance adjustments)
rtol_density = 1e-3      # Relative density accuracy
rtol_pressure = 1e-3     # Relative pressure accuracy

# Helper function to compare moisture content between thermodynamic states
compare_moisture(param_set, a::ThermodynamicState, b::ThermodynamicState) =
    compare_moisture(param_set, a, PhasePartition(param_set, b))

# Compare total moisture for equilibrium states
compare_moisture(param_set, ts::PhaseEquil, q_pt::PhasePartition) =
    getproperty(PhasePartition(param_set, ts), :tot) ≈ getproperty(q_pt, :tot)

# Compare all moisture components for non-equilibrium states
compare_moisture(param_set, ts::PhaseNonEquil, q_pt::PhasePartition) = all((
    getproperty(PhasePartition(param_set, ts), :tot) ≈ getproperty(q_pt, :tot),
    getproperty(PhasePartition(param_set, ts), :liq) ≈ getproperty(q_pt, :liq),
    getproperty(PhasePartition(param_set, ts), :ice) ≈ getproperty(q_pt, :ice),
))

# Function to extract thermodynamic parameters for testing
function extract_thermodynamic_parameters(param_set)
    FT = eltype(param_set)

    # Gas constants
    _R_d = FT(TP.R_d(param_set))           # Dry air gas constant
    _Rv_over_Rd = FT(TP.Rv_over_Rd(param_set))  # Vapor/dry air gas constant ratio
    _R_v = FT(TP.R_v(param_set))           # Water vapor gas constant

    # Heat capacities
    _cp_d = FT(TP.cp_d(param_set))         # Dry air specific heat at constant pressure
    _cp_v = FT(TP.cp_v(param_set))         # Water vapor specific heat at constant pressure
    _cp_l = FT(TP.cp_l(param_set))         # Liquid water specific heat
    _cp_i = FT(TP.cp_i(param_set))         # Ice specific heat
    _cv_d = FT(TP.cv_d(param_set))         # Dry air specific heat at constant volume
    _cv_v = FT(TP.cv_v(param_set))         # Water vapor specific heat at constant volume
    _cv_l = FT(TP.cv_l(param_set))         # Liquid water specific heat at constant volume
    _cv_i = FT(TP.cv_i(param_set))         # Ice specific heat at constant volume

    # Reference temperatures and energies
    _T_0 = FT(TP.T_0(param_set))           # Reference temperature
    _e_int_v0 = FT(TP.e_int_v0(param_set)) # Reference vapor internal energy
    _e_int_i0 = FT(TP.e_int_i0(param_set)) # Reference ice internal energy

    # Latent heats
    _LH_v0 = FT(TP.LH_v0(param_set))       # Latent heat of vaporization
    _LH_s0 = FT(TP.LH_s0(param_set))       # Latent heat of sublimation
    _LH_f0 = FT(TP.LH_f0(param_set))       # Latent heat of fusion

    # Triple point and phase transition temperatures
    _press_triple = FT(TP.press_triple(param_set))  # Triple point pressure
    _T_triple = FT(TP.T_triple(param_set))          # Triple point temperature
    _T_freeze = FT(TP.T_freeze(param_set))          # Freezing temperature
    _T_icenuc = FT(TP.T_icenuc(param_set))          # Homogeneous ice nucleation temperature

    # Temperature and pressure ranges
    _T_min = FT(TP.T_min(param_set))       # Minimum temperature
    _T_max = FT(TP.T_max(param_set))       # Maximum temperature
    _p_ref_theta = FT(TP.p_ref_theta(param_set))  # Reference pressure for potential temperature
    _kappa_d = FT(TP.kappa_d(param_set))   # Dry air adiabatic exponent

    return (
        _R_d,
        _Rv_over_Rd,
        _R_v,
        _cp_d,
        _cp_v,
        _cp_l,
        _cp_i,
        _cv_d,
        _cv_v,
        _cv_l,
        _cv_i,
        _T_0,
        _e_int_v0,
        _e_int_i0,
        _LH_v0,
        _LH_s0,
        _LH_f0,
        _press_triple,
        _T_triple,
        _T_freeze,
        _T_icenuc,
        _T_min,
        _T_max,
        _p_ref_theta,
        _kappa_d,
    )
end
