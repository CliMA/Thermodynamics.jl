"""
# Common setup for functional Thermodynamics.jl tests

This file contains common imports, parameter sets, and tolerances for tests that
exercise the functional API.
"""

using Test
using Thermodynamics
import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import RootSolvers as RS
import ForwardDiff

# For generating physically plausible input arrays (no ThermodynamicState types needed)
using .TestedProfiles

# Test both Float32 and Float64
array_types = [Array{Float32}, Array{Float64}]

# Reference parameter sets
param_set_Float64 = TP.ThermodynamicsParameters(Float64)
param_set_Float32 = TP.ThermodynamicsParameters(Float32)

# Saturation adjustment tolerance (relative change of temperature between consecutive iterations)
rtol_temperature = 1e-6

# Generic tolerances used throughout the functional tests
atol_temperature = 1e-1
rtol_temperature_fd = 5e-4


