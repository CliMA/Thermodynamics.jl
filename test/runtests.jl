#=
In the event that CLIMA's version of ParameterSet diverges
from the one defined in MoistThermodynamics's test suite,
we expose MoistThermodynamics's test suite within the package
itself so that we can test MoistThermodynamics with CLIMA's
ParameterSet.
=#

import MoistThermodynamics
struct ParameterSet end
include("ParameterSet.jl")
using MoistThermodynamics.Parameters

# include("testdata.jl")

param_set = ParameterSet()

MoistThermodynamics.TestMoistThermo.test_moist_thermo(param_set)
