# Thermodynamics.jl Test Suite
# 
# This file runs the complete test suite for the Thermodynamics package.
# Tests are organized into the following categories:
# 
# 1. Aqua tests: Code quality and performance checks (disabled due to system issues)
# 2. TemperatureProfiles: Physical consistency of temperature profiles
# 3. Core thermodynamic tests: Split into functional categories
# 4. Data tests: Robustness testing with realistic atmospheric data
# 5. GPU tests: GPU compatibility and kernel testing (run separately)
#
# To run GPU tests: julia --project=. test/runtests_gpu.jl CuArray

# Run aqua tests for code quality and performance (fail fast)
include("aqua.jl")

# Run temperature profile tests for physical consistency
include("TemperatureProfiles.jl")
include("TestedProfiles.jl")

# Functional tests
include("test_common.jl")

# Run core thermodynamic function tests (split into functional categories)

include("dry_adiabatic_processes.jl")
include("correctness.jl")
include("exceptions.jl")
include("saturation_adjustment.jl")
include("convergence_saturation_adjustment.jl")
include("default_saturation_adjustment.jl")
include("type_stability.jl")
include("optimization_tests.jl")
include("ad_tests.jl")
include("liquid_fraction_ramp_tests.jl")

using Documenter
@info "Running doctests..."
doctest(Thermodynamics)
