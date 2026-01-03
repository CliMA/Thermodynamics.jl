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

# Run temperature profile tests for physical consistency
include("TemperatureProfiles.jl")
include("TestedProfiles.jl")

# Include common parameters, functions, and imports for all tests
include("test_common.jl")

# Run core thermodynamic function tests (split into functional categories)
include("dry_adiabatic_processes.jl")
include("correctness.jl")
include("default_behavior_accuracy.jl")
include("exceptions.jl")
include("constructor_consistency.jl")
include("type_stability.jl")
include("dry_limit.jl")
include("miscellaneous.jl")
include("data_tests.jl")

# Run aqua tests for code quality and performance
include("aqua.jl")

# Clean up any temporary files created during testing
rm(joinpath(@__DIR__, "logfilepath_Float32.toml"); force = true)
rm(joinpath(@__DIR__, "logfilepath_Float64.toml"); force = true)

