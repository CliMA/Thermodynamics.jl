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

# Functional tests
include("test_common.jl")

# Run core thermodynamic function tests (split into functional categories)
include("dry_adiabatic_processes.jl")
include("correctness.jl")
include("exceptions.jl")
include("saturation_adjustment.jl")
include("convergence_saturation_adjustment.jl")
include("type_stability.jl")

const _include_deprecated =
    lowercase(get(ENV, "THERMODYNAMICS_INCLUDE_DEPRECATED", "true")) ∉ ("0", "false", "no")

if _include_deprecated
    # Deprecated tests (PhasePartition and thermodynamic state methods in `src/depr_*`); soon to be removed
    include("depr_test_common.jl")
    include("depr_default_behavior_accuracy.jl")
    include("depr_constructor_consistency.jl")
    include("depr_type_stability.jl")
    include("depr_dry_limit.jl")
    include("depr_miscellaneous.jl")
    include("depr_data_tests.jl")
    include("depr_exceptions.jl")
end

# Run aqua tests for code quality and performance
include("aqua.jl")
