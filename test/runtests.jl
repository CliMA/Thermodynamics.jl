if !haskey(ENV, "BUILDKITE")
    import Pkg
    Pkg.develop(Pkg.PackageSpec(; path = dirname(@__DIR__)))
end

# read parameters needed for tests
import CLIMAParameters
full_parameter_set = CLIMAParameters.create_parameter_dict(dict_type = "alias")


#include("parameter_tests.jl")
include("TemperatureProfiles.jl")
include("relations.jl")
