if !haskey(ENV, "BUILDKITE")
    import Pkg
    Pkg.develop(Pkg.PackageSpec(; path = dirname(@__DIR__)))
end

# read parameters needed for tests
import CLIMAParameters
full_parameter_set =
    CLIMAParameters.create_parameter_struct(dict_type = "alias")


include("TemperatureProfiles.jl")
include("relations.jl")
