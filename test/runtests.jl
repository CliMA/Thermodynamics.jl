if !haskey(ENV, "BUILDKITE")
    import Pkg
    Pkg.develop(Pkg.PackageSpec(; path = dirname(@__DIR__)))
end

# read parameters needed for tests
using TOML
planet_parse = TOML.parsefile(joinpath(@__DIR__, "planet_parameters.toml"))
param_dict = Dict{Symbol, Float64}()
for (key, val) in planet_parse
    # In the future - we will use the full names,
    # param_dict[Symbol(key)] = val["value"]

    # for now we use the aliases
    param_dict[Symbol(val["alias"])] = val["value"]
end
full_parameter_set = (; param_dict...)

include("parameter_tests.jl")
include("TemperatureProfiles.jl")
include("relations.jl")
