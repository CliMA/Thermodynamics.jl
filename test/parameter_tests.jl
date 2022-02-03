

using Test


### read parameters needed for tests
using TOML
planet_parse = TOML.parsefile(joinpath(@__DIR__, "planet_parameters.toml"))
param_dict = Dict{Symbol, Any}()
for (key, val) in planet_parse
    # In the future - we will use the full names,
    # param_dict[Symbol(key)] = val["value"]

    # for now we use the aliases
    param_dict[Symbol(val["alias"])] = val["value"]
end
full_parameter_set = (; param_dict...)
universal_constant_aliases = [
    :gas_constant,
    :light_speed,
    :h_Planck,
    :k_Boltzmann,
    :Stefan,
    :astro_unit,
    :avogad,
]



#get CLIMAParameters to check consistency
using CLIMAParameters
using CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set_cpp = EarthParameterSet()

using Thermodynamics


@testset "parameter read tests" begin

    #check parameter_file agrees with current CP defaults
    for (k, v) in zip(keys(full_parameter_set), full_parameter_set)
        if ~(k in universal_constant_aliases)
            @test (v ≈ @eval $k(param_set_cpp))
        else
            @test (v ≈ @eval $k())
        end
    end

    #check new local implementation agrees with CP
    param_set_by_alias =
        Thermodynamics.ThermodynamicsParameters(full_parameter_set)
    for fn in fieldnames(typeof(param_set_by_alias))
        v = getfield(param_set_by_alias, fn)
        if ~(fn in universal_constant_aliases)
            @test (v ≈ @eval $fn(param_set_cpp))
        else
            @test (v ≈ @eval $fn())
        end
    end
end
