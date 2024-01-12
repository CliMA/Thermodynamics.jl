include("aqua.jl")
include("TemperatureProfiles.jl")
include("relations.jl")

rm(joinpath(@__DIR__, "logfilepath_Float32.toml"); force = true)
rm(joinpath(@__DIR__, "logfilepath_Float64.toml"); force = true)
