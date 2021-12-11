if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end

include("TemperatureProfiles.jl")
include("relations.jl")
include("runtests_gpu.jl")
