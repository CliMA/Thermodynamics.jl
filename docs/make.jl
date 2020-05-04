
rm(joinpath(@__DIR__, "Manifest.toml"), force = true) # Remove local Manifest

push!(LOAD_PATH, joinpath(@__DIR__, "..", "env", "Plots")) # add Plots env
push!(LOAD_PATH, joinpath(@__DIR__, ".."))                 # add MoistThermodynamics env

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate(; verbose = true)

using MoistThermodynamics, Documenter

pages = Any[
    "Home" => "index.md",
    "Installation" => "Installation.md",
    "API" => "API.md",
    "How-to-guide" => "HowToGuide.md",
    "Tested profiles" => "TestedProfiles.md",
]

mathengine = MathJax(Dict(
    :TeX => Dict(
        :equationNumbers => Dict(:autoNumber => "AMS"),
        :Macros => Dict(),
    ),
))

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1,
)

makedocs(
    sitename = "MoistThermodynamics.jl",
    format = format,
    clean = true,
    modules = [Documenter, MoistThermodynamics],
    pages = pages,
)

deploydocs(
    repo = "github.com/climate-machine/MoistThermodynamics.jl.git",
    target = "build",
    push_preview = true,
)
