
rm(joinpath(@__DIR__, "Manifest.toml"), force = true)       # Remove local Manifest.toml
rm(joinpath(@__DIR__, "..", "Manifest.toml"), force = true) # Remove local Manifest.toml

push!(LOAD_PATH, joinpath(@__DIR__, "..", "env", "Plots")) # add Plots env

# Avoiding having to add Thermodynamics deps to docs/ environment:
push!(LOAD_PATH, joinpath(@__DIR__, ".."))                 # add Thermodynamics env

using Pkg
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))

using Thermodynamics, Documenter

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
    sitename = "Thermodynamics.jl",
    format = format,
    clean = true,
    modules = [Documenter, Thermodynamics],
    pages = pages,
)

deploydocs(
    repo = "github.com/CliMA/Thermodynamics.jl.git",
    target = "build",
    push_preview = true,
)
