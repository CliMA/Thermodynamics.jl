rm(joinpath(@__DIR__, "Manifest.toml"), force = true)       # Remove local Manifest.toml
rm(joinpath(@__DIR__, "..", "Manifest.toml"), force = true) # Remove local Manifest.toml

# Avoiding having to add deps to docs/ environment:
env_viz = joinpath(@__DIR__, "..", "env", "viz")
env_doc = @__DIR__

using Pkg
push!(LOAD_PATH, env_viz); Pkg.activate(env_viz); Pkg.instantiate(; verbose=true)
push!(LOAD_PATH, env_doc); Pkg.activate(env_doc); Pkg.instantiate(; verbose=true)

cd(joinpath(@__DIR__, "..")) do
    Pkg.develop(PackageSpec(path="."))
    Pkg.activate(pwd())
    Pkg.instantiate(; verbose=true)
end

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
