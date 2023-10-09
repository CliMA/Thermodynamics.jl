using Thermodynamics, Documenter
using DocumenterCitations

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"

bib = CitationBibliography(joinpath(@__DIR__, "bibliography.bib"))

pages = Any[
    "Home" => "index.md",
    "Installation" => "Installation.md",
    "API" => "API.md",
    "How-to-guide" => "HowToGuide.md",
    "Tested profiles" => "TestedProfiles.md",
    "Temperature profiles" => "TemperatureProfiles.md",
    "Developer docs" => "DevDocs.md",
    "Clausius Clapeyron relation" => "Clausius_Clapeyron.md",
    "Thermodynamics overview" => "Formulation.md",
    "References" => "References.md",
]

mathengine = MathJax(
    Dict(
        :TeX => Dict(
            :equationNumbers => Dict(:autoNumber => "AMS"),
            :Macros => Dict(),
        ),
    ),
)

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1,
)

makedocs(;
    plugins = [bib],
    sitename = "Thermodynamics.jl",
    format = format,
    checkdocs = :exports,
    clean = true,
    doctest = true,
    modules = [Thermodynamics],
    pages = pages,
)

deploydocs(
    repo = "github.com/CliMA/Thermodynamics.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)
