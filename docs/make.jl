using Thermodynamics, Documenter

pages = Any[
    "Home" => "index.md",
    "Installation" => "Installation.md",
    "API" => "API.md",
    "How-to-guide" => "HowToGuide.md",
    "Tested profiles" => "TestedProfiles.md",
    "Temperature profiles" => "TemperatureProfiles.md",
    "Developer docs" => "DevDocs.md",
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
    strict = true,
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
