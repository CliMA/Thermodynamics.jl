using Thermodynamics
using Documenter, DocumenterCitations, Literate, Printf

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"

bib = CitationBibliography(joinpath(@__DIR__, "bibliography.bib"))

pages = Any[
    "Home" => "index.md",
    "Mathematical Formulation" => "Formulation.md",
    "API Reference" => "API.md",
    "How-To Guide" => "HowToGuide.md",
    "Temperature Profiles" => "TemperatureProfiles.md",
    "Tested Profiles" => "TestedProfiles.md",
    "Clausius-Clapeyron Validation" => "Clausius_Clapeyron.md",
    "Saturation Adjustment Convergence" => "SaturationAdjustmentConvergence.md",
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

const MiB = 2^20

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1,
    size_threshold = 2MiB,
    assets = ["src/assets/custom.css"],
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

@info "Clean up temporary .jld2 and .nc output created by doctests or literated examples..."

"""
    recursive_find(directory, pattern)

Return list of filepaths within `directory` that contains the `pattern::Regex`.
"""
recursive_find(directory, pattern) =
    mapreduce(vcat, walkdir(directory)) do (root, dirs, files)
        joinpath.(root, filter(contains(pattern), files))
    end

files = []
for pattern in [r"\.jld2", r"\.nc"]
    global files = vcat(files, recursive_find(@__DIR__, pattern))
end

for file in files
    rm(file)
end

deploydocs(
    repo = "github.com/CliMA/Thermodynamics.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)
