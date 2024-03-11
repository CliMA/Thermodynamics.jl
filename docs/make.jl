using Thermodynamics
using Documenter, DocumenterCitations, Literate, Printf

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"

bib = CitationBibliography(joinpath(@__DIR__, "bibliography.bib"))

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR = joinpath(@__DIR__, "src/literated")

example_scripts = ["density_from_temperature_pressure_humidity.jl"]

if isdir(OUTPUT_DIR)
    rm(OUTPUT_DIR)
end

mkdir(OUTPUT_DIR)

cp(
    joinpath(EXAMPLES_DIR, "JRA55_atmospheric_state_Jan_1_1991.jld2"),
    joinpath(OUTPUT_DIR, "JRA55_atmospheric_state_Jan_1_1991.jld2"),
)

for example in example_scripts
    example_filepath = joinpath(EXAMPLES_DIR, example)

    withenv("JULIA_DEBUG" => "Literate") do
        start_time = time_ns()
        Literate.markdown(
            example_filepath,
            OUTPUT_DIR;
            flavor = Literate.DocumenterFlavor(),
            execute = true,
        )
        elapsed = 1e-9 * (time_ns() - start_time)
        @info @sprintf("%s example took %s seconds to build.", example, elapsed)
    end
end

example_pages = [
    "Density from temperature, pressure, and humidity" => "literated/density_from_temperature_pressure_humidity.md",
]

pages = Any[
    "Home" => "index.md",
    "Installation" => "Installation.md",
    "API" => "API.md",
    "How-to-guide" => "HowToGuide.md",
    "Examples" => example_pages,
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

const MiB = 2^20

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1,
    size_threshold = 2 * MiB,
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
