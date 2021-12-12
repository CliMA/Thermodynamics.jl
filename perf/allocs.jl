if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
import Thermodynamics
import RootSolvers
import CLIMAParameters
import Coverage
import Plots

# Packages to monitor

mod_dir(x) = dirname(dirname(pathof(x)))
all_dirs_to_monitor = [
    ".",
    mod_dir(Thermodynamics),
    mod_dir(RootSolvers),
    mod_dir(CLIMAParameters),
]

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"

constructors = ["ρeq", "pθq", "pTq"]

allocs = Dict()
for constructor in constructors
    ENV["ALLOCATION_CONSTRUCTOR"] = constructor
    run(`julia --project=test/ --track-allocation=all perf/alloc_per_constructor.jl`)

    allocs[constructor] = Coverage.analyze_malloc(all_dirs_to_monitor)

    # Clean up files
    for d in all_dirs_to_monitor
        all_files = [
            joinpath(root, f)
            for (root, dirs, files) in Base.Filesystem.walkdir(d)
            for f in files
        ]
        all_mem_files = filter(x -> endswith(x, ".mem"), all_files)
        for f in all_mem_files
            rm(f)
        end
    end
end

@info "Post-processing allocations"

function plot_allocs(constructor, allocs_per_case, n_unique_bytes)
    p = Plots.plot()
    @info "Allocations for $constructor"

    function filename_only(fn)
        fn = first(split(fn, ".jl")) * ".jl"
        splitby = "central/scratch/climaci/thermodynamics-ci/"
        if occursin(splitby, fn)
            fn = last(split(fn, splitby))
        end
        splitby = "/thermodynamics-ci/"
        if occursin(splitby, fn)
            fn = joinpath("Thermodynamics.jl", last(split(fn, splitby)))
        end
        splitby = "/depot/cpu/packages/"
        if occursin(splitby, fn)
            fn = last(split(fn, splitby))
        end
        return fn
    end
    function compile_pkg(fn, linenumber)
        c1 = endswith(filename_only(fn), "Thermodynamics.jl")
        c2 = linenumber == 1
        return c1 && c2
    end

    filter!(x -> x.bytes ≠ 0, allocs_per_case)
    filter!(x -> !compile_pkg(x.filename, x.linenumber), allocs_per_case)

    for alloc in allocs_per_case
        println(alloc)
    end
    println("Number of allocating sites: $(length(allocs_per_case))")
    case_bytes = getproperty.(allocs_per_case, :bytes)[end:-1:1]
    case_filename = getproperty.(allocs_per_case, :filename)[end:-1:1]
    case_linenumber = getproperty.(allocs_per_case, :linenumber)[end:-1:1]
    all_bytes = Int[]
    filenames = String[]
    linenumbers = Int[]
    loc_ids = String[]
    for (bytes, filename, linenumber) in
        zip(case_bytes, case_filename, case_linenumber)
        compile_pkg(filename, linenumber) && continue # Skip loading module
        loc_id = "$(filename_only(filename))" * "$linenumber"
        if !(bytes in all_bytes) && !(loc_id in loc_ids)
            push!(all_bytes, bytes)
            push!(filenames, filename)
            push!(linenumbers, linenumber)
            push!(loc_ids, loc_id)
            if length(all_bytes) ≥ n_unique_bytes
                break
            end
        end
    end

    all_bytes = all_bytes ./ 10^3
    if isempty(all_bytes)
        @info "$constructor: 0 allocations! 🎉"
        return nothing
    end
    max_bytes = maximum(all_bytes)
    @info "$constructor: $all_bytes"
    xtick_name(filename, linenumber) = "$filename, line number: $linenumber"
    markershape = (:square, :hexagon, :circle, :star, :utriangle, :dtriangle)
    for (bytes, filename, linenumber) in zip(all_bytes, filenames, linenumbers)
        Plots.plot!(
            [0],
            [bytes];
            seriestype = :scatter,
            label = xtick_name(filename_only(filename), linenumber),
            markershape = markershape[1],
            markersize = 1 + bytes / max_bytes * 10,
        )
        markershape = (markershape[end], markershape[1:(end - 1)]...)
    end
    p1 = Plots.plot!(ylabel = "Allocations (KB)", title = constructor)
    subset_allocs_per_case = collect(Base.Iterators.Reverse(allocs_per_case))
    p2 = Plots.plot(
        1:length(subset_allocs_per_case),
        getproperty.(subset_allocs_per_case, :bytes) ./ 1000;
        xlabel = "i-th allocating line (truncated and sorted)",
        ylabel = "Allocations (KB)",
        markershape = :circle,
    )
    Plots.plot(p1, p2, layout = Plots.grid(2, 1))
    Plots.savefig(joinpath(folder, "allocations_$constructor.png"))
end

folder = "perf/allocations_output"
mkpath(folder)

@info "Allocated bytes for constructor broadcasted across tested profiles"
for constructor in constructors
    plot_allocs(constructor, allocs[constructor], 10)
end
