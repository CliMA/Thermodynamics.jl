using Test
using Thermodynamics
using Aqua

@testset "Aqua tests (performance)" begin
    # Test for unbound type parameters in method signatures
    # This prevents issues like https://github.com/JuliaLang/julia/issues/29393
    # where unbound type parameters can cause compilation problems
    ua = Aqua.detect_unbound_args_recursively(Thermodynamics)
    @test length(ua) == 0

    # Test for method ambiguities across dependencies
    # Method ambiguities can cause compilation delays and unexpected behavior
    # See: https://github.com/SciML/OrdinaryDiffEq.jl/issues/1750
    ambs = Aqua.detect_ambiguities(Thermodynamics; recursive = true)

    # Filter ambiguities to only include those involving Thermodynamics methods
    # This helps focus on ambiguities that are actually relevant to our package
    pkg_match(pkgname, pkdir::Nothing) = false
    pkg_match(pkgname, pkdir::AbstractString) = occursin(pkgname, pkdir)
    filter!(x -> pkg_match("Thermodynamics", pkgdir(last(x).module)), ambs)

    # Display any remaining ambiguities for debugging
    # The goal is to have zero ambiguities to minimize compilation latency
    for method_ambiguity in ambs
        @show method_ambiguity
    end
    @test length(ambs) ≤ 0
end

@testset "Aqua tests - remaining" begin
    # Run comprehensive Aqua tests for code quality
    # This includes tests for:
    # - Project structure and dependencies
    # - Type stability
    # - Package compatibility
    # - Documentation coverage
    # Note: We exclude ambiguities and unbound_args since they're tested above
    Aqua.test_all(Thermodynamics; ambiguities = false, unbound_args = false)
end
