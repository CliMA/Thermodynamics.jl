"""
# Miscellaneous Test Suite

This file contains various miscellaneous tests including ProfileSet Iterator, Base.zero, T_guess, and data collection.
"""

using Test
using Thermodynamics
import Thermodynamics as TD
import Thermodynamics.Parameters as TP
using Thermodynamics.TestedProfiles

# Saturation adjustment tolerance (relative change of temperature between consecutive iterations)
rtol_temperature = 1e-4

@testset "Thermodynamics - ProfileSet Iterator" begin
    ArrayType = Array{Float64}
    FT = eltype(ArrayType)
    param_set = TP.ThermodynamicsParameters(FT)
    profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
    (; T, q_pt, z, phase_type) = profiles
    @test all(z .≈ (nt.z for nt in profiles))
    @test all(T .≈ (nt.T for nt in profiles))
    @test all(getproperty.(q_pt, :tot) .≈ (nt.q_pt.tot for nt in profiles))
    @test all(phase_type .== (nt.phase_type for nt in profiles))
end

@testset "Base.zero" begin
    FT = Float32
    @test zero(PhasePartition{FT}).tot == 0
    @test zero(PhaseDry{FT}).ρ == 0
    @test zero(PhaseEquil{FT}).ρ == 0
    @test zero(PhaseNonEquil{FT}).ρ == 0
end

@testset "Thermodynamics - Test T_guess" begin
    ArrayType = Array{Float64}
    FT = eltype(ArrayType)
    param_set = TP.ThermodynamicsParameters(FT)
    profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
    (; p, ρ, e_int, h, θ_liq_ice, q_tot, T, phase_type) = profiles
    T_guess = T .+ (FT(0.2) .* randn(FT, length(T)))
    args = (q_tot, 40, FT(rtol_temperature))
    ts =
        PhaseEquil_ρeq.(param_set, ρ, e_int, args..., RS.NewtonsMethod, T_guess)
    ts = PhaseEquil_ρθq.(param_set, ρ, θ_liq_ice, args..., T_guess)
    ts = PhaseEquil_peq.(param_set, p, e_int, args..., RS.SecantMethod, T_guess)
    ts = PhaseEquil_phq.(param_set, p, h, args..., RS.SecantMethod, T_guess)
    ts =
        PhaseEquil_ρpq.(
            param_set,
            ρ,
            p,
            q_tot,
            true,
            40,
            FT(rtol_temperature),
            RS.NewtonsMethodAD,
            T_guess,
        )
    ts =
        PhaseEquil_pθq.(
            param_set,
            p,
            θ_liq_ice,
            args...,
            RS.SecantMethod,
            T_guess,
        )
end

TD.solution_type() = RS.VerboseSolution()
@testset "Test Data Collection" begin
    ArrayType = Array{Float64}
    FT = eltype(ArrayType)
    param_set = TP.ThermodynamicsParameters(FT)
    profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
    (; ρ, e_int, q_tot) = profiles
    ts = PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot)
    data = TD.DataCollection.get_data()
    TD.DataCollection.print_summary(data)
    TD.DataCollection.reset_stats()
end
TD.solution_type() = RS.CompactSolution() 