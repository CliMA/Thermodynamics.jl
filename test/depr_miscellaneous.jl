"""
# Miscellaneous Test Suite

This file contains various miscellaneous tests including ProfileSet Iterator, Base.zero, T_guess, and data collection.
"""

@testset "Thermodynamics - Miscellaneous" begin
    @testset "ProfileSet Iterator" begin
        ArrayType = Array{Float64}
        FT = eltype(ArrayType)
        param_set = FT == Float64 ? param_set_Float64 : param_set_Float32
        profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
        (; T, z) = profiles
        @test all(z .≈ (nt.z for nt in profiles))
        @test all(T .≈ (nt.T for nt in profiles))
    end

    @testset "Base.zero" begin
        FT = Float32
        @test zero(PhasePartition{FT}).tot == 0
        @test zero(PhaseDry{FT}).ρ == 0
        @test zero(PhaseEquil{FT}).ρ == 0
        @test zero(PhaseNonEquil{FT}).ρ == 0
    end

    @testset "T_guess" begin
        ArrayType = Array{Float64}
        FT = eltype(ArrayType)
        param_set = FT == Float64 ? param_set_Float64 : param_set_Float32
        profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
        (; p, ρ, e_int, h, θ_li, q_tot, T) = profiles
        T_guess = T .+ (FT(0.2) .* randn(FT, length(T)))
        args = (q_tot, 40, FT(rtol_temperature))
        ts =
            PhaseEquil_ρeq.(
                param_set,
                ρ,
                e_int,
                args...,
                RS.NewtonsMethod,
                T_guess,
            )
        ts =
            PhaseEquil_ρeq.(
                param_set,
                ρ,
                e_int,
                args...,
                RS.NewtonsMethodAD,
                T_guess,
            )
        ts = PhaseEquil_ρθq.(param_set, ρ, θ_li, args..., T_guess)
        ts =
            PhaseEquil_peq.(
                param_set,
                p,
                e_int,
                args...,
                RS.SecantMethod,
                T_guess,
            )
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
                θ_li,
                args...,
                RS.SecantMethod,
                T_guess,
            )
    end

    @testset "Data Collection" begin
        TD.solution_type() = RS.VerboseSolution()
        ArrayType = Array{Float64}
        FT = eltype(ArrayType)
        param_set = FT == Float64 ? param_set_Float64 : param_set_Float32
        profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
        (; ρ, e_int, q_tot) = profiles
        ts = PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot)
        data = TD.DataCollection.get_data()
        TD.DataCollection.print_summary(data)
        TD.DataCollection.reset_stats()
        TD.solution_type() = RS.CompactSolution()
    end
end
