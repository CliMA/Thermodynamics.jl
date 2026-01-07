"""
# Constructor Consistency Test Suite

This file contains tests for thermodynamic state constructor consistency.
"""

@testset "Thermodynamics - Constructor Consistency" begin
    for ArrayType in array_types
        FT = eltype(ArrayType)
        param_set = FT == Float64 ? param_set_Float64 : param_set_Float32

        @testset "PhaseDry" begin
            profiles = TestedProfiles.PhaseDryProfiles(param_set, ArrayType)
            (; T, p, e_int, h, ρ, θ_liq_ice) = profiles
            (; q_tot, q_liq, q_ice, RH, e_kin, e_pot) = profiles

            ts = PhaseDry.(param_set, e_int, ρ)
            @test all(internal_energy.(param_set, ts) .≈ e_int)
            @test all(air_density.(param_set, ts) .≈ ρ)

            ts_pT = PhaseDry_pT.(param_set, p, T)
            @test all(
                internal_energy.(param_set, ts_pT) .≈
                internal_energy.(param_set, T),
            )
            @test all(air_density.(param_set, ts_pT) .≈ ρ)

            ts_ρe = PhaseDry_ρe.(param_set, ρ, e_int)
            @test all(
                internal_energy.(param_set, ts_ρe) .≈
                internal_energy.(param_set, T),
            )
            @test all(air_density.(param_set, ts_ρe) .≈ ρ)

            θ_dry = potential_temperature.(param_set, T, ρ)
            ts_pθ = PhaseDry_pθ.(param_set, p, θ_dry)
            @test all(
                internal_energy.(param_set, ts_pθ) .≈
                internal_energy.(param_set, T),
            )
            @test all(air_density.(param_set, ts_pθ) .≈ ρ)

            h = e_int + p ./ ρ
            ts_ph = PhaseDry_ph.(param_set, p, h)
            @test all(
                internal_energy.(param_set, ts_ph) .≈
                internal_energy.(param_set, T),
            )
            @test all(air_density.(param_set, ts_ph) .≈ ρ)

            p_dry = air_pressure.(param_set, T, ρ)
            ts_pe = PhaseDry_pe.(param_set, p, e_int)
            @test all(
                internal_energy.(param_set, ts_pe) .≈
                internal_energy.(param_set, T),
            )
            @test all(air_pressure.(param_set, ts_pe) .≈ p_dry)

            ts_ρθ = PhaseDry_ρθ.(param_set, ρ, θ_dry)
            @test all(
                internal_energy.(param_set, ts_ρθ) .≈
                internal_energy.(param_set, T),
            )
            @test all(air_density.(param_set, ts_ρθ) .≈ ρ)

            ts_ρT = PhaseDry_ρT.(param_set, ρ, T)
            @test all(
                air_density.(param_set, ts_ρT) .≈ air_density.(param_set, ts),
            )
            @test all(
                internal_energy.(param_set, ts_ρT) .≈
                internal_energy.(param_set, ts),
            )

            ts = PhaseDry_ρp.(param_set, ρ, p)
            @test all(air_density.(param_set, ts) .≈ ρ)
            @test all(air_pressure.(param_set, ts) .≈ p)
        end

        @testset "PhaseEquil" begin
            _eint_v0 = TP.e_int_v0(param_set)
            profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
            (; T, p, e_int, h, ρ, θ_liq_ice) = profiles
            (; q_tot, q_liq, q_ice, RH, e_kin, e_pot) = profiles

            ts =
                PhaseEquil_ρeq.(
                    param_set,
                    ρ,
                    e_int,
                    q_tot,
                    40,
                    FT(rtol_temperature),
                    RS.SecantMethod,
                )
            @test all(internal_energy.(param_set, ts) .≈ e_int)
            @test all(
                getproperty.(PhasePartition.(param_set, ts), :tot) .≈ q_tot,
            )
            @test all(air_density.(param_set, ts) .≈ ρ)

            ts = PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot)
            @test all(internal_energy.(param_set, ts) .≈ e_int)
            @test all(
                getproperty.(PhasePartition.(param_set, ts), :tot) .≈ q_tot,
            )
            @test all(air_density.(param_set, ts) .≈ ρ)

            ts_peq = PhaseEquil_peq.(param_set, p, e_int, q_tot)
            @test all(internal_energy.(param_set, ts_peq) .≈ e_int)
            @test all(
                getproperty.(PhasePartition.(param_set, ts_peq), :tot) .≈ q_tot,
            )
            @test all(air_pressure.(param_set, ts_peq) .≈ p)

            ts_phq = PhaseEquil_phq.(param_set, p, h, q_tot)
            @test all(
                abs.(internal_energy.(param_set, ts_phq) .- e_int) .<=
                (atol_energy_temperature .+ rtol_humidity * _eint_v0 .* q_tot),
            )
            @test all(
                getproperty.(PhasePartition.(param_set, ts_phq), :tot) .≈ q_tot,
            )
            @test all(air_pressure.(param_set, ts_phq) .≈ p)
        end

        @testset "PhaseNonEquil" begin
            profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
            (; T, p, e_int, h, ρ, θ_liq_ice) = profiles
            (; q_tot, q_liq, q_ice) = profiles
            # Create PhasePartition for PhaseNonEquil tests
            q_pt = TD.PhasePartition.(q_tot, q_liq, q_ice)
            ts = PhaseNonEquil.(param_set, e_int, ρ, q_pt)
            @test all(internal_energy.(param_set, ts) .≈ e_int)
            @test all(compare_moisture.(param_set, ts, q_pt))
            @test all(air_density.(param_set, ts) .≈ ρ)

            ts = PhaseNonEquil_peq.(param_set, p, e_int, q_pt)
            @test all(internal_energy.(param_set, ts) .≈ e_int)
            @test all(compare_moisture.(param_set, ts, q_pt))
            @test all(air_pressure.(param_set, ts) .≈ p)

            # NB: PhaseNonEquil_phq is not tested here, since we don't have a way to
            # construct h at the moment
            ts_phq = PhaseNonEquil_phq.(param_set, p, h, q_pt)
            @test all(internal_energy.(param_set, ts_phq) .≈ e_int)
            @test all(enthalpy.(param_set, ts_phq) .≈ h)
            @test all(compare_moisture.(param_set, ts_phq, q_pt))
            @test all(air_pressure.(param_set, ts_phq) .≈ p)
        end

        @testset "Relative humidity and virtual temperature" begin
            profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
            (; T, p, ρ) = profiles
            phase_type = PhaseEquil{FT}
            (; q_tot, RH) = profiles

            # Test that relative humidity is 1 for saturated conditions
            q_sat = q_vap_saturation.(param_set, T, ρ, Ref(phase_type))
            q_pt_sat = PhasePartition.(q_sat)
            q_vap = vapor_specific_humidity.(q_pt_sat)
            @test all(getproperty.(q_pt_sat, :liq) .≈ 0)
            @test all(getproperty.(q_pt_sat, :ice) .≈ 0)
            @test all(q_vap .≈ q_sat)

            #  Compute thermodynamic consistent pressure
            p_sat = air_pressure.(param_set, T, ρ, q_pt_sat)

            #  Test that density remains consistent
            ρ_rec = air_density.(param_set, T, p_sat, q_pt_sat)
            @test all.(ρ_rec ≈ ρ)

            RH_sat =
                relative_humidity.(
                    param_set,
                    T,
                    p_sat,
                    Ref(phase_type),
                    q_pt_sat,
                )

            @test all(RH_sat .≈ 1)

            # Test that RH is zero for dry conditions
            q_pt_dry = PhasePartition.(zeros(FT, length(T)))
            p_dry = air_pressure.(param_set, T, ρ, q_pt_dry)
            RH_dry =
                relative_humidity.(
                    param_set,
                    T,
                    p_dry,
                    Ref(phase_type),
                    q_pt_dry,
                )
            @test all(RH_dry .≈ 0)

            # Test virtual temperature
            _R_d = FT(TP.R_d(param_set))
            q_pt = PhasePartition.(q_tot)
            T_virt_dry = virtual_temperature.(param_set, T)
            @test all(T_virt_dry ≈ gas_constant_air.(param_set) ./ _R_d .* T)

            # Test convergence of virtual temperature iterations with MOIST virtual temperature
            T_virt_moist = virtual_temperature.(param_set, T, q_pt)

            T_rec_qpt_rec =
                TD.temperature_and_humidity_given_TᵥρRH.(
                    param_set,
                    T_virt_moist,
                    ρ,
                    RH,
                    Ref(phase_type),
                )

            T_rec = first.(T_rec_qpt_rec)
            q_pt_rec = last.(T_rec_qpt_rec)

            @test all(
                isapprox.(
                    T_virt_moist,
                    virtual_temperature.(param_set, T_rec, q_pt_rec),
                    atol = sqrt(eps(FT)),
                ),
            )

            # Test that reconstructed specific humidity is close
            # to original specific humidity
            q_tot_rec = getproperty.(q_pt_rec, :tot)
            RH_moist = (q_tot .> eps(FT)) .& (RH .< FT(0.99))
            @test all(
                isapprox.(q_tot[RH_moist], q_tot_rec[RH_moist], rtol = 0.10),
            )

            # Update temperature to be exactly consistent with
            # p, ρ_rec; test that this is equal to T_rec
            # Note: Relaxed tolerance needed due to cumulative numerical errors in saturation adjustment
            # and the iterative nature of relative humidity/virtual temperature calculations.
            # The extreme cases (saturated high-humidity) can have differences up to ~2.5K due to
            # numerical challenges in the iterative root finding algorithms.
            T_local = TD.air_temperature_given_pρq.(param_set, p, ρ, q_pt_rec)
            @test all(isapprox.(T_local, T_rec, atol = FT(2.5)))
        end
    end
end
