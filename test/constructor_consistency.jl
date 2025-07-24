"""
# Constructor Consistency Test Suite

This file contains tests for thermodynamic state constructor consistency.
"""

@testset "Thermodynamics - Constructor Consistency" begin

    # Make sure `ThermodynamicState` arguments are returned unchanged

    for ArrayType in array_types
        FT = eltype(ArrayType)
        param_set = TP.ThermodynamicsParameters(FT)
        _eint_v0 = TP.e_int_v0(param_set)

        profiles = TestedProfiles.PhaseDryProfiles(param_set, ArrayType)
        (; T, p, e_int, h, ρ, θ_liq_ice, phase_type) = profiles
        (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = profiles

        # PhaseDry
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

        θ_dry = dry_pottemp.(param_set, T, ρ)
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
        @test all(air_density.(param_set, ts_ρT) .≈ air_density.(param_set, ts))
        @test all(
            internal_energy.(param_set, ts_ρT) .≈
            internal_energy.(param_set, ts),
        )


        ts = PhaseDry_ρp.(param_set, ρ, p)
        @test all(air_density.(param_set, ts) .≈ ρ)
        @test all(air_pressure.(param_set, ts) .≈ p)
        e_tot_proposed =
            TD.total_energy_given_ρp.(param_set, ρ, p, e_kin, e_pot)
        @test all(total_energy.(param_set, ts, e_kin, e_pot) .≈ e_tot_proposed)


        profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
        (; T, p, e_int, h, ρ, θ_liq_ice, phase_type) = profiles
        (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = profiles

        # PhaseEquil
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
        @test all(getproperty.(PhasePartition.(param_set, ts), :tot) .≈ q_tot)
        @test all(air_density.(param_set, ts) .≈ ρ)

        ts = PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot)
        @test all(internal_energy.(param_set, ts) .≈ e_int)
        @test all(getproperty.(PhasePartition.(param_set, ts), :tot) .≈ q_tot)
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

        # ts_phq is challenged by extremely humid conditions (e.g., q ≈ 0.16). 
        # If it fails, this lists where
        errors = abs.(internal_energy.(param_set, ts_phq) .- e_int)
        tolerances =
            atol_energy_temperature .+ rtol_humidity * _eint_v0 .* q_tot
        failing = findall(errors .> tolerances)
        # if !isempty(failing)
        #     @info "PhaseEquil_phq energy test failures" n=length(failing)
        #     @show p[failing]
        #     @show h[failing]
        #     @show q_tot[failing]
        #     @show errors[failing]
        #     @show tolerances[failing]
        # end

        @test all(
            getproperty.(PhasePartition.(param_set, ts_phq), :tot) .≈ q_tot,
        )
        @test all(air_pressure.(param_set, ts_phq) .≈ p)

        ts_pθq = PhaseEquil_pθq.(param_set, p, θ_liq_ice, q_tot)
        @test all(air_pressure.(param_set, ts_pθq) .≈ p)

        @test all(
            isapprox.(
                liquid_ice_pottemp.(param_set, ts_pθq),
                θ_liq_ice,
                atol = 2.25 * atol_temperature,  # Larger tol because of condensate effects
            ),
        )
        @test all(
            getproperty.(PhasePartition.(param_set, ts_pθq), :tot) .≈ q_tot,
        )

        ts = PhaseEquil_ρpq.(param_set, ρ, p, q_tot, true)
        @test all(air_density.(param_set, ts) .≈ ρ)
        @test all(air_pressure.(param_set, ts) .≈ p)
        @test all(getproperty.(PhasePartition.(param_set, ts), :tot) .≈ q_tot)

        # Test against total_energy_given_ρp when not iterating
        ts = PhaseEquil_ρpq.(param_set, ρ, p, q_tot, false)
        e_tot_proposed =
            TD.total_energy_given_ρp.(
                param_set,
                ρ,
                p,
                e_kin,
                e_pot,
                PhasePartition.(q_tot),
            )
        @test all(total_energy.(param_set, ts, e_kin, e_pot) .≈ e_tot_proposed)

        # PhaseNonEquil
        ts = PhaseNonEquil.(param_set, e_int, ρ, q_pt)
        @test all(internal_energy.(param_set, ts) .≈ e_int)
        @test all(compare_moisture.(param_set, ts, q_pt))
        @test all(air_density.(param_set, ts) .≈ ρ)

        ts = PhaseNonEquil_peq.(param_set, p, e_int, q_pt)
        @test all(internal_energy.(param_set, ts) .≈ e_int)
        @test all(compare_moisture.(param_set, ts, q_pt))
        @test all(air_pressure.(param_set, ts) .≈ p)

        ts_phq = PhaseNonEquil_phq.(param_set, p, h, q_pt)
        @test all(internal_energy.(param_set, ts_phq) .≈ e_int)
        @test all(specific_enthalpy.(param_set, ts_phq) .≈ h)
        @test all(compare_moisture.(param_set, ts_phq, q_pt))
        @test all(air_pressure.(param_set, ts_phq) .≈ p)

        # TD.air_temperature_given_pθq-liquid_ice_pottemp inverse
        θ_liq_ice_ =
            TD.liquid_ice_pottemp_given_pressure.(param_set, T, p, q_pt)
        @test all(
            TD.air_temperature_given_pθq.(param_set, p, θ_liq_ice_, q_pt) .≈ T,
        )

        # liquid_ice_pottemp-TD.air_temperature_given_pθq inverse
        T = TD.air_temperature_given_pθq.(param_set, p, θ_liq_ice, q_pt)
        @test all(
            TD.liquid_ice_pottemp_given_pressure.(param_set, T, p, q_pt) .≈
            θ_liq_ice,
        )

        # Accurate but expensive `PhaseNonEquil_ρθq` constructor (Non-linear temperature from θ_liq_ice)
        T_non_linear =
            TD.air_temperature_given_ρθq_nonlinear.(
                param_set,
                ρ,
                θ_liq_ice,
                20,
                RS.ResidualTolerance(FT(5e-5)),
                q_pt,
            )
        T_expansion =
            TD.air_temperature_given_ρθq.(param_set, ρ, θ_liq_ice, q_pt)
        @test all(isapprox.(T_non_linear, T_expansion, atol = atol_temperature))
        e_int_ = internal_energy.(param_set, T_non_linear, q_pt)
        ts = PhaseNonEquil.(param_set, e_int_, ρ, q_pt)
        @test all(T_non_linear .≈ air_temperature.(param_set, ts))
        @test all(
            isapprox(
                θ_liq_ice,
                liquid_ice_pottemp.(param_set, ts),
                atol = atol_temperature,
            ),
        )

        # PhaseEquil_ρθq
        ts =
            PhaseEquil_ρθq.(
                param_set,
                ρ,
                θ_liq_ice,
                q_tot,
                45,
                FT(rtol_temperature),
            )
        @test all(
            isapprox.(
                liquid_ice_pottemp.(param_set, ts),
                θ_liq_ice,
                atol = atol_temperature,
            ),
        )
        @test all(
            isapprox.(air_density.(param_set, ts), ρ, rtol = rtol_density),
        )
        @test all(getproperty.(PhasePartition.(param_set, ts), :tot) .≈ q_tot)

        # The PhaseEquil_pθq constructor
        # passes the consistency test within sufficient physical precision.
        # However, it fails to satisfy the consistency test within machine
        # precision for the input pressure.

        # PhaseEquil_pθq
        ts =
            PhaseEquil_pθq.(
                param_set,
                p,
                θ_liq_ice,
                q_tot,
                35,
                FT(rtol_temperature),
            )
        @test all(
            isapprox.(
                liquid_ice_pottemp.(param_set, ts),
                θ_liq_ice,
                atol = 2.25 * atol_temperature,
            ),
        )
        @test all(compare_moisture.(param_set, ts, q_pt))
        @test all(air_pressure.(param_set, ts) .≈ p)

        # PhaseNonEquil_pθq
        ts = PhaseNonEquil_pθq.(param_set, p, θ_liq_ice, q_pt)
        @test all(liquid_ice_pottemp.(param_set, ts) .≈ θ_liq_ice)
        @test all(air_pressure.(param_set, ts) .≈ p)
        @test all(compare_moisture.(param_set, ts, q_pt))

        ts = PhaseNonEquil_ρpq.(param_set, ρ, p, q_pt)
        @test all(air_density.(param_set, ts) .≈ ρ)
        @test all(air_pressure.(param_set, ts) .≈ p)
        @test all(
            getproperty.(PhasePartition.(param_set, ts), :tot) .≈
            getproperty.(q_pt, :tot),
        )
        @test all(
            getproperty.(PhasePartition.(param_set, ts), :liq) .≈
            getproperty.(q_pt, :liq),
        )
        @test all(
            getproperty.(PhasePartition.(param_set, ts), :ice) .≈
            getproperty.(q_pt, :ice),
        )
        e_tot_proposed =
            TD.total_energy_given_ρp.(param_set, ρ, p, e_kin, e_pot, q_pt)
        @test all(total_energy.(param_set, ts, e_kin, e_pot) .≈ e_tot_proposed)

        # PhaseNonEquil_ρθq
        ts =
            PhaseNonEquil_ρθq.(
                param_set,
                ρ,
                θ_liq_ice,
                q_pt,
                5,
                FT(rtol_temperature),
            )
        @test all(
            isapprox.(
                θ_liq_ice,
                liquid_ice_pottemp.(param_set, ts),
                atol = atol_temperature,
            ),
        )
        @test all(air_density.(param_set, ts) .≈ ρ)
        @test all(compare_moisture.(param_set, ts, q_pt))

        profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
        (; T, p, e_int, ρ, θ_liq_ice, phase_type) = profiles
        (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = profiles

        # Test that relative humidity is 1 for saturated conditions
        q_sat = q_vap_saturation.(param_set, T, ρ, Ref(phase_type))
        q_pt_sat = PhasePartition.(q_sat)
        q_vap = vapor_specific_humidity.(q_pt_sat)
        @test all(getproperty.(q_pt_sat, :liq) .≈ 0)
        @test all(getproperty.(q_pt_sat, :ice) .≈ 0)
        @test all(q_vap .≈ q_sat)

        # Compute thermodynamic consistent pressure
        p_sat = air_pressure.(param_set, T, ρ, q_pt_sat)

        # Test that density remains consistent
        ρ_rec = air_density.(param_set, T, p_sat, q_pt_sat)
        @test all.(ρ_rec ≈ ρ)

        RH_sat =
            relative_humidity.(param_set, T, p_sat, Ref(phase_type), q_pt_sat)

        @test all(RH_sat .≈ 1)

        # Test that RH is zero for dry conditions
        q_pt_dry = PhasePartition.(zeros(FT, length(T)))
        p_dry = air_pressure.(param_set, T, ρ, q_pt_dry)
        RH_dry =
            relative_humidity.(param_set, T, p_dry, Ref(phase_type), q_pt_dry)
        @test all(RH_dry .≈ 0)

        # Test virtual temperature 
        _R_d = FT(TP.R_d(param_set))
        T_virt = virtual_temperature.(param_set, T, q_pt)
        @test all(T_virt ≈ gas_constant_air.(param_set, q_pt) ./ _R_d .* T)

        T_rec_qpt_rec =
            TD.temperature_and_humidity_given_TᵥρRH.(
                param_set,
                T_virt,
                ρ,
                RH,
                Ref(phase_type),
            )

        T_rec = first.(T_rec_qpt_rec)
        q_pt_rec = last.(T_rec_qpt_rec)

        # Test convergence of virtual temperature iterations
        @test all(
            isapprox.(
                T_virt,
                virtual_temperature.(param_set, T_rec, q_pt_rec),
                atol = sqrt(eps(FT)),
            ),
        )

        # Test that reconstructed specific humidity is close
        # to original specific humidity
        q_tot_rec = getproperty.(q_pt_rec, :tot)
        RH_moist = q_tot .> eps(FT)
        @test all(isapprox.(q_tot[RH_moist], q_tot_rec[RH_moist], rtol = 5e-2))

        # Update temperature to be exactly consistent with
        # p, ρ, q_pt_rec; test that this is equal to T_rec
        T_local =
            TD.air_temperature_from_ideal_gas_law.(param_set, p, ρ, q_pt_rec)
        @test all(isapprox.(T_local, T_rec, atol = 2 * sqrt(eps(FT))))
    end
end
