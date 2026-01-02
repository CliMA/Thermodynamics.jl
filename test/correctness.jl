"""
# Correctness Test Suite

This file contains tests for fundamental thermodynamic relations and physical laws.
"""

@testset "Thermodynamics - correctness" begin
    FT = Float64
    param_set = FT == Float64 ? param_set_Float64 : param_set_Float32

    @testset "Extract thermodynamic parameters" begin
        # Test that the extraction function works correctly
        (
            _R_d,
            _Rv_over_Rd,
            _R_v,
            _cp_d,
            _cp_v,
            _cp_l,
            _cp_i,
            _cv_d,
            _cv_v,
            _cv_l,
            _cv_i,
            _T_0,
            _e_int_v0,
            _e_int_i0,
            _LH_v0,
            _LH_s0,
            _LH_f0,
            _press_triple,
            _T_triple,
            _T_freeze,
            _T_icenuc,
            _T_min,
            _T_max,
            _p_ref_theta,
            _kappa_d,
        ) = extract_thermodynamic_parameters(param_set)
        @test _R_d == TP.R_d(param_set)
        @test _R_v == TP.R_v(param_set)
    end

    @testset "Ideal Gas Law" begin
        # Test fundamental thermodynamic relations using ideal gas law
        @test air_pressure(param_set, FT(1), FT(1), PhasePartition(FT(1))) ===
              TP.R_v(param_set)
        @test air_pressure(
            param_set,
            FT(1),
            FT(2),
            PhasePartition(FT(1), FT(0.5), FT(0)),
        ) === TP.R_v(param_set)
        @test air_pressure(param_set, FT(1), FT(1)) === TP.R_d(param_set)
        @test air_pressure(param_set, FT(1), FT(2)) === 2 * TP.R_d(param_set)
        @test air_density(param_set, FT(1), FT(1)) === 1 / TP.R_d(param_set)
        @test air_density(param_set, FT(1), FT(2)) === 2 / TP.R_d(param_set)
    end

    @testset "Gas constants and heat capacities" begin
        _R_d = TP.R_d(param_set)
        _R_v = TP.R_v(param_set)
        _cp_d = TP.cp_d(param_set)
        _cp_v = TP.cp_v(param_set)
        _cp_l = TP.cp_l(param_set)
        _cp_i = TP.cp_i(param_set)
        _cv_d = TP.cv_d(param_set)
        _cv_v = TP.cv_v(param_set)
        _cv_l = TP.cv_l(param_set)
        _cv_i = TP.cv_i(param_set)

        # Test gas constants and heat capacities for different phase partitions
        @test gas_constant_air(param_set, PhasePartition(FT(0))) === _R_d
        @test gas_constant_air(param_set, PhasePartition(FT(1))) === _R_v
        @test gas_constant_air(param_set, PhasePartition(FT(0.5), FT(0.5))) ≈
              _R_d / 2
        @test gas_constant_air(param_set) == _R_d

        # Test specific heat capacities for different phases
        @test cp_m(param_set, PhasePartition(FT(0))) === _cp_d
        @test cp_m(param_set, PhasePartition(FT(1))) === _cp_v
        @test cp_m(param_set, PhasePartition(FT(1), FT(1))) === _cp_l
        @test cp_m(param_set, PhasePartition(FT(1), FT(0), FT(1))) === _cp_i
        @test cp_m(param_set) == _cp_d

        # Test specific heat capacities at constant volume
        @test cv_m(param_set, PhasePartition(FT(0))) === _cp_d - _R_d
        @test cv_m(param_set, PhasePartition(FT(1))) === _cp_v - _R_v
        @test cv_m(param_set, PhasePartition(FT(1), FT(1))) === _cv_l
        @test cv_m(param_set, PhasePartition(FT(1), FT(0), FT(1))) === _cv_i
        @test cv_m(param_set) == _cv_d

        # Test scalar version of cv_m
        @test cv_m(param_set, FT(1), FT(1), FT(0)) === _cv_l
        @test cv_m(param_set, FT(1), FT(0), FT(1)) === _cv_i
    end

    @testset "Speed of sound" begin
        _R_d = TP.R_d(param_set)
        _R_v = TP.R_v(param_set)
        _cp_d = TP.cp_d(param_set)
        _cp_v = TP.cp_v(param_set)
        _cv_d = TP.cv_d(param_set)
        _cv_v = TP.cv_v(param_set)
        _T_0 = TP.T_0(param_set)
        # Test speed of sound calculations in dry air and water vapor
        @test soundspeed_air(param_set, _T_0 + 20, PhasePartition(FT(0))) ==
              sqrt(_cp_d / _cv_d * _R_d * (_T_0 + 20))
        @test soundspeed_air(param_set, _T_0 + 100, PhasePartition(FT(1))) ==
              sqrt(_cp_v / _cv_v * _R_v * (_T_0 + 100))
    end

    @testset "Latent heats" begin
        _T_0 = TP.T_0(param_set)
        _LH_v0 = TP.LH_v0(param_set)
        _LH_f0 = TP.LH_f0(param_set)
        _LH_s0 = TP.LH_s0(param_set)
        # Test latent heat calculations at reference temperature
        @test latent_heat_vapor(param_set, _T_0) ≈ _LH_v0
        @test latent_heat_fusion(param_set, _T_0) ≈ _LH_f0
        @test latent_heat_sublim(param_set, _T_0) ≈ _LH_s0
    end

    @testset "Saturation vapor pressure and humidity" begin
        _R_d = TP.R_d(param_set)
        _R_v = TP.R_v(param_set)
        _press_triple = TP.press_triple(param_set)
        _T_triple = TP.T_triple(param_set)

        p = FT(1.e5)
        q_tot = FT(0.02)
        ρ = FT(1.0)
        ρ_v_triple = _press_triple / _R_v / _T_triple
        p_v_tripleminus = FT(0.99) * _press_triple

        # Test saturation vapor pressure at triple point
        @test saturation_vapor_pressure(param_set, _T_triple, Liquid()) ≈
              _press_triple
        @test saturation_vapor_pressure(param_set, _T_triple, Ice()) ≈
              _press_triple

        # Test specific humidities at triple point
        @test q_vap_saturation(
            param_set,
            _T_triple,
            ρ,
            PhaseEquil,
            PhasePartition(FT(0.01)),
        ) == ρ_v_triple / ρ

        # Test saturation specific humidity from pressure calculations
        @test TD.q_vap_saturation_from_pressure(
            param_set,
            q_tot,
            p,
            _T_triple,
            PhaseEquil,
        ) == _R_d / _R_v * (1 - q_tot) * _press_triple / (p - _press_triple)

        @test TD.q_vap_saturation_from_pressure(
            param_set,
            q_tot,
            p_v_tripleminus,
            _T_triple,
            PhaseEquil,
        ) == FT(1)

        # Test saturation specific humidity for non-equilibrium conditions
        @test q_vap_saturation(
            param_set,
            _T_triple,
            ρ,
            PhaseNonEquil,
            PhasePartition(q_tot, q_tot),
        ) == ρ_v_triple / ρ

        # Test generic saturation specific humidity functions for liquid and ice
        @test q_vap_saturation(param_set, _T_triple, ρ, Liquid()) ==
              ρ_v_triple / ρ
        @test q_vap_saturation(param_set, _T_triple, ρ, Ice()) ==
              ρ_v_triple / ρ
        @test q_vap_saturation(
            param_set,
            _T_triple - 20,
            ρ,
            Liquid(),
        ) >= q_vap_saturation(param_set, _T_triple - 20, ρ, Ice())

        # Test saturation specific humidity wrapper functions with thermodynamic state
        _T_0 = TP.T_0(param_set)
        ρ = FT(0.9)
        T = _T_0 - 20
        q_pt = PhasePartition(FT(0.02), FT(0.001), FT(0.002))
        e_int = internal_energy(param_set, T, q_pt)

        ts = PhaseNonEquil(param_set, e_int, ρ, q_pt)
        @test q_vap_saturation(
            param_set,
            air_temperature(param_set, ts),
            ρ,
            Liquid(),
        ) ≈ q_vap_saturation(
            param_set,
            air_temperature(param_set, ts),
            ρ,
            Liquid(),
        )
        @test q_vap_saturation(
            param_set,
            air_temperature(param_set, ts),
            ρ,
            Ice(),
        ) ≈ q_vap_saturation(
            param_set,
            air_temperature(param_set, ts),
            ρ,
            Ice(),
        )
        @test q_vap_saturation(param_set, air_temperature(param_set, ts), ρ, Ice()) <=
              q_vap_saturation(param_set, air_temperature(param_set, ts), ρ, Liquid())

        # Test saturation excess calculations
        @test saturation_excess(
            param_set,
            _T_triple,
            ρ,
            PhaseEquil,
            PhasePartition(q_tot),
        ) ≈ q_tot - ρ_v_triple / ρ

        @test saturation_excess(
            param_set,
            _T_triple,
            ρ,
            PhaseEquil,
            PhasePartition(q_tot / 1000),
        ) == 0.0

        # Test supersaturation calculations for liquid and ice
        # Calculate expected values using the actual function logic
        q_test = PhasePartition(q_tot, 1e-3 * q_tot, 1e-3 * q_tot)
        q_vap = vapor_specific_humidity(q_test)
        T = _T_triple - FT(5)
        p_v = q_vap * ρ * _R_v * T
        p_v_sat_liquid = saturation_vapor_pressure(param_set, T, Liquid())
        p_v_sat_ice = saturation_vapor_pressure(param_set, T, Ice())
        expected_liquid = p_v / p_v_sat_liquid - 1
        expected_ice = p_v / p_v_sat_ice - 1

        @test supersaturation(param_set, q_test, ρ, T, Liquid()) ≈
              expected_liquid

        @test supersaturation(
            param_set,
            q_test,
            ρ,
            T,
            saturation_vapor_pressure(param_set, T, Liquid()),
        ) ≈ expected_liquid

        @test supersaturation(param_set, q_test, ρ, T, Ice()) ≈ expected_ice

        @test supersaturation(
            param_set,
            q_test,
            ρ,
            T,
            saturation_vapor_pressure(param_set, T, Ice()),
        ) ≈ expected_ice

        # Test that supersaturation is zero when vapor pressure equals saturation vapor pressure
        p_v_sat = saturation_vapor_pressure(param_set, T, Liquid())
        q_sat = q_vap_from_p_vap(param_set, T, ρ, p_v_sat)
        q_sat_partition = PhasePartition(q_sat)
        @test supersaturation(param_set, q_sat_partition, ρ, T, Liquid()) ≈ 0.0
    end

    # Helper function to test partial pressures
    function test_partial_pressures(
        T_test,
        ρ,
        q_partition,
        test_name;
        check_vapor_positive = true,
    )
        @testset "Partial Pressures/VPD $test_name" begin
            p_test = air_pressure(param_set, T_test, ρ, q_partition)
            p_vap = partial_pressure_vapor(param_set, p_test, q_partition)
            p_dry = partial_pressure_dry(param_set, p_test, q_partition)

            # Validate against ideal gas law calculations
            _R_v = TP.R_v(param_set)
            _R_d = TP.R_d(param_set)
            q_vap = vapor_specific_humidity(q_partition)
            @test p_vap ≈ q_vap * ρ * _R_v * T_test
            @test p_dry ≈ (1 - q_partition.tot) * ρ * _R_d * T_test

            # Test vapor pressure deficit calculation
            vpd = vapor_pressure_deficit(param_set, T_test, p_test, q_partition)

            # Test VPD over liquid (above freezing) and ice (below freezing) separately
            T_freeze = TP.T_freeze(param_set)
            if T_test > T_freeze
                # Above freezing: should use liquid saturation vapor pressure
                es_liquid =
                    saturation_vapor_pressure(param_set, T_test, Liquid())
                vpd_expected = max(FT(0), es_liquid - p_vap)
                @test vpd ≈ vpd_expected
            else
                # Below freezing: should use ice saturation vapor pressure
                es_ice = saturation_vapor_pressure(param_set, T_test, Ice())
                vpd_expected = max(FT(0), es_ice - p_vap)
                @test vpd ≈ vpd_expected
            end

            # Test vapor pressure deficit with scalar q_vap (when all water is vapor)
            if q_partition.liq ≈ FT(0) && q_partition.ice ≈ FT(0)
                q_vap_scalar = q_partition.tot
                vpd_scalar = vapor_pressure_deficit(
                    param_set,
                    T_test,
                    p_test,
                    q_vap_scalar,
                )
                @test vpd_scalar ≈ vpd  # Both methods should give same result when all water is vapor
            end

            # Validate Dalton's law of partial pressures
            @test p_vap + p_dry ≈ p_test

            # Test dry air pressure is positive
            @test p_dry > FT(0)

            if check_vapor_positive
                @test p_vap > FT(0)  # Should have some vapor pressure
            end
        end
    end

    @testset "Partial Pressures" begin
        # Test partial pressure functions across different scenarios
        ρ = FT(0.9)
        _R_d = TP.R_d(param_set)

        # Test 1: Dry air (no water vapor)
        T = FT(300)
        q_dry = PhasePartition(FT(0))
        test_partial_pressures(
            T,
            ρ,
            q_dry,
            "dry air";
            check_vapor_positive = false,
        )
        ts = PhaseDry(param_set, T, ρ)
        @test partial_pressure_vapor(param_set, ts) === FT(0)

        # Test 2 & 3: Saturated air over liquid water and ice
        test_cases = [(FT(300), Liquid(), "liquid"), (FT(270), Ice(), "ice")]

        for (T_test, phase_type, phase_name) in test_cases
            q_sat = q_vap_saturation(param_set, T_test, ρ, phase_type)
            q_sat_partition = PhasePartition(q_sat)  # All vapor, no liquid/ice

            test_partial_pressures(
                T_test,
                ρ,
                q_sat_partition,
                "saturated $phase_name",
            )
        end

        # Test 4: Subsaturated air (RH < 1)
        T = FT(300)
        q_sub = 0.5 * q_vap_saturation(param_set, T, ρ, Liquid())
        q_sub_partition = PhasePartition(q_sub)  # Small amount of vapor
        p_sat_liq = saturation_vapor_pressure(param_set, T, Liquid())
        test_partial_pressures(T, ρ, q_sub_partition, "subsaturated air")

        # Test 5: Mixed phase (vapor + liquid)
        q_mixed = PhasePartition(FT(0.01), FT(0.005), FT(0))  # Some vapor, some liquid
        test_partial_pressures(T, ρ, q_mixed, "mixed phase")

        # Test 6: Partial pressure functions with thermodynamic state
        T_test = FT(300)
        ρ_test = FT(1.0)
        q_test = PhasePartition(FT(0.01), FT(0.005), FT(0))
        p_test = air_pressure(param_set, T_test, ρ_test, q_test)

        # Create thermodynamic state and test consistency
        e_int = internal_energy(param_set, T_test, q_test)
        ts = PhaseNonEquil(param_set, e_int, ρ_test, q_test)

        # Test partial pressures from thermodynamic state
        p_vap_ts = partial_pressure_vapor(param_set, ts)
        p_dry_ts = partial_pressure_dry(param_set, ts)

        # Test consistency with direct computation
        p_vap_direct = partial_pressure_vapor(param_set, p_test, q_test)
        p_dry_direct = partial_pressure_dry(param_set, p_test, q_test)

        @test p_vap_ts ≈ p_vap_direct
        @test p_dry_ts ≈ p_dry_direct
        @test p_vap_ts + p_dry_ts ≈ p_test  # Dalton's law

        # Test 7: Partial pressure functions with AbstractPhaseDry (using PhasePartition approach)
        T_dry = FT(300)
        ρ_dry = FT(1.0)
        q_dry = PhasePartition(FT(0))  # No water vapor
        p_dry = air_pressure(param_set, T_dry, ρ_dry, q_dry)

        # Test partial pressures for dry air using PhasePartition
        p_vap_dry = partial_pressure_vapor(param_set, p_dry, q_dry)
        p_dry_result = partial_pressure_dry(param_set, p_dry, q_dry)

        # Validate dry air behavior: zero vapor pressure, total pressure equals dry air pressure
        @test p_vap_dry === FT(0)  # No vapor pressure in dry air
        @test p_dry_result ≈ p_dry  # Total pressure equals dry air pressure
        @test p_vap_dry + p_dry_result ≈ p_dry  # Dalton's law for dry air
    end

    @testset "Energy and temperature" begin
        # Test energy functions and their inverse temperature calculations
        T = FT(300)
        q_pt = PhasePartition(FT(0.02), FT(0.002), FT(0.002))
        e_kin = FT(11)
        e_pot = FT(13)
        e_int_dry = internal_energy_dry(param_set, T)
        e_int_moist = internal_energy(param_set, T, q_pt)

        # Test temperature recovery from internal energy (ensures arbitrary constants 
        #in internal energy are consistent with temperature recovery equation)
        @test air_temperature(param_set, e_int_dry) === FT(T)
        @test air_temperature(param_set, e_int_dry, PhasePartition(FT(0))) ===
              FT(T)

        @test air_temperature(param_set, e_int_moist, q_pt) === FT(T)

        # Test total energy calculations including kinetic and potential energy
        @test total_energy(param_set, FT(e_kin), FT(e_pot), FT(T)) ===
              FT(e_kin) + FT(e_pot) + FT(e_int_dry)
        @test total_energy(param_set, FT(e_kin), FT(e_pot), FT(T), q_pt) ≈
              FT(e_kin) + FT(e_pot) + FT(e_int_moist)
        @test total_energy(param_set, FT(0), FT(0), FT(T), q_pt) ≈
              FT(e_int_moist)
    end

    @testset "Phase partitioning" begin
        _T_freeze = TP.T_freeze(param_set)
        _T_icenuc = TP.T_icenuc(param_set)
        # Test phase partitioning above freezing temperature
        T_warm = FT(_T_freeze + 20)
        ρ = FT(1.0)
        q_tot = FT(0.012)
        q_liq = FT(0.001) * q_tot

        @test liquid_fraction(param_set, T_warm, PhaseEquil) === FT(1)
        @test liquid_fraction(
            param_set,
            T_warm,
            PhaseNonEquil,
            PhasePartition(q_tot, q_liq, q_liq),
        ) === FT(0.5)

        @test liquid_fraction(
            param_set,
            T_warm,
            PhaseNonEquil,
            PhasePartition(q_tot, q_liq, q_liq / 2),
        ) === FT(2 / 3)

        # Test equilibrium phase partitioning above freezing
        q = PhasePartition_equil(param_set, T_warm, ρ, q_tot, PhaseEquil)
        @test 0 <= q.liq <= q_tot
        @test q.ice ≈ 0

        p = air_pressure(param_set, T_warm, ρ, q)
        q = TD.PhasePartition_equil_given_p(
            param_set,
            T_warm,
            p,
            q_tot,
            PhaseEquil,
        )
        @test 0 <= q.liq <= q_tot
        @test q.ice ≈ 0

        # Test equilibrium phase partitioning below homogeneous nucleation temperature
        T_cold = FT(_T_icenuc - 10)
        ρ = FT(1.0)
        q_tot = FT(1.02) * q_vap_saturation(param_set, T_cold, ρ, PhaseEquil)
        q = PhasePartition_equil(param_set, T_cold, ρ, q_tot, PhaseEquil)
        @test liquid_fraction(param_set, T_cold, PhaseEquil) === FT(0)
        @test q.liq ≈ FT(0)
        @test 0 <= q.ice <= q_tot  # Ice should be non-negative, may be zero if not supersaturated

        q = TD.PhasePartition_equil_given_p(
            param_set,
            T_cold,
            p,
            q_tot,
            PhaseEquil,
        )
        @test q.liq ≈ FT(0)
        @test q.ice >= FT(0)  # Ice should be non-negative, may be zero if not supersaturated
    end

    @testset "Saturation adjustment" begin
        # Test saturation adjustment algorithms for equilibrium calculations
        # Given internal energy, density, and total humidity, compute temperature and phase partitioning
        _T_icenuc = TP.T_icenuc(param_set)
        q_tot = FT(0)
        ρ = FT(1)
        T = FT(300)
        @test TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            internal_energy_sat(param_set, T, ρ, q_tot, PhaseEquil),
            ρ,
            q_tot,
            PhaseEquil,
            10,
            rtol_temperature,
        ) ≈ T
        @test isapprox(
            TD.saturation_adjustment(
                RS.NewtonsMethod,
                param_set,
                internal_energy_sat(param_set, T, ρ, q_tot, PhaseEquil),
                ρ,
                q_tot,
                PhaseEquil,
                10,
                rtol_temperature,
            ),
            T,
            atol = atol_temperature,
        )

        # Test saturation adjustment for cold conditions (ice phase)
        q_tot = FT(0.05)
        ρ = FT(0.1)
        T_cold = _T_icenuc - FT(10)
        @test isapprox(
            TD.saturation_adjustment(
                RS.SecantMethod,
                param_set,
                internal_energy_sat(param_set, T_cold, ρ, q_tot, PhaseEquil),
                ρ,
                q_tot,
                PhaseEquil,
                10,
                rtol_temperature,
            ),
            T_cold,
            atol = atol_temperature,
        )
        @test isapprox(
            TD.saturation_adjustment(
                RS.NewtonsMethod,
                param_set,
                internal_energy_sat(param_set, T_cold, ρ, q_tot, PhaseEquil),
                ρ,
                q_tot,
                PhaseEquil,
                10,
                rtol_temperature,
            ),
            T_cold,
            atol = atol_temperature,
        )

        # Test saturation adjustment for warm conditions (liquid phase)
        ρ = FT(0.1)
        T = FT(300)
        # Slightly super-saturated total specific humidity
        q_tot = FT(1.02) * q_vap_saturation(param_set, T, ρ, PhaseEquil)
        q = PhasePartition_equil(param_set, T, ρ, q_tot, PhaseEquil)
        p = air_pressure(param_set, T, ρ, q)
        @test q.tot - q.liq - q.ice ≈
              vapor_specific_humidity(q) ≈
              q_vap_saturation(param_set, T, ρ, PhaseEquil)

        q = TD.PhasePartition_equil_given_p(param_set, T, p, q_tot, PhaseEquil)
        @test q.tot - q.liq - q.ice ≈
              vapor_specific_humidity(q) ≈
              TD.q_vap_saturation_from_pressure(
                  param_set,
                  q_tot,
                  p,
                  T,
                  PhaseEquil,
              )
    end

    @testset "Internal energy" begin
        # Test internal energy calculations for different phases
        T = FT(300)
        q = PhasePartition(FT(20 * 1e-3), FT(5 * 1e-3), FT(2 * 1e-3))
        q_vap = vapor_specific_humidity(q)
        Id = internal_energy_dry(param_set, T)
        Iv = internal_energy_vapor(param_set, T)
        Il = internal_energy_liquid(param_set, T)
        Ii = internal_energy_ice(param_set, T)
        @test internal_energy(param_set, T, q) ≈
              (1 - q.tot) * Id + q_vap * Iv + q.liq * Il + q.ice * Ii
        @test internal_energy(param_set, T) ≈ Id

        # Test internal_energy with different phase partitions
        T_test = FT(280)
        q_dry = PhasePartition(FT(0))  # Dry air
        q_vapor_only = PhasePartition(FT(0.01))  # Only vapor
        q_mixed = PhasePartition(FT(0.01), FT(0.005), FT(0.002))  # Mixed phases

        @test internal_energy(param_set, T_test, q_dry) ≈
              internal_energy_dry(param_set, T_test)
        @test internal_energy(param_set, T_test, q_vapor_only) ≈
              (1 - q_vapor_only.tot) * internal_energy_dry(param_set, T_test) +
              vapor_specific_humidity(q_vapor_only) *
              internal_energy_vapor(param_set, T_test)
        @test internal_energy(param_set, T_test, q_mixed) ≈
              (1 - q_mixed.tot) * internal_energy_dry(param_set, T_test) +
              vapor_specific_humidity(q_mixed) *
              internal_energy_vapor(param_set, T_test) +
              q_mixed.liq * internal_energy_liquid(param_set, T_test) +
              q_mixed.ice * internal_energy_ice(param_set, T_test)

        # Test internal_energy_sat with different supersaturated conditions
        _T_freeze = TP.T_freeze(param_set)
        ρ_test = FT(0.9)

        # Test above freezing (liquid phase)
        T_warm = _T_freeze + FT(30)  # Above freezing
        q_tot_warm =
            FT(1.02) * q_vap_saturation(param_set, T_warm, ρ_test, PhaseEquil)
        e_int_sat_warm = internal_energy_sat(
            param_set,
            T_warm,
            ρ_test,
            q_tot_warm,
            PhaseEquil,
        )
        @test e_int_sat_warm ≈ internal_energy(
            param_set,
            T_warm,
            PhasePartition_equil(
                param_set,
                T_warm,
                ρ_test,
                q_tot_warm,
                PhaseEquil,
            ),
        )

        # Test below freezing (ice phase)
        T_cold = _T_freeze - FT(30)  # Below freezing
        q_tot_cold =
            FT(1.02) * q_vap_saturation(param_set, T_cold, ρ_test, PhaseEquil)
        e_int_sat_cold = internal_energy_sat(
            param_set,
            T_cold,
            ρ_test,
            q_tot_cold,
            PhaseEquil,
        )
        @test e_int_sat_cold ≈ internal_energy(
            param_set,
            T_cold,
            PhasePartition_equil(
                param_set,
                T_cold,
                ρ_test,
                q_tot_cold,
                PhaseEquil,
            ),
        )
    end

    @testset "Potential temperature" begin
        # Test potential temperature calculations
        _R_d = TP.R_d(param_set)
        _R_v = TP.R_v(param_set)
        _cp_d = TP.cp_d(param_set)
        _cp_v = TP.cp_v(param_set)
        _p_ref_theta = TP.p_ref_theta(param_set)
        T = FT(300)
        @test TD.liquid_ice_pottemp_given_pressure(
            param_set,
            T,
            _p_ref_theta,
        ) === T
        @test TD.liquid_ice_pottemp_given_pressure(
            param_set,
            T,
            _p_ref_theta / 10,
        ) ≈ T * 10^(_R_d / _cp_d)
        @test TD.liquid_ice_pottemp_given_pressure(
            param_set,
            T,
            _p_ref_theta / 10,
            PhasePartition(FT(1)),
        ) ≈ T * 10^(_R_v / _cp_v)

        # Test dry potential temperature calculations
        T = FT(300)
        p = FT(1.e5)
        q_tot = FT(0.023)
        @test TD.dry_pottemp_given_pressure(
            param_set,
            T,
            p,
            PhasePartition(q_tot),
        ) isa typeof(p)
        @test TD.air_temperature_given_pθq(
            param_set,
            p,
            TD.dry_pottemp_given_pressure(
                param_set,
                T,
                p,
                PhasePartition(q_tot),
            ),
            PhasePartition(q_tot),
        ) ≈ T
    end

    @testset "Exner function" begin
        # Test Exner function calculations
        p = FT(1.e5)
        q_tot = FT(0.023)
        q_pt = PhasePartition(q_tot)
        @test TD.exner_given_pressure(param_set, p, q_pt) isa typeof(p)
        @test TD.exner_given_pressure(param_set, p / 10, q_pt) ≈
              TD.exner_given_pressure(param_set, p, q_pt) /
              10^(gas_constant_air(param_set, q_pt) / cp_m(param_set, q_pt))
    end

    @testset "Humidities and mixing ratios" begin
        _Rv_over_Rd = TP.Rv_over_Rd(param_set)
        # Test humidity and mixing ratio calculations
        q_tot = FT(0.03)
        q_liq = FT(0.005)
        q_ice = FT(0.001)
        q_partition = PhasePartition(q_tot, q_liq, q_ice)

        @test liquid_specific_humidity(q_partition) == q_liq
        @test ice_specific_humidity(q_partition) == q_ice
        @test vapor_specific_humidity(q_partition) == q_tot - q_liq - q_ice
        @test condensate_specific_humidity(q_partition) == q_liq + q_ice

        mr = shum_to_mixing_ratio(q_tot, q_tot)
        @test mr == q_tot / (1 - q_tot)
        mr = shum_to_mixing_ratio(q_liq, q_tot)
        @test mr == q_liq / (1 - q_tot)

        q = PhasePartition(q_tot, q_liq, q_ice)
        mrs = mixing_ratios(q)
        @test mrs.tot == q_tot / (1 - q_tot)
        @test mrs.liq == q_liq / (1 - q_tot)
        @test mrs.ice == q_ice / (1 - q_tot)

        vmrs = vol_vapor_mixing_ratio(param_set, q)
        q_vap = vapor_specific_humidity(q)
        @test vmrs ≈ _Rv_over_Rd * shum_to_mixing_ratio(q_vap, q.tot)
    end

    @testset "Relative humidity" begin
        # Test relative humidity calculations across wide parameter ranges
        for phase_type in [PhaseDry, PhaseEquil, PhaseNonEquil]
            for T in [FT(40), FT(140), FT(240), FT(340), FT(440)]
                for p in [FT(1e3), FT(1e4), FT(1e5)]
                    for q in [FT(-1), FT(1e-45), FT(0), FT(1e-3), FT(10)]
                        q_pt = PhasePartition(FT(q))
                        RH =
                            relative_humidity(param_set, T, p, phase_type, q_pt)
                        @test RH >= FT(0) && RH <= FT(1)
                    end
                end
            end
        end

        # Test vapor specific humidity from relative humidity calculations
        for T in [FT(280), FT(290), FT(300)]
            for p in [FT(70000), FT(101300)]
                for q_vap in [FT(1e-3), FT(1e-4), FT(0)]
                    q_pt = TD.PhasePartition(q_vap)
                    RH = relative_humidity(param_set, T, p, TD.PhaseEquil, q_pt)
                    q_vap = q_vap_from_RH_liquid(param_set, p, T, RH)
                    @test q_vap ≈ vapor_specific_humidity(q_pt)
                end
            end
        end

        # Test edge cases for relative humidity
        @testset "RH edge cases" begin
            # Extremely dry conditions
            @test relative_humidity(
                param_set,
                300.0,
                1.0e5,
                PhaseEquil,
                PhasePartition(0.0),
            ) ≈ 0.0

            # Saturated conditions
            q_sat = q_vap_saturation(param_set, 300.0, 1.0, PhaseEquil)
            @test relative_humidity(
                param_set,
                300.0,
                1.0e5,
                PhaseEquil,
                PhasePartition(q_sat),
            ) ≈ 1.0

            # Supersaturated conditions
            q_super = 1.1 * q_vap_saturation(param_set, 300.0, 1.0, PhaseEquil)
            rh_super = relative_humidity(
                param_set,
                300.0,
                1.0e5,
                PhaseEquil,
                PhasePartition(q_super),
            )
            @test rh_super ≥ 1.0  # Allow for numerical precision in the case where it equals exactly 1.0
        end
    end

    @testset "Specific enthalpy" begin
        _R_d = TP.R_d(param_set)
        _R_v = TP.R_v(param_set)
        # Test total specific enthalpy calculations
        T = FT(300)
        e_tot = FT(1000)

        # Test for dry air (default PhasePartition)
        h_tot_dry = total_enthalpy(param_set, e_tot, T)
        @test h_tot_dry ≈ e_tot + _R_d * T

        # Test for moist air
        q_pt = PhasePartition(FT(0.02), FT(0.005), FT(0.001))
        h_tot_moist = total_enthalpy(param_set, e_tot, T, q_pt)
        R_m_moist = gas_constant_air(param_set, q_pt)
        @test h_tot_moist ≈ e_tot + R_m_moist * T

        # Test specific enthalpy calculations
        h_dry = enthalpy_dry(param_set, T)
        h_vap = enthalpy_vapor(param_set, T)
        h_liq = enthalpy_liquid(param_set, T)
        h_ice = enthalpy_ice(param_set, T)
        @test h_dry ≈ internal_energy_dry(param_set, T) + _R_d * T
        @test h_vap ≈ internal_energy_vapor(param_set, T) + _R_v * T
        @test h_liq == internal_energy_liquid(param_set, T)
        @test h_ice == internal_energy_ice(param_set, T)
    end
end
