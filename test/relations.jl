"""
# Thermodynamics Relations Test Suite

This file contains comprehensive tests for the Thermodynamics.jl package, covering:

## Test Categories:
- **Correctness**: Consistency and correctness of fundamental thermodynamic relations and physical laws
- **Thermodynamic States**: State constructors and consistency checks
- **Performance**: Type stability and computational efficiency

## Key Features Tested:
- Ideal gas law and thermodynamic consistency
- Gas constants and heat capacities for different phases
- Saturation vapor pressure and humidity relations
- Partial pressure calculations (vapor and dry air)
- Internal energy calculations for dry, vapor, liquid, and ice phases
- Thermodynamic state constructors and their consistency
- Phase partitioning and equilibrium calculations
- Performance across different floating-point types (Float32, Float64)

## Test Structure:
- Uses parameter sets for different thermodynamic configurations
- Tests both scalar and array-based calculations
- Includes tolerance-based comparisons for floating-point accuracy
- Covers edge cases with unusual thermodynamic states
"""

using Test

using Random
import RootSolvers as RS
using LinearAlgebra
import ForwardDiff

using Thermodynamics
import Thermodynamics as TD
import Thermodynamics.Parameters as TP
using Thermodynamics.TemperatureProfiles
using Thermodynamics.TestedProfiles
import ClimaParams as CP

# Tolerances for tested quantities:
param_set_Float64 = TP.ThermodynamicsParameters(Float64)
atol_temperature = 5e-1
atol_energy = TP.cv_d(param_set_Float64) * atol_temperature
rtol_temperature = 1e-1
rtol_density = rtol_temperature
rtol_pressure = 1e-1
rtol_energy = 1e-1

# Test both Float32 and Float64 for type stability
array_types = [Array{Float32}, Array{Float64}]

include("data_tests.jl")

# Helper function to compare moisture content between thermodynamic states
compare_moisture(param_set, a::ThermodynamicState, b::ThermodynamicState) =
    compare_moisture(param_set, a, PhasePartition(param_set, b))

# Compare total moisture for equilibrium states
compare_moisture(param_set, ts::PhaseEquil, q_pt::PhasePartition) =
    getproperty(PhasePartition(param_set, ts), :tot) ≈ getproperty(q_pt, :tot)

# Compare all moisture components for non-equilibrium states
compare_moisture(param_set, ts::PhaseNonEquil, q_pt::PhasePartition) = all((
    getproperty(PhasePartition(param_set, ts), :tot) ≈ getproperty(q_pt, :tot),
    getproperty(PhasePartition(param_set, ts), :liq) ≈ getproperty(q_pt, :liq),
    getproperty(PhasePartition(param_set, ts), :ice) ≈ getproperty(q_pt, :ice),
))

@testset "Thermodynamics - isentropic processes" begin
    for ArrayType in array_types
        FT = eltype(ArrayType)
        param_set = TP.ThermodynamicsParameters(FT)

        # Extract thermodynamic parameters for testing
        # Gas constants
        _R_d = FT(TP.R_d(param_set))           # Dry air gas constant
        _Rv_over_Rd = FT(TP.Rv_over_Rd(param_set))  # Vapor/dry air gas constant ratio
        _R_v = FT(TP.R_v(param_set))           # Water vapor gas constant

        # Heat capacities
        _cp_d = FT(TP.cp_d(param_set))         # Dry air specific heat at constant pressure
        _cp_v = FT(TP.cp_v(param_set))         # Water vapor specific heat at constant pressure
        _cp_l = FT(TP.cp_l(param_set))         # Liquid water specific heat
        _cp_i = FT(TP.cp_i(param_set))         # Ice specific heat
        _cv_d = FT(TP.cv_d(param_set))         # Dry air specific heat at constant volume
        _cv_v = FT(TP.cv_v(param_set))         # Water vapor specific heat at constant volume
        _cv_l = FT(TP.cv_l(param_set))         # Liquid water specific heat at constant volume
        _cv_i = FT(TP.cv_i(param_set))         # Ice specific heat at constant volume

        # Reference temperatures and energies
        _T_0 = FT(TP.T_0(param_set))           # Reference temperature
        _e_int_v0 = FT(TP.e_int_v0(param_set)) # Reference vapor internal energy
        _e_int_i0 = FT(TP.e_int_i0(param_set)) # Reference ice internal energy

        # Latent heats
        _LH_v0 = FT(TP.LH_v0(param_set))       # Latent heat of vaporization
        _LH_s0 = FT(TP.LH_s0(param_set))       # Latent heat of sublimation
        _LH_f0 = FT(TP.LH_f0(param_set))       # Latent heat of fusion

        # Triple point and phase transition temperatures
        _press_triple = FT(TP.press_triple(param_set))  # Triple point pressure
        _T_triple = FT(TP.T_triple(param_set))          # Triple point temperature
        _T_freeze = FT(TP.T_freeze(param_set))          # Freezing temperature
        _T_icenuc = FT(TP.T_icenuc(param_set))          # Homogeneous ice nucleation temperature

        # Temperature and pressure ranges
        _T_min = FT(TP.T_min(param_set))       # Minimum temperature
        _T_max = FT(TP.T_max(param_set))       # Maximum temperature
        _p_ref_theta = FT(TP.p_ref_theta(param_set))  # Reference pressure for potential temperature
        _kappa_d = FT(TP.kappa_d(param_set))   # Dry air adiabatic exponent

        profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
        (; T, p, e_int, ρ, θ_liq_ice, phase_type) = profiles
        (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = profiles

        # Test ideal gas law consistency across all profiles
        T_idgl = TD.air_temperature_from_ideal_gas_law.(param_set, p, ρ, q_pt)
        @test all(T .≈ T_idgl)

        Φ = FT(1)
        Random.seed!(15)
        perturbation = FT(0.1) * rand(FT, length(T))

        # Test adiabatic processes with perturbed ambient conditions
        T∞, p∞ = T .* perturbation, p .* perturbation
        @test air_temperature.(param_set, p, θ_liq_ice, DryAdiabaticProcess()) ≈
              (p ./ _p_ref_theta) .^ (_R_d / _cp_d) .* θ_liq_ice
        @test TD.air_pressure_given_θ.(
            param_set,
            θ_liq_ice,
            Φ,
            DryAdiabaticProcess(),
        ) ≈ _p_ref_theta .* (1 .- Φ ./ (θ_liq_ice .* _cp_d)) .^ (_cp_d / _R_d)
        @test air_pressure.(param_set, T, T∞, p∞, DryAdiabaticProcess()) ≈
              p∞ .* (T ./ T∞) .^ (FT(1) / _kappa_d)
    end
end

@testset "Thermodynamics - correctness" begin
    FT = Float64
    param_set = TP.ThermodynamicsParameters(FT)
    _R_d = FT(TP.R_d(param_set))
    _Rv_over_Rd = FT(TP.Rv_over_Rd(param_set))
    _cp_d = FT(TP.cp_d(param_set))
    _cp_v = FT(TP.cp_v(param_set))
    _cp_l = FT(TP.cp_l(param_set))
    _cv_d = FT(TP.cv_d(param_set))
    _cv_v = FT(TP.cv_v(param_set))
    _cv_l = FT(TP.cv_l(param_set))
    _cv_i = FT(TP.cv_i(param_set))
    _T_0 = FT(TP.T_0(param_set))
    _e_int_v0 = FT(TP.e_int_v0(param_set))
    _e_int_i0 = FT(TP.e_int_i0(param_set))
    _LH_v0 = FT(TP.LH_v0(param_set))
    _LH_s0 = FT(TP.LH_s0(param_set))
    _cp_i = FT(TP.cp_i(param_set))
    _LH_f0 = FT(TP.LH_f0(param_set))
    _press_triple = FT(TP.press_triple(param_set))
    _R_v = FT(TP.R_v(param_set))
    _T_triple = FT(TP.T_triple(param_set))
    _T_freeze = FT(TP.T_freeze(param_set))
    _T_min = FT(TP.T_min(param_set))
    _p_ref_theta = FT(TP.p_ref_theta(param_set))
    _T_max = FT(TP.T_max(param_set))
    _kappa_d = FT(TP.kappa_d(param_set))
    _T_icenuc = FT(TP.T_icenuc(param_set))

    # Test fundamental thermodynamic relations using ideal gas law
    @test air_pressure(param_set, FT(1), FT(1), PhasePartition(FT(1))) === _R_v
    @test air_pressure(
        param_set,
        FT(1),
        FT(2),
        PhasePartition(FT(1), FT(0.5), FT(0)),
    ) === _R_v
    @test air_pressure(param_set, FT(1), FT(1)) === _R_d
    @test air_pressure(param_set, FT(1), FT(2)) === 2 * _R_d
    @test air_density(param_set, FT(1), FT(1)) === 1 / _R_d
    @test air_density(param_set, FT(1), FT(2)) === 2 / _R_d

    # Test gas constants and heat capacities for different phase partitions
    @test gas_constant_air(param_set, PhasePartition(FT(0))) === _R_d
    @test gas_constant_air(param_set, PhasePartition(FT(1))) === _R_v
    @test gas_constant_air(param_set, PhasePartition(FT(0.5), FT(0.5))) ≈
          _R_d / 2
    @test gas_constant_air(param_set, FT) == _R_d

    # Test backward compatibility for renamed parameters
    @test TP.molmass_ratio(param_set) === _Rv_over_Rd

    # Test specific heat capacities for different phases
    @test cp_m(param_set, PhasePartition(FT(0))) === _cp_d
    @test cp_m(param_set, PhasePartition(FT(1))) === _cp_v
    @test cp_m(param_set, PhasePartition(FT(1), FT(1))) === _cp_l
    @test cp_m(param_set, PhasePartition(FT(1), FT(0), FT(1))) === _cp_i
    @test cp_m(param_set, FT) == _cp_d

    # Test specific heat capacities at constant volume
    @test cv_m(param_set, PhasePartition(FT(0))) === _cp_d - _R_d
    @test cv_m(param_set, PhasePartition(FT(1))) === _cp_v - _R_v
    @test cv_m(param_set, PhasePartition(FT(1), FT(1))) === _cv_l
    @test cv_m(param_set, PhasePartition(FT(1), FT(0), FT(1))) === _cv_i
    @test cv_m(param_set, FT) == _cv_d

    # Test speed of sound calculations for different phases
    @test soundspeed_air(param_set, _T_0 + 20, PhasePartition(FT(0))) ==
          sqrt(_cp_d / _cv_d * _R_d * (_T_0 + 20))
    @test soundspeed_air(param_set, _T_0 + 100, PhasePartition(FT(1))) ==
          sqrt(_cp_v / _cv_v * _R_v * (_T_0 + 100))

    # Test latent heat calculations at reference temperature
    @test latent_heat_vapor(param_set, _T_0) ≈ _LH_v0
    @test latent_heat_fusion(param_set, _T_0) ≈ _LH_f0
    @test latent_heat_sublim(param_set, _T_0) ≈ _LH_s0

    # Test saturation vapor pressure and humidity calculations
    p = FT(1.e5)
    q_tot = FT(0.23)
    ρ = FT(1.0)
    ρ_v_triple = _press_triple / _R_v / _T_triple
    p_v_tripleminus = FT(0.99) * _press_triple

    # Test saturation vapor pressure at triple point
    @test saturation_vapor_pressure(param_set, _T_triple, Liquid()) ≈
          _press_triple
    @test saturation_vapor_pressure(param_set, _T_triple, Ice()) ≈ _press_triple

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

    # Test saturation specific humidityfor non-equilibrium conditions
    @test q_vap_saturation(
        param_set,
        _T_triple,
        ρ,
        PhaseNonEquil,
        PhasePartition(q_tot, q_tot),
    ) == ρ_v_triple / ρ

    # Test generic saturation specific humidity functions for liquid and ice
    @test q_vap_saturation_generic(param_set, _T_triple, ρ, Liquid()) ==
          ρ_v_triple / ρ
    @test q_vap_saturation_generic(param_set, _T_triple, ρ, Ice()) ==
          ρ_v_triple / ρ
    @test q_vap_saturation_generic(param_set, _T_triple - 20, ρ, Liquid()) >=
          q_vap_saturation_generic(param_set, _T_triple - 20, ρ, Ice())

    # Test saturation specific humidity wrapper functions with thermodynamic state
    ρ = FT(1)
    ρu = FT[1, 2, 3]
    ρe = FT(1100) - _R_d * _T_0
    e_pot = FT(93)
    e_int = internal_energy(ρ, ρe, ρu, e_pot)
    q_pt = PhasePartition(FT(0.02), FT(0.002), FT(0.002))
    ts = PhaseNonEquil(param_set, e_int, ρ, q_pt)
    @test q_vap_saturation_generic(
        param_set,
        air_temperature(param_set, ts),
        ρ,
        Liquid(),
    ) ≈ q_vap_saturation_liquid(param_set, ts)
    @test q_vap_saturation_generic(
        param_set,
        air_temperature(param_set, ts),
        ρ,
        Ice(),
    ) ≈ q_vap_saturation_ice(param_set, ts)
    @test q_vap_saturation_ice(param_set, ts) <=
          q_vap_saturation_liquid(param_set, ts)

    # Test saturation excess calculations
    @test saturation_excess(
        param_set,
        _T_triple,
        ρ,
        PhaseEquil,
        PhasePartition(q_tot),
    ) == q_tot - ρ_v_triple / ρ
    @test saturation_excess(
        param_set,
        _T_triple,
        ρ,
        PhaseEquil,
        PhasePartition(q_tot / 1000),
    ) == 0.0

    # Test supersaturation calculations for liquid and ice
    @test supersaturation(
        param_set,
        PhasePartition(q_tot, 1e-3 * q_tot, 1e-3 * q_tot),
        ρ,
        _T_triple,
        Liquid(),
    ) ≈ 0.998 * q_tot / ρ_v_triple / ρ - 1
    @test supersaturation(
        param_set,
        PhasePartition(q_tot, 1e-3 * q_tot, 1e-3 * q_tot),
        ρ,
        _T_triple,
        saturation_vapor_pressure(param_set, _T_triple, Liquid()),
    ) ≈ 0.998 * q_tot / ρ_v_triple / ρ - 1

    @test supersaturation(
        param_set,
        PhasePartition(q_tot, 1e-3 * q_tot, 1e-3 * q_tot),
        ρ,
        _T_triple,
        Ice(),
    ) ≈ 0.998 * q_tot / ρ_v_triple / ρ - 1
    @test supersaturation(
        param_set,
        PhasePartition(q_tot, 1e-3 * q_tot, 1e-3 * q_tot),
        ρ,
        _T_triple,
        saturation_vapor_pressure(param_set, _T_triple, Ice()),
    ) ≈ 0.998 * q_tot / ρ_v_triple / ρ - 1

    # Helper function to test partial pressures
    function test_partial_pressures(
        T_test,
        ρ,
        q_partition,
        test_name;
        check_vapor_positive = true,
    )

        p_test = air_pressure(param_set, T_test, ρ, q_partition)
        p_vap = partial_pressure_vapor(param_set, p_test, q_partition)
        p_dry = partial_pressure_dry(param_set, p_test, q_partition)

        # Validate against ideal gas law calculations
        q_vap = vapor_specific_humidity(q_partition)
        @test p_vap ≈ q_vap * ρ * _R_v * T_test
        @test p_dry ≈ (1 - q_partition.tot) * ρ * _R_d * T_test

        # Test vapor pressure deficit calculation
        vpd = vapor_pressure_deficit_liquid(
            param_set,
            T_test,
            p_test,
            q_partition,
        )
        es = saturation_vapor_pressure(param_set, T_test, Liquid())
        @test vpd ≈ max(FT(0), es - p_vap)  # Account for ReLU function

        # Test vapor pressure deficit with scalar q_vap (when all water is vapor)
        if q_partition.liq ≈ FT(0) && q_partition.ice ≈ FT(0)
            q_vap_scalar = q_partition.tot
            vpd_scalar = vapor_pressure_deficit_liquid(
                param_set,
                T_test,
                p_test,
                q_vap_scalar,
            )
            @test vpd_scalar ≈ vpd  # Both methods should give same result when all water is vapor
        end

        # Validate Dalton's law of partial pressures
        @test p_vap + p_dry ≈ p_test

        if check_vapor_positive
            @test p_vap > FT(0)  # Should have some vapor pressure
        end
    end

    # Test partial pressure functions across different scenarios
    ρ = FT(1)

    # Test 1: Dry air (no water vapor)
    T = FT(300)
    q_dry = PhasePartition(FT(0))
    test_partial_pressures(T, ρ, q_dry, "dry air"; check_vapor_positive = false)

    # Test 2 & 3: Saturated air over liquid water and ice
    test_cases = [(FT(300), Liquid(), "liquid"), (FT(270), Ice(), "ice")]

    for (T_test, phase_type, phase_name) in test_cases
        q_sat = q_vap_saturation_generic(param_set, T_test, ρ, phase_type)
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
    q_sub = 0.5 * q_vap_saturation_generic(param_set, T, ρ, Liquid())
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

    # Test 8: Partial pressure functions with some vapor present
    T_moist = FT(300)
    ρ_moist = FT(1.0)
    q_moist = PhasePartition(FT(0.01))  # Some water vapor
    p_moist = air_pressure(param_set, T_moist, ρ_moist, q_moist)

    # Test partial pressures for moist air
    p_vap_moist = partial_pressure_vapor(param_set, p_moist, q_moist)
    p_dry_moist = partial_pressure_dry(param_set, p_moist, q_moist)

    # Validate moist air behavior: positive vapor pressure, reduced dry air pressure
    @test p_vap_moist > FT(0)  # Should have some vapor pressure
    @test p_dry_moist < p_moist  # Dry air pressure should be less than total
    @test p_vap_moist + p_dry_moist ≈ p_moist  # Dalton's law for moist air
    @test p_dry_moist > FT(0)  # Dry air pressure should still be positive

    # Test energy functions and their inverse temperature calculations
    T = FT(300)
    q_pt = PhasePartition(FT(0.02), FT(0.002), FT(0.002))
    e_kin = FT(11)
    e_pot = FT(13)
    e_int_dry = internal_energy_dry(param_set, T)
    e_int_moist = internal_energy(param_set, T, q_pt)

    # Test temperature recovery from internal energy
    @test air_temperature(param_set, e_int_dry) === FT(T)
    @test air_temperature(param_set, e_int_dry, PhasePartition(FT(0))) === FT(T)

    @test air_temperature(param_set, e_int_moist, q_pt) === FT(T)

    # Test total energy calculations including kinetic and potential energy
    @test total_energy(param_set, FT(e_kin), FT(e_pot), FT(T)) ===
          FT(e_kin) + FT(e_pot) + FT(e_int_dry)
    @test total_energy(param_set, FT(e_kin), FT(e_pot), FT(T), q_pt) ≈
          FT(e_kin) + FT(e_pot) + FT(e_int_moist)
    @test total_energy(param_set, FT(0), FT(0), FT(T), q_pt) ≈ FT(e_int_moist)

    # Test phase partitioning above freezing temperature
    T_warm = FT(_T_freeze + 20)
    ρ = FT(1.0)
    q_tot = FT(0.21)
    q_liq = FT(0.01) * q_tot

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
    @test 0 < q.liq <= q_tot
    @test q.ice ≈ 0

    p = air_pressure(param_set, T_warm, ρ, q)
    q = TD.PhasePartition_equil_given_p(param_set, T_warm, p, q_tot, PhaseEquil)
    @test 0 < q.liq <= q_tot
    @test q.ice ≈ 0

    # Test equilibrium phase partitioning below homogeneous nucleation temperature
    T_cold = FT(_T_icenuc - 10)
    ρ = FT(1.0)
    q_tot = FT(1.02) * q_vap_saturation(param_set, T_cold, ρ, PhaseEquil)
    q = PhasePartition_equil(param_set, T_cold, ρ, q_tot, PhaseEquil)
    @test liquid_fraction(param_set, T_cold, PhaseEquil) === FT(0)
    @test q.liq ≈ FT(0)
    @test q.ice >= FT(0)  # Ice should be non-negative, may be zero if not supersaturated

    q = TD.PhasePartition_equil_given_p(param_set, T_cold, p, q_tot, PhaseEquil)
    @test q.liq ≈ FT(0)
    @test q.ice >= FT(0)  # Ice should be non-negative, may be zero if not supersaturated

    # Test saturation adjustment algorithms for equilibrium calculations
    # Given internal energy, density, and total humidity, compute temperature and phase partitioning
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
    @test abs(
        TD.saturation_adjustment(
            RS.NewtonsMethod,
            param_set,
            internal_energy_sat(param_set, T, ρ, q_tot, PhaseEquil),
            ρ,
            q_tot,
            PhaseEquil,
            10,
            rtol_temperature,
        ) - T,
    ) < rtol_temperature

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
        rtol = rtol_temperature,
    )
    @test abs(
        TD.saturation_adjustment(
            RS.NewtonsMethod,
            param_set,
            internal_energy_sat(param_set, T_cold, ρ, q_tot, PhaseEquil),
            ρ,
            q_tot,
            PhaseEquil,
            10,
            rtol_temperature,
        ) - T_cold,
    ) < rtol_temperature

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
          TD.q_vap_saturation_from_pressure(param_set, q_tot, p, T, PhaseEquil)

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
    ρ_test = FT(1.0)

    # Test above freezing (liquid phase)
    T_warm = _T_freeze + FT(30)  # Above freezing
    q_tot_warm =
        FT(1.02) * q_vap_saturation(param_set, T_warm, ρ_test, PhaseEquil)
    e_int_sat_warm =
        internal_energy_sat(param_set, T_warm, ρ_test, q_tot_warm, PhaseEquil)
    @test e_int_sat_warm ≈ internal_energy(
        param_set,
        T_warm,
        PhasePartition_equil(param_set, T_warm, ρ_test, q_tot_warm, PhaseEquil),
    )

    # Test below freezing (ice phase)
    T_cold = _T_freeze - FT(30)  # Below freezing
    q_tot_cold =
        FT(1.02) * q_vap_saturation(param_set, T_cold, ρ_test, PhaseEquil)
    e_int_sat_cold =
        internal_energy_sat(param_set, T_cold, ρ_test, q_tot_cold, PhaseEquil)
    @test e_int_sat_cold ≈ internal_energy(
        param_set,
        T_cold,
        PhasePartition_equil(param_set, T_cold, ρ_test, q_tot_cold, PhaseEquil),
    )

    # Test potential temperature calculations
    T = FT(300)
    @test TD.liquid_ice_pottemp_given_pressure(param_set, T, _p_ref_theta) === T
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
    q_tot = FT(0.23)
    @test TD.dry_pottemp_given_pressure(
        param_set,
        T,
        p,
        PhasePartition(q_tot),
    ) isa typeof(p)
    @test TD.air_temperature_given_pθq(
        param_set,
        p,
        TD.dry_pottemp_given_pressure(param_set, T, p, PhasePartition(q_tot)),
        PhasePartition(q_tot),
    ) ≈ T

    # Test Exner function calculations
    p = FT(1.e5)
    q_tot = FT(0.23)
    q_pt = PhasePartition(q_tot)
    @test TD.exner_given_pressure(param_set, p, q_pt) isa typeof(p)
    @test TD.exner_given_pressure(param_set, p / 10, q_pt) ≈
          TD.exner_given_pressure(param_set, p, q_pt) /
          10^(gas_constant_air(param_set, q_pt) / cp_m(param_set, q_pt))

    # Test humidity and mixing ratio calculations
    q_tot = 0.1
    q_liq = 0.05
    q_ice = 0.01
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

    # Test relative humidity calculations across wide parameter ranges
    for phase_type in [PhaseDry, PhaseEquil, PhaseNonEquil]
        for T in [FT(40), FT(140), FT(240), FT(340), FT(440)]
            for p in [FT(1e3), FT(1e4), FT(1e5)]
                for q in [FT(-1), FT(1e-45), FT(0), FT(1e-3), FT(10)]
                    q_pt = PhasePartition(FT(q))
                    RH = relative_humidity(param_set, T, p, phase_type, q_pt)
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
end

@testset "Thermodynamics - default behavior accuracy" begin
    # Input arguments should be accurate within machine precision
    # Temperature is approximated via saturation adjustment, and should be within a physical tolerance

    or(a, b) = a || b
    for ArrayType in array_types
        FT = eltype(ArrayType)
        param_set = TP.ThermodynamicsParameters(FT)
        profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
        (; T, p, e_int, ρ, θ_liq_ice, phase_type) = profiles
        (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = profiles

        RH_sat_mask = or.(RH .> 1, RH .≈ 1)
        RH_unsat_mask = .!or.(RH .> 1, RH .≈ 1)
        ts = PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot)
        @test all(saturated.(param_set, ts[RH_sat_mask]))
        @test !any(saturated.(param_set, ts[RH_unsat_mask]))

        # Test Clausius Clapeyron relation
        k = findfirst(q -> q > 0.01, q_tot) # test for one value with q_tot above some threshhold
        ts_sol = TD.PhaseEquil_ρTq(param_set, ρ[k], T[k], q_tot[k])

        function q_vap_sat(_T::FT) where {FT}
            _ρ = TD.air_density(param_set, ts_sol)
            _q_tot = TD.total_specific_humidity(param_set, ts_sol)
            _phase_type = PhaseEquil{FT}
            _q_pt = PhasePartition_equil(
                param_set,
                _T,
                oftype(_T, _ρ),
                oftype(_T, _q_tot),
                _phase_type,
            )
            return TD.q_vap_saturation(
                param_set,
                _T,
                oftype(_T, _ρ),
                _phase_type,
                _q_pt,
            )
        end

        function ∂q_vap_sat_∂T_vs_T(_T::FT) where {FT}
            _ρ = TD.air_density(param_set, ts_sol)
            _q_tot = TD.total_specific_humidity(param_set, ts_sol)
            _phase_type = PhaseEquil{FT}
            _λ = TD.liquid_fraction(param_set, ts_sol)
            _q_pt = TD.PhasePartition_equil(
                param_set,
                _T,
                oftype(_T, _ρ),
                oftype(_T, _q_tot),
                _phase_type,
            )
            _q_vap_sat =
                TD.q_vap_saturation(param_set, _T, _ρ, _phase_type, _q_pt)
            return TD.∂q_vap_sat_∂T(
                param_set,
                oftype(_T, _λ),
                _T,
                oftype(_T, _q_vap_sat),
            )
        end

        ∂q_vap_sat_∂T_fd = _T -> ForwardDiff.derivative(q_vap_sat, _T)
        @test all(
            isapprox.(
                log.(∂q_vap_sat_∂T_fd.(T)),
                log.(∂q_vap_sat_∂T_vs_T.(T));
                rtol = 2e-2,
            ),
        )

        # PhaseEquil (freezing)
        _T_freeze = FT(TP.T_freeze(param_set))
        e_int_upper =
            internal_energy_sat.(
                param_set,
                Ref(_T_freeze + sqrt(eps(FT))),
                ρ,
                q_tot,
                phase_type,
            )
        e_int_lower =
            internal_energy_sat.(
                param_set,
                Ref(_T_freeze - sqrt(eps(FT))),
                ρ,
                q_tot,
                phase_type,
            )
        _e_int = (e_int_upper .+ e_int_lower) / 2
        ts = PhaseEquil_ρeq.(param_set, ρ, _e_int, q_tot)
        @test all(
            isapprox.(
                air_temperature.(param_set, ts),
                Ref(_T_freeze),
                rtol = rtol_temperature,
            ),
        )

        # Args needs to be in sync with PhaseEquil:
        ts =
            PhaseEquil_ρeq.(
                param_set,
                ρ,
                _e_int,
                q_tot,
                8,
                FT(1e-1),
                RS.SecantMethod,
            )
        @test all(
            isapprox.(
                air_temperature.(param_set, ts),
                Ref(_T_freeze),
                rtol = rtol_temperature,
            ),
        )

        # PhaseEquil
        ts_exact = PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot, 100, FT(1e-6))
        ts = PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot)
        @test all(
            isapprox.(
                T,
                air_temperature.(param_set, ts),
                rtol = rtol_temperature,
            ),
        )

        # Should be machine accurate (because ts contains `e_int`,`ρ`,`q_tot`):
        @test all(compare_moisture.(param_set, ts, ts_exact))
        @test all(
            internal_energy.(param_set, ts) .≈
            internal_energy.(param_set, ts_exact),
        )
        @test all(
            air_density.(param_set, ts) .≈ air_density.(param_set, ts_exact),
        )
        # Approximate (temperature must be computed via saturation adjustment):
        @test all(
            isapprox.(
                air_pressure.(param_set, ts),
                air_pressure.(param_set, ts_exact),
                rtol = rtol_pressure,
            ),
        )
        @test all(
            isapprox.(
                air_temperature.(param_set, ts),
                air_temperature.(param_set, ts_exact),
                rtol = rtol_temperature,
            ),
        )

        dry_mask = abs.(q_tot .- 0) .< eps(FT)
        q_dry = q_pt[dry_mask]
        @test all(
            condensate.(q_pt) .==
            getproperty.(q_pt, :liq) .+ getproperty.(q_pt, :ice),
        )
        @test all(has_condensate.(q_dry) .== false)

        e_tot = total_energy.(param_set, ts, e_kin, e_pot)
        _cp_d = FT(TP.cp_d(param_set))
        @test all(
            specific_enthalpy.(param_set, ts) .≈
            e_int .+
            gas_constant_air.(param_set, ts) .* air_temperature.(param_set, ts),
        )
        @test all(
            total_specific_enthalpy.(param_set, ts, e_tot) .≈
            specific_enthalpy.(param_set, ts) .+ e_kin .+ e_pot,
        )
        @test all(
            moist_static_energy.(param_set, ts, e_pot) .≈
            specific_enthalpy.(param_set, ts) .+ e_pot,
        )
        @test all(
            moist_static_energy.(param_set, ts, e_pot) .≈
            total_specific_enthalpy.(param_set, ts, e_tot) .- e_kin,
        )
        @test all(
            virtual_dry_static_energy.(param_set, ts, e_pot) .≈
            _cp_d .* virtual_temperature.(param_set, ts) .+ e_pot,
        )

        # PhaseEquil
        ts_exact =
            PhaseEquil_ρeq.(
                param_set,
                ρ,
                e_int,
                q_tot,
                100,
                FT(1e-6),
                RS.SecantMethod,
            )
        ts =
            PhaseEquil_ρeq.(
                param_set,
                ρ,
                e_int,
                q_tot,
                35,
                FT(1e-4),
                RS.SecantMethod,
            ) # Needs to be in sync with default
        # Should be machine accurate (because ts contains `e_int`,`ρ`,`q_tot`):
        @test all(compare_moisture.(param_set, ts, ts_exact))
        @test all(
            internal_energy.(param_set, ts) .≈
            internal_energy.(param_set, ts_exact),
        )
        @test all(
            air_density.(param_set, ts) .≈ air_density.(param_set, ts_exact),
        )
        # Approximate (temperature must be computed via saturation adjustment):
        @test all(
            isapprox.(
                air_pressure.(param_set, ts),
                air_pressure.(param_set, ts_exact),
                rtol = rtol_pressure,
            ),
        )
        @test all(
            isapprox.(
                air_temperature.(param_set, ts),
                air_temperature.(param_set, ts_exact),
                rtol = rtol_temperature,
            ),
        )

        # PhaseEquil_ρθq
        ts_exact = PhaseEquil_ρθq.(param_set, ρ, θ_liq_ice, q_tot, 45, FT(1e-6))
        ts = PhaseEquil_ρθq.(param_set, ρ, θ_liq_ice, q_tot)
        # Should be machine accurate:
        @test all(
            air_density.(param_set, ts) .≈ air_density.(param_set, ts_exact),
        )
        @test all(compare_moisture.(param_set, ts, ts_exact))
        # Approximate (temperature must be computed via saturation adjustment):
        @test all(
            isapprox.(
                internal_energy.(param_set, ts),
                internal_energy.(param_set, ts_exact),
                atol = atol_energy,
            ),
        )
        @test all(
            isapprox.(
                liquid_ice_pottemp.(param_set, ts),
                liquid_ice_pottemp.(param_set, ts_exact),
                rtol = rtol_temperature,
            ),
        )
        @test all(
            isapprox.(
                air_temperature.(param_set, ts),
                air_temperature.(param_set, ts_exact),
                rtol = rtol_temperature,
            ),
        )

        # PhaseEquil_pθq
        ts_exact = PhaseEquil_pθq.(param_set, p, θ_liq_ice, q_tot, 40, FT(1e-6))
        ts = PhaseEquil_pθq.(param_set, p, θ_liq_ice, q_tot)

        ts =
            PhaseEquil_pθq.(
                param_set,
                p,
                θ_liq_ice,
                q_tot,
                40,
                FT(1e-4),
                RS.RegulaFalsiMethod,
            )
        # Should be machine accurate:
        @test all(compare_moisture.(param_set, ts, ts_exact))
        # Approximate (temperature must be computed via saturation adjustment):
        @test all(
            isapprox.(
                air_density.(param_set, ts),
                air_density.(param_set, ts_exact),
                rtol = rtol_density,
            ),
        )
        @test all(
            isapprox.(
                internal_energy.(param_set, ts),
                internal_energy.(param_set, ts_exact),
                atol = atol_energy,
            ),
        )
        @test all(
            isapprox.(
                liquid_ice_pottemp.(param_set, ts),
                liquid_ice_pottemp.(param_set, ts_exact),
                rtol = rtol_temperature,
            ),
        )
        @test all(
            isapprox.(
                air_temperature.(param_set, ts),
                air_temperature.(param_set, ts_exact),
                rtol = rtol_temperature,
            ),
        )

        # PhaseEquil_pθq (freezing)
        _T_freeze = FT(TP.T_freeze(param_set))

        function θ_liq_ice_closure(T, p, q_tot)
            q_pt_closure = TD.PhasePartition_equil_given_p(
                param_set,
                T,
                p,
                q_tot,
                phase_type,
            )
            θ_liq_ice_close = TD.liquid_ice_pottemp_given_pressure(
                param_set,
                T,
                p,
                q_pt_closure,
            )
            return θ_liq_ice_close
        end

        θ_liq_ice_upper =
            θ_liq_ice_closure.(Ref(_T_freeze + sqrt(eps(FT))), p, q_tot)
        θ_liq_ice_lower =
            θ_liq_ice_closure.(Ref(_T_freeze - sqrt(eps(FT))), p, q_tot)
        θ_liq_ice_mid = (θ_liq_ice_upper .+ θ_liq_ice_lower) ./ 2

        ts_lower = PhaseEquil_pθq.(param_set, p, θ_liq_ice_lower, q_tot)
        ts_upper = PhaseEquil_pθq.(param_set, p, θ_liq_ice_upper, q_tot)
        ts_mid = PhaseEquil_pθq.(param_set, p, θ_liq_ice_mid, q_tot)

        @test count(air_temperature.(param_set, ts_lower) .== Ref(_T_freeze)) ≥
              217
        @test count(air_temperature.(param_set, ts_upper) .== Ref(_T_freeze)) ≥
              217
        @test count(air_temperature.(param_set, ts_mid) .== Ref(_T_freeze)) ≥
              1296
        # we should do this instead, but we're failing because some inputs are bad
        # E.g. p ~ 110_000 Pa, q_tot ~ 0.16, which results in negative θ_liq_ice
        # This means that we should probably update our tested profiles.
        # @test all(air_temperature.(param_set, ts_lower) .== Ref(_T_freeze))
        # @test all(air_temperature.(param_set, ts_upper) .== Ref(_T_freeze))
        # @test all(air_temperature.(param_set, ts_mid) .== Ref(_T_freeze))

        # @show ρ, θ_liq_ice, q_pt
        # PhaseNonEquil_ρθq
        ts_exact =
            PhaseNonEquil_ρθq.(param_set, ρ, θ_liq_ice, q_pt, 40, FT(1e-6))
        ts = PhaseNonEquil_ρθq.(param_set, ρ, θ_liq_ice, q_pt)
        # Should be machine accurate:
        @test all(compare_moisture.(param_set, ts, ts_exact))
        @test all(
            air_density.(param_set, ts) .≈ air_density.(param_set, ts_exact),
        )
        # Approximate (temperature must be computed via non-linear solve):
        @test all(
            isapprox.(
                internal_energy.(param_set, ts),
                internal_energy.(param_set, ts_exact),
                atol = atol_energy,
            ),
        )
        @test all(
            isapprox.(
                liquid_ice_pottemp.(param_set, ts),
                liquid_ice_pottemp.(param_set, ts_exact),
                rtol = rtol_temperature,
            ),
        )
        @test all(
            isapprox.(
                air_temperature.(param_set, ts),
                air_temperature.(param_set, ts_exact),
                rtol = rtol_temperature,
            ),
        )

    end

end

@testset "Thermodynamics - exceptions on failed convergence" begin

    ArrayType = Array{Float64}
    FT = eltype(ArrayType)
    param_set = TP.ThermodynamicsParameters(FT)
    profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
    (; T, p, e_int, ρ, θ_liq_ice, phase_type) = profiles
    (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = profiles

    @test_throws ErrorException TD.saturation_adjustment.(
        RS.NewtonsMethod,
        param_set,
        e_int,
        ρ,
        q_tot,
        Ref(phase_type),
        2,
        FT(1e-10),
    )

    @test_throws ErrorException TD.saturation_adjustment.(
        RS.SecantMethod,
        param_set,
        e_int,
        ρ,
        q_tot,
        Ref(phase_type),
        2,
        FT(1e-10),
    )

    @test_throws ErrorException TD.saturation_adjustment_given_peq.(
        RS.SecantMethod,
        param_set,
        p,
        e_int,
        q_tot,
        Ref(phase_type),
        2,
        FT(1e-10),
    )

    T_virt = T # should not matter: testing for non-convergence
    @test_throws ErrorException TD.temperature_and_humidity_given_TᵥρRH.(
        param_set,
        T_virt,
        ρ,
        RH,
        Ref(phase_type),
        2,
        RS.ResidualTolerance(FT(1e-10)),
    )

    @test_throws ErrorException TD.air_temperature_given_ρθq_nonlinear.(
        param_set,
        ρ,
        θ_liq_ice,
        2,
        RS.ResidualTolerance(FT(1e-10)),
        q_pt,
    )

    @test_throws ErrorException TD.saturation_adjustment_given_ρθq.(
        param_set,
        ρ,
        θ_liq_ice,
        q_tot,
        Ref(phase_type),
        2,
        RS.ResidualTolerance(FT(1e-10)),
    )

    @test_throws ErrorException TD.saturation_adjustment_given_pθq.(
        RS.SecantMethod,
        param_set,
        p,
        θ_liq_ice,
        q_tot,
        Ref(phase_type),
        2,
        FT(1e-10),
    )

    @test_throws ErrorException TD.saturation_adjustment_given_pθq.(
        RS.NewtonsMethodAD,
        param_set,
        p,
        θ_liq_ice,
        q_tot,
        Ref(phase_type),
        2,
        FT(1e-10),
    )

    @test_throws ErrorException TD.saturation_adjustment_ρpq.(
        RS.NewtonsMethodAD,
        param_set,
        ρ,
        p,
        q_tot,
        Ref(phase_type),
        2,
        FT(1e-10),
    )

end

@testset "Thermodynamics - constructor consistency" begin

    # Make sure `ThermodynamicState` arguments are returned unchanged

    for ArrayType in array_types
        FT = eltype(ArrayType)
        param_set = TP.ThermodynamicsParameters(FT)

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
                FT(1e-1),
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
        # TODO: should this pass?
        # @test all(isapprox.(internal_energy.(param_set, ts_phq), e_int, atol = atol_energy)) # fails
        # @show maximum(abs.(internal_energy.(param_set, ts_phq) .- e_int)) # ~2258.34
        # @show maximum(abs.(internal_energy.(param_set, ts_phq) .- e_int) ./ e_int * 100) # ~0.502
        @test all(
            isapprox.(
                internal_energy.(param_set, ts_phq),
                e_int,
                rtol = rtol_energy,
            ),
        )
        @test all(
            getproperty.(PhasePartition.(param_set, ts_phq), :tot) .≈ q_tot,
        )
        @test all(air_pressure.(param_set, ts_phq) .≈ p)

        ts_pθq = PhaseEquil_pθq.(param_set, p, θ_liq_ice, q_tot)
        @test all(air_pressure.(param_set, ts_pθq) .≈ p)
        # TODO: Run some tests to make sure that this decreses with
        # decreasing temperature_tol (and increasing maxiter)
        # @show maximum(abs.(liquid_ice_pottemp.(param_set, ts_pθq) .- θ_liq_ice))
        @test all(
            isapprox.(
                liquid_ice_pottemp.(param_set, ts_pθq),
                θ_liq_ice,
                rtol = rtol_temperature,
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
        @test all(isapprox.(T_non_linear, T_expansion, rtol = rtol_temperature))
        e_int_ = internal_energy.(param_set, T_non_linear, q_pt)
        ts = PhaseNonEquil.(param_set, e_int_, ρ, q_pt)
        @test all(T_non_linear .≈ air_temperature.(param_set, ts))
        @test all(
            isapprox(
                θ_liq_ice,
                liquid_ice_pottemp.(param_set, ts),
                rtol = rtol_temperature,
            ),
        )

        # PhaseEquil_ρθq
        ts = PhaseEquil_ρθq.(param_set, ρ, θ_liq_ice, q_tot, 45, FT(1e-3))
        @test all(
            isapprox.(
                liquid_ice_pottemp.(param_set, ts),
                θ_liq_ice,
                rtol = rtol_temperature,
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
        ts = PhaseEquil_pθq.(param_set, p, θ_liq_ice, q_tot, 35, FT(1e-3))
        @test all(
            isapprox.(
                liquid_ice_pottemp.(param_set, ts),
                θ_liq_ice,
                rtol = rtol_temperature,
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
        ts = PhaseNonEquil_ρθq.(param_set, ρ, θ_liq_ice, q_pt, 5, FT(1e-3))
        @test all(
            isapprox.(
                θ_liq_ice,
                liquid_ice_pottemp.(param_set, ts),
                rtol = rtol_temperature,
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


@testset "Thermodynamics - type-stability" begin

    # NOTE: `Float32` saturation adjustment tends to have more difficulty
    # with converging to the same tolerances as `Float64`, so they're relaxed here.
    ArrayType = Array{Float32}
    FT = eltype(ArrayType)
    param_set = TP.ThermodynamicsParameters(FT)

    profiles = TestedProfiles.PhaseDryProfiles(param_set, ArrayType)
    (; T, p, e_int, ρ, θ_liq_ice, phase_type) = profiles
    (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = profiles

    θ_dry = dry_pottemp.(param_set, T, ρ)
    ts_dry = PhaseDry.(param_set, e_int, ρ)
    ts_dry_ρp = PhaseDry_ρp.(param_set, ρ, p)
    ts_dry_pT = PhaseDry_pT.(param_set, p, T)
    ts_dry_ρθ = PhaseDry_ρθ.(param_set, ρ, θ_dry)
    ts_dry_pθ = PhaseDry_pθ.(param_set, p, θ_dry)

    profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
    (; T, p, e_int, ρ, θ_liq_ice, phase_type) = profiles
    (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = profiles

    ρu = FT[1.0, 2.0, 3.0]
    @test typeof.(internal_energy.(ρ, ρ .* e_int, Ref(ρu), e_pot)) ==
          typeof.(e_int)

    ts_eq = PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot, 15, FT(1e-4))
    e_tot = total_energy.(param_set, ts_eq, e_kin, e_pot)

    ts_T =
        PhaseEquil_ρTq.(
            param_set,
            air_density.(param_set, ts_eq),
            air_temperature.(param_set, ts_eq),
            q_tot,
        )
    ts_Tp =
        PhaseEquil_pTq.(
            param_set,
            air_pressure.(param_set, ts_eq),
            air_temperature.(param_set, ts_eq),
            q_tot,
        )
    ts_ρp =
        PhaseEquil_ρpq.(
            param_set,
            air_density.(param_set, ts_eq),
            air_pressure.(param_set, ts_eq),
            q_tot,
        )

    @test all(
        air_temperature.(param_set, ts_T) .≈ air_temperature.(param_set, ts_Tp),
    )
    @test all(air_pressure.(param_set, ts_T) .≈ air_pressure.(param_set, ts_Tp))
    @test all(
        total_specific_humidity.(param_set, ts_T) .≈
        total_specific_humidity.(param_set, ts_Tp),
    )

    ts_neq = PhaseNonEquil.(param_set, e_int, ρ, q_pt)
    ts_ρT_neq = PhaseNonEquil_ρTq.(param_set, ρ, T, q_pt)
    ts_pT_neq = PhaseNonEquil_pTq.(param_set, p, T, q_pt)

    ts_θ_liq_ice_eq =
        PhaseEquil_ρθq.(param_set, ρ, θ_liq_ice, q_tot, 45, FT(1e-4))
    ts_θ_liq_ice_eq_p =
        PhaseEquil_pθq.(param_set, p, θ_liq_ice, q_tot, 40, FT(1e-4))
    ts_θ_liq_ice_neq = PhaseNonEquil_ρθq.(param_set, ρ, θ_liq_ice, q_pt)
    ts_θ_liq_ice_neq_p = PhaseNonEquil_pθq.(param_set, p, θ_liq_ice, q_pt)

    for ts in (
        ts_dry,
        ts_dry_ρp,
        ts_dry_pT,
        ts_dry_ρθ,
        ts_dry_pθ,
        ts_eq,
        ts_T,
        ts_Tp,
        ts_ρp,
        ts_neq,
        ts_ρT_neq,
        ts_pT_neq,
        ts_θ_liq_ice_eq,
        ts_θ_liq_ice_eq_p,
        ts_θ_liq_ice_neq,
        ts_θ_liq_ice_neq_p,
    )
        @test typeof.(soundspeed_air.(param_set, ts)) == typeof.(e_int)
        @test typeof.(gas_constant_air.(param_set, ts)) == typeof.(e_int)
        @test typeof.(specific_enthalpy.(param_set, ts)) == typeof.(e_int)
        @test typeof.(vapor_specific_humidity.(param_set, ts)) == typeof.(e_int)
        @test typeof.(relative_humidity.(param_set, ts)) == typeof.(e_int)
        @test typeof.(air_pressure.(param_set, ts)) == typeof.(e_int)
        @test typeof.(air_density.(param_set, ts)) == typeof.(e_int)
        @test typeof.(total_specific_humidity.(param_set, ts)) == typeof.(e_int)
        @test typeof.(liquid_specific_humidity.(param_set, ts)) ==
              typeof.(e_int)
        @test typeof.(ice_specific_humidity.(param_set, ts)) == typeof.(e_int)
        @test typeof.(cp_m.(param_set, ts)) == typeof.(e_int)
        @test typeof.(cv_m.(param_set, ts)) == typeof.(e_int)
        @test typeof.(air_temperature.(param_set, ts)) == typeof.(e_int)
        @test typeof.(internal_energy_sat.(param_set, ts)) == typeof.(e_int)
        @test typeof.(internal_energy.(param_set, ts)) == typeof.(e_int)
        @test typeof.(internal_energy_dry.(param_set, ts)) == typeof.(e_int)
        @test typeof.(internal_energy_vapor.(param_set, ts)) == typeof.(e_int)
        @test typeof.(internal_energy_liquid.(param_set, ts)) == typeof.(e_int)
        @test typeof.(internal_energy_ice.(param_set, ts)) == typeof.(e_int)
        @test typeof.(latent_heat_vapor.(param_set, ts)) == typeof.(e_int)
        @test typeof.(latent_heat_sublim.(param_set, ts)) == typeof.(e_int)
        @test typeof.(latent_heat_fusion.(param_set, ts)) == typeof.(e_int)
        @test typeof.(q_vap_saturation.(param_set, ts)) == typeof.(e_int)
        @test typeof.(q_vap_saturation_liquid.(param_set, ts)) == typeof.(e_int)
        @test typeof.(q_vap_saturation_ice.(param_set, ts)) == typeof.(e_int)
        @test typeof.(saturation_excess.(param_set, ts)) == typeof.(e_int)
        @test typeof.(liquid_fraction.(param_set, ts)) == typeof.(e_int)
        @test typeof.(liquid_ice_pottemp.(param_set, ts)) == typeof.(e_int)
        @test typeof.(dry_pottemp.(param_set, ts)) == typeof.(e_int)
        @test typeof.(exner.(param_set, ts)) == typeof.(e_int)
        @test typeof.(liquid_ice_pottemp_sat.(param_set, ts)) == typeof.(e_int)
        @test typeof.(specific_volume.(param_set, ts)) == typeof.(e_int)
        @test typeof.(supersaturation.(param_set, ts, Ice())) == typeof.(e_int)
        @test typeof.(supersaturation.(param_set, ts, Liquid())) ==
              typeof.(e_int)
        @test typeof.(virtual_pottemp.(param_set, ts)) == typeof.(e_int)
        @test typeof.(specific_entropy.(param_set, ts)) == typeof.(e_int)
        @test eltype.(gas_constants.(param_set, ts)) == typeof.(e_int)

        @test typeof.(total_specific_enthalpy.(param_set, ts, e_tot)) ==
              typeof.(e_int)
        @test typeof.(moist_static_energy.(param_set, ts, e_pot)) ==
              typeof.(e_int)
        @test typeof.(getproperty.(PhasePartition.(param_set, ts), :tot)) ==
              typeof.(e_int)
        @test typeof.(virtual_dry_static_energy.(param_set, ts, e_pot)) ==
              typeof.(e_int)
    end
    @test typeof(
        q_vap_from_RH_liquid(param_set, FT(100000), FT(275), FT(0.8)),
    ) == FT
end

@testset "Thermodynamics - dry limit" begin

    ArrayType = Array{Float64}
    FT = eltype(ArrayType)
    param_set = TP.ThermodynamicsParameters(FT)
    profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
    (; T, p, e_int, ρ, θ_liq_ice, phase_type) = profiles
    (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = profiles

    # PhasePartition test is noisy, so do this only once:
    ts_dry = PhaseDry(param_set, first(e_int), first(ρ))
    ts_eq =
        PhaseEquil_ρeq(param_set, first(ρ), first(e_int), typeof(first(ρ))(0))
    @test PhasePartition(param_set, ts_eq).tot ≈
          PhasePartition(param_set, ts_dry).tot
    @test PhasePartition(param_set, ts_eq).liq ≈
          PhasePartition(param_set, ts_dry).liq
    @test PhasePartition(param_set, ts_eq).ice ≈
          PhasePartition(param_set, ts_dry).ice

    @test mixing_ratios(param_set, ts_eq).tot ≈
          mixing_ratios(param_set, ts_dry).tot
    @test mixing_ratios(param_set, ts_eq).liq ≈
          mixing_ratios(param_set, ts_dry).liq
    @test mixing_ratios(param_set, ts_eq).ice ≈
          mixing_ratios(param_set, ts_dry).ice

    ts_dry = PhaseDry.(param_set, e_int, ρ)
    ts_eq = PhaseEquil_ρeq.(param_set, ρ, e_int, q_tot .* 0)

    @test all(
        gas_constant_air.(param_set, ts_eq) .≈
        gas_constant_air.(param_set, ts_dry),
    )
    @test all(
        relative_humidity.(param_set, ts_eq) .≈
        relative_humidity.(param_set, ts_dry),
    )
    @test all(
        air_pressure.(param_set, ts_eq) .≈ air_pressure.(param_set, ts_dry),
    )
    @test all(air_density.(param_set, ts_eq) .≈ air_density.(param_set, ts_dry))
    @test all(
        specific_volume.(param_set, ts_eq) .≈
        specific_volume.(param_set, ts_dry),
    )
    @test all(
        total_specific_humidity.(param_set, ts_eq) .≈
        total_specific_humidity.(param_set, ts_dry),
    )
    @test all(
        liquid_specific_humidity.(param_set, ts_eq) .≈
        liquid_specific_humidity.(param_set, ts_dry),
    )
    @test all(
        ice_specific_humidity.(param_set, ts_eq) .≈
        ice_specific_humidity.(param_set, ts_dry),
    )
    @test all(cp_m.(param_set, ts_eq) .≈ cp_m.(param_set, ts_dry))
    @test all(cv_m.(param_set, ts_eq) .≈ cv_m.(param_set, ts_dry))
    @test all(
        air_temperature.(param_set, ts_eq) .≈
        air_temperature.(param_set, ts_dry),
    )
    @test all(
        internal_energy.(param_set, ts_eq) .≈
        internal_energy.(param_set, ts_dry),
    )
    @test all(
        internal_energy_sat.(param_set, ts_eq) .≈
        internal_energy_sat.(param_set, ts_dry),
    )
    @test all(
        internal_energy_dry.(param_set, ts_eq) .≈
        internal_energy_dry.(param_set, ts_dry),
    )
    @test all(
        internal_energy_vapor.(param_set, ts_eq) .≈
        internal_energy_vapor.(param_set, ts_dry),
    )
    @test all(
        internal_energy_liquid.(param_set, ts_eq) .≈
        internal_energy_liquid.(param_set, ts_dry),
    )
    @test all(
        internal_energy_ice.(param_set, ts_eq) .≈
        internal_energy_ice.(param_set, ts_dry),
    )
    @test all(
        soundspeed_air.(param_set, ts_eq) .≈ soundspeed_air.(param_set, ts_dry),
    )
    @test all(
        supersaturation.(param_set, ts_eq, Ice()) .≈
        supersaturation.(param_set, ts_dry, Ice()),
    )
    @test all(
        supersaturation.(param_set, ts_eq, Liquid()) .≈
        supersaturation.(param_set, ts_dry, Liquid()),
    )
    @test all(
        latent_heat_vapor.(param_set, ts_eq) .≈
        latent_heat_vapor.(param_set, ts_dry),
    )
    @test all(
        latent_heat_sublim.(param_set, ts_eq) .≈
        latent_heat_sublim.(param_set, ts_dry),
    )
    @test all(
        latent_heat_fusion.(param_set, ts_eq) .≈
        latent_heat_fusion.(param_set, ts_dry),
    )
    @test all(
        q_vap_saturation.(param_set, ts_eq) .≈
        q_vap_saturation.(param_set, ts_dry),
    )
    @test all(
        q_vap_saturation_liquid.(param_set, ts_eq) .≈
        q_vap_saturation_liquid.(param_set, ts_dry),
    )
    @test all(
        q_vap_saturation_ice.(param_set, ts_eq) .≈
        q_vap_saturation_ice.(param_set, ts_dry),
    )
    @test all(
        saturation_excess.(param_set, ts_eq) .≈
        saturation_excess.(param_set, ts_dry),
    )
    @test all(
        liquid_fraction.(param_set, ts_eq) .≈
        liquid_fraction.(param_set, ts_dry),
    )
    @test all(
        liquid_ice_pottemp.(param_set, ts_eq) .≈
        liquid_ice_pottemp.(param_set, ts_dry),
    )
    @test all(dry_pottemp.(param_set, ts_eq) .≈ dry_pottemp.(param_set, ts_dry))
    @test all(
        virtual_pottemp.(param_set, ts_eq) .≈
        virtual_pottemp.(param_set, ts_dry),
    )
    @test all(
        specific_entropy.(param_set, ts_eq) .≈
        specific_entropy.(param_set, ts_dry),
    )
    @test all(
        liquid_ice_pottemp_sat.(param_set, ts_eq) .≈
        liquid_ice_pottemp_sat.(param_set, ts_dry),
    )
    @test all(exner.(param_set, ts_eq) .≈ exner.(param_set, ts_dry))

    @test all(
        saturation_vapor_pressure.(param_set, ts_eq, Ice()) .≈
        saturation_vapor_pressure.(param_set, ts_dry, Ice()),
    )
    @test all(
        saturation_vapor_pressure.(param_set, ts_eq, Liquid()) .≈
        saturation_vapor_pressure.(param_set, ts_dry, Liquid()),
    )
    @test all(
        first.(gas_constants.(param_set, ts_eq)) ≈
        first.(gas_constants.(param_set, ts_dry)),
    )
    @test all(
        last.(gas_constants.(param_set, ts_eq)) ≈
        last.(gas_constants.(param_set, ts_dry)),
    )

end

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

@testset "Thermodynamics - test T_guess" begin
    ArrayType = Array{Float64}
    FT = eltype(ArrayType)
    param_set = TP.ThermodynamicsParameters(FT)
    profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
    (; p, ρ, e_int, h, θ_liq_ice, q_tot, T, phase_type) = profiles
    T_guess = T .+ (FT(0.2) .* randn(FT, length(T)))
    args = (q_tot, 40, FT(1e-1))
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
            FT(1e-4),
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
@testset "Test data collection" begin
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
