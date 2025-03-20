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

array_types = [Array{Float32}, Array{Float64}]

include("data_tests.jl")

compare_moisture(param_set, a::ThermodynamicState, b::ThermodynamicState) =
    compare_moisture(param_set, a, PhasePartition(param_set, b))

compare_moisture(param_set, ts::PhaseEquil, q_pt::PhasePartition) =
    getproperty(PhasePartition(param_set, ts), :tot) ≈ getproperty(q_pt, :tot)

compare_moisture(param_set, ts::PhaseNonEquil, q_pt::PhasePartition) = all((
    getproperty(PhasePartition(param_set, ts), :tot) ≈ getproperty(q_pt, :tot),
    getproperty(PhasePartition(param_set, ts), :liq) ≈ getproperty(q_pt, :liq),
    getproperty(PhasePartition(param_set, ts), :ice) ≈ getproperty(q_pt, :ice),
))

@testset "Thermodynamics - isentropic processes" begin
    for ArrayType in array_types
        FT = eltype(ArrayType)
        param_set = TP.ThermodynamicsParameters(FT)

        _R_d = FT(TP.R_d(param_set))
        _molmass_ratio = FT(TP.molmass_ratio(param_set))
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

        profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
        (; T, p, e_int, ρ, θ_liq_ice, phase_type) = profiles
        (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = profiles

        # Test state for thermodynamic consistency (with ideal gas law)
        T_idgl = TD.air_temperature_from_ideal_gas_law.(param_set, p, ρ, q_pt)
        @test all(T .≈ T_idgl)

        Φ = FT(1)
        Random.seed!(15)
        perturbation = FT(0.1) * rand(FT, length(T))

        # TODO: Use reasonable values for ambient temperature/pressure
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
    _molmass_ratio = FT(TP.molmass_ratio(param_set))
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

    # ideal gas law
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

    # gas constants and heat capacities
    @test gas_constant_air(param_set, PhasePartition(FT(0))) === _R_d
    @test gas_constant_air(param_set, PhasePartition(FT(1))) === _R_v
    @test gas_constant_air(param_set, PhasePartition(FT(0.5), FT(0.5))) ≈
          _R_d / 2
    @test gas_constant_air(param_set, FT) == _R_d

    @test cp_m(param_set, PhasePartition(FT(0))) === _cp_d
    @test cp_m(param_set, PhasePartition(FT(1))) === _cp_v
    @test cp_m(param_set, PhasePartition(FT(1), FT(1))) === _cp_l
    @test cp_m(param_set, PhasePartition(FT(1), FT(0), FT(1))) === _cp_i
    @test cp_m(param_set, FT) == _cp_d

    @test cv_m(param_set, PhasePartition(FT(0))) === _cp_d - _R_d
    @test cv_m(param_set, PhasePartition(FT(1))) === _cp_v - _R_v
    @test cv_m(param_set, PhasePartition(FT(1), FT(1))) === _cv_l
    @test cv_m(param_set, PhasePartition(FT(1), FT(0), FT(1))) === _cv_i
    @test cv_m(param_set, FT) == _cv_d

    # speed of sound
    @test soundspeed_air(param_set, _T_0 + 20, PhasePartition(FT(0))) ==
          sqrt(_cp_d / _cv_d * _R_d * (_T_0 + 20))
    @test soundspeed_air(param_set, _T_0 + 100, PhasePartition(FT(1))) ==
          sqrt(_cp_v / _cv_v * _R_v * (_T_0 + 100))

    # specific latent heats
    @test latent_heat_vapor(param_set, _T_0) ≈ _LH_v0
    @test latent_heat_fusion(param_set, _T_0) ≈ _LH_f0
    @test latent_heat_sublim(param_set, _T_0) ≈ _LH_s0

    # saturation vapor pressure and specific humidity
    p = FT(1.e5)
    q_tot = FT(0.23)
    ρ = FT(1.0)
    ρ_v_triple = _press_triple / _R_v / _T_triple
    p_v_tripleminus = FT(0.99) * _press_triple
    @test saturation_vapor_pressure(param_set, _T_triple, Liquid()) ≈
          _press_triple
    @test saturation_vapor_pressure(param_set, _T_triple, Ice()) ≈ _press_triple

    phase_type = PhaseDry
    @test q_vap_saturation(
        param_set,
        _T_triple,
        ρ,
        phase_type,
        PhasePartition(FT(0)),
    ) == ρ_v_triple / ρ

    @test TD.q_vap_saturation_from_pressure(
        param_set,
        q_tot,
        p,
        _T_triple,
        phase_type,
    ) == _R_d / _R_v * (1 - q_tot) * _press_triple / (p - _press_triple)

    @test TD.q_vap_saturation_from_pressure(
        param_set,
        q_tot,
        p_v_tripleminus,
        _T_triple,
        phase_type,
    ) == FT(1)

    phase_type = PhaseNonEquil
    @test q_vap_saturation(
        param_set,
        _T_triple,
        ρ,
        phase_type,
        PhasePartition(q_tot, q_tot),
    ) == ρ_v_triple / ρ

    @test q_vap_saturation_generic(param_set, _T_triple, ρ, Liquid()) ==
          ρ_v_triple / ρ
    @test q_vap_saturation_generic(param_set, _T_triple, ρ, Ice()) ==
          ρ_v_triple / ρ
    @test q_vap_saturation_generic(param_set, _T_triple - 20, ρ, Liquid()) >=
          q_vap_saturation_generic(param_set, _T_triple - 20, ρ, Ice())

    # test the wrapper for q_vap_saturation over liquid water and ice
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

    phase_type = PhaseDry
    @test saturation_excess(
        param_set,
        _T_triple,
        ρ,
        phase_type,
        PhasePartition(q_tot),
    ) == q_tot - ρ_v_triple / ρ
    @test saturation_excess(
        param_set,
        _T_triple,
        ρ,
        phase_type,
        PhasePartition(q_tot / 1000),
    ) == 0.0

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

    # energy functions and inverse (temperature)
    T = FT(300)
    e_kin = FT(11)
    e_pot = FT(13)
    @test air_temperature(param_set, _cv_d * (T - _T_0) - _R_d * _T_0) === FT(T)
    @test air_temperature(
        param_set,
        _cv_d * (T - _T_0) - _R_d * _T_0,
        PhasePartition(FT(0)),
    ) === FT(T)

    @test air_temperature(
        param_set,
        cv_m(param_set, PhasePartition(FT(0))) * (T - _T_0) - _R_d * _T_0,
        PhasePartition(FT(0)),
    ) === FT(T)
    @test air_temperature(
        param_set,
        cv_m(param_set, PhasePartition(FT(q_tot))) * (T - _T_0) -
        (1 - q_tot) * _R_d * _T_0 + q_tot * _e_int_v0,
        PhasePartition(q_tot),
    ) ≈ FT(T)

    @test total_energy(param_set, FT(e_kin), FT(e_pot), _T_0) ===
          FT(e_kin) + FT(e_pot) - _R_d * _T_0
    @test total_energy(param_set, FT(e_kin), FT(e_pot), FT(T)) ≈
          FT(e_kin) + FT(e_pot) + _cv_d * (T - _T_0) - _R_d * _T_0
    @test total_energy(param_set, FT(0), FT(0), _T_0, PhasePartition(q_tot)) ≈
          q_tot * _e_int_v0 - (1 - q_tot) * _R_d * _T_0

    # phase partitioning in equilibrium
    q_liq = FT(0.1)
    T = FT(_T_icenuc - 10)
    ρ = FT(1.0)
    q_tot = FT(0.21)
    phase_type = PhaseDry
    @test liquid_fraction(param_set, T, phase_type) === FT(0)
    phase_type = PhaseNonEquil
    @test liquid_fraction(
        param_set,
        T,
        phase_type,
        PhasePartition(q_tot, q_liq, q_liq),
    ) === FT(0.5)
    phase_type = PhaseDry
    q = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    @test q.liq ≈ FT(0)
    @test 0 < q.ice <= q_tot

    q = TD.PhasePartition_equil_given_p(param_set, T, p, q_tot, phase_type)
    @test q.liq ≈ FT(0)
    @test 0 < q.ice <= q_tot

    T = FT(_T_freeze + 10)
    ρ = FT(0.1)
    q_tot = FT(0.60)
    @test liquid_fraction(param_set, T, phase_type) === FT(1)
    phase_type = PhaseNonEquil
    @test liquid_fraction(
        param_set,
        T,
        phase_type,
        PhasePartition(q_tot, q_liq, q_liq / 2),
    ) === FT(2 / 3)
    phase_type = PhaseDry
    q = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    @test 0 < q.liq <= q_tot
    @test q.ice ≈ 0

    q = TD.PhasePartition_equil_given_p(param_set, T, p, q_tot, phase_type)
    @test 0 < q.liq <= q_tot
    @test q.ice ≈ 0

    # saturation adjustment in equilibrium (i.e., given the thermodynamic
    # variables E_int, p, q_tot, compute the temperature and partitioning of the phases
    q_tot = FT(0)
    ρ = FT(1)
    phase_type = PhaseEquil
    @test TD.saturation_adjustment(
        RS.SecantMethod,
        param_set,
        internal_energy_sat(param_set, 300.0, ρ, q_tot, phase_type),
        ρ,
        q_tot,
        phase_type,
        10,
        rtol_temperature,
    ) ≈ 300.0
    @test abs(
        TD.saturation_adjustment(
            RS.NewtonsMethod,
            param_set,
            internal_energy_sat(param_set, 300.0, ρ, q_tot, phase_type),
            ρ,
            q_tot,
            phase_type,
            10,
            rtol_temperature,
        ) - 300.0,
    ) < rtol_temperature

    q_tot = FT(0.21)
    ρ = FT(0.1)
    @test isapprox(
        TD.saturation_adjustment(
            RS.SecantMethod,
            param_set,
            internal_energy_sat(param_set, 200.0, ρ, q_tot, phase_type),
            ρ,
            q_tot,
            phase_type,
            10,
            rtol_temperature,
        ),
        200.0,
        rtol = rtol_temperature,
    )
    @test abs(
        TD.saturation_adjustment(
            RS.NewtonsMethod,
            param_set,
            internal_energy_sat(param_set, 200.0, ρ, q_tot, phase_type),
            ρ,
            q_tot,
            phase_type,
            10,
            rtol_temperature,
        ) - 200.0,
    ) < rtol_temperature
    q = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    @test q.tot - q.liq - q.ice ≈
          vapor_specific_humidity(q) ≈
          q_vap_saturation(param_set, T, ρ, phase_type)

    q = TD.PhasePartition_equil_given_p(param_set, T, p, q_tot, phase_type)
    @test q.tot - q.liq - q.ice ≈
          vapor_specific_humidity(q) ≈
          TD.q_vap_saturation_from_pressure(param_set, q_tot, p, T, phase_type)

    ρ = FT(1)
    ρu = FT[1, 2, 3]
    ρe = FT(1100)
    e_pot = FT(93)
    @test internal_energy(ρ, ρe, ρu, e_pot) ≈ 1000.0

    # internal energies for dry, vapor, liquid and ice
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

    # potential temperatures
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

    # dry potential temperatures. FIXME: add correctness tests
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

    # Exner function. FIXME: add correctness tests
    p = FT(1.e5)
    q_tot = FT(0.23)
    @test TD.exner_given_pressure(param_set, p, PhasePartition(q_tot)) isa
          typeof(p)

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
    @test vmrs ≈ _molmass_ratio * shum_to_mixing_ratio(q_vap, q.tot)

    # Sanity check for PhaseEquil_pTRH constructor
    # No humidity
    p = FT(1.e5); T=FT(300.) ; RH=FT(0.)
    ts_pTRH = PhaseEquil_pTRH(param_set, p, T, RH)
    q_tot_expected = FT(0)
    @test q_tot_expected ≈ TD.total_specific_humidity(param_set, ts_pTRH)

    # q at Saturation
    p = FT(1000); T=300 ; RH=1
    ts_pTRH = PhaseEquil_pTRH(param_set, p, T, RH)
    saturation_vapor_pressure = saturation_vapor_pressure(param_set, T, Liquid())
    q_tot_expected = FT(0.05559498223324131)
    @test q_tot_expected ≈ TD.total_specific_humidity(param_set, ts_pTRH)
    
    # Relative humidity sanity checks
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

        # test Clausius Clapeyron relation
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

        # We can't pass on this yet, due to https://github.com/CliMA/ClimateMachine.jl/issues/263
        # ts_pθq = PhaseEquil_pθq.(param_set, p, θ_liq_ice, q_tot, 40, FT(1e-3), RS.NewtonsMethodAD)
        # @test all(air_pressure.(param_set, ts_pθq) .≈ p)
        # @test all(isapprox.(
        #     liquid_ice_pottemp.(param_set, ts_pθq),
        #     θ_liq_ice,
        #     rtol = rtol_temperature,
        # ))
        # @test all(getproperty.(PhasePartition.(param_set, ts_pθq), :tot) .≈ q_tot)

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
        # passes the consistency test within sufficient physical precision,
        # however, it fails to satisfy the consistency test within machine
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

        # TODO: Add this test back in
        @test all(RH_sat .≈ 1)

        # Test that RH is zero for dry conditions
        q_pt_dry = PhasePartition.(zeros(FT, length(T)))
        p_dry = air_pressure.(param_set, T, ρ, q_pt_dry)
        RH_dry =
            relative_humidity.(param_set, T, p_dry, Ref(phase_type), q_pt_dry)
        @test all(RH_dry .≈ 0)


        # Test virtual temperature and inverse functions:
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
