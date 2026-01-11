"""
Convergence + correctness tests for the functional `saturation_adjustment` API.

This suite has two parts:
1) **Correctness**: Uses equilibrium thermodynamic profiles from `TestedProfiles.jl`
   (which provide realistic, internally consistent reference values) and checks that
   each IndepVars formulation converges and recovers the reference temperature and
   phase partitioning.
2) **Non-redundant coverage**: A small **smoke** convergence test that constructs
   inputs from additional `TemperatureProfiles` implementations (`DryAdiabaticProfile`,
   `IsothermalProfile`) to ensure those profile types do not trigger solver pathologies.

Additionally, for `ρeq` we time multiple RootSolvers methods (with default hyperparameters)
and log elapsed time (no hard timing thresholds).
"""

const TDTP_SA = TD.TemperatureProfiles

@testset "Thermodynamics - saturation_adjustment convergence+correctness" begin
    for FT in (Float32, Float64)
        param_set = FT == Float64 ? param_set_Float64 : param_set_Float32

        maxiter = 40
        tol = FT(rtol_temperature)

        approx_tight(a, b) =
            isapprox(a, b; rtol = FT(100) * eps(FT), atol = FT(100) * eps(FT))

        function check_partition(Tsol, ρ_eff, q0, ql, qi)
            (ql_exp, qi_exp) = TD.condensate_partition(param_set, Tsol, ρ_eff, q0)
            @test approx_tight(ql, ql_exp)
            @test approx_tight(qi, qi_exp)
            return nothing
        end

        @testset "TestedProfiles equilibrium columns ($FT)" begin
            profiles = TestedProfiles.PhaseEquilProfiles(param_set, Array{FT})
            (; T, p, ρ, q_tot) = profiles

            # Sample ~120 points across the full profile grid.
            idxs = unique(round.(Int, range(1, length(T), length = 250)))

            # For each sampled point, compute reference targets for each formulation.
            function targets(i::Int)
                T0 = T[i]
                p0 = p[i]
                ρ0 = ρ[i]
                q0 = q_tot[i]

                # ρ-based targets
                (q_liq_ρ, q_ice_ρ) = TD.condensate_partition(param_set, T0, ρ0, q0)
                e_int_ρ = TD.internal_energy_sat(param_set, T0, ρ0, q0)
                p_ρ = TD.air_pressure(param_set, T0, ρ0, q0, q_liq_ρ, q_ice_ρ)
                θ_ρ = TD.liquid_ice_pottemp(param_set, T0, ρ0, q0, q_liq_ρ, q_ice_ρ)

                # p-based targets (matches pe/ph/pθ_li internals: ρ(T) computed from (p,T,q_tot))
                ρ_p = TD.air_density(param_set, T0, p0, q0)
                (q_liq_p, q_ice_p) = TD.condensate_partition(param_set, T0, ρ_p, q0)
                e_int_p = TD.internal_energy_sat(param_set, T0, ρ_p, q0)
                h_p = TD.enthalpy_sat(param_set, T0, ρ_p, q0)
                θ_p =
                    TD.liquid_ice_pottemp_given_pressure(
                        param_set,
                        T0,
                        p0,
                        q0,
                        q_liq_p,
                        q_ice_p,
                    )

                return (;
                    T0,
                    p0,
                    ρ0,
                    q0,
                    e_int_ρ,
                    p_ρ,
                    θ_ρ,
                    e_int_p,
                    h_p,
                    θ_p,
                )
            end

            @testset "All IndepVars recover reference state ($FT)" begin
                for i in idxs
                    inp = targets(i)

                    # ρeq
                    let (; T, q_liq, q_ice, converged) = TD.saturation_adjustment(
                            RS.SecantMethod,
                            param_set,
                            TD.ρe(),
                            inp.ρ0,
                            inp.e_int_ρ,
                            inp.q0,
                            maxiter,
                            tol,
                        )
                        @test isfinite(T)
                        @test converged
                        @test isapprox(
                            T,
                            inp.T0;
                            atol = FT(atol_temperature),
                            rtol = FT(0),
                        )
                        check_partition(T, inp.ρ0, inp.q0, q_liq, q_ice)
                    end

                    # peq
                    let (; T, q_liq, q_ice, converged) = TD.saturation_adjustment(
                            RS.SecantMethod,
                            param_set,
                            TD.pe(),
                            inp.p0,
                            inp.e_int_p,
                            inp.q0,
                            maxiter,
                            tol,
                        )
                        @test isfinite(T)
                        @test converged
                        @test isapprox(
                            T,
                            inp.T0;
                            atol = FT(atol_temperature),
                            rtol = FT(0),
                        )
                        ρ_eff = TD.air_density(param_set, T, inp.p0, inp.q0)
                        check_partition(T, ρ_eff, inp.q0, q_liq, q_ice)
                    end

                    # phq
                    let (; T, q_liq, q_ice, converged) = TD.saturation_adjustment(
                            RS.SecantMethod,
                            param_set,
                            TD.ph(),
                            inp.p0,
                            inp.h_p,
                            inp.q0,
                            maxiter,
                            tol,
                        )
                        @test isfinite(T)
                        @test converged
                        @test isapprox(
                            T,
                            inp.T0;
                            atol = FT(atol_temperature),
                            rtol = FT(0),
                        )
                        ρ_eff = TD.air_density(param_set, T, inp.p0, inp.q0)
                        check_partition(T, ρ_eff, inp.q0, q_liq, q_ice)
                    end

                    # pρ (uses ρ-based equilibrium; p target is p_ρ)
                    let (; T, q_liq, q_ice, converged) = TD.saturation_adjustment(
                            RS.SecantMethod,
                            param_set,
                            TD.pρ(),
                            inp.p_ρ,
                            inp.ρ0,
                            inp.q0,
                            maxiter,
                            tol,
                        )
                        @test isfinite(T)
                        @test converged
                        @test isapprox(
                            T,
                            inp.T0;
                            atol = FT(atol_temperature),
                            rtol = FT(0),
                        )
                        check_partition(T, inp.ρ0, inp.q0, q_liq, q_ice)
                    end

                    # pθ_li (p-based)
                    let (; T, q_liq, q_ice, converged) = TD.saturation_adjustment(
                            RS.SecantMethod,
                            param_set,
                            TD.pθ_li(),
                            inp.p0,
                            inp.θ_p,
                            inp.q0,
                            maxiter,
                            tol,
                        )
                        @test isfinite(T)
                        @test converged
                        @test isapprox(
                            T,
                            inp.T0;
                            atol = FT(atol_temperature),
                            rtol = FT(0),
                        )
                        ρ_eff = TD.air_density(param_set, T, inp.p0, inp.q0)
                        check_partition(T, ρ_eff, inp.q0, q_liq, q_ice)
                    end

                    # ρθ_li (ρ-based)
                    let (; T, q_liq, q_ice, converged) = TD.saturation_adjustment(
                            RS.SecantMethod,
                            param_set,
                            TD.ρθ_li(),
                            inp.ρ0,
                            inp.θ_ρ,
                            inp.q0,
                            maxiter,
                            tol,
                        )
                        @test isfinite(T)
                        @test converged
                        @test isapprox(
                            T,
                            inp.T0;
                            atol = FT(atol_temperature),
                            rtol = FT(0),
                        )
                        check_partition(T, inp.ρ0, inp.q0, q_liq, q_ice)
                    end
                end

                # Accuracy check for forced_fixed_iters (GPU optimization)
                # We expect better than 0.1 K accuracy with maxiter=3.
                @testset "forced_fixed_iters accuracy ($FT)" begin
                    fixed_maxiter = 3
                    fixed_tol = FT(1e-4) # Ignored but passed
                    max_err = FT(0)
                    chunk_errs = FT[]
                    for i in idxs
                        inp = targets(i)

                        # Only applicable to ρe formulation
                        res = TD.saturation_adjustment(
                            RS.NewtonsMethod,
                            param_set,
                            TD.ρe(),
                            inp.ρ0,
                            inp.e_int_ρ,
                            inp.q0,
                            fixed_maxiter,
                            fixed_tol,
                            nothing, # T_guess (ignored)
                            true, # forced_fixed_iters
                        )

                        @test res.converged
                        err = abs(res.T - inp.T0)
                        max_err = max(max_err, err)
                        push!(chunk_errs, err)

                        if err > FT(0.1) # Outliers appear to occur, if at all, at unrealistically high temperatures/humidities
                            @info "Outlier detected" FT i err T_ref = inp.T0 p = inp.p0 q_tot =
                                inp.q0 ρ = inp.ρ0
                        end

                        check_partition(res.T, inp.ρ0, inp.q0, res.q_liq, res.q_ice)
                    end
                    @info "Max error forced_fixed_iters (maxiter=$fixed_maxiter)" FT max_err

                    # Allow a few outliers to exceed strict tolerance
                    @test max_err < FT(0.5)

                    # Check percentage of points with error < 0.1 K
                    n_strict = count(e -> e < FT(0.1), chunk_errs)
                    pct_strict = n_strict / length(chunk_errs)
                    @info "Percentage of points with error < 0.1 K" FT pct_strict
                    @test pct_strict > 0.98
                end
            end

            @testset "ρeq solver variants are timed ($FT)" begin
                solvers =
                    (RS.SecantMethod, RS.BrentsMethod, RS.NewtonsMethod, RS.NewtonsMethodAD)
                timing_idxs = idxs[1:min(end, 200)]

                for solver in solvers
                    # Warmup compilation for this solver
                    let inp = targets(timing_idxs[1])
                        TD.saturation_adjustment(
                            solver,
                            param_set,
                            TD.ρe(),
                            inp.ρ0,
                            inp.e_int_ρ,
                            inp.q0,
                            maxiter,
                            tol,
                        )
                    end

                    elapsed = @elapsed begin
                        for i in timing_idxs
                            inp = targets(i)
                            res = TD.saturation_adjustment(
                                solver,
                                param_set,
                                TD.ρe(),
                                inp.ρ0,
                                inp.e_int_ρ,
                                inp.q0,
                                maxiter,
                                tol,
                            )
                            @test isfinite(res.T)
                            @test res.q_liq ≥ zero(FT) && res.q_ice ≥ zero(FT)
                        end
                    end

                    @info "ρeq saturation_adjustment timing (TestedProfiles)" FT solver elapsed_s =
                        elapsed ncases = length(timing_idxs)
                    @test isfinite(elapsed) && elapsed ≥ 0
                end
            end
        end

        @testset "TemperatureProfiles smoke convergence ($FT)" begin
            R_d = FT(TP.R_d(param_set))
            R_v = FT(TP.R_v(param_set))
            T_freeze = FT(TP.T_freeze(param_set))

            # Minimal set of profile types not covered by TestedProfiles directly.
            profiles = (
                TDTP_SA.DryAdiabaticProfile{FT}(param_set),
                TDTP_SA.IsothermalProfile(param_set, FT),
            )

            z_vals = FT.([0, 5e3, 15e3])
            relsat_vals = FT.([0.6, 1.01])

            for profile in profiles, z in z_vals, relsat in relsat_vals
                T0, p0 = profile(param_set, z)
                phase = ifelse(T0 > T_freeze, TD.Liquid(), TD.Ice())
                q_vap_sat_p = TD.q_vap_from_RH(param_set, p0, T0, one(FT), phase)
                q0 = min(FT(0.3), max(zero(FT), relsat * q_vap_sat_p))
                q_vap = min(q0, q_vap_sat_p)
                R_m = R_d * (one(FT) - q0) + R_v * q_vap
                ρ0 = p0 / (R_m * T0)

                # ρeq target
                e_int_ρ = TD.internal_energy_sat(param_set, T0, ρ0, q0)

                # p-based targets computed using the mapping the implementation uses
                ρ_p = TD.air_density(param_set, T0, p0, q0)
                e_int_p = TD.internal_energy_sat(param_set, T0, ρ_p, q0)
                h_p = TD.enthalpy_sat(param_set, T0, ρ_p, q0)

                # Quick smoke: ensure each IndepVars converges and returns a self-consistent partition
                let (; T, q_liq, q_ice) = TD.saturation_adjustment(
                        RS.SecantMethod,
                        param_set,
                        TD.ρe(),
                        ρ0,
                        e_int_ρ,
                        q0,
                        maxiter,
                        tol,
                    )
                    @test isfinite(T)
                    check_partition(T, ρ0, q0, q_liq, q_ice)
                end

                let (; T, q_liq, q_ice) = TD.saturation_adjustment(
                        RS.SecantMethod,
                        param_set,
                        TD.pe(),
                        p0,
                        e_int_p,
                        q0,
                        maxiter,
                        tol,
                    )
                    @test isfinite(T)
                    ρ_eff = TD.air_density(param_set, T, p0, q0)
                    check_partition(T, ρ_eff, q0, q_liq, q_ice)
                end

                let (; T, q_liq, q_ice) = TD.saturation_adjustment(
                        RS.SecantMethod,
                        param_set,
                        TD.ph(),
                        p0,
                        h_p,
                        q0,
                        maxiter,
                        tol,
                    )
                    @test isfinite(T)
                    ρ_eff = TD.air_density(param_set, T, p0, q0)
                    check_partition(T, ρ_eff, q0, q_liq, q_ice)
                end
            end
        end
    end
end
