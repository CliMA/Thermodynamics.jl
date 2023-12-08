import CLIMAParameters as CP

function get_parameter_set(::Type{FT}) where {FT}
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    aliases = string.(fieldnames(Parameters.ThermodynamicsParameters))
    param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    return Parameters.ThermodynamicsParameters{FT}(; param_pairs...)
end

parameter_set(::Type{Float64}) = get_parameter_set(Float64)
parameter_set(::Type{Float32}) = get_parameter_set(Float32)

@setup_workload begin
    @compile_workload begin
        for ArrayType in [Array{Float32}, Array{Float64}]
            FT = eltype(ArrayType)
            param_set = parameter_set(FT)
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

            ρu = FT[1, 2, 3]
            typeof.(internal_energy.(ρ, ρ .* e_int, Ref(ρu), e_pot)) ==
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

            all(
                air_temperature.(param_set, ts_T) .≈ air_temperature.(param_set, ts_Tp),
            )
            all(air_pressure.(param_set, ts_T) .≈ air_pressure.(param_set, ts_Tp))
            all(
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
                typeof.(soundspeed_air.(param_set, ts)) == typeof.(e_int)
                typeof.(gas_constant_air.(param_set, ts)) == typeof.(e_int)
                typeof.(specific_enthalpy.(param_set, ts)) == typeof.(e_int)
                typeof.(vapor_specific_humidity.(param_set, ts)) == typeof.(e_int)
                typeof.(relative_humidity.(param_set, ts)) == typeof.(e_int)
                typeof.(air_pressure.(param_set, ts)) == typeof.(e_int)
                typeof.(air_density.(param_set, ts)) == typeof.(e_int)
                typeof.(total_specific_humidity.(param_set, ts)) == typeof.(e_int)
                typeof.(liquid_specific_humidity.(param_set, ts)) ==
                      typeof.(e_int)
                typeof.(ice_specific_humidity.(param_set, ts)) == typeof.(e_int)
                typeof.(cp_m.(param_set, ts)) == typeof.(e_int)
                typeof.(cv_m.(param_set, ts)) == typeof.(e_int)
                typeof.(air_temperature.(param_set, ts)) == typeof.(e_int)
                typeof.(internal_energy_sat.(param_set, ts)) == typeof.(e_int)
                typeof.(internal_energy.(param_set, ts)) == typeof.(e_int)
                typeof.(internal_energy_dry.(param_set, ts)) == typeof.(e_int)
                typeof.(internal_energy_vapor.(param_set, ts)) == typeof.(e_int)
                typeof.(internal_energy_liquid.(param_set, ts)) == typeof.(e_int)
                typeof.(internal_energy_ice.(param_set, ts)) == typeof.(e_int)
                typeof.(latent_heat_vapor.(param_set, ts)) == typeof.(e_int)
                typeof.(latent_heat_sublim.(param_set, ts)) == typeof.(e_int)
                typeof.(latent_heat_fusion.(param_set, ts)) == typeof.(e_int)
                typeof.(q_vap_saturation.(param_set, ts)) == typeof.(e_int)
                typeof.(q_vap_saturation_liquid.(param_set, ts)) == typeof.(e_int)
                typeof.(q_vap_saturation_ice.(param_set, ts)) == typeof.(e_int)
                typeof.(saturation_excess.(param_set, ts)) == typeof.(e_int)
                typeof.(liquid_fraction.(param_set, ts)) == typeof.(e_int)
                typeof.(liquid_ice_pottemp.(param_set, ts)) == typeof.(e_int)
                typeof.(dry_pottemp.(param_set, ts)) == typeof.(e_int)
                typeof.(exner.(param_set, ts)) == typeof.(e_int)
                typeof.(liquid_ice_pottemp_sat.(param_set, ts)) == typeof.(e_int)
                typeof.(specific_volume.(param_set, ts)) == typeof.(e_int)
                typeof.(supersaturation.(param_set, ts, Ice())) == typeof.(e_int)
                typeof.(supersaturation.(param_set, ts, Liquid())) ==
                      typeof.(e_int)
                typeof.(virtual_pottemp.(param_set, ts)) == typeof.(e_int)
                typeof.(specific_entropy.(param_set, ts)) == typeof.(e_int)
                eltype.(gas_constants.(param_set, ts)) == typeof.(e_int)

                typeof.(total_specific_enthalpy.(param_set, ts, e_tot)) ==
                      typeof.(e_int)
                typeof.(moist_static_energy.(param_set, ts, e_pot)) ==
                      typeof.(e_int)
                typeof.(getproperty.(PhasePartition.(param_set, ts), :tot)) ==
                      typeof.(e_int)
                typeof.(virtual_dry_static_energy.(param_set, ts, e_pot)) ==
                      typeof.(e_int)
            end


            profiles = TestedProfiles.PhaseDryProfiles(param_set, ArrayType)
            for prof in profiles
                (; T, p, e_int, ρ, θ_liq_ice, phase_type) = prof
                (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = prof

                θ_dry = dry_pottemp(param_set, T, ρ)
                ts_dry = PhaseDry(param_set, e_int, ρ)
                ts_dry_ρp = PhaseDry_ρp(param_set, ρ, p)
                ts_dry_pT = PhaseDry_pT(param_set, p, T)
                ts_dry_ρθ = PhaseDry_ρθ(param_set, ρ, θ_dry)
                ts_dry_pθ = PhaseDry_pθ(param_set, p, θ_dry)
            end

            profiles = TestedProfiles.PhaseEquilProfiles(param_set, ArrayType)
            for prof in profiles
                (; T, p, e_int, ρ, θ_liq_ice, phase_type) = prof
                (; q_tot, q_liq, q_ice, q_pt, RH, e_kin, e_pot) = prof

                ρu = FT[1, 2, 3]
                typeof(internal_energy(ρ, ρ * e_int, ρu, e_pot)) ==
                      typeof(e_int)

                ts_eq = PhaseEquil_ρeq(param_set, ρ, e_int, q_tot, 15, FT(1e-4))
                e_tot = total_energy(param_set, ts_eq, e_kin, e_pot)

                ts_T =
                    PhaseEquil_ρTq(
                        param_set,
                        air_density(param_set, ts_eq),
                        air_temperature(param_set, ts_eq),
                        q_tot,
                    )
                ts_Tp =
                    PhaseEquil_pTq(
                        param_set,
                        air_pressure(param_set, ts_eq),
                        air_temperature(param_set, ts_eq),
                        q_tot,
                    )

                ts_ρp =
                    PhaseEquil_ρpq(
                        param_set,
                        air_density(param_set, ts_eq),
                        air_pressure(param_set, ts_eq),
                        q_tot,
                    )

                all(
                    air_temperature(param_set, ts_T) ≈ air_temperature(param_set, ts_Tp),
                )
                all(air_pressure(param_set, ts_T) ≈ air_pressure(param_set, ts_Tp))
                all(
                    total_specific_humidity(param_set, ts_T) ≈
                    total_specific_humidity(param_set, ts_Tp),
                )

                ts_neq = PhaseNonEquil(param_set, e_int, ρ, q_pt)
                ts_ρT_neq = PhaseNonEquil_ρTq(param_set, ρ, T, q_pt)
                ts_pT_neq = PhaseNonEquil_pTq(param_set, p, T, q_pt)

                ts_θ_liq_ice_eq =
                    PhaseEquil_ρθq(param_set, ρ, θ_liq_ice, q_tot, 45, FT(1e-4))
                ts_θ_liq_ice_eq_p =
                    PhaseEquil_pθq(param_set, p, θ_liq_ice, q_tot, 40, FT(1e-4))
                ts_θ_liq_ice_neq = PhaseNonEquil_ρθq(param_set, ρ, θ_liq_ice, q_pt)
                ts_θ_liq_ice_neq_p = PhaseNonEquil_pθq(param_set, p, θ_liq_ice, q_pt)

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
                    typeof(soundspeed_air(param_set, ts)) == typeof(e_int)
                    typeof(gas_constant_air(param_set, ts)) == typeof(e_int)
                    typeof(specific_enthalpy(param_set, ts)) == typeof(e_int)
                    typeof(vapor_specific_humidity(param_set, ts)) == typeof(e_int)
                    typeof(relative_humidity(param_set, ts)) == typeof(e_int)
                    typeof(air_pressure(param_set, ts)) == typeof(e_int)
                    typeof(air_density(param_set, ts)) == typeof(e_int)
                    typeof(total_specific_humidity(param_set, ts)) == typeof(e_int)
                    typeof(liquid_specific_humidity(param_set, ts)) ==
                          typeof(e_int)
                    typeof(ice_specific_humidity(param_set, ts)) == typeof(e_int)
                    typeof(cp_m(param_set, ts)) == typeof(e_int)
                    typeof(cv_m(param_set, ts)) == typeof(e_int)
                    typeof(air_temperature(param_set, ts)) == typeof(e_int)
                    typeof(internal_energy_sat(param_set, ts)) == typeof(e_int)
                    typeof(internal_energy(param_set, ts)) == typeof(e_int)
                    typeof(internal_energy_dry(param_set, ts)) == typeof(e_int)
                    typeof(internal_energy_vapor(param_set, ts)) == typeof(e_int)
                    typeof(internal_energy_liquid(param_set, ts)) == typeof(e_int)
                    typeof(internal_energy_ice(param_set, ts)) == typeof(e_int)
                    typeof(latent_heat_vapor(param_set, ts)) == typeof(e_int)
                    typeof(latent_heat_sublim(param_set, ts)) == typeof(e_int)
                    typeof(latent_heat_fusion(param_set, ts)) == typeof(e_int)
                    typeof(q_vap_saturation(param_set, ts)) == typeof(e_int)
                    typeof(q_vap_saturation_liquid(param_set, ts)) == typeof(e_int)
                    typeof(q_vap_saturation_ice(param_set, ts)) == typeof(e_int)
                    typeof(saturation_excess(param_set, ts)) == typeof(e_int)
                    typeof(liquid_fraction(param_set, ts)) == typeof(e_int)
                    typeof(liquid_ice_pottemp(param_set, ts)) == typeof(e_int)
                    typeof(dry_pottemp(param_set, ts)) == typeof(e_int)
                    typeof(exner(param_set, ts)) == typeof(e_int)
                    typeof(liquid_ice_pottemp_sat(param_set, ts)) == typeof(e_int)
                    typeof(specific_volume(param_set, ts)) == typeof(e_int)
                    typeof(supersaturation(param_set, ts, Ice())) == typeof(e_int)
                    typeof(supersaturation(param_set, ts, Liquid())) ==
                          typeof(e_int)
                    typeof(virtual_pottemp(param_set, ts)) == typeof(e_int)
                    typeof(specific_entropy(param_set, ts)) == typeof(e_int)
                    eltype(gas_constants(param_set, ts)) == typeof(e_int)

                    typeof(total_specific_enthalpy(param_set, ts, e_tot)) ==
                          typeof(e_int)
                    typeof(moist_static_energy(param_set, ts, e_pot)) ==
                          typeof(e_int)
                    typeof(getproperty(PhasePartition(param_set, ts), :tot)) ==
                          typeof(e_int)
                    typeof(virtual_dry_static_energy(param_set, ts, e_pot)) ==
                          typeof(e_int)
                end
            end


        end

    end
end
