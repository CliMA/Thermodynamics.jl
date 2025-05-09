@time "measure precompile" begin
    import Thermodynamics as TD
    params = TD.Parameters.ThermodynamicsParameters{Float64}(;
        T_0 = 273.16,
        MSLP = 101325,
        p_ref_theta = 101325,
        cp_v = 1859,
        cp_l = 4181,
        cp_i = 2100,
        LH_v0 = 2500800,
        LH_s0 = 2834400,
        press_triple = 611.657,
        T_triple = 273.16,
        T_freeze = 273.15,
        T_min = 220,
        T_max = 1000,
        entropy_reference_temperature = 298.15,
        entropy_dry_air = 6864.8,
        entropy_water_vapor = 10513.6,
        kappa_d = 0.28571428571,
        gas_constant = 8.3144598,
        molmass_dryair = 0.02897,
        molmass_water = 0.01801528,
        T_surf_ref = 290,
        T_min_ref = 220,
        grav = 9.8,
        T_icenuc = 233,
        pow_icenuc = 1,
    );
    ts = TD.PhaseEquil_ρθq(params, 1.0, 374.0, 0.01);
    T = TD.air_temperature(params, ts)
end
nothing
