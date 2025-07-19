import Random

"""
    sample_range(; param_set, e_int_range, ρ_range, q_tot_range, n_samples)

Generate random thermodynamic states within specified ranges to test constructor robustness.
This helps identify potential convergence issues with the saturation adjustment algorithm.

# Arguments
- `param_set`: Thermodynamic parameter set
- `e_int_range`: Range for internal energy values
- `ρ_range`: Range for density values  
- `q_tot_range`: Range for total specific humidity values
- `n_samples`: Number of random samples to generate

# Returns
- Nothing, but creates thermodynamic states to test constructor behavior
"""
function sample_range(; param_set, e_int_range, ρ_range, q_tot_range, n_samples)
    for i in 1:n_samples
        e_int = Random.rand(e_int_range)
        ρ = Random.rand(ρ_range)
        q_tot = Random.rand(q_tot_range)

        # Test PhaseEquil_ρeq constructor with saturation adjustment
        # This is a critical test as it involves iterative convergence
        # Note: fails with 3 iterations, but succeeds with 4
        ts = PhaseEquil_ρeq(param_set, ρ, e_int, q_tot, 4)
    end
end

@testset "Data tests" begin
    FT = Float64
    param_set = TP.ThermodynamicsParameters(FT)

    # Set random seed for reproducible tests
    Random.seed!(1234)

    # Test constructor robustness with realistic atmospheric ranges
    # These ranges represent typical atmospheric conditions
    sample_range(;
        param_set,
        e_int_range = (28311.801716:30981.514836),  # Internal energy range (J/kg)
        ρ_range = (1.124755:1.129586),              # Density range (kg/m³)
        q_tot_range = (0.011897:0.013305),          # Total humidity range (kg/kg)
        n_samples = 11_000,                         # Large sample size for robustness
    )
end

@testset "pθq data-driven tests" begin
    param_set = TP.ThermodynamicsParameters(Float64)

    # Known problematic cases that cause convergence issues
    # These are cases where the saturation adjustment algorithm fails
    #! format: off
    pθq_broken = [
        # (; p = , θ_liq_ice = , q_tot = ),
        # (; p = , θ_liq_ice = , q_tot = ),
        # (; p = , θ_liq_ice = , q_tot = ),
        # (; p = , θ_liq_ice = , q_tot = ),
        # (; p = , θ_liq_ice = , q_tot = ),
    ]
    #! format: on

    # Configuration for saturation adjustment
    # config = (50, 1e-3, RootSolvers.RegulaFalsiMethod)
    # config = (50, 1e-3, RootSolvers.NewtonMethodAD)
    config = ()

    # Test known problematic cases
    if !isempty(pθq_broken)
        p_arr = getproperty.(pθq_broken, :p)
        θ_liq_ice_arr = getproperty.(pθq_broken, :θ_liq_ice)
        q_tot_arr = getproperty.(pθq_broken, :q_tot)

        for (p, θ_liq_ice, q_tot) in zip(p_arr, θ_liq_ice_arr, q_tot_arr)
            @test_broken begin
                ts = PhaseEquil_pθq(param_set, p, θ_liq_ice, q_tot, config...)
                air_pressure(param_set, ts) == p
            end
        end
    end

    # Test cases that should work correctly
    # These represent realistic atmospheric conditions that should converge
    #! format: off
    pθq = [
        (; p = 82307.30319719888, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 88357.42589002676, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 81090.35731696963, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 75461.41343839701, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 71158.23329080557, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 68032.73180723468, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 59906.16044902247, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 51740.77281055945, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 50056.7894012541, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 48806.70567178993, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
        (; p = 27564.73889538213, θ_liq_ice = 296.7074342326652, q_tot = 0.01894019026929829,),
    ]
    #! format: on

    # Test working cases
    if !isempty(pθq)
        p_arr = getproperty.(pθq, :p)
        θ_liq_ice_arr = getproperty.(pθq, :θ_liq_ice)
        q_tot_arr = getproperty.(pθq, :q_tot)

        for (p, θ_liq_ice, q_tot) in zip(p_arr, θ_liq_ice_arr, q_tot_arr)
            @test begin
                ts = PhaseEquil_pθq(param_set, p, θ_liq_ice, q_tot, config...)
                air_pressure(param_set, ts) == p
            end
        end
    end
end
