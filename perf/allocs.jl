import Thermodynamics
import RootSolvers
import ClimaParams
import Profile, ProfileCanvas

include("common.jl")
output_dir = joinpath(pkgdir(Thermodynamics), "perf", "allocations_output")
mkpath(output_dir)
thermo_state_ρeq(param_set, ρ, e_int, q_tot)
@show @allocated thermo_state_ρeq(param_set, ρ, e_int, q_tot)
Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate = 1 thermo_state_ρeq(
    param_set,
    ρ,
    e_int,
    q_tot,
)
results = Profile.Allocs.fetch()
Profile.Allocs.clear()
profile = ProfileCanvas.view_allocs(results)
ProfileCanvas.html_file(joinpath(output_dir, "allocs_ρeq.html"), profile)

thermo_state_pθq(param_set, p, θ_liq_ice, q_tot)
@show @allocated thermo_state_pθq(param_set, p, θ_liq_ice, q_tot)
Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate = 1 thermo_state_pθq(
    param_set,
    p,
    θ_liq_ice,
    q_tot,
)
results = Profile.Allocs.fetch()
Profile.Allocs.clear()
profile = ProfileCanvas.view_allocs(results)
ProfileCanvas.html_file(joinpath(output_dir, "allocs_pθq.html"), profile)

thermo_state_pTq(param_set, p, T, q_tot)
@show @allocated thermo_state_pTq(param_set, p, T, q_tot)
Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate = 1 thermo_state_pTq(param_set, p, T, q_tot)
results = Profile.Allocs.fetch()
Profile.Allocs.clear()
profile = ProfileCanvas.view_allocs(results)
ProfileCanvas.html_file(joinpath(output_dir, "allocs_pTq.html"), profile)
