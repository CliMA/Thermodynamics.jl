# Launch with `julia --project --track-allocation=user`
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
import Profile

include("common.jl")
constructor = ENV["ALLOCATION_CONSTRUCTOR"]
@info "Recording allocations for $constructor"

thermo_state_map = Dict(
    "ρeq" => thermo_state_ρeq,
    "pθq" => thermo_state_pθq,
    "pTq" => thermo_state_pTq,
)
thermo_state = thermo_state_map[constructor]

thermo_state() # compile first
Profile.clear_malloc_data()
thermo_state()

# Quit julia (which generates .mem files), then call
#=
import Coverage
allocs = Coverage.analyze_malloc("src")
=#
