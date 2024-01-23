# Launch with `julia --project --track-allocation=user`
import Profile

include("common.jl")
constructor = ENV["ALLOCATION_CONSTRUCTOR"]

thermo_state_map = Dict(
    "ρeq" => thermo_state_ρeq,
    "pθq" => thermo_state_pθq,
    "pTq" => thermo_state_pTq,
)
thermo_state = thermo_state_map[constructor]

thermo_state() # compile first
Profile.clear_malloc_data()
thermo_state()
