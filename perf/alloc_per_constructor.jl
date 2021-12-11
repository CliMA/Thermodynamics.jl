# Launch with `julia --project --track-allocation=user`
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
import Profile

include("common.jl")
constructor = ENV["ALLOCATION_CONSTRUCTOR"]
@info "Recording allocations for $constructor"


if constructor == "ρeq"
    thermo_state = thermo_state_ρeq_newton
else
end
if constructor == "pθq"
    thermo_state = thermo_state_pθq
else
end

thermo_state() # compile first
Profile.clear_malloc_data()
thermo_state()

# Quit julia (which generates .mem files), then call
#=
import Coverage
allocs = Coverage.analyze_malloc("src")
=#
