"""
    DataCollection

This module is designed to help judge the accuracy and
performance for a particular formulation, tolerance, and or
solver configuration, by providing tools to collect various
statistics when `saturation_adjustment` is called.

## Example:
```
import Thermodynamics as TD
import RootSolvers as RS

function do_work()
    # Calls TD.saturation_adjustment()..., possibly many times
end

TD.solution_type() = RS.VerboseSolution()
do_work()
TD.DataCollection.print_summary()
```
!!! warn
    This data collection was designed for unthreaded single processor
    runs, and may not work correctly for threaded / multi-processor runs.
"""
module DataCollection

import RootSolvers
const RS = RootSolvers

# Stats to collect
const ref_max_iter = Ref{Int}(0)
const ref_call_counter = Ref{Int}(0)
const ref_converged_counter = Ref{Int}(0)
const ref_non_converged_counter = Ref{Int}(0)
const ref_iter_performed = Ref{Int}(0)

@inline log_meta(sol::RS.CompactSolutionResults) = nothing

function log_meta(sol::RS.VerboseSolutionResults)
    ref_max_iter[] = max(ref_max_iter[], sol.iter_performed)
    if sol.converged
        ref_converged_counter[] += 1
    else
        ref_non_converged_counter[] += 1
    end
    ref_call_counter[] += 1
    ref_iter_performed[] += sol.iter_performed
    return nothing
end

function reset_stats()
    ref_max_iter[] = 0
    ref_call_counter[] = 0
    ref_converged_counter[] = 0
    ref_non_converged_counter[] = 0
    ref_iter_performed[] = 0
    return nothing
end

function get_data()
    max_iter = ref_max_iter[]
    call_counter = ref_call_counter[]
    converged_counter = ref_converged_counter[]
    non_converged_counter = ref_non_converged_counter[]
    iter_performed = ref_iter_performed[]
    return (;
        max_iter,
        call_counter,
        iter_performed,
        converged_counter,
        non_converged_counter,
    )
end

function print_summary()
    return print_summary(get_data())
end

function print_summary(data)
    max_iter = data.max_iter
    call_counter = data.call_counter
    iter_performed = data.iter_performed
    converged_counter = data.converged_counter
    non_converged_counter = data.non_converged_counter
    mean_iter_performed = ifelse(call_counter > 0, iter_performed / call_counter, NaN)
    @info "Thermodynamics `saturation_adjustment` statistics:" max_iter call_counter mean_iter_performed converged_counter non_converged_counter
end

end # module
