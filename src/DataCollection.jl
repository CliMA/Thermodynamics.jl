"""
    DataCollection

This module is designed to help judge the accuracy and
performance for a particular formulation, tolerance, and or
solver configuration, by providing tools to collect various
statistics when Thermodynamic `saturation_adjustment` is called.

## Example:
```
import Thermodynamics as TD
import RootSolvers as RS

function do_work()
    # Calls TD.PhaseEquil_ρeq(VerboseLogger(WarningLogger()), )..., possibly many times
end

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
    return nothing
end

function get_data()
    max_iter = ref_max_iter[]
    call_counter = ref_call_counter[]
    converged_counter = ref_converged_counter[]
    non_converged_counter = ref_non_converged_counter[]
    return (; max_iter, call_counter, converged_counter, non_converged_counter)
end

function print_summary(data)
    max_iter = data.max_iter
    call_counter = data.call_counter
    converged_counter = data.converged_counter
    non_converged_counter = data.non_converged_counter
    average_max_iter = max_iter / call_counter
    @info "Thermodynamics `saturation_adjustment` statistics:" max_iter call_counter average_max_iter converged_counter non_converged_counter
end

end # module
