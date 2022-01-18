import Pkg
Pkg.develop(path = ".")
import Thermodynamics
import RootSolvers
import CLIMAParameters
import ReportMetrics
const RM = ReportMetrics
td_dir = RM.mod_dir(Thermodynamics)

dirs_to_monitor = [
    ".",
    RM.mod_dir(Thermodynamics),
    RM.mod_dir(RootSolvers),
    RM.mod_dir(CLIMAParameters),
]

for constructor in ["ρeq", "pθq", "pTq"]
    ENV["ALLOCATION_CONSTRUCTOR"] = constructor
    RM.report_allocs(;
        job_name = constructor,
        run_cmd = `$(Base.julia_cmd()) --project=perf/ --track-allocation=all perf/alloc_per_constructor.jl`,
        dirs_to_monitor = dirs_to_monitor,
        process_filename = x -> last(split(x, "packages/")),
    )
end
