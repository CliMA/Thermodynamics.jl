import Thermodynamics
import RootSolvers
import CLIMAParameters
import ReportMetrics

dirs_to_monitor =
    [".", pkgdir(Thermodynamics), pkgdir(RootSolvers), pkgdir(CLIMAParameters)]

for constructor in ["ρeq", "pθq", "pTq"]
    ENV["ALLOCATION_CONSTRUCTOR"] = constructor
    ReportMetrics.report_allocs(;
        job_name = constructor,
        run_cmd = `$(Base.julia_cmd()) --project=perf/ --track-allocation=all perf/alloc_per_constructor.jl`,
        dirs_to_monitor = dirs_to_monitor,
        process_filename = x -> last(split(x, "packages/")),
    )
end
