"""
    AbstractThermodynamicsLogger

A Thermodynamics logger can be passed to
thermo state constructors that (may) perform
saturation adjustment. This logger will statically
decide if a warning statement is printed, or
if the function should error on non-convergence.
"""
abstract type AbstractThermodynamicsLogger end
const ATL = AbstractThermodynamicsLogger

"""
    WarningLogger

Warn when saturation adjustment did not converge.
"""
struct WarningLogger <: AbstractThermodynamicsLogger end

"""
    ErrorLogger

Error (but do not warn) when saturation adjustment
did not converge.

This is really only useful in our development docs.
"""
struct ErrorLogger <: AbstractThermodynamicsLogger end

"""
    WarnAndErrorLogger

Warn and error (default) when saturation adjustment did not converge.
"""
struct WarnAndErrorLogger <: AbstractThermodynamicsLogger end

"""
    VerboseLogger

A verbose logger, for when we use RootSolvers `VerboseSolution` type.
"""
struct VerboseLogger{L} <: AbstractThermodynamicsLogger
    logger::L
end

"""
    NullLogger

Do nothing when saturation adjustment did not converge.
"""
struct NullLogger <: AbstractThermodynamicsLogger end

warn_msg(::WarningLogger) = true
warn_msg(::WarnAndErrorLogger) = true
warn_msg(l::VerboseLogger) = warn_msg(l.logger)
warn_msg(::AbstractThermodynamicsLogger) = false

error_msg(::ErrorLogger) = true
error_msg(::WarnAndErrorLogger) = true
error_msg(l::VerboseLogger) = error_msg(l.logger)
error_msg(::AbstractThermodynamicsLogger) = false

Base.broadcastable(x::AbstractThermodynamicsLogger) = tuple(x)
