# Auxiliary functions for Thermodynamics.jl

"""
    ReLU(x)

Internal function. Rectified Linear Unit: returns `max(0, x)`.
"""
@inline ReLU(x) = max(zero(x), x)

"""
    fast_power(x, y)

Internal function. Fast power function using `exp(y * log(x))`.

This is faster than Julia's `^` operator for bases very close to 1. Julia's `^`
dispatches through a general code path that handles complex numbers and edge cases;
for real positive `x` near 1 (as occurs in the Clausius-Clapeyron computation
where `x = T / T_triple ≈ 1`), `exp(y * log(x))` avoids those branches and is
significantly faster.

Note: requires `x > 0` since it uses `log(x)`.
"""
@inline fast_power(x, y) = exp(y * log(x))

"""
    ϵ_numerics(FT)

Smallest acceptable number that is different than zero.
"""
@inline ϵ_numerics(FT) = sqrt(floatmin(FT))

"""
    T_positive_floor(FT)

Internal function. Smallest temperature at which the saturation functions and their
derivatives are numerically evaluable.

Every quantity the saturation-adjustment iteration evaluates is finite for any strictly
positive temperature: the saturation vapor pressure underflows to zero, and the `1/T`,
`1/T²` and `log(T)` factors stay finite. At exactly zero they produce NaN, and below zero
the logarithms throw. `sqrt(eps(FT))` (≈ 1.5e-8 K in `Float64`, ≈ 3.5e-4 K in `Float32`)
is far below any physical temperature yet keeps `L/(R_v T²)` finite in `Float32`, unlike
`ϵ_numerics`, which is small enough to overflow that quotient.

This bounds the *search* of the solvers, never their answer: unsaturated states return
their exact temperature, and no saturated state has a solution this cold.
"""
@inline T_positive_floor(::Type{FT}) where {FT} = sqrt(eps(FT))

# Integer arguments reach this when a caller passes an integer literal `0` for a condensate
# specific humidity and the element type is taken from that argument. `floatmin` has no
# integer method, so without this the call is a MethodError rather than a no-op regularization.
@inline ϵ_numerics(::Type{<:Integer}) = 0
