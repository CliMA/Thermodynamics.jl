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
ϵ_numerics(FT) = sqrt(floatmin(FT))
ϵ_numerics(::Type{<:Integer}) = 0
