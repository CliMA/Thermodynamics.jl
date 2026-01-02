# Auxiliary functions for Thermodynamics.jl

"""
    ReLU(x)

Internal function. Rectified Linear Unit: returns `max(0, x)`.
"""
@inline ReLU(x) = max(zero(x), x)

"""
    fast_power(x, y)

Internal function. Fast power function using `exp(y * log(x))`.

This provides better performance than the standard `^` operator for bases very
close to 1, avoiding slow paths in the standard library.

Note: requires `x > 0` since it uses `log(x)`.

See: https://stackoverflow.com/questions/14687665/very-slow-stdpow-for-bases-very-close-to-1
"""
@inline fast_power(x, y) = exp(y * log(x))
