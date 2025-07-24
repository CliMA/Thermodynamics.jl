# Auxiliary functions for Thermodynamics.jl

"""
    ReLU(x)

Rectified Linear Unit function. Returns the maximum of zero and the input value.

# Arguments
- `x`: Input value

# Returns
- Maximum of zero and `x`

# Examples
```julia
ReLU(3.0)  # returns 3.0
ReLU(-2.0)  # returns 0.0
ReLU(0.0)   # returns 0.0
```
"""
@inline ReLU(x) = max(zero(x), x)


"""
    fast_power(x, y)

Fast power function using exponential and logarithm for bases very close to 1.
This function provides better performance than the standard `^` operator for cases
where the base is very close to 1, avoiding slow paths in the standard library.

# Arguments
- `x`: Base value (should be close to 1 for optimal performance)
- `y`: Exponent value

# Returns
- `x^y` calculated using `exp(y * log(x))`

# Notes
This function is specifically optimized for cases where the base `x` is very close to 1,
as the standard `^` operator can have performance issues in these scenarios.
See: https://stackoverflow.com/questions/14687665/very-slow-stdpow-for-bases-very-close-to-1

# Examples
```julia
fast_power(1.001, 2.5)  # More efficient than 1.001^2.5
fast_power(0.999, 3.0)  # Handles bases close to 1 efficiently
```
"""
@inline fast_power(x, y) = exp(y * log(x))