using Test
using Thermodynamics
import Thermodynamics as TD
const RS = TD.RS
using StaticArrays
using InteractiveUtils
import ClimaParams as CP

# Define parameters
param_set = TD.Parameters.ThermodynamicsParameters(Float32)

# Define inputs
p = 101325.0f0
e_int = 200000.0f0
q_tot = 0.01f0
maxiter = 10
tol = 1.0f-4
T_guess = 290.0f0

# Method
M = RS.SecantMethod

println("Inspecting sa_numerical_method(pe)...")
@code_warntype TD.sa_numerical_method(M, param_set, TD.pe(), p, e_int, q_tot, T_guess)
