# Author: Ivan Bioli (https://github.com/IvanBioli)

import Pkg
# Pkg.activate("elementifinitiunipv_pkg") # Uncomment this line if you want to activate the package environment
using Revise
using LinearAlgebra
using SparseArrays
using Plots
using LaTeXStrings

# Define the exact solution and problem data
u(x) = cos(π * x)            # Exact solution u(x) = cos(πx)
f(x) = π^2 * cos(π * x)      # Right-hand side f(x) = π²cos(πx)

# Boundary conditions (Dirichlet)
g0, g1 = u(0), u(1)          # Dirichlet boundary conditions u(0) = g0 and u(1) = g1

# Function to solve the Poisson equation using finite differences
function poisson1d(N, f, g_a, g_b)

    ############### ADD CODE HERE ###############

end

# Function to compute the error between the exact and approximate solution
function compute_error(u_exact, u_approx, x)
    # Compute the error in the infinity norm
    error = ############### ADD CODE HERE ###############
    return error
end

# Mesh sizes to test different grid resolutions
Nvals = [10, 20, 40, 80, 160, 320]  # Testing for N = 10, 100, 1000, ...
mesh_sizes = Float64[]                        # Store mesh sizes
errors = Float64[]                            # Store errors

# Loop over different values of N
for N in Nvals
    h = 1 / N  # Mesh size
    push!(mesh_sizes, h)
    
    # Compute the approximate solution using poisson1d
    uh = poisson1d(N, f, g0, g1)
    
    # Compute the error
    x = LinRange(0, 1, N+1)  # Grid points from 0 to 1 with N intervals
    err = compute_error(u, uh, x)
    push!(errors, err)
end

# Plot error vs mesh size h (log-log scale)
plt = plot(mesh_sizes, errors, xscale=:log10, yscale=:log10, marker=:o, label="FD Method",
    xlabel=L"h", ylabel=L"$\max_i |u(x_i) - u_i|$", title="Finite Differences Error")
plot!(mesh_sizes, mesh_sizes.^2, xscale=:log10, yscale=:log10, linestyle=:dash, label=L"\mathcal{O}(h^2)")  # Reference line
# Show and save the figure
savefig("figures/ex01_1.pdf")
display(plt)
