# Author: Ivan Bioli (https://github.com/IvanBioli)

import Pkg
Pkg.activate("elementifinitiunipv_pkg") # Uncomment this line if you want to activate the package environment
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

    h = 1/N  # Mesh size
    x = LinRange(0, 1, N+1)  # Grid points from 0 to 1 with N intervals
    # Construct the finite difference matrix A
    # Construct the sparse matrix A using spdiagm
    A = (1 / h^2) * spdiagm(
        -1 => -ones(N-2),  # Below diagonal
        0 => 2 * ones(N-1), # Main diagonal
        1 => -ones(N-2)     # Above diagonal
    )
    # Construct the right-hand side vector b (excluding boundaries)
    b = f.(x[2:end-1])  # Values of f(x) at internal points
    b[1] += g_a / h^2   # Apply Dirichlet boundary condition at x=0
    b[end] += g_b / h^2  # Apply Dirichlet boundary condition at x=1
    
    # Solve the linear system
    u_in = A \ b
    # Add boundary values
    u = [g_a; u_in; g_b]
    return u

end

# Function to compute the error between the exact and approximate solution
function compute_error(u_exact, u_approx, x)
    # Compute the error in the infinity norm
    error = maximum(abs.(u_exact.(x) .- u_approx))
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
