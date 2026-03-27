# Author: Ivan Bioli (https://github.com/IvanBioli)

import Pkg
Pkg.activate("elementifinitiunipv_pkg") # Uncomment this line if you want to activate the package environment
using Revise
includet("../../modules/Meshing.jl")
using Plots

function add_triangle!(p1, p2, p3)
    x = [p1[1], p2[1], p3[1], p1[1]]  # x-coordinates (closed loop)
    y = [p1[2], p2[2], p3[2], p1[2]]  # y-coordinates (closed loop)
    plot!(x, y, label="", linewidth=1, color=:black)
end

# Function to plot the mesh
function plot_mesh(T, p)
    fig = plot() # Initialize an empty plot
    for j in eachindex(axes(T, 2))
        v1 = p[:,T[1,j]]
        v2 = p[:,T[2,j]]
        v3 = p[:,T[3,j]]
        add_triangle!(v1, v2, v3)
    end
    plot!(fig, aspect_ratio=1)
    return fig 
end

# Create the mesh
h = 0.05
out_file = mesh_square(h)
T, p = get_nodes_connectivity(out_file)

# Plot the mesh
fig = plot_mesh(T, p)
savefig("figures/ex02_1.pdf")
display(fig)

# Compare with Meshes.jl
import Meshes
mesh = to_Meshes(T, p)
Meshes.viz(mesh, showsegments = true)
