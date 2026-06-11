# Author: Ivan Bioli (https://github.com/IvanBioli)
begin
    import Pkg
    Pkg.activate("elementifinitiunipv_pkg")
    using Revise

    # Load the necessary files
    includet("../../modules/Meshing.jl")
    includet("../../modules/Quadrature_adv.jl")
    includet("../../modules/Assembly.jl")
    includet("../../modules/Assembly_mixed.jl")

    import Meshes
    using Plots
    using LaTeXStrings
    import GridapGmsh
end

##################################################################
begin
begin
    # Define data
    p(x) = 0.25 * (1 - x[1]^2 - x[2]^2)
    ∇p(x) = [-0.5 * x[1], -0.5 * x[2]]
    Δp(x) = -1
    k(x) = 1
    ∇k(x) = [0.0, 0]

    μ(x) = 1 ./ k(x)
    u(x) = - k(x) * ∇p(x)
    f(x) = -Δp(x) * k(x) + ∇p(x) ⋅ ∇k(x)
    divu(x) = f(x)
end


begin
msh_sizes = 2 .^ range(-3, -5, 6)
L2_error_p = similar(msh_sizes)
H1div_error_u = similar(msh_sizes)
for (i, h) in enumerate(msh_sizes)
    # Build the mesh with mesh-size h    
    out_file = mesh_circle(h)
    edges2nodes, elems2edges, elems2orientation = get_edges_info(out_file)
    msh = Mesh(get_nodes_connectivity(out_file)...)
    set_edges_info!(msh, edges2nodes, elems2edges, elems2orientation)

    # Assemble the system
    initialize_assembly_mixed!(msh)
    local_assembler(Ae, Be, fe, msh, cell_index) = darcy_assemble_local_mixed!(Ae, Be, fe, msh, cell_index, f, μ)
    K, rhs = assemble_global_mixed(msh, local_assembler)

    # Solve the linear system
    x = K \ rhs
    uh, ph = x[1:get_nedges(msh)], x[get_nedges(msh)+1:end]

    # Compute error
    L2_error_p[i] = L2error_mixed_p(p, ph, msh, Q2_ref)
    H1div_error_u[i] = H1diverror_mixed_u(u, divu, uh, msh, Q2_ref)
end
end

import Plots
begin
    plt = plot(msh_sizes, L2_error_p, xaxis=:log, yaxis=:log, label=L"${\|\|p-p_h \|\|}_{L^2}$", marker=:circle)
    plot!(msh_sizes, L2_error_p[1] * msh_sizes ./ msh_sizes[1], xaxis=:log, yaxis=:log, linestyle=:dash, label=L"$\mathcal{O}(h)$")
    xlabel!(plt, L"$h$")
    ylabel!(plt, "Errore")
    plot!(legend=:bottomright)
    savefig(plt, "./figures_julia/ex08_1_L2error.pdf")
    plot(plt) # Show plot
end

begin
    plt = plot(msh_sizes, H1div_error_u, xaxis=:log, yaxis=:log, label=L"${\|\|u-u_h \|\|}_{H^1(\mathrm{div})}$", marker=:circle)
    plot!(msh_sizes, H1div_error_u[1] * msh_sizes ./ msh_sizes[1], xaxis=:log, yaxis=:log, linestyle=:dash, label=L"$\mathcal{O}(h)$")
    xlabel!(plt, L"$h$")
    ylabel!(plt, "Errore")
    plot!(legend=:bottomright)
    savefig(plt, "./figures_julia/ex08_1_H1error.pdf")
    plot(plt) # Show plot
end
end


##################################################################
begin
begin
    # Define data
    p(x) = sin(2 * pi*x[1]) * sin(pi*x[2])
    ∇p(x) = [2 * pi * cos(2 * pi*x[1]) * sin(pi*x[2]), pi * sin(2 * pi*x[1]) * cos(pi*x[2])]
    Δp(x) = - 5 * pi^2 * p(x) 
    k(x) = 2
    ∇k(x) = [0.0, 0]

    μ(x) = 1 ./ k(x)
    u(x) = - k(x) * ∇p(x)
    f(x) = -Δp(x) * k(x) + ∇p(x) ⋅ ∇k(x)
    divu(x) = f(x)
end


begin
msh_sizes = 2 .^ range(-3, -5, 6)
L2_error_p = similar(msh_sizes)
H1div_error_u = similar(msh_sizes)
for (i, h) in enumerate(msh_sizes)
    # Build the mesh with mesh-size h    
    out_file = mesh_square(h)
    edges2nodes, elems2edges, elems2orientation = get_edges_info(out_file)
    msh = Mesh(get_nodes_connectivity(out_file)...)
    set_edges_info!(msh, edges2nodes, elems2edges, elems2orientation)

    # Assemble the system
    initialize_assembly_mixed!(msh)
    local_assembler(Ae, Be, fe, msh, cell_index) = darcy_assemble_local_mixed!(Ae, Be, fe, msh, cell_index, f, μ)
    K, rhs = assemble_global_mixed(msh, local_assembler)

    # Solve the linear system
    x = K \ rhs
    uh, ph = x[1:get_nedges(msh)], x[get_nedges(msh)+1:end]

    # Compute error
    L2_error_p[i] = L2error_mixed_p(p, ph, msh, Q2_ref)
    H1div_error_u[i] = H1diverror_mixed_u(u, divu, uh, msh, Q2_ref)
end
end

import Plots
begin
    plt = plot(msh_sizes, L2_error_p, xaxis=:log, yaxis=:log, label=L"${\|\|p-p_h \|\|}_{L^2}$", marker=:circle)
    plot!(msh_sizes, L2_error_p[1] * msh_sizes ./ msh_sizes[1], xaxis=:log, yaxis=:log, linestyle=:dash, label=L"$\mathcal{O}(h)$")
    xlabel!(plt, L"$h$")
    ylabel!(plt, "Errore")
    plot!(legend=:bottomright)
    savefig(plt, "./figures_julia/ex08_2_L2error.pdf")
    plot(plt) # Show plot
end

begin
    plt = plot(msh_sizes, H1div_error_u, xaxis=:log, yaxis=:log, label=L"${\|\|u-u_h \|\|}_{H^1(\mathrm{div})}$", marker=:circle)
    plot!(msh_sizes, H1div_error_u[1] * msh_sizes ./ msh_sizes[1], xaxis=:log, yaxis=:log, linestyle=:dash, label=L"$\mathcal{O}(h)$")
    xlabel!(plt, L"$h$")
    ylabel!(plt, "Errore")
    plot!(legend=:bottomright)
    savefig(plt, "./figures_julia/ex08_2_H1error.pdf")
    plot(plt) # Show plot
end
end

