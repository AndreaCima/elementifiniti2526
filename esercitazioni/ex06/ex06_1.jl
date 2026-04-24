# Author: Ivan Bioli (https://github.com/IvanBioli)
begin
    import Pkg
    Pkg.activate("elementifinitiunipv_pkg")
    using Revise

    # Load the necessary files
    includet("../../modules/Meshing.jl")
    includet("../../modules/Quadrature_adv.jl")
    includet("../../modules/Assembly.jl")

    import Meshes
    using Plots
    using LaTeXStrings
    import GridapGmsh
end


################################# k(x, y) = 1 #################################
begin
    # Define data
    f(x) = exp(2*x[1] + x[2])
    k(x) = 1
    g(x) = 0

    # Choose mesh-size
    h = 0.1

    # Build the mesh with mesh-size h    
    out_file = mesh_square(h; v1=[-1,-1], v2=[1,1])
    T, p = get_nodes_connectivity(out_file)
    msh = Mesh(T, p)
    # Get Dirichlet dofs
    bnd_tags, bnd_coords = get_boundary_nodes(out_file; labels = ["left"])
    set_dirichletdofs!(msh, bnd_tags)

    # Assemble the system
    initialize_assembly!(msh)
    local_assembler(Ke, fe, msh, cell_index) = darcy_assemble_local!(Ke, fe, msh, cell_index, f, k)
    A, b = assemble_global(msh, local_assembler)
    # Condense system imposing Dirichlet BCs
    A_cond, b_cond, uh = impose_dirichlet(A, b, g, msh)
    # Solve linear system
    uh[get_freedofs(msh)] = A_cond \ b_cond

    # Plot the solution
    plt = plot_flat(msh, uh; plot_msh = true)
    savefig(plt, "./figures_julia/ex06_1_solution_flat.pdf") # Save plot

    plt = plot_surf(msh, uh; plot_msh = true)
    PlotlyJS.savefig(plt, "./figures_julia/ex06_1_solution_surf.pdf")
    display(plt) # Save plot
end

############################## k(x, y) = 1+x^2+y^2 ##############################
begin
    # Define data
    f(x) = exp(2*x[1] + x[2])
    k(x) = 1 + x[1]^2 + x[2]^2
    g(x) = 0
    
    # Choose mesh-size
    h = 0.1
    
    # Build the mesh with mesh-size h    
    out_file = mesh_square(h; v1=[-1,-1], v2=[1,1])
    T, p = get_nodes_connectivity(out_file)
    msh = Mesh(T, p)
    # Get Dirichlet dofs
    bnd_tags, bnd_coords = get_boundary_nodes(out_file; labels = ["left"])
    set_dirichletdofs!(msh, bnd_tags)
    
    # Assemble the system
    initialize_assembly!(msh)
    local_assembler(Ke, fe, msh, cell_index) = darcy_assemble_local!(Ke, fe, msh, cell_index, f, k)
    A, b = assemble_global(msh, local_assembler)
    # Condense system imposing Dirichlet BCs
    A_cond, b_cond, uh = impose_dirichlet(A, b, g, msh)
    # Solve linear system
    uh[get_freedofs(msh)] = A_cond \ b_cond
    
    # Plot the solution
    plt = plot_flat(msh, uh; plot_msh = true)
    savefig(plt, "./figures_julia/ex06_2_solution_flat.pdf") # Save plot
    
    plt = plot_surf(msh, uh; plot_msh = true)
    PlotlyJS.savefig(plt, "./figures_julia/ex06_2_solution_surf.pdf")
    display(plt) # Save plot
end
