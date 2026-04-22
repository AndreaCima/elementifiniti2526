# Author: Ivan Bioli (https://github.com/IvanBioli)
# Inspired by code written by Jochen Hinz (https://github.com/JochenHinz) for MATH-451 @ EPFL

using Memoize
using SparseArrays

"""
    initialize_assembly!(mesh::Mesh)

Initialize the assembly process for the given mesh by computing the necessary geometric quantities.

# Arguments
- `mesh::Mesh`: The mesh object for which the assembly is initialized.
"""
function initialize_assembly!(mesh::Mesh)
    get_Bk!(mesh)
    get_detBk!(mesh)
    get_invBk!(mesh)
end

########################### GLOBAL ASSEMBLER ########################### 
"""
    assemble_global(mesh::Mesh, local_assembler!)

Assemble the global stiffness matrix and force vector for the given mesh using the provided local assembler function.

# Arguments
- `mesh::Mesh`: The mesh object.
- `local_assembler!`: A function that assembles the local stiffness matrix and force vector.

# Returns
- `K::SparseMatrixCSC`: The global stiffness matrix.
- `f::Vector`: The global force vector.
"""
function assemble_global(mesh::Mesh, local_assembler!)
    T = mesh.T
    p = mesh.p
    n_points = size(p,2)
    n_triangles = size(T, 2)
    rows = []
    cols = []
    data = Float64[]
    
    rows_f = []
    data_f = Float64[]

    A_loc = zeros(3,3)
    f_loc = zeros(3)
    for k in 1:n_triangles
        local_assembler!(A_loc, f_loc, mesh, k) 
        indices = T[:, k]
    
        for i in 1:3
            i_glob = indices[i]
            for j in 1:3
                j_glob = indices[j]
                push!(rows, i_glob)
                push!(cols, j_glob)
                push!(data, A_loc[i,j])
            end
            push!(rows_f, i_glob)
            push!(data_f, f_loc[i])
        end    
    end
    A_glob = sparse(rows, cols, data, n_points, n_points)
    F_glob = Matrix(sparse(rows_f, ones(size(rows_f)), data_f))
    return A_glob, F_glob
end

"""
    impose_dirichlet(A, b, g, mesh)

Impose Dirichlet boundary conditions on the system.

# Arguments
- `A`: The global stiffness matrix.
- `b`: The global force vector.
- `g`: The Dirichlet boundary condition function.
- `mesh::Mesh`: The mesh object.

# Returns
- `A_cond`: The modified stiffness matrix with Dirichlet conditions imposed.
- `b_cond`: The modified force vector with Dirichlet conditions imposed.
- `uh`: The solution vector with Dirichlet conditions applied.
"""
function impose_dirichlet(A, b, g, mesh)
    F = mesh.freedofs
    D = mesh.dirichletdofs
    T = mesh.T
    p = mesh.p

    pD = p[:, T[D]] # coordinate dei punti di bordo

    # pD_g = pD

    # for col in eachindex(pD)
    #     pD_g[:, col] = g.(pD[:, col])
    # end

    A_cond = A[F, F]
    b_cond = b[F] - A[F, D]*g.(eachcol(pD))

    uh = zeros(size(p, 2))
    u_F = A_cond \ b_cond
    u_D = g.(eachcol(pD))

    uh[F] = u_F
    uh[D] = u_D

    return A_cond, b_cond, uh

end


########################### POISSON PROBLEM ###########################
"""
    shapef_2DLFE(quadrule::TriQuad)

Compute the shape functions for the Poisson problem.

# Arguments
- `quadrule::TriQuad`: The quadrature rule.

# Returns
- `shapef`: The shape functions evaluated at the quadrature points.
"""
function shapef_2DLFE(quadrule::TriQuad) # aggiungere @memoize prima di function una volta controllato che il codice funzioni
    # restituisce una matrice 3 x n_quad in cui nella colonna i valoto le tre funzioni di base nell'i-esimo punto di quadratura
    PQ = quadrule.points
    shapef = zeros(3, size(PQ, 2))
    shapef[1, :] = 1 .- PQ[1, :] .- PQ[2, :] # prima funzione di base
    shapef[2, :] = PQ[1, :]  # seconda funzione di base
    shapef[3, :] = PQ[2, :]  # terza funzione di base
    return shapef
end

"""
    ∇shapef_2DLFE(quadrule::TriQuad)

Compute the gradients of the shape functions for the Poisson problem.

# Arguments
- `quadrule::TriQuad`: The quadrature rule.

# Returns
- `∇shapef`: The gradients of the shape functions evaluated at the quadrature points.
"""
function ∇shapef_2DLFE(quadrule::TriQuad) # aggiungere @memoize prima di function una volta controllato che il codice funzioni
    # il risultato dipende solo dal numero di punti della quadratura e non dalla loro espressione
    # restituisce un tensore 2 x n_funzioni di base(3) x numero punti di quadratura 
    PQ = quadrule.points
    face = [-1 1 0; -1 0 1]
    ∇shapef = repeat(face, 1, 1, size(PQ, 2))
    return ∇shapef
end

"""
    poisson_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f)

Assemble the local stiffness matrix and force vector for the Poisson problem.

# Arguments
- `Ke::Matrix`: The local stiffness matrix to be assembled.
- `fe::Vector`: The local force vector to be assembled.
- `mesh::Mesh`: The mesh object.
- `cell_index::Integer`: The index of the current cell.
- `f`: The source term function.

# Returns
- `Ke`: The assembled local stiffness matrix.
- `fe`: The assembled local force vector.
"""
function poisson_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f)
    B, a = get_Bk!(mesh)
    detB = get_detBk!(mesh)
    invB = get_invBk!(mesh)

    Bk = B[:, :, cell_index]
    ak = a[:, cell_index]
    detBk = detB[cell_index]
    invBk = invB[:, :, cell_index]

    Q_matrix = Q0_ref
    phi_grad = ∇shapef_2DLFE(Q_matrix) # per la stiffness
    # phi_val_matrix = shapef_2DLFE(Q_matrix)
    W_matrix = Q_matrix.weights

    Q_vector = Q2_ref
    phi_val_vector = shapef_2DLFE(Q_vector)
    q_val = Q_vector.points # per valutare la funzione f
    W_vector = Q_vector.weights
    
    # Inizializzo le nuove Ke e f a zero
    fill!(Ke, 0)
    fill!(fe, 0)

    for i=1:3
        for j=1:3
            # fisso la funzione di base e prendo tutte le coordinate dei punti e tutti i punti di quadratura
            fattore_∇phi_i = transpose(invBk)*phi_grad[:, i, :]
            fattore_∇phi_j = transpose(invBk)*phi_grad[:, j, :]

            for s in 1:size(Q_matrix.points, 2) # sommo su tutti i punti di quadratura
                Ke[i, j] += dot(fattore_∇phi_i[:, s], fattore_∇phi_j[:, s])*detBk*W_matrix[s]
            end
        
        end
        for t in 1:size(Q_vector.points, 2) # sommo su tutti i nodi di quadratura
            fe[i] += f(Bk*q_val[:, t] + ak)*phi_val_vector[i, t] * detBk * W_vector[t]
        end

    end    


end