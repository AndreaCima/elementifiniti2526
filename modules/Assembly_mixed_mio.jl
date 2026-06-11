# Author: Ivan Bioli (https://github.com/IvanBioli)
# Inspired by code written by Jochen Hinz (https://github.com/JochenHinz) for MATH-451 @ EPFL

using Memoize
using SparseArrays

"""
    initialize_assembly_mixed!(mesh::Mesh)

Initializes the assembly process for a mixed finite element method on the given `mesh`.

This function performs the following steps:
- Computes the element transformation matrices by calling `get_Bk!`.
- Calculates the determinants of the transformation matrices by calling `get_detBk!`.

# Arguments
- `mesh::Mesh`: The mesh data structure on which the assembly is to be initialized.

# Side Effects
Modifies the internal state of `mesh` by updating its transformation matrices and their determinants.

# See Also
- [`get_Bk!`](@ref)
- [`get_detBk!`](@ref)
"""
function initialize_assembly_mixed!(mesh::Mesh)
    get_Bk!(mesh)
    get_detBk!(mesh)
end

########################### GLOBAL ASSEMBLER ########################### 
"""
    assemble_global_mixed(mesh::Mesh, local_assembler!)

Assembles the global mixed finite element system matrices and right-hand side vector for a given mesh.

# Arguments
- `mesh::Mesh`: The mesh data structure containing information about elements, edges, and connectivity.
- `local_assembler!`: A function that computes the local element matrices and force vector for a given cell. It should have the signature `local_assembler!(Ae, Be, fe, mesh, cell_index)`.

# Returns
- `K`: The global system matrix assembled as a block matrix, where the upper-left block corresponds to the velocity matrix, the upper-right and lower-left blocks correspond to the coupling (pressure) matrices, and the lower-right block is a zero matrix.
- `b`: The global right-hand side vector, concatenating zeros for velocity DOFs and the assembled force vector for pressure DOFs.
"""
function assemble_global_mixed(mesh::Mesh, local_assembler!)

    ###########################################################################
    # matrice [ A  tB ]
    #         [ B   0 ]  è una matrice (#(lati totali) + N_tri) x (#(lati totali) + N_tri)
    # f = vettori degli integrali di f sui triangoli 
    ########################################################################### 
    T = mesh.T;
    p = mesh.p;
    E = mesh.elems2edges # matrice 3 x Ntri in cui ogni colonna contiene gli indici dei tre lati di un elemento 
    Ntri = size(T,2)
    Npoints = size(p,2)
    Nedges = size(mesh.edges2nodes, 2) # mesh.edges2nodes è una matrice 2 x Nedges 

    A_loc = zeros(3,3)
    B_loc = zeros(1,3)
    fe_loc = [0.0]

    #############################################################################
    # versione con matrici piene

    # A_glob = zeros(Nedges, Nedges)
    # B_glob = zeros(Ntri, Nedges)
    # O = zeros(Ntri, Ntri)
    # b = zeros(Nedges + Ntri)

    # for k = 1:Ntri
    #     local_assembler!(A_loc, B_loc, fe_loc, mesh, k)
    #     indici = E[:, k] # dobbiamo ragionare sui vertici
       
    #     b[Nedges + k] = fe_loc[1]
    #     A_glob[indici, indici] .+= A_loc
    #     B_glob[k,indici] = B_loc


    # end

    # println("-------------------------------")
    # println("lati: ", Nedges)
    # println("triangoli: ", Ntri)
    # println("rows A: ", size(A_glob, 1), "\t cols A: ", size(A_glob, 2))
    # println("rows B: ", size(B_glob, 1), "\t cols B: ", size(B_glob, 2))
    # # println("rows K: ", size(K, 1), "\t cols K: ", size(K, 2))
    # println("lunghezza b: ", length(b))
    # # println("rango di K: ", rank(K))
    # println("rango di A: ", rank(A_glob))
    # println("rango di B: ", rank(B_glob))
    # # println("nodi di Dirichlet: ", D)
    # println("-------------------------------")
    #######################################################################################

    # dim A = Nedges x Nedges
    rows_A = []
    cols_A = []
    data_A = Float64[]

    # dim B = Ntri x Nedges
    rows_B = []
    cols_B = []
    data_B = Float64[]

    # dim f = Ntri, a cui poi aggiungo all'inizio Nedges zeri
    rows_f = []
    data_f = Float64[]

    A_loc = zeros(3,3)
    B_loc = zeros(1, 3)
    f_loc = [0.0]


    for k in 1:Ntri
        local_assembler!(A_loc, B_loc, f_loc, mesh, k) 
        indices = E[:, k]
        append!(rows_B, fill(k, 3)) # indice di colonna, per ogni triangolo ci sono tre valori che poi calcolo sotto

        push!(rows_f, k)
        push!(data_f, f_loc[1]) # è uno scalare
    
        for i in 1:3
            i_glob = indices[i]
            push!(cols_B, i_glob)
            push!(data_B, B_loc[1, i])

            for j in 1:3
                j_glob = indices[j]
                push!(rows_A, i_glob)
                push!(cols_A, j_glob)
                push!(data_A, A_loc[i,j])
            end
            
        end    
    end

    A_glob = sparse(rows_A, cols_A, data_A, Nedges, Nedges)
    B_glob = sparse(rows_B, cols_B, data_B, Ntri, Nedges)
    F_glob = Matrix(sparse(rows_f, ones(size(rows_f)), data_f))
    
    # println("-------------------------------")
    # println(typeof(rows_B))
    # println(size(rows_B))
    # println(typeof(cols_B))
    # println(size(rows_B))
    # println(typeof(data_B))
    # println(size(B_loc))
    # println("-------------------------------")
    
    

    K = [A_glob transpose(B_glob); B_glob spzeros(Ntri, Ntri)]

    b = [zeros(Nedges); F_glob]





    return K, b


end


"""
    shapef_2D_RT0FE(quadrule::TriQuad)

Compute the Raviart-Thomas RT0 vector-valued shape functions on the reference triangle at the given quadrature points.

# Arguments
- `quadrule::TriQuad`: A quadrature rule object for triangles, containing the quadrature points as a 2×n matrix.

# Returns
- `shapef::Array{Float64,3}`: A 3D array of size (2, 3, n), where:
    - The first dimension (2) corresponds to the vector components (x and y).
    - The second dimension (3) corresponds to the three RT0 basis functions.
    - The third dimension (n) corresponds to the number of quadrature points.

# Details
The RT0 basis functions on the reference triangle are:
- `f1(x, y) = [x; y - 1]`
- `f2(x, y) = [x; y]`
- `f3(x, y) = [x - 1; y]`

The function evaluates these basis functions at each quadrature point and returns them in a single array for efficient use in finite element assembly.

# Memoization
The function is memoized to cache results for repeated calls with the same quadrature rule.
"""
# FIXME: PUT MEMOIZE BACK AFTER IMPLEMENTATION
# @memoize function shapef_2D_RT0FE(quadrule::TriQuad)
function shapef_2D_RT0FE(quadrule::TriQuad)
    PQ = quadrule.points
    n = size(PQ, 2)
    # è un tensore 2 x 3 x n_quad, ovvero (dimensione spazio) x (numero funzioni di base) x (numero punti di quadratura)
    shapef_RT = zeros(2, 3, n)

    shapef_RT[1, 1, :] = reshape(PQ[1, :], 1, 1, n)
    shapef_RT[1, 2, :] = reshape(PQ[1, :], 1, 1, n)
    shapef_RT[1, 3, :] = reshape(PQ[1, :] .- 1, 1, 1, n)

    shapef_RT[2, 1, :] = reshape(PQ[2, :] .- 1, 1, 1, n)
    shapef_RT[2, 2, :] = reshape(PQ[2, :], 1, 1, n)
    shapef_RT[2, 3, :] = reshape(PQ[2, :], 1, 1, n)

    return shapef_RT
end


"""
    divshapef_2D_RT0FE(quadrule::TriQuad)

Compute the divergence of the Raviart-Thomas RT0 vector-valued basis functions on a 2D triangle at the given quadrature points.

# Arguments
- `quadrule::TriQuad`: A quadrature rule object containing the quadrature points for integration over a triangle.

# Returns
- `divshapef::Array{Int,3}`: A `(1, 3, n)` array, where `n` is the number of quadrature points. Each entry contains the divergence of the three RT0 basis functions, which is constant and equal to 2 for each function.

# Notes
- The RT0 basis functions on the reference triangle are:
    - `f₁(x, y) = [x; y - 1]`
    - `f₂(x, y) = [x; y]`
    - `f₃(x, y) = [x - 1; y]`
- The divergence of each basis function is constant and equal to 2.
- The result is repeated for each quadrature point.
"""
# FIXME: PUT MEMOIZE BACK AFTER IMPLEMENTATION
# @memoize function divshapef_2D_RT0FE(quadrule::TriQuad)
function divshapef_2D_RT0FE(quadrule::TriQuad)
    PQ = quadrule.points
    n = size(PQ, 2)

    Δshapef = 2*ones(1, 3, n)

    return Δshapef
     
end

"""
    darcy_assemble_local_mixed!(Ae::Matrix, Be::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f, μ)

Assembles the local matrices and vector for the mixed finite element formulation of the Darcy problem on a single cell.

# Arguments
- `Ae::Matrix`: Local stiffness matrix (to be filled in-place).
- `Be::Matrix`: Local divergence matrix (to be filled in-place).
- `fe::Vector`: Local right-hand side vector (to be filled in-place).
- `mesh::Mesh`: Mesh data structure containing geometry and connectivity.
- `cell_index::Integer`: Index of the current cell in the mesh.
- `f`: Function representing the source term, evaluated as `f(x)` at a point `x`.
- `μ`: Function representing the permeability coefficient, evaluated as `μ(x)` at a point `x`.

# Description
This function computes the local contributions to the global system for the mixed finite element discretization of the Darcy problem. It uses RT0 (Raviart-Thomas of lowest order) basis functions and Piola transformation for mapping reference shape functions to the physical element. The function performs numerical integration using a quadrature rule, and assembles the local matrices and vector by looping over quadrature points and basis functions.

# Notes
- The function assumes that the shape functions are oriented consistently using `elems2orientation`.
- The local matrices and vector are reset to zero at the beginning of the function.
- The function modifies `Ae`, `Be`, and `fe` in-place and returns them for convenience.

# Returns
- `(Ae, Be, fe)`: The assembled local stiffness matrix, divergence matrix, and right-hand side vector.
"""
########################### DARCY PROBLEM ###########################
function darcy_assemble_local_mixed!(Ae::Matrix, Be::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f, μ)
    B, a = get_Bk!(mesh)
    detB = get_detBk!(mesh)
    invB = get_invBk!(mesh)
    Bk = B[:, :, cell_index]
    ak = a[:, cell_index]
    detBk = detB[cell_index]
    invBk = invB[:, :, cell_index]

    quadrule = Q2_ref

    pe = Bk * quadrule.points .+ ak # punti trasformati 

    shapef = shapef_2D_RT0FE(quadrule)
    divshapef = divshapef_2D_RT0FE(quadrule)

    fill!(Ae, 0.0)
    fill!(Be, 0.0)
    fill!(fe, 0.0)

    for p in eachindex(axes(pe, 2)) # sommo sui punti di quadratura
        dΩ = quadrule.weights[p]* detBk
        fe .-= f(pe[:, p]) * dΩ

        for i = 1:3
            sign_i = mesh.elems2orientation[i, cell_index]

            phi_i = (sign_i/detBk) * Bk * shapef[:, i, p]
            div_phi_i = (sign_i/detBk) * divshapef[1, i, p]

            Be[1, i] -= div_phi_i * dΩ 

            for j = 1:3
                sign_j = mesh.elems2orientation[j, cell_index]
                phi_j = (sign_j/detBk) * Bk * shapef[:, j, p]

                Ae[i, j] += (μ(pe[:, p]) * phi_i) ⋅ phi_j * dΩ
            end
        end
    end
    # println("-------------------------------")

    # println("rows A: ", size(Ae, 1), "\t cols A: ", size(Ae, 2))
    # println("rows B: ", size(Be, 1), "\t cols B: ", size(Be, 2))
    # println("lunghezza F: ",  length(fe))
    
    # println("rango di A: ", rank(Ae))
    # println("rango di B: ", rank(Be))
    # println("-------------------------------")
    return Ae, Be, fe

end


########################## DEFINE FUNCTIONS TO COMPUTE ERROR ##########################
"""
    L2error_mixed_p(p::Function, ph::Vector, mesh::Mesh, ref_quad::TriQuad) -> Float64

Compute the L² error between an exact function `p` and its discrete approximation `ph` over a given mesh using quadrature.

# Arguments
- `p::Function`: The exact function to be evaluated.
- `ph::Vector`: Vector of discrete approximations (one per element).
- `mesh::Mesh`: The mesh structure containing element connectivity and geometry.
- `ref_quad::TriQuad`: Quadrature rule on the reference triangle, providing points and weights.

# Returns
- `Float64`: The computed L² error over the mesh.
"""
function L2error_mixed_p(p::Function, ph::Vector, mesh::Mesh, ref_quad::TriQuad)
    WQ = ref_quad.weights;
    PQ = ref_quad.points;
    Bk, ak = get_Bk!(mesh);
    detBk = get_detBk!(mesh);
    T = mesh.T;
    Npoints = size(mesh.p, 2)
    Ntri = size(mesh.T, 2)

    L2error = 0.0 

    for t in 1:Ntri # sommo gli errori sui singoli triangoli

        a = ak[:, t]
        B = Bk[:,:, t]
        detB = detBk[t]

        # println(ph)
        
        ph_t = ph[t]

        # devo calcolare l'area di |p - ph|^2 in cui ph è costante per ogni triangolo

        pe = B * PQ .+ a # punti trasformati
        p_ev = p.(eachcol(pe))

        L2error += (p_ev .- ph_t).^2 ⋅ WQ * detB
    end

    L2error = sqrt(L2error)
    return L2error



end

"""
    H1diverror_mixed_u(u::Function, divu::Function, uh::Vector, mesh::Mesh, ref_quad::TriQuad) -> Float64

Compute the H(div) error norm between the exact solution `(u, divu)` and the finite element solution `uh` 
for a mixed finite element method using Raviart-Thomas (RT0) elements on a triangular mesh.

# Arguments
- `u::Function`: The exact vector-valued solution function, accepting a point and returning a 2D vector.
- `divu::Function`: The exact divergence of the solution, accepting a point and returning a scalar.
- `uh::Vector`: The vector of degrees of freedom for the finite element solution (RT0 coefficients).
- `mesh::Mesh`: The mesh data structure, containing element connectivity, geometry, and orientation.
- `ref_quad::TriQuad`: Quadrature rule on the reference triangle, providing points and weights.

# Returns
- `Float64`: The combined H(div) error norm, i.e., `sqrt(∫|u - uh|^2) + sqrt(∫|divu - div(uh)|^2)` over the domain.
"""
function H1diverror_mixed_u(u::Function, divu::Function, uh::Vector, mesh::Mesh, ref_quad::TriQuad)
    ###########################################################################
    ############################ ADD CODE HERE ################################
    ########################################################################### 
end