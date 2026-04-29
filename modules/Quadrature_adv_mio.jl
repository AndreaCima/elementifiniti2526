# Author: Ivan Bioli (https://github.com/IvanBioli)
# Inspired by code written by Jochen Hinz (https://github.com/JochenHinz) for MATH-451 @ EPFL

"""
    struct TriQuad

A structure representing a triangular quadrature rule.

# Fields
- `name::String`: The name of the quadrature rule.
- `order::Integer`: The order of the quadrature rule.
- `points::Matrix`: The quadrature points.
- `weights::Array`: The quadrature weights.
"""
struct TriQuad
    name::String
    order::Integer
    points::Matrix
    weights::Array
end

M0 = zeros(2, 1)
M0[1, 1] = 1/3
M0[2, 1] = 1/3

Q0_ref = TriQuad("Q0", 2, M0, [1/2])
Q1_ref = TriQuad("Q1", 2, [0 1 0; 0 0 1], [1/6 1/6 1/6])
Q2_ref = TriQuad("Q2", 3, [1/2 0 1/2; 0 1/2 1/2], [1/6 1/6 1/6])

"""
    Quadrature(u, mesh::Mesh, ref_quad::TriQuad)

Perform numerical integration of a function over a mesh using a given quadrature rule.

# Arguments
- `u`: The function to be integrated.
- `mesh::Mesh`: The mesh object.
- `ref_quad::TriQuad`: The reference quadrature rule.

# Returns
- `I_approx::Float64`: The approximate integral of the function over the mesh.
"""
function Quadrature(u, mesh::Mesh, ref_quad::TriQuad)
    I = 0.0
    WQ = ref_quad.weights
    PQ = ref_quad.points
    q = length(WQ)
    Bk, ak = get_Bk!(mesh)
    detBk = get_detBk!(mesh)
    Tri = mesh.T

    for i in eachindex(axes(Tri, 2))
        T = 0
        B = Bk[:, :, i]
        a = ak[:, i]
        d = detBk[i]
        for j in 1:q
            p = PQ[:, j]
            w = WQ[j]
            pT = a + B * p

            T += w * u(pT)
        end
        I += d*T
    end
    return I
end

# Evaluation of a function
"""
    eval_u(u::Function, points_elem::Matrix, mesh::Mesh, tri_idx::Integer, quadrule::TriQuad)

Evaluate a function at given points within an element.

# Arguments
- `u::Function`: The function to be evaluated.
- `points_elem::Matrix`: The points at which to evaluate the function.
- `mesh::Mesh`: The mesh object (ignored).
- `tri_idx::Integer`: The index of the current element (ignored).
- `quadrule::TriQuad`: The quadrature rule (ignored).

# Returns
- `u_evals::Matrix`: The evaluated function values at the given points.
"""
function eval_u(u::Function, points_elem::Matrix, mesh::Mesh, tri_idx::Integer, quadrule::TriQuad)
    return u.(eachcol(points_elem))
end

"""
    eval_u(uh::Vector, points_elem::Matrix, mesh::Mesh, tri_idx::Integer, quadrule::TriQuad)

Evaluate a linear finite element solution at given quadrature points within an element.

# Arguments
- `uh::Vector`: The finite element solution vector.
- `points_elem::Matrix`: The points at which to evaluate the solution (ignored).
- `mesh::Mesh`: The mesh object.
- `tri_idx::Integer`: The index of the current element.
- `quadrule::TriQuad`: The quadrature rule.

# Returns
- `uh_evals::Matrix`: The evaluated solution values at the given points.
"""
function eval_u(uh::Vector, points_elem::Matrix, mesh::Mesh, tri_idx::Integer, quadrule::TriQuad)
    ###########################################################################
    ####################### PUT YOUR CODE HERE ################################
    ###########################################################################
end

"""
    L2error(u::Function, uh::Vector, mesh::Mesh, ref_quad::TriQuad)

Compute the L2 error between a function and a finite element solution over a mesh.

# Arguments
- `u::Function`: The exact solution function.
- `uh::Vector`: The finite element solution vector.
- `mesh::Mesh`: The mesh object.
- `ref_quad::TriQuad`: The reference quadrature rule.

# Returns
- `L2_error::Float64`: The L2 error between the exact solution and the finite element solution.
"""
function L2error(u::Function, uh::Vector, mesh::Mesh, ref_quad::TriQuad)

    WQ = ref_quad.weights;
    PQ = ref_quad.points;
    Bk, ak = get_Bk!(mesh);
    detBk = get_detBk!(mesh);
    T = mesh.T;

    # valutazioni delle funzioni di base sui punti di quadratura di riferimento
    φ_val = shapef_2DLFE(ref_quad);

    L2error = 0.0 

    for t in size(T,2)

        a = ak[:, t]
        B = Bk[:,:, t]
        detB = detBk[t]
        
        uh_t = uh[T[:,t]]

        # ciclo sui punti di quadratura
        for j = 1:size(WQ,2)
            q = PQ[:, j] # j-esimo punto di quadratura

            
            φ_j =  φ_val[:,j] 
            uh_val_j = dot(uh_t, φ_j) 


            L2error += detB * WQ[j] * (abs(   u( B*q + a) -  uh_val_j ))^2 
        end 


    end

    L2error = sqrt(L2error);
    return L2error


end

"""
    H1semierror(∇u::Function, uh::Vector, mesh::Mesh, ref_quad::TriQuad)

Compute the H1 semi-norm error between the gradient of a function and a finite element solution over a mesh.

# Arguments
- `∇u::Function`: The gradient of the exact solution function.
- `uh::Vector`: The finite element solution vector.
- `mesh::Mesh`: The mesh object.
- `ref_quad::TriQuad`: The reference quadrature rule.

# Returns
- `H1_semi_error::Float64`: The H1 semi-norm error between the gradient of the exact solution and the finite element solution.
"""
function H1semierror(∇u::Function, uh::Vector, mesh::Mesh, ref_quad::TriQuad)

    WQ = ref_quad.weights;
    PQ = ref_quad.points;
    Bk, ak = get_Bk!(mesh);
    detBk = get_detBk!(mesh);
    invBk = get_invBk!(mesh);
    T = mesh.T;

    # valutazioni dei gradienti delle funzioni di base sui punti di quadratura di riferimento
    ∇φ_val = ∇shapef_2DLFE(ref_quad)

    H1semierror = 0.0
    # Ciclo sui triangoli (calcola gli integrali localmente e li somma)
    for t in size(T,2)

        # richiama le quantità sull'elemento t fissato
        a = ak[:, t]
        B = Bk[:,:, t]
        invB = invBk[:, :, t]
        detB = detBk[t]
        
        # valori della funzione approssimata uh sui tre vertici del triangolo t fissato
        uh_t = uh[T[:,t]]

        # ciclo sui punti di quadratura
        for j = 1:size(WQ,2)
            q = PQ[:, j] # j-esimo punto di quadratura sul triangolo di riferimento

            ∇uh_j = zeros(2)
            for i =1:3
                ∇uh_j += uh_t[i] * ( (invB)'* ∇φ_val[:,i,j]   )
            end 



            H1semierror += detB* WQ[j] * (norm(   ∇u( B*q + a) -  ∇uh_j ))^2 
        end 


    end

    H1semierror = sqrt(H1semierror);
    return H1semierror

end