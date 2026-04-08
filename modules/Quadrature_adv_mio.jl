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
        B = Bk[i]
        a = ak[i]
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