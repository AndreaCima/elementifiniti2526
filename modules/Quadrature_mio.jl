# Author: Ivan Bioli (https://github.com/IvanBioli)

"""
    triarea(V1, V2, V3)

Calculate the area of a triangle given its vertices.

# Arguments
- `V1`: The first vertex of the triangle.
- `V2`: The second vertex of the triangle.
- `V3`: The third vertex of the triangle.

# Returns
- `area::Float64`: The area of the triangle.
"""
function triarea(V1, V2, V3)
    l1 = norm(V2 - V3);
    l2 = norm(V1 - V3);
    l3 = norm(V1 - V2);
    p = (l1 + l2 + l3) / 2;
    area = sqrt(p * (p - l1) * (p - l2) * (p - l3));
    return area
end

"""
    Q0(p, T, u)

Perform numerical integration using the Q0 quadrature rule (i.e., baricenter formula) over a mesh.
This quadrature rule has order 1.

# Arguments
- `p::Matrix`: The coordinates of the mesh nodes.
- `T::Matrix`: The connectivity matrix of the mesh elements.
- `u::Function`: The function to be integrated.

# Returns
- `I_approx::Float64`: The approximate integral of the function over the mesh.
"""
function Q0(p, T, u)
    I_approx = 0;
    for i in eachindex(axes(T, 2))
        V1 = p[:, T[1, i]];
        V2 = p[:, T[2, i]];
        V3 = p[:, T[3, i]];
        area = triarea(V1, V2, V3);
        bt = (V1 + V2 + V3)/3;
        I_approx += area * u(bt);

    end
    return I_approx
end

"""
    Q1(p, T, u)

Perform numerical integration using the Q1 quadrature rule (i.e., vertex formula) over a mesh.
This quadrature rule has order 1.

# Arguments
- `p::Matrix`: The coordinates of the mesh nodes.
- `T::Matrix`: The connectivity matrix of the mesh elements.
- `u::Function`: The function to be integrated.

# Returns
- `I_approx::Float64`: The approximate integral of the function over the mesh.
"""
function Q1(p, T, u)
    I_approx = 0;
    for i in eachindex(axes(T, 2))
        V1 = p[:, T[1, i]];
        V2 = p[:, T[2, i]];
        V3 = p[:, T[3, i]];
        area = triarea(V1, V2, V3);
        I_approx += area * (u(V1) + u(V2) + u(V3))/3;

    end
    return I_approx
end

"""
    Q2(p, T, u)

Perform numerical integration using the Q2 quadrature rule (i.e., midpoints rule) over a mesh.
This quadrature rule has order 2.

# Arguments
- `p::Matrix`: The coordinates of the mesh nodes.
- `T::Matrix`: The connectivity matrix of the mesh elements.
- `u::Function`: The function to be integrated.

# Returns
- `I_approx::Float64`: The approximate integral of the function over the mesh.
"""
function Q2(p, T, u)
    I_approx = 0;
    for i in eachindex(axes(T, 2))
        V1 = p[:, T[1, i]];
        V2 = p[:, T[2, i]];
        V3 = p[:, T[3, i]];
        area = triarea(V1, V2, V3);
        m1 = (V2 + V3) / 2;
        m2 = (V1 + V3) / 2;
        m3 = (V1 + V2) / 2;
        I_approx += area * (u(m1) + u(m2) + u(m3))/3;
    end
    return I_approx
end