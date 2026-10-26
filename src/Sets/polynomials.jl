using Polyhedra

"""
    struct PolySet{T, B, U, P} <: AbstractSet{U}
        degree::Int
        p::P
    end

Set ``\\{\\, x \\mid p(x) \\le 1 \\,\\}`` where `p` is a polynomial
of degree at most `degree`.
"""
struct PolySet{T, B, U, P} <: AbstractSet{U}
    degree::Int
    p::P
end

perspective_variable(::PolySet) = nothing

"""
    struct ConvexPolySet{T, B, U, P, C} <: AbstractSet{U}
        degree::Int
        p::P
        convexity_proof::C # may be nothing after applying LinearMap
    end

Set ``\\{\\, x \\mid p(x) \\le 1 \\,\\}`` where `p` is a homogeneous polynomial
of degree `degree`.
"""
struct ConvexPolySet{T, B, U, P, C} <: AbstractSet{U}
    degree::Int
    p::P
    convexity_proof::C # may be nothing after applying LinearMap
end

perspective_variable(::ConvexPolySet) = nothing
function gauge1(set::ConvexPolySet)
    return set.p
end

convexity_proof(set::ConvexPolySet) = set.convexity_proof

"""
    struct ConvexPolynomialSet{T, B, U, P, V} <: AbstractSet{U}
        degree::Int
        q::P
        z::V
        x::Vector{V}
    end

Set ``\\{\\, (z, x) \\mid p(z, x) \\le 0 \\,\\}`` or
``H \\{\\, (z, x) \\mid q(z, x) \\le z^{\\texttt{degree}} \\,\\}`` where `p` and
`q` are homogeneous polynomials of degree `degree` and `H` is a householder
matrix.
"""
struct ConvexPolynomialSet{T, B, U, P, V} <: AbstractSet{U}
    degree::Int
    q::P
    z::V
    x::Vector{V}
end
perspective_gauge0(set) = set.q - set.z^set.degree
perspective_variable(set::ConvexPolynomialSet) = set.z
space_variables(set::ConvexPolynomialSet) = set.x
