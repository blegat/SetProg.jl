using Polyhedra
using SumOfSquares

"""
    struct PolynomialSublevelSetAtOrigin{T, P<:AbstractPolynomial{T}}
        degree::Int
        p::P
    end

Set ``\\{\\, x \\mid p(x) \\le 1 \\,\\}`` where `p` is a homogeneous polynomial
of degree `degree`.
"""
struct PolynomialSublevelSetAtOrigin{T} <: AbstractSet{T}
    degree::Int
    p::GramMatrix{T, DynamicPolynomials.Monomial{true},
                     DynamicPolynomials.MonomialVector{true}}
end

function space_variables(set::PolynomialSublevelSetAtOrigin)
    return variables(set.p)
end

"""
    struct ConvexPolynomialSublevelSetAtOrigin{T, P<:AbstractPolynomial{T}}
        degree::Int
        p::P
    end

Set ``\\{\\, x \\mid p(x) \\le 1 \\,\\}`` where `p` is a homogeneous polynomial
of degree `degree`.
"""
struct ConvexPolynomialSublevelSetAtOrigin{T} <: AbstractSet{T}
    degree::Int
    p::GramMatrix{T, DynamicPolynomials.Monomial{true},
                     DynamicPolynomials.MonomialVector{true}}
    convexity_proof::Union{Nothing, SumOfSquares.SymMatrix{T}} # may be nothing after applying LinearMap
end

function space_variables(set::ConvexPolynomialSublevelSetAtOrigin)
    return variables(set.p)
end
function dimension(set::ConvexPolynomialSublevelSetAtOrigin)
    return nvariables(set.p)
end
function gauge1(set::ConvexPolynomialSublevelSetAtOrigin)
    return set.p
end
function zero_eliminate(set::ConvexPolynomialSublevelSetAtOrigin, I)
    vars = space_variables(set)[I]
    K = findall(mono -> all(var -> iszero(degree(mono, var)), vars),
                set.p.x)
    M = set.p.Q[K, K]
    Q = SumOfSquares.SymMatrix([M[i, j] for j in 1:length(K) for i in 1:j],
                               length(K))
    monos = set.p.x[K]
    J = setdiff(1:dimension(set), I)
    monos = DynamicPolynomials.MonomialVector(monos.vars[J],
                                             Vector{Int}[z[J] for z in monos.Z])
    p = SumOfSquares.GramMatrix(Q, monos)
    return ConvexPolynomialSublevelSetAtOrigin(set.degree, p, nothing)
end

convexity_proof(set::ConvexPolynomialSublevelSetAtOrigin) = set.convexity_proof

function scaling_function(set::Union{PolynomialSublevelSetAtOrigin,
                                     ConvexPolynomialSublevelSetAtOrigin})
    # We convert the GramMatrix to a polynomial to avoid having to do the
    # conversion for every substitution.
    p = polynomial(set.p)
    vars = variables(p)
    @assert length(vars) == 2
    vx, vy = vars
    return (x, y) -> p(vx => x, vy => y)^(1 / set.degree)
end

"""
    struct ConvexPolynomialSet{T} <: AbstractSet{T}
        degree::Int
        q::GramMatrix{T, DynamicPolynomials.Monomial{true},
                         DynamicPolynomials.MonomialVector{true}}
        z::SpaceVariable
        x::Vector{SpaceVariable}
    end

Set ``\\{\\, (z, x) \\mid p(z, x) \\le 0 \\,\\}`` or
``H \\{\\, (z, x) \\mid q(z, x) \\le z^{\\texttt{degree}} \\,\\}`` where `p` and
`q` are homogeneous polynomials of degree `degree` and `H` is a householder
matrix.
"""
struct ConvexPolynomialSet{T} <: AbstractSet{T}
    degree::Int
    q::GramMatrix{T, DynamicPolynomials.Monomial{true},
                     DynamicPolynomials.MonomialVector{true}}
    z::SpaceVariable
    x::Vector{SpaceVariable}
end
perspective_gauge0(set) = set.q - set.z^set.degree
perspective_variable(set::ConvexPolynomialSet) = set.z
space_variables(set::ConvexPolynomialSet) = set.x
function gauge1(set::ConvexPolynomialSet{T}) where T
    return subs(set.q, perspective_variable(set) => one(T))
end
