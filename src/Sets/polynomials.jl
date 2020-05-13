using Polyhedra
using SumOfSquares

"""
    struct PolynomialSublevelSetAtOrigin{T, B, U} <: AbstractSet{U}
        degree::Int
        p::GramMatrix{T, B, U}
    end

Set ``\\{\\, x \\mid p(x) \\le 1 \\,\\}`` where `p` is a homogeneous polynomial
of degree `degree`.
"""
struct PolynomialSublevelSetAtOrigin{T, B, U} <: AbstractSet{U}
    degree::Int
    p::GramMatrix{T, B, U}
end

function space_variables(set::PolynomialSublevelSetAtOrigin)
    return variables(set.p)
end

"""
    struct ConvexPolynomialSublevelSetAtOrigin{T, B, U} <: AbstractSet{U}
        degree::Int
        p::GramMatrix{T, B, U}
        convexity_proof::Union{Nothing, SumOfSquares.SymMatrix{T}} # may be nothing after applying LinearMap
    end

Set ``\\{\\, x \\mid p(x) \\le 1 \\,\\}`` where `p` is a homogeneous polynomial
of degree `degree`.
"""
struct ConvexPolynomialSublevelSetAtOrigin{T, B, U} <: AbstractSet{U}
    degree::Int
    p::GramMatrix{T, B, U}
    convexity_proof::Union{Nothing, SumOfSquares.SymMatrix{T}} # may be nothing after applying LinearMap
end
function ConvexPolynomialSublevelSetAtOrigin(
    degree::Int,
    p::GramMatrix{T, B, U},
    convexity_proof::SumOfSquares.SymMatrix{T}) where {T, B, U}
    return ConvexPolynomialSublevelSetAtOrigin{T, B, U}(degree, p, convexity_proof)
end
function ConvexPolynomialSublevelSetAtOrigin(
    degree::Int,
    p::GramMatrix{S},
    convexity_proof::SumOfSquares.SymMatrix{T}) where {S, T}
    V = promote_type(S, T)
    _convert(mat) = SumOfSquares.SymMatrix(convert(Vector{U}, mat.Q), mat.n)
    return ConvexPolynomialSublevelSetAtOrigin(
        degree, GramMatrix(_convert(p.Q), p.basis), _convert(convexity_proof))
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
                set.p.basis.monomials)
    M = set.p.Q[K, K]
    Q = SumOfSquares.SymMatrix([M[i, j] for j in 1:length(K) for i in 1:j],
                               length(K))
    monos = set.p.basis.monomials[K]
    J = setdiff(1:dimension(set), I)
    monos = DynamicPolynomials.MonomialVector(monos.vars[J],
                                             Vector{Int}[z[J] for z in monos.Z])
    p = SumOfSquares.GramMatrix(Q, MB.MonomialBasis(monos))
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
    struct ConvexPolynomialSet{T, U} <: AbstractSet{U}
        degree::Int
        q::GramMatrix{T, MonoBasis, u}
        z::SpaceVariable
        x::Vector{SpaceVariable}
    end

Set ``\\{\\, (z, x) \\mid p(z, x) \\le 0 \\,\\}`` or
``H \\{\\, (z, x) \\mid q(z, x) \\le z^{\\texttt{degree}} \\,\\}`` where `p` and
`q` are homogeneous polynomials of degree `degree` and `H` is a householder
matrix.
"""
struct ConvexPolynomialSet{T, B, U} <: AbstractSet{U}
    degree::Int
    q::GramMatrix{T, B, U}
    z::SpaceVariable
    x::Vector{SpaceVariable}
end
perspective_gauge0(set) = set.q - set.z^set.degree
perspective_variable(set::ConvexPolynomialSet) = set.z
space_variables(set::ConvexPolynomialSet) = set.x
function gauge1(set::ConvexPolynomialSet{T}) where T
    return subs(set.q, perspective_variable(set) => one(T))
end
