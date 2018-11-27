using Polyhedra
using SumOfSquares

"""
    struct DualPolynomialSet{T} <: AbstractSet{T}
        degree::Int
        p::DynamicPolynomials.Polynomial{true, T}
    end

Sets whose perspective-dual is ``\\{\\, (z, x) \\mid p(z, x) \\le 0 \\,\\}``
where `p` is a homogeneous polynomial of degree `degree`.
"""
struct DualPolynomialSet{T} <: AbstractSet{T}
    degree::Int
    p::DynamicPolynomials.Polynomial{true, T}
end

"""
    struct PolarConvexPolynomialSublevelSetAtOrigin{T} <: AbstractSet{T}
        degree::Int
        p::MatPolynomial{T, DynamicPolynomials.Monomial{true},
                         DynamicPolynomials.MonomialVector{true}}
        convexity_proof::Union{Nothing, SumOfSquares.SymMatrix{T}} # may be nothing after applying LinearMap
    end

Set whose polar is ``\\{\\, x \\mid p(x) \\le 1 \\,\\}`` where `p` is a
homogeneous polynomial of degree `degree`.
"""
struct PolarConvexPolynomialSublevelSetAtOrigin{T} <: AbstractSet{T}
    degree::Int
    p::MatPolynomial{T, DynamicPolynomials.Monomial{true},
                     DynamicPolynomials.MonomialVector{true}}
    convexity_proof::Union{Nothing, SumOfSquares.SymMatrix{T}} # may be nothing after applying LinearMap
end

"""
    struct DualConvexPolynomialCone{T, U} <: AbstractSet{T}
        degree::Int
        p::MatPolynomial{T, DynamicPolynomials.Monomial{true},
                            DynamicPolynomials.MonomialVector{true}}
        q::DynamicPolynomials.Polynomial{true, U}
    end

Set whose dual is ``\\{\\, (z, x) \\mid p(z, x) \\le 0 \\,\\}`` or
``\\{\\, (z, x) \\mid q(z, x) \\le z^{\\texttt{degree}} \\,\\}`` where `p` and
`q` are homogeneous polynomials of degree `degree`.
"""
struct DualConvexPolynomialCone{T, U} <: AbstractSet{T}
    degree::Int
    q::MatPolynomial{T, DynamicPolynomials.Monomial{true},
                        DynamicPolynomials.MonomialVector{true}}
    p::DynamicPolynomials.Polynomial{true, U}
    z::DynamicPolynomials.PolyVar{true}
    x::Vector{DynamicPolynomials.PolyVar{true}}
end
function DualConvexPolynomialCone(degree::Integer, q::MatPolynomial, z, x)
    return DualConvexPolynomialCone(degree, q, q - z^degree, z, x)
end
dimension(d::DualConvexPolynomialCone) = length(d.x)

"""
    dual_contour(f::Function, nhalfspaces::Int, T::Type)

Return a polytope of `nhalfspaces` halfspaces defined by normal vectors of
equally spaced angles for the polar of the 1-sublevel set of the homogeneous
function `f(x, y)`.
"""
function dual_contour(f::Function, nhalfspaces::Int, T::Type)
    αs = range(0, stop=2π, length=nhalfspaces)
    h = hrep(Polyhedra.HyperPlane{T, Vector{T}}[],
             Polyhedra.HalfSpace{T, Vector{T}}[], d=2)
    for (i, α) in enumerate(range(0, stop=2π - 2π/nhalfspaces, length=nhalfspaces))
        a = cos(α)
        b = sin(α)
        r = f(a, b)
        # f is homogeneous so f(a/r, b/r) = 1 so the halfspace is
        # a*x/r + b*y/r ≤ 1 or equivalently a*x + b*y ≤ r
        intersect!(h, HalfSpace([a, b], r))
    end
    return polyhedron(h)
end

@recipe function f(set::PolarConvexPolynomialSublevelSetAtOrigin{T}; npoints=64) where T
    seriestype --> :shape
    legend --> false
    dual_contour(scaling_function(set), npoints, T)
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
    p::MatPolynomial{T, DynamicPolynomials.Monomial{true},
                     DynamicPolynomials.MonomialVector{true}}
    convexity_proof::Union{Nothing, SumOfSquares.SymMatrix{T}} # may be nothing after applying LinearMap
end

@recipe function f(set::ConvexPolynomialSublevelSetAtOrigin; npoints=64)
    seriestype --> :shape
    legend --> false
    primal_contour(scaling_function(set), npoints)
end

function scaling_function(set::Union{PolarConvexPolynomialSublevelSetAtOrigin,
                                     ConvexPolynomialSublevelSetAtOrigin})
    # We convert the MatPolynomial to a polynomial to avoid having to do the
    # conversion for every substitution.
    p = polynomial(set.p)
    vars = variables(p)
    @assert length(vars) == 2
    vx, vy = vars
    return (x, y) -> p(vx => x, vy => y)^(1 / set.degree)
end

function polar(set::ConvexPolynomialSublevelSetAtOrigin)
    return PolarConvexPolynomialSublevelSetAtOrigin( set.degree, set.p,
                                              set.convexity_proof)
end
function polar(set::PolarConvexPolynomialSublevelSetAtOrigin)
    return ConvexPolynomialSublevelSetAtOrigin(set.degree, set.p, set.convexity_proof)
end
