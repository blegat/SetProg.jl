using Polyhedra
using SumOfSquares

"""
    struct PolarPolynomialSublevelSet{T, P<:AbstractPolynomial{T}}
        degree::Int
        p::P
    end

Set whose polar is ``\\{\\, x \\mid p(x) \\le 1 \\,\\}`` where `p` is a
homogeneous polynomial of degree `degree`.
"""
struct PolarPolynomialSublevelSetAtOrigin{T} <: AbstractSet{T}
    degree::Int
    p::MatPolynomial{T, DynamicPolynomials.Monomial{true},
                     DynamicPolynomials.MonomialVector{true}}
    convexity_proof::Union{Nothing, SumOfSquares.SymMatrix{T}}
end

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

@recipe function f(set::PolarPolynomialSublevelSetAtOrigin{T}; npoints=64) where T
    seriestype --> :shape
    legend --> false
    dual_contour(scaling_function(set), npoints, T)
end


"""
    struct PolarPolynomialSublevelSetAtOrigin{T, P<:AbstractPolynomial{T}}
        degree::Int
        p::P
    end

Set ``\\{\\, x \\mid p(x) \\le 1 \\,\\}`` where `p` is a homogeneous polynomial
of degree `degree`.
"""
struct PolynomialSublevelSetAtOrigin{T} <: AbstractSet{T}
    degree::Int
    p::MatPolynomial{T, DynamicPolynomials.Monomial{true},
                     DynamicPolynomials.MonomialVector{true}}
    convexity_proof::Union{Nothing, SumOfSquares.SymMatrix{T}}
end

@recipe function f(set::PolynomialSublevelSetAtOrigin; npoints=64)
    seriestype --> :shape
    legend --> false
    primal_contour(scaling_function(set), npoints)
end

function scaling_function(set::Union{PolarPolynomialSublevelSetAtOrigin,
                                     PolynomialSublevelSetAtOrigin})
    # We convert the MatPolynomial to a polynomial to avoid having to do the
    # conversion for every substitution.
    p = polynomial(set.p)
    vars = variables(p)
    @assert length(vars) == 2
    vx, vy = vars
    return (x, y) -> p(vx => x, vy => y)^(1 / set.degree)
end

function polar(set::PolynomialSublevelSetAtOrigin)
    return PolarPolynomialSublevelSetAtOrigin( set.degree, set.p,
                                              set.convexity_proof)
end
function polar(set::PolarPolynomialSublevelSetAtOrigin)
    return PolynomialSublevelSetAtOrigin(set.degree, set.p, set.convexity_proof)
end
