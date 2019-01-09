using Polyhedra
using SumOfSquares

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

function space_variables(set::ConvexPolynomialSublevelSetAtOrigin)
    return variables(set.p)
end
function dimension(set::ConvexPolynomialSublevelSetAtOrigin)
    return nvariables(set.p)
end
function gauge1(set::ConvexPolynomialSublevelSetAtOrigin)
    return set.p
end

convexity_proof(set::ConvexPolynomialSublevelSetAtOrigin) = set.convexity_proof

@recipe function f(set::ConvexPolynomialSublevelSetAtOrigin; npoints=64)
    seriestype --> :shape
    legend --> false
    primal_contour(scaling_function(set), npoints)
end

function scaling_function(set::ConvexPolynomialSublevelSetAtOrigin)
    # We convert the MatPolynomial to a polynomial to avoid having to do the
    # conversion for every substitution.
    p = polynomial(set.p)
    vars = variables(p)
    @assert length(vars) == 2
    vx, vy = vars
    return (x, y) -> p(vx => x, vy => y)^(1 / set.degree)
end

"""
    dual_contour(f::Function, nhalfspaces::Int, T::Type)

Return a polytope of `nhalfspaces` halfspaces defined by normal vectors of
equally spaced angles for the polar of the 1-sublevel set of the homogeneous
function `f(x, y)`.
"""
function dual_contour(f::Function, nhalfspaces::Int, ::Type{T},
                      point::Vector{T} = [0.0, 0.0],
                      x_axis::Vector{T} = [1.0, 0.0],
                      y_axis::Vector{T} = [0.0, 1.0],
                      cone = false) where T
    h = hrep(Polyhedra.HyperPlane{T, Vector{T}}[],
             Polyhedra.HalfSpace{T, Vector{T}}[], d=length(x_axis))
    for α in range(0, stop=2π - 2π/nhalfspaces, length=nhalfspaces)
        ray = x_axis * cos(α) + y_axis * sin(α)
        λ = f(ray...)
        # We have f(ray/λ) = 1 so the halfspace is
        # (point + ray / λ) ⋅ x ≤ 1 for non-cone
        # (point + ray / λ) ⋅ x ≥ 0 for coen
        a = point + ray / λ
        intersect!(h, HalfSpace(cone ? -a : a, cone ? zero(T) : one(T)))
    end
    return polyhedron(h)
end



@recipe function f(set::PolarOf{ConvexPolynomialSublevelSetAtOrigin{T}}; npoints=64) where T
    seriestype --> :shape
    legend --> false
    dual_contour(scaling_function(polar(set)), npoints, T)
end

"""
    struct ConvexPolynomialSet{T} <: AbstractSet{T}
        degree::Int
        q::MatPolynomial{T, DynamicPolynomials.Monomial{true},
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
    q::MatPolynomial{T, DynamicPolynomials.Monomial{true},
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
