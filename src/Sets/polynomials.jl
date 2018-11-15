using MultivariatePolynomials

abstract type AbstractPolySet{T} <: AbstractSet{T} end

"""
    struct PolarPolynomialSublevelSet{T, P<:AbstractPolynomial{T}}
        degree::Int
        p::P
    end

Set whose polar is ``\\{\\, x \\mid p(x) \\le 1 \\,\\}`` where `p` is a
homogeneous polynomial of degree `degree`.
"""
struct PolarPolynomialSublevelSet{T, P<:AbstractPolynomial{T}}
    degree::Int
    p::P
end

"""
    dual_contour(f::Function, nhalfspaces::Int)

Return a polytope of `nhalfspaces` halfspaces defined by normal vectors of
equally spaced angles for the polar of the 1-sublevel set of the homogeneous
function `f(x, y)`.
"""
function primal_contour(f::Function, npoints::Int, T::Type)
    αs = range(0, stop=2π, length=npoints)
    h = hrep(HyperPlane{T, Vector{T}}[], HalfSpace{T, Vector{T}}[], d=2)
    for (i, α) in enumerate(range(0, stop=2π - 2π/npoints, length=npoints))
        a = cos(α)
        b = sin(α)
        r = f(a, b)
        # f is homogeneous so f(a/r, b/r) = 1 so the halfspace is
        # a*x/r + b*y/r ≤ 1 or equivalently a*x + b*y ≤ r
        intersect!(h, HalfSpace([a, b], r))
    end
    return polyhedron(h)
end

@recipe function f(set::PolarPolynomialSublevelSet{T}) where T
    vars = variables(set)
    @assert length(vars) == 2
    vx, vy = vars
    seriestype --> :shape
    legend --> false
    dual_contour((x, y) -> set.p(vx => x, vy => y)^(1 / set.degree),
                 256)
end


"""
    struct PolarPolynomialSublevelSet{T, P<:AbstractPolynomial{T}}
        degree::Int
        p::P
    end

Set ``\\{\\, x \\mid p(x) \\le 1 \\,\\}`` where `p` is a homogeneous polynomial
of degree `degree`.
"""
struct PolynomialSublevelSet{T, P<:AbstractPolynomial{T}}
    degree::Int
    p::P
end

@recipe function f(set::PolynomialSublevelSet; npoints::Int=64)
    vars = variables(set)
    @assert length(vars) == 2
    vx, vy = vars
    seriestype --> :shape
    legend --> false
    primal_contour((x, y) -> set.p(vx => x, vy => y)^(1 / set.degree),
                   npoints)
end
