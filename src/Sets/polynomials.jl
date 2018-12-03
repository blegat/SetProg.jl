using Polyhedra
using SumOfSquares

"""
    struct DualPolynomialSet{T} <: AbstractSet{T}
        degree::Int
        p::DynamicPolynomials.Polynomial{true, T}
        h::Vector{Float64}
        z::SpaceVariable
        x::Vector{SpaceVariable}
    end

Sets whose perspective-dual is ``\\{\\, (z, x) \\mid p(z, x) \\le 0 \\,\\}``
where `p` is a homogeneous polynomial of degree `degree`.
"""
struct DualPolynomialSet{T} <: AbstractSet{T}
    degree::Int
    p::DynamicPolynomials.Polynomial{true, T}
    h::Vector{Float64}
    z::SpaceVariable
    x::Vector{SpaceVariable}
end

space_variables(set::DualPolynomialSet) = set.x

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

function space_variables(set::Union{ConvexPolynomialSublevelSetAtOrigin,
                                    PolarConvexPolynomialSublevelSetAtOrigin})
    return variables(set.p)
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
    return PolarConvexPolynomialSublevelSetAtOrigin(set.degree, set.p,
                                                    set.convexity_proof)
end
function polar(set::PolarConvexPolynomialSublevelSetAtOrigin)
    return ConvexPolynomialSublevelSetAtOrigin(set.degree, set.p,
                                               set.convexity_proof)
end

"""
    struct DualConvexPolynomialCone{T, U} <: AbstractSet{T}
        degree::Int
        p::MatPolynomial{T, DynamicPolynomials.Monomial{true},
                            DynamicPolynomials.MonomialVector{true}}
        q::DynamicPolynomials.Polynomial{true, U}
    end

Set whose dual is ``\\{\\, (z, x) \\mid p(z, x) \\le 0 \\,\\}`` or
``H \\{\\, (z, x) \\mid q(z, x) \\le z^{\\texttt{degree}} \\,\\}`` where `p` and
`q` are homogeneous polynomials of degree `degree`.
"""
struct DualConvexPolynomialCone{T, U} <: AbstractSet{T}
    degree::Int
    q::MatPolynomial{T, DynamicPolynomials.Monomial{true},
                        DynamicPolynomials.MonomialVector{true}}
    p::DynamicPolynomials.Polynomial{true, U}
    h::Vector{Float64}
    H::Matrix{Float64}
    z::SpaceVariable
    x::Vector{SpaceVariable}
end
function DualConvexPolynomialCone(degree::Integer, q::MatPolynomial, h::Vector,
                                  z, x)
    H = _householder(h)
    y = [z; x]
    p = (q - z^degree)(y => H * y)
    return DualConvexPolynomialCone(degree, q, p, h, H, z, x)
end
dimension(d::DualConvexPolynomialCone) = length(d.x)
space_variables(set::DualConvexPolynomialCone) = set.x
function Polyhedra.project(set::DualConvexPolynomialCone, I)
    project(set, [I])
end
function Polyhedra.project(set::DualConvexPolynomialCone{T},
                           I::AbstractVector) where T
    J = setdiff(1:dimension(set), I)
    DualPolynomialSet(set.degree, subs(set.p, set.x[J] => zeros(T, length(J))),
                      set.h[I], set.z, set.x[I])
end

function scaling_function(set::Union{DualConvexPolynomialCone,
                                     DualPolynomialSet})
    @assert length(set.x) == 2
    vars = [set.z; set.x]
    # z is a halfspace of the primal so a ray of the dual
    z = [1.0, 0.0, 0.0]
    in_set(Δ::Vector) = set.p(vars => z + Δ) < 0
    @assert in_set(zeros(3))
    return (Δz, Δx, Δy) -> begin
        Δ = [Δz, Δx, Δy]
        _in_set(λ::Real) = in_set(Δ * λ)
        λ = 1.0
        while _in_set(λ)
            if λ > 1e10
                error("Error in plotting : the `InteriorPoint` seems to be on the boundary")
            end
            λ *= 2
        end
        λmin = 0.0
        λmax = λ
        # Binary search. Invariant: in_set(λmin) and !in_set(λmax)
        while abs(λmin - λmax) > 1e-8
            λ = (λmin + λmax) / 2
            if _in_set(λ)
                λmin = λ
            else
                λmax = λ
            end
        end
        λ = (λmin + λmax) / 2
        return 1 / λ
    end
end

@recipe function f(set::Union{DualConvexPolynomialCone{T},
                              DualPolynomialSet{T}}; npoints=64) where T
    seriestype --> :shape
    legend --> false
    # z is a halfspace of the primal so a ray of the dual
    z = [1.0, 0.0, 0.0]
    h1, h2 = set.h
    # a is a ray of the primal so a halfspace of the dual
    a = [1, h1, h2]
    b = [h1, -1, 0]
    @assert abs(dot(a, b)) < 1e-8
    c = [h2 * (1 - h1^2) / (1 + h1^2), h1*h2 / (1 + h1^2), -1]
    @assert abs(dot(b, c)) < 1e-8
    @assert abs(dot(a, c)) < 1e-8
    polyhedron = dual_contour(scaling_function(set), npoints, T,
                              z, b, c, true)
    # We fix z to 1.0 and eliminate it, this is cheap for H-rep
    pp = fixandeliminate(polyhedron, 1, 1.0)
    pp
end
