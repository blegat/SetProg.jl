"""
    primal_contour(f::Function, npoints::Int)

Return `npoints` points with equally spaced angles of the 1-sublevel set of the
homogeneous function `f(x, y)`.
"""
function primal_contour(f::Function, npoints::Int)
    x = Vector{Float64}(undef, npoints)
    y = Vector{Float64}(undef, npoints)
    for (i, α) in enumerate(range(0, stop=2π - 2π/npoints, length=npoints))
        x0 = cos(α)
        y0 = sin(α)
        r = f(x0, y0)
        # f is homogeneous so f(x0/r, y0/r) = 1
        x[i] = x0 / r
        y[i] = y0 / r
    end
    return x, y
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

function Polyhedra.planar_contour(ell::PerspectiveDualOrPolarOrNot{<:AbstractEllipsoid};
                           kws...)
    return Polyhedra.planar_contour(ellipsoid(ell); kws...)
end

function Polyhedra.planar_contour(ell::EllipsoidAtOrigin; npoints=64)
    @assert dimension(ell) == 2
    Q = ell.Q
    return primal_contour((x, y) -> sqrt(x^2 * Q[1, 1] + 2x*y * Q[1, 2] + y^2 * Q[2, 2]),
                          npoints)
end

function Polyhedra.planar_contour(set::Union{PolynomialSublevelSetAtOrigin,
                                             ConvexPolynomialSublevelSetAtOrigin};
                                  npoints=64)
    return primal_contour(scaling_function(set), npoints)
end

function Polyhedra.planar_contour(set::PolarOf{ConvexPolynomialSublevelSetAtOrigin{T}};
                 npoints=64) where T
    return Polyhedra.planar_contour(dual_contour(scaling_function(polar(set)),
                                                 npoints, T))
end

function Polyhedra.planar_contour(set::PerspectiveDual{T, <:Householder};
                           npoints=64) where T
    @assert dimension(set) == 2
    # z is a halfspace of the primal so a ray of the dual
    z = [1.0, 0.0, 0.0]
    h1, h2 = set.set.h
    # a is a ray of the primal so a halfspace of the dual
    a = [1, h1, h2]
    b = [h1, -1, 0]
    @assert abs(dot(a, b)) < 1e-8
    c = [h2 / (1 + h1^2), h1*h2 / (1 + h1^2), -1]
    @assert abs(dot(b, c)) < 1e-8
    @assert abs(dot(a, c)) < 1e-8
    polyhedron = dual_contour(scaling_function(set), npoints, T,
                              z, b, c, true)
    # We fix z to 1.0 and eliminate it, this is cheap for H-rep
    return Polyhedra.planar_contour(fixandeliminate(polyhedron, 1, 1.0))
end

function Polyhedra.planar_contour(t::Translation; kws...)
    @assert dimension(t) == 2
    x, y = Polyhedra.planar_contour(t.set; kws...)
    return x .+ t.c[1], y .+ t.c[2]
end

@recipe function f(set::AbstractSet; npoints=64)
    @assert dimension(set) == 2
    seriestype --> :shape
    legend --> false
    Polyhedra.planar_contour(set; npoints=npoints)
end
