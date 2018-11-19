abstract type AbstractEllipsoid{T} <: AbstractSet{T} end
dimension(ell::AbstractEllipsoid) = LinearAlgebra.checksquare(ell.Q)
struct Ellipsoid{T} <: AbstractEllipsoid{T}
    Q::Symmetric{T, Matrix{T}}
    c::Vector{T}
end

struct EllipsoidAtOrigin{T} <: AbstractEllipsoid{T}
    Q::Symmetric{T, Matrix{T}}
end
function convert(::Type{Ellipsoid{T}}, ell::EllipsoidAtOrigin{T}) where T
    Ellipsoid(ell.Q, zeros(T, dimension(ell)))
end

struct PolarEllipsoidAtOrigin{T} <: AbstractEllipsoid{T}
    Q::Symmetric{T, Matrix{T}}
end
function convert(::Type{Ellipsoid{T}}, ell::PolarEllipsoidAtOrigin{T}) where T
    convert(Ellipsoid{T}, convert(EllipsoidAtOrigin{T}, ell))
end
function convert(::Type{EllipsoidAtOrigin{T}},
                 ell::PolarEllipsoidAtOrigin{T}) where T
    EllipsoidAtOrigin(inv(ell.Q))
end

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

@recipe function f(aell::AbstractEllipsoid{T}; npoints=64) where T
    @assert dimension(aell) == 2
    ell = convert(Ellipsoid{T}, aell)
    seriestype --> :shape
    legend --> false
    Q = ell.Q
    primal_contour((x, y) -> sqrt(x^2 * Q[1, 1] + 2x*y * Q[1, 2] + y^2 * Q[2, 2]),
                   npoints)
end

function polar(ell::EllipsoidAtOrigin)
    return PolarEllipsoidAtOrigin(ell.Q)
end
function polar(ell::PolarEllipsoidAtOrigin)
    return EllipsoidAtOrigin(ell.Q)
end
