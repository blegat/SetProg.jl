module Sets
using LinearAlgebra

abstract type AbstractSet{T} end

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

using RecipesBase
@recipe function f(aell::AbstractEllipsoid{T}) where T
    @assert dimension(aell) == 2
    ell = convert(Ellipsoid{T}, aell)
    αs = range(0, stop=2π, length=1024)
    ps = [[cos(α), sin(α)] for α in αs]
    r = [sqrt(dot(p, ell.Q * p)) for p in ps]
    seriestype --> :shape
    legend --> false
    ell.c[1] .+ cos.(αs) ./ r, ell.c[2] .+ sin.(αs) ./ r
end
end
