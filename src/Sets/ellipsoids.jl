abstract type AbstractEllipsoid{T} <: AbstractSet{T} end
dimension(ell::AbstractEllipsoid) = LinearAlgebra.checksquare(ell.Q)

struct Ellipsoid{T} <: AbstractEllipsoid{T}
    Q::Symmetric{T, Matrix{T}}
    center::Vector{T}
end

struct EllipsoidAtOrigin{T} <: AbstractEllipsoid{T}
    Q::Symmetric{T, Matrix{T}}
end
function Base.convert(::Type{Ellipsoid{T}}, ell::EllipsoidAtOrigin{T}) where T
    Ellipsoid(ell.Q, zeros(T, dimension(ell)))
end

struct PolarEllipsoidAtOrigin{T} <: AbstractEllipsoid{T}
    Q::Symmetric{T, Matrix{T}}
end
function Base.convert(::Type{Ellipsoid{T}}, ell::PolarEllipsoidAtOrigin{T}) where T
    convert(Ellipsoid{T}, convert(EllipsoidAtOrigin{T}, ell))
end
function Base.convert(::Type{EllipsoidAtOrigin{T}},
                 ell::PolarEllipsoidAtOrigin{T}) where T
    EllipsoidAtOrigin(inv(ell.Q))
end

struct LiftedEllipsoid{T}
    P::Matrix{T}
end
dimension(ell::LiftedEllipsoid) = LinearAlgebra.checksquare(ell.P) - 1

function LiftedEllipsoid(ell::Ellipsoid)
    md = ell.Q * ell.c
    δ = ell.c' * md-1
    d = -md
    D = ell.Q
    P = [δ d'
         d D]
    LiftedEllipsoid(P)
end

function Base.convert(::Type{Ellipsoid{T}}, ell::LiftedEllipsoid) where T
    convert(Ellipsoid{T}, Ellipsoid(ell))
end
function Bbβλ(P)
    n = LinearAlgebra.checksquare(P) - 1
    ix = 1 .+ (1:n)
    β = P[1, 1]
    b = P[1, ix]
    B = P[ix, ix]
    λ = dot(b, B \ b) - β
    @assert λ >= 0
    B, b, β, λ
end
function Ellipsoid(ell::LiftedEllipsoid)
    # P is
    # λ * [c'Qc-1  -c'Q
    #         -Qc   Q]
    # Let P be [β b'; b B]
    # We have
    # β = λ c'Qc - λ
    # b = -λ Qc <=> Q^{-1/2}b = -λ Q^{1/2} c
    # hence
    # λ c'Qc = β + λ
    # λ^2 c'Qc = b'Q^{-1}b = λ b'B^{-1}b <=> λ c'Qc = b'B^{-1}b
    # Hence λ = b'B^{-1}b - β
    B, b, β, λ = Bbβλ(ell.P)
    c = -(B \ b)
    Q = B / λ
    Ellipsoid(Symmetric(Q), c)
end

"""
    householder(x)

Householder reflection
```math
I - 2 v v^T / (v^T v)
```
It is symmetric and orthogonal.
"""
function householder(x)
    y = copy(x)
    t = LinearAlgebra.reflector!(y)
    v = [1; y[2:end]]
    I - t * v * v'
end
_householder(h) = householder([1.0; h]) # We add 1, for perspective variable z

function _HPH(D, d, δ, H)
    P = [δ d'
         d D]
    HPH = H * P * H
end

abstract type DualQuadCone{T, S} <: AbstractEllipsoid{T} end

"""
    struct CenterDualQuadCone{T}
        p::DynamicPolynomials.Polynomial{true}
        Q::Symmetric{T, Matrix{T}}
        h::Vector{Float64} # h is an center-like point
        H::Matrix{Float64}
    end

```math
p(y) =
y^\\top H^\\top
\\begin{bmatrix}
  -1 & 0\\
   0 & Q
\\end{bmatrix}
H y
```
where ``y = (z, x)``.
"""
struct CenterDualQuadCone{T, S} <: DualQuadCone{T, S}
    p::DynamicPolynomials.Polynomial{true}
    Q::Symmetric{T, Matrix{T}}
    h::Vector{Float64} # h is an center-like point
    H::Matrix{Float64}
end
function CenterDualQuadCone(Q::Symmetric, y, h::Vector{Float64})
    H = _householder(point)
    p = y' * _HPH(Q, zeros(length(h)), -1.0, H) * y
    CenterDualQuadCone(p, Q, h, H)
end
_HPH(q::CenterDualQuadCone) = _HPH(q.Q, zeros(size(q.Q, 1)), -1.0, q.H)
samecenter(q1::CenterDualQuadCone, q2::CenterDualQuadCone) = q1.h == q2.h

"""
    struct InteriorDualQuadCone{T}
        p::DynamicPolynomials.Polynomial{true}
        Q::Symmetric{T, Matrix{T}}
        b::Vector{S}
        β::S
        h::Vector{Float64} # h is an interior point
        H::Matrix{Float64}
    end

```math
p =
y^\\top H^\\top
\\begin{bmatrix}
  \\beta & b^\\top\\
  b & Q
\\end{bmatrix}
H y
```
"""
struct InteriorDualQuadCone{T, S} <: DualQuadCone{T, S}
    p::DynamicPolynomials.Polynomial{true, S}
    Q::Symmetric{T, Matrix{T}}
    b::Vector{T}
    β::T
    h::Vector{Float64} # h is an interior point
    H::Matrix{Float64}
end
function InteriorDualQuadCone(Q::Symmetric, b::Vector, β, y, h::Vector{Float64})
    H = _householder(h)
    p = y' * _HPH(Q, b, β, H) * y
    InteriorDualQuadCone(p, Q, b, β, h, H)
end
_HPH(q::InteriorDualQuadCone) = _HPH(q.Q, q.b, q.β, q.H)
samecenter(::InteriorDualQuadCone, ::InteriorDualQuadCone) = false

function Base.convert(::Type{LiftedEllipsoid{T}}, qc::DualQuadCone) where T
    LiftedEllipsoid{T}(inv(_HPH(qc)))
end
function Base.convert(::Type{Ellipsoid{T}}, qc::DualQuadCone) where T
    convert(Ellipsoid{T}, convert(LiftedEllipsoid{T}, qc))
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
    x, y = primal_contour((x, y) -> sqrt(x^2 * Q[1, 1] + 2x*y * Q[1, 2] + y^2 * Q[2, 2]),
                          npoints)
    ell.center[1] .+ x, ell.center[2] .+ y
end

function polar(ell::EllipsoidAtOrigin)
    return PolarEllipsoidAtOrigin(ell.Q)
end
function polar(ell::PolarEllipsoidAtOrigin)
    return EllipsoidAtOrigin(ell.Q)
end
