abstract type AbstractEllipsoid{T} <: AbstractSet{T} end
dimension(ell::AbstractEllipsoid) = LinearAlgebra.checksquare(ell.Q)

struct Ellipsoid{T} <: AbstractEllipsoid{T}
    Q::Symmetric{T, Matrix{T}}
    center::Vector{T}
end
Ellipsoid(ell::Ellipsoid) = ell
function Polyhedra.project(ell::Ellipsoid, I)
    homogeneous = EllipsoidAtOrigin(project(EllipsoidAtOrigin(ell.Q), I))
    return Ellipsoid(homogeneous.Q, ell.center[I])
end

"""
    struct EllipsoidAtOrigin{T} <: AbstractEllipsoid{T}
        Q::Symmetric{T, Matrix{T}}
    end
"""
struct EllipsoidAtOrigin{T} <: AbstractEllipsoid{T}
    Q::Symmetric{T, Matrix{T}}
end
function Ellipsoid(ell::EllipsoidAtOrigin)
    Ellipsoid(ell.Q, zeros(eltype(ell.Q), dimension(ell)))
end
function Polyhedra.project(ell::EllipsoidAtOrigin, I)
    return project(polar_representation(ell), I)
end
convexity_proof(ell::EllipsoidAtOrigin) = ell.Q

function Ellipsoid(ell::PolarOf{<:EllipsoidAtOrigin})
    Ellipsoid(polar_representation(ell))
end
function polar_representation(ell::PolarOf{<:EllipsoidAtOrigin})
    EllipsoidAtOrigin(inv(ell.set.Q))
end
function polar_representation(ell::EllipsoidAtOrigin)
    polar(EllipsoidAtOrigin(inv(ell.Q)))
end

function Polyhedra.project(ell::PolarOf{<:EllipsoidAtOrigin}, I)
    return polar(EllipsoidAtOrigin(Symmetric(ell.set.Q[I, I])))
end

struct LiftedEllipsoid{T}
    P::Matrix{T}
end
dimension(ell::LiftedEllipsoid) = LinearAlgebra.checksquare(ell.P) - 1

function perspective_variables(ell::Union{Ellipsoid, EllipsoidAtOrigin,
                                          LiftedEllipsoid})
    return nothing
end
function space_variables(ell::Union{Ellipsoid, EllipsoidAtOrigin,
                                    LiftedEllipsoid})
    return nothing
end

function LiftedEllipsoid(ell::Ellipsoid)
    md = ell.Q * ell.c
    δ = ell.c' * md-1
    d = -md
    D = ell.Q
    P = [δ d'
         d D]
    LiftedEllipsoid(P)
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

abstract type PerspectiveEllipsoid{T, S} <: AbstractEllipsoid{T} end

# The first variable is the perspective variable z
perspective_variable(ell::PerspectiveEllipsoid) = variables(ell.p)[1]
space_variables(ell::PerspectiveEllipsoid) = variables(ell.p)[2:end]
convexity_proof(ell::PerspectiveEllipsoid) = ell.Q

"""
    struct PerspectiveCenterEllipsoid{T}
        p::DynamicPolynomials.Polynomial{true}
        Q::Symmetric{T, Matrix{T}}
        h::Vector{Float64} # h is an center-like point
        H::Matrix{Float64}
    end

Set ``\\{\\, (z, x) \\mid p(z, x) \\le 0 \\,\\}`` where `p` is a
quadratic forms given by:
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
struct PerspectiveCenterEllipsoid{T, S} <: PerspectiveEllipsoid{T, S}
    p::DynamicPolynomials.Polynomial{true}
    Q::Symmetric{T, Matrix{T}}
    h::Vector{Float64} # h is an center-like point
    H::Matrix{Float64}
end
function PerspectiveCenterEllipsoid(Q::Symmetric, y, h::Vector)
    H = _householder(h)
    p = y' * _HPH(Q, zeros(length(h)), -1.0, H) * y
    PerspectiveCenterEllipsoid(p, Q, h, H)
end
_HPH(q::PerspectiveCenterEllipsoid) = _HPH(q.Q, zeros(size(q.Q, 1)), -1.0, q.H)
samecenter(q1::PerspectiveCenterEllipsoid, q2::PerspectiveCenterEllipsoid) = q1.h == q2.h

"""
    struct PerspectiveInteriorEllipsoid{T}
        p::DynamicPolynomials.Polynomial{true}
        Q::Symmetric{T, Matrix{T}}
        b::Vector{S}
        β::S
        h::Vector{Float64} # h is an interior point
        H::Matrix{Float64}
    end

Set whose dual is ``\\{\\, (z, x) \\mid p(z, x) \\le 0 \\,\\}`` where `p` is a
quadratic forms given by:
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
struct PerspectiveInteriorEllipsoid{T, S} <: PerspectiveEllipsoid{T, S}
    p::DynamicPolynomials.Polynomial{true, S}
    Q::Symmetric{T, Matrix{T}}
    b::Vector{T}
    β::T
    h::Vector{Float64} # h is an interior point
    H::Matrix{Float64}
end
function PerspectiveInteriorEllipsoid(Q::Symmetric, b::Vector, β, y, h::Vector)
    H = _householder(h)
    p = y' * _HPH(Q, b, β, H) * y
    PerspectiveInteriorEllipsoid(p, Q, b, β, h, H)
end
_HPH(q::PerspectiveInteriorEllipsoid) = _HPH(q.Q, q.b, q.β, q.H)
samecenter(::PerspectiveInteriorEllipsoid, ::PerspectiveInteriorEllipsoid) = false

function LiftedEllipsoid(qc::PerspectiveDualOf{<:PerspectiveEllipsoid})
    return LiftedEllipsoid(inv(_HPH(perspective_dual(qc))))
end
function Ellipsoid(qc::PerspectiveDualOf{<:PerspectiveEllipsoid})
    return Ellipsoid(LiftedEllipsoid(qc))
end
function Polyhedra.project(ell::PerspectiveDualOf{<:PerspectiveEllipsoid}, I)
    return project(Ellipsoid(ell), I)
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

@recipe function f(aell::PolarOrNot{<:AbstractEllipsoid}; npoints=64)
    @assert dimension(aell) == 2
    ell = Ellipsoid(aell)
    seriestype --> :shape
    legend --> false
    Q = ell.Q
    x, y = primal_contour((x, y) -> sqrt(x^2 * Q[1, 1] + 2x*y * Q[1, 2] + y^2 * Q[2, 2]),
                          npoints)
    ell.center[1] .+ x, ell.center[2] .+ y
end

function PerspectiveInteriorEllipsoid(ell::LiftedEllipsoid)
    Pd = inv(le.P)
    H = SetProg.Sets._householder(h[state])
    HPdH = H * Pd * H
    # HPdH is not like a solution what would be obtained by solving the program
    # since the λ computed for unlifting it is maybe not one.
    # Therefore, the S-procedure's λ for the constraints will be different.
    B, b, β, λ = Bbβλ(HPdH)
    ps[state] = y' * _HPH(B/λ, b/λ, β/λ, H) * y
    error("TODO: LiftedEllipsoid -> PerspectiveInteriorEllipsoid")
end
