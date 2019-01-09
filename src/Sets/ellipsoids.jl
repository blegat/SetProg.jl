abstract type AbstractEllipsoid{T} <: AbstractSet{T} end
dimension(ell::AbstractEllipsoid) = LinearAlgebra.checksquare(ell.Q)

struct HyperSphere <: AbstractEllipsoid{Bool}
    dim::Int
end

struct Ellipsoid{T} <: AbstractEllipsoid{T}
    Q::Symmetric{T, Matrix{T}}
    center::Vector{T}
end
Ellipsoid(ell::Ellipsoid) = ell
function Polyhedra.project(ell::Ellipsoid, I)
    homogeneous = polar_representation(project(EllipsoidAtOrigin(ell.Q), I))
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


function _HPH(D, d, δ)
    P = [δ d'
         d D]
end

_HPH(ell::EllipsoidAtOrigin) = _HPH(ell.Q, zeros(size(q.Q, 1)), -1.0)

"""
    struct ShiftedEllipsoid{T}
        Q::Symmetric{T, Matrix{T}}
        b::Vector{S}
        β::S
    end

Set ``\\{\\, x \\mid x^\\top Q x + 2 b^\\top x + \\beta \\le 0 \\,\\}``.
"""
struct ShiftedEllipsoid{T} <: AbstractEllipsoid{T}
    Q::Symmetric{T, Matrix{T}}
    b::Vector{T}
    β::T
end
_HPH(q::ShiftedEllipsoid) = _HPH(q.Q, q.b, q.β)
convexity_proof(ell::ShiftedEllipsoid) = ell.Q

function LiftedEllipsoid(qc::HouseDualOf{<:AbstractEllipsoid})
    return LiftedEllipsoid(inv(_HPH(perspective_dual(qc))))
end
function Ellipsoid(qc::HouseDualOf{<:AbstractEllipsoid})
    return Ellipsoid(LiftedEllipsoid(qc))
end
function Polyhedra.project(ell::HouseDualOf{<:AbstractEllipsoid}, I)
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

@recipe function f(aell::PerspectiveDualOrPolarOrNot{<:AbstractEllipsoid};
                   npoints=64)
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
    H = _householder(h[state])
    HPdH = H * Pd * H
    # HPdH is not like a solution what would be obtained by solving the program
    # since the λ computed for unlifting it is maybe not one.
    # Therefore, the S-procedure's λ for the constraints will be different.
    B, b, β, λ = Bbβλ(HPdH)
    ps[state] = y' * H * _HPH(B/λ, b/λ, β/λ) * H * y
    error("TODO: LiftedEllipsoid -> PerspectiveInteriorEllipsoid")
end
