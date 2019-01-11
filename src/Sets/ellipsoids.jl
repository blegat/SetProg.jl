abstract type AbstractEllipsoid{T} <: AbstractSet{T} end
dimension(ell::AbstractEllipsoid) = LinearAlgebra.checksquare(ell.Q)

struct HyperSphere <: AbstractEllipsoid{Bool}
    dim::Int
end

"""
    struct EllipsoidAtOrigin{T} <: AbstractEllipsoid{T}
        Q::Symmetric{T, Matrix{T}}
    end
"""
struct EllipsoidAtOrigin{T} <: AbstractEllipsoid{T}
    Q::Symmetric{T, Matrix{T}}
end
ellipsoid(ell::EllipsoidAtOrigin) = ell
function Polyhedra.project(ell::EllipsoidAtOrigin, I)
    return project(polar_representation(ell), I)
end
convexity_proof(ell::EllipsoidAtOrigin) = ell.Q

function ellipsoid(ell::PolarOf{<:EllipsoidAtOrigin})
    ellipsoid(polar_representation(ell))
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

function perspective_variables(ell::Union{EllipsoidAtOrigin, LiftedEllipsoid})
    return nothing
end
function space_variables(ell::Union{EllipsoidAtOrigin, LiftedEllipsoid})
    return nothing
end

function LiftedEllipsoid(t::Translation{<:EllipsoidAtOrigin})
    ell = t.set
    md = ell.Q * t.c
    δ = t.c' * md-1
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
function ellipsoid(ell::LiftedEllipsoid)
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
    Translation(EllipsoidAtOrigin(Symmetric(Q)), c)
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
function ellipsoid(qc::HouseDualOf{<:AbstractEllipsoid})
    return ellipsoid(LiftedEllipsoid(qc))
end
function Polyhedra.project(ell::HouseDualOf{<:AbstractEllipsoid},
                           I::AbstractVector)
    return project(ellipsoid(ell), I)
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
