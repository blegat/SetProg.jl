# A^{-1} * S
struct LinearPreImage{S, T, MT <: AbstractMatrix{T}}
    set::S
    A::MT
end

# S + c
struct Translation{S, T, VT <: AbstractVector{T}}
    set::S
    c::VT
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

struct Householder{T, S <: AbstractSet{T}, U} <: AbstractSet{T}
    set::S
    p::DynamicPolynomials.Polynomial{true, U}
    h::Vector{Float64}
    z::SpaceVariable
    x::Vector{SpaceVariable}
end
perspective_gauge0(set::Householder) = set.p
perspective_variable(set::Householder) = set.z
space_variables(set::Householder) = set.x
convexity_proof(set::Householder) = convexity_proof(set.set)

function Polyhedra.project(set::Householder,
                           I)
    project(set, [I])
end
function Polyhedra.project(set::PerspectiveDualOf{Householder{T, S, U}},
                           I::AbstractVector) where {T, S, U}
    J = setdiff(1:dimension(set), I)
    dual = perspective_dual(set)
    p = subs(dual.p,
             dual.x[J] => zeros(T, length(J)))
    proj = Householder(UnknownSet{T}(), p, dual.h[I], dual.z, dual.x[I])
    return perspective_dual(proj)
end

@recipe function f(set::PerspectiveDual{T, <:Householder};
                   npoints=64) where T
    seriestype --> :shape
    legend --> false
    # z is a halfspace of the primal so a ray of the dual
    z = [1.0, 0.0, 0.0]
    h1, h2 = set.set.h
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

function _HPH(set::Householder)
    H = _householder(set.h)
    return H * _HPH(set.set) * H
end

const HouseDualOf{S, T, U} = PerspectiveDualOf{Householder{T, S, U}}
